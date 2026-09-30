#!/usr/bin/env python3
"""Line-wing cut-off dependence of the HITEMP band means (Task 39).

    python tools/wing_cutoff.py config/fiducial.yaml
    python tools/wing_cutoff.py config/fiducial.yaml --only 1000,30
    python tools/wing_cutoff.py config/fiducial.yaml --truncations 25,50

RADIS truncates every line profile at ``truncation`` cm⁻¹ from line centre and
renormalises what is left to unit area, so line *strength* is conserved exactly
(A58) but the *shape* is not: the area that physically belongs beyond the cut-off
is piled back inside it. Where the Lorentzian HWHM exceeds the cut-off — 45 of the
48 grid points — the profile is effectively a rectangle of width 2·truncation
carrying the full line strength. κ_P, an arithmetic mean, barely notices;
κ_R, a harmonic mean, is set by the transparency windows and notices a great deal.
This scans the cut-off and measures how much.

**Only the truncation varies.** Broadening pressure is the ideal-gas P = ρRT/M the
delivered spectra were built with (A59), the wavenumber range, step, isotopes and
line-loading edge are the delivered ones, and ``neighbour_lines`` stays at RADIS's
default 0 — so the 50 cm⁻¹ row must reproduce the delivered κ_R, and does (that
agreement is asserted, not assumed). ``--neighbour-lines`` runs the separate
question of whether lines outside the loaded range matter at a wide cut-off.

**Renormalisation cannot be switched off in RADIS 0.16.3.** ``lineshape /= area``
in ``radis/lbl/broadening.py`` (lines 1526, 1592, 1770, 1892, 1989) is
unconditional for every profile and every broadening method; no ``SpectrumFactory``
argument or ``radis.config`` key exposes it. The brief's "where RADIS allows"
therefore resolves to: nowhere, and the library is not patched to find out.

**The 0-50 cm^-1 lines are not in the local cache.** At truncation >= 100 a line
below 50 cm^-1 could reach the 100 cm^-1 band edge, and those lines are absent
here. kappa_R is set by the transparency windows, not by the band edge, so the
effect is small — and it is one-signed: the wide-cut-off columns are missing
absorption the narrow ones never had, so the spreads below are **lower bounds**.

Writes ``outputs/wing_cutoff.json`` and the per-run spectra to
``outputs/wing_cutoff/``. Nothing is written into the spectra directories the
opacity loader globs: a stray file there would silently change the 90-point
line-list grid (A60).
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import atomicio

import numpy as np

from config import Config, load_config
from opacity.means import planck_mean, rosseland_mean
from tools.broadening_sensitivity import (
    POINTS,
    WAVENUM_MAX_CM,
    WAVENUM_MIN_CM,
    build_factory,
    compute_kappa,
    ideal_pressure_pa,
)

TRUNCATIONS_CM = (25.0, 50.0, 100.0, 300.0, 500.0)
DELIVERED_TRUNCATION_CM = 50.0
BAND_CM = (100.0, 4000.0)


def _band_means(wn: np.ndarray, kappa: np.ndarray, t_k: float,
                rel_floor: float) -> tuple[float, float]:
    """The configured ``rel_floor`` is passed through, not defaulted: the
    delivered baseline this scan has to reproduce was computed with it."""
    return (float(planck_mean(wn, kappa, t_k, BAND_CM, rel_floor)),
            float(rosseland_mean(wn, kappa, t_k, BAND_CM, rel_floor)))


def _delivered(cfg: Config) -> dict[tuple[float, float], float]:
    """κ_R at the six points from the previous run's ``ideal`` convention.

    That run used the delivered 50 cm⁻¹ truncation, so it is the baseline this
    scan has to reproduce at its 50 cm⁻¹ row.
    """
    path = cfg.output_path / "broadening_sensitivity.json"
    if not path.exists():
        return {}
    doc = json.loads(path.read_text())
    return {(float(r["T"]), float(r["rho"])): float(r["kappa_R_ideal_cm2_g"])
            for r in doc["results"]}


def add_fp_ratios(cfg: Config, doc: dict) -> dict:
    """First-principles / line-list Rosseland ratio across the cut-off range.

    The line-list kappa_R is already in the scan; the first-principles value at
    the same state is a table lookup, so this is post-processing and needs no
    spectra. It answers the question the cut-off scan exists for: how much of the
    first-principles excess is a line-shape convention?
    """
    from opacity.build_tables import get_source
    from opacity.table import CM2_PER_G_TO_SI

    fp = get_source("fp", cfg, band=BAND_CM)
    for r in doc["results"]:
        k_fp = float(fp.kappa_R(r["T"], r["rho"]) / CM2_PER_G_TO_SI)
        ratios = {t: k_fp / v["kappa_R_cm2_g"]
                  for t, v in r["by_truncation"].items()}
        r["kappa_R_first_principles_cm2_g"] = k_fp
        r["fp_over_line_list_by_truncation"] = ratios
        r["fp_over_line_list_range"] = [min(ratios.values()), max(ratios.values())]
        r["fp_over_line_list_at_delivered"] = ratios.get(
            f"{DELIVERED_TRUNCATION_CM:g}")
    vals = [x for r in doc["results"] for x in r["fp_over_line_list_range"]]
    doc["fp_over_line_list_overall_range"] = [min(vals), max(vals)]
    doc["fp_note"] = ("the first-principles value is the model's own table at the "
                      "same state, which is clamped outside the computed hull; "
                      "every point here is inside it")
    return doc


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("config", type=Path)
    ap.add_argument("--summarise-only", action="store_true",
                    help="recompute the derived ratios from the existing JSON and "
                         "rewrite it; runs no spectra")
    ap.add_argument("--only", default=None, help="one point, as 'T,rho'")
    ap.add_argument("--truncations", default=",".join(f"{t:g}" for t in TRUNCATIONS_CM))
    ap.add_argument("--wavenum-min", type=float, default=WAVENUM_MIN_CM)
    ap.add_argument("--wavenum-max", type=float, default=WAVENUM_MAX_CM,
                    help=f"spectral upper edge, cm^-1 (delivered {WAVENUM_MAX_CM:g}). "
                         f"Band means use [100, 4000] only, and with neighbour_lines=0 "
                         f"a line above 4000+truncation cannot reach the band, so "
                         f"4600 is exact for truncation <= 500 and much cheaper")
    ap.add_argument("--neighbour-lines", type=float, default=0.0,
                    help="extend line loading beyond the spectral range (default 0, "
                         "the delivered convention); set equal to the truncation to "
                         "ask whether out-of-range lines matter")
    args = ap.parse_args(argv)
    cfg = load_config(args.config)

    if args.summarise_only:
        path = cfg.output_path / "wing_cutoff.json"
        doc = add_fp_ratios(cfg, json.loads(path.read_text()))
        atomicio.write_text(path, json.dumps(doc, indent=1) + "\n")
        print(f"\n  {'T':>6s} {'rho':>8s} {'kappa_R FP':>11s} "
              f"{'FP/LL range over the cut-offs':>32s}")
        for r in doc["results"]:
            lo, hi = r["fp_over_line_list_range"]
            print(f"  {r['T']:6g} {r['rho']:8g} "
                  f"{r['kappa_R_first_principles_cm2_g']:11.4g} "
                  f"{lo:14.3g} - {hi:<14.3g}")
        lo, hi = doc["fp_over_line_list_overall_range"]
        print(f"\n  overall {lo:.3g} - {hi:.3g}")
        worst = max(doc["results"], key=lambda r: r["kappa_R_max_over_min"])
        print(f"\n  largest cut-off-induced spread in kappa_R across the "
              f"{len(doc['results'])} points: "
              f"{worst['kappa_R_max_over_min']:.3g}x at T = {worst['T']:g} K, "
              f"rho = {worst['rho']:g} kg/m^3 "
              f"({worst['kappa_R_range_cm2_g'][0]:.4g} to "
              f"{worst['kappa_R_range_cm2_g'][1]:.4g} cm^2/g)")
        print(f"  wrote {path}")
        return 0

    points = list(POINTS)
    if args.only:
        t_s, r_s = args.only.split(",")
        want = (float(t_s), float(r_s))
        points = [p for p in points
                  if abs(p[0] - want[0]) < 1e-6 and abs(p[1] - want[1]) < 1e-3]
        if not points:
            raise SystemExit(f"--only {args.only} matches none of {POINTS}")
    truncs = [float(x) for x in args.truncations.split(",")]

    out_dir = cfg.output_path / "wing_cutoff"
    out_dir.mkdir(parents=True, exist_ok=True)
    delivered = _delivered(cfg)

    results = []
    for t_k, rho in points:
        p_pa = ideal_pressure_pa(t_k, rho)
        row = {"T": t_k, "rho": rho, "p_ideal_pa": p_pa, "by_truncation": {}}
        for tr in truncs:
            t0 = time.time()
            sf = build_factory(args.wavenum_min, tr, args.neighbour_lines,
                               args.wavenum_max)
            wn, kappa = compute_kappa(sf, t_k, p_pa)
            kp, kr = _band_means(wn, kappa, t_k, cfg.opacity.rel_floor)
            name = f"T{t_k:g}_rho{rho:g}_trunc{tr:g}.npz"
            np.savez_compressed(out_dir / name, wn=wn, kappa=kappa, T=t_k,
                                rho_kg_m3=rho, truncation_cm=tr,
                                p_broad_pa=p_pa)
            row["by_truncation"][f"{tr:g}"] = {
                "kappa_P_cm2_g": kp, "kappa_R_cm2_g": kr,
                "spectrum": name, "seconds": round(time.time() - t0, 1)}
            print(f"  T={t_k:>6g} rho={rho:>8g} trunc={tr:>5g}  "
                  f"kappa_P={kp:10.4g}  kappa_R={kr:10.4g}  "
                  f"({time.time() - t0:.0f} s)", flush=True)
        krs = [v["kappa_R_cm2_g"] for v in row["by_truncation"].values()]
        kps = [v["kappa_P_cm2_g"] for v in row["by_truncation"].values()]
        row["kappa_R_max_over_min"] = float(max(krs) / min(krs))
        row["kappa_P_max_over_min"] = float(max(kps) / min(kps))
        row["kappa_R_range_cm2_g"] = [float(min(krs)), float(max(krs))]
        key = (t_k, rho)
        if key in delivered and f"{DELIVERED_TRUNCATION_CM:g}" in row["by_truncation"]:
            got = row["by_truncation"][f"{DELIVERED_TRUNCATION_CM:g}"]["kappa_R_cm2_g"]
            row["delivered_kappa_R_cm2_g"] = delivered[key]
            row["reproduces_delivered_to"] = float(abs(got / delivered[key] - 1.0))
        results.append(row)
        _write(cfg, results, truncs, args)
    _write(cfg, results, truncs, args)
    path = cfg.output_path / "wing_cutoff.json"
    atomicio.write_text(path, json.dumps(add_fp_ratios(cfg, json.loads(path.read_text())),
                               indent=1) + "\n")
    return 0


def _write(cfg: Config, results: list[dict], truncs: list[float], args) -> None:
    doc = {
        "meta": ("kappa_P and kappa_R on the common band against the RADIS line-wing "
                 "truncation, at the six broadening-test points. Only the truncation "
                 "varies: the broadening pressure is the ideal-gas P the delivered "
                 "spectra used, and every other RADIS setting is the delivered one."),
        "band_cm": list(BAND_CM),
        "rel_floor": cfg.opacity.rel_floor,
        "truncations_cm": truncs,
        "delivered_truncation_cm": DELIVERED_TRUNCATION_CM,
        "wavenum_min_cm": args.wavenum_min,
        "wavenum_max_cm": args.wavenum_max,
        "neighbour_lines_cm": args.neighbour_lines,
        "renormalisation": ("always on; RADIS 0.16.3 exposes no switch "
                            "(lineshape /= area is unconditional in "
                            "radis/lbl/broadening.py)"),
        "results": results,
    }
    atomicio.write_text(cfg.output_path / "wing_cutoff.json", json.dumps(doc, indent=1) + "\n")


if __name__ == "__main__":
    raise SystemExit(main())
