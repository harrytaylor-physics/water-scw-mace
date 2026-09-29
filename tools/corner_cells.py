#!/usr/bin/env python3
"""Merge the four corner cells into the first-principles grid (Tasks 42-43).

    python tools/corner_cells.py config/fiducial.yaml --before   # snapshot first
    python tools/corner_cells.py config/fiducial.yaml --check    # Task 42
    python tools/corner_cells.py config/fiducial.yaml --stamp    # Task 43 merge

Four new cells at (700, 10), (700, 30), (700, 60) and (800, 30) fill the gap
between the cold-dense rows (650/800 K at rho >= 100) and the dilute rows
(>= 1000 K at 10-60), which is exactly where six of the 108 line-list boundaries
fell outside the hull and were clamped.

**The checks come before the marker, not after.** ``units_marker`` is what stops
a legacy factor-1000 spectrum loading (A50), so an unstamped file is refused by
``opacity.loaders.load_spectrum`` by design. ``--check`` therefore reads the npz
directly, the way ``tools/text_checks._load`` does, computes band-limited means
itself, and only ``--stamp`` writes the marker — and only if every check passed.

The embedded ``planck``/``rosseland`` scalars are **not** used: they are computed
over each file's own full range (1.1-6000 cm^-1), not the common band, and A4
records that band-limited means are always recomputed from the spectra.
"""

from __future__ import annotations

import argparse
import csv
import json
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np

from config import Config, load_config

SRC_DIR = Path("~/subneptune_data/fp_corner").expanduser()
BAND_CM = (100.0, 4000.0)
# (T, rho) and the source subdirectory each arrives in.
CELLS = (
    (700.0, 10.0, "T700_rho10_05fs"),
    (700.0, 30.0, "T700_rho30_05fs"),
    (700.0, 60.0, "T700_rho60_05fs"),
    (800.0, 30.0, "T800_rho30_05fs"),
)
# The continuity check: a cell that differs from its 1000 K neighbour at the same
# density by more than this in kappa_P is the signature of the A50 units error.
MAX_NEIGHBOUR_FACTOR = 10.0
UNITS_MARKER_VALUE = "cm2_per_g:A74"
SNAPSHOT_DIR = "pre52"
SNAPSHOT_FILES = ("manuscript_numbers.json", "crossover_v3.json",
                  "S4_domain_scan_extended.csv", "S8_rcb_flip.csv",
                  "rho_scaling.json", "text_checks.json")


def _read_npz(path: Path) -> dict:
    with np.load(path) as z:
        out = {}
        for k in z.files:
            v = z[k]
            if v.shape == ():
                out[k] = float(v) if v.dtype.kind in "fiub" else str(v)
            else:
                out[k] = np.asarray(v, float)
        out["_keys"] = list(z.files)
    return out


def _band_means(cfg: Config, wn, kappa, t_k):
    from opacity.means import planck_mean, rosseland_mean

    order = np.argsort(wn)
    wn, kappa = wn[order], kappa[order]
    return (float(planck_mean(wn, kappa, t_k, BAND_CM, cfg.opacity.rel_floor)),
            float(rosseland_mean(wn, kappa, t_k, BAND_CM, cfg.opacity.rel_floor)))


def _existing(cfg: Config, t_k: float, rho: float):
    """Band means of an existing production cell, or None."""
    p = cfg.fp_opacity_path / f"T{t_k:g}_rho{rho:g}_opacity.npz"
    if not p.exists():
        return None
    d = _read_npz(p)
    return _band_means(cfg, d["wn"], d["kappa"], d["T"])


def model_state(cfg: Config, t_k: float, rho: float) -> dict:
    """What the current table returns at (T, rho), and whether it is clamped."""
    from opacity.build_tables import get_source
    from opacity.table import CM2_PER_G_TO_SI

    src = get_source("fp", cfg, band=BAND_CM)
    clamped = bool(np.asarray(src.table.clamped_mask(t_k, rho)).ravel()[0])
    kp = float(src.table.kappa_planck(t_k, rho) / CM2_PER_G_TO_SI)
    kr = float(src.table.kappa_rosseland(t_k, rho) / CM2_PER_G_TO_SI)
    out = {"kappa_P_cm2_g": kp, "kappa_R_cm2_g": kr,
           "provenance": "clamped" if clamped else "interpolated",
           "n_points": int(src.table.n_points)}
    if clamped:
        tg = np.asarray(src.table.temperature_k, float)
        rg = np.asarray(src.table.rho_si, float)
        lin = src.table._lin
        scale, offset = np.asarray(lin.scale), np.asarray(lin.offset)
        q = np.array([(t_k - offset[0]) / scale[0],
                      (np.log10(rho) - offset[1]) / scale[1]])
        g = np.column_stack([(tg - offset[0]) / scale[0],
                             (np.log10(rg) - offset[1]) / scale[1]])
        i = int(np.argmin(((g - q) ** 2).sum(axis=1)))
        out["clamped_to"] = {"T_K": float(tg[i]), "rho_kg_m3": float(rg[i])}
    return out


def triangle(cfg: Config, t_k: float, rho: float, source: str = "fp"):
    """The interpolation simplex's vertices at (T, rho), or None if clamped."""
    from opacity.build_tables import get_source

    src = get_source(source, cfg, band=BAND_CM)
    tab = src.table
    lin = tab._lin
    scale, offset = np.asarray(lin.scale), np.asarray(lin.offset)
    q = np.array([[(t_k - offset[0]) / scale[0],
                   (np.log10(rho) - offset[1]) / scale[1]]])
    s = int(lin.tri.find_simplex(q)[0])
    if s < 0:
        return None
    verts = lin.tri.simplices[s]
    tg = np.asarray(tab.temperature_k, float)
    rg = np.asarray(tab.rho_si, float)
    return sorted([[float(tg[v]), float(rg[v])] for v in verts])


def snapshot(cfg: Config) -> dict:
    """Everything the 'before' half of the summary needs, taken while the table
    still has 48 points. outputs/ is regenerated in place, so without this the
    before state is gone — it is gitignored, so git holds no copy of it."""
    out_dir = cfg.output_path / SNAPSHOT_DIR
    out_dir.mkdir(parents=True, exist_ok=True)
    for name in SNAPSHOT_FILES:
        p = cfg.output_path / name
        if p.exists():
            shutil.copy2(p, out_dir / name)
    grid = cfg.fp_opacity_path / "grid.json"
    if grid.exists():
        shutil.copy2(grid, out_dir / "grid.json")

    from figures.ed_items import boundary_table

    doc = {"meta": ("the 48-point state, taken before the four corner cells were "
                    "merged, so the summary can quote before and after"),
           "n_points": int(model_state(cfg, 1000.0, 100.0)["n_points"]),
           "new_cells": {}, "boundaries": {}, "fiducial": {}}
    for t_k, rho, _ in CELLS:
        doc["new_cells"][f"{t_k:g}_{rho:g}"] = model_state(cfg, t_k, rho)

    bnds = boundary_table(cfg)
    doc["counts_by_category"] = {}
    for b in bnds:
        doc["counts_by_category"][b["category"]] = \
            doc["counts_by_category"].get(b["category"], 0) + 1
    between = [b for b in bnds if b["category"] == "between computed points"]
    doc["between_computed_points"] = [
        {**{k: b[k] for k in ("T", "rho", "P", "t_eq_k", "t_int_k",
                              "mass_core_mearth", "wmf")},
         **model_state(cfg, b["T"], b["rho"]),
         "triangle": triangle(cfg, b["T"], b["rho"])}
        for b in between]

    flip_path = cfg.output_path / "S8_rcb_flip.csv"
    if flip_path.exists():
        rows = list(csv.DictReader(flip_path.open()))
        doc["strongest"] = [{
            "case": r["case"], "T": float(r["t_k"]), "rho": float(r["rho_si"]),
            "kappa_r_fp_cm2_g": float(r["kappa_r_fp_cm2_g"]),
            "kappa_r_hitemp_cm2_g": float(r["kappa_r_hitemp_cm2_g"]),
            "fp_over_hitemp": float(r["fp_over_hitemp"]),
            "threshold_factor": float(r["threshold_factor"]),
            "fp_clamped_at_rcb": r["fp_clamped_at_rcb"] == "True",
            "triangle": triangle(cfg, float(r["t_k"]), float(r["rho_si"])),
        } for r in rows]
    (out_dir / "state.json").write_text(json.dumps(doc, indent=1) + "\n")
    return doc


def check(cfg: Config) -> dict:
    """Task 42. Returns the table; ``ok`` is False if any cell fails."""
    rows, ok = [], True
    for t_k, rho, sub in CELLS:
        path = SRC_DIR / sub / "opacity_v3.npz"
        if not path.exists():
            rows.append({"T": t_k, "rho": rho, "error": f"{path} missing"})
            ok = False
            continue
        d = _read_npz(path)
        kp, kr = _band_means(cfg, d["wn"], d["kappa"], d["T"])
        r = {"T": t_k, "rho": rho, "file": str(path), "keys": d["_keys"],
             "kappa_P_cm2_g": kp, "kappa_R_cm2_g": kr,
             "embedded_full_range": {"planck": d["planck"],
                                     "rosseland": d["rosseland"]},
             "N": int(d["N"]), "dt_fs": d["dt_fs"], "smooth_cm": d["smooth_cm"],
             "has_units_marker": "units_fixed" in d["_keys"],
             "wn_span_cm": [float(d["wn"].min()), float(d["wn"].max())]}
        # units sanity, independent of the neighbour test: kappa = alpha/rho and
        # kappa/sigma is the same constant as every production file.
        r["alpha_over_rho_matches_kappa"] = bool(np.allclose(
            d["alpha"] / (rho / 1000.0), d["kappa"], rtol=1e-9))
        r["kappa_positive_finite"] = bool(np.all(np.isfinite(d["kappa"]))
                                          and np.all(d["kappa"] > 0))
        r["covers_band"] = bool(d["wn"].min() <= BAND_CM[0]
                                and d["wn"].max() >= BAND_CM[1])
        ref = _read_npz(cfg.fp_opacity_path / "T1000_rho10_opacity.npz")
        r["kappa_over_sigma_matches_production"] = bool(np.isclose(
            float(np.median(d["kappa"] / d["sigma"])),
            float(np.median(ref["kappa"] / ref["sigma"])), rtol=1e-9))

        # (i) the same density at 1000 K, and (700, 30) for the 800 K cell
        neigh = {}
        n1000 = _existing(cfg, 1000.0, rho)
        if n1000:
            neigh["1000K_same_rho"] = {
                "kappa_P_cm2_g": n1000[0], "kappa_R_cm2_g": n1000[1],
                "ratio_P": kp / n1000[0], "ratio_R": kr / n1000[1]}
        r["neighbours"] = neigh
        # (ii) what the current model returns there
        r["model_before"] = model_state(cfg, t_k, rho)
        r["model_before"]["ratio_P_new_over_model"] = \
            kp / r["model_before"]["kappa_P_cm2_g"]
        r["model_before"]["ratio_R_new_over_model"] = \
            kr / r["model_before"]["kappa_R_cm2_g"]

        fac = (max(r["neighbours"]["1000K_same_rho"]["ratio_P"],
                   1.0 / r["neighbours"]["1000K_same_rho"]["ratio_P"])
               if neigh else float("inf"))
        r["neighbour_factor_P"] = fac
        r["passes"] = bool(fac <= MAX_NEIGHBOUR_FACTOR
                           and r["alpha_over_rho_matches_kappa"]
                           and r["kappa_positive_finite"] and r["covers_band"]
                           and r["kappa_over_sigma_matches_production"])
        ok = ok and r["passes"]
        rows.append(r)

    # the extra comparison the brief asks for at (800, 30)
    for r in rows:
        if r.get("T") == 800.0 and r.get("rho") == 30.0:
            other = next((x for x in rows if x.get("T") == 700.0
                          and x.get("rho") == 30.0), None)
            if other:
                r["neighbours"]["700K_same_rho"] = {
                    "kappa_P_cm2_g": other["kappa_P_cm2_g"],
                    "kappa_R_cm2_g": other["kappa_R_cm2_g"],
                    "ratio_P": r["kappa_P_cm2_g"] / other["kappa_P_cm2_g"],
                    "ratio_R": r["kappa_R_cm2_g"] / other["kappa_R_cm2_g"]}
    return {"ok": ok, "rows": rows, "band_cm": list(BAND_CM),
            "max_neighbour_factor": MAX_NEIGHBOUR_FACTOR}


def stamp(cfg: Config, checked: dict) -> dict:
    """Copy the four cells into fp_opacity/ with the units marker, and rewrite
    grid.json. Refuses unless every check passed."""
    if not checked["ok"]:
        raise SystemExit("checks did not pass; nothing stamped")
    written = []
    for t_k, rho, sub in CELLS:
        src = SRC_DIR / sub / "opacity_v3.npz"
        dst = cfg.fp_opacity_path / f"T{t_k:g}_rho{rho:g}_opacity.npz"
        with np.load(src) as z:
            payload = {k: z[k] for k in z.files}
        payload["units_fixed"] = np.array(UNITS_MARKER_VALUE)
        np.savez_compressed(dst, **payload)
        written.append(str(dst))
    paths = sorted(cfg.fp_opacity_path.glob("T*_opacity.npz"))
    grid = []
    for p in paths:
        with np.load(p) as z:
            t, r = float(z["T"]), float(z["rho_kg_m3"])
        grid.append({"T": int(t) if t == int(t) else t,
                     "rho": int(r) if r == int(r) else r})
    grid.sort(key=lambda g: (g["T"], g["rho"]))
    (cfg.fp_opacity_path / "grid.json").write_text(json.dumps(grid, indent=1) + "\n")
    return {"written": written, "n_points": len(grid),
            "grid_json": str(cfg.fp_opacity_path / "grid.json")}


def crossover_rows(cfg: Config) -> dict:
    """The four new crossover rows, with the line list taken two ways.

    ``crossover_v3_rows`` reads the line list off its own triangulated table at
    every one of the 52 points — that is the existing convention and the one ED
    Table 1 reports, so it is not changed here for four rows. The brief asks for
    the line list *at exactly these four states*, so the direct spectra are
    computed as well and the two are printed side by side. Three of the four sit
    on the 90-point line-list grid, so for those the comparison is a free check
    that the table interpolation is doing what it claims.
    """
    from opacity.crossover import crossover_v3_rows
    from opacity.mtckd import apply_continuum
    from tools.text_checks import _load

    rows = {(round(r["T"], 6), round(r["rho"], 6)): r
            for r in crossover_v3_rows(cfg, BAND_CM)}
    out = {"band_cm": list(BAND_CM), "rows": [],
           "line_list_convention": (
               "table: the triangulated line-list table at the point, as every "
               "one of the 52 crossover rows uses; direct: the line-list spectrum "
               "computed at exactly that state")}
    for t_k, rho, _ in CELLS:
        row = dict(rows[(round(t_k, 6), round(rho, 6))])
        direct = {}
        try:
            wn, k = _load(cfg, "hitemp", t_k, rho)
            kp, kr = _band_means(cfg, wn, k, t_k)
            direct["kappa_P_line_cm2_g"] = kp
            direct["kappa_R_line_cm2_g"] = kr
            kc = apply_continuum("mtckd_self", wn, k, t_k, rho, cfg.mtckd_ref_path)
            kpc, krc = _band_means(cfg, wn, kc, t_k)
            direct["kappa_P_line_mtckd_cm2_g"] = kpc
            direct["kappa_R_line_mtckd_cm2_g"] = krc
            direct["ratio_planck_bare"] = row["kappa_P_fp_cm2_g"] / kp
            direct["ratio_rosseland_bare"] = row["kappa_R_fp_cm2_g"] / kr
            direct["ratio_planck_mtckd"] = row["kappa_P_fp_cm2_g"] / kpc
            direct["ratio_rosseland_mtckd"] = row["kappa_R_fp_cm2_g"] / krc
            direct["mtckd_share_kappa_R"] = 1.0 - kr / krc
            direct["table_over_direct_kappa_R"] = row["kappa_R_line_cm2_g"] / kr
            direct["table_over_direct_kappa_P"] = row["kappa_P_line_cm2_g"] / kp
        except FileNotFoundError as exc:  # noqa: BLE001
            direct["error"] = str(exc).split("\n")[0]
        row["line_list_direct"] = direct
        out["rows"].append(row)
    return out


def _print_crossover(doc: dict) -> None:
    print(f"\n  TASK 43 — the four new crossover rows "
          f"(band {doc['band_cm'][0]:g}-{doc['band_cm'][1]:g} cm^-1)\n")
    hdr = (f"  {'T':>5s} {'rho':>5s} {'kappa_P FP':>10s} {'kappa_R FP':>10s} "
           f"{'P/bare':>8s} {'R/bare':>8s} {'P/+MT':>8s} {'R/+MT':>8s} "
           f"{'MT share':>9s} {'LL':>6s}")
    print(hdr); print("  " + "-" * (len(hdr) - 2))
    for r in doc["rows"]:
        print(f"  {r['T']:5g} {r['rho']:5g} {r['kappa_P_fp_cm2_g']:10.4g} "
              f"{r['kappa_R_fp_cm2_g']:10.4g} {r['ratio_planck_bare']:8.3f} "
              f"{r['ratio_rosseland_bare']:8.3f} {r['ratio_planck_mtckd']:8.3f} "
              f"{r['ratio_rosseland_mtckd']:8.3f} "
              f"{r['mtckd_share_kappa_R']:9.4f} "
              f"{'grid' if r['ll_on_grid'] else 'interp':>6s}")
    print(f"\n  the same rows with the line list taken directly at the state:")
    for r in doc["rows"]:
        d = r["line_list_direct"]
        if "error" in d:
            print(f"    ({r['T']:g}, {r['rho']:g})  {d['error']}")
            continue
        print(f"    ({r['T']:g}, {r['rho']:g})  P/bare {d['ratio_planck_bare']:7.3f} "
              f"R/bare {d['ratio_rosseland_bare']:7.3f} "
              f"P/+MT {d['ratio_planck_mtckd']:7.3f} "
              f"R/+MT {d['ratio_rosseland_mtckd']:7.3f} "
              f"MT share {d['mtckd_share_kappa_R']:7.4f}   "
              f"table/direct kappa_R x{d['table_over_direct_kappa_R']:.4f}")


def _print_check(doc: dict) -> None:
    print(f"\n  TASK 42 — four corner cells, band {BAND_CM[0]:g}-{BAND_CM[1]:g} cm^-1\n")
    hdr = (f"  {'T':>5s} {'rho':>5s} {'kappa_P':>9s} {'kappa_R':>9s} "
           f"{'vs 1000 K':>18s} {'model before':>24s} {'pass':>5s}")
    print(hdr); print("  " + "-" * (len(hdr) - 2))
    for r in doc["rows"]:
        if "error" in r:
            print(f"  {r['T']:5g} {r['rho']:5g}  {r['error']}")
            continue
        n = r["neighbours"].get("1000K_same_rho")
        nb = (f"P x{n['ratio_P']:.2f} R x{n['ratio_R']:.2f}" if n else "n/a")
        m = r["model_before"]
        mb = (f"{m['provenance'][:6]} P x{m['ratio_P_new_over_model']:.2f} "
              f"R x{m['ratio_R_new_over_model']:.2f}")
        print(f"  {r['T']:5g} {r['rho']:5g} {r['kappa_P_cm2_g']:9.4g} "
              f"{r['kappa_R_cm2_g']:9.4g} {nb:>18s} {mb:>24s} "
              f"{'yes' if r['passes'] else 'NO':>5s}")
    for r in doc["rows"]:
        extra = r.get("neighbours", {}).get("700K_same_rho")
        if extra:
            print(f"\n  (800, 30) against (700, 30): kappa_P x{extra['ratio_P']:.3f}, "
                  f"kappa_R x{extra['ratio_R']:.3f}")
    for r in doc["rows"]:
        if "error" in r:
            continue
        m = r["model_before"]
        clamp = (f" (clamped to T={m['clamped_to']['T_K']:g}, "
                 f"rho={m['clamped_to']['rho_kg_m3']:g})" if "clamped_to" in m else "")
        print(f"    ({r['T']:g}, {r['rho']:g}) the 48-point model returned "
              f"kappa_P {m['kappa_P_cm2_g']:.4g}, kappa_R {m['kappa_R_cm2_g']:.4g}, "
              f"{m['provenance']}{clamp}")
    print(f"\n  units marker present already: "
          f"{[r.get('has_units_marker') for r in doc['rows']]}")
    print(f"  every cell within a factor {MAX_NEIGHBOUR_FACTOR:g} of its 1000 K "
          f"neighbour: {all(r.get('passes') for r in doc['rows'])}")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("config", type=Path)
    ap.add_argument("--before", action="store_true")
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--stamp", action="store_true")
    ap.add_argument("--crossover", action="store_true",
                    help="print the four new crossover rows (after merging)")
    args = ap.parse_args(argv)
    cfg = load_config(args.config)

    if args.crossover:
        doc = crossover_rows(cfg)
        (cfg.output_path / "corner_cells_crossover.json").write_text(
            json.dumps(doc, indent=1) + "\n")
        _print_crossover(doc)
        return 0

    if args.before:
        doc = snapshot(cfg)
        print(f"\n  snapshot of the {doc['n_points']}-point state -> "
              f"{cfg.output_path / SNAPSHOT_DIR}")
        print(f"    categories: {doc['counts_by_category']}")
        print(f"    between computed points: "
              f"{len(doc['between_computed_points'])} boundaries")
        return 0

    doc = check(cfg)
    (cfg.output_path / "corner_cells_check.json").write_text(
        json.dumps(doc, indent=1, default=str) + "\n")
    _print_check(doc)
    if args.stamp:
        if not doc["ok"]:
            print("\n  REFUSED: a cell failed its checks; nothing was stamped.")
            return 1
        res = stamp(cfg, doc)
        print(f"\n  stamped {len(res['written'])} cells; grid is now "
              f"{res['n_points']} points")
    return 0 if doc["ok"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
