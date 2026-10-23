#!/usr/bin/env python3
"""Task 22 / A61: where the dilute FP intensity excess lives, band by band.

    python tools/band_diagnostic.py config/fiducial.yaml

At ρ = 10-60 kg/m³ the first-principles Planck mean exceeds bare HITEMP by
1.32-1.53× while the Rosseland ratio falls to 1.08. A band-integrated mean cannot
say whether that excess is spread over the whole spectrum or concentrated in one
band, and the answer attributes it to a different part of the calculation. This
splits it at the six dilute points (1000, 1500, 2000 K) × (10, 30 kg/m³).

Both sources are read as delivered — FP from ``fp_opacity/``, bare HITEMP from
``hitemp_opacity/`` at its own ideal-gas broadening pressure — and integrated
over identical band limits on their own wavenumber grids (see ``opacity/bands.py``
for the integrand and the band definitions).

Writes ``outputs/band_diagnostic.json``.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import atomicio

import numpy as np

from config import Config, load_config
from opacity.bands import BANDS, band_absorption
from opacity.registry import get_spec

POINTS = [(t, rho) for t in (1000.0, 1500.0, 2000.0) for rho in (10.0, 30.0)]


def _load(cfg: Config, source: str, t_k: float, rho: float) -> tuple[np.ndarray, np.ndarray]:
    """(wn, kappa) for one delivered spectrum, matched by its stored (T, ρ)."""
    spec = get_spec(source)
    for path in sorted(cfg.spectra_dir(spec.subdir).glob(spec.glob)):
        with np.load(path) as d:
            if abs(float(d["T"]) - t_k) > 0.5:
                continue
            if abs(float(d["rho_kg_m3"]) - rho) > max(1e-3, 1e-3 * rho):
                continue
            return np.asarray(d["wn"], float), np.asarray(d["kappa"], float)
    raise FileNotFoundError(
        f"no {source} spectrum at T={t_k:g} K, rho={rho:g} kg/m^3 in "
        f"{cfg.spectra_dir(spec.subdir)}"
    )


def rows(cfg: Config, fp_name: str = "fp", line_list: str = "hitemp") -> list[dict]:
    out = []
    for t_k, rho in POINTS:
        wn_fp, k_fp = _load(cfg, fp_name, t_k, rho)
        wn_ll, k_ll = _load(cfg, line_list, t_k, rho)
        a_fp = band_absorption(wn_fp, k_fp, t_k)
        a_ll = band_absorption(wn_ll, k_ll, t_k)
        out.append({
            "T": t_k,
            "rho": rho,
            "A_fp": a_fp,
            "A_line": a_ll,
            "ratio": {b: a_fp[b] / a_ll[b] for b in BANDS},
            "stretch_over_bend_fp": a_fp["stretch"] / a_fp["bend"],
            "stretch_over_bend_line": a_ll["stretch"] / a_ll["bend"],
        })
    return out


def _verdict(rows: list[dict]) -> tuple[str, bool]:
    """The reading the diagnostic exists to produce, plus whether it is monotone.

    The brief offered three outcomes — uniform (dipole surface), stretch-confined
    (dynamics + harmonic quantum correction), or window-dominated (collision
    induced). The data can also show a *tilt*: a ratio that rises monotonically
    with band centre and crosses unity partway. That is a stronger and more
    specific statement than "stretch-confined", so it is detected separately
    rather than collapsed into one of the three.
    """
    core = ("libration", "bend", "stretch")
    ordered = ("libration", "bend", "window", "stretch")  # ascending band centre
    spread = max(
        max(r["ratio"][b] for b in core) / min(r["ratio"][b] for b in core) for r in rows
    )
    stretch_lead = min(
        r["ratio"]["stretch"] / max(r["ratio"]["libration"], r["ratio"]["bend"])
        for r in rows
    )
    monotone = all(
        all(r["ratio"][a] < r["ratio"][b] for a, b in zip(ordered, ordered[1:]))
        for r in rows
    )
    if spread < 1.25:
        return ("near-uniform across libration/bend/stretch -> an overall intensity "
                "scale, i.e. the DIPOLE SURFACE"), monotone
    if stretch_lead > 1.25:
        base = ("concentrated in the stretch -> the CLASSICAL DYNAMICS + HARMONIC "
                "QUANTUM CORRECTION, which act on the high-frequency modes")
        if monotone:
            base += (". Stronger than that: the ratio rises MONOTONICALLY with band "
                     "centre at every point and crosses unity between the bend and "
                     "the window, so it is a frequency-dependent tilt -- a deficit "
                     "at low frequency as well as an excess at high -- not an "
                     "excess bolted onto an otherwise matching spectrum")
        return base, monotone
    return ("neither uniform nor stretch-confined -> no single attribution; report "
            "the per-band spread and do not assign a cause"), monotone


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("config", type=Path)
    ap.add_argument("--line-list", default="hitemp")
    args = ap.parse_args(argv)
    cfg = load_config(args.config)

    data = rows(cfg, line_list=args.line_list)
    band_names = list(BANDS)

    print(f"\n  TASK 22 — band-integrated absorption A = ∫ ν(1−e^(−hcν/kT)) κ(ν) dν,")
    print(f"  first-principles / {args.line_list.upper()} (bare, as delivered), at the "
          f"six dilute points.")
    print(f"  Bands [cm⁻¹]: " + ", ".join(f"{b} {BANDS[b][0]:g}-{BANDS[b][1]:g}"
                                          for b in band_names))
    head = (f"\n  {'T [K]':>6s} {'rho':>6s} | "
            + " ".join(f"{b:>10s}" for b in band_names))
    print(head)
    print("  " + "-" * (len(head) - 3))
    for r in data:
        print(f"  {r['T']:>6.0f} {r['rho']:>6.0f} | "
              + " ".join(f"{r['ratio'][b]:>10.3f}" for b in band_names))

    print(f"\n  stretch/bend ratio WITHIN each source (independent of any overall "
          f"intensity scale):")
    head = (f"  {'T [K]':>6s} {'rho':>6s} | {'FP':>10s} {args.line_list.upper():>10s} "
            f"| {'FP / LL':>9s}")
    print(head)
    print("  " + "-" * (len(head) - 2))
    for r in data:
        print(f"  {r['T']:>6.0f} {r['rho']:>6.0f} | {r['stretch_over_bend_fp']:>10.3f} "
              f"{r['stretch_over_bend_line']:>10.3f} | "
              f"{r['stretch_over_bend_fp'] / r['stretch_over_bend_line']:>9.3f}")

    core = ("libration", "bend", "stretch")
    lo = min(min(r["ratio"][b] for b in core) for r in data)
    hi = max(max(r["ratio"][b] for b in core) for r in data)
    wlo = min(r["ratio"]["window"] for r in data)
    whi = max(r["ratio"]["window"] for r in data)
    verdict, monotone = _verdict(data)
    print(f"\n  libration/bend/stretch ratios span {lo:.3f}-{hi:.3f}; "
          f"window spans {wlo:.3f}-{whi:.3f}.")
    print(f"  READING: {verdict}")
    if monotone:
        order = ("libration", "bend", "window", "stretch")
        print(f"  The monotonic rise, point by point ({' < '.join(order)}):")
        for r in data:
            print(f"    T={r['T']:>5.0f} rho={r['rho']:>3.0f}: "
                  + "  ".join(f"{r['ratio'][b]:.3f}" for b in order))
    print(f"  The window column is the collision-induced contribution — a line list "
          f"has almost\n  nothing between the bend and stretch fundamentals, so a "
          f"large ratio there is expected\n  and is not evidence about the dipole "
          f"surface either way.")

    out_dir = cfg.output_path
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / "band_diagnostic.json"
    atomicio.write_text(path, json.dumps({
        "meta": (
            "Band-integrated absorption A = int nu (1 - exp(-hc nu / kT)) kappa(nu) dnu "
            "for first-principles vs bare line-list spectra at the six dilute points "
            "(A61). Each source is integrated on its own delivered wavenumber grid over "
            "identical band limits; no radiative weighting enters, so the ratio is a "
            "comparison of integrated band intensity."
        ),
        "bands_cm": {b: list(v) for b, v in BANDS.items()},
        "line_list": args.line_list,
        "core_bands_ratio_range": [lo, hi],
        "window_ratio_range": [wlo, whi],
        "verdict": verdict,
        "ratio_monotone_in_band_centre": monotone,
        "results": data,
    }, indent=1) + "\n")
    print(f"\n  wrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
