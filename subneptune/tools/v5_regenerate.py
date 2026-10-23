#!/usr/bin/env python3
"""Regenerate the first-principles spectra as v5 and carry their errors (A83).

    python tools/v5_regenerate.py config/fiducial.yaml --write

Output: ``outputs/v5_regenerate.json``, spectra in ``<data_root>/fp_opacity_v5``.

The estimator is A81's, settled: a **Hann-tapered periodogram** of the
mean-removed cell dipole with ``sum(w^2)`` normalisation, the harmonic quantum
correction entering **once** (as the cancellation A79 derived), and the
**continuous** ``omega^2`` prefactor. A82 measured the discrete-derivative form
``[2 sin(w dt/2)/dt]^2`` as an alternative and it moves kappa_R by -0.71% to
-1.26%; that form is a finite-difference diagnostic and is not used here.

The smoothing width is **sd = 10 cm^-1 at every cell**, replacing v3's 20 (and
15 at the four 1500 K points). A82 chose it: it is the widest width whose
noise-corrected kappa_R is stationary to 1.0%, and its smoothing bias is
0.7-1.5% at five of the six cells tested. At (1000 K, 10 kg m^-3) the bias is
4.4%, which A82 could not reduce because the statistical uncertainty there is
8.0% -- the 30 ps record, not the width, is the limit.

Since that uncertainty is now the dominant systematic in the dilute corner, each
v5 file carries it: kappa_R and kappa_P are recomputed from the two 15 ps halves
of the record and stored as ``kappa_R_sigma`` and ``kappa_P_sigma``, with
``sigma_30 = |a - b| / 2`` (the halves' difference has variance 2 sigma_15^2 and
sigma_30 = sigma_15 / sqrt(2), so the spread over two is the full-record sigma).
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np

from config import Config, load_config
from opacity.dipole_spectrum import spectrum_from_dipoles

DIPOLE_DIR = Path("ed_data/arc_dipoles")
V5_SUBDIR = "fp_opacity_v5"
SMOOTH_CM = 10.0
UNITS_MARKER = "units_fixed"
UNITS_VALUE = "cm2_per_g:A83"
VERSION_VALUE = "v5:tapered-periodogram:sd10"
ESTIMATOR = ("Hann-tapered periodogram, sum(w^2) normalisation, single harmonic "
             "correction, continuous omega^2; sigma = beta w^2 I2/(6 eps0 c N)")


def _series_dir(root: Path, t_k: float, rho: float) -> Path:
    for suffix in ("_05fs", ""):
        d = root / DIPOLE_DIR / f"T{t_k:g}_rho{rho:g}{suffix}"
        if (d / "dipoles.npy").is_file():
            return d
    raise FileNotFoundError(f"no dipole series for T={t_k:g}, rho={rho:g}")


def regenerate(cfg: Config, root: Path, write: bool) -> dict:
    dest = cfg.data_root_path / V5_SUBDIR
    if write:
        dest.mkdir(parents=True, exist_ok=True)
    rows = []
    for path in sorted((cfg.data_root_path / "fp_opacity").glob("T*_opacity.npz")):
        with np.load(path) as z:
            t_k, rho = float(z["T"]), float(z["rho_kg_m3"])
            n_mol, dt_fs = int(z["N"]), float(z["dt_fs"])
            v3_smooth = float(z["smooth_cm"])
            v3_planck, v3_ross = float(z["planck"]), float(z["rosseland"])
        m = np.load(_series_dir(root, t_k, rho) / "dipoles.npy")
        full = spectrum_from_dipoles(m, t_k, dt_fs, n_mol, rho / 1000.0,
                                     smooth_sigma_cm=SMOOTH_CM)
        half = len(m) // 2
        a = spectrum_from_dipoles(m[:half], t_k, dt_fs, n_mol, rho / 1000.0,
                                  smooth_sigma_cm=SMOOTH_CM)
        b = spectrum_from_dipoles(m[half:], t_k, dt_fs, n_mol, rho / 1000.0,
                                  smooth_sigma_cm=SMOOTH_CM)
        kr_sigma = abs(a["rosseland"] - b["rosseland"]) / 2.0
        kp_sigma = abs(a["planck"] - b["planck"]) / 2.0
        rows.append({
            "T_K": t_k, "rho_kg_m3": rho, "N": n_mol, "dt_fs": dt_fs,
            "n_frames": int(m.shape[0]),
            "kappa_R": full["rosseland"], "kappa_R_sigma": kr_sigma,
            "kappa_R_sigma_frac": kr_sigma / full["rosseland"],
            "kappa_P": full["planck"], "kappa_P_sigma": kp_sigma,
            "kappa_P_sigma_frac": kp_sigma / full["planck"],
            "v3_kappa_R": v3_ross, "v3_kappa_P": v3_planck,
            "v5_over_v3_kappa_R": full["rosseland"] / v3_ross,
            "v5_over_v3_kappa_P": full["planck"] / v3_planck,
            "v3_smooth_cm": v3_smooth})
        if write:
            np.savez(dest / f"T{t_k:g}_rho{rho:g}_opacity.npz",
                     wn=full["wn"], sigma=full["sigma"], sigma_raw=full["sigma_raw"],
                     alpha=full["alpha"], kappa=full["kappa"],
                     planck=full["planck"], rosseland=full["rosseland"],
                     kappa_R_sigma=kr_sigma, kappa_P_sigma=kp_sigma,
                     T=t_k, rho_kg_m3=rho, N=n_mol, dt_fs=dt_fs,
                     smooth_cm=SMOOTH_CM,
                     **{UNITS_MARKER: np.array(UNITS_VALUE)},
                     version=np.array(VERSION_VALUE),
                     estimator=np.array(ESTIMATOR))
    if write:
        # grid.json carries over; the v3 crossover file is v3-derived and is
        # recomputed downstream, so it is not copied under a loadable name.
        src = cfg.data_root_path / "fp_opacity" / "grid.json"
        if src.exists():
            (dest / "grid.json").write_text(src.read_text())
    krf = np.array([r["kappa_R_sigma_frac"] for r in rows])
    kpf = np.array([r["kappa_P_sigma_frac"] for r in rows])
    return {"meta": __doc__.split("\n\n")[0], "smooth_cm": SMOOTH_CM,
            "estimator": ESTIMATOR, "units_marker": {UNITS_MARKER: UNITS_VALUE},
            "version": VERSION_VALUE, "n_cells": len(rows),
            "written": bool(write), "directory": str(dest), "rows": rows,
            "kappa_R_sigma_frac": {"min": float(krf.min()), "max": float(krf.max()),
                                   "median": float(np.median(krf))},
            "kappa_P_sigma_frac": {"min": float(kpf.min()), "max": float(kpf.max()),
                                   "median": float(np.median(kpf))},
            "v5_over_v3_kappa_R": {
                "min": float(min(r["v5_over_v3_kappa_R"] for r in rows)),
                "max": float(max(r["v5_over_v3_kappa_R"] for r in rows)),
                "median": float(np.median([r["v5_over_v3_kappa_R"] for r in rows]))},
            "v5_over_v3_kappa_P": {
                "min": float(min(r["v5_over_v3_kappa_P"] for r in rows)),
                "max": float(max(r["v5_over_v3_kappa_P"] for r in rows)),
                "median": float(np.median([r["v5_over_v3_kappa_P"] for r in rows]))}}


def _print(doc: dict) -> None:
    print(f"\n  v5: {doc['n_cells']} cells, sd = {doc['smooth_cm']:g} cm^-1, "
          f"{'written to ' + doc['directory'] if doc['written'] else 'not written'}")
    k = doc["kappa_R_sigma_frac"]
    p = doc["kappa_P_sigma_frac"]
    print(f"    statistical uncertainty  kappa_R {100*k['min']:.1f}-{100*k['max']:.1f}% "
          f"(median {100*k['median']:.1f}%);  kappa_P {100*p['min']:.2f}-"
          f"{100*p['max']:.2f}% (median {100*p['median']:.2f}%)")
    r = doc["v5_over_v3_kappa_R"]
    q = doc["v5_over_v3_kappa_P"]
    print(f"    v5/v3  kappa_R {r['min']:.3f}-{r['max']:.3f} (median {r['median']:.3f});"
          f"  kappa_P {q['min']:.3f}-{q['max']:.3f} (median {q['median']:.3f})")
    worst = sorted(doc["rows"], key=lambda x: -x["kappa_R_sigma_frac"])[:5]
    print(f"\n    least certain cells:")
    for w in worst:
        print(f"      T={w['T_K']:6g} rho={w['rho_kg_m3']:9g}  "
              f"kappa_R = {w['kappa_R']:10.3f} +- {100*w['kappa_R_sigma_frac']:.1f}%")
    print()


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("config", type=Path)
    ap.add_argument("--write", action="store_true")
    args = ap.parse_args(argv)
    cfg = load_config(args.config)
    doc = regenerate(cfg, Path(__file__).resolve().parents[1], args.write)
    out = cfg.output_path / "v5_regenerate.json"
    out.write_text(json.dumps(doc, indent=1) + "\n")
    _print(doc)
    print(f"  wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
