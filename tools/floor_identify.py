#!/usr/bin/env python3
"""Identify the dilute noise floor before any ARC work (A81).

    python tools/floor_identify.py config/fiducial.yaml

Existing dipole series, this machine only. Output: ``outputs/floor_identify.json``.

A80 held propagation because the A75 floor, recomputed through the fixed
estimator, was only 1.8x below the transparency window at (1000 K, 10). A81 asks
what that floor actually is, from the dipole series themselves.

**(1) It is not discontinuities.** Across all 52 cells and 3.1 million
consecutive-frame steps there is **not one** outlier above 10 sigma, and the
largest |dM| anywhere is 0.20 of the box length against the ~0.4 L a wrapped
atomic charge would give -- a 450-sigma event at the dilute cells, and there are
none. The dipole model is not being evaluated across a periodic image.

**(2) It is white in the dipole's time derivative, not in the dipole.** The
apparent curvature of the high-frequency tail is an artefact of writing
sigma ~ omega^2 S with the *continuous* omega^2: the series is sampled at 0.5 fs
and the discrete transfer function is [2 sin(omega dt/2)/dt]^2. Under that, the
tail is flat to within 0.05 in log-log slope from 6000 cm^-1 to Nyquist. Flat in
the discrete sigma means S ~ omega^-2, and since the physical
sigma = beta omega^2 S / (6 eps0 c N), **the floor is a constant in sigma** --
the same level at 2550 cm^-1 as at 30000.

That is the step signature in the brief's sense, and (1) says there are no large
steps, so it is many small ones: dM stays strongly autocorrelated (lag-1 0.91 to
0.98, i.e. smooth) with only 0.2-4% of its power in the flat tail. The residue
is intrinsic dipole-model evaluation noise, ~4000-25000x larger than float32
quantisation, and the production dynamics are NVE VelocityVerlet
(``run_prod_md_05fs.py`` line 50), so it is not thermostat noise either.

**A75 and A52 both assumed the floor was white in M**, which gives
sigma ~ nu^2 and so a floor 3.4x smaller at 2550 than at the 3900-5500 cm^-1
window where it was estimated. It is flat instead, which is why their margins
were optimistic.

**(3) Subtracting it.** The level is measured per cell as the lowest 2000-cm^-1
binned median of sigma between 6000 cm^-1 and Nyquist -- conservative, since at
hot cells real absorption is still decaying there -- with the binned standard
error as its uncertainty (0.4-4.1%). It is then subtracted as a constant.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
from scipy.ndimage import gaussian_filter1d

import opacity.dipole_spectrum as ds
from config import Config, load_config
from opacity.dipole_spectrum import K_B, two_sided_power_spectrum

DIPOLE_DIR = Path("ed_data/arc_dipoles")
NA, MW = 6.02214076e23, 18.01528
OUTLIER_SIGMA = 10.0
FLOOR_BAND_CM = (6000.0, 32000.0)
FLOOR_BIN_CM = 2000.0
WINDOW_CM = (1900.0, 3100.0)
MARK_CM = 2550.0
MARGIN_REQUIRED = 10.0


def _cells(cfg: Config, root: Path):
    """(T, rho, N, dt_fs, smooth_cm, series path) for every delivered cell."""
    groups: dict[tuple[float, float], list[str]] = {}
    for d in sorted(os.listdir(root / DIPOLE_DIR)):
        m = re.match(r"T(\d+(?:\.\d+)?)_rho([\d.]+?)(_05fs)?$", d)
        if m:
            groups.setdefault((float(m.group(1)), float(m.group(2))), []).append(d)
    out = []
    for path in sorted((cfg.data_root_path / "fp_opacity").glob("T*_opacity.npz")):
        with np.load(path) as z:
            t_k, rho = float(z["T"]), float(z["rho_kg_m3"])
            meta = (int(z["N"]), float(z["dt_fs"]), float(z["smooth_cm"]))
        cands = groups.get((t_k, rho), [])
        # A80's manifest resolved every cell to the _05fs 60002-frame series
        pick = next((c for c in cands if c.endswith("_05fs")), cands[0] if cands else None)
        if pick:
            out.append((t_k, rho, *meta, root / DIPOLE_DIR / pick / "dipoles.npy"))
    return out


def discontinuity_scan(cells) -> dict:
    """(1) Consecutive-frame steps against the smooth expectation."""
    rows, total_out, total_steps = [], 0, 0
    for t_k, rho, n_mol, dt_fs, _, path in cells:
        m = np.load(path).astype(float)
        box_a = (n_mol * MW / NA / (rho * 1e-3)) ** (1.0 / 3.0) * 1e8
        step = np.linalg.norm(np.diff(m, axis=0), axis=1)
        med = float(np.median(step))
        mad = float(np.median(np.abs(step - med)) * 1.4826)
        n_out = int((step > med + OUTLIER_SIGMA * mad).sum())
        rows.append({"T_K": t_k, "rho_kg_m3": rho, "N": n_mol,
                     "box_length_A": box_a, "sigma_step": mad,
                     "n_outliers_10sigma": n_out, "n_steps": int(step.size),
                     "max_step": float(step.max()),
                     "max_step_over_box_length": float(step.max() / box_a),
                     "max_step_in_sigma": float((step.max() - med) / mad)})
        total_out += n_out
        total_steps += int(step.size)
    return {"rows": rows, "total_outliers": total_out, "total_steps": total_steps,
            "outlier_rate_per_1e4": 1e4 * total_out / total_steps,
            "max_step_over_box_length": max(r["max_step_over_box_length"] for r in rows),
            "wrap_step_would_be": ("~0.4-0.8 of the box length, for an atomic "
                                   "charge of 0.4-0.8 e"),
            "verdict": ("no discontinuities" if total_out == 0
                        else f"{total_out} outliers")}


def spectral_shape(cells) -> dict:
    """(2) The tail exponent, with the continuous and the discrete prefactor."""
    rows = []
    for t_k, rho, n_mol, dt_fs, _, path in cells:
        dt_s = dt_fs * 1e-15
        omega, itilde = two_sided_power_spectrum(np.load(path), dt_s)
        wn = omega / (2.0 * np.pi) / ds.C_CM
        cont = omega ** 2 * itilde
        disc = (2.0 * np.sin(omega * dt_s / 2.0) / dt_s) ** 2 * itilde
        rec = {"T_K": t_k, "rho_kg_m3": rho, "slope_by_band": {}}
        for lo, hi in ((6000, 10000), (10000, 16000), (16000, 24000), (24000, 32000)):
            m = (wn >= lo) & (wn < hi)
            rec["slope_by_band"][f"{lo}-{hi}"] = {
                "continuous_m": float(np.polyfit(np.log(wn[m]), np.log(cont[m]), 1)[0]),
                "discrete_m": float(np.polyfit(np.log(wn[m]), np.log(disc[m]), 1)[0])}
        rows.append(rec)
    hi_disc = [r["slope_by_band"]["24000-32000"]["discrete_m"] for r in rows]
    hi_cont = [r["slope_by_band"]["24000-32000"]["continuous_m"] for r in rows]
    return {"rows": rows,
            "near_nyquist_discrete_m": {"median": float(np.median(hi_disc)),
                                        "min": float(min(hi_disc)),
                                        "max": float(max(hi_disc))},
            "near_nyquist_continuous_m": {"median": float(np.median(hi_cont)),
                                          "min": float(min(hi_cont)),
                                          "max": float(max(hi_cont))},
            "interpretation": ("flat under the discrete transfer function => "
                               "S ~ omega^-2 => the physical sigma floor is a "
                               "constant in nu")}


def smoothness(cells) -> dict:
    """Is the residue a few large steps or many small ones? dM's own ACF."""
    rows = []
    for t_k, rho, n_mol, dt_fs, _, path in cells:
        if rho > 60 and not (t_k == 3000.0 and rho > 2000):
            continue
        m = np.load(path).astype(np.float64)
        dm = np.diff(m, axis=0)
        dm -= dm.mean(0)
        a1 = float(np.mean([np.corrcoef(dm[:-1, k], dm[1:, k])[0, 1] for k in range(3)]))
        dt_s = dt_fs * 1e-15
        omega, itilde = two_sided_power_spectrum(m, dt_s)
        wn = omega / (2.0 * np.pi) / ds.C_CM
        sd = (2.0 * np.sin(omega * dt_s / 2.0) / dt_s) ** 2 * itilde
        hi = (wn >= 24000) & (wn <= 32000)
        rows.append({"T_K": t_k, "rho_kg_m3": rho, "acf1_of_dM": a1,
                     "white_fraction_of_dM_power":
                         float(np.median(sd[hi]) / np.mean(sd[wn > 0]))})
    return {"rows": rows,
            "acf1_range": [min(r["acf1_of_dM"] for r in rows),
                           max(r["acf1_of_dM"] for r in rows)],
            "verdict": ("dM is smooth (lag-1 ACF near 1), so the residue is many "
                        "small steps rather than a few large ones")}


def subtract_and_regate(cells) -> dict:
    """(3) Measure the flat floor, subtract it, and rerun gate (iii)."""
    bins = np.arange(FLOOR_BAND_CM[0], FLOOR_BAND_CM[1] + 1, FLOOR_BIN_CM)
    rows = []
    for t_k, rho, n_mol, dt_fs, smooth, path in cells:
        dt_s = dt_fs * 1e-15
        omega, itilde = two_sided_power_spectrum(np.load(path), dt_s)
        wn = omega / (2.0 * np.pi) / ds.C_CM
        sig = ((1.0 / (K_B * t_k)) * omega ** 2 * itilde
               / (6.0 * ds.EPS0 * ds.C_VAC * n_mol) * 1e4)
        meds, errs = [], []
        for lo, hi in zip(bins[:-1], bins[1:]):
            m = (wn >= lo) & (wn < hi)
            meds.append(float(np.median(sig[m])))
            errs.append(float(np.std(sig[m]) / np.sqrt(m.sum()) * 1.2533))
        i = int(np.argmin(meds))
        floor, ferr = meds[i], errs[i]
        keep = (wn > 0.0) & (wn <= 6000.0)
        grid = np.arange(max(wn[keep].min(), 0.0), 6001.0, 1.0)
        raw = np.nan_to_num(np.interp(grid, wn[keep], sig[keep],
                                      left=np.nan, right=np.nan), nan=0.0)
        before = gaussian_filter1d(raw, smooth)
        after = gaussian_filter1d(np.clip(raw - floor, 0.0, None), smooth)
        _, k0 = ds.sigma_to_alpha_kappa(before, rho / 1000.0)
        _, k1 = ds.sigma_to_alpha_kappa(after, rho / 1000.0)
        _, kf = ds.sigma_to_alpha_kappa(np.full_like(grid, floor), rho / 1000.0)
        p0, r0 = ds.planck_and_rosseland(grid, k0, t_k)
        p1, r1 = ds.planck_and_rosseland(grid, np.maximum(k1, 1e-30), t_k)
        sel = (grid >= WINDOW_CM[0]) & (grid <= WINDOW_CM[1])
        rows.append({
            "T_K": t_k, "rho_kg_m3": rho,
            "floor_sigma_cm2": floor, "floor_rel_uncertainty": ferr / floor,
            "floor_over_window_sigma": float(floor / np.median(raw[sel])),
            "margin_before": float(np.min(k0[sel]) / float(np.median(kf[sel]))),
            "margin_after": float(np.min(k1[sel]) / float(np.median(kf[sel]))),
            "kappa_R_before": r0, "kappa_R_after": r1, "kappa_R_ratio": r1 / r0,
            "kappa_P_ratio": p1 / p0})
    dilute = [r for r in rows if r["rho_kg_m3"] <= 60
              or (r["T_K"] == 800.0 and r["rho_kg_m3"] == 100.0)]
    by_band = {}
    for lo, hi, lab in ((0, 60, "<=60"), (60, 300, "60-300"), (300, 1e9, ">300")):
        s = [r for r in rows if lo < r["rho_kg_m3"] <= hi]
        if s:
            by_band[lab] = {
                "n": len(s),
                "median_floor_over_window": float(np.median(
                    [r["floor_over_window_sigma"] for r in s])),
                "median_kappa_R_ratio": float(np.median([r["kappa_R_ratio"] for r in s])),
                "median_kappa_P_ratio": float(np.median([r["kappa_P_ratio"] for r in s]))}
    return {"rows": rows, "by_density_band": by_band,
            "worst_margin_before": float(min(r["margin_before"] for r in dilute)),
            "worst_margin_after": float(min(r["margin_after"] for r in dilute)),
            "floor_uncertainty_range": [float(min(r["floor_rel_uncertainty"] for r in rows)),
                                        float(max(r["floor_rel_uncertainty"] for r in rows))],
            "n_cells_losing_over_20pc_of_kappa_R":
                int(sum(1 for r in rows if r["kappa_R_ratio"] < 0.8)),
            "gate_iii_passes": bool(min(r["margin_after"] for r in dilute)
                                    >= MARGIN_REQUIRED)}


def sha256_manifest(root: Path) -> dict:
    """A tracked fingerprint of the bundle, which is itself gitignored (A80)."""
    out = {}
    for d in sorted(os.listdir(root / DIPOLE_DIR)):
        f = root / DIPOLE_DIR / d / "dipoles.npy"
        if not f.is_file():
            continue
        h = hashlib.sha256()
        with open(f, "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                h.update(chunk)
        out[d] = {"sha256": h.hexdigest(), "bytes": f.stat().st_size}
    return out


def run(cfg: Config) -> dict:
    root = Path(__file__).resolve().parents[1]
    cells = _cells(cfg, root)
    return {"meta": __doc__.split("\n\n")[0], "n_cells": len(cells),
            "discontinuities": discontinuity_scan(cells),
            "spectral_shape": spectral_shape(cells),
            "smoothness": smoothness(cells),
            "subtraction": subtract_and_regate(cells)}


def _print(doc: dict) -> None:
    d = doc["discontinuities"]
    print(f"\n  (1) Discontinuities: {d['total_outliers']} outliers above "
          f"{OUTLIER_SIGMA:g} sigma in {d['total_steps']:,} steps across "
          f"{doc['n_cells']} cells")
    print(f"      largest step anywhere = {d['max_step_over_box_length']:.3f} of the "
          f"box length; a wrapped charge would give {d['wrap_step_would_be']}")
    s = doc["spectral_shape"]
    print(f"\n  (2) Tail exponent near Nyquist: continuous prefactor "
          f"{s['near_nyquist_continuous_m']['median']:+.2f}, discrete "
          f"{s['near_nyquist_discrete_m']['median']:+.2f} "
          f"[{s['near_nyquist_discrete_m']['min']:+.2f}, "
          f"{s['near_nyquist_discrete_m']['max']:+.2f}]")
    print(f"      {s['interpretation']}")
    sm = doc["smoothness"]
    print(f"      dM lag-1 ACF {sm['acf1_range'][0]:.2f}-{sm['acf1_range'][1]:.2f}: "
          f"{sm['verdict']}")
    t = doc["subtraction"]
    print(f"\n  (3) Floor subtracted (uncertainty "
          f"{100*t['floor_uncertainty_range'][0]:.1f}-"
          f"{100*t['floor_uncertainty_range'][1]:.1f}%)")
    print(f"      {'rho band':>10s}{'n':>4s}{'floor/window':>14s}"
          f"{'kappa_R after/before':>22s}{'kappa_P':>10s}")
    for lab, v in t["by_density_band"].items():
        print(f"      {lab:>10s}{v['n']:4d}{v['median_floor_over_window']:14.2f}"
              f"{v['median_kappa_R_ratio']:22.3f}{v['median_kappa_P_ratio']:10.3f}")
    print(f"      gate (iii) margin: {t['worst_margin_before']:.1f}x before, "
          f"{t['worst_margin_after']:.1f}x after  ->  "
          f"{'PASS' if t['gate_iii_passes'] else 'FAIL'}")
    print(f"      cells losing >20% of kappa_R: "
          f"{t['n_cells_losing_over_20pc_of_kappa_R']} of {doc['n_cells']}\n")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("config", type=Path)
    ap.add_argument("--manifest", type=Path, default=None,
                    help="write a SHA-256 manifest of the dipole bundle here")
    args = ap.parse_args(argv)
    cfg = load_config(args.config)
    root = Path(__file__).resolve().parents[1]
    doc = run(cfg)
    path = cfg.output_path / "floor_identify.json"
    path.write_text(json.dumps(doc, indent=1) + "\n")
    _print(doc)
    if args.manifest:
        args.manifest.write_text(json.dumps(sha256_manifest(root), indent=1) + "\n")
        print(f"  wrote {args.manifest}")
    print(f"  wrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
