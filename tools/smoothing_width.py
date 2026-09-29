#!/usr/bin/env python3
"""Choose the smoothing width for the first-principles spectra (A82).

    python tools/smoothing_width.py config/fiducial.yaml

Existing dipole series. Output: ``outputs/smoothing_width.json``.

A81 removed the estimator's spectral leakage, which had been filling the
transparency windows. With the windows genuinely deep, kappa_R -- a harmonic
mean, so set by the minima -- became sensitive to the Gaussian smoothing width,
and A81 found it did not converge: 13.2% between sd 10 and 40 and still moving
2.6% between 5 and 10. Which width is right is a physics question, and this tool
answers it by separating the two things that move kappa_R in opposite directions.

  * **Smoothing bias.** Convolving a spectrum whose windows are two or three
    decades deep fills them from the band edges, which raises kappa_R. It grows
    with the width and is a property of the true spectrum, not of the noise. It
    is measured on the Welch-8 estimate, which has ~3x the degrees of freedom of
    the full-record periodogram and so is nearly noise-free, relative to its own
    sd = 2 value.
  * **Noise bias.** A chi^2_nu spectral estimate has `E[1/S_hat] = (1/S)
    nu/(nu-2)`, and kappa_R integrates `1/kappa`, so noise biases kappa_R *low*
    by `(nu-2)/nu`. It shrinks with the width, because smoothing raises nu. The
    correction applied here is `nu/(nu-2)` with

        nu = 6 * 2 s sqrt(pi) / (1.5 * dnu_fft)

    -- 6 for three Cartesian components of a chi^2_2 periodogram, the Gaussian
    kernel's effective bin count for width `s`, and 1.5 for the Hann taper's
    noise-equivalent bandwidth, which correlates adjacent bins.

The statistical uncertainty is taken empirically from the two 15 ps halves of
each record rather than from the nu above, so the analytic nu is used only for
the bias it predicts, never for an error bar.

The estimator is A81's throughout: Hann-tapered periodogram, ``sum(w^2)``
normalisation, one harmonic correction, and the **discrete-derivative**
prefactor `[2 sin(w dt/2)/dt]^2` in place of `w^2` as the brief specifies -- a
choice worth 1.2% at 4000 cm^-1 and less below, quantified per cell here.
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
from scipy.ndimage import gaussian_filter1d

import opacity.dipole_spectrum as ds
from config import Config, load_config

DIPOLE_DIR = Path("ed_data/arc_dipoles")
GRID = np.arange(1.0, 6001.0, 1.0)
WIDTHS = (2.0, 5.0, 7.5, 10.0, 15.0, 20.0, 30.0, 40.0)
REPORT_WIDTHS = (2.0, 5.0, 7.5, 10.0, 15.0, 20.0, 30.0, 40.0)
# Six cells spanning the grid; (1000, 10) and (700, 60) are required by the brief.
CELLS = ((1000.0, 10.0), (700.0, 60.0), (1500.0, 300.0),
         (2000.0, 1200.0), (2500.0, 500.0), (3000.0, 2228.71))
STATIONARY_TOL = 0.02
BIAS_TOL = 0.02
HANN_NEB = 1.5          # Hann taper noise-equivalent bandwidth, in bins
WELCH_SEGMENTS = 8


def _series(root: Path, t_k: float, rho: float):
    for suffix in ("_05fs", ""):
        d = root / DIPOLE_DIR / f"T{t_k:g}_rho{rho:g}{suffix}"
        if (d / "dipoles.npy").is_file():
            return np.load(d / "dipoles.npy")
    raise FileNotFoundError((t_k, rho))


def sigma_on_grid(m_eang, t_k, dt_s, n_molecules, n_segments=1, discrete=True):
    """A81's estimator, optionally Welch-averaged over overlapping segments."""
    m = np.asarray(m_eang, float) * ds.EANG_TO_DEBYE * ds.DEBYE_TO_CM
    m = m - m.mean(0)
    n = len(m)
    if n_segments == 1:
        spans = [(0, n)]
    else:
        seg = int(2 * n / (n_segments + 1))
        spans = [(s, s + seg) for s in range(0, n - seg + 1, seg // 2)]
    acc, count = None, 0
    for a, b in spans:
        s = m[a:b] - m[a:b].mean(0)
        length = len(s)
        w = np.hanning(length)
        ft = np.fft.rfft(s * w[:, None], axis=0)
        p = dt_s * np.sum(np.abs(ft) ** 2, axis=1) / np.sum(w ** 2)
        acc = p if acc is None else acc + p
        count += 1
    itilde = acc / count
    length = spans[0][1] - spans[0][0]
    omega = 2.0 * np.pi * np.fft.rfftfreq(length, d=dt_s)
    wn = omega / (2.0 * np.pi) / ds.C_CM
    pre = (2.0 * np.sin(omega * dt_s / 2.0) / dt_s) ** 2 if discrete else omega ** 2
    sig = ((1.0 / (ds.K_B * t_k)) * pre * itilde
           / (6.0 * ds.EPS0 * ds.C_VAC * n_molecules) * 1e4)
    keep = (wn > 0.0) & (wn <= 6000.0)
    resolution = float(wn[1] - wn[0])
    return (np.nan_to_num(np.interp(GRID, wn[keep], sig[keep],
                                    left=np.nan, right=np.nan), nan=0.0),
            resolution, len(spans))


def kappa_r(raw, sd, t_k, rho):
    _, kap = ds.sigma_to_alpha_kappa(gaussian_filter1d(raw, sd), rho / 1000.0)
    return ds.planck_and_rosseland(GRID, kap, t_k)[1]


def dof(sd_cm, resolution_cm, n_averaged=1):
    """Effective chi^2 degrees of freedom of the smoothed tapered periodogram."""
    n_eff = 2.0 * sd_cm * np.sqrt(np.pi) / (HANN_NEB * resolution_cm)
    return 6.0 * n_eff * n_averaged


def study(root: Path, cfg: Config) -> dict:
    out = {"widths_cm": list(REPORT_WIDTHS), "cells": []}
    for t_k, rho in CELLS:
        path = cfg.data_root_path / "fp_opacity" / f"T{t_k:g}_rho{rho:g}_opacity.npz"
        with np.load(path) as z:
            n_mol, dt_fs = int(z["N"]), float(z["dt_fs"])
        dt_s = dt_fs * 1e-15
        m = _series(root, t_k, rho)
        full, res, _ = sigma_on_grid(m, t_k, dt_s, n_mol)
        cont, _, _ = sigma_on_grid(m, t_k, dt_s, n_mol, discrete=False)
        welch, res_w, n_seg = sigma_on_grid(m, t_k, dt_s, n_mol,
                                            n_segments=WELCH_SEGMENTS)
        half = len(m) // 2
        h1, _, _ = sigma_on_grid(m[:half], t_k, dt_s, n_mol)
        h2, _, _ = sigma_on_grid(m[half:], t_k, dt_s, n_mol)
        # Two bias data. Welch-8 is nearly noise-free but has 5 cm^-1 intrinsic
        # resolution, so it cannot certify a 2% bias for widths near 5 -- its
        # sd = 2 value is already smoothed to an effective ~5.4 cm^-1. The
        # full-record noise-corrected value at sd = 2 is properly resolved
        # (1.11 cm^-1) at the cost of being noisy, and is the datum used for the
        # verdict; the Welch one is kept as a cross-check.
        ref_w = kappa_r(welch, WIDTHS[0], t_k, rho)
        rec = {"T_K": t_k, "rho_kg_m3": rho, "N": n_mol,
               "fft_resolution_cm": res, "welch_resolution_cm": res_w,
               "welch_n_segments": n_seg, "rows": []}
        for sd in REPORT_WIDTHS:
            nu = dof(sd, res)
            k_raw = kappa_r(full, sd, t_k, rho)
            k_cor = k_raw * nu / (nu - 2.0)
            k_w = kappa_r(welch, sd, t_k, rho)
            a, b = kappa_r(h1, sd, t_k, rho), kappa_r(h2, sd, t_k, rho)
            rec["rows"].append({
                "sd_cm": sd, "nu": nu,
                "kappa_R_as_delivered": k_raw,
                "noise_correction": nu / (nu - 2.0),
                "kappa_R_noise_corrected": k_cor,
                "kappa_R_welch": k_w,
                "smoothing_bias": k_w / ref_w - 1.0,
                "half_records": [a, b],
                "half_spread_frac": abs(a - b) / (a + b) * 2.0,
                "discrete_vs_continuous_prefactor":
                    kappa_r(full, sd, t_k, rho) / kappa_r(cont, sd, t_k, rho) - 1.0,
            })
        nu2 = dof(WIDTHS[0], res)
        ref_full = kappa_r(full, WIDTHS[0], t_k, rho) * nu2 / (nu2 - 2.0)
        rec["bias_datum_full_record_sd2"] = ref_full
        for r in rec["rows"]:
            r["smoothing_bias_vs_full_record"] = (
                r["kappa_R_noise_corrected"] / ref_full - 1.0)
            # sigma on the full record from two 15 ps halves: their difference
            # has variance 2 sigma_15^2 and sigma_30 = sigma_15/sqrt(2), so
            # sigma_30/kappa = spread/2.
            r["statistical_uncertainty_frac"] = r["half_spread_frac"] / 2.0
        for i, r in enumerate(rec["rows"]):
            r["step_from_previous"] = (
                None if i == 0
                else r["kappa_R_noise_corrected"]
                / rec["rows"][i - 1]["kappa_R_noise_corrected"] - 1.0)
        out["cells"].append(rec)
    return out


def choose(doc: dict) -> dict:
    """The width where the noise-corrected kappa_R is stationary and the
    smoothing bias is small, at every cell."""
    per_width = {}
    for sd in REPORT_WIDTHS:
        steps, biases, stat, failing = [], [], [], []
        for c in doc["cells"]:
            r = next(x for x in c["rows"] if x["sd_cm"] == sd)
            if r["step_from_previous"] is not None:
                steps.append(abs(r["step_from_previous"]))
            b = abs(r["smoothing_bias_vs_full_record"])
            biases.append(b)
            stat.append(r["statistical_uncertainty_frac"])
            if b >= BIAS_TOL or (r["step_from_previous"] is not None
                                 and abs(r["step_from_previous"]) >= STATIONARY_TOL):
                failing.append([c["T_K"], c["rho_kg_m3"]])
        per_width[sd] = {
            "worst_step_from_previous": max(steps) if steps else None,
            "worst_smoothing_bias": max(biases),
            "worst_statistical_uncertainty": max(stat),
            "n_cells_passing": len(doc["cells"]) - len(failing),
            "cells_failing": failing,
            "stationary": bool(steps and max(steps) < STATIONARY_TOL),
            "bias_ok": bool(max(biases) < BIAS_TOL)}
    # sd = 2 is the bias datum, so its own bias is 0 by construction and it
    # cannot be a candidate; the choice is made among the widths above it.
    ok = [sd for sd, v in per_width.items()
          if sd > 2.0 and v["stationary"] and v["bias_ok"]]
    return {"per_width": per_width, "satisfying_both": ok,
            "chosen": (max(ok) if ok else None),
            "stationary_tol": STATIONARY_TOL, "bias_tol": BIAS_TOL}


def run(cfg: Config) -> dict:
    root = Path(__file__).resolve().parents[1]
    doc = study(root, cfg)
    doc["choice"] = choose(doc)
    return doc


def _print(doc: dict) -> None:
    for c in doc["cells"]:
        print(f"\n  T = {c['T_K']:g} K, rho = {c['rho_kg_m3']:g} "
              f"(resolution {c['fft_resolution_cm']:.2f} cm^-1)")
        print(f"    {'sd':>5s}{'nu':>8s}{'kR':>9s}{'kR corr':>9s}{'step':>8s}"
              f"{'kR Welch':>10s}{'sm. bias':>10s}{'half spread':>13s}")
        for r in c["rows"]:
            st = "" if r["step_from_previous"] is None else f"{100*r['step_from_previous']:+7.1f}"
            print(f"    {r['sd_cm']:5g}{r['nu']:8.0f}{r['kappa_R_as_delivered']:9.3f}"
                  f"{r['kappa_R_noise_corrected']:9.3f}{st:>8s}"
                  f"{r['kappa_R_welch']:10.3f}{100*r['smoothing_bias']:9.1f}%"
                  f"{100*r['half_spread_frac']:12.1f}%")
    ch = doc["choice"]
    print(f"\n  Width selection (stationary < {100*ch['stationary_tol']:.0f}%, "
          f"bias < {100*ch['bias_tol']:.0f}%), worst over the six cells")
    print(f"    {'sd':>5s}{'worst step':>12s}{'worst bias':>12s}{'worst stat':>12s}"
          f"{'cells ok':>10s}{'both':>7s}")
    for sd, v in ch["per_width"].items():
        st = "-" if v["worst_step_from_previous"] is None else f"{100*v['worst_step_from_previous']:.1f}%"
        print(f"    {sd:5g}{st:>12s}{100*v['worst_smoothing_bias']:11.1f}%"
              f"{100*v['worst_statistical_uncertainty']:11.1f}%"
              f"{v['n_cells_passing']:6d}/{len(doc['cells'])}"
              f"{'yes' if v['stationary'] and v['bias_ok'] else 'no':>7s}")
    print(f"\n  satisfying both: {ch['satisfying_both'] or 'NONE'}   "
          f"chosen: {ch['chosen']}\n")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("config", type=Path)
    args = ap.parse_args(argv)
    cfg = load_config(args.config)
    doc = run(cfg)
    path = cfg.output_path / "smoothing_width.json"
    path.write_text(json.dumps(doc, indent=1) + "\n")
    _print(doc)
    print(f"  wrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
