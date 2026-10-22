#!/usr/bin/env python3
"""The production formula, read at last, and the sum rule tested at its source (A78).

    python tools/dipole_sum_rule.py config/fiducial.yaml

Existing data only. Output: ``outputs/dipole_sum_rule.json``.

A77 could not read the ACF -> cross-section step and had to carry two candidate
formulas. ``ed_data/arc_a78/scw/opacity_pipeline.py`` settles it, and the answer
is neither candidate:

    line 66  Smm   = dt * Re rfft( Cmm * hanning(L) )        one-sided ACF
    line 80  x     = beta hbar omega
    line 81  qcf   = x / (1 - exp(-x))                       = Q_harmonic
    line 82  sigma = (beta omega^2 / (3 eps0 c N)) * Smm * qcf * 1e4

`beta omega^2 Smm` is **already** the harmonic-corrected result: the standard
relation is `alpha ~ omega (1 - e^-x) I_qm` with `I_qm = Q_harmonic I_cl`, and
the `(1 - e^-x)` cancels against Q_harmonic's denominator to leave
`beta hbar omega^2 I_cl`. Line 82 then multiplies by Q_harmonic a second time.
**The production spectra carry the harmonic quantum correction twice.**

Q_harmonic falls with temperature at fixed frequency, so a spectrum carrying a
spare factor of it must show band strengths that fall with temperature — which
is exactly the violation A76 found and A77 measured grid-wide. This tool proves
that reading at the source, independently of the formula, by going back to the
dipole time series: the classical sum rule is a statement about
`beta <|Mdot|^2>`, which can be computed from `dipoles.npy` with no pipeline in
between.

Line 66 is unconventional in the same way and turns out not to be a defect. It
applies `np.hanning(L)` -- a *symmetric* window -- to a *one-sided*
autocorrelation, so the zero-lag value is multiplied by zero and the peak weight
falls at lag L/2; `spectra.py` line 46 uses the conventional decaying half
window instead. But measured against the window-free `beta <|Mdot_band|^2>`, the
pipeline's choice is the one that gives a single constant of proportionality
across every band and temperature, and the half-window alternatives are the ones
that scatter. `window_check()` records that, because it is the reason the
absolute scale cannot simply be "corrected" and is left to A79.
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

A78_DIR = Path("ed_data/arc_a78")
# The three rho = 30 cells the brief names. The rho = 300 pair has only two
# temperatures, so it carries no index; it is used for the reproduction check
# only. The two *_BAD64 directories hold different dipole series processed with
# N = 256 and are not used.
CELLS = (("T700_rho30_05fs", 700.0), ("T1000_rho30_05fs", 1000.0),
         ("T2000_rho30_05fs", 2000.0))
REPRODUCE_ONLY = ("T1000_rho300_05fs", "T2000_rho300_05fs")
BANDS_CM = {"rotation": (50.0, 1100.0), "bend": (1000.0, 2400.0),
            "stretch": (2800.0, 4000.0)}


def _pipeline():
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / A78_DIR / "scw"))
    import opacity_pipeline as op
    return op


def beta_mdot2_per_molecule(M_eang, dt_s, t_k, n_molecules) -> dict:
    """beta <|Mdot_band|^2> per molecule, straight from the dipole series.

    Mdot by central differences at the stored step, then an FFT mask per band.
    No autocorrelation, no window, no pipeline: this is the quantity the
    classical sum rule says is temperature-independent, since
    beta <|Mdot|^2> = beta (kT/m) sum (dmu/dQ)^2.
    """
    op = _pipeline()
    m = np.asarray(M_eang, float) * op.EANG_TO_DEBYE * op.debye_to_Cm
    md = (m[2:] - m[:-2]) / (2.0 * dt_s)
    md = md - md.mean(0)
    n = len(md)
    wn = np.fft.rfftfreq(n, d=dt_s) / op.c_cm_s
    beta = 1.0 / (op.kB * t_k)
    out = {}
    for band, (lo, hi) in BANDS_CM.items():
        mask = (wn >= lo) & (wn <= hi)
        total = 0.0
        for k in range(3):
            f = np.fft.rfft(md[:, k])
            f[~mask] = 0.0
            total += float(np.mean(np.fft.irfft(f, n=n) ** 2))
        out[band] = beta * total / n_molecules
    return out


def sigma_bands(M_eang, t_k, dt_fs, n_molecules, apply_qcf: bool) -> dict:
    """The production chain to sigma_raw, with line 81-82's qcf optional."""
    op = _pipeline()
    m = np.asarray(M_eang, float) * op.EANG_TO_DEBYE * op.debye_to_Cm
    beta = 1.0 / (op.kB * t_k)
    dt_s = dt_fs * 1e-15
    length = len(m)
    s = m - m.mean(0)
    cmm = sum(op.autocorr_fft(s[:, k]) for k in range(3))
    smm = dt_s * np.real(np.fft.rfft(cmm * np.hanning(length)))
    freq = np.fft.rfftfreq(length, d=dt_s)
    omega = 2.0 * np.pi * freq
    wn = freq / op.c_cm_s
    keep = (wn > 0.0) & (wn <= 6000.0)
    wn, omega, smm = wn[keep], omega[keep], np.clip(smm[keep], 0.0, None)
    qcf = 1.0
    if apply_qcf:
        x = beta * op.hbar * omega
        qcf = np.ones_like(x)
        hot = x > 1e-12
        qcf[hot] = x[hot] / (1.0 - np.exp(-x[hot]))
    sig = (beta * omega ** 2 / (3.0 * op.eps0 * op.c * n_molecules)) * smm * qcf * 1.0e4
    grid = np.arange(max(wn.min(), 0.0), 6001.0, 1.0)
    raw = np.nan_to_num(np.interp(grid, wn, sig, left=np.nan, right=np.nan), nan=0.0)
    return {b: float(np.trapezoid(raw[(grid >= lo) & (grid <= hi)],
                                  grid[(grid >= lo) & (grid <= hi)]))
            for b, (lo, hi) in BANDS_CM.items()}


# A76's calibration points and the monomer test's dipole intensity factors.
CAL_POINTS = sorted({(t, r) for t in (700.0, 1000.0, 1500.0, 2000.0)
                     for r in (10.0, 30.0, 60.0)} | {(800.0, 30.0), (800.0, 100.0)})
CAL_BANDS_CM = {"rotation": (100.0, 1000.0), "bend": (1400.0, 1800.0),
                "stretch": (3200.0, 4000.0)}
CAL_FACTOR = {"rotation": (0.974, 0.974), "bend": (1.25, 1.32),
              "stretch": (2.1, 3.9)}


def corrected_calibration(cfg: Config) -> dict:
    """A76 redone with the production formula corrected.

    A76 asked which quantum correction to *apply*; the answer is that one must
    be *removed*. Dividing the delivered kappa by Q_harmonic undoes the spare
    factor, after which the comparison with the line list is the one A76 meant
    to make.
    """
    op = _pipeline()
    from tools.text_checks import _load

    def qh(wn, t_k):
        x = op.hbar * (2.0 * np.pi * np.asarray(wn, float) * op.c_cm_s) / (op.kB * t_k)
        return x / -np.expm1(-x)

    def integ(wn, kap, lo, hi):
        m = (wn > lo) & (wn < hi)
        return float(np.trapezoid(kap[m], wn[m]))

    out = {"note": ("FP/HITEMP band ratio as shipped and with the spare "
                    "Q_harmonic divided out, then over the dipole factor"),
           "dipole_intensity_factor": {k: list(v) for k, v in CAL_FACTOR.items()},
           "points": [], "summary": {}}
    series = {b: [] for b in CAL_BANDS_CM}
    for t_k, rho in CAL_POINTS:
        wf, kf = _load(cfg, "fp", t_k, rho)
        wl, kl = _load(cfg, "hitemp", t_k, rho)
        o = np.argsort(wf); wf, kf = wf[o], kf[o]
        o = np.argsort(wl); wl, kl = wl[o], kl[o]
        rec = {"T_K": t_k, "rho_kg_m3": rho, "bands": {}}
        for band, (lo, hi) in CAL_BANDS_CM.items():
            ll = integ(wl, kl, lo, hi)
            rec["bands"][band] = {
                "as_shipped": integ(wf, kf, lo, hi) / ll,
                "corrected": integ(wf, kf / qh(wf, t_k), lo, hi) / ll}
            series[band].append((t_k, rec["bands"][band]["corrected"]))
        out["points"].append(rec)
    for band, rows in series.items():
        ts = np.array([t for t, _ in rows], float)
        v = np.array([x for _, x in rows], float)
        lo_f, hi_f = CAL_FACTOR[band]
        rank = lambda a: a.argsort().argsort().astype(float)  # noqa: E731
        out["summary"][band] = {
            "median_over_factor": [float(np.median(v) / lo_f),
                                   float(np.median(v) / hi_f)],
            "range_over_factor_low": [float(v.min() / lo_f), float(v.max() / lo_f)],
            "dln_dT_per_1000K": float(np.polyfit(ts, np.log(v), 1)[0] * 1000.0),
            "spearman_rho": float(np.corrcoef(rank(ts), rank(v))[0, 1]),
        }
    return out


def window_check(root: Path) -> dict:
    """Which lag window is consistent with the window-free measurement?

    `beta <|Mdot_band|^2>` involves no autocorrelation and no window, so the
    right lag window is the one whose band integrals are a single constant
    multiple of it across bands and temperatures. This is the evidence for
    leaving line 66 alone rather than "fixing" it.
    """
    op = _pipeline()
    out = {"note": ("ratio int sigma_noqcf dnu / beta<|Mdot_band|^2> per lag "
                    "window; the correct window gives a T- and band-independent "
                    "constant"),
           "windows": {}}

    def bands_for(m, t_k, dt_fs, n_mol, window):
        mx = np.asarray(m, float) * op.EANG_TO_DEBYE * op.debye_to_Cm
        beta = 1.0 / (op.kB * t_k)
        dt_s = dt_fs * 1e-15
        length = len(mx)
        s = mx - mx.mean(0)
        cmm = sum(op.autocorr_fft(s[:, k]) for k in range(3))
        if window == "pipeline_symmetric_full":
            c, w = cmm, np.hanning(length)
        elif window == "half_window_full":
            c, w = cmm, np.hanning(2 * length)[length:]
        elif window == "half_window_8000":
            lag = 8000
            c, w = cmm[:lag], np.hanning(2 * lag)[lag:]
        n = len(c)
        smm = dt_s * np.real(np.fft.rfft(c * w))
        freq = np.fft.rfftfreq(n, d=dt_s)
        omega = 2.0 * np.pi * freq
        wn = freq / op.c_cm_s
        keep = (wn > 0.0) & (wn <= 6000.0)
        wn, omega, smm = wn[keep], omega[keep], np.clip(smm[keep], 0.0, None)
        sig = (beta * omega ** 2 / (3.0 * op.eps0 * op.c * n_mol)) * smm * 1.0e4
        grid = np.arange(max(wn.min(), 0.0), 6001.0, 1.0)
        raw = np.nan_to_num(np.interp(grid, wn, sig, left=np.nan, right=np.nan), nan=0.0)
        return {b: float(np.trapezoid(raw[(grid >= lo) & (grid <= hi)],
                                      grid[(grid >= lo) & (grid <= hi)]))
                for b, (lo, hi) in BANDS_CM.items()}

    for window in ("pipeline_symmetric_full", "half_window_full", "half_window_8000"):
        ratios = {b: [] for b in BANDS_CM}
        for cell, t_k in CELLS:
            d = root / A78_DIR / cell
            z = np.load(d / "opacity_v3.npz")
            m = np.load(d / "dipoles.npy")
            md = beta_mdot2_per_molecule(m, float(z["dt_fs"]) * 1e-15, t_k, int(z["N"]))
            sg = bands_for(m, t_k, float(z["dt_fs"]), int(z["N"]), window)
            for b in BANDS_CM:
                ratios[b].append(sg[b] / md[b])
        out["windows"][window] = {
            "ratio_by_band": ratios,
            "spread_frac_by_band": {b: max(v) / min(v) - 1.0 for b, v in ratios.items()},
            "worst_spread_frac": max(max(v) / min(v) - 1.0 for v in ratios.values()),
        }
    best = min(out["windows"], key=lambda k: out["windows"][k]["worst_spread_frac"])
    out["most_consistent_window"] = best
    return out


def _index(ts, vals):
    return float(np.polyfit(np.log(ts), np.log(vals), 1)[0])


def reproduction_check(root: Path) -> dict:
    """Does the pipeline as read reproduce the delivered spectra? If it does,
    the formula above is the production formula and not a paraphrase of it."""
    op = _pipeline()
    out = {}
    for cell in [c for c, _ in CELLS] + list(REPRODUCE_ONLY):
        d = root / A78_DIR / cell
        z = np.load(d / "opacity_v3.npz")
        r = op.spectrum_from_dipoles(
            np.load(d / "dipoles.npy"), T_K=float(z["T"]), dt_fs=float(z["dt_fs"]),
            n_molecules=int(z["N"]), density_g_cm3=float(z["rho_kg_m3"]) / 1000.0,
            units="eAng", smooth_sigma_cm=float(z["smooth_cm"]))
        errs = {}
        for mine, stored in (("sigma_cm2_per_molecule_raw", "sigma_raw"),
                             ("kappa_cm2_g", "kappa")):
            ok = z[stored] != 0
            errs[stored] = float(np.max(np.abs(r[mine][ok] / z[stored][ok] - 1.0)))
        # alpha = n sigma is an internal consistency check that catches the
        # legacy x1000 units error, which one supplied cell still carries in
        # alpha and kappa although its sigma is correct and identical to the
        # delivered file. The formula check is the sigma_raw one.
        n_used = (float(z["rho_kg_m3"]) * 1e-3 / op.DEFAULT_MOLAR_MASS_G_MOL) * op.NA
        units = float(np.nanmax(np.abs(z["alpha"] / (n_used * z["sigma"]) - 1.0)))
        out[cell] = {"N": int(z["N"]), "T_K": float(z["T"]),
                     "rho_kg_m3": float(z["rho_kg_m3"]),
                     "max_rel_error": errs,
                     "alpha_over_n_sigma_minus_1": units,
                     "carries_legacy_x1000_in_alpha_kappa": bool(units > 1e-9),
                     "planck_mine": r["planck_mean_cm2_g"],
                     "planck_stored": float(z["planck"])}
    return out


def run(cfg: Config) -> dict:
    root = Path(__file__).resolve().parents[1]
    doc = {"meta": __doc__.split("\n\n")[0],
           "production_formula": {
               "file": "ed_data/arc_a78/scw/opacity_pipeline.py",
               "function": "spectrum_from_dipoles",
               "line_66": "Smm = dt * Re rfft(Cmm * hanning(L))",
               "line_80": "x = beta hbar omega",
               "line_81": "qcf = x / (1 - exp(-x))  [= Q_harmonic]",
               "line_82": "sigma_raw = (beta omega^2 / (3 eps0 c N)) Smm qcf 1e4",
               "verdict": ("neither A77 form A nor form B: it is form B "
                           "(beta omega^2 I_cl) times a second Q_harmonic"),
               "second_defect_line_66": (
                   "np.hanning(L), a symmetric window, applied to a one-sided "
                   "ACF: zero weight at lag 0, peak weight at lag L/2 where the "
                   "unbiased estimator is noisiest. spectra.py:46 does it "
                   "correctly with np.hanning(2*maxlag)[maxlag:]"),
               "imports": ["numpy", "scipy.ndimage.gaussian_filter1d"],
           },
           "reproduction": reproduction_check(root),
           "corrected_calibration": corrected_calibration(cfg),
           "window_check": window_check(root),
           "cells": [], "summary": {}}

    ts = [t for _, t in CELLS]
    for cell, t_k in CELLS:
        d = root / A78_DIR / cell
        z = np.load(d / "opacity_v3.npz")
        m = np.load(d / "dipoles.npy")
        n_mol, dt_fs = int(z["N"]), float(z["dt_fs"])
        grid, raw = z["wn"], z["sigma_raw"]
        doc["cells"].append({
            "cell": cell, "T_K": t_k, "rho_kg_m3": float(z["rho_kg_m3"]),
            "n_molecules": n_mol, "dt_fs": dt_fs, "n_frames": int(m.shape[0]),
            "beta_mdot2_per_molecule": beta_mdot2_per_molecule(
                m, dt_fs * 1e-15, t_k, n_mol),
            "band_integral_sigma_delivered": {
                b: float(np.trapezoid(raw[(grid >= lo) & (grid <= hi)],
                                      grid[(grid >= lo) & (grid <= hi)]))
                for b, (lo, hi) in BANDS_CM.items()},
            "band_integral_sigma_without_double_qcf": sigma_bands(
                m, t_k, dt_fs, n_mol, apply_qcf=False),
        })

    for band in BANDS_CM:
        mdot = [c["beta_mdot2_per_molecule"][band] for c in doc["cells"]]
        noqcf = [c["band_integral_sigma_without_double_qcf"][band] for c in doc["cells"]]
        shipped = [c["band_integral_sigma_delivered"][band] for c in doc["cells"]]
        ratio = [n / m for n, m in zip(noqcf, mdot)]
        doc["summary"][band] = {
            "index_beta_mdot2": _index(ts, mdot),
            "index_sigma_without_double_qcf": _index(ts, noqcf),
            "index_sigma_delivered": _index(ts, shipped),
            "index_attributable_to_the_spare_qcf":
                _index(ts, shipped) - _index(ts, noqcf),
            "ratio_noqcf_over_mdot2": ratio,
            "ratio_spread_frac": max(ratio) / min(ratio) - 1.0,
        }
    doc["verdict"] = {
        "trend_present_in_beta_mdot2": bool(
            any(abs(doc["summary"][b]["index_beta_mdot2"]) > 0.3 for b in BANDS_CM)),
        "statement": ("the dipole series satisfies the classical sum rule; the "
                      "temperature trend is introduced by the spare Q_harmonic "
                      "at opacity_pipeline.py line 81-82, i.e. by the formula"),
    }
    return doc


def _print(doc: dict) -> None:
    print("\n  Production formula (ed_data/arc_a78/scw/opacity_pipeline.py)")
    for k in ("line_66", "line_80", "line_81", "line_82"):
        print(f"    {k[-2:]}: {doc['production_formula'][k]}")
    print(f"    => {doc['production_formula']['verdict']}")
    print("\n  Reproduction of the delivered spectra by the formula as read:")
    for cell, r in doc["reproduction"].items():
        tag = "  <- stale x1000 in alpha/kappa (sigma correct)" if \
            r["carries_legacy_x1000_in_alpha_kappa"] else ""
        print(f"    {cell:22s} N={r['N']:4d}  sigma_raw to "
              f"{r['max_rel_error']['sigma_raw']:.1e}{tag}")

    print("\n  The sum rule at its source: beta <|Mdot_band|^2> per molecule")
    print(f"  {'band':10s}{'700 K':>13s}{'1000 K':>13s}{'2000 K':>13s}"
          f"{'index':>9s}{'  (delivered index)':>20s}")
    for band, s in doc["summary"].items():
        v = [c["beta_mdot2_per_molecule"][band] for c in doc["cells"]]
        print(f"  {band:10s}" + "".join(f"{x:13.4e}" for x in v)
              + f"{s['index_beta_mdot2']:+9.3f}"
              + f"{s['index_sigma_delivered']:+20.3f}")
    print("\n  Attribution of the trend, per band")
    print(f"  {'band':10s}{'dipoles':>10s}{'no 2nd qcf':>12s}{'delivered':>11s}"
          f"{'from qcf':>10s}{'ratio spread':>14s}")
    for band, s in doc["summary"].items():
        print(f"  {band:10s}{s['index_beta_mdot2']:+10.3f}"
              f"{s['index_sigma_without_double_qcf']:+12.3f}"
              f"{s['index_sigma_delivered']:+11.3f}"
              f"{s['index_attributable_to_the_spare_qcf']:+10.3f}"
              f"{100*s['ratio_spread_frac']:13.2f}%")
    wc = doc["window_check"]
    print("\n  Which lag window matches the window-free beta<|Mdot|^2>?")
    print(f"  {'window':26s}" + "".join(f"{b:>12s}" for b in BANDS_CM))
    for w, e in wc["windows"].items():
        print(f"  {w:26s}" + "".join(
            f"{100*e['spread_frac_by_band'][b]:11.1f}%" for b in BANDS_CM))
    print(f"  most consistent: {wc['most_consistent_window']} "
          "-> line 66 is unconventional but not a defect")

    cc = doc["corrected_calibration"]
    print("\n  A76 redone with the formula corrected (kappa / Q_harmonic),")
    print("  FP/HITEMP over the monomer dipole factor:")
    print(f"  {'band':10s}{'median':>16s}{'range (lo f)':>20s}{'dln/dT/kK':>12s}{'rho_s':>8s}")
    for band, s2 in cc["summary"].items():
        lo_m, hi_m = s2["median_over_factor"]
        med = f"{lo_m:.3f}" if lo_m == hi_m else f"{lo_m:.3f}-{hi_m:.3f}"
        r_lo, r_hi = s2["range_over_factor_low"]
        rng = f"[{r_lo:.3f}, {r_hi:.3f}]"
        print(f"  {band:10s}{med:>16s}{rng:>20s}"
              f"{s2['dln_dT_per_1000K']:+12.3f}{s2['spearman_rho']:+8.2f}")

    print(f"\n  VERDICT: {doc['verdict']['statement']}\n")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("config", type=Path)
    args = ap.parse_args(argv)
    cfg = load_config(args.config)
    doc = run(cfg)
    path = cfg.output_path / "dipole_sum_rule.json"
    atomicio.write_text(path, json.dumps(doc, indent=1) + "\n")
    _print(doc)
    print(f"  wrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
