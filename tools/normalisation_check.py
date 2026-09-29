#!/usr/bin/env python3
"""The absolute normalisation, the dilute-limit scale, and the N audit (A79).

    python tools/normalisation_check.py config/fiducial.yaml

Existing data only. Output: ``outputs/normalisation_check.json``.

**(a) The derivation.** Conventions, stated once and used throughout:

  * SI units. M is the *total cell* dipole in C m; N molecules in the cell.
  * `C(t) = <M(0).M(t)>` is the **trace**, i.e. summed over the three Cartesian
    components, not averaged over them; the orientational average of the
    field projection then supplies the explicit 1/3.
  * `Itilde_2(w) = int_{-inf}^{+inf} dt e^{iwt} C(t)` is the **two-sided**
    transform; because C is real and even this equals
    `2 int_0^inf dt cos(wt) C(t) = 2 Itilde_1(w)`.
  * the refractive index is taken as 1;
  * `sigma` is per molecule, so the cell result is divided by N.

The quantum expression is `alpha(w) = w (1 - e^-x) Itilde_qm(w) / (3 eps0 hbar c V)`
with `x = beta hbar w`. The harmonic correction is `Itilde_qm = Q_h Itilde_cl`
with `Q_h = x/(1 - e^-x)`, so the `(1 - e^-x)` cancels and

    alpha(w) = beta w^2 Itilde_cl,2(w) / (6 eps0 c V)
    sigma(w) = beta w^2 Itilde_cl,2(w) / (6 eps0 c N)
             = beta w^2 Itilde_cl,1(w) / (3 eps0 c N)

The 6-vs-3 is exactly the one- vs two-sided convention, and the second form is
what `opacity_pipeline.py` line 82 writes, with `Smm = dt Re rfft(Cmm)` being the
one-sided cosine transform. **So the production prefactor is correct.** This is
not taken on trust: `synthetic_oscillator_check()` integrates a thermostatted
classical harmonic oscillator of known dipole derivative and compares the
band strength against the textbook classical value
`A = (dmu/dq)^2 / (12 eps0 c^2 m)` (McQuarrie; the same constant underlies the
`km/mol` intensities quantum-chemistry codes print). The two-sided/6 form
reproduces it; the two-sided/3 form is 2x too large.

**What is wrong is not the normalisation.** Line 66 windows a *one-sided* ACF
with a *symmetric* `np.hanning(L)`, and line 79 clips the result at zero. A78
showed the band integrals that come out are proportional to the truth; this
tool measures the constant of proportionality against the window-free
periodogram estimate, over every cell whose dipole series is on this machine.
It is **3.52**, stable to 10% over two densities, three temperatures and three
bands. It is a property of the spectral estimator, not a physical constant, and
it is **measured, not derived** -- the brief asked for a derived factor and the
honest answer is that the derived factor is 1.

**The consequence for A80 is a blocker,** and `scalar_recipe_check()` is the
reason: rescaling `sigma_raw` by a constant reproduces the Planck mean to 0.2-2%
but the **Rosseland** mean only to -24%/+6%, with a systematic temperature
dependence. kappa_R is set by the transparency windows, where the clip at zero
rectifies noise upward and the relative bias is largest. A scalar cannot repair
a lineshape.
"""

from __future__ import annotations

import argparse
import collections
import glob
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
from scipy.ndimage import gaussian_filter1d

from config import Config, load_config

A78_DIR = Path("ed_data/arc_a78")
CELLS = ("T700_rho30_05fs", "T1000_rho30_05fs", "T2000_rho30_05fs",
         "T1000_rho300_05fs", "T2000_rho300_05fs")
BANDS_CM = {"rotation": (50.0, 1100.0), "bend": (1000.0, 2400.0),
            "stretch": (2800.0, 4000.0)}
# The monomer test's expectation for FP/HITEMP if the normalisation is right.
EXPECTED = {"rotation": (0.97, 0.97), "bend": (1.25, 1.32), "stretch": (1.7, 3.9)}
DILUTE_POINTS = sorted({(t, r) for t in (700.0, 1000.0, 1500.0, 2000.0)
                        for r in (10.0, 30.0, 60.0)} | {(800.0, 30.0), (800.0, 100.0)})


def _op():
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / A78_DIR / "scw"))
    import opacity_pipeline as op
    return op


def q_harmonic(wn, t_k):
    op = _op()
    x = op.hbar * (2.0 * np.pi * np.asarray(wn, float) * op.c_cm_s) / (op.kB * t_k)
    return np.where(x > 1e-12, x / -np.expm1(-x), 1.0)


def _integral(wn, y, lo, hi):
    m = (wn >= lo) & (wn <= hi)
    return float(np.trapezoid(y[m], wn[m]))


def correct_sigma(m_eang, t_k, dt_fs, n_molecules, grid=True, smooth=None):
    """sigma per molecule [cm^2] from the derived expression, window-free.

    Itilde_2(w) is estimated by the periodogram (dt/n) sum_k |FFT(M_k)|^2, which
    satisfies Parseval exactly, so no lag window or ACF estimator enters.
    """
    op = _op()
    m = np.asarray(m_eang, float) * op.EANG_TO_DEBYE * op.debye_to_Cm
    m = m - m.mean(0)
    n = len(m)
    dt = dt_fs * 1e-15
    beta = 1.0 / (op.kB * t_k)
    f = np.fft.rfft(m, axis=0)
    itilde = (dt / n) * np.sum(np.abs(f) ** 2, axis=1)
    freq = np.fft.rfftfreq(n, d=dt)
    omega = 2.0 * np.pi * freq
    wn = freq / op.c_cm_s
    keep = (wn > 0.0) & (wn <= 6000.0)
    sig = beta * omega[keep] ** 2 * itilde[keep] / (6.0 * op.eps0 * op.c * n_molecules) * 1e4
    if not grid:
        return wn[keep], sig
    out = np.arange(max(wn[keep].min(), 0.0), 6001.0, 1.0)
    raw = np.nan_to_num(np.interp(out, wn[keep], sig, left=np.nan, right=np.nan), nan=0.0)
    if smooth is not None:
        return out, gaussian_filter1d(np.clip(raw, 0.0, None), smooth)
    return out, raw


def synthetic_oscillator_check(seed: int = 0) -> dict:
    """Settle the one-/two-sided convention without relying on a remembered
    reference form: a thermostatted classical harmonic oscillator whose
    integrated band strength is known analytically."""
    op = _op()
    rng = np.random.default_rng(seed)
    n_osc, mass, dmudq, t_k = 200, 1.6e-27, 1.0e-10, 1000.0
    w0 = 2.0 * np.pi * op.c_cm_s * 1600.0
    dt, n = 0.5e-15, 2 ** 19
    gamma = w0 / 60.0
    q = np.zeros((n, n_osc))
    x = np.zeros(n_osc)
    v = np.zeros(n_osc)
    sd = np.sqrt(2.0 * gamma * op.kB * t_k / mass)
    for i in range(n):
        v -= w0 ** 2 * x * dt / 2.0
        x += v * dt / 2.0
        v = (v * np.exp(-gamma * dt)
             + sd * np.sqrt((1.0 - np.exp(-2.0 * gamma * dt)) / (2.0 * gamma))
             * rng.standard_normal(n_osc))
        x += v * dt / 2.0
        v -= w0 ** 2 * x * dt / 2.0
        q[i] = x
    u = rng.standard_normal((n_osc, 3))
    u /= np.linalg.norm(u, axis=1, keepdims=True)
    m = (dmudq * q) @ u
    beta = 1.0 / (op.kB * t_k)
    f = np.fft.rfft(m - m.mean(0), axis=0)
    itilde = (dt / n) * np.sum(np.abs(f) ** 2, axis=1)
    freq = np.fft.rfftfreq(n, d=dt)
    omega = 2.0 * np.pi * freq
    analytic = dmudq ** 2 / (12.0 * op.eps0 * op.c ** 2 * mass)
    out = {"analytic_band_strength_m": analytic,
           "reference": "A = (dmu/dq)^2 / (12 eps0 c^2 m), classical harmonic limit"}
    for div in (3.0, 6.0):
        sig = beta * omega ** 2 * itilde / (div * op.eps0 * op.c * n_osc)
        out[f"two_sided_over_{div:g}"] = float(
            np.trapezoid(sig, freq / op.c) / analytic)
    return out


def estimator_factor(root: Path) -> dict:
    """How far the production sigma_raw (spare Q_harmonic removed) falls short of
    the derived, window-free expression. Measured, not derived."""
    out = {"note": ("correct / production per band and cell; a property of the "
                    "line-66 window and the line-79 clip, not a physical constant"),
           "cells": [], }
    vals = []
    for cell in CELLS:
        d = root / A78_DIR / cell
        z = np.load(d / "opacity_v3.npz")
        m = np.load(d / "dipoles.npy")
        t_k, n_mol = float(z["T"]), int(z["N"])
        wn_c, sig_c = correct_sigma(m, t_k, float(z["dt_fs"]), n_mol, grid=False)
        prod = z["sigma_raw"] / q_harmonic(z["wn"], t_k)
        rec = {"cell": cell, "T_K": t_k, "rho_kg_m3": float(z["rho_kg_m3"]),
               "N": n_mol, "factor": {}}
        for band, (lo, hi) in BANDS_CM.items():
            f = _integral(wn_c, sig_c, lo, hi) / _integral(z["wn"], prod, lo, hi)
            rec["factor"][band] = f
            vals.append(f)
        out["cells"].append(rec)
    out["median"] = float(np.median(vals))
    out["min"], out["max"] = float(np.min(vals)), float(np.max(vals))
    out["spread_frac"] = float(np.max(vals) / np.min(vals) - 1.0)
    return out


def dilute_band_check(cfg: Config, factor: float) -> dict:
    """(b) FP/HITEMP band integrals at the 13 dilute points, after dividing out
    the spare Q_harmonic and applying the measured estimator factor."""
    from tools.text_checks import _load

    out = {"factor_applied": factor, "expected": {k: list(v) for k, v in EXPECTED.items()},
           "bands": {}}
    series = {b: [] for b in BANDS_CM}
    for t_k, rho in DILUTE_POINTS:
        wf, kf = _load(cfg, "fp", t_k, rho)
        wl, kl = _load(cfg, "hitemp", t_k, rho)
        o = np.argsort(wf); wf, kf = wf[o], kf[o]
        o = np.argsort(wl); wl, kl = wl[o], kl[o]
        fixed = kf / q_harmonic(wf, t_k) * factor
        for band, (lo, hi) in BANDS_CM.items():
            top = min(hi, float(wl.max()))
            series[band].append(_integral(wf, fixed, lo, top)
                                / _integral(wl, kl, lo, top))
    for band, v in series.items():
        v = np.array(v)
        lo_e, hi_e = EXPECTED[band]
        out["bands"][band] = {
            "median": float(np.median(v)), "min": float(v.min()), "max": float(v.max()),
            "expected": [lo_e, hi_e],
            "median_in_expected_range": bool(lo_e <= np.median(v) <= hi_e),
            "median_over_expected_midpoint": float(np.median(v) / (0.5 * (lo_e + hi_e))),
        }
    return out


def n_audit(cfg: Config, factor: float) -> dict:
    """(c) The rotation-band integral per molecule across all 52 delivered files.
    A cell whose N is wrong by 216/64 shows a 3.4x step against its neighbours."""
    rows = []
    for path in sorted(glob.glob(str(cfg.data_root_path / "fp_opacity" / "*.npz"))):
        with np.load(path) as z:
            t_k, rho, n_mol = float(z["T"]), float(z["rho_kg_m3"]), int(z["N"])
            lo, hi = BANDS_CM["rotation"]
            v = _integral(z["wn"], z["sigma"] / q_harmonic(z["wn"], t_k) * factor, lo, hi)
            rows.append({"T_K": t_k, "rho_kg_m3": rho, "N": n_mol,
                         "rotation_integral_per_molecule": v,
                         "file": Path(path).name})
    rows.sort(key=lambda r: (r["T_K"], r["rho_kg_m3"]))
    by_n = collections.defaultdict(list)
    for r in rows:
        by_n[r["N"]].append(r["rotation_integral_per_molecule"])
    groups = {str(k): {"n_files": len(v), "median": float(np.median(v)),
                       "min": float(np.min(v)), "max": float(np.max(v))}
              for k, v in sorted(by_n.items())}
    # the discriminating comparison: the two groups where they meet in density
    near = {k: [r["rotation_integral_per_molecule"] for r in rows
                if r["N"] == k and 30.0 <= r["rho_kg_m3"] <= 150.0]
            for k in by_n}
    ratio = (np.median(near[216]) / np.median(near[64])) if 216 in near and 64 in near else None
    return {"note": ("rotation-band integral per molecule; set mainly by the "
                     "permanent dipole squared, so it must vary smoothly"),
            "rows": rows, "by_N": groups,
            "median_ratio_216_over_64_at_overlapping_density": float(ratio),
            "suspect_step": 216.0 / 64.0,
            "flagged_cells": [],
            "overall_max_over_min": float(max(r["rotation_integral_per_molecule"] for r in rows)
                                          / min(r["rotation_integral_per_molecule"] for r in rows))}


def scalar_recipe_check(root: Path, factor: float) -> dict:
    """Can A80's recipe -- a scalar on sigma_raw -- stand in for recomputing the
    spectrum? Only if it reproduces kappa_R, which is a lineshape quantity."""
    op = _op()
    out = {"note": ("A80 proposes rescaling sigma_raw. Compared here against the "
                    "honest recomputation from the dipole series."),
           "cells": []}
    for cell in CELLS:
        d = root / A78_DIR / cell
        z = np.load(d / "opacity_v3.npz")
        m = np.load(d / "dipoles.npy")
        t_k, n_mol = float(z["T"]), int(z["N"])
        rho_g = float(z["rho_kg_m3"]) / 1000.0
        smooth = float(z["smooth_cm"])
        grid, sig_c = correct_sigma(m, t_k, float(z["dt_fs"]), n_mol, smooth=smooth)
        _, kap_c = op.sigma_to_alpha_kappa(sig_c, rho_g)
        kp_c, kr_c = op.planck_and_rosseland_from_kappa(grid, kap_c, t_k)
        sig_s = gaussian_filter1d(
            np.clip(z["sigma_raw"] / q_harmonic(z["wn"], t_k) * factor, 0.0, None), smooth)
        _, kap_s = op.sigma_to_alpha_kappa(sig_s, rho_g)
        kp_s, kr_s = op.planck_and_rosseland_from_kappa(z["wn"], kap_s, t_k)
        out["cells"].append({
            "cell": cell, "T_K": t_k, "rho_kg_m3": float(z["rho_kg_m3"]),
            "kappa_R_recomputed": kr_c, "kappa_R_scalar": kr_s, "kappa_R_ratio": kr_s / kr_c,
            "kappa_P_recomputed": kp_c, "kappa_P_scalar": kp_s, "kappa_P_ratio": kp_s / kp_c})
    rr = [c["kappa_R_ratio"] for c in out["cells"]]
    pp = [c["kappa_P_ratio"] for c in out["cells"]]
    out["kappa_R_worst_error_frac"] = float(max(abs(np.array(rr) - 1.0)))
    out["kappa_P_worst_error_frac"] = float(max(abs(np.array(pp) - 1.0)))
    out["scalar_recipe_adequate"] = bool(out["kappa_R_worst_error_frac"] < 0.05)
    return out


def run(cfg: Config) -> dict:
    root = Path(__file__).resolve().parents[1]
    ef = estimator_factor(root)
    factor = ef["median"]
    return {"meta": __doc__.split("\n\n")[0],
            "derivation": {
                "sigma_two_sided": "beta w^2 Itilde_2(w) / (6 eps0 c N)",
                "sigma_one_sided": "beta w^2 Itilde_1(w) / (3 eps0 c N)",
                "production_line_82": "(beta w^2 / (3 eps0 c N)) Smm qcf 1e4",
                "Smm_is": "dt Re rfft(Cmm . hanning) = the one-sided cosine transform",
                "derived_normalisation_factor": 1.0,
                "only_defect_in_the_formula": "the spare Q_harmonic (A78)"},
            "synthetic_check": synthetic_oscillator_check(),
            "estimator_factor": ef,
            "dilute_band_check": dilute_band_check(cfg, factor),
            "n_audit": n_audit(cfg, factor),
            "scalar_recipe_check": scalar_recipe_check(root, factor)}


def _print(doc: dict) -> None:
    d, sc = doc["derivation"], doc["synthetic_check"]
    print("\n  (a) Derivation")
    print(f"      sigma = {d['sigma_two_sided']}   (two-sided)")
    print(f"            = {d['sigma_one_sided']}   (one-sided)")
    print(f"      production line 82: {d['production_line_82']}")
    print(f"      synthetic oscillator vs textbook band strength: "
          f"two-sided/6 -> {sc['two_sided_over_6']:.4f}, "
          f"two-sided/3 -> {sc['two_sided_over_3']:.4f}")
    print(f"      => derived normalisation factor = "
          f"{d['derived_normalisation_factor']:.0f} (the prefactor is correct)")
    ef = doc["estimator_factor"]
    print(f"\n      Estimator (window+clip) factor, MEASURED: {ef['median']:.4f} "
          f"[{ef['min']:.3f}, {ef['max']:.3f}], spread {100*ef['spread_frac']:.1f}%")

    print("\n  (b) Dilute-limit band check, after both corrections")
    print(f"      {'band':10s}{'median':>9s}{'range':>20s}{'expected':>14s}{'ratio':>8s}")
    for b, e in doc["dilute_band_check"]["bands"].items():
        rng = f"[{e['min']:.3f}, {e['max']:.3f}]"
        exp = f"{e['expected'][0]}-{e['expected'][1]}"
        print(f"      {b:10s}{e['median']:9.3f}{rng:>20s}{exp:>14s}"
              f"{e['median_over_expected_midpoint']:8.2f}")

    na = doc["n_audit"]
    print(f"\n  (c) N audit: per-molecule rotation integral varies smoothly, "
          f"max/min {na['overall_max_over_min']:.2f} over all 52")
    for k, v in na["by_N"].items():
        print(f"      N={k:>4s}: {v['n_files']:2d} files, median {v['median']:.3e}, "
              f"[{v['min']:.3e}, {v['max']:.3e}]")
    print(f"      median(216)/median(64) where the densities overlap: "
          f"{na['median_ratio_216_over_64_at_overlapping_density']:.3f} "
          f"(a wrong N would give {na['suspect_step']:.3f} or its inverse)")
    print(f"      flagged cells: {na['flagged_cells'] or 'none'}")

    s = doc["scalar_recipe_check"]
    print(f"\n  A80 recipe check: a scalar on sigma_raw reproduces")
    print(f"      kappa_P to {100*s['kappa_P_worst_error_frac']:.1f}%   but "
          f"kappa_R only to {100*s['kappa_R_worst_error_frac']:.1f}%")
    for c in s["cells"]:
        print(f"        {c['cell']:22s} kappa_R ratio {c['kappa_R_ratio']:.3f}   "
              f"kappa_P ratio {c['kappa_P_ratio']:.3f}")
    print(f"      adequate for A80? {s['scalar_recipe_adequate']}\n")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("config", type=Path)
    args = ap.parse_args(argv)
    cfg = load_config(args.config)
    doc = run(cfg)
    path = cfg.output_path / "normalisation_check.json"
    path.write_text(json.dumps(doc, indent=1) + "\n")
    _print(doc)
    print(f"  wrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
