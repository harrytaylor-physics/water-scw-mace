#!/usr/bin/env python3
"""Regenerate the first-principles spectra with the fixed estimator (A80).

    python tools/v4_regenerate.py config/fiducial.yaml            # gates only
    python tools/v4_regenerate.py config/fiducial.yaml --write    # + write v4

Output: ``outputs/v4_regenerate.json``.

A78 found the harmonic quantum correction applied twice; A79 showed the
normalisation was right and the remaining factor of 3.52 was the spectral
estimator (a symmetric Hanning on a one-sided ACF, then a clip at zero), and
that no scalar can repair it because kappa_R is a lineshape quantity.
``opacity.dipole_spectrum`` replaces the estimator outright: a window-free
periodogram, the harmonic correction entering once, no clip, and everything
downstream of the estimator kept exactly as delivered.

**The brief gates this before propagating, and one gate fails.** The A75 dilute
noise floor, recomputed through the fixed estimator, is no longer far below the
transparency-window opacity: at the cold dilute points it is within a factor of
a few, and at (1000 K, 10 kg m^-3) it is 1.8x below the deepest point of the
window and 55% of it. The v3 clip at zero had been suppressing the windows,
which is the same defect that made kappa_R unrepairable by a scalar; removing it
raises the windows and exposes how much of what is left is dipole-model noise.
So this tool stops at the gates: it does not switch the registry and it reruns
nothing downstream. ``--write`` puts the v4 spectra in their own directory so
the result can be inspected, which changes no delivered product.
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

import opacity.dipole_spectrum as ds
from config import Config, load_config
from opacity.dipole_spectrum import K_B, spectrum_from_dipoles, two_sided_power_spectrum

DIPOLE_DIR = Path("ed_data/arc_dipoles")
V4_SUBDIR = "fp_opacity_v4"
UNITS_MARKER = "units_fixed"
UNITS_VALUE = "cm2_per_g:A80"
VERSION_VALUE = "v4:fixed-estimator"

BANDS_CM = {"rotation": (50.0, 1100.0), "bend": (1000.0, 2400.0),
            "stretch": (2800.0, 4000.0)}
SMOOTH_SENSITIVITY = ((700.0, 30.0), (1000.0, 100.0), (1500.0, 300.0),
                      (2000.0, 1200.0), (2500.0, 500.0), (3000.0, 2228.71))
# A75's uncorrelated dipole-residual power fraction, from the lag test.
R2_UNCORRELATED = 7.26e-5
WINDOW_CM = (1900.0, 3100.0)
MARK_CM = 2550.0
FLOOR_MARGIN_REQUIRED = 10.0


def manifest(cfg: Config, root: Path) -> dict:
    """Match every delivered cell to its dipole series, by reproduction.

    Fifteen tags have two directories -- an older ~30002-frame run and an
    ``_05fs`` 60002-frame one. Rather than trust the names, each candidate is
    pushed through the *unfixed* production backend and the one that reproduces
    the delivered sigma_raw wins.
    """
    sys.path.insert(0, str(root / "ed_data/arc_a78/scw"))
    import opacity_pipeline as prod

    groups: dict[tuple[float, float], list[str]] = {}
    for d in sorted(os.listdir(root / DIPOLE_DIR)):
        m = re.match(r"T(\d+(?:\.\d+)?)_rho([\d.]+?)(_05fs)?$", d)
        if m:
            groups.setdefault((float(m.group(1)), float(m.group(2))), []).append(d)

    rows, missing = [], []
    for path in sorted((cfg.data_root_path / "fp_opacity").glob("T*_opacity.npz")):
        with np.load(path) as z:
            t_k, rho = float(z["T"]), float(z["rho_kg_m3"])
            n_mol, dt_fs = int(z["N"]), float(z["dt_fs"])
            smooth = float(z["smooth_cm"])
            stored = z["sigma_raw"].copy()
        cands = groups.get((t_k, rho), [])
        if not cands:
            missing.append({"T_K": t_k, "rho_kg_m3": rho})
            continue
        best = None
        for d in cands:
            m = np.load(root / DIPOLE_DIR / d / "dipoles.npy")
            r = prod.spectrum_from_dipoles(m, T_K=t_k, dt_fs=dt_fs, n_molecules=n_mol,
                                           density_g_cm3=rho / 1000.0, units="eAng",
                                           smooth_sigma_cm=smooth)
            got = r["sigma_cm2_per_molecule_raw"]
            if got.shape != stored.shape:
                continue
            ok = stored != 0
            err = float(np.max(np.abs(got[ok] / stored[ok] - 1.0)))
            if best is None or err < best[1]:
                best = (d, err)
        rows.append({"T_K": t_k, "rho_kg_m3": rho, "N": n_mol, "dt_fs": dt_fs,
                     "smooth_cm": smooth, "series": best[0],
                     "reproduction_max_rel_error": best[1],
                     "n_candidates": len(cands)})
    return {"n_delivered": len(rows) + len(missing), "n_matched": len(rows),
            "missing": missing, "rows": rows,
            "worst_reproduction_error": max(r["reproduction_max_rel_error"] for r in rows),
            "note": ("all matches resolve to the _05fs 60002-frame series; the "
                     "15 bare directories are older runs and are unused")}


def _series(root: Path, man: dict, t_k: float, rho: float):
    for r in man["rows"]:
        if r["T_K"] == t_k and r["rho_kg_m3"] == rho:
            return np.load(root / DIPOLE_DIR / r["series"] / "dipoles.npy"), r
    raise KeyError((t_k, rho))


def gate_oscillator() -> dict:
    """(2)(i) The A79 harmonic-oscillator check, through the v4 estimator."""
    rng = np.random.default_rng(0)
    n_osc, mass, dmudq, t_k = 200, 1.6e-27, 1.0e-10, 1000.0
    w0 = 2.0 * np.pi * ds.C_CM * 1600.0
    dt, n = 0.5e-15, 2 ** 19
    gamma = w0 / 60.0
    x = np.zeros(n_osc); v = np.zeros(n_osc); q = np.zeros((n, n_osc))
    sd = np.sqrt(2.0 * gamma * K_B * t_k / mass)
    for i in range(n):
        v -= w0 ** 2 * x * dt / 2.0; x += v * dt / 2.0
        v = (v * np.exp(-gamma * dt) + sd
             * np.sqrt((1.0 - np.exp(-2.0 * gamma * dt)) / (2.0 * gamma))
             * rng.standard_normal(n_osc))
        x += v * dt / 2.0; v -= w0 ** 2 * x * dt / 2.0
        q[i] = x
    u = rng.standard_normal((n_osc, 3)); u /= np.linalg.norm(u, axis=1, keepdims=True)
    m_cm = (dmudq * q) @ u
    omega, itilde = two_sided_power_spectrum(m_cm / (ds.EANG_TO_DEBYE * ds.DEBYE_TO_CM), dt)
    sig = (1.0 / (K_B * t_k)) * omega ** 2 * itilde / (6.0 * ds.EPS0 * ds.C_VAC * n_osc)
    numeric = float(np.trapezoid(sig, omega / (2.0 * np.pi) / ds.C_VAC))
    analytic = dmudq ** 2 / (12.0 * ds.EPS0 * ds.C_VAC ** 2 * mass)
    ratio = numeric / analytic
    return {"numeric_over_analytic": ratio, "tolerance": 0.03,
            "passes": bool(abs(ratio - 1.0) < 0.03)}


def gate_band_integrals(root: Path, man: dict) -> dict:
    """(2)(ii) Each band integral equals beta<|Mdot_band|^2> to 2%.

    Mdot is taken **spectrally** (multiply by i.omega), not by central
    differences: the central-difference estimator carries a known
    [sin(w dt)/(w dt)]^2 bias worth 5% at the top of the stretch band, which is
    a property of the test and not of the estimator under test.
    """
    pred_const = 1e4 / (12.0 * ds.EPS0 * ds.C_VAC * ds.C_CM)
    out = {"derived_constant": pred_const, "tolerance": 0.02, "cells": []}
    worst = 0.0
    for t_k in (700.0, 1000.0, 2000.0):
        m, meta = _series(root, man, t_k, 30.0)
        r = spectrum_from_dipoles(m, t_k, meta["dt_fs"], meta["N"], 0.03,
                                  smooth_sigma_cm=meta["smooth_cm"])
        mx = np.asarray(m, float) * ds.EANG_TO_DEBYE * ds.DEBYE_TO_CM
        mx = mx - mx.mean(0)
        n = len(mx)
        dts = meta["dt_fs"] * 1e-15
        omega = 2.0 * np.pi * np.fft.rfftfreq(n, d=dts)
        wn = omega / (2.0 * np.pi) / ds.C_CM
        rec = {"T_K": t_k, "rho_kg_m3": 30.0, "ratio": {}}
        for band, (lo, hi) in BANDS_CM.items():
            mask = (wn >= lo) & (wn <= hi)
            total = 0.0
            for k in range(3):
                f = np.fft.rfft(mx[:, k]) * (1j * omega)
                f[~mask] = 0.0
                total += float(np.mean(np.fft.irfft(f, n=n) ** 2))
            mdot = (1.0 / (K_B * t_k)) * total / meta["N"]
            sel = (r["wn"] >= lo) & (r["wn"] <= hi)
            got = float(np.trapezoid(r["sigma_raw"][sel], r["wn"][sel]))
            rec["ratio"][band] = got / (mdot * pred_const)
            worst = max(worst, abs(rec["ratio"][band] - 1.0))
        out["cells"].append(rec)
    out["worst_deviation"] = worst
    out["passes"] = bool(worst < 0.02)
    return out


def gate_noise_floor(root: Path, man: dict) -> dict:
    """(2)(iii) The A75 floor under the fixed estimator, against the window.

    A white residual carrying a fraction R2 of the dipole power has a flat
    two-sided spectrum Itilde = R2 <|M|^2> dt, which goes through the *same*
    estimator as the signal. Compared pointwise -- at 2550 cm^-1, which is where
    A75 quoted its margins, and at the deepest point of the 1900-3100 window.
    """
    out = {"r2": R2_UNCORRELATED, "required_margin": FLOOR_MARGIN_REQUIRED,
           "points": []}
    for row in man["rows"]:
        t_k, rho = row["T_K"], row["rho_kg_m3"]
        if rho > 60.0 and not (t_k == 800.0 and rho == 100.0):
            continue
        m, meta = _series(root, man, t_k, rho)
        r = spectrum_from_dipoles(m, t_k, meta["dt_fs"], meta["N"], rho / 1000.0,
                                  smooth_sigma_cm=meta["smooth_cm"])
        mx = np.asarray(m, float) * ds.EANG_TO_DEBYE * ds.DEBYE_TO_CM
        mx = mx - mx.mean(0)
        var = float(np.sum(np.mean(mx ** 2, axis=0)))
        itilde = R2_UNCORRELATED * var * meta["dt_fs"] * 1e-15
        wn = r["wn"]
        omega = 2.0 * np.pi * wn * ds.C_CM
        sig_f = ((1.0 / (K_B * t_k)) * omega ** 2 * itilde
                 / (6.0 * ds.EPS0 * ds.C_VAC * meta["N"]) * 1e4)
        _, kap_f = ds.sigma_to_alpha_kappa(sig_f, rho / 1000.0)
        sel = (wn >= WINDOW_CM[0]) & (wn <= WINDOW_CM[1])
        i_mark = int(np.argmin(np.abs(wn - MARK_CM)))
        i_min = int(np.where(sel)[0][np.argmin(r["kappa"][sel])])
        # The lag-test r2 was measured at (800, 61) and (1000, 130). The cells
        # here are more dilute and noisier, so the spectrum's own plateau above
        # every water band is the better estimate of their floor -- an upper
        # bound at hot points, where real hot-band absorption lives up there,
        # but close to pure noise at the cold dilute ones.
        hi_m = (wn >= 4500.0) & (wn <= 5800.0)
        implied = R2_UNCORRELATED * float(np.median(r["sigma"][hi_m])
                                          / np.median(sig_f[hi_m]))
        out["points"].append({
            "T_K": t_k, "rho_kg_m3": rho,
            "margin_at_2550": float(r["kappa"][i_mark] / kap_f[i_mark]),
            "margin_at_window_minimum": float(r["kappa"][i_min] / kap_f[i_min]),
            "noise_fraction_at_window_minimum": float(kap_f[i_min] / r["kappa"][i_min]),
            "r2_implied_by_plateau": implied,
            "margin_at_window_minimum_using_plateau_r2": float(
                r["kappa"][i_min] / (kap_f[i_min] * implied / R2_UNCORRELATED))})
    a = [p["margin_at_2550"] for p in out["points"]]
    b = [p["margin_at_window_minimum"] for p in out["points"]]
    c = [p["margin_at_window_minimum_using_plateau_r2"] for p in out["points"]]
    out["worst_margin_using_plateau_r2"] = float(min(c))
    out["worst_margin_at_2550"] = float(min(a))
    out["worst_margin_at_window_minimum"] = float(min(b))
    out["passes"] = bool(min(a) >= FLOOR_MARGIN_REQUIRED
                         and min(b) >= FLOOR_MARGIN_REQUIRED)
    return out


def smoothing_sensitivity(root: Path, man: dict) -> dict:
    """kappa_R and kappa_P at sd = 10 and 40 against the delivered choice."""
    out = {"points": []}
    for t_k, rho in SMOOTH_SENSITIVITY:
        m, meta = _series(root, man, t_k, rho)
        rec = {"T_K": t_k, "rho_kg_m3": rho, "sd_delivered": meta["smooth_cm"]}
        for sd in (10.0, meta["smooth_cm"], 40.0):
            r = spectrum_from_dipoles(m, t_k, meta["dt_fs"], meta["N"], rho / 1000.0,
                                      smooth_sigma_cm=sd)
            rec[f"kappa_R_sd{sd:g}"] = r["rosseland"]
            rec[f"kappa_P_sd{sd:g}"] = r["planck"]
        rec["kappa_R_ratio_10_over_40"] = rec["kappa_R_sd10"] / rec["kappa_R_sd40"]
        rec["kappa_P_ratio_10_over_40"] = rec["kappa_P_sd10"] / rec["kappa_P_sd40"]
        out["points"].append(rec)
    out["kappa_R_worst_spread_frac"] = float(
        max(abs(p["kappa_R_ratio_10_over_40"] - 1.0) for p in out["points"]))
    out["kappa_P_worst_spread_frac"] = float(
        max(abs(p["kappa_P_ratio_10_over_40"] - 1.0) for p in out["points"]))
    return out


def write_v4(cfg: Config, root: Path, man: dict) -> dict:
    """Write the v4 spectra to their own directory. The registry is untouched."""
    dest = cfg.data_root_path / V4_SUBDIR
    dest.mkdir(parents=True, exist_ok=True)
    written = []
    for row in man["rows"]:
        t_k, rho = row["T_K"], row["rho_kg_m3"]
        m, meta = _series(root, man, t_k, rho)
        r = spectrum_from_dipoles(m, t_k, meta["dt_fs"], meta["N"], rho / 1000.0,
                                  smooth_sigma_cm=meta["smooth_cm"])
        out = dest / f"T{t_k:g}_rho{rho:g}_opacity.npz"
        np.savez(out, wn=r["wn"], sigma=r["sigma"], sigma_raw=r["sigma_raw"],
                 alpha=r["alpha"], kappa=r["kappa"], planck=r["planck"],
                 rosseland=r["rosseland"], T=t_k, rho_kg_m3=rho, N=meta["N"],
                 dt_fs=meta["dt_fs"], smooth_cm=meta["smooth_cm"],
                 **{UNITS_MARKER: np.array(UNITS_VALUE)},
                 version=np.array(VERSION_VALUE),
                 estimator=np.array("periodogram; sigma = beta w^2 I2/(6 eps0 c N); "
                                    "single harmonic correction; no zero clip"))
        written.append(out.name)
    # grid.json is a list of (T, rho) and carries over. crossover_hitemp.json
    # holds v3-DERIVED crossover points and would be silently wrong against v4,
    # so it is copied under a name nothing loads.
    src = cfg.data_root_path / "fp_opacity" / "grid.json"
    if src.exists():
        (dest / "grid.json").write_text(src.read_text())
    src = cfg.data_root_path / "fp_opacity" / "crossover_hitemp.json"
    if src.exists():
        (dest / "crossover_hitemp_v3_DO_NOT_USE.json").write_text(src.read_text())
    stale = dest / "crossover_hitemp.json"
    if stale.exists():
        stale.unlink()
    return {"directory": str(dest), "n_written": len(written),
            "units_marker": {UNITS_MARKER: UNITS_VALUE}, "version": VERSION_VALUE,
            "registry_switched": False,
            "note": "written for inspection only; nothing downstream reads this"}


def run(cfg: Config, write: bool) -> dict:
    root = Path(__file__).resolve().parents[1]
    man = manifest(cfg, root)
    doc = {"meta": __doc__.split("\n\n")[0], "manifest": man,
           "gates": {"oscillator": gate_oscillator(),
                     "band_integrals": gate_band_integrals(root, man),
                     "noise_floor": gate_noise_floor(root, man)},
           "smoothing_sensitivity": smoothing_sensitivity(root, man)}
    doc["gates_all_pass"] = all(g["passes"] for g in doc["gates"].values())
    doc["propagated"] = False
    if write:
        doc["v4_written"] = write_v4(cfg, root, man)
    return doc


def _print(doc: dict) -> None:
    m = doc["manifest"]
    print(f"\n  Manifest: {m['n_matched']}/{m['n_delivered']} delivered cells matched; "
          f"missing {len(m['missing'])}; worst reproduction "
          f"{m['worst_reproduction_error']:.1e}")
    g = doc["gates"]
    print(f"\n  Gates before propagating")
    o = g["oscillator"]
    print(f"    (i)   oscillator      {o['numeric_over_analytic']:.4f}  "
          f"{'PASS' if o['passes'] else 'FAIL'}")
    b = g["band_integrals"]
    print(f"    (ii)  band integrals  worst deviation "
          f"{100*b['worst_deviation']:.2f}%  {'PASS' if b['passes'] else 'FAIL'}")
    f = g["noise_floor"]
    print(f"    (iii) noise floor     worst margin {f['worst_margin_at_2550']:.1f}x "
          f"at 2550, {f['worst_margin_at_window_minimum']:.1f}x at the window "
          f"minimum  {'PASS' if f['passes'] else 'FAIL'}")
    print(f"          with the plateau-implied r2 instead: "
          f"{f['worst_margin_using_plateau_r2']:.1f}x")
    if not f["passes"]:
        bad = sorted(f["points"], key=lambda p: p["margin_at_window_minimum"])[:4]
        for p in bad:
            print(f"          T={p['T_K']:.0f} rho={p['rho_kg_m3']:g}: "
                  f"{p['margin_at_window_minimum']:.1f}x, "
                  f"{100*p['noise_fraction_at_window_minimum']:.0f}% of the window "
                  f"opacity is dipole noise")
    s = doc["smoothing_sensitivity"]
    print(f"\n  Smoothing sensitivity (sd 10 vs 40): kappa_R "
          f"{100*s['kappa_R_worst_spread_frac']:.1f}%, kappa_P "
          f"{100*s['kappa_P_worst_spread_frac']:.1f}%")
    print(f"\n  gates all pass: {doc['gates_all_pass']}; propagated: {doc['propagated']}")
    if "v4_written" in doc:
        w = doc["v4_written"]
        print(f"  wrote {w['n_written']} v4 spectra to {w['directory']} "
              f"(registry NOT switched)")
    print()


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("config", type=Path)
    ap.add_argument("--write", action="store_true",
                    help="write the v4 spectra to their own directory")
    args = ap.parse_args(argv)
    cfg = load_config(args.config)
    doc = run(cfg, args.write)
    path = cfg.output_path / "v4_regenerate.json"
    path.write_text(json.dumps(doc, indent=1) + "\n")
    _print(doc)
    print(f"  wrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
