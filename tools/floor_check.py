#!/usr/bin/env python3
"""Is there a noise floor under the dilute first-principles spectra?

    python tools/floor_check.py config/fiducial.yaml

Existing data only; no new physics is run and no existing file is written.
Output: ``outputs/floor_check.json``.

(a) Four spectral quantities at the 13 dilute points and at every vertex of the
    ten strongest boundaries' interpolation triangles, for both sources.
(b) The floor implied by the *uncorrelated* part of the dipole residual, taken
    from ``ed_data/arc/scw/lagtest_v2.json`` rather than from the A52 amplitude
    ratio, through the production normalisation.
(c) The same kappa_R with the lineshape quantum correction changed.

**What the delivered spectra actually carry**, read from the code rather than
the prose: ``ed_data/arc/water-scw-mace/convergence/04_spectra/
spectrum_convergence.py`` line 28 is

    alpha = w * (1 - np.exp(-HBAR*w/(KB*T))) * I

with ``I`` the classical dipole-autocorrelation spectrum. That is the
detailed-balance (stimulated-emission) factor and **no quantum correction on the
lineshape** — Q = 1 in the Ramirez et al. (2004) taxonomy. The ARC prose calls
the w(1-exp(-x)) factor itself "the quantum correction", which is where the name
"harmonic" in the request comes from; the operation asked for is unambiguous
either way, and both readings are reported below:

  Q_standard(x) = 2 / (1 + exp(-x))        -> 1 as x -> 0
  Q_harmonic(x) = x / (1 - exp(-x))
  Q_schofield(x) = exp(x/2)                 x = hc nu / kT

``as_delivered`` multiplies by Q itself (the delivered lineshape is classical);
``relative_to_harmonic`` multiplies by Q/Q_harmonic (the request's literal
reading, that the delivered spectra already carry the harmonic correction). The
three coincide as x -> 0 and separate in the stretch band.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import atomicio

import numpy as np

from config import Config, load_config

BAND_CM = (100.0, 4000.0)
MARK_CM = 2550.0
GAP_CM = (1000.0, 1200.0)
WINDOW_CM = (1900.0, 3100.0)
LAGTEST = Path("ed_data/arc/scw/lagtest_v2.json")
FLOOR_TEMPERATURES_K = (700.0, 1000.0)
DILUTE_RHO = (10.0, 30.0, 60.0)
NEAR_FLOOR_FACTOR = 2.0


def _points(cfg: Config) -> list[dict]:
    """The 13 dilute points, plus every vertex of the ten triangles."""
    from tools.corner_cells import triangle

    dilute = sorted({(t, r) for t in (700.0, 1000.0, 1500.0, 2000.0)
                     for r in DILUTE_RHO} | {(800.0, 30.0)})
    verts: set[tuple[float, float]] = set()
    for row in csv.DictReader((cfg.output_path / "S8_rcb_flip.csv").open()):
        tri = triangle(cfg, float(row["t_k"]), float(row["rho_si"]))
        if tri:
            verts.update((round(a, 6), round(b, 6)) for a, b in tri)
    out = []
    for t, r in sorted(set(dilute) | verts):
        out.append({"T_K": t, "rho_kg_m3": r,
                    "is_dilute": (t, r) in set(dilute),
                    "is_triangle_vertex": (round(t, 6), round(r, 6)) in verts})
    return out


def _x(wn, t_k):
    from constants import C_LIGHT, H_PLANCK, K_B
    return H_PLANCK * (C_LIGHT * 100.0) * np.asarray(wn, float) / (K_B * t_k)


def qcf(name: str, wn, t_k):
    """Quantum correction factor on the lineshape."""
    x = _x(wn, t_k)
    if name == "classical":
        return np.ones_like(x)
    if name == "standard":
        return 2.0 / (1.0 + np.exp(-x))
    if name == "harmonic":
        return x / -np.expm1(-x)
    if name == "schofield":
        return np.exp(0.5 * x)
    raise ValueError(name)


def _spectral(cfg: Config, source: str, t_k: float, rho: float) -> dict | None:
    from opacity.means import planck_mean, rosseland_mean
    from tools.text_checks import _load

    try:
        wn, kap = _load(cfg, source, t_k, rho)
    except FileNotFoundError:
        return None
    order = np.argsort(wn)
    wn, kap = np.asarray(wn, float)[order], np.asarray(kap, float)[order]
    m = (wn >= BAND_CM[0]) & (wn <= BAND_CM[1])
    g = (wn >= GAP_CM[0]) & (wn <= GAP_CM[1])
    rel = cfg.opacity.rel_floor
    extra = {}
    if source == "fp":
        path = cfg.fp_opacity_path / f"T{t_k:g}_rho{rho:g}_opacity.npz"
        with np.load(path) as z:
            extra["smooth_cm"] = float(z["smooth_cm"])
            extra["n_molecules"] = int(z["N"])
        # the minimum is of the SMOOTHED spectrum; check it is a resolved
        # feature and not the kernel's edge
        near = (wn >= 2600.0) & (wn <= 2800.0)
        extra["n_samples_2600_2800"] = int(near.sum())
    return {**extra,
        "kappa_R_cm2_g": float(rosseland_mean(wn, kap, t_k, BAND_CM, rel)),
        "kappa_P_cm2_g": float(planck_mean(wn, kap, t_k, BAND_CM, rel)),
        "kappa_at_2550_cm2_g": float(np.interp(MARK_CM, wn, kap)),
        "kappa_mean_1000_1200_cm2_g": float(np.trapezoid(kap[g], wn[g])
                                            / (wn[g][-1] - wn[g][0])),
        "kappa_min_in_band_cm2_g": float(kap[m].min()),
        "kappa_min_at_wavenumber_cm": float(wn[m][int(np.argmin(kap[m]))]),
        "kappa_max_in_band_cm2_g": float(kap[m].max()),
        "rel_floor_times_max_cm2_g": float(rel * kap[m].max()),
        "n_samples_in_band": int(m.sum()),
        "_wn": wn, "_kappa": kap,
    }


def lagtest_ratios() -> dict:
    """r^2 = (uncorrelated residual power) / (signal power), per lag-test cell.

    The residual is only partly white: its autocorrelation at the first lag is
    0.94-0.96 in the two cold cells, so 4-6% of its variance is uncorrelated and
    the rest is structure the spectrum resolves. Only the uncorrelated part
    spreads flat and behaves as a floor.
    """
    doc = json.loads(LAGTEST.read_text())
    out = {}
    for cell, v in doc.items():
        acf1 = float(v["acf_resid"][1])
        resid, signal = float(v["resid_mD"]), float(v["signal_mD"])
        out[cell] = {
            "resid_mD": resid, "signal_mD": signal,
            "acf_resid_lag1": acf1,
            "uncorrelated_fraction": 1.0 - acf1,
            "residual_variance_mD2": resid ** 2,
            "uncorrelated_variance_mD2": resid ** 2 * (1.0 - acf1),
            "signal_variance_mD2": signal ** 2,
            "r2_uncorrelated": resid ** 2 * (1.0 - acf1) / signal ** 2,
            "r2_total_residual": (resid / signal) ** 2,
        }
    return out


def floor_for(cfg: Config, t_k: float, rho: float, r2: float) -> dict | None:
    """The floor curve for one spectrum at a given r^2, production normalisation.

    Reuses ``figures.plot_s7_dilute_floor.floor_curve`` — the same construction
    A52 used — with r^2 replaced by the lag-test estimate. The floor inherits the
    w(1-exp(-x)) shape, so it is not flat in kappa.
    """
    from figures.plot_s7_dilute_floor import floor_curve, nyquist_wavenumber_cm

    path = cfg.fp_opacity_path / f"T{t_k:g}_rho{rho:g}_opacity.npz"
    if not path.exists():
        return None
    with np.load(path) as z:
        wn = np.asarray(z["wn"], float)
        kap = np.asarray(z["kappa"], float)
        dt_fs = float(z["dt_fs"])
    order = np.argsort(wn)
    wn, kap = wn[order], kap[order]
    nyq = nyquist_wavenumber_cm(dt_fs)
    floor = floor_curve(wn, kap, t_k, dt_fs, float(np.sqrt(r2)), nyq)
    m = (wn >= BAND_CM[0]) & (wn <= BAND_CM[1])
    g = (wn >= GAP_CM[0]) & (wn <= GAP_CM[1])
    w = (wn >= WINDOW_CM[0]) & (wn <= WINDOW_CM[1])
    return {
        "r2": r2, "dt_fs": dt_fs, "nyquist_cm": nyq,
        "floor_at_2550_cm2_g": float(np.interp(MARK_CM, wn, floor)),
        "floor_mean_1000_1200_cm2_g": float(np.trapezoid(floor[g], wn[g])
                                            / (wn[g][-1] - wn[g][0])),
        "floor_mean_window_cm2_g": float(np.trapezoid(floor[w], wn[w])
                                         / (wn[w][-1] - wn[w][0])),
        "floor_min_in_band_cm2_g": float(floor[m].min()),
        "floor_max_in_band_cm2_g": float(floor[m].max()),
        "wavenumber_fraction_above_2x_floor": float(
            np.trapezoid((kap[m] > NEAR_FLOOR_FACTOR * floor[m]).astype(float),
                         wn[m]) / (wn[m][-1] - wn[m][0])),
    }


def qcf_variants(cfg: Config, t_k: float, rho: float, fp: dict,
                 ll: dict | None) -> dict:
    """kappa_R under each lineshape correction, both readings of 'undo'."""
    from opacity.means import rosseland_mean

    rel = cfg.opacity.rel_floor
    wn, kap = fp["_wn"], fp["_kappa"]
    q_harm = qcf("harmonic", wn, t_k)
    out = {"as_delivered": {}, "relative_to_harmonic": {}}
    for name in ("classical", "standard", "harmonic", "schofield"):
        q = qcf(name, wn, t_k)
        out["as_delivered"][name] = float(
            rosseland_mean(wn, kap * q, t_k, BAND_CM, rel))
        out["relative_to_harmonic"][name] = float(
            rosseland_mean(wn, kap * q / q_harm, t_k, BAND_CM, rel))
    out["delivered_kappa_R_cm2_g"] = fp["kappa_R_cm2_g"]
    if ll is not None:
        for mode in ("as_delivered", "relative_to_harmonic"):
            out[f"fp_over_hitemp_{mode}"] = {
                k: v / ll["kappa_R_cm2_g"] for k, v in out[mode].items()}
        out["fp_over_hitemp_delivered"] = (fp["kappa_R_cm2_g"]
                                           / ll["kappa_R_cm2_g"])
    return out


def run(cfg: Config) -> dict:
    from opacity.build_tables import get_source
    from opacity.table import CM2_PER_G_TO_SI

    lag = lagtest_ratios()
    cold = [v["r2_uncorrelated"] for k, v in lag.items()
            if k in ("snip_T800_rho61", "snip_T1000_rho130")]
    r2_cold = float(np.mean(cold))
    uncorr_cold = float(np.mean([v["uncorrelated_fraction"] for k, v in lag.items()
                                 if k in ("snip_T800_rho61", "snip_T1000_rho130")]))
    from figures.plot_s7_dilute_floor import (
        DIPOLE_FLOOR_MDEBYE,
        DIPOLE_SIGNAL_RMS_MDEBYE,
    )
    r2_a52 = (DIPOLE_FLOOR_MDEBYE / DIPOLE_SIGNAL_RMS_MDEBYE) ** 2

    ll_table = get_source("hitemp", cfg, band=BAND_CM)
    doc = {
        "meta": __doc__.split("\n\n")[0],
        "band_cm": list(BAND_CM), "mark_cm": MARK_CM,
        "gap_cm": list(GAP_CM), "window_cm": list(WINDOW_CM),
        "rel_floor": cfg.opacity.rel_floor,
        "delivered_prefactor": ("alpha = w (1 - exp(-hc w / kT)) I_cl, from "
                                "ed_data/arc/water-scw-mace/convergence/04_spectra/"
                                "spectrum_convergence.py "
                                "line 28: detailed balance, no lineshape QCF"),
        "smoothing_cm": "read per point from each npz (smooth_cm)",
        "lagtest": {"file": str(LAGTEST), "cells": lag,
                    "r2_cold_cells_mean": r2_cold,
                    "cold_cells": ["snip_T800_rho61", "snip_T1000_rho130"],
                    "r2_a52_amplitude_ratio": r2_a52,
                    "uncorrelated_fraction_cold_mean": uncorr_cold,
                    "r2_a52_times_uncorrelated_fraction": r2_a52 * uncorr_cold,
                    "a52_over_lagtest": r2_a52 / r2_cold,
                    "amplitude_ratio_lagtest_cold": float(np.mean(
                        [v["r2_total_residual"] for k, v in lag.items()
                         if k in ("snip_T800_rho61", "snip_T1000_rho130")])),
                    "note": ("A52's 5.75/12.4 mDebye/atom and the lag test's "
                             "resid_mD/signal_mD are not stated to be the same "
                             "quantity; both are reported rather than reconciled")},
        "points": [],
    }
    for p in _points(cfg):
        t_k, rho = p["T_K"], p["rho_kg_m3"]
        fp = _spectral(cfg, "fp", t_k, rho)
        ll = _spectral(cfg, "hitemp", t_k, rho)
        rec = dict(p)
        rec["first_principles"] = {k: v for k, v in fp.items()
                                   if not k.startswith("_")}
        if ll is None:
            rec["hitemp"] = {
                "unavailable": ("no line-list spectrum at this state: T = 650 K "
                                "is off the 90-point grid's temperature set and "
                                "rho = 150 off its density set; no new runs"),
                "kappa_R_table_interpolated_cm2_g":
                    float(ll_table.kappa_R(t_k, rho) / CM2_PER_G_TO_SI)}
        else:
            rec["hitemp"] = {k: v for k, v in ll.items() if not k.startswith("_")}
        rec["floor_lagtest"] = floor_for(cfg, t_k, rho, r2_cold)
        rec["floor_a52"] = floor_for(cfg, t_k, rho, r2_a52)
        # A52's amplitude ratio with only the uncorrelated fraction applied:
        # isolates which of the two differences matters.
        rec["floor_a52_uncorrelated_only"] = floor_for(
            cfg, t_k, rho, r2_a52 * uncorr_cold)
        for tag in ("floor_lagtest", "floor_a52", "floor_a52_uncorrelated_only"):
            f = rec[tag]
            if not f:
                continue
            f["kappa_2550_over_floor"] = (fp["kappa_at_2550_cm2_g"]
                                          / f["floor_at_2550_cm2_g"])
            f["kappa_gap_over_floor"] = (fp["kappa_mean_1000_1200_cm2_g"]
                                         / f["floor_mean_1000_1200_cm2_g"])
            f["kappa_R_over_floor_window"] = (fp["kappa_R_cm2_g"]
                                              / f["floor_mean_window_cm2_g"])
            f["kappa_R_within_factor_2_of_floor"] = bool(
                f["kappa_R_over_floor_window"] < NEAR_FLOOR_FACTOR)
        rec["quantum_correction"] = qcf_variants(cfg, t_k, rho, fp, ll)
        doc["points"].append(rec)

    # (b)'s explicit 700/1000 K block, at each dilute density
    doc["floor_by_temperature"] = {}
    for t_k in FLOOR_TEMPERATURES_K:
        for rho in DILUTE_RHO:
            f = floor_for(cfg, t_k, rho, r2_cold)
            if f is None:
                continue
            fp = _spectral(cfg, "fp", t_k, rho)
            doc["floor_by_temperature"][f"{t_k:g}_{rho:g}"] = {
                "floor_window_cm2_g": f["floor_mean_window_cm2_g"],
                "floor_at_2550_cm2_g": f["floor_at_2550_cm2_g"],
                "kappa_at_2550_cm2_g": fp["kappa_at_2550_cm2_g"],
                "kappa_R_cm2_g": fp["kappa_R_cm2_g"],
                "kappa_2550_over_floor": (fp["kappa_at_2550_cm2_g"]
                                          / f["floor_at_2550_cm2_g"]),
                "kappa_R_over_floor_window": (fp["kappa_R_cm2_g"]
                                              / f["floor_mean_window_cm2_g"])}
    return doc


def _print(doc: dict) -> None:
    print(f"\n  (a) {len(doc['points'])} points — first principles | HITEMP\n")
    hdr = (f"  {'T':>5s} {'rho':>6s} {'kR':>9s} {'k2550':>9s} {'k(gap)':>9s} "
           f"{'min k':>9s} | {'kR':>9s} {'k2550':>9s} {'k(gap)':>9s} {'min k':>9s}")
    print(hdr); print("  " + "-" * (len(hdr) - 2))
    for p in doc["points"]:
        f, h = p["first_principles"], p["hitemp"]
        if "unavailable" in h:
            hh = (f"{h['kappa_R_table_interpolated_cm2_g']:9.4g} "
                  f"{'(table)':>9s} {'-':>9s} {'-':>9s}")
        else:
            hh = (f"{h['kappa_R_cm2_g']:9.4g} {h['kappa_at_2550_cm2_g']:9.4g} "
                  f"{h['kappa_mean_1000_1200_cm2_g']:9.4g} "
                  f"{h['kappa_min_in_band_cm2_g']:9.4g}")
        print(f"  {p['T_K']:5g} {p['rho_kg_m3']:6g} {f['kappa_R_cm2_g']:9.4g} "
              f"{f['kappa_at_2550_cm2_g']:9.4g} "
              f"{f['kappa_mean_1000_1200_cm2_g']:9.4g} "
              f"{f['kappa_min_in_band_cm2_g']:9.4g} | {hh}")

    lg = doc["lagtest"]
    print(f"\n  (b) uncorrelated residual from {lg['file']}\n")
    print(f"  {'cell':<20s} {'resid':>8s} {'signal':>8s} {'ACF1':>7s} "
          f"{'1-ACF1':>7s} {'r2_uncorr':>11s} {'r2_total':>9s}")
    for k, v in lg["cells"].items():
        print(f"  {k:<20s} {v['resid_mD']:8.3f} {v['signal_mD']:8.2f} "
              f"{v['acf_resid_lag1']:7.4f} {v['uncorrelated_fraction']:7.4f} "
              f"{v['r2_uncorrelated']:11.4g} {v['r2_total_residual']:9.4g}")
    print(f"\n  cold cells mean r2 = {lg['r2_cold_cells_mean']:.4g}; "
          f"A52's amplitude-ratio r2 = {lg['r2_a52_amplitude_ratio']:.4g} "
          f"({lg['a52_over_lagtest']:.0f}x larger)")
    print(f"\n  floor at the dilute densities (lag-test r2):")
    print(f"  {'T':>5s} {'rho':>5s} {'floor(window)':>14s} {'floor(2550)':>12s} "
          f"{'k(2550)':>10s} {'k/floor':>9s} {'kR/floor':>9s}")
    for key, v in doc["floor_by_temperature"].items():
        t, r = key.split("_")
        print(f"  {float(t):5g} {float(r):5g} {v['floor_window_cm2_g']:14.4g} "
              f"{v['floor_at_2550_cm2_g']:12.4g} {v['kappa_at_2550_cm2_g']:10.4g} "
              f"{v['kappa_2550_over_floor']:9.1f} "
              f"{v['kappa_R_over_floor_window']:9.1f}")

    print(f"\n  (c) kappa_R under each lineshape correction "
          f"(delivered = classical lineshape)\n")
    print(f"  {'T':>5s} {'rho':>6s} {'classical':>10s} {'standard':>9s} "
          f"{'harmonic':>9s} {'schofield':>10s} | {'FP/LL cl':>9s} "
          f"{'std':>7s} {'harm':>7s} {'schof':>7s}")
    for p in doc["points"]:
        q = p["quantum_correction"]
        d = q["as_delivered"]
        r = q.get("fp_over_hitemp_as_delivered")
        rr = (f"{r['classical']:9.2f} {r['standard']:7.2f} {r['harmonic']:7.2f} "
              f"{r['schofield']:7.2f}" if r else f"{'-':>9s} {'-':>7s} "
              f"{'-':>7s} {'-':>7s}")
        print(f"  {p['T_K']:5g} {p['rho_kg_m3']:6g} {d['classical']:10.4g} "
              f"{d['standard']:9.4g} {d['harmonic']:9.4g} {d['schofield']:10.4g} "
              f"| {rr}")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("config", type=Path)
    args = ap.parse_args(argv)
    cfg = load_config(args.config)
    doc = run(cfg)
    path = cfg.output_path / "floor_check.json"
    atomicio.write_text(path, json.dumps(doc, indent=1) + "\n")
    _print(doc)
    print(f"\n  wrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
