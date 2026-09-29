#!/usr/bin/env python3
"""Calibrate the lineshape quantum correction against the line list (A76).

    python tools/qcf_calibration.py config/fiducial.yaml

Existing data only; no new physics is run. Output: ``outputs/qcf_calibration.json``.

The delivered first-principles spectra carry alpha = w (1 - exp(-hc w / kT)) I_cl
— the detailed-balance factor with **no quantum correction on the lineshape**
(``ed_data/arc/water-scw-mace/convergence/04_spectra/
spectrum_convergence.py`` line 28; A75). Which of the
three standard corrections belongs there is not decidable from the spectra alone,
but it is decidable against the line list, because the line list is quantum by
construction and the dipole model's own band-intensity error is known
independently from the monomer test.

So: integrate kappa(nu) over each band for both sources, form FP/HITEMP per band
under each correction, and divide out the monomer test's intensity factor. What
is left should be 1, with no trend in temperature, for the correct choice — the
rotation and bend bands are the diagnostic because their intensity factors are
measured tightly (0.974; 1.25-1.32) while the stretch factor spans 2.1-3.9 and
cannot discriminate.

The correction is applied to the first-principles side only. HITEMP is a
line-by-line quantum calculation and takes none.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np

from config import Config, load_config

# The 13 points at rho <= 60, plus (800, 100).
POINTS = sorted({(t, r) for t in (700.0, 1000.0, 1500.0, 2000.0)
                 for r in (10.0, 30.0, 60.0)} | {(800.0, 30.0), (800.0, 100.0)})
BANDS_CM = {"rotation": (100.0, 1000.0),
            "bend": (1400.0, 1800.0),
            "stretch": (3200.0, 4000.0)}
CORRECTIONS = ("delivered", "standard", "harmonic", "schofield")

# Band intensity factor of the dipole model, from the monomer test (supplied).
# The bend pair is the square of the measured derivative ratios 1.12 and 1.15
# (intensity goes as |dmu/dQ|^2), which is the one internal check available on
# these numbers; the rotation and stretch values are taken as given.
INTENSITY_FACTOR = {"rotation": (0.974, 0.974),
                    "bend": (1.25, 1.32),
                    "stretch": (2.1, 3.9)}
INTENSITY_SOURCE = {
    "rotation": "monomer test, supplied: 0.974",
    "bend": "monomer test, supplied: 1.25-1.32 (= 1.12^2 and 1.15^2, the "
            "measured bend derivative ratios at r_OH = 0.957 and 1.00 A)",
    "stretch": "monomer test, supplied: 2.1 at equilibrium to 3.9 at 1.00 A",
}
# Which bands can discriminate: a factor known to a few per cent can, one
# spanning a factor of 1.9 cannot.
DIAGNOSTIC_BANDS = ("rotation", "bend")


def qcf(name: str, wn, t_k):
    """Lineshape quantum correction factor. x = hc nu / kT."""
    from constants import C_LIGHT, H_PLANCK, K_B

    x = H_PLANCK * (C_LIGHT * 100.0) * np.asarray(wn, float) / (K_B * t_k)
    if name == "delivered":       # Q = 1: the classical lineshape, as shipped
        return np.ones_like(x)
    if name == "standard":
        return 2.0 / (1.0 + np.exp(-x))
    if name == "harmonic":
        return x / -np.expm1(-x)
    if name == "schofield":
        return np.exp(0.5 * x)
    raise ValueError(name)


def _spectrum(cfg: Config, source: str, t_k: float, rho: float):
    from tools.text_checks import _load

    wn, kap = _load(cfg, source, t_k, rho)
    order = np.argsort(wn)
    return np.asarray(wn, float)[order], np.asarray(kap, float)[order]


def _band_integral(wn, kap, lo, hi) -> float:
    """int kappa dnu over [lo, hi], with the edges interpolated so the limits do
    not depend on either source's sampling."""
    m = (wn > lo) & (wn < hi)
    w = np.concatenate([[lo], wn[m], [hi]])
    k = np.concatenate([[float(np.interp(lo, wn, kap))], kap[m],
                        [float(np.interp(hi, wn, kap))]])
    return float(np.trapezoid(k, w))


def _trend(temps, values, densities=None) -> dict:
    """Trend of a corrected ratio against T, in both linear and log space.

    Ratios are multiplicative, so d ln(ratio)/dT is the honest measure; the
    linear slope is kept because it is what the eye reads off the table. The
    fixed-density sub-series separates a temperature trend from a density trend,
    since the 14 points span rho = 10-100 as well as T = 700-2000 K.
    """
    t = np.asarray(temps, float)
    v = np.asarray(values, float)
    ok = np.isfinite(v) & (v > 0)
    t, v = t[ok], v[ok]
    if t.size < 3 or np.ptp(t) == 0:
        return {"n": int(t.size)}
    slope, intercept = np.polyfit(t, v, 1)
    resid = v - (slope * t + intercept)
    ss_tot = float(((v - v.mean()) ** 2).sum())
    log_slope = float(np.polyfit(t, np.log(v), 1)[0])
    rt, rv = t.argsort().argsort().astype(float), v.argsort().argsort().astype(float)
    out = {"n": int(t.size),
           "slope_per_1000K": float(slope * 1000.0),
           "slope_per_1000K_as_fraction_of_median": float(
               slope * 1000.0 / np.median(v)),
           "dlnratio_dT_per_1000K": log_slope * 1000.0,
           "ratio_hottest_over_coldest": float(v[t.argmax()] / v[t.argmin()]),
           "spearman_rho": float(np.corrcoef(rt, rv)[0, 1]),
           "r_squared": float(1.0 - (resid ** 2).sum() / ss_tot) if ss_tot else None}
    if densities is not None:
        d = np.asarray(densities, float)[ok]
        # same trend using only rho = 30, which spans the full temperature range
        m = d == 30.0
        if m.sum() >= 3 and np.ptp(t[m]) > 0:
            out["at_fixed_rho_30"] = {
                "n": int(m.sum()),
                "dlnratio_dT_per_1000K": float(
                    np.polyfit(t[m], np.log(v[m]), 1)[0] * 1000.0),
                "ratio_hottest_over_coldest": float(
                    v[m][t[m].argmax()] / v[m][t[m].argmin()])}
        # and the density trend at fixed T, pooled over the four isotherms
        slopes = []
        for t0 in np.unique(t):
            mm = (t == t0) & (d > 0)
            if mm.sum() >= 3 and np.ptp(d[mm]) > 0:
                slopes.append(float(np.polyfit(np.log(d[mm]), np.log(v[mm]), 1)[0]))
        if slopes:
            out["dlnratio_dlnrho_at_fixed_T"] = {
                "median": float(np.median(slopes)),
                "min": float(np.min(slopes)), "max": float(np.max(slopes)),
                "n_isotherms": len(slopes)}
    return out


def run(cfg: Config) -> dict:
    doc = {
        "meta": __doc__.split("\n\n")[0],
        "integrand": "int kappa(nu) dnu over each band, edges interpolated",
        "bands_cm": {k: list(v) for k, v in BANDS_CM.items()},
        "corrections": list(CORRECTIONS),
        "correction_applied_to": ("the first-principles side only; HITEMP is a "
                                  "line-by-line quantum calculation"),
        "delivered_form": ("alpha = w (1 - exp(-hc w / kT)) I_cl, i.e. Q = 1 on "
                           "the lineshape (A75)"),
        "intensity_factor": {k: list(v) for k, v in INTENSITY_FACTOR.items()},
        "intensity_factor_source": INTENSITY_SOURCE,
        "diagnostic_bands": list(DIAGNOSTIC_BANDS),
        "points": [],
    }
    for t_k, rho in POINTS:
        wn_f, k_f = _spectrum(cfg, "fp", t_k, rho)
        wn_l, k_l = _spectrum(cfg, "hitemp", t_k, rho)
        rec = {"T_K": t_k, "rho_kg_m3": rho, "bands": {}}
        for band, (lo, hi) in BANDS_CM.items():
            ll = _band_integral(wn_l, k_l, lo, hi)
            entry = {"hitemp_integral_cm_g": ll, "first_principles": {},
                     "ratio": {}, "corrected_ratio": {}}
            for name in CORRECTIONS:
                fp = _band_integral(wn_f, k_f * qcf(name, wn_f, t_k), lo, hi)
                entry["first_principles"][name] = fp
                entry["ratio"][name] = fp / ll
                f_lo, f_hi = INTENSITY_FACTOR[band]
                entry["corrected_ratio"][name] = {
                    "at_factor_low": fp / ll / f_lo,
                    "at_factor_high": fp / ll / f_hi,
                }
            rec["bands"][band] = entry
        doc["points"].append(rec)

    # summary: median, range and temperature trend of the corrected ratio
    doc["summary"] = {}
    temps = [p["T_K"] for p in doc["points"]]
    rhos = [p["rho_kg_m3"] for p in doc["points"]]
    for band in BANDS_CM:
        doc["summary"][band] = {}
        for name in CORRECTIONS:
            block = {}
            for end in ("at_factor_low", "at_factor_high"):
                vals = [p["bands"][band]["corrected_ratio"][name][end]
                        for p in doc["points"]]
                block[end] = {
                    "median": float(np.median(vals)),
                    "min": float(np.min(vals)), "max": float(np.max(vals)),
                    "log10_spread_dex": float(np.log10(np.max(vals)
                                                       / np.min(vals))),
                    "distance_of_median_from_1": float(abs(np.median(vals) - 1.0)),
                    "trend_with_T": _trend(temps, vals, rhos),
                }
            doc["summary"][band][name] = block

    # The band-strength sum rule. Q_harmonic x (1 - exp(-x)) = x, so the
    # harmonic-corrected form is alpha ~ beta w^2 I_cl -- the plain classical
    # absorption formula. For a classical trajectory
    #     int beta w^2 I_cl dw ~ beta <|mu_dot|^2> = beta (kT/m) sum (dmu/dQ)^2
    # is independent of temperature, and the same sum rule holds for the line
    # list's net (absorption minus stimulated emission) band strength. So
    # int kappa_harmonic dnu over a window wide enough to hold the whole band
    # must be flat in T on both sides. Any slope is upstream of the lineshape
    # and cannot be repaired by any choice of quantum correction factor.
    doc["sum_rule"] = {
        "statement": ("int kappa_harmonic dnu over a wide window is "
                      "T-independent for a classical trajectory (index p = 0)"),
        "windows_cm": {"rotation": [50.0, 1100.0], "bend": [1000.0, 2400.0],
                       "stretch": [2800.0, 4000.0], "full": [100.0, 4000.0]},
        "index_p_of_T": {},
    }
    for tag, (lo, hi) in doc["sum_rule"]["windows_cm"].items():
        entry = {"first_principles_by_rho": {}, "hitemp_at_rho_30": None}
        for rho0 in (10.0, 30.0, 60.0):
            ts = [t for t, r in POINTS if r == rho0]
            vals = []
            for t in ts:
                wf, kf = _spectrum(cfg, "fp", t, rho0)
                vals.append(_band_integral(wf, kf * qcf("harmonic", wf, t), lo, hi))
            entry["first_principles_by_rho"][f"{rho0:g}"] = float(
                np.polyfit(np.log(ts), np.log(vals), 1)[0])
        ts = [t for t, r in POINTS if r == 30.0]
        vals = []
        for t in ts:
            wl, kl = _spectrum(cfg, "hitemp", t, 30.0)
            vals.append(_band_integral(wl, kl, lo, hi))
        entry["hitemp_at_rho_30"] = float(np.polyfit(np.log(ts), np.log(vals), 1)[0])
        doc["sum_rule"]["index_p_of_T"][tag] = entry

    # Robustness: is the bend trend a window artefact? The first-principles bend
    # band sits 5-26 cm^-1 above the line list's, so a fixed 1400-1800 window
    # could in principle catch a drifting fraction of it. Widening the window
    # settles it.
    doc["bend_window_robustness"] = {
        "windows_cm": [[1400.0, 1800.0], [1200.0, 2100.0], [1000.0, 2400.0]],
        "note": ("delivered-form ratio in each window, with the band centroid "
                 "and the fraction of the 1000-2400 integral inside 1400-1800"),
        "points": [],
    }
    for t_k, rho in POINTS:
        wf, kf = _spectrum(cfg, "fp", t_k, rho)
        wl, kl = _spectrum(cfg, "hitemp", t_k, rho)

        def centroid(w, k, lo=1400.0, hi=1800.0):
            m = (w >= lo) & (w <= hi)
            return float(np.trapezoid(w[m] * k[m], w[m]) / np.trapezoid(k[m], w[m]))

        doc["bend_window_robustness"]["points"].append({
            "T_K": t_k, "rho_kg_m3": rho,
            "ratio_by_window": [_band_integral(wf, kf, lo, hi)
                                / _band_integral(wl, kl, lo, hi)
                                for lo, hi in
                                doc["bend_window_robustness"]["windows_cm"]],
            "centroid_cm": {"first_principles": centroid(wf, kf),
                            "hitemp": centroid(wl, kl)},
            "fraction_of_1000_2400_inside_window": {
                "first_principles": _band_integral(wf, kf, 1400, 1800)
                                    / _band_integral(wf, kf, 1000, 2400),
                "hitemp": _band_integral(wl, kl, 1400, 1800)
                          / _band_integral(wl, kl, 1000, 2400)},
        })
    temps_r = [p["T_K"] for p in doc["bend_window_robustness"]["points"]]
    doc["bend_window_robustness"]["trend_by_window"] = [
        _trend(temps_r, [p["ratio_by_window"][i]
                         for p in doc["bend_window_robustness"]["points"]])
        for i in range(3)]

    # Neither side carries MT_CKD here: the loader returns the stored kappa, so
    # the line list is lines only while the first-principles spectrum, being bulk
    # MD, contains the continuum intrinsically. That asymmetry flatters the
    # first-principles side, so adding the continuum to the denominator can only
    # deepen a deficit. It matters in the rotation band at cold, dense states.
    doc["continuum_asymmetry"] = {
        "note": ("MT_CKD self-continuum added to the HITEMP rotation integral; "
                 "the first-principles spectrum already contains it"),
        "points": [],
    }
    from opacity.mtckd import load_reference, self_continuum_kappa
    ref = load_reference(cfg.mtckd_ref_path)
    lo, hi = BANDS_CM["rotation"]
    for t_k, rho in POINTS:
        wl, kl = _spectrum(cfg, "hitemp", t_k, rho)
        wf, kf = _spectrum(cfg, "fp", t_k, rho)
        ll = _band_integral(wl, kl, lo, hi)
        lc = _band_integral(wl, kl + self_continuum_kappa(wl, t_k, rho, ref), lo, hi)
        fp = _band_integral(wf, kf, lo, hi)
        doc["continuum_asymmetry"]["points"].append({
            "T_K": t_k, "rho_kg_m3": rho,
            "continuum_share_of_line_integral": (lc - ll) / ll,
            "ratio_lines_only": fp / ll,
            "ratio_with_continuum": fp / lc})

    # the verdict, over the diagnostic bands only
    doc["verdict"] = {}
    for band in DIAGNOSTIC_BANDS:
        rows = []
        for name in CORRECTIONS:
            b = doc["summary"][band][name]["at_factor_low"]
            b2 = doc["summary"][band][name]["at_factor_high"]
            # score at the midpoint of the factor range
            med = 0.5 * (b["median"] + b2["median"])
            tr = b["trend_with_T"].get("dlnratio_dT_per_1000K", 0.0)
            rows.append({"correction": name, "median_corrected_ratio": med,
                         "distance_from_1_in_dex": abs(float(np.log10(med))),
                         "dlnratio_dT_per_1000K": tr,
                         "ratio_hottest_over_coldest":
                             b["trend_with_T"].get("ratio_hottest_over_coldest"),
                         "spearman_rho": b["trend_with_T"].get("spearman_rho")})
        doc["verdict"][band] = {
            "rows": rows,
            "closest_median_to_1":
                min(rows, key=lambda r: r["distance_from_1_in_dex"])["correction"],
            # a "closest" is only a finding if something is actually close
            "any_within_0p15_dex_of_1": [
                r["correction"] for r in rows
                if r["distance_from_1_in_dex"] <= 0.15],
            "all_corrections_leave_a_trend": all(
                abs(r["dlnratio_dT_per_1000K"]) > 0.05 for r in rows),
            "flattest_in_T":
                min(rows, key=lambda r: abs(r["dlnratio_dT_per_1000K"]))["correction"],
        }
    return doc


def _print(doc: dict) -> None:
    print("\n  FP/HITEMP band ratio, int kappa dnu, per correction\n")
    for band in BANDS_CM:
        print(f"  --- {band} {doc['bands_cm'][band]} cm^-1 "
              f"(intensity factor {doc['intensity_factor'][band]}) ---")
        print(f"  {'T':>5s} {'rho':>5s} " +
              " ".join(f"{c:>10s}" for c in CORRECTIONS))
        for p in doc["points"]:
            r = p["bands"][band]["ratio"]
            print(f"  {p['T_K']:5g} {p['rho_kg_m3']:5g} " +
                  " ".join(f"{r[c]:10.3f}" for c in CORRECTIONS))
        print()

    print("\n  Corrected ratio = (FP/HITEMP) / intensity factor\n")
    for band in BANDS_CM:
        print(f"  --- {band} ---")
        print(f"  {'correction':<15s} {'median':>9s} {'min':>9s} {'max':>9s} "
              f"{'dex':>6s} {'dlnr/dT/kK':>11s} {'hot/cold':>9s} {'rho_s':>6s}")
        for name in CORRECTIONS:
            for end, tag in (("at_factor_low", "lo"), ("at_factor_high", "hi")):
                b = doc["summary"][band][name][end]
                t = b["trend_with_T"]
                same = (doc["intensity_factor"][band][0]
                        == doc["intensity_factor"][band][1])
                if same and tag == "hi":
                    continue
                label = name if same else f"{name} ({tag})"
                print(f"  {label:<15s} {b['median']:9.3f} {b['min']:9.3f} "
                      f"{b['max']:9.3f} {b['log10_spread_dex']:6.2f} "
                      f"{t.get('dlnratio_dT_per_1000K', float('nan')):+11.3f} "
                      f"{t.get('ratio_hottest_over_coldest', float('nan')):9.2f} "
                      f"{t.get('spearman_rho', float('nan')):6.2f}")
        print()

    print("  Band-strength sum rule: index p of int kappa_harmonic dnu ~ T^p")
    print("  (classical trajectory and line list both require p = 0)")
    print(f"  {'window':<10s}{'FP rho=10':>10s}{'rho=30':>9s}{'rho=60':>9s}{'HITEMP':>9s}")
    for tag, e in doc["sum_rule"]["index_p_of_T"].items():
        b = e["first_principles_by_rho"]
        print(f"  {tag:<10s}{b['10']:10.2f}{b['30']:9.2f}{b['60']:9.2f}"
              f"{e['hitemp_at_rho_30']:9.2f}")
    print()

    r = doc["bend_window_robustness"]
    print("  Bend trend vs window (delivered form): dln(ratio)/dT per 1000 K")
    for w, t in zip(r["windows_cm"], r["trend_by_window"]):
        print(f"    {w[0]:.0f}-{w[1]:.0f} cm^-1: {t['dlnratio_dT_per_1000K']:+.3f}"
              f"   hot/cold {t['ratio_hottest_over_coldest']:.2f}")
    print()

    for band, v in doc["verdict"].items():
        near = v["any_within_0p15_dex_of_1"]
        head = (f"median closest to 1 = {v['closest_median_to_1']}" if near
                else "no correction brings the median within 0.15 dex of 1 "
                     "(only flatness discriminates here)")
        print(f"  VERDICT ({band}): {head}; flattest in T = {v['flattest_in_T']}")
        for r in v["rows"]:
            print(f"    {r['correction']:<10s} median {r['median_corrected_ratio']:7.3f} "
                  f"({r['distance_from_1_in_dex']:.3f} dex from 1)  "
                  f"dln/dT {r['dlnratio_dT_per_1000K']:+.3f} /1000K  "
                  f"hot/cold {r['ratio_hottest_over_coldest']:.2f}  "
                  f"rho_s {r['spearman_rho']:+.2f}")
        print()


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("config", type=Path)
    args = ap.parse_args(argv)
    cfg = load_config(args.config)
    doc = run(cfg)
    path = cfg.output_path / "qcf_calibration.json"
    path.write_text(json.dumps(doc, indent=1) + "\n")
    _print(doc)
    print(f"  wrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
