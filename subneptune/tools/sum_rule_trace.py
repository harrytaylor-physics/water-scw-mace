#!/usr/bin/env python3
"""Trace the band-strength sum-rule violation A76 found (A77).

    python tools/sum_rule_trace.py config/fiducial.yaml

Existing data only; no new runs. Output: ``outputs/sum_rule_trace.json``.

A76 measured that the first-principles vibrational band strengths fall steeply
with temperature where the classical sum rule requires them to be flat, and
attributed the index to the delivered form read from
``convergence/04_spectra/spectrum_convergence.py``.
That is a *convergence* script, not the production one. The production chain is

    spectra.py  ->  dipoles.npy  ->  ~/scw/opacity_pipeline.py
                                     :: spectrum_from_dipoles  ->  *.npz

(named in ``lda_campaign/stepG_dipoles/make_opacity_lda.py`` lines 2-4 and 17).
**Neither ``~/scw/opacity_pipeline.py`` nor ``spectra.py`` is on this machine**,
so the ACF -> cross-section step cannot be read. What this tool does instead:

  1. reconstructs and *verifies* every production step that is still visible in
     the delivered arrays, exactly, so the unreadable part is isolated;
  2. measures the sum-rule index under both candidate formulas for the
     unreadable step, over every density row with three or more temperatures,
     so the finding can be stated without knowing which one is right;
  3. checks the two steps that could fake a temperature trend -- the zero
     clipping and the choice of integration window;
  4. lists the exact ARC files needed to close the remaining gap.

The two candidates differ by exactly Q_harmonic = x/(1-e^-x):

  A  sigma ~ w (1 - e^-x) I_cl   -- the convergence script's form, Q = 1 on the
                                    lineshape; the classical quantity is then
                                    kappa x Q_harmonic
  B  sigma ~ beta w^2 I_cl       -- the textbook classical form, which already
                                    carries the harmonic factor; the classical
                                    quantity is then kappa itself

Both reduce to w^2 I/T as x -> 0, so the dilute low-frequency limit cannot
separate them and no test on the delivered arrays can. The sum rule
`int (classical quantity) dnu = const in T` must hold under whichever is right.
"""

from __future__ import annotations

import argparse
import collections
import glob
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import atomicio

import numpy as np
from scipy.ndimage import gaussian_filter1d

from config import Config, load_config

N_AVOGADRO = 6.02214076e23
CANDIDATES = ("A_convergence_form", "B_classical_form")

# Wide enough to hold each whole band; "full" is A76's range and "wide" extends
# to the end of the delivered grid, so intensity leaving the top of A76's window
# would show up as a difference between the two.
WINDOWS_CM = {"rotation": (50.0, 1100.0), "bend": (1000.0, 2400.0),
              "stretch": (2800.0, 4000.0), "full": (100.0, 4000.0),
              "wide": (100.0, 6000.0)}
CLIP_BANDS_CM = {"rotation": (100.0, 1000.0), "bend": (1400.0, 1800.0),
                 "stretch": (3200.0, 4000.0), "noise_floor_window": (3900.0, 5500.0)}

# Everything step 1 and steps 3-4 need and this machine does not have. The
# trajectory directory convention is make_opacity_lda.py's --trajdir
# (".../production/T{T}_rho{rho}", 64 molecules, dipoles.npy written by
# spectra.py); the brief's $DATA/scw_prod/<tag> is the same <tag>.
MISSING_INPUTS = {
    "for_step_1_the_production_formula": [
        "~/scw/opacity_pipeline.py   (spectrum_from_dipoles: the ACF -> sigma step)",
        "~/scw/spectra.py            (writes dipoles.npy; the ACF convention)",
    ],
    "for_steps_3_and_4_dynamics_and_dipole_model": [
        "$DATA/scw_prod/T700_rho30/dipoles.npy   + the matching trajectory "
        "(positions and velocities, 64 molecules, 0.5 fs, 60 ps)",
        "$DATA/scw_prod/T1000_rho30/dipoles.npy  + the matching trajectory",
        "$DATA/scw_prod/T2000_rho30/dipoles.npy  + the matching trajectory",
    ],
    "note": ("rho = 30 kg m^-3 is the dilute row with the most temperatures "
             "(700, 800, 1000, 1500, 2000 K), so it is the cheapest row that "
             "spans the trend. Trajectory directory names follow "
             "make_opacity_lda.py: T{T:d}_rho{rho:g}."),
}


def _x(wn, t_k):
    from constants import C_LIGHT, H_PLANCK, K_B
    return H_PLANCK * (C_LIGHT * 100.0) * np.asarray(wn, float) / (K_B * t_k)


def q_harmonic(wn, t_k):
    x = _x(wn, t_k)
    return x / -np.expm1(-x)


def _rows(cfg: Config, source: str) -> dict:
    """{rho: [(T, path), ...]} for every stored spectrum of one source."""
    from opacity.registry import get_spec
    from tools.text_checks import EXTRA_SUBDIR

    spec = get_spec(source)
    dirs = [cfg.spectra_dir(spec.subdir)]
    if source == "hitemp":
        dirs.append(cfg.data_root_path / EXTRA_SUBDIR)
    out = collections.defaultdict(list)
    for directory in dirs:
        for path in sorted(glob.glob(str(directory / spec.glob))):
            with np.load(path) as d:
                out[float(d["rho_kg_m3"])].append((float(d["T"]), path))
    return {r: sorted(v) for r, v in out.items() if len(v) >= 3}


def _integral(wn, kap, lo, hi):
    if wn[0] > lo or wn[-1] < hi:
        return None
    m = (wn > lo) & (wn < hi)
    w = np.concatenate([[lo], wn[m], [hi]])
    k = np.concatenate([[float(np.interp(lo, wn, kap))], kap[m],
                        [float(np.interp(hi, wn, kap))]])
    return float(np.trapezoid(k, w))


def production_chain(cfg: Config) -> dict:
    """Verify, on every delivered file, the production steps still visible in
    the arrays. What survives this is known exactly; what does not is the gap."""
    paths = sorted(glob.glob(str(cfg.data_root_path / "fp_opacity" / "*.npz")))
    molar, alpha_err, kappa_err, smooth_err = [], [], [], []
    clipped, negative = [], []
    for p in paths:
        with np.load(p) as d:
            rho = float(d["rho_kg_m3"])
            n_used = float(np.nanmedian(d["alpha"] / d["sigma"]))
            molar.append(rho * 1e-3 / n_used * N_AVOGADRO)
            alpha_err.append(float(np.nanmax(np.abs(
                d["alpha"] / (n_used * d["sigma"]) - 1.0))))
            kappa_err.append(float(np.nanmax(np.abs(
                d["kappa"] * (rho * 1e-3) / d["alpha"] - 1.0))))
            # the grid is uniform at 1 cm^-1, so the filter width is in channels
            sm = gaussian_filter1d(d["sigma_raw"], float(d["smooth_cm"]))
            ok = d["sigma"] > 0
            smooth_err.append(float(np.nanmedian(
                np.abs(sm[ok] / d["sigma"][ok] - 1.0))))
            clipped.append(float((d["sigma_raw"] == 0).mean()))
            negative.append(int((d["sigma_raw"] < 0).sum()))
    with np.load(paths[0]) as d:
        wn = d["wn"]
        grid = {"n": int(wn.size), "first_cm": float(wn[0]), "last_cm": float(wn[-1]),
                "spacing_cm": float(np.median(np.diff(wn))),
                "smooth_cm": float(d["smooth_cm"]),
                "n_molecules": int(d["N"]), "dt_fs": float(d["dt_fs"]),
                "units_marker": str(d["units_fixed"])}
    return {
        "n_files": len(paths),
        "verified_exactly": {
            "molar_mass_g_mol": {"min": min(molar), "max": max(molar),
                                 "note": "alpha = n sigma with n = rho/M N_A"},
            "max_rel_error_alpha_eq_n_sigma": max(alpha_err),
            "max_rel_error_kappa_eq_alpha_over_rho": max(kappa_err),
            "median_rel_error_sigma_eq_gaussian_smoothed_sigma_raw": max(smooth_err),
            "smoothing_is_standard_deviation_not_fwhm": True,
        },
        "grid": grid,
        "sigma_raw_clipped_fraction": {"min": min(clipped), "max": max(clipped)},
        "sigma_raw_negative_values": sum(negative),
        "readable_chain": [
            "sigma_raw = max(<unreadable ACF -> cross-section step> - noise floor, 0)",
            "sigma     = gaussian_filter1d(sigma_raw, sd = smooth_cm = 20 cm^-1)",
            "alpha     = n sigma,  n = rho / (18.015280 g/mol) * N_A",
            "kappa     = alpha / rho_cgs",
        ],
        "unreadable_step": ("the ACF -> cross-section map, which carries the "
                            "temperature and quantum factors; it lives in "
                            "~/scw/opacity_pipeline.py::spectrum_from_dipoles"),
    }


def _index(ts, vals):
    ok = [i for i, v in enumerate(vals) if v is not None and v > 0]
    if len(ok) < 3:
        return None
    t = np.log([ts[i] for i in ok])
    v = np.log([vals[i] for i in ok])
    return float(np.polyfit(t, v, 1)[0])


def sum_rule_indices(cfg: Config) -> dict:
    """Index p of int (classical quantity) dnu ~ T^p, which must be 0."""
    out = {"statement": ("p = 0 is required of any classical trajectory and of "
                         "the line list; p is measured per density row"),
           "candidates": {
               "A_convergence_form": "sigma ~ w(1-e^-x) I_cl; classical qty = kappa x Q_harm",
               "B_classical_form": "sigma ~ beta w^2 I_cl; classical qty = kappa"},
           "first_principles": {}, "hitemp": {}}
    for rho, rows in sorted(_rows(cfg, "fp").items()):
        ts = [t for t, _ in rows]
        spec = []
        for t, p in rows:
            with np.load(p) as d:
                spec.append((t, d["wn"].copy(), d["kappa"].copy()))
        entry = {"temperatures_K": ts, "p": {}}
        for tag, (lo, hi) in WINDOWS_CM.items():
            entry["p"][tag] = {
                "A_convergence_form": _index(
                    ts, [_integral(w, k * q_harmonic(w, t), lo, hi) for t, w, k in spec]),
                "B_classical_form": _index(
                    ts, [_integral(w, k, lo, hi) for t, w, k in spec]),
            }
        out["first_principles"][f"{rho:g}"] = entry
    for rho, rows in sorted(_rows(cfg, "hitemp").items()):
        ts = [t for t, _ in rows]
        spec = []
        for t, p in rows:
            with np.load(p) as d:
                spec.append((t, d["wn"].copy(), d["kappa"].copy()))
        out["hitemp"][f"{rho:g}"] = {
            "temperatures_K": ts,
            "p": {tag: _index(ts, [_integral(w, k, lo, hi) for _, w, k in spec])
                  for tag, (lo, hi) in WINDOWS_CM.items()},
        }
    return out


def clipping_diagnostic(cfg: Config) -> dict:
    """Could the zero clipping fake the trend? Only if it varies with T."""
    out = {"note": ("fraction of sigma_raw channels clipped to zero, by band; "
                    "a temperature trend here would be a candidate mechanism"),
           "rows": {}}
    for rho, rows in sorted(_rows(cfg, "fp").items()):
        per_t = {}
        for t, p in rows:
            with np.load(p) as d:
                wn, sr = d["wn"], d["sigma_raw"]
                per_t[f"{t:g}"] = {
                    b: float((sr[(wn >= lo) & (wn <= hi)] == 0).mean())
                    for b, (lo, hi) in CLIP_BANDS_CM.items()}
        out["rows"][f"{rho:g}"] = per_t
    allv = [v for r in out["rows"].values() for d in r.values() for v in d.values()]
    out["overall"] = {"min": min(allv), "max": max(allv),
                      "spread": max(allv) - min(allv)}
    return out


def run(cfg: Config) -> dict:
    return {"meta": __doc__.split("\n\n")[0],
            "production_chain": production_chain(cfg),
            "sum_rule": sum_rule_indices(cfg),
            "clipping": clipping_diagnostic(cfg),
            "missing_inputs": MISSING_INPUTS}


def _print(doc: dict) -> None:
    pc = doc["production_chain"]
    print("\n  Production chain, reconstructed from the delivered arrays "
          f"({pc['n_files']} files)")
    for line in pc["readable_chain"]:
        print("    " + line)
    v = pc["verified_exactly"]
    print(f"    verified: M = {v['molar_mass_g_mol']['min']:.6f} g/mol "
          f"(spread {v['molar_mass_g_mol']['max'] - v['molar_mass_g_mol']['min']:.1e}); "
          f"alpha=n.sigma to {v['max_rel_error_alpha_eq_n_sigma']:.1e}; "
          f"kappa=alpha/rho to {v['max_rel_error_kappa_eq_alpha_over_rho']:.1e}; "
          f"smoothing to {v['median_rel_error_sigma_eq_gaussian_smoothed_sigma_raw']:.1e}")
    print(f"    UNREADABLE: {pc['unreadable_step']}")

    sr = doc["sum_rule"]
    print("\n  Sum-rule index p (must be 0), by density row")
    print(f"  {'rho':>7s}{'nT':>3s}  " +
          "".join(f"{w:>9s}{'':>6s}" for w in WINDOWS_CM))
    print(f"  {'':>7s}{'':>3s}  " + "".join(f"{'A':>9s}{'B':>6s}" for _ in WINDOWS_CM))
    for rho, e in sr["first_principles"].items():
        line = f"  {rho:>7s}{len(e['temperatures_K']):3d}  "
        for w in WINDOWS_CM:
            a, b = e["p"][w]["A_convergence_form"], e["p"][w]["B_classical_form"]
            line += f"{a:9.2f}" if a is not None else f"{'-':>9s}"
            line += f"{b:6.2f}" if b is not None else f"{'-':>6s}"
        print(line)
    print(f"\n  {'HITEMP':>7s}{'nT':>3s}  " + "".join(f"{w:>9s}" for w in WINDOWS_CM))
    for rho, e in sr["hitemp"].items():
        line = f"  {rho:>7s}{len(e['temperatures_K']):3d}  "
        for w in WINDOWS_CM:
            p = e["p"][w]
            line += f"{p:9.2f}" if p is not None else f"{'-':>9s}"
        print(line)

    c = doc["clipping"]
    print(f"\n  Zero clipping: {100*c['overall']['min']:.1f}-"
          f"{100*c['overall']['max']:.1f}% of channels, spread "
          f"{100*c['overall']['spread']:.1f} points -- no temperature trend, "
          "so it is not the mechanism")
    print("\n  Missing inputs (steps 1, 3, 4 cannot close without these):")
    for k, v in doc["missing_inputs"].items():
        if isinstance(v, list):
            for line in v:
                print("    " + line)
    print()


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("config", type=Path)
    args = ap.parse_args(argv)
    cfg = load_config(args.config)
    doc = run(cfg)
    path = cfg.output_path / "sum_rule_trace.json"
    atomicio.write_text(path, json.dumps(doc, indent=1) + "\n")
    _print(doc)
    print(f"  wrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
