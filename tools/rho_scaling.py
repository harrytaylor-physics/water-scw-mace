#!/usr/bin/env python3
"""Does the first-principles absorption grow faster than rho^2? (Task 38)

    python tools/rho_scaling.py config/fiducial.yaml

Every continuum model the literature offers — MT_CKD, far-wing theory, the
laboratory self-continuum coefficients — is a **pair** quantity: two molecules
interacting, so absorption per unit volume goes as the square of the number
density. Write alpha = rho * kappa for the absorption coefficient; then

    alpha / rho^2  =  kappa / rho

is flat in rho for anything that is exactly pair-wise, rises where three-body and
higher clusters contribute, and falls towards 1/rho where isolated lines dominate
(alpha proportional to rho, the dilute limit).

This computes that quantity along every isotherm carrying at least three
densities, at three places in the spectrum:

  ``window_2550``      kappa at 2550 cm^-1, the transparency window the
                       manuscript quotes (same wavenumber as tools/text_checks.py);
  ``window_1000_1200`` mean kappa over 1000-1200 cm^-1, the gap between the
                       rotation and bend bands;
  ``rosseland``        the Rosseland mean on the common band, which is what the
                       structure model actually uses.

and the local power-law index d ln alpha / d ln rho between adjacent densities.

**The line-list sources are run on their own grid, not the first-principles one.**
The 90-point line-list grid is a clean T x rho product (rho = 0.1 ... 1000) and
shares six isotherms with the 48 computed points, but only 20 of the 48 states
have a line-list spectrum at the same (T, rho). Interpolating one onto the other
would smooth the very derivative being measured, so each source is differenced
along its own densities and the shared isotherms are compared.

Writes ``outputs/rho_scaling.json``.
"""

from __future__ import annotations

import argparse
import collections
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np

from config import Config, load_config
from opacity.build_tables import compute_point_means
from opacity.means import rosseland_mean
from tools.text_checks import MARK_WN, _load

BAND_CM = (100.0, 4000.0)
GAP_WINDOW_CM = (1000.0, 1200.0)
MARK_HALFWIDTH_CM = 50.0  # for the robustness average around MARK_WN
MIN_DENSITIES = 3
SOURCES = ("fp", "hitemp", "hitemp_mtckd")

# The most closely spaced pair anywhere on the grid (Delta log10 rho = 0.046).
# Because the index is a difference divided by that spacing, it is the pair where
# point-to-point scatter is most amplified — which is what makes it the scatter
# estimator the brief asks for.
SCATTER_PAIR = (3000.0, 729.344, 810.382)


def _run_parameters(cfg: Config, t_k: float, rho_si: float) -> dict:
    """Molecule count and spectral smoothing width of one first-principles run.

    Neither is uniform over the 48 points — N is 216 at rho <= 60 on the 1000,
    1500 and 2000 K isotherms and 64 everywhere else, and the smoothing kernel is
    15 cm^-1 at four of the 1500 K points and 20 cm^-1 at the rest. A narrower
    kernel lets a transparency window sit deeper, so a kernel that changes along
    an isotherm changes the window value for a reason that is not density. Both
    are recorded per point and checked per isotherm.
    """
    for path in sorted(cfg.fp_opacity_path.glob("T*_opacity.npz")):
        with np.load(path) as d:
            if (abs(float(d["T"]) - t_k) < 0.5
                    and abs(float(d["rho_kg_m3"]) - rho_si) <= max(1e-3, 1e-3 * rho_si)):
                return {"n_molecules": int(d["N"]),
                        "smooth_cm": float(d["smooth_cm"]),
                        "dt_fs": float(d["dt_fs"])}
    return {}


def _diagnostics(cfg: Config, source: str, t_k: float, rho_si: float) -> dict:
    """The three spectral diagnostics at one state, in cm^2/g."""
    wn, kappa = _load(cfg, source, t_k, rho_si)
    if source == "hitemp_mtckd":
        from opacity.mtckd import apply_continuum
        kappa = apply_continuum("mtckd_self", wn, kappa, t_k, rho_si,
                                cfg.mtckd_ref_path)
    order = np.argsort(wn)
    wn, kappa = wn[order], kappa[order]
    m = (wn >= GAP_WINDOW_CM[0]) & (wn <= GAP_WINDOW_CM[1])
    mm = (wn >= MARK_WN - MARK_HALFWIDTH_CM) & (wn <= MARK_WN + MARK_HALFWIDTH_CM)
    return {
        "window_2550": float(np.interp(MARK_WN, wn, kappa)),
        "window_2550_pm50_mean": float(np.trapezoid(kappa[mm], wn[mm])
                                       / (wn[mm][-1] - wn[mm][0])),
        "window_1000_1200": float(np.trapezoid(kappa[m], wn[m])
                                  / (wn[m][-1] - wn[m][0])),
        "rosseland": float(rosseland_mean(wn, kappa, t_k, BAND_CM,
                                          cfg.opacity.rel_floor)),
    }


def _isotherms(cfg: Config, source: str) -> dict[float, list[float]]:
    t, r, _, _ = compute_point_means(source if source != "hitemp_mtckd" else "hitemp",
                                     cfg, BAND_CM)
    g: dict[float, list[float]] = collections.defaultdict(list)
    for a, b in zip(t, r):
        g[float(a)].append(float(b))
    return {k: sorted(v) for k, v in sorted(g.items()) if len(v) >= MIN_DENSITIES}


def _index(rho_a: float, q_a: float, rho_b: float, q_b: float) -> float:
    """d ln alpha / d ln rho between two states, from q = alpha/rho^2.

    alpha = q * rho^2, so d ln alpha / d ln rho = 2 + d ln q / d ln rho.
    """
    return 2.0 + float(np.log(q_b / q_a) / np.log(rho_b / rho_a))


def run(cfg: Config) -> dict:
    keys = ("window_2550", "window_2550_pm50_mean", "window_1000_1200", "rosseland")
    out: dict = {
        "meta": ("alpha/rho^2 = kappa/rho against rho along each isotherm with at "
                 "least three densities, and the local power-law index "
                 "d ln alpha / d ln rho between adjacent densities. Flat in "
                 "alpha/rho^2 (index 2) is exact pair-wise absorption; index 1 is "
                 "the dilute line limit."),
        "band_cm": list(BAND_CM),
        "gap_window_cm": list(GAP_WINDOW_CM),
        "mark_wavenumber_cm": MARK_WN,
        "units": {"q": "cm^5 g^-2 (kappa in cm^2/g divided by rho in g/cm^3)",
                  "kappa": "cm^2/g", "rho": "kg/m^3"},
        "min_densities_per_isotherm": MIN_DENSITIES,
        "sources": {},
    }
    for source in SOURCES:
        per_iso = {}
        for t_k, rhos in _isotherms(cfg, source).items():
            rows = []
            for rho in rhos:
                try:
                    d = _diagnostics(cfg, source, t_k, rho)
                except FileNotFoundError:
                    continue
                rho_cgs = rho / 1000.0
                row = {"rho_kg_m3": rho,
                       "kappa": {k: d[k] for k in keys},
                       "q": {k: d[k] / rho_cgs for k in keys}}
                if source == "fp":
                    row["run"] = _run_parameters(cfg, t_k, rho)
                rows.append(row)
            if len(rows) < MIN_DENSITIES:
                continue
            idx = []
            for a, b in zip(rows, rows[1:]):
                idx.append({
                    "rho_lo_kg_m3": a["rho_kg_m3"], "rho_hi_kg_m3": b["rho_kg_m3"],
                    "rho_geom_kg_m3": float(np.sqrt(a["rho_kg_m3"] * b["rho_kg_m3"])),
                    "dln_rho": float(np.log(b["rho_kg_m3"] / a["rho_kg_m3"])),
                    "index": {k: _index(a["rho_kg_m3"], a["q"][k],
                                        b["rho_kg_m3"], b["q"][k]) for k in keys}})
            entry = {"rows": rows, "index": idx,
                     "index_assigned_to": "geometric mean of the pair"}
            if source == "fp":
                dense = [r for r in rows if r["rho_kg_m3"] >= DENSE_RHO_MIN]
                ns = sorted({r["run"].get("n_molecules") for r in dense})
                sm = sorted({r["run"].get("smooth_cm") for r in dense})
                entry["dense_run_parameters"] = {
                    "n_molecules": ns, "smooth_cm": sm,
                    "homogeneous": len(ns) == 1 and len(sm) == 1}
            per_iso[f"{t_k:g}"] = entry
        out["sources"][source] = per_iso
    out["scatter"] = _scatter(out)
    out["scatter_pooled"] = _scatter_pooled(out)
    out["breakdown"] = _breakdown(out)
    out["fitted_slope_dense"] = _fitted_slope(out)
    out["breakdown_map"] = _lambda_map(cfg)
    return out


# Hard-sphere collision diameter for water; the value the brief specifies.
COLLISION_DIAMETER_M = 2.8e-10
# The two published tables in use and the pressure each is documented to.
TABLE_LIMITS_BAR = {"Freedman et al. 2014": 300.0, "Kempton et al. 2023": 1000.0}
LIMIT_ISOTHERMS_K = (1000.0, 2000.0)


def _lambda_map(cfg: Config) -> dict:
    """Lambda = mean free path / collision diameter, and where the tables stop.

    lambda_mfp = 1/(sqrt(2) n sigma) with sigma = pi d^2, so

        Lambda = lambda_mfp / d = 1 / (sqrt(2) pi d^3 n),  n = rho / m_H2O,

    which depends on rho and not at all on T. Lambda = 1 and Lambda = 0.1 are
    therefore vertical lines on the (T, rho) plane, not contours — a fact worth
    stating rather than drawing round.
    """
    from constants import MMW_H2O, N_A

    m_kg = MMW_H2O / N_A
    denom = np.sqrt(2.0) * np.pi * COLLISION_DIAMETER_M ** 3 / m_kg
    out = {
        "definition": ("Lambda = lambda_mfp / d, lambda_mfp = 1/(sqrt(2) n sigma), "
                       "sigma = pi d^2, d = 2.8 Angstrom, n = rho / m_H2O"),
        "collision_diameter_m": COLLISION_DIAMETER_M,
        "depends_on": "density only; Lambda is independent of temperature",
        "rho_at_lambda_kg_m3": {"1": float(1.0 / denom), "0.1": float(10.0 / denom)},
        "table_limits": {},
    }
    eos = None
    try:
        from eos.aqua import AquaEOS
        eos = AquaEOS.from_config(cfg)
    except Exception as exc:  # noqa: BLE001
        out["table_limits_error"] = f"AQUA unavailable: {exc}"
    if eos is not None:
        for name, p_bar in TABLE_LIMITS_BAR.items():
            per_t = {}
            for t_k in LIMIT_ISOTHERMS_K:
                st = eos.state(p_bar * 1.0e5, t_k)
                per_t[f"{t_k:g}"] = {"rho_kg_m3": float(st.rho),
                                     "lambda": float(1.0 / (denom * st.rho))}
            out["table_limits"][name] = {"pressure_bar": p_bar, "by_isotherm": per_t}
    return out


DENSE_RHO_MIN = 100.0


def _fitted_slope(doc: dict) -> dict:
    """Least-squares d ln alpha / d ln rho over rho >= 100, with a standard error.

    Part (b) asks for adjacent-point indices, and those are what ``breakdown``
    reports. But an adjacent difference divides by a small spacing and throws away
    every other point on the isotherm, so at this scatter level no single pair can
    say much. The same data fitted as a line over the dense end uses all of them,
    and its standard error follows from the same sigma_ln_alpha.
    """
    out = {"rho_min_kg_m3": DENSE_RHO_MIN, "per_source": {}}
    for source, per_src in doc["sources"].items():
        per_iso = {}
        for t_s, iso in per_src.items():
            rows = [r for r in iso["rows"] if r["rho_kg_m3"] >= DENSE_RHO_MIN]
            if len(rows) < 3:
                continue
            for key, sig in ((k, v["sigma_ln_alpha"])
                             for k, v in doc["scatter"]["per_diagnostic"].items()):
                x = np.log(np.array([r["rho_kg_m3"] for r in rows]))
                y = np.log(np.array([r["q"][key] * (r["rho_kg_m3"] / 1000.0) ** 2
                                     for r in rows]))
                slope, _ = np.polyfit(x, y, 1)
                se = float(sig / np.sqrt(((x - x.mean()) ** 2).sum()))
                per_iso.setdefault(t_s, {})[key] = {
                    "slope": float(slope), "standard_error": se,
                    "n_points": len(rows),
                    "sigma_above_2": float((slope - 2.0) / se) if se > 0 else None}
        out["per_source"][source] = per_iso
    return out


def _scatter_pooled(doc: dict) -> dict:
    """sigma_ln_alpha pooled over every first-principles isotherm with >= 4 points.

    The brief's estimator is one pair, so it carries one degree of freedom; this
    is the same quantity with about twenty, as a check that the pair is not an
    outlier in either direction. Each isotherm is fitted with a quadratic in
    ln rho and the residuals pooled.
    """
    out = {"per_diagnostic": {}, "method": ("residuals of a quadratic in ln rho, "
                                            "pooled over first-principles isotherms "
                                            "with >= 4 densities")}
    for key in doc["scatter"]["per_diagnostic"]:
        res, dof = [], 0
        for iso in doc["sources"]["fp"].values():
            rhos = np.array([r["rho_kg_m3"] for r in iso["rows"]])
            if rhos.size < 4:
                continue
            y = np.log(np.array([r["q"][key] for r in iso["rows"]])
                       * (rhos / 1000.0) ** 2)
            x = np.log(rhos)
            res.append(y - np.polyval(np.polyfit(x, y, 2), x))
            dof += rhos.size - 3
        r = np.concatenate(res)
        out["per_diagnostic"][key] = {
            "sigma_ln_alpha": float(np.sqrt((r ** 2).sum() / dof)),
            "dof": int(dof), "n_isotherms": len(res)}
    return out


def _scatter(doc: dict) -> dict:
    """Index noise, set by the most closely spaced pair on the grid.

    The pair's own index is compared with a smooth quadratic in ln rho fitted to
    ln alpha along the whole isotherm. The difference, times the pair's ln rho
    spacing, is the noise in Delta ln alpha; halved by sqrt(2) it is the
    per-point sigma, and the index noise at any other spacing follows as
    sqrt(2) sigma / Delta ln rho.
    """
    t_k, rho_a, rho_b = SCATTER_PAIR
    iso = doc["sources"]["fp"][f"{t_k:g}"]
    rhos = np.array([r["rho_kg_m3"] for r in iso["rows"]])
    out = {"pair": [t_k, rho_a, rho_b], "per_diagnostic": {}}
    for key in ("window_2550", "window_2550_pm50_mean", "window_1000_1200",
                "rosseland"):
        q = np.array([r["q"][key] for r in iso["rows"]])
        ln_alpha = np.log(q * (rhos / 1000.0) ** 2)
        coef = np.polyfit(np.log(rhos), ln_alpha, 2)
        smooth_index = float(np.polyval(np.polyder(coef),
                                        np.log(np.sqrt(rho_a * rho_b))))
        pair = next(i for i in iso["index"]
                    if abs(i["rho_lo_kg_m3"] - rho_a) < 1e-3)
        d_ln_rho = pair["dln_rho"]
        sigma_ln_alpha = abs(pair["index"][key] - smooth_index) * d_ln_rho / np.sqrt(2)
        out["per_diagnostic"][key] = {
            "pair_index": pair["index"][key],
            "smooth_index_at_pair": smooth_index,
            "dln_rho_pair": d_ln_rho,
            "sigma_ln_alpha": float(sigma_ln_alpha),
            "sigma_index_at_pair_spacing": float(np.sqrt(2) * sigma_ln_alpha
                                                 / d_ln_rho)}
    return out


def _breakdown(doc: dict) -> dict:
    """Per isotherm, the first density whose index exceeds 2 by more than noise."""
    out = {}
    for key, per in doc["scatter"]["per_diagnostic"].items():
        sig_ln = per["sigma_ln_alpha"]
        per_iso = {}
        for t_s, iso in doc["sources"]["fp"].items():
            hit = None
            for i in iso["index"]:
                sigma_here = float(np.sqrt(2) * sig_ln / i["dln_rho"])
                if i["index"][key] - 2.0 > sigma_here:
                    hit = {"rho_hi_kg_m3": i["rho_hi_kg_m3"],
                           "rho_geom_kg_m3": i["rho_geom_kg_m3"],
                           "index": i["index"][key],
                           "sigma_index_here": sigma_here}
                    break
            per_iso[t_s] = hit  # None => not within the grid
        out[key] = per_iso
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("config", type=Path)
    args = ap.parse_args(argv)
    cfg = load_config(args.config)
    doc = run(cfg)
    path = cfg.output_path / "rho_scaling.json"
    path.write_text(json.dumps(doc, indent=1) + "\n")

    print("\n  RHO SCALING — index d ln alpha / d ln rho (2 = exact pair-wise)\n")
    for src, per in doc["sources"].items():
        print(f"  {src}")
        for t_s, iso in per.items():
            sk = " ".join(f"{i['rho_geom_kg_m3']:>7.1f}:{i['index']['rosseland']:>5.2f}"
                          for i in iso["index"])
            print(f"    {t_s:>5s} K  {sk}")
    print("\n  scatter (Rosseland):",
          f"sigma_index at the {SCATTER_PAIR[1]:g}/{SCATTER_PAIR[2]:g} spacing = "
          f"{doc['scatter']['per_diagnostic']['rosseland']['sigma_index_at_pair_spacing']:.3g}")
    print(f"\n  wrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
