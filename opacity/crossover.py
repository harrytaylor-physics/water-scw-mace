"""FP vs line-list crossover metrics on the merged 48-point grid (A53).

``crossover_hitemp.json`` holds the owner's independent RADIS analysis, which
covers 27 of the 48 first-principles states. This module recomputes the same
Planck/Rosseland ratios through *our* pipeline — the identical band-limited table
builder used for every source (SPEC §3), so the differential is exact — over all
48 points, which are treated throughout as one grid.

How the line list is evaluated at an FP point. The grid is not a subset of
HITEMP's (T = 650 and ρ = 60/150/500 have no HITEMP twin), so exact
(T, ρ) matching — what ``figures/plot_fx_opacity_compare.py`` does — cannot
reach the new points. The line list is therefore evaluated through its own
``GridOpacitySource``, i.e. triangulated interpolation in (T, log ρ) on log κ,
which is exactly how the structure model reads it. HITEMP spans
T = 500-3000 K, ρ = 0.1-1000 kg/m³, so every FP point is interior in T and all
but the nine with ρ > 1000 are interior in ρ; those nine keep the existing
handling (nearest-neighbour clamp) and are flagged ``ll_clamped`` in every
output rather than silently averaged in.

Agreement with the independent analysis, at the 27 shared states (A53):

* both coordinates on HITEMP's own grid  -> within **0.7%** on both means;
* T interpolated (T = 2250, 2750)        -> within **3%**;
* ρ clamped (ρ > 1000)                   -> up to **30%** on the Rosseland mean.

The first line is the real cross-check and it reproduces A25's ~1% claim.

The third line was originally attributed to the density clamp alone. A59 measured
the clamp against direct RADIS runs at the true density — **-1.62%** at
ρ/ρ_ceiling = 1.20 and **-5.89%** at 2.23, reproducing those points' deviations to
better than 0.3 percentage points — so the clamp is a **<=6%** effect, not a 30%
one. The two 20-30% outliers sit at T = 2250 and 2750 K, the only temperatures off
HITEMP's own T grid, so they are T-interpolated as well as ρ-clamped; the residual
is unexplained and this repository cannot separate an interpolation error from a
difference in the delivered analysis at those states.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

from config import Config
from opacity.build_tables import compute_point_means, get_source
from opacity.mtckd import (
    LBLRTM_LINE_BOUND_CM,
    decomposition_is_well_posed,
    lorentz_hwhm_cm,
)
from opacity.table import CM2_PER_G_TO_SI

DILUTE_RHO_MAX = 60.0  # matches figures/plot_s7_dilute_floor.py (A52)


def _reference(cfg: Config) -> dict:
    """The owner's RADIS metrics, keyed by (T, ρ), or {} if absent."""
    path = cfg.fp_opacity_path / "crossover_hitemp.json"
    if not path.exists():
        return {}
    data = json.loads(path.read_text())["results"]
    return {(float(r["T"]), float(r["rho"])): r for r in data}


def _floor_fractions(cfg: Config) -> dict:
    """Per-point window floor fractions from the S7 run (A52), or {} if not run."""
    path = cfg.output_path / "S7_dilute_floor.csv"
    if not path.exists():
        return {}
    out = {}
    with path.open() as fh:
        for row in csv.DictReader(fh):
            key = (float(row["temperature_k"]), float(row["rho_si"]))
            out[key] = {
                "frac_above_window": float(row["frac_above_window"]),
                "frac_above_cores": float(row["frac_above_cores"]),
                "rosseland_window_share": float(row["rosseland_window_share"]),
            }
    return out


def _match(reference: dict, t_k: float, rho: float) -> dict | None:
    """Reference entry for (T, ρ), tolerating the JSON's rounded coordinates."""
    for (t_ref, rho_ref), entry in reference.items():
        if abs(t_ref - t_k) < 0.5 and abs(rho_ref - rho) <= max(1e-3, 1e-3 * rho):
            return entry
    return None


def crossover_rows(
    cfg: Config,
    band_cm: tuple[float, float],
    fp_name: str = "fp",
    line_list: str = "hitemp",
) -> list[dict]:
    """One row per FP grid point: both sources' means, their ratios, provenance.

    κ values are in cm²/g (what ``compute_point_means`` returns); the line list is
    converted back from the source's SI so both sides are quoted in one unit.
    """
    t_fp, rho_fp, kr_fp, kp_fp = compute_point_means(fp_name, cfg, band_cm)
    order = np.lexsort((rho_fp, t_fp))
    t_fp, rho_fp, kr_fp, kp_fp = (a[order] for a in (t_fp, rho_fp, kr_fp, kp_fp))

    ll = get_source(line_list, cfg, band=band_cm)
    kr_ll = np.asarray(ll.kappa_R(t_fp, rho_fp)) / CM2_PER_G_TO_SI
    kp_ll = np.asarray(ll.kappa_P(t_fp, rho_fp)) / CM2_PER_G_TO_SI
    clamped = np.asarray(ll.clamped_mask(t_fp, rho_fp)).ravel()

    t_grid = np.unique(np.asarray(ll.table.temperature_k, float))
    rho_grid = np.unique(np.asarray(ll.table.rho_si, float))
    reference = _reference(cfg)
    floors = _floor_fractions(cfg)

    rows = []
    for i in range(t_fp.size):
        t_k, rho = float(t_fp[i]), float(rho_fp[i])
        on_t = bool(np.any(np.isclose(t_grid, t_k, rtol=1e-9)))
        on_rho = bool(np.any(np.isclose(rho_grid, rho, rtol=1e-9)))
        row = {
            "T": t_k,
            "rho": rho,
            "kappa_R_fp_cm2_g": float(kr_fp[i]),
            "kappa_P_fp_cm2_g": float(kp_fp[i]),
            f"kappa_R_{line_list}_cm2_g": float(kr_ll[i]),
            f"kappa_P_{line_list}_cm2_g": float(kp_ll[i]),
            "M2_planck_ratio": float(kp_fp[i] / kp_ll[i]),
            "M2_rosseland_ratio": float(kr_fp[i] / kr_ll[i]),
            # Spectral contrast per source: κ_P/κ_R is large exactly when the
            # band has deep transparency windows (the harmonic mean collapses
            # while the arithmetic one does not). Comparing the two sources'
            # contrast says whether a ratio is an opacity difference or a
            # window-depth difference (A53).
            "contrast_fp": float(kp_fp[i] / kr_fp[i]),
            f"contrast_{line_list}": float(kp_ll[i] / kr_ll[i]),
            "ll_clamped": bool(clamped[i]),
            "ll_on_grid": on_t and on_rho,
            "is_dilute": rho <= DILUTE_RHO_MAX,
        }
        ref = _match(reference, t_k, rho)
        if ref is not None:
            row["ref_planck_ratio"] = float(ref["M2_planck_ratio"])
            row["ref_rosseland_ratio"] = float(ref["M2_rosseland_ratio"])
            row["dev_planck_pct"] = 100.0 * (
                row["M2_planck_ratio"] / row["ref_planck_ratio"] - 1.0)
            row["dev_rosseland_pct"] = 100.0 * (
                row["M2_rosseland_ratio"] / row["ref_rosseland_ratio"] - 1.0)
        floor = floors.get((t_k, rho))
        if floor is not None:
            row.update({f"floor_{k}": v for k, v in floor.items()})
        rows.append(row)
    return rows


def write_json(cfg: Config, rows: list[dict], band_cm: tuple[float, float],
               line_list: str = "hitemp", name: str = "crossover_hitemp_v2.json") -> Path:
    """Write the v2 metrics to ``outputs/`` (a computed product, not delivered data)."""
    out_dir = cfg.output_path
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / name
    payload = {
        "meta": (
            f"FP vs {line_list.upper()} band-limited means recomputed through the "
            f"subneptune table builder on the merged 48-point FP grid (A53). Band "
            f"[{band_cm[0]:g}, {band_cm[1]:g}] cm^-1, rel_floor {cfg.opacity.rel_floor:g}. "
            f"{line_list.upper()} is evaluated by triangulated interpolation in "
            f"(T, log rho) on log kappa — the same call the structure model makes. "
            f"Points with ll_clamped=true are outside its rho range (>1000 kg/m^3) "
            f"and hold kappa at the nearest grid value."
        ),
        "band_cm": list(band_cm),
        "rel_floor": cfg.opacity.rel_floor,
        "line_list": line_list,
        "n_points": len(rows),
        "results": rows,
    }
    path.write_text(json.dumps(payload, indent=1) + "\n")
    return path


def agreement_summary(rows: list[dict]) -> dict:
    """Max |deviation| from the independent analysis, split by how κ_LL was got."""
    groups = {"on_grid": [], "interpolated": [], "clamped": []}
    for r in rows:
        if "dev_planck_pct" not in r:
            continue
        key = ("clamped" if r["ll_clamped"]
               else "on_grid" if r["ll_on_grid"] else "interpolated")
        groups[key].append(r)
    out = {}
    for key, these in groups.items():
        out[key] = {
            "n": len(these),
            "max_abs_dev_planck_pct": max(
                (abs(r["dev_planck_pct"]) for r in these), default=0.0),
            "max_abs_dev_rosseland_pct": max(
                (abs(r["dev_rosseland_pct"]) for r in these), default=0.0),
        }
    return out


# --------------------------------------------------------------------------
# v3: both references — bare line list, and line list + MT_CKD self-continuum
# --------------------------------------------------------------------------
# The 48 points are ONE grid. They are grouped by where they sit in (T, rho),
# because the caveats differ by regime, not by when a point was computed (A62):
#   "dilute"     — rho <= 60, where A52 found the transparency windows that carry
#                  kappa_R sit at or below the dipole model's noise floor;
#   "cold_dense" — T <= 1000 K and rho <= 500, the regime the structural claim
#                  lives in and where the MT_CKD continuum dominates (A58);
#   "dense"      — the rest: hotter, denser, and where the line list's own
#                  broadening convention starts to dominate (A59).
GROUPS = ("dilute", "cold_dense", "dense")
COLD_DENSE_T_MAX = 1000.0
COLD_DENSE_RHO_MAX = 500.0


def _ideal_pressure_pa(temperature_k: float, rho_si: float) -> float:
    """Ideal-gas P at the cell's density — the broadening pressure the delivered
    HITEMP spectra were built at (A17, A59)."""
    from constants import MMW_H2O

    return float(rho_si) * 8.314462618 * float(temperature_k) / MMW_H2O


def _group_of(row: dict) -> str:
    if row["is_dilute"]:
        return "dilute"
    if row["T"] <= COLD_DENSE_T_MAX and row["rho"] <= COLD_DENSE_RHO_MAX:
        return "cold_dense"
    return "dense"


def crossover_v3_rows(
    cfg: Config,
    band_cm: tuple[float, float],
    fp_name: str = "fp",
    line_list: str = "hitemp",
    continuum_source: str = "hitemp_mtckd",
) -> list[dict]:
    """One row per FP point with BOTH references: bare line list and line+MT_CKD.

    Everything is computed by the identical code path — same band, same
    ``compute_point_means``, same triangulated (T, log rho) interpolation of the
    line-list table onto the FP points — so the two ratio columns differ only by
    the continuum (A58). kappa values are cm^2/g.

    ``mtckd_share_kappa_R`` is 1 - kappa_R(line) / kappa_R(line+continuum). In
    the grey limit kappa_R(line+cont) = kappa_l + kappa_c and this reduces
    exactly to continuum / (line + continuum); away from it, it is the fraction
    of the combined Rosseland mean the continuum is responsible for, which is the
    quantity that matters for the structure model.
    """
    t_fp, rho_fp, kr_fp, kp_fp = compute_point_means(fp_name, cfg, band_cm)
    order = np.lexsort((rho_fp, t_fp))
    t_fp, rho_fp, kr_fp, kp_fp = (a[order] for a in (t_fp, rho_fp, kr_fp, kp_fp))

    ll = get_source(line_list, cfg, band=band_cm)
    mt = get_source(continuum_source, cfg, band=band_cm)
    kr_ll = np.asarray(ll.kappa_R(t_fp, rho_fp)) / CM2_PER_G_TO_SI
    kp_ll = np.asarray(ll.kappa_P(t_fp, rho_fp)) / CM2_PER_G_TO_SI
    kr_mt = np.asarray(mt.kappa_R(t_fp, rho_fp)) / CM2_PER_G_TO_SI
    kp_mt = np.asarray(mt.kappa_P(t_fp, rho_fp)) / CM2_PER_G_TO_SI
    clamped = np.asarray(ll.clamped_mask(t_fp, rho_fp)).ravel()

    t_grid = np.unique(np.asarray(ll.table.temperature_k, float))
    rho_grid = np.unique(np.asarray(ll.table.rho_si, float))
    rho_ceiling = float(rho_grid.max())

    rows = []
    for i in range(t_fp.size):
        t_k, rho = float(t_fp[i]), float(rho_fp[i])
        row = {
            "T": t_k,
            "rho": rho,
            "kappa_R_fp_cm2_g": float(kr_fp[i]),
            "kappa_P_fp_cm2_g": float(kp_fp[i]),
            "kappa_R_line_cm2_g": float(kr_ll[i]),
            "kappa_P_line_cm2_g": float(kp_ll[i]),
            "kappa_R_line_mtckd_cm2_g": float(kr_mt[i]),
            "kappa_P_line_mtckd_cm2_g": float(kp_mt[i]),
            # FP / bare line list (the manuscript's present Fig. 2b)
            "ratio_planck_bare": float(kp_fp[i] / kp_ll[i]),
            "ratio_rosseland_bare": float(kr_fp[i] / kr_ll[i]),
            # FP / (line list + MT_CKD self-continuum)
            "ratio_planck_mtckd": float(kp_fp[i] / kp_mt[i]),
            "ratio_rosseland_mtckd": float(kr_fp[i] / kr_mt[i]),
            # what fraction of the combined mean the continuum carries
            "mtckd_share_kappa_R": float(1.0 - kr_ll[i] / kr_mt[i]),
            "mtckd_share_kappa_P": float(1.0 - kp_ll[i] / kp_mt[i]),
            "ll_clamped": bool(clamped[i]),
            "ll_on_grid": bool(np.any(np.isclose(t_grid, t_k, rtol=1e-9))
                               and np.any(np.isclose(rho_grid, rho, rtol=1e-9))),
            "is_dilute": rho <= DILUTE_RHO_MAX,
            # A rho-clamped point holds kappa at the ceiling. The continuum's
            # ABSORPTION is quadratic in density, but its MASS OPACITY is linear
            # (kappa = alpha/rho), and it is kappa that is tabulated — so the
            # continuum's contribution is under-counted by rho/rho_ceiling, not
            # by its square (A58).
            "continuum_underscale": float(max(1.0, rho / rho_ceiling)),
            # Whether "line list + continuum" is a decomposition of anything here.
            # MT_CKD is the residual beyond a line-by-line calculation truncated at
            # 25 cm^-1 with the plinth removed; where the Lorentz half-width exceeds
            # that bound, the split is not well posed and the delivered spectra
            # (RADIS, truncated at 50 cm^-1, no plinth subtraction) double-count the
            # 25-50 cm^-1 region on top. See opacity/mtckd.py.
            "lorentz_hwhm_cm": lorentz_hwhm_cm(t_k, _ideal_pressure_pa(t_k, rho)),
            "decomposition_well_posed": decomposition_is_well_posed(
                t_k, _ideal_pressure_pa(t_k, rho)),
        }
        row["group"] = _group_of(row)
        rows.append(row)
    return rows


def summarize_v3(rows: list[dict]) -> dict:
    """Per-group min/max of each ratio, plus the sub-unity continuum-inclusive set."""
    keys = ("ratio_planck_bare", "ratio_rosseland_bare",
            "ratio_planck_mtckd", "ratio_rosseland_mtckd", "mtckd_share_kappa_R")
    out = {"groups": {}, "all": {}}
    for label, these in [(g, [r for r in rows if r["group"] == g]) for g in GROUPS] + [
        ("all", rows)
    ]:
        stats = {"n": len(these)}
        for key in keys:
            vals = [r[key] for r in these]
            stats[key] = {"min": min(vals), "max": max(vals)} if vals else None
        if label == "all":
            out["all"] = stats
        else:
            out["groups"][label] = stats
    out["sub_unity_mtckd"] = [
        {k: r[k] for k in ("T", "rho", "group", "ratio_planck_mtckd",
                           "ratio_rosseland_mtckd", "mtckd_share_kappa_R", "ll_clamped")}
        for r in rows
        if r["ratio_planck_mtckd"] < 1.0 or r["ratio_rosseland_mtckd"] < 1.0
    ]
    return out


def write_v3_json(
    cfg: Config, rows: list[dict], band_cm: tuple[float, float],
    line_list: str = "hitemp", continuum_source: str = "hitemp_mtckd",
    name: str = "crossover_v3.json",
) -> Path:
    """Write the two-reference metrics + per-group summary to ``outputs/``."""
    out_dir = cfg.output_path
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / name
    payload = {
        "meta": (
            f"FP vs {line_list.upper()} on the 48-point FP grid, against TWO references "
            f"(A58): the bare line list, and the line list plus the MT_CKD 4.3 water "
            f"SELF-continuum evaluated at each spectrum's own (T, rho). Both are "
            f"band-limited means over [{band_cm[0]:g}, {band_cm[1]:g}] cm^-1 computed by "
            f"the identical code path, with the line-list table read by the same "
            f"triangulated (T, log rho) interpolation the structure model uses, so the "
            f"two ratio columns differ ONLY by the continuum. Foreign continuum is "
            f"identically zero: it carries (1 - x_H2O) and the envelope is pure steam. "
            f"The continuum is EXTRAPOLATED far outside MT_CKD's fitted range in both T "
            f"(fitted ~250-350 K, used 500-3000 K) and density (used up to ~1350x the "
            f"1013 mb / 296 K reference). mtckd_share_kappa_R = 1 - kappa_R(line) / "
            f"kappa_R(line+continuum), which reduces to continuum/(line+continuum) in "
            f"the grey limit. Points with ll_clamped=true sit above the line list's rho "
            f"ceiling; kappa there is held at the ceiling, so the continuum's "
            f"contribution is UNDER-counted by continuum_underscale = rho/rho_ceiling "
            f"(linear, because the tabulated quantity is the MASS opacity kappa = "
            f"alpha/rho, and alpha is the quadratic one). decomposition_well_posed "
            f"is false where the self-broadened Lorentz half-width exceeds LBLRTM's "
            f"{LBLRTM_LINE_BOUND_CM:g} cm^-1 line-by-line bound, which is the "
            f"convention MT_CKD is defined as the residual to: there 'line list + "
            f"continuum' is not a decomposition of anything. The delivered spectra "
            f"(RADIS, truncated at 50 cm^-1, profile RENORMALISED to unit area, no "
            f"plinth subtraction) disagree with MT_CKD's convention about the "
            f"25-50 cm^-1 region, and additionally fold the >50 cm^-1 wing back into "
            f"it. The mismatch is real; its SIGN is not determined, because "
            f"renormalisation conserves line strength (A59)."
        ),
        "band_cm": list(band_cm),
        "rel_floor": cfg.opacity.rel_floor,
        "line_list": line_list,
        "continuum_source": continuum_source,
        "mtckd_reference": str(cfg.mtckd_ref_path),
        "n_points": len(rows),
        "summary": summarize_v3(rows),
        "results": rows,
    }
    path.write_text(json.dumps(payload, indent=1) + "\n")
    return path
