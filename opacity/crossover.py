"""FP vs line-list crossover metrics on the merged 48-point grid (A53).

``crossover_hitemp.json`` holds the owner's independent RADIS analysis at the 27
legacy FP points. This module recomputes the same M2 Planck/Rosseland ratios
through *our* pipeline — the identical band-limited table builder used for every
source (SPEC §3), so the differential is exact — and extends them to all 48
points now that the extension campaign has landed.

How the line list is evaluated at an FP point. The extension grid is not a
subset of HITEMP's (T = 650 and ρ = 60/150/500 have no HITEMP twin), so exact
(T, ρ) matching — what ``figures/plot_fx_opacity_compare.py`` does — cannot
reach the new points. The line list is therefore evaluated through its own
``GridOpacitySource``, i.e. triangulated interpolation in (T, log ρ) on log κ,
which is exactly how the structure model reads it. HITEMP spans
T = 500-3000 K, ρ = 0.1-1000 kg/m³, so every FP point is interior in T and all
but the nine with ρ > 1000 are interior in ρ; those nine keep the existing
handling (nearest-neighbour clamp) and are flagged ``ll_clamped`` in every
output rather than silently averaged in.

Agreement with the independent analysis, at the 27 legacy points (A53):

* both coordinates on HITEMP's own grid  -> within **0.7%** on both means;
* T interpolated (T = 2250, 2750)        -> within **3%**;
* ρ clamped (ρ > 1000)                   -> up to **30%** on the Rosseland mean,
  because RADIS computed those at the true density while we hold κ at ρ = 1000.

The first line is the real cross-check and it reproduces A25's ~1% claim.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

from config import Config
from opacity.build_tables import compute_point_means, get_source
from opacity.table import CM2_PER_G_TO_SI

DILUTE_RHO_MAX = 60.0  # matches figures/plot_s7_dilute_floor.py (A52)
LEGACY_GRID_RHO_MIN = 108.796  # the 27-point grid's ρ floor: everything below is new


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
            "is_extension_point": rho < LEGACY_GRID_RHO_MIN or t_k < 1000.0,
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
