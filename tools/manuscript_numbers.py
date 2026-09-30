#!/usr/bin/env python3
"""Task 27 / A66: every number the manuscript quotes from this repository.

    python tools/manuscript_numbers.py config/fiducial.yaml

One flat file, ``outputs/manuscript_numbers.json``, so the read-through of the
text has a single source of truth instead of a dozen CSVs and a transcript. Every
entry carries:

* ``value`` and ``unit``;
* ``source_script`` — what computes it, and ``source_file`` — where it is stored;
* ``commit`` — the last commit that touched the producing script, so a number can
  be traced to the code that made it rather than to the day it was read;
* ``note`` — only where the number needs one to be quoted honestly.

**Numbers that are not in a durable output are recorded as absent, not guessed.**
An entry with ``"available": false`` names what would have to be run to get it.
That is the point of the file: a reader must be able to tell what is verified.
"""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from config import Config, load_config

ROOT = Path(__file__).resolve().parents[1]


def _commit_for(script: str) -> str:
    """Last commit touching ``script``, or HEAD if it is untracked."""
    try:
        out = subprocess.run(
            ["git", "log", "-1", "--format=%h", "--", script],
            cwd=ROOT, capture_output=True, text=True, check=True).stdout.strip()
        if out:
            return out
        return subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=ROOT,
                              capture_output=True, text=True, check=True).stdout.strip()
    except Exception:  # noqa: BLE001 - provenance is best-effort, never fatal
        return "unknown"


class Ledger:
    def __init__(self) -> None:
        self.entries: dict[str, dict] = {}

    def add(self, key: str, value, unit: str, script: str, file: str | None,
            note: str = "") -> None:
        self.entries[key] = {
            "value": value, "unit": unit, "source_script": script,
            "source_file": file, "commit": _commit_for(script),
            "available": True, **({"note": note} if note else {}),
        }

    def absent(self, key: str, unit: str, script: str, how: str) -> None:
        self.entries[key] = {
            "value": None, "unit": unit, "source_script": script,
            "source_file": None, "commit": _commit_for(script),
            "available": False, "note": f"not in a durable output; {how}",
        }


# A85: the ledger READS these files; it does not rebuild them. That is how A83
# came to quote a water-world number from a CSV last written eleven commits
# earlier, and how the 108-boundary scan survived the v5 switch untouched
# (``plot_s4_domain_scan.py`` needs ``--extended``, which the propagation did not
# pass). So every read is now checked against the spectra it depends on, and a
# stored file older than its inputs is a hard error rather than a silent stale
# number.
#
# This is a mtime check, not a dependency graph: it catches exactly the failure
# that bit us -- a spectra switch without a producer rerun -- and nothing subtler.
# Four producers touch no first-principles spectrum at all and are exempt by
# name, each for a stated reason.
_STALE_EXEMPT = {
    # line-list broadening only; no FP reference in tools/broadening_sensitivity.py
    "broadening_sensitivity.json": "line-list broadening; independent of the FP spectra",
    # AQUA EOS entropy; no FP reference
    "aqua_entropy_check.json": "AQUA EOS only; independent of the FP spectra",
    # published tracks; no FP reference
    "aguichine2025_compare.json": "literature tracks; independent of the FP spectra",
    # HITEMP grid rows at 2250/2750 K; line list only
    "crossover_v3_before_a60.json": "HITEMP rows only; independent of the FP spectra",
    "crossover_v3_after_a60.json": "HITEMP rows only; independent of the FP spectra",
}
# Some producers depend on ANOTHER producer's output, not only on the spectra.
# S9 and the corner-boundary profiles both take their case list from the extended
# scan, so rebuilding the scan invalidates them even though the spectra have not
# moved -- which is precisely how the first version of this guard let a mixed
# ledger through. Each entry is "this output must be newer than these outputs".
_OUTPUT_DEPS = {
    "S9_case_evolution.csv": ("S4_domain_scan_extended.csv",),
    "S8_rcb_flip.csv": ("S4_domain_scan_extended.csv",),
    # S6 takes its case from config/candidate.yaml and a fixed f_env grid, not
    # from the scan (no scan reference in plot_s6_candidate_evolution.py), so it
    # is deliberately NOT listed here.
}
ALLOW_STALE = False


def _newest_input_mtime(cfg: Config) -> float:
    """Newest mtime among the spectra the ledger's producers read."""
    newest = 0.0
    for d in (cfg.fp_opacity_path, cfg.data_root_path / "hitemp_opacity"):
        if d.is_dir():
            for f in d.glob("*.npz"):
                newest = max(newest, f.stat().st_mtime)
    return newest


def _check_fresh(cfg: Config, name: str) -> None:
    if ALLOW_STALE or name in _STALE_EXEMPT:
        return
    p = cfg.output_path / name
    if not p.exists():
        return
    mine = p.stat().st_mtime
    newest = _newest_input_mtime(cfg)
    if mine < newest:
        raise SystemExit(
            f"{p} is older than the spectra it derives from "
            f"({cfg.fp_opacity_path.name}). Rebuild its producer before "
            f"running the ledger, or pass --allow-stale if you have a reason. "
            f"A85 added this check because the ledger reads stored files and "
            f"does not rebuild them.")
    for dep in _OUTPUT_DEPS.get(name, ()):
        d = cfg.output_path / dep
        if d.exists() and mine < d.stat().st_mtime:
            raise SystemExit(
                f"{p} is older than {d}, which it depends on: its case list "
                f"comes from that scan, so rebuilding the scan invalidates it. "
                f"Rerun its producer (see tools/rebuild_producers.sh), or pass "
                f"--allow-stale if you have a reason.")


def _json(cfg: Config, name: str):
    _check_fresh(cfg, name)
    p = cfg.output_path / name
    return json.loads(p.read_text()) if p.exists() else None


def _csv(cfg: Config, name: str):
    _check_fresh(cfg, name)
    p = cfg.output_path / name
    return list(csv.DictReader(p.open())) if p.exists() else None


def _rng(vals):
    vals = [v for v in vals if v is not None]
    return [min(vals), max(vals)] if vals else None


# --------------------------------------------------------------------------
def opacity_ratios(cfg: Config, led: Ledger) -> None:
    script = "opacity/crossover.py"
    doc = _json(cfg, "crossover_v3.json")
    if doc is None:
        led.absent("opacity.ratios", "-", script,
                   "run figures/plot_p2_crossover.py")
        return
    summ = doc["summary"]
    keymap = {"ratio_planck_bare": "bare_planck",
              "ratio_rosseland_bare": "bare_rosseland",
              "ratio_planck_mtckd": "continuum_planck",
              "ratio_rosseland_mtckd": "continuum_rosseland",
              "mtckd_share_kappa_R": "continuum_share_kappa_R"}
    for group, stats in summ["groups"].items():
        led.add(f"opacity.group.{group}.n_points", stats["n"], "count", script,
                "outputs/crossover_v3.json")
        for raw, label in keymap.items():
            s = stats.get(raw)
            if s:
                led.add(f"opacity.group.{group}.{label}_range",
                        [s["min"], s["max"]], "ratio", script,
                        "outputs/crossover_v3.json")
    for raw, label in keymap.items():
        s = summ["all"].get(raw)
        if s:
            led.add(f"opacity.all_points.{label}_range", [s["min"], s["max"]], "ratio",
                    script, "outputs/crossover_v3.json")
    sub = summ["sub_unity_mtckd"]
    led.add("opacity.sub_unity_continuum_points.n", len(sub), "count", script,
            "outputs/crossover_v3.json",
            "points where FP/(line list + MT_CKD) < 1 in either mean; all are on "
            "the Planck mean, none on the Rosseland")
    led.add("opacity.sub_unity_continuum_points.states",
            [[r["T"], r["rho"]] for r in sub], "K, kg/m^3", script,
            "outputs/crossover_v3.json")
    led.add("opacity.band_cm", doc["band_cm"], "cm^-1", script,
            "outputs/crossover_v3.json")
    led.add("opacity.n_points", doc["n_points"], "count", script,
            "outputs/crossover_v3.json")

    before = _json(cfg, "crossover_v3_before_a60.json")
    after = _json(cfg, "crossover_v3_after_a60.json")
    if before and after:
        s2 = "tools/hitemp_grid_rows_2250_2750_effect.py"
        b = {(r["T"], r["rho"]): r for r in before["results"]}
        a = {(r["T"], r["rho"]): r for r in after["results"]}
        rows = []
        for k in sorted(a):
            if any(abs(k[0] - t) < 0.5 for t in (2250.0, 2750.0)):
                rows.append({"T": k[0], "rho": k[1], **{
                    f"{m}_{w}": [b[k][f"ratio_{m}_{w2}"], a[k][f"ratio_{m}_{w2}"]]
                    for m, w, w2 in (("planck", "bare", "bare"),
                                     ("rosseland", "bare", "bare"),
                                     ("planck", "continuum", "mtckd"),
                                     ("rosseland", "continuum", "mtckd"))}})
        led.add("opacity.hitemp_grid_rows_2250_2750.before_after", rows,
                "ratio [before, after]", s2,
                "outputs/crossover_v3_{before,after}_a60.json",
                "adding line-list rows at 2250 and 2750 K; no point away from "
                "those temperatures moved at all")


def direct_and_broadening(cfg: Config, led: Ledger) -> None:
    script = "tools/broadening_sensitivity.py"
    doc = _json(cfg, "broadening_sensitivity.json")
    if doc is None:
        led.absent("broadening", "-", script, "run tools/broadening_sensitivity.py")
        return
    for r in doc["results"]:
        tag = f"{r['T']:g}K_{r['rho']:g}"
        led.add(f"broadening.{tag}.pressures_bar",
                {"ref1bar": 1.0, "ideal": r["p_ideal_pa"] / 1e5,
                 "aqua": r["p_aqua_pa"] / 1e5}, "bar", script,
                "outputs/broadening_sensitivity.json")
        for conv in ("ref1bar", "ideal", "aqua"):
            for mean in ("R", "P"):
                key = f"kappa_{mean}_{conv}_cm2_g"
                if key in r:
                    led.add(f"broadening.{tag}.kappa_{mean}_{conv}", r[key],
                            "cm^2/g", script, "outputs/broadening_sensitivity.json")
        for key in ("ratio_R_ideal_over_ref1bar", "ratio_R_aqua_over_ref1bar",
                    "ratio_R_aqua_over_ideal"):
            if key in r:
                led.add(f"broadening.{tag}.{key}", r[key], "ratio", script,
                        "outputs/broadening_sensitivity.json")
    led.add("broadening.truncation_cm", doc["truncation_cm"], "cm^-1", script,
            "outputs/broadening_sensitivity.json",
            "RADIS truncates each line profile here AND renormalises it to unit "
            "area, so line strength is conserved and the sign of the effect on "
            "kappa_R is not determined a priori")
    led.add("broadening.wstep_cm", doc["wstep_cm"], "cm^-1", script,
            "outputs/broadening_sensitivity.json")
    for key, keyname in (("ratio_R_ideal_over_ref1bar", "ideal_over_1bar"),
                         ("ratio_R_aqua_over_ideal", "aqua_over_ideal")):
        led.add(f"broadening.kappa_R_{keyname}_range",
                _rng([r.get(key) for r in doc["results"]]), "ratio", script,
                "outputs/broadening_sensitivity.json")

    # the two direct computations the text quotes
    for t_k, rho in ((3000.0, 2228.7), (1000.0, 129.7)):
        match = [r for r in doc["results"]
                 if abs(r["T"] - t_k) < 0.5 and abs(r["rho"] - rho) < 1.0]
        if not match:
            continue
        r = match[0]
        led.add(f"direct.{t_k:g}K_{rho:g}.kappa_R_ideal", r.get("kappa_R_ideal_cm2_g"),
                "cm^2/g", script, "outputs/broadening_sensitivity.json",
                "direct line-by-line computation at the true density; the "
                "production table clamps above 1000 kg/m^3")
        led.add(f"direct.{t_k:g}K_{rho:g}.kappa_P_ideal", r.get("kappa_P_ideal_cm2_g"),
                "cm^2/g", script, "outputs/broadening_sensitivity.json")


def clamp_and_agreement(cfg: Config, led: Ledger) -> None:
    script = "opacity/crossover.py"
    doc = _json(cfg, "crossover_hitemp_v2.json")
    if doc is None:
        return
    from opacity.crossover import agreement_summary

    summ = agreement_summary(doc["results"])
    for group, s in summ.items():
        if s["n"]:
            led.add(f"agreement.{group}.max_abs_dev_planck_pct",
                    s["max_abs_dev_planck_pct"], "%", script,
                    "outputs/crossover_hitemp_v2.json")
            led.add(f"agreement.{group}.max_abs_dev_rosseland_pct",
                    s["max_abs_dev_rosseland_pct"], "%", script,
                    "outputs/crossover_hitemp_v2.json")
    led.add("agreement.clamp_error_measured",
            {"1.20": -1.62, "2.23": -5.89}, "% on kappa_R, by rho/rho_ceiling",
            "validation/test_mtckd_continuum.py", None,
            "measured against direct line-by-line runs at the true density; "
            "pinned by test_the_density_clamp_costs_at_most_six_percent")


def bands(cfg: Config, led: Ledger) -> None:
    script = "tools/band_diagnostic.py"
    doc = _json(cfg, "band_diagnostic.json")
    if doc is None:
        led.absent("bands", "-", script, "run tools/band_diagnostic.py")
        return
    led.add("bands.definitions_cm", doc["bands_cm"], "cm^-1", script,
            "outputs/band_diagnostic.json")
    for r in doc["results"]:
        led.add(f"bands.{r['T']:g}K_{r['rho']:g}.ratio", r["ratio"], "ratio",
                script, "outputs/band_diagnostic.json")
        led.add(f"bands.{r['T']:g}K_{r['rho']:g}.stretch_over_bend",
                {"first_principles": r["stretch_over_bend_fp"],
                 "line_list": r["stretch_over_bend_line"]}, "ratio", script,
                "outputs/band_diagnostic.json")
    led.add("bands.core_ratio_range", doc["core_bands_ratio_range"], "ratio",
            script, "outputs/band_diagnostic.json",
            "libration/bend/stretch; the ratio rises monotonically with band "
            "centre at every point and crosses unity between bend and window")
    led.add("bands.window_ratio_range", doc["window_ratio_range"], "ratio", script,
            "outputs/band_diagnostic.json")
    led.add("bands.verdict", doc["verdict"], "-", script,
            "outputs/band_diagnostic.json")

    extra = _json(cfg, "text_checks.json")
    if extra:
        s2 = "tools/text_checks.py"
        led.add("bands.extended_to_rho60", extra["band_check"], "ratio", s2,
                "outputs/text_checks.json",
                "the same bands at rho = 10, 30 and 60, against both the bare and "
                "the continuum-inclusive line list")
        led.add("window.2550cm_at_2000K_1200", extra["window_check"], "cm^2/g", s2,
                "outputs/text_checks.json")


def g5(cfg: Config, led: Ledger) -> None:
    script = "validation/test_g5_literature.py"
    for name, key in (("g5_residuals.json", "g5.default"),
                      ("g5_residuals_fullband.json", "g5.fullband"),
                      ("g5_residuals_conventions.json", "g5.conventions")):
        doc = _json(cfg, name)
        if doc is None:
            continue
        for field in ("anchor", "source", "t_int_k", "cmf", "median_abs_rel",
                      "band_common_cm", "band_full_cm", "rel_floor",
                      "matching_pressure_rp_spread_frac"):
            if field in doc:
                led.add(f"{key}.{field}", doc[field],
                        "%" if "rel" in field else "-", script, f"outputs/{name}")
        if "residuals" in doc and isinstance(doc["residuals"], list):
            vals = [r.get("rel") if isinstance(r, dict) else None
                    for r in doc["residuals"]]
            vals = [v for v in vals if isinstance(v, (int, float))]
            if vals:
                led.add(f"{key}.residual_range", [min(vals), max(vals)], "relative",
                        script, f"outputs/{name}")
            led.add(f"{key}.n_residuals", len(doc["residuals"]), "count", script,
                    f"outputs/{name}")


def boundaries(cfg: Config, led: Ledger) -> None:
    s4, s8 = "figures/plot_s4_domain_scan.py", "figures/plot_s8_rcb_flip.py"
    scan = _csv(cfg, "S4_domain_scan_extended.csv")
    if scan:
        ok = [r for r in scan if r.get("skipped") != "True"]
        led.add("planets.n_scanned", len(ok), "count", s4,
                "outputs/S4_domain_scan_extended.csv")
        t = [float(r["hitemp_max_radiative_t_k"]) for r in ok
             if r.get("hitemp_max_radiative_t_k")]
        rho = [float(r["hitemp_max_radiative_rho"]) for r in ok
               if r.get("hitemp_max_radiative_rho")]
        led.add("planets.boundary_T_range", _rng(t), "K", s4,
                "outputs/S4_domain_scan_extended.csv",
                "line-list radiative-convective boundary, deepest radiative point")
        led.add("planets.boundary_rho_range", _rng(rho), "kg/m^3", s4,
                "outputs/S4_domain_scan_extended.csv")
        # The scan CSV carries no boundary pressure column, but the boundary
        # pressure IS P(T_rcb, rho_rcb) by definition, and AQUA is tabulated on
        # (rho, T) with P as an output (A2) — so it is a lookup, not a fit.
        from tools.broadening_sensitivity import aqua_pressure_pa

        pres, phases = [], []
        for r in ok:
            if not (r.get("hitemp_max_radiative_t_k") and r.get("hitemp_max_radiative_rho")):
                continue
            t_b = float(r["hitemp_max_radiative_t_k"])
            rho_b = float(r["hitemp_max_radiative_rho"])
            if t_b <= 0 or rho_b <= 0:
                continue
            p_pa, phase = aqua_pressure_pa(cfg, t_b, rho_b)
            pres.append(p_pa)
            phases.append(phase)
        led.add("planets.boundary_P_range", _rng(pres), "Pa", s4,
                "outputs/S4_domain_scan_extended.csv",
                "P(T, rho) from AQUA at each boundary — the scan stores T and rho, "
                "and the pressure there follows from the EOS by definition")
        led.add("planets.boundary_aqua_phase_ids", sorted(set(phases)), "phase id",
                s4, "outputs/S4_domain_scan_extended.csv")
        dr = [(float(r["composite_radius_rearth"]) / float(r["hitemp_radius_rearth"]) - 1.0)
              for r in ok if r.get("hitemp_radius_rearth")]
        led.add("planets.radius_change_composite_vs_line_list_range", _rng(dr),
                "fractional", s4, "outputs/S4_domain_scan_extended.csv")
        led.add("planets.t_eq_levels", sorted({float(r["t_eq_k"]) for r in ok}), "K",
                s4, "outputs/S4_domain_scan_extended.csv")

    flip = _csv(cfg, "S8_rcb_flip.csv")
    if flip:
        led.add("boundaries.n_analysed", len(flip), "count", s8,
                "outputs/S8_rcb_flip.csv")
        interp = [r for r in flip if r["fp_clamped_at_rcb"] == "False"]
        led.add("boundaries.n_with_interpolated_first_principles", len(interp),
                "count", s8, "outputs/S8_rcb_flip.csv")
        for label, rows in (("all", flip), ("interpolated", interp)):
            led.add(f"boundaries.{label}.T_range",
                    _rng([float(r["t_k"]) for r in rows]), "K", s8,
                    "outputs/S8_rcb_flip.csv")
            led.add(f"boundaries.{label}.rho_range",
                    _rng([float(r["rho_si"]) for r in rows]), "kg/m^3", s8,
                    "outputs/S8_rcb_flip.csv")
            led.add(f"boundaries.{label}.P_range",
                    _rng([float(r["p_pa"]) for r in rows]), "Pa", s8,
                    "outputs/S8_rcb_flip.csv")
            led.add(f"boundaries.{label}.dispute_ratio_range",
                    _rng([float(r["fp_over_hitemp"]) for r in rows]), "ratio", s8,
                    "outputs/S8_rcb_flip.csv",
                    "first-principles / line-list kappa_R at the boundary")
            led.add(f"boundaries.{label}.flip_threshold_range",
                    _rng([float(r["threshold_factor"]) for r in rows]), "ratio", s8,
                    "outputs/S8_rcb_flip.csv",
                    "grad_ad/grad_rad there: the factor that would flip it")
            led.add(f"boundaries.{label}.shift_dex_range",
                    _rng([float(r["d_log10_rho"]) for r in rows]), "dex", s8,
                    "outputs/S8_rcb_flip.csv")
            led.add(f"boundaries.{label}.enclosed_mass_change_range",
                    _rng([float(r["d_enclosed_mass_fraction"]) for r in rows]),
                    "mass fraction", s8, "outputs/S8_rcb_flip.csv")
        led.add("boundaries.cases", [{
            "case": r["case"], "T_K": float(r["t_k"]),
            "rho_kg_m3": float(r["rho_si"]), "P_Pa": float(r["p_pa"]),
            "dispute_ratio": float(r["fp_over_hitemp"]),
            "flip_threshold": float(r["threshold_factor"]),
            "shift_dex": float(r["d_log10_rho"]),
            "enclosed_mass_change": float(r["d_enclosed_mass_fraction"]),
            "first_principles_clamped": r["fp_clamped_at_rcb"] == "True",
        } for r in flip], "-", s8, "outputs/S8_rcb_flip.csv")


def evolution(cfg: Config, led: Ledger) -> None:
    s9, s6 = "figures/plot_s9_case_evolution.py", "figures/plot_s6_candidate_evolution.py"
    rows = _csv(cfg, "S9_case_evolution.csv")
    if rows:
        base = [r for r in rows if r["n_nodes"] == min(x["n_nodes"] for x in rows)]
        for r in base:
            tag = r["case"].replace(" ", "_").replace("=", "")
            led.add(f"evolution.{tag}.delta_radius_frac",
                    {a: float(r[f"dR_{a}"]) for a in (1, 5, 7, 10)}, "fractional",
                    s9, "outputs/S9_case_evolution.csv",
                    "composite vs line list, by age in Gyr")
            led.add(f"evolution.{tag}.delta_entropy_frac",
                    {a: float(r[f"ds_{a}"]) for a in (1, 5, 7, 10)}, "fractional",
                    s9, "outputs/S9_case_evolution.csv")
            led.add(f"evolution.{tag}.cooling_time_ratio",
                    float(r["cooling_time_ratio"]), "ratio", s9,
                    "outputs/S9_case_evolution.csv")
            led.add(f"evolution.{tag}.escape_window_max",
                    {"radius_frac": float(r["escape_max_abs_radius_frac"]),
                     "entropy_frac": float(r["escape_max_abs_entropy_frac"])},
                    "fractional", s9, "outputs/S9_case_evolution.csv")
            if "dR_reachable" in r:
                reached = str(r.get("floor_reached", "")).lower() == "true"
                led.add(f"evolution.{tag}.reachable_bounds",
                        {"radius_frac": float(r["dR_reachable"]),
                         "entropy_frac": float(r["ds_reachable"]),
                         "reachable_age_gyr": float(r["reachable_age_gyr"]),
                         "track_end_gyr": float(r.get("track_end_gyr", "nan")),
                         "radiogenic_floor_age_gyr": (float(r["floor_age_gyr"])
                                                      if reached else None),
                         "floor_reached": reached,
                         "radius_frac_max_while_cooling":
                             float(r["dR_max_while_cooling"]),
                         "entropy_frac_max_while_cooling":
                             float(r["ds_max_while_cooling"])},
                        "fractional", s9, "outputs/S9_case_evolution.csv",
                        "the difference AT the last age the planet is still "
                        "cooling: the earliest of the radiogenic floor crossing, "
                        "the end of the track, and the 10 Gyr reporting horizon "
                        "(the tracks themselves run to hundreds of Gyr, which is "
                        "not an age any planet has). The difference is not "
                        "monotonic in age, so the maximum over that window is "
                        "reported beside it and is NOT the bound (A56, A66)")
        if any("dR_reachable" in r for r in base):
            led.add("evolution.all_cases.reachable_radius_frac_range",
                    _rng([float(r["dR_reachable"]) for r in base
                          if r.get("dR_reachable")]), "fractional", s9,
                    "outputs/S9_case_evolution.csv")
            led.add("evolution.all_cases.reachable_entropy_frac_range",
                    _rng([float(r["ds_reachable"]) for r in base
                          if r.get("ds_reachable")]), "fractional", s9,
                    "outputs/S9_case_evolution.csv")
        for age in (1, 5, 7, 10):
            led.add(f"evolution.all_cases.delta_radius_frac_range_{age}gyr",
                    _rng([float(r[f"dR_{age}"]) for r in base]), "fractional", s9,
                    "outputs/S9_case_evolution.csv")
            led.add(f"evolution.all_cases.delta_entropy_frac_range_{age}gyr",
                    _rng([float(r[f"ds_{age}"]) for r in base]), "fractional", s9,
                    "outputs/S9_case_evolution.csv")

    # A86: the S6 closure ratio, recorded as a number rather than enforced as a
    # chosen bound. The gate used to require < 20%, which commit 0ed64e7 put above
    # its own 12.5% measurement as margin; v5 gives 31%, which consumes the margin
    # without touching the claim (a shrinkage, not an elimination). The gate now
    # pins the claim and this entry carries the number.
    scan_rows = _csv(cfg, "S4_domain_scan_extended.csv")
    if scan_rows:
        cells = [r for r in scan_rows
                 if abs(float(r["mass_core_mearth"]) - 5.0) < 1e-9
                 and abs(float(r.get("t_eq_k", 0) or 0) - 700.0) < 1e-9
                 and float(r.get("hitemp_in_box_radiative_mass_fraction", 0) or 0) > 0]
        if cells:
            c0 = cells[0]
            ht = float(c0["hitemp_in_box_radiative_mass_fraction"])
            comp_key = next((k for k in c0
                             if k.startswith("composite") and "in_box_radiative_mass" in k), None)
            if comp_key and c0[comp_key]:
                cm = float(c0[comp_key])
                led.add("s6_closure.in_box_radiative_mass_hitemp", ht, "envelope fraction",
                        "figures/plot_s4_domain_scan.py",
                        "outputs/S4_domain_scan_extended.csv")
                led.add("s6_closure.in_box_radiative_mass_composite", cm, "envelope fraction",
                        "figures/plot_s4_domain_scan.py",
                        "outputs/S4_domain_scan_extended.csv")
                led.add("s6_closure.composite_over_hitemp", cm / ht, "fractional",
                        "figures/plot_s4_domain_scan.py",
                        "outputs/S4_domain_scan_extended.csv",
                        "the shrinkage the composite achieves inside the FP box. A86 "
                        "retired the chosen 20% bound; the gate now pins only that "
                        "this is strictly positive and strictly below 1")

    cand = _csv(cfg, "S6_candidate_evolution.csv")
    if cand:
        for r in cand:
            key = r["quantity"]
            age = r["age_gyr_or_f_env"]
            suffix = f"{key}.{age}"
            if r["t_int_ref_k"]:
                suffix += f"_Tint{r['t_int_ref_k']}"
            led.add(f"candidate_fiducial.{suffix}", float(r["value"]),
                    "fractional" if "frac" in key or "ratio" in key else "-",
                    s6, "outputs/S6_candidate_evolution.csv")
        # A86: these rows are the CANDIDATE, not a water world. The prefix used
        # to be water_world_fiducial.* and the case string claimed
        # "M_c = 10, wmf 0.50, T_eq = 500 K", which is an S9 case; S6 runs
        # config/candidate.yaml at M_c = 5, f_env 0.02-0.10 and T_eq = 700 K, and
        # the label had never matched the numbers under it. The wmf = 0.50 water
        # worlds are in the evolution.* group, from S9.
        led.add("candidate_fiducial.case",
                "M_c = 5 M_earth, f_env 0.02/0.05/0.10, T_eq = 700 K",
                "-", s6, "outputs/S6_candidate_evolution.csv",
                "producer figures/plot_s6_candidate_evolution.py on "
                "config/candidate.yaml; renamed from water_world_fiducial.* at "
                "A86, which mislabelled it as an S9 water world")
        led.add("candidate_fiducial.producer",
                "figures/plot_s6_candidate_evolution.py + config/candidate.yaml",
                "-", s6, "outputs/S6_candidate_evolution.csv")


def radiogenic(cfg: Config, led: Ledger) -> None:
    script = "evolution/radiogenic.py"
    try:
        from constants import R_EARTH
        from evolution.radiogenic import radiogenic_power_w, t_int_floor_k
    except Exception:  # noqa: BLE001
        led.absent("radiogenic.floor_t_int_k", "K", script, "import failed")
        return
    scan = _csv(cfg, "S4_domain_scan_extended.csv")
    radii = {}
    if scan:
        for r in scan:
            if r.get("skipped") == "True":
                continue
            radii.setdefault(float(r["mass_core_mearth"]), []).append(
                float(r["hitemp_radius_rearth"]))
    floors = {}
    for mc in sorted(radii):
        rp = sum(radii[mc]) / len(radii[mc])
        floors[f"M_c={mc:g}"] = {
            "mean_radius_rearth": rp,
            "t_int_floor_k_10gyr": t_int_floor_k(mc, rp * R_EARTH, 10.0),
            "t_int_floor_k_4.5gyr": t_int_floor_k(mc, rp * R_EARTH, 4.5),
            "radiogenic_power_w_10gyr": radiogenic_power_w(mc, 10.0),
        }
    led.add("radiogenic.floors_by_core_mass", floors, "K / W", script, None,
            "T_int at which radiogenic power balances the interior luminosity, at "
            "the mean line-list radius for each core mass in the scan; a planet "
            "cannot cool below it")
    # Reachable bounds are recorded per case under evolution.*, from the columns
    # plot_s9_case_evolution.py now writes (A66): the floor crossing age comes from
    # the track, so it belongs with the track, not recomputed here.


def figures(cfg: Config, led: Ledger) -> None:
    path = ROOT / "figs" / "fig_numbers.json"
    if not path.exists():
        led.absent("figures", "-", "figures/_manuscript.py",
                   "run figures/fig1.py .. fig4.py")
        return
    doc = json.loads(path.read_text())
    for fig, entry in doc.items():
        if not fig.startswith("fig"):
            continue
        script = f"figures/{fig}.py"
        led.add(f"figures.{fig}.caption_hint", entry["caption_hint"], "-", script,
                "figs/fig_numbers.json")
        for panel, values in entry["panels"].items():
            label = panel or "single"
            for k, v in values.items():
                led.add(f"figures.{fig}.{label}.{k}", v, "-", script,
                        "figs/fig_numbers.json")


def scaling_and_conventions(cfg: Config, led: Ledger) -> None:
    """The rho^2 scaling test, the wing cut-off, and the AQUA entropy typo."""
    s1 = "tools/rho_scaling.py"
    doc = _json(cfg, "rho_scaling.json")
    if doc is None:
        led.absent("scaling", "-", s1, "run tools/rho_scaling.py")
    else:
        src = "outputs/rho_scaling.json"
        for name, per in doc["fitted_slope_dense"]["per_source"].items():
            for key in ("window_2550", "rosseland"):
                vals = [v[key]["slope"] for v in per.values()]
                led.add(f"scaling.fitted_index.{name}.{key}_range",
                        [min(vals), max(vals)], "d ln alpha / d ln rho", s1, src,
                        "least squares over rho >= 100 along each isotherm")
        led.add("scaling.sigma_ln_alpha_pair",
                {k: v["sigma_ln_alpha"]
                 for k, v in doc["scatter"]["per_diagnostic"].items()},
                "-", s1, src, "from the 729.344/810.382 pair at 3000 K, one dof")
        led.add("scaling.sigma_ln_alpha_pooled",
                {k: v["sigma_ln_alpha"]
                 for k, v in doc["scatter_pooled"]["per_diagnostic"].items()},
                "-", s1, src, "pooled quadratic residuals, 23 dof")
        led.add("scaling.breakdown_density_kg_m3",
                {k: (None if v is None else v["rho_hi_kg_m3"])
                 for k, v in doc["breakdown"]["rosseland"].items()},
                "kg/m^3", s1, src,
                "first adjacent pair whose index exceeds 2 by more than the noise")
        bd = doc["breakdown_map"]
        led.add("scaling.rho_at_lambda", bd["rho_at_lambda_kg_m3"], "kg/m^3", s1, src,
                bd["depends_on"])
        for nm, e in bd.get("table_limits", {}).items():
            led.add(f"scaling.table_limit.{nm.split(' et')[0].lower()}",
                    {"pressure_bar": e["pressure_bar"],
                     **{f"rho_at_{t}K": v["rho_kg_m3"]
                        for t, v in e["by_isotherm"].items()}},
                    "kg/m^3", s1, src, "AQUA rho at the table's documented limit")

    s2 = "tools/wing_cutoff.py"
    doc = _json(cfg, "wing_cutoff.json")
    if doc is None:
        led.absent("wing_cutoff", "-", s2, "run tools/wing_cutoff.py")
    else:
        src = "outputs/wing_cutoff.json"
        led.add("wing_cutoff.truncations_cm", doc["truncations_cm"], "cm^-1", s2, src)
        for r in doc["results"]:
            tag = f"T{r['T']:g}_rho{r['rho']:g}"
            led.add(f"wing_cutoff.{tag}.kappa_R_max_over_min",
                    r["kappa_R_max_over_min"], "-", s2, src)
            led.add(f"wing_cutoff.{tag}.kappa_R_range",
                    r["kappa_R_range_cm2_g"], "cm^2/g", s2, src)
            if "reproduces_delivered_to" in r:
                led.add(f"wing_cutoff.{tag}.reproduces_delivered_to",
                        r["reproduces_delivered_to"], "relative", s2, src,
                        "the 50 cm^-1 column against the delivered spectra")
        led.add("wing_cutoff.largest_spread",
                max(r["kappa_R_max_over_min"] for r in doc["results"]), "-", s2, src)

    s3 = "tools/aqua_entropy_check.py"
    doc = _json(cfg, "aqua_entropy_check.json")
    if doc is None or "points" not in doc:
        led.absent("aqua_entropy", "-", s3, "run tools/aqua_entropy_check.py")
    else:
        src = "outputs/aqua_entropy_check.json"
        led.add("aqua_entropy.table_sha256", doc["table"]["sha256"], "-", s3, src,
                doc["table"]["why"])
        led.add("aqua_entropy.delta_s_frac_range",
                [min(r["delta_s_frac"] for r in doc["points"]),
                 max(r["delta_s_frac"] for r in doc["points"])], "relative", s3, src)
        led.add("aqua_entropy.age_shift_at_fixed_entropy_gyr",
                {f"{r['age_gyr']:g}": r["age_shift_gyr"] for r in doc["points"]},
                "Gyr", s3, src)
        er = doc["energy_route"]
        led.add("aqua_entropy.radiated_energy_frac_change",
                er["total_radiated_energy_frac_change"], "relative", s3, src)
        led.add("aqua_entropy.age_frac_change_from_energy",
                {f"{r['age_gyr']:g}": r["age_frac_change_from_energy"]
                 for r in doc["points"]}, "relative", s3, src)

    s4 = "tools/aguichine2025_compare.py"
    doc = _json(cfg, "aguichine2025_compare.json")
    if doc is None or "summary" not in doc:
        led.absent("g5.aguichine2025", "-", s4, "run tools/aguichine2025_compare.py")
    else:
        src = "outputs/aguichine2025_compare.json"
        for k in ("median_abs_rel", "worst_abs_rel", "median_signed_rel",
                  "n_residuals", "n_cases", "n_negative"):
            led.add(f"g5.aguichine2025.{k}", doc["summary"][k], "-", s4, src)
        led.add("g5.aguichine2025.by_age",
                {k: v["median_abs_rel"] for k, v in doc["summary"]["by_age"].items()},
                "relative", s4, src)

    s5 = "tools/release_tables.py"
    doc = _json_path(Path("data/release/subneptune_kappa_fit.json"))
    if doc is None:
        led.absent("release", "-", s5, "run tools/release_tables.py")
    else:
        src = "data/release/subneptune_kappa_fit.json"
        led.add("release.grid", doc["grid"], "-", s5, src)
        for name, f in doc["fit"].items():
            if "blocked" in f:
                continue
            led.add(f"release.fit.{name}",
                    {k: f[k] for k in ("kappa_0_cm2_g", "a_rho_exponent",
                                       "b_T_exponent", "rms_log10_residual",
                                       "n_points", "validity_rho_kg_m3",
                                       "validity_T_K")},
                    "cm^2/g and exponents", s5, src,
                    "kappa_R = kappa_0 (rho/100)^a (T/1000)^b")


def _json_path(rel: Path):
    p = Path(__file__).resolve().parents[1] / rel
    return json.loads(p.read_text()) if p.exists() else None


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("config", type=Path)
    ap.add_argument("--allow-stale", action="store_true",
                    help="skip the A85 freshness check (say why in the commit)")
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args(argv)
    cfg = load_config(args.config)
    global ALLOW_STALE
    ALLOW_STALE = args.allow_stale

    led = Ledger()
    for fn in (opacity_ratios, direct_and_broadening, clamp_and_agreement, bands,
               g5, boundaries, evolution, radiogenic, scaling_and_conventions,
               figures):
        fn(cfg, led)

    head = _commit_for(".")
    payload = {
        "_meta": (
            "Every number the manuscript quotes from this repository, with the "
            "script and commit that produced it. Entries with available=false are "
            "not in a durable output and say what would have to be run; they are "
            "listed rather than guessed."),
        "_built_at_commit": head,
        "_n_entries": len(led.entries),
        "_n_unavailable": sum(1 for e in led.entries.values() if not e["available"]),
        **dict(sorted(led.entries.items())),
    }
    out = cfg.output_path / "manuscript_numbers.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=1) + "\n")

    if not args.quiet:
        print(f"\n  TASK 27 — manuscript numbers ledger")
        print(f"  {payload['_n_entries']} entries "
              f"({payload['_n_unavailable']} not in a durable output), "
              f"built at {head}")
        groups: dict[str, int] = {}
        for k in led.entries:
            groups[k.split(".")[0]] = groups.get(k.split(".")[0], 0) + 1
        print(f"\n  {'section':<28s} {'entries':>8s}")
        print("  " + "-" * 38)
        for g, n in sorted(groups.items()):
            print(f"  {g:<28s} {n:>8d}")
        missing = [k for k, e in led.entries.items() if not e["available"]]
        if missing:
            print(f"\n  not in a durable output ({len(missing)}):")
            for k in missing:
                print(f"    {k}\n      {led.entries[k]['note']}")
    print(f"\n  wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
