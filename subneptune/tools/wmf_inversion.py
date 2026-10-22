#!/usr/bin/env python3
"""What the opacity choice costs an inferred water fraction (A88 task 4, A89).

    python tools/wmf_inversion.py config/fiducial.yaml

A radius is what a transit measures; a water fraction is what somebody infers from
it with an interior model. If the two opacities give different radii for the same
planet, then fitting an observed radius with the line list rather than with the
first-principles composite returns a different water fraction, and that difference
is the part of the opacity question an observer actually feels.

Method, at fixed core mass and fixed age. ``R_comp`` and ``R_LL`` are the radii the
two opacities give at the case's own water fraction ``w0``. An observer fitting the
line list to the composite's radius infers ``w_LL`` with
``R_LL(w_LL) = R_comp(w0)``; linearising about ``w0``,

    dwmf = w0 - w_LL = -(R_comp - R_LL) / (dR/dwmf)

so a NEGATIVE ``dwmf`` means the line list would have overstated the water
fraction and adopting the first-principles composite lowers it. ``dR/dwmf`` is a
central difference of the line-list radius over neighbouring water fractions on
the same core, which is why each case carries a small ``wmf`` grid.

The gradient is used rather than interpolating a bracket because ``np.interp``
clamps at the ends of its grid, and in A88 that silently reported a difference of
exactly zero where the true answer lay outside the grid. A shift that lands
outside the tested range is reported as such and nothing is extrapolated.

Writes ``outputs/wmf_inversion.json``.
"""

from __future__ import annotations

import argparse
import dataclasses
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import atomicio
from config import Config, load_config
from evolution.track import age_grid_problem, cooling_track

AGES_GYR = (1.0, 5.0, 10.0)
# Below this the shift is solver noise, not a displacement: dR between the two
# tracks is ~1e-6 R_E where they agree and dR/dwmf is ~2, so a genuine shift
# cannot be smaller than ~1e-6. Reporting such a row as "outside the tested
# range" because w0 + 1e-7 exceeds 0.50 would be an artefact dressed as a finding.
NOISE_FLOOR_DWMF = 1.0e-5
GRID = (12.0, 450.0, 44)
TESTED_WMF = (0.05, 0.20, 0.50)

# (key, label, core mass, w0, T_eq, wmf grid for the gradient, candidate config)
CASES = (
    ("water_world_mc10", "water world M_c=10", 10.0, 0.50, 500.0, (0.20, 0.30, 0.50), False),
    ("water_world_mc15", "water world M_c=15", 15.0, 0.50, 500.0, (0.20, 0.30, 0.50), False),
    ("candidate", "candidate M_c=5", 5.0, 0.05, 700.0, (0.05, 0.10, 0.20), True),
)


def _at(age, val, target):
    order = np.argsort(age)
    a, v = np.asarray(age)[order], np.asarray(val)[order]
    if not (a.min() <= target <= a.max()):
        return None
    return float(np.interp(target, a, v))


def _radius_at_ages(cfg, source, t_eq, grid):
    tr = cooling_track(cfg, source, t_eq, t_int_grid_k=grid)
    return ({a: _at(tr.age_gyr, tr.radius_p_rearth, a) for a in AGES_GYR},
            age_grid_problem(tr), float(tr.energy_residual))


def run(cfg: Config, cfg_candidate: Config) -> dict:
    grid = np.geomspace(*GRID[:2], int(GRID[2]))
    doc = {"_meta": __doc__.split("\n\n")[0], "ages_gyr": list(AGES_GYR),
           "tested_wmf_range": list(TESTED_WMF),
           "sign_convention": ("dwmf = w0 - w_LL = -(R_comp - R_LL)/(dR/dwmf); "
                               "negative means the line list overstated the water "
                               "fraction"),
           "cases": {}}
    for key, label, mcore, w0, t_eq, wgrid, is_cand in CASES:
        base = cfg_candidate if is_cand else cfg
        def _cfg(w):
            return dataclasses.replace(
                base, core=dataclasses.replace(base.core, mass_mearth=mcore),
                envelope=dataclasses.replace(base.envelope, f_env=w))
        rec = {"label": label, "core_mass_mearth": mcore, "w0": w0, "t_eq_k": t_eq,
               "config": "config/candidate.yaml" if is_cand else "config/fiducial.yaml",
               "wmf_grid": list(wgrid), "failed_points": [], "rows": []}
        r_ll_grid, notes = {}, {}
        for w in wgrid:
            try:
                r, prob, res = _radius_at_ages(_cfg(w), "hitemp", t_eq, grid)
            except Exception as exc:                      # a grid point that will not build
                rec["failed_points"].append({"wmf": w, "error": str(exc)[:200]})
                continue
            r_ll_grid[w] = r
            notes[w] = {"age_grid_problem": prob, "energy_residual": res}
        try:
            r_comp, prob_c, res_c = _radius_at_ages(_cfg(w0), "composite", t_eq, grid)
        except Exception as exc:
            rec["blocked"] = f"composite track failed: {str(exc)[:200]}"
            doc["cases"][key] = rec
            continue
        rec["composite_note"] = {"age_grid_problem": prob_c, "energy_residual": res_c}
        rec["line_list_notes"] = notes
        have = sorted(r_ll_grid)
        for age in AGES_GYR:
            row = {"age_gyr": age}
            rc = r_comp.get(age)
            rl = r_ll_grid.get(w0, {}).get(age)
            if rc is None or rl is None:
                row["blocked"] = "one track does not reach this age"
                rec["rows"].append(row); continue
            # central difference where the neighbours exist, one-sided otherwise
            pts = [(w, r_ll_grid[w][age]) for w in have if r_ll_grid[w][age] is not None]
            if len(pts) < 2:
                row["blocked"] = "fewer than two usable water fractions"
                rec["rows"].append(row); continue
            ws = np.array([p[0] for p in pts]); rs = np.array([p[1] for p in pts])
            i0 = int(np.argmin(abs(ws - w0)))
            grad = float(np.gradient(rs, ws)[i0])
            # w0 is the end of its own grid in every case here, where np.gradient
            # falls back to a three-point one-sided rule. The plain two-point
            # secant to the neighbour is recorded beside it so the sensitivity of
            # the answer to that choice is visible rather than buried.
            j = i0 - 1 if i0 > 0 else i0 + 1
            secant = float((rs[i0] - rs[j]) / (ws[i0] - ws[j]))
            d = -(rc - rl) / grad
            w_new = w0 + d
            tol = 1e-9
            row.update({
                "radius_rearth": {"composite": rc, "hitemp": rl, "delta": rc - rl},
                "line_list_radius_by_wmf": {f"{w:g}": r for w, r in pts},
                "dR_dwmf_rearth": grad,
                "dR_dwmf_rearth_secant": secant,
                "dwmf": d, "dwmf_percent": 100.0 * d / w0,
                "dwmf_secant": -(rc - rl) / secant,
                "implied_wmf": w_new,
                "below_noise_floor": bool(abs(d) < NOISE_FLOOR_DWMF),
                "inside_tested_range": bool(
                    abs(d) < NOISE_FLOOR_DWMF
                    or min(TESTED_WMF) - tol <= w_new <= max(TESTED_WMF) + tol),
            })
            rec["rows"].append(row)
        doc["cases"][key] = rec
    return doc


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("config", type=Path)
    ap.add_argument("--candidate-config", type=Path, default=Path("config/candidate.yaml"))
    args = ap.parse_args(argv)
    cfg = load_config(args.config)
    doc = run(cfg, load_config(args.candidate_config))
    path = cfg.output_path / "wmf_inversion.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    atomicio.write_text(path, json.dumps(doc, indent=1) + "\n")
    print(f"\n  WMF INVERSION -> {path}")
    for key, rec in doc["cases"].items():
        print(f"\n  {rec['label']}  w0 = {rec['w0']:g}  [{rec['config']}]")
        if rec.get("blocked"):
            print(f"    BLOCKED: {rec['blocked']}"); continue
        for f in rec["failed_points"]:
            print(f"    wmf = {f['wmf']:g} did not build: {f['error'][:90]}")
        print(f"    {'age(Gyr)':>8} {'dR (R_E)':>10} {'dR/dwmf':>9} {'dwmf':>10} "
              f"{'dwmf %':>8} {'implied wmf':>12}  range")
        for row in rec["rows"]:
            if row.get("blocked"):
                print(f"    {row['age_gyr']:8.0f}   {row['blocked']}"); continue
            print(f"    {row['age_gyr']:8.0f} {row['radius_rearth']['delta']:+10.5f}"
                  f" {row['dR_dwmf_rearth']:9.3f} {row['dwmf']:+10.5f}"
                  f" {row['dwmf_percent']:+8.3f} {row['implied_wmf']:12.5f}"
                  f"  {'no shift' if row['below_noise_floor'] else ('inside' if row['inside_tested_range'] else 'OUTSIDE tested range')}"
                  f"   (secant dwmf {row['dwmf_secant']:+.5f})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
