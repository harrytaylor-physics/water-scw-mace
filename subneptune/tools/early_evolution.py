#!/usr/bin/env python3
"""Radius and entropy over the first 100 Myr, line list against composite (A89).

    python tools/early_evolution.py config/fiducial.yaml

The late-time cooling ratios (S9, ED Fig. 8) say the two opacities give different
ages at the same radius. This asks the complementary question: how early does the
difference start? For the two water worlds and the candidate it reports R_p and
the interior entropy at four ages inside the first 100 Myr under each opacity,
and then the age at which the radii first separate at all.

The answer is that nothing separates in the first 100 Myr, and the reason is in
the composite's own definition rather than in the spectra. The composite blends
first-principles opacity in only above ``opacity.blend_window_kg_m3``; while the
radiative-convective boundary sits below that window the composite *evaluates to
the line list*, so the two tracks are identical to machine precision by
construction. This tool therefore also prints the RCB density against the window,
which is what turns an exact zero from a suspicious result into an explained one.

Two cautions on reading the output. Age zero is the top of the T_int grid, not a
formation state; and the oldest ages on a grid taken down to 12 K are an
extrapolation of the cold end, where the age integral diverges -- they are printed
for completeness and nothing is quoted from them.

Writes ``outputs/early_evolution.json``.
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

AGES_GYR = (0.01, 0.03, 0.05, 0.10)
RCB_PROBE_GYR = (0.01, 0.03, 0.1, 0.3, 1.0, 3.0)
SEPARATION_REARTH = 1.0e-4
GRID = (12.0, 450.0, 70)

CASES = (
    ("water_world_mc10", "water world M_c=10 wmf=0.50 T_eq=500", 10.0, 0.50, 500.0, False),
    ("water_world_mc15", "water world M_c=15 wmf=0.50 T_eq=500", 15.0, 0.50, 500.0, False),
    ("candidate", "candidate M_c=5 f_env=0.05 T_eq=700", 5.0, 0.05, 700.0, True),
)


def _at(age, val, target):
    """Interpolate, returning None outside the range -- np.interp clamps silently,
    which in A88 turned a grid-floor artefact into a difference of exactly zero."""
    order = np.argsort(age)
    a, v = np.asarray(age)[order], np.asarray(val)[order]
    if not (a.min() <= target <= a.max()):
        return None
    return float(np.interp(target, a, v))


def run(cfg: Config, cfg_candidate: Config) -> dict:
    grid = np.geomspace(*GRID[:2], int(GRID[2]))
    out = {"_meta": __doc__.split("\n\n")[0],
           "separation_threshold_rearth": SEPARATION_REARTH,
           "t_int_grid_k": [GRID[0], GRID[1]], "n_nodes": int(GRID[2]),
           "cases": {}}
    for key, label, mcore, f_env, t_eq, is_cand in CASES:
        base = cfg_candidate if is_cand else cfg
        c = dataclasses.replace(
            base,
            core=dataclasses.replace(base.core, mass_mearth=mcore),
            envelope=dataclasses.replace(base.envelope, f_env=f_env))
        window = list(c.opacity.blend_window_kg_m3)
        tr = {s: cooling_track(c, s, t_eq, t_int_grid_k=grid)
              for s in ("hitemp", "composite")}
        ll, co = tr["hitemp"], tr["composite"]
        rec = {
            "label": label, "config": "config/candidate.yaml" if is_cand
            else "config/fiducial.yaml",
            "core_mass_mearth": mcore, "f_env": f_env, "t_eq_k": t_eq,
            "blend_window_kg_m3": window,
            "age_grid_problem": {s: age_grid_problem(t) for s, t in tr.items()},
            "energy_residual": {s: float(t.energy_residual) for s, t in tr.items()},
            "ages_gyr": list(AGES_GYR), "rows": [],
        }
        for a in AGES_GYR:
            r_l, r_c = _at(ll.age_gyr, ll.radius_p_rearth, a), _at(co.age_gyr, co.radius_p_rearth, a)
            s_l, s_c = _at(ll.age_gyr, ll.entropy_j_kg_k, a), _at(co.age_gyr, co.entropy_j_kg_k, a)
            if None in (r_l, r_c, s_l, s_c):
                rec["rows"].append({"age_gyr": a, "outside_track": True})
                continue
            rec["rows"].append({
                "age_gyr": a,
                "radius_rearth": {"hitemp": r_l, "composite": r_c, "delta": r_c - r_l},
                "entropy_j_kg_k": {"hitemp": s_l, "composite": s_c, "delta": s_c - s_l},
                "rcb_rho_kg_m3_composite": _at(co.age_gyr, co.rcb_rho_si, a)})
        rec["rcb_against_window"] = [
            {"age_gyr": a, "rho_rcb_kg_m3": _at(co.age_gyr, co.rcb_rho_si, a)}
            for a in RCB_PROBE_GYR if _at(co.age_gyr, co.rcb_rho_si, a) is not None]
        # where the radii first part company at all
        fine = np.geomspace(0.002, min(ll.age_gyr.max(), co.age_gyr.max()) * 0.999, 400)
        rec["first_separation"] = None
        for a in fine:
            r_l, r_c = _at(ll.age_gyr, ll.radius_p_rearth, a), _at(co.age_gyr, co.radius_p_rearth, a)
            if r_l is None or r_c is None:
                continue
            if abs(r_c - r_l) > SEPARATION_REARTH:
                rec["first_separation"] = {
                    "age_gyr": float(a), "delta_rearth": r_c - r_l,
                    "rho_rcb_kg_m3": _at(co.age_gyr, co.rcb_rho_si, a)}
                break
        out["cases"][key] = rec
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("config", type=Path)
    ap.add_argument("--candidate-config", type=Path,
                    default=Path("config/candidate.yaml"))
    args = ap.parse_args(argv)
    cfg = load_config(args.config)
    doc = run(cfg, load_config(args.candidate_config))
    path = cfg.output_path / "early_evolution.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    atomicio.write_text(path, json.dumps(doc, indent=1) + "\n")

    print(f"\n  EARLY EVOLUTION -> {path}")
    for key, rec in doc["cases"].items():
        print(f"\n  {rec['label']}   [{rec['config']}]  "
              f"blend window {rec['blend_window_kg_m3']} kg/m^3")
        print(f"    {'age(Myr)':>9} {'dR (R_E)':>11} {'dS (J/kg/K)':>13} {'rho_rcb':>9}")
        for row in rec["rows"]:
            if row.get("outside_track"):
                print(f"    {row['age_gyr']*1000:9.0f}   outside the track")
                continue
            print(f"    {row['age_gyr']*1000:9.0f} {row['radius_rearth']['delta']:+11.6f}"
                  f" {row['entropy_j_kg_k']['delta']:+13.3f}"
                  f" {row['rcb_rho_kg_m3_composite']:9.3f}")
        fs = rec["first_separation"]
        print(f"    first |dR| > {SEPARATION_REARTH:g} R_E: "
              + (f"{fs['age_gyr']*1000:.0f} Myr (rho_rcb = {fs['rho_rcb_kg_m3']:.2f})"
                 if fs else "never on this track"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
