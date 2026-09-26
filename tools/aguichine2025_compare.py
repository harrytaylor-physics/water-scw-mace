#!/usr/bin/env python3
"""Radius-age tracks against Aguichine et al. (2025) (Task 40b).

    python tools/aguichine2025_compare.py config/fiducial.yaml
    python tools/aguichine2025_compare.py config/fiducial.yaml --cases 4

The 2021 anchor the G5 gate uses is a *static* mass-radius relation digitized
from a figure. The 2025 Steam Worlds Evolution release
(Zenodo 10.5281/zenodo.15043384, arXiv 2412.17945) is a machine-readable grid of
*evolutionary tracks* — 1440 of them — so it can be matched age for age instead
of at one assumed interior state. That is a strictly better anchor, and this
measures whether the coverage supports making it the gate's reference.

**Conventions, matched rather than assumed.**

* Radius: their R_1mibar is the radius at 1 microbar = 0.1 Pa, which is the
  level the G5 anchor already reads (``anchor.radius_pressure_bar = 1e-6``), so
  the primary comparison needs no conversion at all. Their R_20mbar and R_b
  (1000 bar) are carried too; for a 20 M_E water world R_20mbar/R_b is 1.013,
  so the convention is worth about 1%, not tens of percent.
* Mass: they tabulate the planet mass M_p; this model is parameterized by the
  core mass, so M_c = (1 - WMF) M_p.
* Core: "Earth-like core (32.5% iron-sulfur alloy, 67.5% mantle)", so cmf =
  0.325 — not the cmf = 0 the 2021 no-metallic-core case needed.
* Star: the G rows. Stellar type enters only their atmospheric radiative
  transfer; both are recorded and the M-G spread is reported as a systematic.
* Age zero: they state the planet forms at t_0 = 1 Myr and age is host-star age.
  This model's age is measured from its own hot start, so ages are compared
  directly and the 1 Myr offset is negligible at 1-10 Gyr (0.1% at 1 Gyr).

Writes ``outputs/aguichine2025_compare.json``.
"""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import replace
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np

from config import Config, load_config

REL = Path("aguichine2025") / "Aguichine2024_SWE_all.dat"
DOI = "10.5281/zenodo.15043384"
COLS = ("ST", "WMF", "T_eq", "M_p", "Age", "R_b", "R_20mbar", "R_1mibar",
        "T_b", "T_int", "s_env", "I_MR2", "L_env", "L_rad", "L_core", "L_grav",
        "u_tot")
AGES_GYR = (1.0, 5.0, 10.0)
THEIR_CMF = 0.325
RADIUS_LEVEL_BAR = 1.0e-6  # their R_1mibar, and the G5 anchor's own level
T_INT_RANGE_K = (18.0, 200.0)
N_NODES = 20
# Matched cases: their own grid values, spanning T_eq, water fraction and mass.
CASES = tuple({"star": "G", "wmf": w, "t_eq_k": t, "m_p_earth": m}
              for t in (500.0, 600.0, 700.0)
              for w in (0.2, 0.5)
              for m in (9.66586, 20.0))


def load_tracks(cfg: Config) -> dict:
    path = Path(cfg.data_root_path) / REL
    if not path.exists():
        return {}
    out: dict = {}
    for ln in path.open():
        if ln.startswith("#") or not ln.strip():
            continue
        f = ln.split()
        if len(f) < len(COLS):
            continue
        key = (f[0], round(float(f[1]), 6), round(float(f[2]), 3),
               round(float(f[3]), 6))
        out.setdefault(key, []).append([float(x) for x in f[4:]])
    return {k: np.array(sorted(v, key=lambda r: r[0])) for k, v in out.items()}


def their_radius(track: np.ndarray, age_gyr: float, column: str) -> float:
    """Their radius at one age, log-interpolated in age. NaN outside the track."""
    idx = {"R_b": 1, "R_20mbar": 2, "R_1mibar": 3}[column]
    age, r = track[:, 0], track[:, idx]
    if age_gyr < age.min() or age_gyr > age.max():
        return float("nan")
    return float(np.interp(np.log(age_gyr), np.log(age), r))


def our_track(cfg: Config, case: dict, source: str = "composite"):
    from evolution.track import cooling_track

    m_core = (1.0 - case["wmf"]) * case["m_p_earth"]
    c = replace(
        cfg,
        core=replace(cfg.core, mass_mearth=m_core, cmf=THEIR_CMF),
        envelope=replace(cfg.envelope, f_env=case["wmf"]),
        boundary=replace(cfg.boundary, p_top_bar=RADIUS_LEVEL_BAR),
    )
    grid = np.geomspace(*T_INT_RANGE_K, N_NODES)
    return cooling_track(c, source, case["t_eq_k"], t_int_grid_k=grid), m_core


def run(cfg: Config, cases) -> dict:
    tracks = load_tracks(cfg)
    out = {
        "meta": __doc__.split("\n\n")[0],
        "reference": {"doi": DOI, "file": str(REL),
                      "n_tracks": len(tracks),
                      "radius_level_bar": RADIUS_LEVEL_BAR,
                      "their_column": "R_1mibar",
                      "core": f"cmf = {THEIR_CMF} (their Earth-like core)"},
        "ages_gyr": list(AGES_GYR),
        "cases": [],
    }
    if not tracks:
        out["blocked"] = (f"reference tracks not found at {cfg.data_root_path / REL}; "
                          f"download from Zenodo {DOI}")
        return out

    for case in cases:
        key = (case["star"], round(case["wmf"], 6), round(case["t_eq_k"], 3),
               round(case["m_p_earth"], 6))
        if key not in tracks:
            out["cases"].append({**case, "skipped": "no matching reference track"})
            continue
        ref = tracks[key]
        try:
            track, m_core = our_track(cfg, case)
        except Exception as exc:  # noqa: BLE001
            out["cases"].append({**case, "skipped": f"model failed: {exc}"})
            continue
        rows = []
        for a in AGES_GYR:
            ours = float(track.radius_at_age(a))
            theirs = their_radius(ref, a, "R_1mibar")
            rows.append({
                "age_gyr": a, "ours_rearth": ours, "theirs_rearth": theirs,
                "rel": (ours / theirs - 1.0) if np.isfinite(theirs) else None,
                "theirs_R_b_rearth": their_radius(ref, a, "R_b"),
                "theirs_R_20mbar_rearth": their_radius(ref, a, "R_20mbar"),
            })
        # the same case round an M star, as the stellar-type systematic
        m_key = ("M",) + key[1:]
        star_spread = None
        if m_key in tracks:
            g5 = their_radius(ref, 5.0, "R_1mibar")
            m5 = their_radius(tracks[m_key], 5.0, "R_1mibar")
            if np.isfinite(g5) and np.isfinite(m5):
                star_spread = float(m5 / g5 - 1.0)
        out["cases"].append({
            **case, "mass_core_mearth": m_core, "points": rows,
            "reference_star_type_spread_at_5gyr": star_spread,
            "our_age_range_gyr": [float(np.min(track.age_gyr)),
                                  float(np.max(track.age_gyr))],
            "reference_age_range_gyr": [float(ref[0, 0]), float(ref[-1, 0])],
        })

    rel = [p["rel"] for c in out["cases"] for p in c.get("points", [])
           if p.get("rel") is not None]
    if rel:
        a = np.abs(np.array(rel))
        out["summary"] = {
            "n_residuals": int(a.size),
            "n_cases": sum(1 for c in out["cases"] if "points" in c),
            "median_abs_rel": float(np.median(a)),
            "worst_abs_rel": float(a.max()),
            "median_signed_rel": float(np.median(rel)),
            "n_negative": int(sum(1 for x in rel if x < 0)),
            "by_age": {f"{ag:g}": {
                "median_abs_rel": float(np.median(np.abs([
                    p["rel"] for c in out["cases"] for p in c.get("points", [])
                    if p.get("rel") is not None and p["age_gyr"] == ag]))),
                "n": int(sum(1 for c in out["cases"] for p in c.get("points", [])
                             if p.get("rel") is not None and p["age_gyr"] == ag))}
                for ag in AGES_GYR},
        }
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("config", type=Path)
    ap.add_argument("--cases", type=int, default=len(CASES))
    args = ap.parse_args(argv)
    cfg = load_config(args.config)
    doc = run(cfg, CASES[: args.cases])
    path = cfg.output_path / "aguichine2025_compare.json"
    path.write_text(json.dumps(doc, indent=1) + "\n")

    if "blocked" in doc:
        print(f"\n  BLOCKED: {doc['blocked']}")
        return 0
    print(f"\n  AGUICHINE 2025 — {doc['reference']['n_tracks']} reference tracks, "
          f"radius at {RADIUS_LEVEL_BAR:g} bar\n")
    print(f"  {'T_eq':>5s} {'wmf':>5s} {'M_p':>7s} {'M_c':>6s} {'age':>5s} "
          f"{'ours':>7s} {'theirs':>7s} {'rel':>8s}")
    for c in doc["cases"]:
        if "points" not in c:
            print(f"  {c['t_eq_k']:5.0f} {c['wmf']:5.2f} {c['m_p_earth']:7.3f} "
                  f"      -     -       -       -   {c['skipped']}")
            continue
        for p in c["points"]:
            rel = "n/a" if p["rel"] is None else f"{p['rel']:+7.2%}"
            print(f"  {c['t_eq_k']:5.0f} {c['wmf']:5.2f} {c['m_p_earth']:7.3f} "
                  f"{c['mass_core_mearth']:6.2f} {p['age_gyr']:5.0f} "
                  f"{p['ours_rearth']:7.4f} {p['theirs_rearth']:7.4f} {rel:>8s}")
    s = doc.get("summary")
    if s:
        print(f"\n  median |rel| {s['median_abs_rel']:.4f}  worst "
              f"{s['worst_abs_rel']:.4f}  signed median {s['median_signed_rel']:+.4f}"
              f"  ({s['n_negative']}/{s['n_residuals']} negative)")
        for ag, v in s["by_age"].items():
            print(f"    at {ag:>2s} Gyr: median |rel| {v['median_abs_rel']:.4f} "
                  f"(n = {v['n']})")
    print(f"\n  wrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
