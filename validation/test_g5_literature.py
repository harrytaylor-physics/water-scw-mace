"""Gate G5: static M–R for steam worlds vs Aguichine et al. (2021) (SPEC §7).

Compares our R_p at matched (M, f_env, T_eq) to the published steam-world M–R,
requiring agreement to ~5–10%, and writes a residual report. The Aguichine
reference table is not on the machine (offline), so this gate SKIPS until
`validation/data/aguichine2021.csv` is provided (see that directory's README).
Residuals are documented, not tuned away (SPEC §7 G5).
"""

from __future__ import annotations

import csv
import json
from dataclasses import replace
from pathlib import Path

import pytest

from eos import AquaEOS
from structure import build_structure

_REF = Path(__file__).parent / "data" / "aguichine2021.csv"
_TOL = 0.10  # 5–10% target; assert at the looser 10% bound, report all residuals


def _load_reference():
    if not _REF.exists():
        pytest.skip(
            f"G5 reference {_REF.name} not present (offline); provide it to enable "
            f"the Aguichine 2021 comparison — see validation/data/README.md"
        )
    with _REF.open() as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        pytest.skip(f"{_REF.name} is empty")
    return rows


def our_radius(cfg, m_earth, f_env, t_eq, t_int, source, eos):
    """R_p [R⊕] at matched conditions (T_int is the sequence label; use a low value)."""
    cfg2 = replace(
        cfg,
        core=replace(cfg.core, mass_mearth=m_earth),
        envelope=replace(cfg.envelope, f_env=f_env),
    )
    return build_structure(cfg2, t_eq, t_int, source, eos=eos).radius_p_rearth


def test_g5_steam_world_mass_radius(cfg):
    rows = _load_reference()
    eos = AquaEOS.from_config(cfg)
    source = cfg.opacity.source
    t_int = 50.0  # cool, near-static; documented choice

    residuals = []
    for r in rows:
        m, f, teq, r_ref = (
            float(r["m_earth"]), float(r["f_env"]), float(r["t_eq_k"]), float(r["r_earth"])
        )
        r_ours = our_radius(cfg, m, f, teq, t_int, source, eos)
        rel = (r_ours - r_ref) / r_ref
        residuals.append({"m_earth": m, "f_env": f, "t_eq_k": teq,
                          "r_ref": r_ref, "r_ours": r_ours, "rel": rel})

    report = cfg.output_path / "g5_residuals.json"
    report.parent.mkdir(parents=True, exist_ok=True)
    report.write_text(json.dumps({"source": source, "t_int_k": t_int, "residuals": residuals}, indent=2))

    worst = max(abs(x["rel"]) for x in residuals)
    assert worst <= _TOL, (
        f"G5: worst |ΔR_p/R_p| = {worst:.1%} exceeds {_TOL:.0%} (source={source}); "
        f"see {report}"
    )
