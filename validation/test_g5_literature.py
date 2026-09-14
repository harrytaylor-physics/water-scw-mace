"""Gate G5: static M–R for steam worlds vs Aguichine et al. (2021) (SPEC §7).

Anchor = the published mass–radius relation of Aguichine, Mousis, Deleuil &
Marcq (2021), ApJ 914, 84. Their parametric fit is Eq. 29,

    log10 R_p = a·log10(M_p) + b + exp(−d·(log10(M_p) + c)),

with a,b,c,d tabulated per (composition, T_irr) in an ONLINE machine-readable
table not contained in the PDF (AMENDMENTS A32). We therefore anchor to points
digitized from their Figure 8 (Ma19+, no metallic core), stored in
`validation/data/aguichine2021_fig8_points.csv`.

We match (M, WMF, T_irr) under the corrected like-for-like conventions of A43,
adopted as defaults in A46: WMF → f_env (pure-steam envelope); **T_eq = T_irr**,
because their Equation (9) T_irr = T_eff·√(R⋆/2a) already *is* the
full-redistribution equilibrium temperature (`anchor.t_irr_to_t_eq = "identity"`);
the radius read at their **0.1 Pa** transiting level (`anchor.radius_pressure_bar
= 1e-6` bar); a pure-rock core (cmf = 0) for their no-metallic-core case; and a
cool mature interior (T_int = 30 K).

Under those conventions the **median |ΔR_p/R_p| is 0.084** (it was 0.165 with the
pre-A46 T_eq = T_irr/√2 and 1 Pa radius: about half the old residual was a
definition mismatch, not physics). Per SPEC §7 G5 we DOCUMENT the remainder
rather than tuning it away, and three candidate explanations have been measured:

  * **The reference model's layered supercritical-water structure.** Every
    corrected residual is still negative and grows with WMF (−3.9 to −9.9% at
    WMF = 0.2, −6.9 to −17.9% at WMF = 0.5), which is where our single pure-steam
    envelope departs from their ice/liquid/supercritical water layers. This is
    the attribution.
  * **Matching pressure contributes 0.35%.** Scanning the interior matching cap
    over 1–100 bar moves R_p by 0.349% at the fiducial (A43) — two orders of
    magnitude too small.
  * **Band truncation does not reduce it.** Recomputing these residuals with
    full-band HITEMP κ_R instead of the [100, 4000] cm⁻¹ science band moves the
    median the wrong way (A40, re-measured under the corrected conventions in
    A46), so the science band is not hiding the offset either.
"""

from __future__ import annotations

import csv
import json
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from eos import AquaEOS
from structure import build_structure

_REF = Path(__file__).parent / "data" / "aguichine2021_fig8_points.csv"
# Order-of-magnitude agreement bound. Under the corrected conventions (A46) the
# low-WMF (steam-envelope-like) points sit at 3.9–9.9%, inside SPEC's 5–10% ideal;
# the high-WMF points reach ~18% from the documented model-definition differences.
# The gate confirms structural agreement; the residual report (below) carries the
# physics.
_TOL = 0.30
# Median |ΔR_p/R_p| under the A46 default conventions, pinned to 1% so a silent
# change of definition, EOS or opacity default cannot pass unnoticed.
_MEDIAN_A46 = 0.0844


def _load_points():
    if not _REF.exists():
        pytest.skip(f"G5 reference {_REF.name} not present")
    with _REF.open() as fh:
        rows = [r for r in csv.DictReader(l for l in fh if not l.startswith("#"))]
    return rows


def anchor_t_eq(cfg, t_irr_k: float) -> float:
    """The anchor's T_irr expressed as our T_eq, per ``cfg.anchor`` (A43/A46).

    ``identity`` (the default since A46) uses their Equation (9) definition,
    T_irr = T_eff·√(R⋆/2a), which *is* the full-redistribution equilibrium
    temperature our T_eq means; ``sqrt2`` is the pre-A46 convention, which
    treated their T_irr as Guillot's irradiation temperature.
    """
    return float(t_irr_k) if cfg.anchor.t_irr_to_t_eq == "identity" else float(t_irr_k) / np.sqrt(2.0)


def our_radius(cfg, m_earth, wmf, t_irr_k, t_int_k, source, eos):
    """R_p [R⊕] for a steam envelope (f_env=WMF) on a rocky core, matched to
    Aguichine (M, WMF, T_irr) under the ``cfg.anchor`` conventions (A43/A46).

    Defaults since A46: T_eq = T_irr, cmf = 0 (no metallic core), radius read at
    their 0.1 Pa transiting level.
    """
    boundary = cfg.boundary
    if cfg.anchor.radius_pressure_bar is not None:
        boundary = replace(boundary, p_top_bar=cfg.anchor.radius_pressure_bar)
    cfg2 = replace(
        cfg,
        core=replace(cfg.core, mass_mearth=(1.0 - wmf) * m_earth, cmf=cfg.anchor.cmf),
        envelope=replace(cfg.envelope, f_env=wmf),
        boundary=boundary,
    )
    return build_structure(cfg2, anchor_t_eq(cfg, t_irr_k), t_int_k, source, eos=eos,
                           n_sample=400).radius_p_rearth


def test_g5_steam_world_mass_radius(cfg):
    rows = _load_points()
    eos = AquaEOS.from_config(cfg)
    source = "hitemp"  # full radiative-zone coverage; R is EOS-dominated here
    t_int = 30.0

    residuals = []
    for r in rows:
        wmf, t_irr, m, r_ref = (
            float(r["wmf"]), float(r["t_irr_k"]), float(r["m_earth"]), float(r["r_earth"])
        )
        r_ours = our_radius(cfg, m, wmf, t_irr, t_int, source, eos)
        rel = (r_ours - r_ref) / r_ref
        residuals.append({"wmf": wmf, "t_irr_k": t_irr, "m_earth": m,
                          "r_aguichine": r_ref, "r_ours": r_ours, "rel": rel})

    median = float(np.median([abs(x["rel"]) for x in residuals]))
    report = cfg.output_path / "g5_residuals.json"
    report.parent.mkdir(parents=True, exist_ok=True)
    report.write_text(json.dumps({
        "anchor": "Aguichine et al. 2021 ApJ 914 84, Fig 8 (digitized); fit Eq. 29",
        "source": source, "t_int_k": t_int, "cmf": cfg.anchor.cmf,
        "conventions": {
            "t_irr_to_t_eq": cfg.anchor.t_irr_to_t_eq,
            "radius_pressure_bar": (cfg.anchor.radius_pressure_bar
                                    if cfg.anchor.radius_pressure_bar is not None
                                    else cfg.boundary.p_top_bar),
            "note": "A46 defaults: their Eq. (9) T_irr is our T_eq; radius at their "
                    "0.1 Pa transiting level.",
        },
        "median_abs_rel": median,
        "residuals": residuals,
    }, indent=2))

    worst = max(abs(x["rel"]) for x in residuals)
    assert worst <= _TOL, f"G5 worst |ΔR_p/R_p| = {worst:.1%} > {_TOL:.0%}; see {report}"
    # A46: the corrected conventions are the defaults, so the published median is
    # now a gate rather than a reported number.
    assert median == pytest.approx(_MEDIAN_A46, rel=0.01), (
        f"G5 median |ΔR_p/R_p| = {median:.4f}, not the A46 value {_MEDIAN_A46:.4f} "
        f"(1% tolerance); see {report}"
    )
