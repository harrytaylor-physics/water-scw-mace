"""G5 definition audit: the anchor conventions and their evidence (A43).

Gate G5 compares our static M–R against Aguichine et al. (2021). Three
definitions have to agree before any residual can be called physics. The one
that was wrong is the irradiation temperature, and the evidence is arithmetic
the paper itself supplies: their Equation (9) irradiance temperature,
T_irr = T_eff·√(R⋆/2a), reproduces the T_irr values of their own Table 3 for
GJ 9827 b/c/d. Dropping the factor 2 does not — it is 41% high. So their T_irr
*is* the zero-albedo, full-redistribution equilibrium temperature our T_eq means,
and dividing it by √2 (the original A32 assumption) modelled planets one √2 too
cold.

The corrections were opt-in `anchor.*` flags in A43 and became the repository
defaults in A46; these gates pin both the arithmetic and the fact that the
defaults are now the corrected conventions, with the pre-A46 ones still
selectable so the audit can attribute each step.
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest

from config import Anchor, load_config
from config.schema import _validate
from eos import AquaEOS
from figures.plot_s5_g5_audit import PRE_A46
from validation.test_g5_literature import anchor_t_eq, our_radius

# GJ 9827 (K2-135), the system in their Table 3, whose planetary data the paper
# takes from the NASA Exoplanet Archive. Two published stellar parameter sets are
# used below so the conclusion does not hinge on one stellar radius:
#   Rodriguez et al. (2018): T_eff = 4340 K, R⋆ = 0.602 R_sun
#   Niraula et al. (2017):   T_eff = 4255 K, R⋆ = 0.651 R_sun
_STELLAR = {"rodriguez2018": (4340.0, 0.602), "niraula2017": (4255.0, 0.651)}
_R_SUN_AU = 0.00465047
_SEMI_MAJOR_AU = {"b": 0.0189, "c": 0.0394, "d": 0.0562}
# Their Table 3 T_irr, read from the paper.
_THEIR_T_IRR_K = {"b": 1184.0, "c": 820.0, "d": 686.0}


def _t_irr_eq9(a_au: float, stellar: str = "rodriguez2018", include_half: bool = True) -> float:
    """T_eff·√(R⋆/2a) — their Eq. (9). ``include_half=False`` drops the factor 2."""
    t_eff, r_star = _STELLAR[stellar]
    denom = (2.0 if include_half else 1.0) * a_au
    return t_eff * np.sqrt(r_star * _R_SUN_AU / denom)


@pytest.mark.parametrize("planet", sorted(_THEIR_T_IRR_K))
def test_their_t_irr_is_the_equilibrium_temperature(planet):
    """A43: Eq. (9) with the factor 2 reproduces their own Table 3 to <1%."""
    predicted = _t_irr_eq9(_SEMI_MAJOR_AU[planet])
    theirs = _THEIR_T_IRR_K[planet]
    assert predicted == pytest.approx(theirs, rel=0.01), (
        f"GJ 9827 {planet}: Eq. (9) gives {predicted:.0f} K vs their {theirs:.0f} K"
    )
    # Without the factor 2 it is ~41% high, i.e. not what they tabulate.
    assert _t_irr_eq9(_SEMI_MAJOR_AU[planet], include_half=False) > 1.3 * theirs


@pytest.mark.parametrize("planet", sorted(_THEIR_T_IRR_K))
def test_the_factor_2_survives_the_stellar_parameter_choice(planet):
    """The √2 conclusion is not an artifact of one published R⋆/T_eff pair."""
    theirs = _THEIR_T_IRR_K[planet]
    alt = _t_irr_eq9(_SEMI_MAJOR_AU[planet], stellar="niraula2017")
    assert alt == pytest.approx(theirs, rel=0.03)  # ~1.7% high, still not 41%
    assert _t_irr_eq9(_SEMI_MAJOR_AU[planet], stellar="niraula2017",
                      include_half=False) > 1.3 * theirs


def test_anchor_conversion_matches_the_convention(cfg):
    """``identity`` (default since A46) passes their T_irr through; ``sqrt2`` divides."""
    assert anchor_t_eq(cfg, 1000.0) == pytest.approx(1000.0)
    old = replace(cfg, anchor=replace(cfg.anchor, **PRE_A46))
    assert anchor_t_eq(old, 1000.0) == pytest.approx(1000.0 / np.sqrt(2.0))
    # The correction makes the modelled planet hotter by exactly √2.
    assert anchor_t_eq(cfg, 1000.0) / anchor_t_eq(old, 1000.0) == pytest.approx(
        np.sqrt(2.0))


def test_repository_defaults_are_the_corrected_conventions(cfg, repo_root):
    """A46: G5 now runs like-for-like by default."""
    assert cfg.anchor == Anchor()
    assert cfg.anchor.t_irr_to_t_eq == "identity"
    assert cfg.anchor.radius_pressure_bar == pytest.approx(1.0e-6)  # 0.1 Pa
    assert cfg.anchor.cmf == 0.0
    # The candidate config inherits the same anchor defaults.
    assert load_config(repo_root / "config" / "candidate.yaml").anchor == Anchor()


def test_anchor_validation_rejects_nonsense(cfg):
    for bad in ({"t_irr_to_t_eq": "sqrt3"}, {"radius_pressure_bar": 0.0},
                {"radius_pressure_bar": -1.0}, {"cmf": 1.0}):
        with pytest.raises(ValueError):
            _validate(replace(cfg, anchor=replace(cfg.anchor, **bad)))


def test_corrections_move_the_radius_in_the_expected_direction(cfg):
    """Hotter irradiation and a shallower radius level both grow R_p."""
    eos = AquaEOS.from_config(cfg)
    args = dict(m_earth=5.0, wmf=0.2, t_irr_k=1000.0, t_int_k=30.0, source="hitemp", eos=eos)
    old_cfg = replace(cfg, anchor=replace(cfg.anchor, **PRE_A46))
    base = our_radius(old_cfg, **args)
    hotter = our_radius(replace(old_cfg, anchor=replace(old_cfg.anchor,
                                                        t_irr_to_t_eq="identity")), **args)
    higher = our_radius(replace(old_cfg, anchor=replace(old_cfg.anchor,
                                                        radius_pressure_bar=1e-6)), **args)
    assert hotter > base
    assert higher > base
    # The radius level is the smaller effect: one decade of pressure is a couple
    # of scale heights, the √2 in T_eq is a different planet.
    assert (higher - base) < (hotter - base)
