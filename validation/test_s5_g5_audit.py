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

The corrections are opt-in `anchor.*` flags; these gates pin both the arithmetic
and the fact that the repository defaults still reproduce the published numbers.
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest

from config import Anchor, load_config
from config.schema import _validate
from eos import AquaEOS
from validation.test_g5_literature import anchor_t_eq, our_radius

# GJ 9827 (K2-135), the system in their Table 3. Stellar and orbital values from
# the NASA Exoplanet Archive (the table's stated source; Niraula et al. 2017,
# Rodriguez et al. 2018): T_eff = 4340 K, R⋆ = 0.602 R_sun.
_T_EFF_K = 4340.0
_R_STAR_RSUN = 0.602
_R_SUN_AU = 0.00465047
_SEMI_MAJOR_AU = {"b": 0.0189, "c": 0.0394, "d": 0.0562}
# Their Table 3 T_irr, read from the paper.
_THEIR_T_IRR_K = {"b": 1184.0, "c": 820.0, "d": 686.0}


def _t_irr_eq9(a_au: float, include_half: bool = True) -> float:
    """T_eff·√(R⋆/2a) — their Eq. (9). ``include_half=False`` drops the factor 2."""
    denom = (2.0 if include_half else 1.0) * a_au
    return _T_EFF_K * np.sqrt(_R_STAR_RSUN * _R_SUN_AU / denom)


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


def test_anchor_conversion_matches_the_convention(cfg):
    """``sqrt2`` (default) divides; ``identity`` passes their T_irr through."""
    assert anchor_t_eq(cfg, 1000.0) == pytest.approx(1000.0 / np.sqrt(2.0))
    corrected = replace(cfg, anchor=replace(cfg.anchor, t_irr_to_t_eq="identity"))
    assert anchor_t_eq(corrected, 1000.0) == pytest.approx(1000.0)
    # The correction makes the modelled planet hotter by exactly √2.
    assert anchor_t_eq(corrected, 1000.0) / anchor_t_eq(cfg, 1000.0) == pytest.approx(
        np.sqrt(2.0))


def test_repository_defaults_are_unchanged(cfg, repo_root):
    """The audit must not move the published G5 numbers (A32)."""
    assert cfg.anchor == Anchor()
    assert cfg.anchor.t_irr_to_t_eq == "sqrt2"
    assert cfg.anchor.radius_pressure_bar is None
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
    base = our_radius(cfg, **args)
    hotter = our_radius(replace(cfg, anchor=replace(cfg.anchor, t_irr_to_t_eq="identity")),
                        **args)
    higher = our_radius(replace(cfg, anchor=replace(cfg.anchor, radius_pressure_bar=1e-6)),
                        **args)
    assert hotter > base
    assert higher > base
    # The radius level is the smaller effect: one decade of pressure is a couple
    # of scale heights, the √2 in T_eq is a different planet.
    assert (higher - base) < (hotter - base)
