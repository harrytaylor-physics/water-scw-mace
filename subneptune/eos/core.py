"""Inert rocky core mass–radius relation (SPEC §4).

Zeng, Sasselov & Jacobsen (2016), ApJ 819, 127: for a differentiated
rock+iron planet the radius follows the power law

    R / R⊕ = (1.07 − 0.21 · CMF) · (M / M⊕)^0.27

with CMF the core (iron) mass fraction. Earth-like composition is CMF = 0.325
(32.5% Fe core + 67.5% MgSiO₃ mantle) — the default here (AMENDMENTS A9). The
fit is accurate to a few percent for M ≲ 8 M⊕; no core luminosity in v1.
"""

from __future__ import annotations

from constants import M_EARTH, R_EARTH

# Zeng (2016) Earth-like composition.
EARTH_LIKE_CMF = 0.325
_EXPONENT = 0.27


def core_radius_rearth(mass_mearth: float, cmf: float = EARTH_LIKE_CMF) -> float:
    """Rocky-core radius in Earth radii (Zeng 2016)."""
    if mass_mearth <= 0:
        raise ValueError("core mass must be positive")
    if not 0.0 <= cmf <= 1.0:
        raise ValueError("core mass fraction must be in [0, 1]")
    return (1.07 - 0.21 * cmf) * mass_mearth**_EXPONENT


def core_radius_si(mass_mearth: float, cmf: float = EARTH_LIKE_CMF) -> float:
    """Rocky-core radius in metres (Zeng 2016)."""
    return core_radius_rearth(mass_mearth, cmf) * R_EARTH


def core_mass_si(mass_mearth: float) -> float:
    """Core mass in kilograms."""
    return mass_mearth * M_EARTH
