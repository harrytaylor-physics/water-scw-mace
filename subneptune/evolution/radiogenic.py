"""Radiogenic T_int floor for an inert rocky core (A30).

v1 omits core luminosity (SPEC §4), so a model's T_int can fall arbitrarily low.
In reality radioactive decay in an Earth-like core supplies a floor luminosity,
setting a minimum internal temperature the planet cools toward but not below.
This estimator quantifies that floor as a stated sensitivity (not used inside
the structure/evolution solvers — an annotation for figures and analysis).

Radiogenic power: an Earth-like core of mass M_c produces ≈ M_c × H_EARTH_NOW,
with H_EARTH_NOW ≈ 20 TW the present-day (4.5 Gyr) terrestrial radiogenic power
(BSE/geoneutrino estimates 15–25 TW). Decay follows the U/Th/K ensemble; a
single effective e-folding reproduces the ~0.5× decline by 10 Gyr adequately for
a floor whose T ∝ L^¼.
"""

from __future__ import annotations

import numpy as np

from constants import SIGMA_SB

H_EARTH_NOW_W = 20.0e12  # present-day Earth radiogenic power [W] (~15-25 TW)
_AGE_NOW_GYR = 4.5  # age at which H_EARTH_NOW applies
_EFFECTIVE_HALFLIFE_GYR = 6.0  # ensemble effective (U/Th/K); ~0.5× per 6 Gyr


def radiogenic_power_w(mass_core_mearth: float, age_gyr: float = _AGE_NOW_GYR) -> float:
    """Radiogenic luminosity of an Earth-like core [W] at ``age_gyr``."""
    decay = 0.5 ** ((age_gyr - _AGE_NOW_GYR) / _EFFECTIVE_HALFLIFE_GYR)
    return mass_core_mearth * H_EARTH_NOW_W * decay


def t_int_floor_k(mass_core_mearth: float, radius_p_m: float, age_gyr: float = _AGE_NOW_GYR) -> float:
    """Equilibrium T_int where L_int = 4πR_p²σT_int⁴ balances radiogenic power.

    The planet cools toward this floor, not below it, so the FP-relevant regime
    (RCB at ρ≥109, requiring T_int ≲ 25 K) is truncated when the floor exceeds it.
    """
    lum = radiogenic_power_w(mass_core_mearth, age_gyr)
    return float((lum / (4.0 * np.pi * radius_p_m**2 * SIGMA_SB)) ** 0.25)
