"""Radiogenic T_int floor at the candidate core (A45).

The candidate fiducial's structural case sits at T_int = 20 K. Whether a real
planet can be there is decided by the radiogenic floor (A30): an Earth-like core
supplies a luminosity the planet cannot cool below. These gates pin the
estimator's scaling and the statement the manuscript needs — that for a
5 M⊕ Earth-like core the floor is above 20 K over the whole plausible radius
range, at every age from 1 to 10 Gyr.
"""

from __future__ import annotations

import numpy as np
import pytest

from constants import R_EARTH
from evolution.radiogenic import H_EARTH_NOW_W, radiogenic_power_w, t_int_floor_k

CANDIDATE_T_INT_K = 20.0
M_CORE = 5.0
# The candidate's R_p is 1.759 R⊕; the floor falls as R_p grows, so bracketing
# generously keeps the gate about physics rather than about one converged model.
R_P_RANGE_REARTH = (1.5, 2.5)


def test_present_day_power_is_earth_like_per_core_mass():
    assert radiogenic_power_w(1.0, 4.5) == pytest.approx(H_EARTH_NOW_W)
    assert radiogenic_power_w(M_CORE, 4.5) == pytest.approx(M_CORE * H_EARTH_NOW_W)


def test_power_decays_and_the_floor_follows_l_quarter():
    """L declines with age; T_int floor ∝ L^¼ at fixed radius."""
    young = radiogenic_power_w(M_CORE, 1.0)
    now = radiogenic_power_w(M_CORE, 4.5)
    old = radiogenic_power_w(M_CORE, 10.0)
    assert young > now > old

    r_m = 1.76 * R_EARTH
    assert t_int_floor_k(M_CORE, r_m, 1.0) / t_int_floor_k(M_CORE, r_m, 4.5) == (
        pytest.approx((young / now) ** 0.25))
    # Doubling the core mass doubles L and raises the floor by 2^¼.
    assert t_int_floor_k(2 * M_CORE, r_m, 4.5) / t_int_floor_k(M_CORE, r_m, 4.5) == (
        pytest.approx(2.0**0.25))


@pytest.mark.parametrize("age_gyr", [1.0, 4.5, 10.0])
def test_candidate_t_int_is_below_the_floor(age_gyr):
    """A45: T_int = 20 K is radiogenically unreachable for a 5 M⊕ core."""
    for r_rearth in R_P_RANGE_REARTH:
        floor = t_int_floor_k(M_CORE, r_rearth * R_EARTH, age_gyr)
        assert floor > CANDIDATE_T_INT_K, (
            f"at {age_gyr} Gyr and R_p = {r_rearth} R⊕ the floor is {floor:.2f} K, no "
            f"longer above the candidate's T_int = {CANDIDATE_T_INT_K:g} K"
        )


def test_floor_falls_with_radius():
    """A bigger planet radiates the same power at a lower T_int."""
    floors = [t_int_floor_k(M_CORE, r * R_EARTH, 4.5) for r in (1.5, 2.0, 2.5)]
    assert np.all(np.diff(floors) < 0)
