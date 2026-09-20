"""Gates for the reported history/cooling-range diagnostics (A40).

Task-4 numbers the manuscript quotes come from three small pieces of arithmetic
that are easy to get subtly wrong: the escape-window extrema (must include the
window edge, must not run past it), the age at which a track's deepest RCB first
reaches the FP grid floor (interpolated, and ``None`` when never reached), and
the min/max of the F3 τ_KH ratio over the scanned (f_env, T_eq) grid. These are
checked against hand-computable inputs.

The *physical* claim behind the escape-window numbers — that the composite and
line-list histories agree inside the first 100 Myr — is already gated on real
models by ``test_evolution.test_history_comparison_diagnostics``; repeating a
full cooling comparison here would only duplicate that cost.
"""

from __future__ import annotations

import numpy as np
import pytest

from evolution.history import HistoryComparison
from figures.plot_f3_cooling_ratio import F_ENV_GRID, ratio_extremes
from structure import FP_GRID_FLOOR_SI

# The legacy FP ρ floor. Used here as an arbitrary interior threshold so the
# interpolation test stays meaningful whatever the live floor is (A54).
THRESHOLD_SI = 108.796


def _history(age, radius_a, radius_b, entropy_a, entropy_b, rcb_a, rcb_b):
    return HistoryComparison(
        source_a="composite", source_b="hitemp", t_eq_k=500.0, f_env=0.05,
        age_gyr=np.asarray(age, float), t_int_k=np.full(len(age), 100.0),
        radius_a=np.asarray(radius_a, float), radius_b=np.asarray(radius_b, float),
        rcb_rho_a=np.asarray(rcb_a, float), rcb_rho_b=np.asarray(rcb_b, float),
        entropy_a=np.asarray(entropy_a, float), entropy_b=np.asarray(entropy_b, float),
    )


def test_escape_window_extrema_are_confined_to_the_window():
    """Values after the window must not leak into the reported maxima."""
    age = [0.01, 0.05, 0.1, 1.0, 10.0]
    h = _history(age,
                 radius_a=[1.01, 1.02, 1.00, 2.00, 3.00],  # the 2.0/3.0 entries are late
                 radius_b=[1.00, 1.00, 1.00, 1.00, 1.00],
                 entropy_a=[100.0, 100.0, 105.0, 500.0, 900.0],
                 entropy_b=[100.0, 100.0, 100.0, 100.0, 100.0],
                 rcb_a=[1.0] * 5, rcb_b=[1.0] * 5)
    esc = h.escape_window_extrema(0.1)
    assert esc["n_grid_points"] == 3
    assert esc["max_abs_radius_frac"] == pytest.approx(0.02)  # at 0.05 Gyr, not 1.0 Gyr
    assert esc["age_at_max_radius_gyr"] == pytest.approx(0.05)
    assert esc["max_abs_entropy_frac"] == pytest.approx(0.05)  # at the window edge
    assert esc["age_at_max_entropy_gyr"] == pytest.approx(0.1)
    assert esc["age_min_gyr"] == pytest.approx(0.01)


def test_escape_window_includes_the_interpolated_edge():
    """A sparse grid must not miss the window edge (the value is interpolated)."""
    h = _history([0.01, 1.0],
                 radius_a=[1.00, 2.00], radius_b=[1.00, 1.00],
                 entropy_a=[100.0, 100.0], entropy_b=[100.0, 100.0],
                 rcb_a=[1.0, 1.0], rcb_b=[1.0, 1.0])
    esc = h.escape_window_extrema(0.1)
    assert esc["n_grid_points"] == 1  # only the 0.01 Gyr node is inside
    edge = np.interp(0.1, [0.01, 1.0], [0.0, 1.0])  # linear in age, as diff_at_age is
    assert esc["max_abs_radius_frac"] == pytest.approx(edge)
    assert esc["age_at_max_radius_gyr"] == pytest.approx(0.1)


def test_first_age_rcb_above_interpolates_the_crossing():
    """The crossing age is interpolated in log age across the bracket."""
    age = [1.0, 10.0, 100.0]
    h = _history(age,
                 radius_a=[1.0] * 3, radius_b=[1.0] * 3,
                 entropy_a=[1.0] * 3, entropy_b=[1.0] * 3,
                 rcb_a=[10.0, 50.0, 200.0],  # crosses THRESHOLD between 10 and 100 Gyr
                 rcb_b=[10.0, 10.0, 10.0])  # never crosses
    # An explicit threshold, not the module default: A54 moved FP_GRID_FLOOR_SI to
    # 10, which this synthetic track starts *at*, so the default would make the
    # crossing trivial and stop testing the interpolation.
    cross = h.first_age_rcb_above(THRESHOLD_SI)
    w = (THRESHOLD_SI - 50.0) / (200.0 - 50.0)
    expected = np.exp(np.log(10.0) + w * (np.log(100.0) - np.log(10.0)))
    assert cross["composite"] == pytest.approx(expected)
    assert cross["hitemp"] is None
    assert cross["rho_si"] == pytest.approx(THRESHOLD_SI)


def test_first_age_rcb_above_handles_an_already_deep_start():
    h = _history([0.5, 1.0],
                 radius_a=[1.0, 1.0], radius_b=[1.0, 1.0],
                 entropy_a=[1.0, 1.0], entropy_b=[1.0, 1.0],
                 rcb_a=[500.0, 600.0], rcb_b=[1.0, 2.0])
    assert h.first_age_rcb_above(THRESHOLD_SI)["composite"] == pytest.approx(0.5)


def test_f3_ratio_extremes_report_their_grid_location():
    """The quoted F3 range must carry the (T_eq, f_env) it came from."""
    n = len(F_ENV_GRID)
    grids = {
        500.0: {"composite": [1.0 + 0.1 * i for i in range(n)]},
        800.0: {"composite": [0.5 + 0.1 * i for i in range(n)]},
    }
    lo, hi = ratio_extremes(grids, "composite", [500.0, 800.0])
    assert lo == (0.5, 800.0, F_ENV_GRID[0])
    assert hi == (pytest.approx(1.0 + 0.1 * (n - 1)), 500.0, F_ENV_GRID[-1])

    # A single-T_eq scan (the default F3 run) reports that T_eq only.
    lo1, hi1 = ratio_extremes(grids, "composite", [500.0])
    assert lo1[1] == 500.0 and hi1[1] == 500.0
