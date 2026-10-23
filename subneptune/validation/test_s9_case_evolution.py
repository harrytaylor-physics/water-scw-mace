"""Cooling comparison at the strongest cases (A56).

Task 18 runs E2's history comparison at the scan cells the composite changes
most, plus the old fiducial on the new tables. Two things need gating: the
cooling-time-to-entropy arithmetic the headline number rests on (tested
synthetically, so it does not depend on a 30-minute run), and the convergence
discipline — A44 found a 24- vs 44-node age grid flipping the sign of a
sub-percent cooling difference, so a number that does not survive doubling is
not a result.
"""

from __future__ import annotations

import csv

import numpy as np
import pytest

import figures.plot_s9_case_evolution as s9


class _Track:
    """Minimal stand-in for CoolingTrack: entropy falling with age."""

    def __init__(self, age_gyr, entropy):
        self.age_gyr = np.asarray(age_gyr, float)
        self.entropy_j_kg_k = np.asarray(entropy, float)


class _Hist:
    def __init__(self, age, entropy_b, track_a, track_b):
        self.age_gyr = np.asarray(age, float)
        self.entropy_b = np.asarray(entropy_b, float)
        self.track_a = track_a
        self.track_b = track_b


def test_age_at_entropy_interpolates_a_falling_track():
    """Entropy decreases with age; the lookup must handle that orientation."""
    track = _Track([1.0, 2.0, 4.0], [3000.0, 2000.0, 1000.0])
    assert s9._age_at_entropy(track, 2000.0) == pytest.approx(2.0)
    assert s9._age_at_entropy(track, 2500.0) == pytest.approx(1.5)
    assert s9._age_at_entropy(track, 1500.0) == pytest.approx(3.0)
    # Outside the track's range -> NaN, never a silent clamp.
    assert np.isnan(s9._age_at_entropy(track, 5000.0))
    assert np.isnan(s9._age_at_entropy(track, 10.0))


def test_cooling_time_ratio_is_ages_to_a_common_entropy():
    """> 1 means the composite takes longer to reach the line list's entropy."""
    age = [1.0, 2.0, 4.0]
    slow = _Track(age, [3000.0, 2500.0, 2000.0])  # cools slowly
    fast = _Track(age, [3000.0, 2000.0, 1000.0])  # cools fast
    hist = _Hist(age, fast.entropy_j_kg_k, track_a=slow, track_b=fast)

    out = s9.cooling_time_ratio_at_entropy(hist, reference_age_gyr=2.0)
    assert out["reference_entropy_j_kg_k"] == pytest.approx(2000.0)
    assert out["age_hitemp_gyr"] == pytest.approx(2.0)
    assert out["age_composite_gyr"] == pytest.approx(4.0)
    assert out["cooling_time_ratio"] == pytest.approx(2.0)

    # Identical tracks -> exactly 1, whatever the reference age.
    same = s9.cooling_time_ratio_at_entropy(
        _Hist(age, fast.entropy_j_kg_k, track_a=fast, track_b=fast), 3.0)
    assert same["cooling_time_ratio"] == pytest.approx(1.0)


def test_the_reported_ages_are_the_ones_task_18_asks_for():
    assert s9.REPORT_AGES_GYR == (1.0, 5.0, 7.0, 10.0)
    assert s9.N_NODES >= 44  # "the converged evolution grid (44+ nodes)"
    assert s9.ESCAPE_WINDOW_GYR == 0.1
    assert s9.OLD_FIDUCIAL["t_eq_k"] == 500.0 and s9.OLD_FIDUCIAL["t_int_k"] == 20.0


@pytest.fixture(scope="module")
def s9_rows(cfg):
    path = cfg.output_path / "S9_case_evolution.csv"
    if not path.exists():
        pytest.skip("S9 not run; see figures/plot_s9_case_evolution.py")
    with path.open() as fh:
        return list(csv.DictReader(fh))


def test_every_case_was_run_at_both_resolutions(s9_rows):
    """The convergence check is not optional: each case appears twice."""
    by_case = {}
    for r in s9_rows:
        by_case.setdefault(r["case"], set()).add(int(r["n_nodes"]))
    assert by_case, "S9 produced no rows"
    for case, nodes in by_case.items():
        assert nodes == {s9.N_NODES, 2 * s9.N_NODES}, (
            f"{case} was not run at both {s9.N_NODES} and {2 * s9.N_NODES} nodes")
    assert any("old fiducial" in c for c in by_case), (
        "the old fiducial must be rerun on the new tables so the change is recorded")


def test_the_escape_window_is_still_a_clean_null(s9_rows):
    """A40's null, re-measured on the 48-point tables (A56).

    Atmospheric escape sculpts planets in the first 100 Myr, so any
    composite-vs-line-list difference there would be an irreversible imprint.
    A40 found |ΔR/R| ≤ 5.7e-5 % and |Δs/s| ≤ 1.3e-4 %; the enlarged grid and the
    new blend window must not have opened that up.
    """
    for r in s9_rows:
        dr = 100.0 * float(r["escape_max_abs_radius_frac"])
        ds = 100.0 * float(r["escape_max_abs_entropy_frac"])
        assert dr < 1.0e-2, (
            f"{r['case']}: escape-window |ΔR/R| rose to {dr:.3e} % — the null that "
            f"the mass-loss story rests on has opened up; re-report it")
        assert ds < 1.0e-1, f"{r['case']}: escape-window |Δs/s| rose to {ds:.3e} %"
