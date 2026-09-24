"""FP/HITEMP crossover recomputed on the merged 48-point grid (A53).

``opacity/crossover.py`` replaces the delivered ``crossover_hitemp.json`` as the
*source* of the M2 ratios (it can reach the 21 extension points, which RADIS
never covered) while keeping it as the independent check. These gates pin the
agreement with that check, the provenance flags that keep a clamped or
floor-limited point visible, and the spectral-contrast diagnostic that says what
the cold Rosseland ratios are actually measuring.
"""

from __future__ import annotations

import json

import pytest

from opacity.crossover import agreement_summary, crossover_rows, write_json

BAND_CM = (100.0, 4000.0)
N_POINTS = 48
N_CLAMPED = 9  # the nine legacy points with rho > 1000, HITEMP's grid ceiling
N_DILUTE = 9


@pytest.fixture(scope="module")
def rows(cfg):
    return crossover_rows(cfg, BAND_CM)


def test_one_row_per_fp_grid_point(rows):
    assert len(rows) == N_POINTS
    assert all(r["M2_planck_ratio"] > 0 and r["M2_rosseland_ratio"] > 0 for r in rows)
    # Sorted by (T, rho), as the printed table assumes.
    keys = [(r["T"], r["rho"]) for r in rows]
    assert keys == sorted(keys)


def test_hitemp_is_clamped_exactly_where_its_grid_ends(rows):
    """Only rho > 1000 clamps; the whole thin and cold part of the grid is interior.

    The last assertion used to be phrased over the points added by a later
    delivery. The 48 points are one grid now (A62), so it is phrased over where
    the points sit instead — which is what it always meant: nothing dilute and
    nothing cold is outside the line list's range.
    """
    clamped = [r for r in rows if r["ll_clamped"]]
    assert len(clamped) == N_CLAMPED
    assert all(r["rho"] > 1000.0 for r in clamped)
    assert all(not r["ll_clamped"] for r in rows if r["rho"] <= 1000.0)
    assert all(not r["ll_clamped"] for r in rows if r["is_dilute"] or r["T"] <= 1000.0)


def test_agreement_with_the_delivered_analysis(rows):
    """The real cross-check: where both coordinates are on HITEMP's own grid.

    Those points involve no interpolation on the line-list side, so any
    disagreement would be ours. A25 claimed ~1%; this holds it.
    """
    summary = agreement_summary(rows)
    on_grid = summary["on_grid"]
    assert on_grid["n"] >= 3
    assert on_grid["max_abs_dev_planck_pct"] < 1.0
    assert on_grid["max_abs_dev_rosseland_pct"] < 1.0
    # T-interpolated points. This was < 5% until A60 put line-list rows at 2250 and
    # 2750 K, the only two temperatures that needed interpolating; now it is < 1%,
    # and the bound is tightened rather than left where it was.
    assert summary["interpolated"]["max_abs_dev_rosseland_pct"] < 1.0
    # Clamped ones still disagree — RADIS had the true density and we hold kappa at
    # the ceiling — but by the amount A59 MEASURED against a line-by-line run at
    # that density (5.89% at the most extreme point), not by the ~30% A53 reported.
    # That 30% was the two T-interpolated points, and A60 removed it. The assertion
    # used to require > 10%, i.e. it pinned the pathology as expected behaviour.
    clamped = summary["clamped"]["max_abs_dev_rosseland_pct"]
    assert 1.0 < clamped < 6.5


def test_dilute_rows_carry_their_noise_floor_fractions(cfg, rows):
    """A ratio the A52 floor makes unreliable must be visible as one."""
    dilute = [r for r in rows if r["is_dilute"]]
    assert len(dilute) == N_DILUTE
    assert all(r["rho"] <= 60.0 for r in dilute)
    if not (cfg.output_path / "S7_dilute_floor.csv").exists():
        pytest.skip("S7 has not been run; floor columns are optional by design")
    for r in dilute:
        assert r["floor_frac_above_window"] < 0.10
        assert r["floor_rosseland_window_share"] > 0.8


def test_spectral_contrast_explains_the_cold_ratios(rows):
    """The cold Rosseland ratios are a window-depth difference, not opacity.

    κ_P/κ_R measures how deep a band's transparency windows are. HITEMP's rises
    by two orders of magnitude as it cools (sparse, narrow lines); FP's does not
    move, because its windows are filled in. That, not a 100× opacity
    difference, is what the T = 650 K ratios are reporting (A52/A53).
    """
    cold = [r for r in rows if r["T"] < 1000.0]
    hot = [r for r in rows if r["T"] >= 2500.0]
    assert cold and hot
    # FP has no strong window structure anywhere on the grid: its contrast spans
    # 1-33 over all 48 points, against HITEMP's 2.4-567.
    assert max(r["contrast_fp"] for r in rows) < 40.0
    # HITEMP's rises sharply when cold, and FP's does not follow.
    assert min(r["contrast_hitemp"] for r in cold) > 100.0
    assert max(r["contrast_hitemp"] for r in hot) < 10.0
    assert all(r["contrast_hitemp"] / r["contrast_fp"] > 15.0 for r in cold)
    # So the cold Rosseland ratios are huge while the Planck ones stay ~1-2.
    assert min(r["M2_rosseland_ratio"] for r in cold) > 20.0
    assert max(r["M2_planck_ratio"] for r in cold) < 3.0


def test_v2_json_round_trips(cfg, rows, tmp_path):
    """The written artefact carries the band, the flags and every point."""
    path = write_json(cfg, rows, BAND_CM)
    payload = json.loads(path.read_text())
    assert payload["n_points"] == N_POINTS
    assert payload["band_cm"] == list(BAND_CM)
    assert payload["line_list"] == "hitemp"
    assert len(payload["results"]) == N_POINTS
    assert sum(r["ll_clamped"] for r in payload["results"]) == N_CLAMPED
    assert "interpolation" in payload["meta"]
