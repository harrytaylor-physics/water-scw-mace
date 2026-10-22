"""Clamp-provenance gate for the deep-closure claim (A41).

`clamp_fraction_radiative` (A26) is the right reliability measure for a *radiative*
zone, but it says nothing about a layer the composite has turned convective — and
closing the deep line-list radiative window is exactly that. `clamp_fraction_all`
and `structure.diagnostics` expose where the κ doing the closing came from.

These gates pin two things: the accounting is self-consistent, and the current
fiducial's deep closure is recorded as resting on *interpolated* FP data. It did
not always: A41 recorded the opposite, every evaluation clamped below the then
1000 K grid floor. The extension campaign (A50) added T = 650/800 K at
ρ = 100-500 and reversed it (A51). If this starts failing again — a different
fiducial, a changed blend — that is a result, not a regression: re-report it.
"""

from __future__ import annotations

import numpy as np
import pytest

from eos import AquaEOS
from figures._common import comparison_band
from opacity import build_source
from structure import build_structure
from structure.diagnostics import clamp_audit, clamp_reasons

FP_T_FLOOR_K = 650.0  # FP grid's coldest point since the extension campaign (A50)
DEEP_RHO_RANGE = (400.0, 550.0)  # where the line list leaves a radiative sliver (A36)


@pytest.fixture(scope="module")
def eos(cfg):
    return AquaEOS.from_config(cfg)


@pytest.fixture(scope="module")
def audit(cfg, eos):
    band = comparison_band(cfg, ["composite", "hitemp"])
    return clamp_audit(cfg, 500.0, 20.0, "composite", DEEP_RHO_RANGE, eos=eos,
                       n_sample=800, band=band)


def test_clamp_fraction_all_matches_the_profile_mask(audit):
    """The new field is the clamped fraction over every sampled point."""
    res = audit["result"]
    prof = res.profile
    whole = audit["whole_profile"]
    assert whole["n_points"] == prof["rho_si"].size
    assert res.clamp_fraction_all == pytest.approx(whole["n_clamped"] / whole["n_points"])
    assert 0.0 <= res.clamp_fraction_all <= 1.0
    # It is a different measure from the existing two, and both survive (A41).
    assert hasattr(res, "clamp_fraction") and hasattr(res, "clamp_fraction_radiative")


def test_fiducial_deep_closure_now_rests_on_interpolated_fp_data(audit):
    """A51 (was A41, REVERSED): the closure at ρ∈[400,550] is FP data now.

    A41 recorded that every one of the 46 FP evaluations closing this range was a
    nearest-neighbour clamp 300 K below the grid's 1000 K floor, which is why the
    headline was demoted. The extension campaign added T = 650 and 800 K at
    ρ = 100-500 — exactly this region — so the same evaluations now interpolate.
    The gate is inverted accordingly and pins the new state.
    """
    fp = audit["in_range"]["by_source"].get("fp")
    assert fp is not None and fp["n_points"] > 0, (
        "the composite no longer uses FP at ρ∈[400,550] — the closure claim has moved"
    )
    assert fp["n_clamped"] == 0, (
        f"{fp['n_clamped']} of {fp['n_points']} FP evaluations at ρ∈[400,550] are "
        "clamped again: the deep closure has stopped being interpolated FP data. "
        "That is a result — re-report it, do not relax this gate."
    )
    assert fp["n_below_t_floor"] == 0
    # The envelope here is ~766-790 K, above the extension grid's 650 K floor.
    assert min(audit["in_range_temperature_k"]) >= FP_T_FLOOR_K


def test_clamp_reasons_attributes_to_the_used_subsource(cfg):
    """Below the blend window the line list is used; above it, FP."""
    comp = build_source("composite", cfg)
    n = 8
    dilute = clamp_reasons(comp, np.full(n, 1200.0), np.full(n, 10.0))
    assert set(dilute["by_source"]) == {"hitemp"}
    assert dilute["n_clamped"] == 0  # in-grid line-list query

    # (700 K, 500) was clamped below the old 1000 K floor; the extension covers
    # it, so it interpolates now (A51). Below 650 K there is still no FP data.
    dense_warm = clamp_reasons(comp, np.full(n, 700.0), np.full(n, 500.0))
    assert set(dense_warm["by_source"]) == {"fp"}
    assert dense_warm["by_source"]["fp"]["n_clamped"] == 0

    dense_cold = clamp_reasons(comp, np.full(n, 500.0), np.full(n, 500.0))
    assert set(dense_cold["by_source"]) == {"fp"}
    entry = dense_cold["by_source"]["fp"]
    assert entry["n_clamped"] == n
    assert entry["n_below_t_floor"] == n
    assert entry["grid_t_k"][0] == pytest.approx(FP_T_FLOOR_K)

    dense_hot = clamp_reasons(comp, np.full(n, 2000.0), np.full(n, 500.0))
    assert dense_hot["by_source"]["fp"]["n_clamped"] == 0  # inside the FP grid


def test_clamp_fraction_all_is_serialized(cfg, eos):
    """The field travels in the json summary with the other reliability metadata."""
    res = build_structure(cfg, 500.0, 100.0, "hitemp", eos=eos, n_sample=200)
    summary = res.summary()
    assert "clamp_fraction_all" in summary
    assert summary["clamp_fraction_all"] == res.clamp_fraction_all
