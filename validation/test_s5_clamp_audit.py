"""Clamp-provenance gate for the deep-closure claim (A41).

`clamp_fraction_radiative` (A26) is the right reliability measure for a *radiative*
zone, but it says nothing about a layer the composite has turned convective — and
closing the deep line-list radiative window is exactly that. `clamp_fraction_all`
and `structure.diagnostics` expose where the κ doing the closing came from.

These gates pin two things: the accounting is self-consistent, and the current
fiducial's deep closure is recorded as resting on *clamped* (nearest-neighbour)
FP opacity. If the second ever starts failing — an extended FP grid, a different
fiducial — that is a result, not a regression: read the message and re-report.
"""

from __future__ import annotations

import numpy as np
import pytest

from eos import AquaEOS
from figures._common import comparison_band
from opacity import build_source
from structure import build_structure
from structure.diagnostics import clamp_audit, clamp_reasons

FP_T_FLOOR_K = 1000.0  # FP grid's coldest point
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


def test_fiducial_deep_closure_rests_on_clamped_fp_opacity(audit):
    """A41: at T_int=20 K, T_eq=500 K the FP κ closing ρ∈[400,550] is extrapolated."""
    fp = audit["in_range"]["by_source"].get("fp")
    assert fp is not None and fp["n_points"] > 0, (
        "the composite no longer uses FP at ρ∈[400,550] — the closure claim has moved"
    )
    assert fp["n_clamped"] == fp["n_points"], (
        f"only {fp['n_clamped']} of {fp['n_points']} FP evaluations at ρ∈[400,550] are "
        "clamped: some of the deep closure is now interpolated FP data. That is a "
        "result — re-report it, do not relax this gate."
    )
    assert fp["n_below_t_floor"] == fp["n_points"], (
        "the clamps at ρ∈[400,550] are no longer all below the FP grid's T floor"
    )
    assert min(audit["in_range_temperature_k"]) < FP_T_FLOOR_K


def test_clamp_reasons_attributes_to_the_used_subsource(cfg):
    """Below the blend window the line list is used; above it, FP."""
    comp = build_source("composite", cfg)
    n = 8
    dilute = clamp_reasons(comp, np.full(n, 1200.0), np.full(n, 10.0))
    assert set(dilute["by_source"]) == {"hitemp"}
    assert dilute["n_clamped"] == 0  # in-grid line-list query

    dense_cold = clamp_reasons(comp, np.full(n, 700.0), np.full(n, 500.0))
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
