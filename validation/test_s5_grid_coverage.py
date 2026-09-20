"""FP grid-coverage accounting and the campaign specification (A48).

The S5 coverage figure claims, per envelope trajectory, what fraction of the
envelope mass the composite evaluated from real FP data and what fraction it
extrapolated — and then emits a list of (T, ρ) nodes as the specification for
the next opacity campaign. Both are only as good as their bookkeeping:

  * every evaluation gets exactly one category, and the categories partition the
    envelope mass;
  * a clamp inside the grid's (T, ρ) box is distinguished from one outside its
    range — the first can be fixed by filling the scattered grid's gaps (A6), the
    second needs the grid extended, and they are different asks;
  * the emitted node list actually **encloses** the samples it claims to, tested
    with the same triangulated interpolator the model uses, not asserted.
"""

from __future__ import annotations

import numpy as np
import pytest

from eos import AquaEOS
from opacity import get_source
from structure import build_structure
from structure.diagnostics import (
    CATEGORIES,
    campaign_nodes,
    category_mass_fractions,
    classify_evaluations,
)

FP_T_FLOOR_K = 1000.0
FP_RHO_FLOOR_SI = 108.796


@pytest.fixture(scope="module")
def eos(cfg):
    return AquaEOS.from_config(cfg)


@pytest.fixture(scope="module")
def composite(cfg):
    return get_source("composite", cfg)


@pytest.fixture(scope="module")
def old_fiducial(cfg, eos):
    """A41's case: the closure that turned out to be clamped FP opacity."""
    return build_structure(cfg, 500.0, 20.0, "composite", eos=eos, n_sample=600)


def test_every_evaluation_gets_exactly_one_known_category(composite, old_fiducial):
    prof = old_fiducial.profile
    cats = classify_evaluations(composite, prof["temperature_k"], prof["rho_si"])
    assert cats.size == prof["temperature_k"].size
    assert set(np.unique(cats)) <= set(CATEGORIES)


def test_categories_partition_the_envelope_mass(composite, old_fiducial):
    prof = old_fiducial.profile
    cats = classify_evaluations(composite, prof["temperature_k"], prof["rho_si"])
    frac = category_mass_fractions(old_fiducial, cats)
    assert set(frac) == set(CATEGORIES)
    assert sum(frac.values()) == pytest.approx(1.0)
    assert all(0.0 <= v <= 1.0 for v in frac.values())


def test_clamped_categories_agree_with_the_clamp_mask(composite, old_fiducial):
    """The three clamp categories must sum to the source's own clamp count."""
    prof = old_fiducial.profile
    t, rho = prof["temperature_k"], prof["rho_si"]
    cats = classify_evaluations(composite, t, rho)
    clamped = np.asarray(composite.clamped_mask(t, rho)).ravel()
    # dilute-flat points are below a line list's ρ floor: justified, not clamps
    # (A26), and the hull test does not count them either.
    labelled = np.isin(cats.astype(str), ["clamped_t", "clamped_rho", "in_box_out_of_hull"])
    assert np.array_equal(labelled, clamped & (cats.astype(str) != "dilute_flat"))


def test_the_old_fiducial_runs_on_extrapolated_opacity(composite, old_fiducial):
    """A41/A48, restated as coverage — and what the extension changed (A51).

    The cold-clamped tail is essentially gone: the extension campaign put FP data
    at T = 650/800 K over ρ = 100-500, so ``clamped_t`` falls from >1% of the
    envelope mass to ~0.02%. What remains is the *other* ask A48 separated out —
    the scattered grid's own interior gaps — which still carry ~99% of the mass.
    So the old fiducial is still not running on interpolated FP data; the reason
    has changed from "too cold for the grid" to "inside the box, outside the
    hull", which is fillable without extending the grid's range at all.
    """
    prof = old_fiducial.profile
    cats = classify_evaluations(composite, prof["temperature_k"], prof["rho_si"])
    frac = category_mass_fractions(old_fiducial, cats)
    assert frac["interpolated"] < 0.01, (
        f"the old fiducial now interpolates {frac['interpolated']:.3f} of its envelope "
        f"mass — A41/A48's finding has moved")
    assert frac["in_box_out_of_hull"] > 0.9
    # Was > 0.01 before the extension; now a residual tail below the 650 K floor.
    assert frac["clamped_t"] < 0.01
    assert frac["clamped_t"] > 0.0


def test_campaign_nodes_enclose_what_they_claim(cfg, composite, old_fiducial):
    """The emitted lattice nodes must actually put the needed samples in-hull."""
    table = get_source("fp", cfg).table
    prof = old_fiducial.profile
    cats = classify_evaluations(composite, prof["temperature_k"], prof["rho_si"])
    spec = campaign_nodes([(old_fiducial, cats)], table, t_step_k=100.0,
                          log_rho_step=0.25, t_floor_k=FP_T_FLOOR_K,
                          rho_floor_si=FP_RHO_FLOOR_SI)
    assert spec["lattice"] == (100.0, 0.25)
    assert spec["fill_encloses"], "the fill nodes do not enclose the in-box samples"
    # Nodes sit on the lattice, and inside the box they claim to fill.
    for t_k, rho in spec["fill_nodes"]:
        assert t_k % 100.0 == pytest.approx(0.0)
        assert (np.log10(rho) / 0.25) == pytest.approx(round(np.log10(rho) / 0.25))
    if spec["n_in_box_needed"] == 0:
        assert spec["fill_nodes"] == []


def test_extension_nodes_are_reported_separately(cfg, composite, old_fiducial):
    """Cold-dense samples need the grid extended, not filled: different list."""
    table = get_source("fp", cfg).table
    prof = old_fiducial.profile
    cats = classify_evaluations(composite, prof["temperature_k"], prof["rho_si"])
    spec = campaign_nodes([(old_fiducial, cats)], table, t_floor_k=FP_T_FLOOR_K,
                          rho_floor_si=FP_RHO_FLOOR_SI)
    t_lo = float(np.min(table.temperature_k))
    t_step = spec["lattice"][0]
    # This case is the reason the addendum exists: its dense envelope is ~630-760 K.
    assert spec["n_extend_needed"] > 0
    assert spec["extend_nodes"]
    # Bracketing corners: a sample just below the floor brackets to the lattice
    # cell spanning it, so a node may sit up to one step ABOVE the floor (A51 —
    # with the floor now 650 K, the 600-700 K cell straddles it). None may sit
    # further above, or it would not be an extension.
    assert all(t <= t_lo + t_step for t, _ in spec["extend_nodes"])
    assert any(t < t_lo for t, _ in spec["extend_nodes"])
    assert all(rho > 0.0 for _, rho in spec["extend_nodes"])
