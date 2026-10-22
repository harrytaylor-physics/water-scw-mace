"""Opacity-perturbation gate: the FP table's uncertainty does not flip the result (A37).

``opacity.fp_scale`` multiplies the FP table's κ before blending. The two claims
the manuscript rests on must survive a ±20-25% perturbation of it:

  1. RCB closure — no composite case leaves a radiative layer at ρ > 300 kg/m³
     in the cold fiducial (the A28/A35 structural claim);
  2. the sign of the composite−HITEMP interior-entropy bias at 7 Gyr (the A29
     cooling-history claim) is the same for every scale factor.

The gate shares the S2 figure's T_int grid (the age integral needs it) but
samples each profile more coarsely: it checks closure and signs, not publication
numbers.
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest

from config.schema import _validate
from eos import AquaEOS
from figures.plot_s2_opacity_scaling import DEEP_RHO_SI, FP_SCALES, run_scan
from opacity import build_source, clear_cache
from structure import deep_radiative_zones

# E2's T_int grid (A29): the age integral needs ≳20 nodes — on a 14-node grid the
# cumulative-Simpson age overshoots and returns negative ages (A38) — and reusing
# E2's grid makes the gate's numbers directly comparable to the S2 figure. Cost is
# kept down with a coarser profile sampling instead.
_TRACK_GRID_K = np.geomspace(12.0, 450.0, 44)
_AGE_REF_GYR = 7.0


@pytest.fixture(scope="module")
def eos(cfg):
    return AquaEOS.from_config(cfg)


@pytest.fixture(scope="module")
def scan(cfg, eos):
    return run_scan(cfg, eos=eos, n_sample=800, t_int_grid_k=_TRACK_GRID_K,
                    n_sample_track=250, age_ref_gyr=_AGE_REF_GYR)


@pytest.mark.parametrize("scale", FP_SCALES)
def test_rcb_closure_holds_for_every_fp_scale(scan, scale):
    """A37: the deep radiative window stays closed under the κ perturbation."""
    res = scan["structures"][scale]
    zones = deep_radiative_zones(res, DEEP_RHO_SI)
    assert zones == [], (
        f"fp_scale={scale:g} leaves {len(zones)} radiative zone(s) at "
        f"ρ > {DEEP_RHO_SI:g} kg/m³: {zones}; max radiative ρ = "
        f"{res.max_radiative_rho_si:.1f} kg/m³"
    )


def test_entropy_bias_sign_is_scale_independent(scan):
    """A37: composite−HITEMP entropy at 7 Gyr keeps its sign across the scan."""
    deltas = [float(r["delta_entropy_j_kg_k"]) for r in scan["rows"]]
    assert all(d != 0.0 for d in deltas), f"a zero entropy bias is not a sign: {deltas}"
    signs = {np.sign(d) for d in deltas}
    assert len(signs) == 1, (
        f"the composite−HITEMP entropy bias at {_AGE_REF_GYR:.0f} Gyr changes sign "
        f"across fp_scale {FP_SCALES}: {deltas} J/kg/K"
    )
    # Every track must actually reach the reference age (no silent extrapolation).
    for row in scan["rows"]:
        assert float(row["track_age_max_gyr"]) >= _AGE_REF_GYR


def test_fp_scale_multiplies_the_fp_table_only(cfg):
    """The knob scales FP's κ_R and κ_P and leaves the line list untouched."""
    clear_cache()
    base = build_source("fp", cfg)
    scaled = build_source("fp", replace(cfg, opacity=replace(cfg.opacity, fp_scale=1.25)))
    t, rho = 2000.0, 500.0
    assert scaled.kappa_R(t, rho) == pytest.approx(1.25 * base.kappa_R(t, rho), rel=1e-12)
    assert scaled.kappa_P(t, rho) == pytest.approx(1.25 * base.kappa_P(t, rho), rel=1e-12)

    ht = build_source("hitemp", cfg)
    ht_scaled = build_source("hitemp", replace(cfg, opacity=replace(cfg.opacity, fp_scale=1.25)))
    assert ht_scaled.kappa_R(1000.0, 10.0) == pytest.approx(ht.kappa_R(1000.0, 10.0), rel=1e-12)


def test_composite_inherits_fp_scale_in_the_dense_limit(cfg):
    """Above the blend window the composite is pure FP, so it carries the scale."""
    comp = build_source("composite", cfg)
    comp_scaled = build_source("composite",
                               replace(cfg, opacity=replace(cfg.opacity, fp_scale=0.8)))
    dense_rho = 800.0  # well above the window -> pure FP
    assert comp_scaled.kappa_R(2000.0, dense_rho) == pytest.approx(
        0.8 * comp.kappa_R(2000.0, dense_rho), rel=1e-12)
    dilute_rho = 10.0  # well below the window -> pure line list, unaffected
    assert comp_scaled.kappa_R(1000.0, dilute_rho) == pytest.approx(
        comp.kappa_R(1000.0, dilute_rho), rel=1e-12)


def test_fp_scale_must_be_positive(cfg):
    with pytest.raises(ValueError):
        _validate(replace(cfg, opacity=replace(cfg.opacity, fp_scale=0.0)))
    with pytest.raises(ValueError):
        _validate(replace(cfg, opacity=replace(cfg.opacity, fp_scale=-1.0)))
