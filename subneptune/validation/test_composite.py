"""Composite opacity source gates (AMENDMENTS A11, A12).

Blending limits (pure line-list dilute / pure FP dense), log-mean at the seam,
monotone transition, seam-ratio diagnostic, usage-weighted clamp attribution,
and the 'line-list not yet delivered' behaviour of the registered composite.
"""

from __future__ import annotations

import numpy as np
import pytest

from opacity import (
    CompositeOpacitySource,
    ConstantOpacitySource,
    GridOpacitySource,
    OpacityTable,
    build_source,
    compute_point_means,
    resolve_band,
)
from opacity.table import CM2_PER_G_TO_SI


def _comp(lo=30.0, hi=100.0):
    dilute = ConstantOpacitySource(0.5, name="hitemp")  # 5 cm^2/g
    dense = ConstantOpacitySource(5.0, name="fp")  # 50 cm^2/g
    return CompositeOpacitySource(dilute, dense, blend_rho_si=(lo, hi))


def test_blend_limits():
    c = _comp()
    assert c.kappa_R(1500.0, 1.0) == pytest.approx(0.5)  # dilute below window
    assert c.kappa_R(1500.0, 30.0) == pytest.approx(0.5)  # window bottom
    assert c.kappa_R(1500.0, 100.0) == pytest.approx(5.0)  # window top
    assert c.kappa_R(1500.0, 1e4) == pytest.approx(5.0)  # dense above window


def test_seam_is_log_mean():
    c = _comp()
    seam = c.kappa_R(1500.0, c.seam_rho_si)
    assert seam == pytest.approx(10 ** ((np.log10(0.5) + np.log10(5.0)) / 2))
    assert c.seam_rho_si == pytest.approx(np.sqrt(30.0 * 100.0))


def test_blend_is_monotone_in_window():
    c = _comp()
    rho = np.linspace(30.0, 100.0, 20)
    k = c.kappa_R(1500.0, rho)
    assert np.all(np.diff(k) >= -1e-12)  # dense > dilute here, so increasing


def test_seam_ratio_diagnostic():
    c = _comp()
    d = c.seam_ratio([1000.0, 2500.0])
    assert d["kappa_R_fp_over_ll"] == pytest.approx([10.0, 10.0])
    assert d["kappa_P_fp_over_ll"] == pytest.approx([10.0, 10.0])
    assert d["seam_rho_si"] == pytest.approx(np.sqrt(30.0 * 100.0))


def test_clamp_attributed_to_used_source(cfg):
    """In the dilute regime the line-list (in-hull) is used, so the composite is
    reliable even though FP is extrapolated there (down-weighted)."""
    band = resolve_band(cfg)
    t, rho, kr, kp = compute_point_means("fp", cfg, band)
    fp = GridOpacitySource("fp", OpacityTable(t, rho, kr * CM2_PER_G_TO_SI, kp * CM2_PER_G_TO_SI, band))
    dilute = ConstantOpacitySource(0.5, name="hitemp")  # covers everything
    c = CompositeOpacitySource(dilute, fp, blend_rho_si=(30.0, 100.0))

    c.reset_clamp_counter()
    c.kappa_R(np.full(50, 700.0), np.full(50, 1.0))  # dilute -> line-list used
    assert c.clamp_fraction == 0.0  # FP clamps here but is not the used source


def test_blend_window_validation():
    d = ConstantOpacitySource(1.0)
    with pytest.raises(ValueError):
        CompositeOpacitySource(d, d, blend_rho_si=(100.0, 30.0))  # lo > hi
    with pytest.raises(ValueError):
        CompositeOpacitySource(d, d, blend_rho_si=(-1.0, 30.0))  # lo <= 0


def test_dilute_flat_extrapolation(cfg):
    """Below the grid ρ floor, a line-list source holds κ flat in ρ (A26) and
    counts the call as flat_rho, not as a hull clamp."""
    ht = build_source("hitemp", cfg)
    floor = ht.rho_floor_si
    ht.reset_clamp_counter()
    # A grid T (in range) at ρ below the floor -> flat κ in ρ.
    k_floor = ht.kappa_R(1000.0, floor)
    k_below = ht.kappa_R(1000.0, floor * 1e-3)
    assert k_below == pytest.approx(k_floor, rel=1e-9)  # flat in ρ
    assert ht.flat_rho_fraction > 0.0
    assert ht.clamp_fraction == 0.0  # not counted as an unreliable clamp


def test_fp_not_dilute_flat(cfg):
    """FP (condensed-phase) does NOT flat-extrapolate; below its floor it clamps."""
    fp = build_source("fp", cfg)
    assert fp.flat_rho_fraction == 0.0
    fp.reset_clamp_counter()
    fp.kappa_R(2000.0, 1e-4)  # far below FP's ρ floor (~109)
    assert fp.clamp_fraction == 1.0
    assert fp.flat_rho_fraction == 0.0


def test_registered_composite_builds(cfg):
    """The registered composite (HITEMP dilute + FP dense) builds now that HITEMP
    is delivered, and blends across the config density window."""
    comp = build_source("composite", cfg)
    assert isinstance(comp, CompositeOpacitySource)
    assert comp.dilute.name == "hitemp"
    assert comp.dense.name == "fp"
    lo, hi = comp.blend_rho_si
    # Below the window -> dilute (HITEMP); above -> dense (FP).
    assert comp.kappa_R(700.0, lo * 0.5) == pytest.approx(comp.dilute.kappa_R(700.0, lo * 0.5), rel=1e-9)
