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


def test_registered_composite_needs_line_list(cfg):
    """Building the registered composite fails clearly until HITEMP lands."""
    with pytest.raises(FileNotFoundError, match="not yet delivered"):
        build_source("composite", cfg)
