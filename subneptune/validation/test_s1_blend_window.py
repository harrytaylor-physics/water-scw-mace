"""Blend-window gate: the deep radiative zones stay closed (A35).

The headline structural result (A28) is that pure HITEMP opens radiative layers
at ρ ≳ 300 kg/m³ in the cold fiducial (T_int = 20 K) and the FP-informed
composite closes them. The composite's surviving RCB, however, sits *inside* the
blend window, so the claim must not depend on the blend prescription. This gate
runs the S1 scan — windows [100,200], [100,300] (default), [100,500] and a sharp
switch at the FP grid floor — and asserts that no composite case leaves a
radiative layer denser than 300 kg/m³.

The *number* of deep zones and their exact density extent are sampling-sensitive
(the zones are thin slivers straddling the steam→dense-water density jump at
~180-230 bar), so the gate asserts presence for HITEMP and absence for the
composite, never a count. Plus cheap checks that the new config knobs behave:
``blend_mode='linear'`` reproduces the pre-existing smoothstep blend exactly.
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest

from config.schema import _validate
from eos import AquaEOS
from figures.plot_s1_blend_sensitivity import CASES, DEEP_RHO_SI, run_cases
from opacity import ConstantOpacitySource, CompositeOpacitySource
from structure import deep_radiative_zones


@pytest.fixture(scope="module")
def eos(cfg):
    return AquaEOS.from_config(cfg)


@pytest.fixture(scope="module")
def scan(cfg, eos):
    """The S1 scan itself (HITEMP + every composite blend case)."""
    rows, results = run_cases(cfg, eos=eos)
    return results


def test_hitemp_opens_a_deep_radiative_zone(scan):
    """Positive control: without it the closure assertion below is vacuous."""
    zones = deep_radiative_zones(scan["HITEMP"], DEEP_RHO_SI)
    assert zones, (
        f"HITEMP no longer opens a radiative zone at ρ > {DEEP_RHO_SI:g} kg/m³ in the "
        "cold fiducial; the composite-closure claim (A28) has lost its control case"
    )


@pytest.mark.parametrize("label", [c[0] for c in CASES])
def test_deep_radiative_zones_closed_for_every_blend_case(scan, label):
    """A35: no composite blend prescription leaves a radiative layer at ρ > 300."""
    res = scan[label]
    zones = deep_radiative_zones(res, DEEP_RHO_SI)
    assert zones == [], (
        f"blend case {label!r} leaves {len(zones)} radiative zone(s) at "
        f"ρ > {DEEP_RHO_SI:g} kg/m³: {zones}; max radiative ρ = "
        f"{res.max_radiative_rho_si:.1f} kg/m³"
    )


def test_linear_mode_is_the_unchanged_smoothstep_blend():
    """'linear' labels the existing C¹ smoothstep log-blend (A12) — unchanged."""
    dilute = ConstantOpacitySource(0.5, name="hitemp")
    dense = ConstantOpacitySource(5.0, name="fp")
    c = CompositeOpacitySource(dilute, dense, blend_rho_si=(100.0, 300.0))
    assert c.mode == "linear"
    rho = np.array([120.0, 173.205, 250.0])
    x = (np.log10(rho) - np.log10(100.0)) / (np.log10(300.0) - np.log10(100.0))
    w = x * x * (3.0 - 2.0 * x)
    expected = 10.0 ** ((1 - w) * np.log10(0.5) + w * np.log10(5.0))
    assert c.kappa_R(np.full(3, 1500.0), rho) == pytest.approx(expected)


def test_sharp_mode_switches_at_the_window_floor():
    """'sharp' is a step at ρ_lo: pure line-list below, pure FP at and above."""
    dilute = ConstantOpacitySource(0.5, name="hitemp")
    dense = ConstantOpacitySource(5.0, name="fp")
    c = CompositeOpacitySource(dilute, dense, blend_rho_si=(108.796, 108.796), mode="sharp")
    assert c.seam_rho_si == pytest.approx(108.796)
    assert c.kappa_R(1500.0, 108.0) == pytest.approx(0.5)
    assert c.kappa_R(1500.0, 108.796) == pytest.approx(5.0)
    assert c.kappa_R(1500.0, 1000.0) == pytest.approx(5.0)


def test_zero_width_window_requires_sharp_mode(cfg):
    """A zero-width window is meaningless for a blend; config + source both refuse."""
    with pytest.raises(ValueError):
        CompositeOpacitySource(ConstantOpacitySource(1.0), ConstantOpacitySource(1.0),
                               blend_rho_si=(109.0, 109.0), mode="linear")
    bad = replace(cfg, opacity=replace(cfg.opacity, blend_window_kg_m3=(109.0, 109.0)))
    with pytest.raises(ValueError):
        _validate(bad)
    ok = replace(cfg, opacity=replace(cfg.opacity, blend_window_kg_m3=(109.0, 109.0),
                                      blend_mode="sharp"))
    _validate(ok)  # allowed: a sharp switch has no width


def test_blend_window_alias_tracks_the_canonical_field(cfg):
    """``composite_blend_rho`` is a deprecated read-only alias (A35)."""
    c = replace(cfg, opacity=replace(cfg.opacity, blend_window_kg_m3=(30.0, 100.0)))
    assert c.opacity.composite_blend_rho == (30.0, 100.0)
