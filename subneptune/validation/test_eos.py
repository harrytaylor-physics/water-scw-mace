"""EOS gate G2 + AQUA loader/inversion checks (SPEC §4, §7).

G2: the AQUA loader reproduces >= 10 published table values to <= 0.5%. We
round-trip single-phase grid nodes (ρ, s, grad_ad) through the (P,T) inversion,
and additionally anchor to independent physical values (liquid water, ideal-gas
vapor). Also covers the coexistence-plateau handling (AMENDMENTS A8) and the
Zeng (2016) core radius.
"""

from __future__ import annotations

import numpy as np
import pytest

from eos import AquaEOS, core_radius_rearth, core_radius_si
from eos.core import EARTH_LIKE_CMF
from eos.loader import load_aqua_tables


@pytest.fixture(scope="module")
def eos(cfg):
    return AquaEOS.from_config(cfg)


@pytest.fixture(scope="module")
def tables(cfg):
    return load_aqua_tables(str(cfg.aqua_table_path))


def _single_phase_nodes(tables, n, seed=2):
    """Pick n single-phase (phase id > 2) grid nodes spread across the table."""
    rng = np.random.default_rng(seed)
    picks = []
    while len(picks) < n:
        i = int(rng.integers(0, tables.n_rho))
        j = int(rng.integers(0, tables.n_T))
        if tables.phase[i, j] > 2:  # vapor/liquid/supercritical
            picks.append((i, j))
    return picks


# --------------------------------------------------------------------------
# G2 — reproduce published table values to <= 0.5%
# --------------------------------------------------------------------------
def test_g2_single_phase_roundtrip(eos, tables):
    nodes = _single_phase_nodes(tables, 40)
    assert len(nodes) >= 10
    for i, j in nodes:
        p_ij = tables.press[i, j]
        t_j = 10.0 ** tables.logT_grid[j]
        rho_i = 10.0 ** tables.logrho_grid[i]
        st = eos.state(p_ij, t_j)
        assert abs(st.rho - rho_i) / rho_i <= 5e-3
        assert abs(st.grad_ad - tables.ad_grad[i, j]) <= 5e-3 * abs(tables.ad_grad[i, j]) + 1e-6
        s_ij = tables.entropy[i, j]
        assert abs(st.entropy - s_ij) <= 5e-3 * abs(s_ij) + 1.0  # J/kg/K floor


def test_g2_internal_energy_and_sound_roundtrip(eos, tables):
    for i, j in _single_phase_nodes(tables, 15, seed=7):
        p_ij = tables.press[i, j]
        t_j = 10.0 ** tables.logT_grid[j]
        st = eos.state(p_ij, t_j)
        u_ij, c_ij = tables.u[i, j], tables.sound[i, j]
        assert abs(st.internal_energy - u_ij) <= 5e-3 * abs(u_ij) + 1.0
        assert abs(st.sound_speed - c_ij) <= 5e-3 * abs(c_ij) + 1e-3


# --------------------------------------------------------------------------
# Independent physical anchors
# --------------------------------------------------------------------------
def test_liquid_water_density(eos):
    # 1 bar, 300 K -> liquid water ~997 kg/m^3 (IAPWS); off-node interpolation.
    st = eos.state(1e5, 300.0)
    assert st.phase == 4  # liquid
    assert st.rho == pytest.approx(997.0, rel=0.03)


def test_superheated_vapor_is_near_ideal_gas(eos):
    # 1 bar, 500 K: P << P_sat -> vapor; close to ideal gas.
    st = eos.state(1e5, 500.0)
    rho_ideal = 1e5 * 0.018015 / (8.314 * 500.0)
    assert st.phase == 3  # vapor
    assert st.rho == pytest.approx(rho_ideal, rel=0.02)


def test_grad_ad_in_physical_range(eos):
    rng = np.random.default_rng(3)
    p = 10.0 ** rng.uniform(4, 11, size=200)
    t = 10.0 ** rng.uniform(2.2, 3.7, size=200)
    grad = eos.grad_ad(p, t)
    assert np.all(grad >= -1e-9)
    assert np.all(grad <= 1.2)  # (dlnT/dlnP)_S can exceed 1 near transitions


# --------------------------------------------------------------------------
# Coexistence-plateau handling (AMENDMENTS A8)
# --------------------------------------------------------------------------
def test_dome_crossing_selects_branches(eos):
    # At 500 K: below P_sat(~26 bar) vapor (low ρ); above it liquid (high ρ).
    rho_below = eos.rho(2.0e6, 500.0)  # 20 bar
    rho_above = eos.rho(3.0e6, 500.0)  # 30 bar
    assert rho_below < 50.0
    assert rho_above > 500.0


def test_logP_strict_is_strictly_increasing(tables):
    d = np.diff(tables.logP_strict, axis=0)
    assert np.all(d > 0)  # every column strictly increasing -> invertible


# --------------------------------------------------------------------------
# Core mass-radius (Zeng 2016)
# --------------------------------------------------------------------------
def test_core_radius_zeng():
    # Earth: 1 M_earth, Earth-like CMF -> ~1 R_earth.
    assert core_radius_rearth(1.0) == pytest.approx(1.07 - 0.21 * EARTH_LIKE_CMF, rel=1e-12)
    # Fiducial 5 M_earth core.
    assert core_radius_rearth(5.0) == pytest.approx(1.547, rel=1e-2)
    assert core_radius_si(5.0) == pytest.approx(core_radius_rearth(5.0) * 6.3781e6, rel=1e-9)


def test_core_radius_validation():
    with pytest.raises(ValueError):
        core_radius_rearth(-1.0)
    with pytest.raises(ValueError):
        core_radius_rearth(5.0, cmf=1.5)


# --------------------------------------------------------------------------
# Vectorized API shape
# --------------------------------------------------------------------------
def test_vector_api_shapes(eos):
    p = np.array([1e5, 1e6, 1e8, 1e10])
    t = np.array([500.0, 800.0, 1500.0, 2000.0])
    assert eos.rho(p, t).shape == (4,)
    assert eos.grad_ad(p, t).shape == (4,)
    assert np.all(np.isfinite(eos.rho(p, t)))
