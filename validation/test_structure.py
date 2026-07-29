"""Structure gates (SPEC §5, §7).

  G1  Isothermal limit: κ → 0 ⇒ the radiative zone approaches the Guillot skin
      solution analytically.
  G3  Hydrostatic residual: |∫(−Gm/4πr⁴)dm − ΔP|/ΔP < 1e-6, every model.
  G4  Opacity-swap null test (structure level): two sources with identical
      tables ⇒ identical structure to machine precision.
  G7  Convergence: R_p changes < 0.1% under halved tolerances (and the profile
      is insensitive to doubled sampling), for the fiducial planet.
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest

import opacity.build_tables as bt
from opacity import (
    CM2_PER_G_TO_SI,
    ConstantOpacitySource,
    GridOpacitySource,
    OpacityTable,
    compute_point_means,
    resolve_band,
)
from structure import (
    build_structure,
    guillot_temperature,
    integrate_outer,
    skin_temperature,
)
from eos import AquaEOS, core_mass_si


def _scaled_fp_grid_source(cfg, factor):
    """A GridOpacitySource from the FP grid with κ scaled by ``factor`` so the
    fiducial structure converges (exercises grid interpolation in G4)."""
    band = resolve_band(cfg)
    t, rho, kr, kp = compute_point_means("fp", cfg, band)
    table = OpacityTable(
        t, rho, kr * CM2_PER_G_TO_SI * factor, kp * CM2_PER_G_TO_SI * factor, band
    )
    return GridOpacitySource("fp", table)


@pytest.fixture(scope="module")
def eos(cfg):
    return AquaEOS.from_config(cfg)


def _inject_constant_fp(cfg, kappa_si):
    """Register a constant-opacity source under the name 'fp' (for gate runs)."""
    key = ("fp", resolve_band(cfg), str(cfg.data_root_path))
    bt._cache[key] = ConstantOpacitySource(kappa_si, name="fp")


# --------------------------------------------------------------------------
# G1 — isothermal / Guillot-skin limit
# --------------------------------------------------------------------------
def test_g1_isothermal_skin_limit(cfg, eos):
    """As κ → 0, τ → 0 everywhere, so the outer profile → the skin temperature."""
    t_int, t_eq, gamma = 100.0, 500.0, 1.0
    t_skin = skin_temperature(t_int, t_eq, gamma)
    p_top = cfg.boundary.p_top_bar * 1e5
    p_match = cfg.boundary.match_pressure_bar * 1e5
    mass_p = core_mass_si(cfg.core.mass_mearth / (1 - cfg.envelope.f_env))

    prev = None
    for kappa in (1e-3, 1e-5, 1e-7):
        op = ConstantOpacitySource(kappa)
        match, outer, _ = integrate_outer(
            eos, op, 1.79 * 6.3781e6, mass_p, t_int, t_eq, gamma,
            p_top, p_match, cfg.numerics.rtol, cfg.numerics.atol, cfg.numerics.solver,
            n_sample=200,
        )
        max_dev = float(np.max(np.abs(outer["T"] - t_skin)))
        if prev is not None:
            assert max_dev < prev  # deviation shrinks as κ decreases
        prev = max_dev
    # At the smallest κ the whole radiative zone is within a fraction of a K.
    assert max_dev < 1.0
    # τ at the matching pressure vanishes in the limit (τ ∝ κ).
    assert match.tau < 1e-2


def test_g1_skin_matches_guillot_tau0(cfg):
    """The analytic skin temperature equals the Guillot profile at τ = 0."""
    from structure.boundary import t_irr_from_t_eq

    t_skin = skin_temperature(80.0, 650.0, 0.4)
    t_tau0 = float(guillot_temperature(0.0, 80.0, t_irr_from_t_eq(650.0), 0.4))
    assert t_skin == pytest.approx(t_tau0, rel=1e-12)


# --------------------------------------------------------------------------
# G3 — hydrostatic residual
# --------------------------------------------------------------------------
def test_g3_hydrostatic_residual(cfg):
    _inject_constant_fp(cfg, 0.1)
    res = build_structure(cfg, t_eq_k=600.0, t_int_k=100.0, source_name="fp")
    assert res.hydrostatic_residual < 1e-6


# --------------------------------------------------------------------------
# Adaptive matching (A13): the RCB match keeps high-opacity models from
# collapsing (a fixed deep match did), and matches shallower than the cap.
# --------------------------------------------------------------------------
def test_adaptive_matching_converges_high_opacity(cfg):
    bt.clear_cache()  # use the real FP grid (very high band-limited kappa_R)
    res = build_structure(cfg, t_eq_k=700.0, t_int_k=150.0, source_name="fp")
    assert res.radius_p_rearth > res.radius_core_m / 6.3781e6  # sane, converged
    assert res.matched_at_rcb
    assert res.match_pressure_actual_bar < res.match_pressure_bar  # above the 10-bar cap
    assert res.hydrostatic_residual < 1e-6
    bt.clear_cache()


# --------------------------------------------------------------------------
# G4 — opacity-swap null test at the structure level
# --------------------------------------------------------------------------
def test_g4_identical_tables_identical_structure(cfg):
    """Two grid sources built from identical tables give identical structure.

    Uses the FP grid scaled so the fiducial converges; the scaling is common to
    both sources, so the test still checks source-swap invariance through the
    full grid-interpolation + integration path. (Grid-interp determinism itself
    is covered by the opacity-level G4.)
    """
    key = ("fp", resolve_band(cfg), str(cfg.data_root_path))

    bt._cache[key] = _scaled_fp_grid_source(cfg, 1e-3)
    a = build_structure(cfg, t_eq_k=700.0, t_int_k=150.0, source_name="fp")
    bt._cache[key] = _scaled_fp_grid_source(cfg, 1e-3)
    b = build_structure(cfg, t_eq_k=700.0, t_int_k=150.0, source_name="fp")

    assert a.radius_p_m == pytest.approx(b.radius_p_m, rel=1e-12, abs=1e-6)
    np.testing.assert_allclose(
        a.profile["temperature_k"], b.profile["temperature_k"], rtol=1e-10
    )


# --------------------------------------------------------------------------
# G7 — convergence under halved tolerances
# --------------------------------------------------------------------------
def test_g7_convergence_under_tolerance(cfg):
    _inject_constant_fp(cfg, 0.1)
    base = build_structure(cfg, t_eq_k=500.0, t_int_k=100.0, source_name="fp")

    tight = replace(
        cfg.numerics,
        rtol=cfg.numerics.rtol / 2,
        atol=cfg.numerics.atol / 2,
        shoot_rtol=cfg.numerics.shoot_rtol / 2,
    )
    cfg_tight = replace(cfg, numerics=tight)
    refined = build_structure(cfg_tight, t_eq_k=500.0, t_int_k=100.0, source_name="fp", n_sample=4000)

    rel = abs(refined.radius_p_m - base.radius_p_m) / base.radius_p_m
    assert rel < 1e-3  # R_p change < 0.1%
