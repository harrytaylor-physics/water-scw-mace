"""Evolution gate G6 + cooling-track sanity (SPEC §6, §7).

G6  Energy bookkeeping: |ΔE − ∫L dt| / E < 1% over the full track.

Also checks the physical monotonicity of a cooling track (radius contracts and
luminosity falls with age) and that a more opaque source reaches a given radius
later ("stay puffy longer"), the definitive time-integrated version of F3.
"""

from __future__ import annotations

import numpy as np
import pytest

import opacity.build_tables as bt
from eos import AquaEOS
from evolution import cooling_time_ratio, cooling_track
from opacity import ConstantOpacitySource


def _register_constant(cfg, name, kappa_si):
    bt._cache[bt.cache_key(name, cfg, None)] = ConstantOpacitySource(kappa_si, name=name)


@pytest.fixture(scope="module")
def eos(cfg):
    return AquaEOS.from_config(cfg)


# Grids are kept small: build_structure's shooting (~2 s) dominates, so the
# evolution gate uses the coarsest grid that still meets each assertion. G6 needs
# ~n=40 for the energy quadrature; the monotonicity/ratio checks do not.
_G6_GRID = np.geomspace(25.0, 500.0, 40)
_COARSE_GRID = np.geomspace(25.0, 500.0, 12)


def test_g6_energy_bookkeeping(cfg, eos):
    _register_constant(cfg, "const_evo", 0.1)
    tr = cooling_track(cfg, "const_evo", 600.0, eos=eos, t_int_grid_k=_G6_GRID, n_sample=120)
    assert tr.energy_residual < 1e-2  # G6


def test_cooling_track_monotone(cfg, eos):
    _register_constant(cfg, "const_evo", 0.1)
    tr = cooling_track(cfg, "const_evo", 600.0, eos=eos, t_int_grid_k=_COARSE_GRID, n_sample=120)
    order = np.argsort(tr.age_gyr)
    radius = tr.radius_p_rearth[order]
    lum = tr.luminosity_w[order]
    # Age increases as T_int decreases: radius contracts and L falls.
    assert np.all(np.diff(radius) <= 1e-9)
    assert np.all(np.diff(lum) <= 0)
    assert tr.age_gyr.min() == pytest.approx(0.0, abs=1e-12)
    assert tr.radius_at_age(tr.age_gyr.max()) < tr.radius_p_rearth.max()  # contracts with age


def test_history_comparison_diagnostics(cfg, eos):
    """A29: composite-vs-line-list history has a well-defined FP-relevant epoch
    (RCB reaches ρ≥109 only at late/cold ages) and the escape window is clean."""
    from evolution import compare_histories

    g = np.geomspace(15.0, 400.0, 12)  # coarse but reaches the cold epoch
    h = compare_histories(cfg, "composite", "hitemp", 600.0, eos=eos, t_int_grid_k=g)
    # Young ages (shallow RCB) -> the two opacities agree; FP-relevance is late.
    assert abs(h.diff_at_age(0.1)) < 5e-3  # <0.5% at 100 Myr (escape window clean)
    ep = h.fp_relevant_epoch_gyr()
    assert ep is not None and ep[0] > 0.5  # FP-relevant epoch starts well after 0.5 Gyr


def test_more_opaque_stays_puffy_longer(cfg, eos):
    """Definitive Phase-2 statement: at a fixed radius, the more opaque source is
    older (took longer to contract there)."""
    _register_constant(cfg, "evo_hi", 1.0)  # more opaque
    _register_constant(cfg, "evo_lo", 0.1)  # less opaque
    hi = cooling_track(cfg, "evo_hi", 600.0, eos=eos, t_int_grid_k=_COARSE_GRID, n_sample=120)
    lo = cooling_track(cfg, "evo_lo", 600.0, eos=eos, t_int_grid_k=_COARSE_GRID, n_sample=120)
    # Compare at a radius both tracks reach.
    r_common = 0.5 * (max(hi.radius_p_rearth.min(), lo.radius_p_rearth.min())
                      + min(hi.radius_p_rearth.max(), lo.radius_p_rearth.max()))
    assert cooling_time_ratio(hi, lo, r_common) > 1.0
