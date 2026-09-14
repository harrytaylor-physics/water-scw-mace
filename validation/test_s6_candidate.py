"""Candidate-fiducial gate (A44): the demonstration is made of FP data.

At the candidate case (config/candidate.yaml: M_c = 5 M⊕, f_env = 5%,
T_eq = 700 K, T_int = 20 K) the claim is specifically that

  * the line list leaves a radiative layer inside the FP data box
    (T ≥ 1000 K, ρ ≥ 108.796 kg/m³) — without this the demonstration is vacuous;
  * the composite closes it, leaving no in-box radiative layer;
  * and it does so with `clamp_fraction_all` < 1%, i.e. essentially every FP
    evaluation interpolated rather than clamped (contrast A41, where the old
    fiducial's closure was 100% clamped).

Also covers the A38 age-integral check that the candidate's evolution run uses.
"""

from __future__ import annotations

import numpy as np
import pytest

from config import load_config
from eos import AquaEOS
from evolution import age_grid_problem
from evolution.track import CoolingTrack
from figures._common import comparison_band
from structure import FP_GRID_FLOOR_SI, build_structure
from structure.diagnostics import radiative_layers

FP_T_FLOOR_K = 1000.0
T_INT_K = 20.0
CLAMP_LIMIT = 0.01
N_SAMPLE = 1500  # the in-box layer is a few samples thick (A42)


@pytest.fixture(scope="module")
def candidate_cfg(repo_root):
    return load_config(repo_root / "config" / "candidate.yaml")


@pytest.fixture(scope="module")
def eos(cfg):
    return AquaEOS.from_config(cfg)


@pytest.fixture(scope="module")
def models(candidate_cfg, eos):
    band = comparison_band(candidate_cfg, ["composite", "hitemp"])
    t_eq = float(candidate_cfg.irradiation.t_eq_k[0])
    return {
        name: build_structure(candidate_cfg, t_eq, T_INT_K, name, eos=eos,
                              n_sample=N_SAMPLE, band=band)
        for name in ("hitemp", "composite")
    }


def test_candidate_config_is_the_scanned_case(candidate_cfg):
    assert candidate_cfg.core.mass_mearth == 5.0
    assert candidate_cfg.envelope.f_env == 0.05
    assert candidate_cfg.irradiation.t_eq_k == (700.0,)
    assert T_INT_K in candidate_cfg.interior.t_int_k
    # The blend prescription and band are the fiducial's, so the comparison moves
    # only in (T_eq, f_env).
    assert candidate_cfg.opacity.band_cm == (100.0, 4000.0)
    assert candidate_cfg.opacity.blend_window_kg_m3 == (100.0, 300.0)
    assert candidate_cfg.opacity.fp_scale == 1.0


def test_line_list_leaves_an_in_box_radiative_layer(models):
    """Positive control: without this the closure below proves nothing."""
    layers = radiative_layers(models["hitemp"], FP_T_FLOOR_K, FP_GRID_FLOOR_SI)
    assert layers, (
        "HITEMP no longer has a radiative layer inside the FP box at the candidate — "
        "the candidate has stopped being a demonstration and must be re-chosen (A42)"
    )
    assert max(x["t_max_k"] for x in layers) >= FP_T_FLOOR_K
    assert min(x["rho_min_si"] for x in layers) >= FP_GRID_FLOOR_SI


def test_composite_closes_it_with_interpolated_opacity(models):
    """A44: zero in-box radiative layers, and the FP κ used is not extrapolated."""
    res = models["composite"]
    layers = radiative_layers(res, FP_T_FLOOR_K, FP_GRID_FLOOR_SI)
    assert layers == [], (
        f"the composite still has {len(layers)} in-box radiative layer(s) at the "
        f"candidate: {layers}"
    )
    assert res.clamp_fraction_all < CLAMP_LIMIT, (
        f"composite clamp_fraction_all = {res.clamp_fraction_all:.4f} exceeds "
        f"{CLAMP_LIMIT:.0%}: the closure is no longer made of interpolated FP data"
    )
    assert res.clamp_reliable


def test_age_grid_problem_flags_a_bad_age_integral():
    """A38: negative or non-monotone ages are reported, not silently used."""
    t_int = np.array([50.0, 100.0, 200.0])
    good = CoolingTrack(
        source="x", t_eq_k=700.0, f_env=0.05, band_cm=None, t_int_k=t_int,
        entropy_j_kg_k=np.zeros(3), radius_p_rearth=np.ones(3),
        luminosity_w=np.ones(3), energy_j=np.ones(3), rcb_bar=np.ones(3),
        age_gyr=np.array([3.0, 1.0, 0.0]), energy_residual=0.0,
    )
    assert age_grid_problem(good) is None

    negative = CoolingTrack(**{**good.__dict__, "age_gyr": np.array([3.0, -1.0, 0.0])})
    assert "negative" in age_grid_problem(negative)

    non_monotone = CoolingTrack(**{**good.__dict__, "age_gyr": np.array([1.0, 3.0, 0.0])})
    assert "monotone" in age_grid_problem(non_monotone)
