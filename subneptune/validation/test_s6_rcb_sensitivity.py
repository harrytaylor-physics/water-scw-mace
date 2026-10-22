"""How far each RCB is from moving, in units of opacity (A49).

Before any new first-principles opacity exists, the useful question is how big a
change in κ_R a radiative-convective boundary would even notice. The model's own
criterion is ∇_rad > ∇_ad with ∇_rad ∝ κ_R, so at fixed structure the answer is
exactly ∇_ad/∇_rad at the sample in question. These gates pin that arithmetic and
its agreement with the RCBs the structure model reports, so the printed factors
are the model's, not a second implementation of it.
"""

from __future__ import annotations

import numpy as np
import pytest

from eos import AquaEOS
from opacity import get_source
from structure import build_structure
from structure.diagnostics import rcb_opacity_sensitivity


@pytest.fixture(scope="module")
def eos(cfg):
    return AquaEOS.from_config(cfg)


@pytest.fixture(scope="module")
def hitemp_cold(cfg, eos):
    """The A41 fiducial case: one RCB with a thin radiative layer above it."""
    return build_structure(cfg, 500.0, 20.0, "hitemp", eos=eos, n_sample=600)


@pytest.fixture(scope="module")
def sens(cfg, hitemp_cold):
    return rcb_opacity_sensitivity(hitemp_cold, get_source("hitemp", cfg))


def test_rcbs_agree_with_the_structure_model(hitemp_cold, sens):
    """Same boundaries, same geometric-mean convention as StructureResult."""
    assert len(sens["rcbs"]) == len(hitemp_cold.rcb_densities_si)
    for got, rho, p_pa in zip(sens["rcbs"], hitemp_cold.rcb_densities_si,
                              hitemp_cold.rcb_pressures_pa):
        assert got["rho_si"] == pytest.approx(rho)
        assert got["p_pa"] == pytest.approx(p_pa)
        assert got["t_k"] > 0.0
        assert got["kappa_r_si"] > 0.0


def test_the_factor_is_the_schwarzschild_ratio(cfg, sens):
    """κ_R × factor_close makes ∇_rad = ∇_ad at the layer's deepest sample."""
    assert sens["layers"], "the fiducial should have a radiative layer to test"
    for layer in sens["layers"]:
        base = layer["base"]
        # ∇_rad ∝ κ_R at fixed structure, so scaling κ scales ∇_rad exactly.
        assert base["grad_rad"] * base["factor_close"] == pytest.approx(base["grad_ad"])
        # A radiative sample needs MORE opacity to convect.
        assert base["factor_close"] >= 1.0
        assert base["grad_rad"] <= base["grad_ad"]


def test_extending_a_layer_needs_less_opacity(sens):
    """The first convective sample below flips the other way: factor < 1."""
    checked = 0
    for layer in sens["layers"]:
        ext = layer["extend"]
        if ext is None:  # layer runs to the core
            continue
        checked += 1
        assert ext["grad_rad"] > ext["grad_ad"]  # it is convective now
        assert ext["factor_extend"] <= 1.0
        assert ext["grad_rad"] * ext["factor_extend"] == pytest.approx(ext["grad_ad"])
    assert checked, "no layer had a convective sample below it"


def test_the_reported_opacity_is_the_source_evaluated_there(cfg, sens):
    """κ_R is read from the source at the sample, not interpolated off a table."""
    source = get_source("hitemp", cfg)
    for layer in sens["layers"]:
        base = layer["base"]
        assert base["kappa_r_si"] == pytest.approx(
            source.kappa_R(base["t_k"], base["rho_si"]))


def test_a_fully_convective_profile_has_no_layers(cfg, hitemp_cold):
    """Degenerate case: no radiative sample, no boundary, no factor."""
    class AllConvective:
        profile = dict(hitemp_cold.profile)
        profile["is_convective"] = np.ones_like(hitemp_cold.profile["is_convective"],
                                                dtype=bool)

    out = rcb_opacity_sensitivity(AllConvective(), get_source("hitemp", cfg))
    assert out["layers"] == []
    assert out["rcbs"] == []


def test_layer_masses_match_the_radiative_layer_accounting(hitemp_cold, sens):
    """The same mass weights as radiative_layers (A42), so the two agree."""
    from structure.diagnostics import radiative_layers

    layers = radiative_layers(hitemp_cold)
    assert len(layers) == len(sens["layers"])
    for a, b in zip(layers, sens["layers"]):
        assert a["mass_fraction"] == pytest.approx(b["mass_fraction"])
        assert a["t_min_k"] == pytest.approx(b["t_min_k"])
        assert a["rho_max_si"] == pytest.approx(b["rho_max_si"])
