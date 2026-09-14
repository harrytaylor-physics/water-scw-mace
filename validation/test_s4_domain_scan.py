"""Radiative-layer bookkeeping behind the domain scan (A42).

The scan's conclusion — that a line-list radiative zone almost never reaches
inside the FP data box (T ≥ 1000 K, ρ ≥ 108.796 kg/m³) — is only as good as the
layer accounting it rests on. These gates check that accounting: mass weights
that sum to the envelope, layers that partition the radiative samples, and a
box filter that actually restricts.
"""

from __future__ import annotations

import numpy as np
import pytest

from eos import AquaEOS
from structure import build_structure, deep_radiative_zones
from structure.diagnostics import (
    _sample_mass_weights,
    in_box_radiative_mass_fraction,
    radiative_layers,
)

FP_T_FLOOR_K = 1000.0
FP_RHO_FLOOR_SI = 108.796


@pytest.fixture(scope="module")
def eos(cfg):
    return AquaEOS.from_config(cfg)


@pytest.fixture(scope="module")
def hitemp_cold(cfg, eos):
    """The A41 fiducial case: a deep, cold line-list radiative sliver."""
    return build_structure(cfg, 500.0, 20.0, "hitemp", eos=eos, n_sample=600)


def test_sample_mass_weights_sum_to_the_envelope():
    """Central-difference weights tile the mass coordinate exactly."""
    m = np.array([10.0, 8.0, 6.0, 3.0, 1.0])  # top -> core, decreasing
    w = _sample_mass_weights(m)
    assert w.sum() == pytest.approx(m[0] - m[-1])
    assert np.all(w > 0)


def test_layers_partition_the_radiative_samples(hitemp_cold):
    """Every radiative sample belongs to exactly one layer, and mass adds up."""
    prof = hitemp_cold.profile
    layers = radiative_layers(hitemp_cold)
    n_radiative = int((~prof["is_convective"]).sum())
    assert sum(x["n_samples"] for x in layers) == n_radiative

    w = _sample_mass_weights(prof["mass_kg"])
    total = float(w[~prof["is_convective"]].sum() / w.sum())
    assert sum(x["mass_fraction"] for x in layers) == pytest.approx(total)
    assert 0.0 <= total <= 1.0


def test_in_box_fraction_is_a_subset_of_the_radiative_fraction(hitemp_cold):
    in_box = in_box_radiative_mass_fraction(hitemp_cold, FP_T_FLOOR_K, FP_RHO_FLOOR_SI)
    everything = in_box_radiative_mass_fraction(hitemp_cold, 0.0, 0.0)
    assert 0.0 <= in_box <= everything <= 1.0
    assert sum(x["mass_fraction"]
               for x in radiative_layers(hitemp_cold, FP_T_FLOOR_K, FP_RHO_FLOOR_SI)) == (
        pytest.approx(in_box))


def test_box_filter_restricts(hitemp_cold):
    """An unreachable box leaves nothing; the open box leaves every layer."""
    assert radiative_layers(hitemp_cold, t_min_k=1e6) == []
    assert radiative_layers(hitemp_cold, rho_min_si=1e9) == []
    assert len(radiative_layers(hitemp_cold)) >= 1


def test_layers_agree_with_the_deep_zone_diagnostic(hitemp_cold):
    """``radiative_layers`` and ``deep_radiative_zones`` (A35) see the same runs."""
    deep = deep_radiative_zones(hitemp_cold, 300.0)
    layers = radiative_layers(hitemp_cold, rho_min_si=300.0)
    assert len(deep) == len(layers)
    for (lo, hi), layer in zip(deep, layers):
        assert layer["rho_min_si"] == pytest.approx(lo)
        assert layer["rho_max_si"] == pytest.approx(hi)


def test_cold_fiducial_has_no_in_box_radiative_mass(hitemp_cold):
    """A42: at T_eq=500 K, T_int=20 K the deep sliver is far below 1000 K."""
    assert in_box_radiative_mass_fraction(hitemp_cold, FP_T_FLOOR_K, FP_RHO_FLOOR_SI) == 0.0
    hot = max(x["t_max_k"] for x in radiative_layers(hitemp_cold))
    assert hot < FP_T_FLOOR_K, (
        f"the line-list radiative zone now reaches {hot:.0f} K, inside the FP grid — "
        "the A41/A42 diagnosis has changed and must be re-reported"
    )
