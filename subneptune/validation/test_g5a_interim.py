"""Interim analytic anchor G5a (AMENDMENTS A20), pending the Aguichine data.

The real G5 (published steam-world M–R vs Aguichine 2021) needs an external
table that is not on the machine, so it skips. In the meantime these two
*derivable* analytic limits validate the structure integration end-to-end with
no external data:

  G5a-1  Bare-core limit: R_p → R_core (Zeng 2016) as f_env → 0, monotonically
         and from above; the residual is the radiative atmosphere skin.
  G5a-2  Isothermal scale height: the radiative-zone thickness matches the
         closed-form ΔR = H·ln(P_base/P_top), H = k_B T/(μ_H2O_molec · g), to
         ~10% — validating hydrostatic balance + EOS in the ideal-gas zone.

These are analytic limits (like G1), NOT a literature M–R anchor; they do not
replace G5.
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest

from constants import G, K_B, MMW_H2O, N_A, R_EARTH
from eos import AquaEOS, core_radius_si
from figures._common import register_constant_reference
from structure import build_structure


@pytest.fixture(scope="module")
def eos(cfg):
    return AquaEOS.from_config(cfg)


def test_g5a1_bare_core_limit(cfg, eos):
    """R_p → R_core as f_env → 0: monotone, from above, small residual."""
    register_constant_reference(cfg)
    r_core = core_radius_si(cfg.core.mass_mearth) / R_EARTH

    radii = []
    for f in (5e-2, 1e-2, 1e-3, 1e-4):
        c = replace(cfg, envelope=replace(cfg.envelope, f_env=f))
        radii.append(build_structure(c, 600.0, 50.0, "constant", eos=eos).radius_p_rearth)

    assert all(rp > r_core for rp in radii)  # envelope always adds radius
    assert all(a > b for a, b in zip(radii, radii[1:]))  # monotone decreasing
    # At f_env=1e-4 only a thin atmosphere skin remains (< 0.1 R⊕).
    assert radii[-1] - r_core < 0.1


def test_g5a2_isothermal_scale_height(cfg, eos):
    """Radiative-zone thickness ≈ analytic isothermal scale height (~10%)."""
    register_constant_reference(cfg)
    c = replace(cfg, envelope=replace(cfg.envelope, f_env=1e-2))
    res = build_structure(c, 600.0, 30.0, "constant", eos=eos)
    prof = res.profile

    p = prof["pressure_pa"]
    r = prof["radius_m"]
    p_top, t_top = p[0], prof["temperature_k"][0]
    p_base = max(res.rcb_pressures_pa)  # radiative zone runs top → RCB

    g_surf = G * res.mass_p_kg / res.radius_p_m**2
    mu_molecule = MMW_H2O / N_A  # kg per molecule
    scale_height = K_B * t_top / (mu_molecule * g_surf)
    dr_analytic = scale_height * np.log(p_base / p_top)

    i_rcb = int(np.searchsorted(p, p_base))
    r_rcb = r[min(i_rcb, len(r) - 1)]
    dr_measured = r[0] - r_rcb

    assert dr_measured / dr_analytic == pytest.approx(1.0, abs=0.12)
