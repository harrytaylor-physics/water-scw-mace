"""Guillot (2010) irradiated grey outer boundary (SPEC §5).

Guillot, T. (2010), A&A 520, A27. Analytic two-stream temperature profile for
an irradiated grey atmosphere with an internal flux. We use the collimated
form (his eq. 49) with a planet-average incidence angle μ⋆ = 1/√3 (the
Eddington isotropic-equivalent), and interpret the config equilibrium
temperature T_eq as the full-redistribution value, T_irr = √2·T_eq
(σT_eq⁴ = ¼ σT_irr⁴). Derivation, conventions, and the T_eq↔T_irr choice are in
`docs/boundary_guillot.md` (AMENDMENTS A10).

The thermal optical depth τ (from the top, increasing inward) is the argument.
γ = κ_visible / κ_thermal is the config knob (fiducial + sensitivity, SPEC §10).
"""

from __future__ import annotations

import numpy as np

# Planet-averaged incidence angle (Eddington isotropic equivalent).
MU_STAR_PLANET_AVG: float = 1.0 / np.sqrt(3.0)


def t_irr_from_t_eq(t_eq_k: float) -> float:
    """Irradiation temperature from full-redistribution equilibrium temperature.

    σT_eq⁴ = ¼ σT_irr⁴  ⇒  T_irr = 4^(1/4) T_eq = √2 · T_eq.
    """
    return 4.0**0.25 * t_eq_k


def guillot_temperature4(
    tau: np.ndarray,
    t_int_k: float,
    t_irr_k: float,
    gamma: float,
    mu_star: float = MU_STAR_PLANET_AVG,
) -> np.ndarray:
    """T⁴(τ) from Guillot (2010) eq. 49 (collimated, incidence μ⋆).

    Internal-flux term + irradiation term. Returns T⁴ [K⁴].
    """
    tau = np.asarray(tau, dtype=float)
    internal = 0.75 * t_int_k**4 * (tau + 2.0 / 3.0)
    e = np.exp(-gamma * tau / mu_star)
    irr = (
        0.75
        * t_irr_k**4
        * mu_star
        * (2.0 / 3.0 + mu_star / gamma + (gamma / (3.0 * mu_star) - mu_star / gamma) * e)
    )
    return internal + irr


def guillot_temperature(
    tau, t_int_k: float, t_irr_k: float, gamma: float, mu_star: float = MU_STAR_PLANET_AVG
):
    """T(τ) [K] from Guillot (2010)."""
    return guillot_temperature4(tau, t_int_k, t_irr_k, gamma, mu_star) ** 0.25


def skin_temperature(
    t_int_k: float, t_eq_k: float, gamma: float, mu_star: float = MU_STAR_PLANET_AVG
) -> float:
    """Skin temperature T(τ→0): the isothermal limit as κ → 0 (gate G1)."""
    t_irr = t_irr_from_t_eq(t_eq_k)
    return float(guillot_temperature(0.0, t_int_k, t_irr, gamma, mu_star))
