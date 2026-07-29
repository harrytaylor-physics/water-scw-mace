"""Envelope integration: Guillot outer region + interior (SPEC §5).

Two stitched integrations for a trial planet radius R_p:

  Outer (P_top → P_match), independent variable ln P, state (r, m, τ):
    T = Guillot(τ); the thermal optical depth grows as dτ = κ_R dP / g.
  Interior (m_match → M_c), independent variable mass m (SPEC §5 form):
    dP/dm = −Gm/(4πr⁴), dr/dm = 1/(4πr²ρ),
    dlnT/dlnP = min(grad_ad, grad_rad) [Schwarzschild],
    grad_rad = 3 κ_R L P / (16π G m σ T⁴), L = L_int constant.

All SI. Constants from ``constants.py``. Tolerances come from config.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp

from constants import G, SIGMA_SB
from eos.base import EOS
from opacity.source import OpacitySource
from structure.boundary import guillot_temperature, t_irr_from_t_eq


@dataclass
class MatchState:
    """State at the interior/atmosphere matching pressure."""

    pressure_pa: float
    radius_m: float
    mass_kg: float
    temperature_k: float
    tau: float


@dataclass
class EnvelopeProfile:
    """Sampled envelope profile (outer + interior concatenated, top → core)."""

    mass_kg: np.ndarray
    radius_m: np.ndarray
    pressure_pa: np.ndarray
    temperature_k: np.ndarray
    rho_si: np.ndarray
    grad: np.ndarray  # actual dlnT/dlnP used
    grad_ad: np.ndarray
    grad_rad: np.ndarray
    is_convective: np.ndarray  # grad_ad < grad_rad
    two_phase: np.ndarray


def luminosity_int(radius_m: float, t_int_k: float) -> float:
    """Internal luminosity L_int = 4π R_p² σ T_int⁴ (SPEC §5)."""
    return 4.0 * np.pi * radius_m**2 * SIGMA_SB * t_int_k**4


# --------------------------------------------------------------------------
# Outer (Guillot) region
# --------------------------------------------------------------------------
def _outer_rhs(ln_p, y, eos, opacity, t_int, t_irr, gamma):
    r, m, tau = y
    p = np.exp(ln_p)
    temperature = float(guillot_temperature(tau, t_int, t_irr, gamma))
    rho = eos.rho(p, temperature)
    g = G * m / r**2
    kappa = opacity.kappa_R(temperature, rho)
    # Hydrostatic in ln P: dr/dlnP = P·dr/dP = −P/(ρg); dm/dlnP = −4π r⁴ P/(Gm).
    dr = -p / (rho * g)
    dm = -4.0 * np.pi * r**4 * p / (G * m)
    # Grey thermal optical depth: dτ/dP = κ_R/g ⇒ dτ/dlnP = P κ_R/g.
    dtau = p * kappa / g
    return [dr, dm, dtau]


def integrate_outer(
    eos: EOS,
    opacity: OpacitySource,
    radius_top_m: float,
    mass_p_kg: float,
    t_int_k: float,
    t_eq_k: float,
    gamma: float,
    p_top_pa: float,
    p_match_pa: float,
    rtol: float,
    atol: float,
    solver: str,
    n_sample: int | None = None,
) -> tuple[MatchState, dict]:
    """Integrate the Guillot outer region from P_top to P_match.

    If ``n_sample`` is given, the profile is returned on that many log-spaced
    pressure points (for figures); otherwise at the solver's own steps.
    """
    t_irr = t_irr_from_t_eq(t_eq_k)
    y0 = [radius_top_m, mass_p_kg, 0.0]
    t_eval = (
        np.linspace(np.log(p_top_pa), np.log(p_match_pa), n_sample) if n_sample else None
    )
    sol = solve_ivp(
        _outer_rhs,
        (np.log(p_top_pa), np.log(p_match_pa)),
        y0,
        method=solver,
        rtol=rtol,
        atol=atol,
        t_eval=t_eval,
        args=(eos, opacity, t_int_k, t_irr, gamma),
    )
    if not sol.success:
        raise RuntimeError(f"outer integration failed: {sol.message}")
    r_m, m_m, tau_m = sol.y[:, -1]
    T_m = float(guillot_temperature(tau_m, t_int_k, t_irr, gamma))
    ln_p = sol.t
    match = MatchState(pressure_pa=p_match_pa, radius_m=r_m, mass_kg=m_m, temperature_k=T_m, tau=tau_m)

    # sampled outer arrays (top -> match)
    P = np.exp(ln_p)
    T = guillot_temperature(sol.y[2], t_int_k, t_irr, gamma)
    prof = {"P": P, "r": sol.y[0], "m": sol.y[1], "T": np.asarray(T), "tau": sol.y[2]}
    return match, prof


# --------------------------------------------------------------------------
# Interior region
# --------------------------------------------------------------------------
def _interior_rhs(m, y, eos, opacity, luminosity):
    r, p, ln_t = y
    temperature = np.exp(ln_t)
    st = eos.state(p, temperature)
    rho = st.rho
    dP_dm = -G * m / (4.0 * np.pi * r**4)
    dr_dm = 1.0 / (4.0 * np.pi * r**2 * rho)
    kappa = opacity.kappa_R(temperature, rho)
    grad_rad = 3.0 * kappa * luminosity * p / (16.0 * np.pi * G * m * SIGMA_SB * temperature**4)
    grad = min(st.grad_ad, grad_rad)
    dlnT_dm = grad * dP_dm / p
    return [dr_dm, dP_dm, dlnT_dm]


def integrate_interior(
    eos: EOS,
    opacity: OpacitySource,
    match: MatchState,
    mass_core_kg: float,
    radius_core_m: float,
    luminosity: float,
    rtol: float,
    atol: float,
    solver: str,
    n_sample: int | None = None,
) -> tuple[float, dict]:
    """Integrate inward in mass from the matching state to the core mass.

    Returns (r_at_core, sampled interior profile). Raises if the envelope
    collapses (r → 0) before the core mass is reached. If ``n_sample`` is given,
    the profile is returned on a mass grid clustered toward the matching point
    (where pressure varies fastest).
    """
    y0 = [match.radius_m, match.pressure_pa, np.log(match.temperature_k)]

    def collapse(m, y, *a):
        return y[0] - 0.05 * radius_core_m

    collapse.terminal = True
    collapse.direction = -1

    t_eval = None
    if n_sample:
        frac = np.linspace(0.0, 1.0, n_sample) ** 2  # cluster near m_match
        t_eval = match.mass_kg - (match.mass_kg - mass_core_kg) * frac

    sol = solve_ivp(
        _interior_rhs,
        (match.mass_kg, mass_core_kg),
        y0,
        method=solver,
        rtol=rtol,
        atol=atol,
        t_eval=t_eval,
        events=collapse,
        args=(eos, opacity, luminosity),
    )
    if not sol.success:
        raise RuntimeError(f"interior integration failed: {sol.message}")
    if sol.t_events[0].size > 0:
        # collapsed before reaching the core -> steer shooter with r -> small
        return float(sol.y[0, -1]), {"collapsed": True}

    r_core = float(sol.y[0, -1])
    prof = {
        "m": sol.t,
        "r": sol.y[0],
        "P": sol.y[1],
        "T": np.exp(sol.y[2]),
        "collapsed": False,
    }
    return r_core, prof
