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
# Runaway guard (A57)
# --------------------------------------------------------------------------
class _EvalBudget:
    """Cap the RHS evaluations one ``solve_ivp`` call may spend.

    LSODA can grind indefinitely on a stiff or ill-posed segment: it keeps
    halving its step instead of failing, so the call never returns and the whole
    run hangs with no diagnostic. That cost Task 18 two multi-hour runs before it
    was caught (A57). A hang is strictly worse than an error here — every caller
    already handles a failed build (the domain scan skips and logs it, the
    shooting solver re-brackets) — so the budget converts one into the other.

    The cap is generous: a healthy integration uses well under 1000 evaluations,
    so tripping it means the solver is not converging. Note that
    ``structure.shooting`` swallows the error on its warm-start attempt and
    falls back to the bracketed solve, so one bad node costs seconds, not a run.
    """

    __slots__ = ("limit", "used", "where")

    def __init__(self, limit: int, where: str) -> None:
        self.limit = int(limit)
        self.used = 0
        self.where = where

    def spend(self) -> None:
        self.used += 1
        if self.used > self.limit:
            raise RuntimeError(
                f"{self.where}: exceeded {self.limit} right-hand-side evaluations "
                f"without converging — the integration is not making progress. "
                f"Raise numerics.max_rhs_evals only if a legitimate model needs "
                f"more; a hang here used to be silent (A57)."
            )


# --------------------------------------------------------------------------
# Outer (Guillot) region
# --------------------------------------------------------------------------
def _outer_rhs(ln_p, y, eos, opacity, t_int, t_irr, gamma, budget=None):
    if budget is not None:
        budget.spend()
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


def _outer_superadiabatic(ln_p, y, eos, opacity, t_int, t_irr, gamma,
                          budget=None):
    """grad_rad − grad_ad in the outer region (Schwarzschild RCB event, A13).

    L_int is 4π R_p² σ T_int⁴ with R_p the top radius (radius_top_m enters via
    the initial condition; here m ≈ M_p so we use the local m and the top r
    through L is passed implicitly by evaluating grad_rad with the local state).
    """
    r, m, tau = y
    p = np.exp(ln_p)
    temperature = float(guillot_temperature(tau, t_int, t_irr, gamma))
    st = eos.state(p, temperature)
    kappa = opacity.kappa_R(temperature, st.rho)
    lum = luminosity_int(y[0], t_int)  # L uses local radius (≈ R_p in the thin outer layer)
    grad_rad = 3.0 * kappa * lum * p / (16.0 * np.pi * G * m * SIGMA_SB * temperature**4)
    return grad_rad - st.grad_ad


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
    max_rhs_evals: int = 200_000,
) -> tuple[MatchState, dict, bool]:
    """Integrate the Guillot outer region from P_top down toward P_match.

    Adaptive matching (A13): the integration stops at the radiative–convective
    boundary (grad_rad = grad_ad) if that is reached before P_match, so the
    analytic radiative profile is only used where the atmosphere is radiative.
    Returns (match_state, sampled_profile, matched_at_rcb). If ``n_sample`` is
    given the profile is returned on that many log-spaced points (via dense
    output) up to the actual matching pressure.
    """
    t_irr = t_irr_from_t_eq(t_eq_k)
    y0 = [radius_top_m, mass_p_kg, 0.0]
    budget = _EvalBudget(max_rhs_evals, "outer (Guillot) integration")
    args = (eos, opacity, t_int_k, t_irr, gamma)

    event = lambda lp, y, *a: _outer_superadiabatic(lp, y, *a)  # noqa: E731
    event.terminal = True
    event.direction = 1.0  # radiative (−) -> convective (+)

    # If the very top is already convective, match immediately at P_top.
    if _outer_superadiabatic(np.log(p_top_pa), y0, *args) >= 0.0:
        T0 = float(guillot_temperature(0.0, t_int_k, t_irr, gamma))
        match = MatchState(p_top_pa, radius_top_m, mass_p_kg, T0, 0.0)
        prof = {"P": np.array([p_top_pa]), "r": np.array([radius_top_m]),
                "m": np.array([mass_p_kg]), "T": np.array([T0]), "tau": np.array([0.0])}
        return match, prof, True

    sol = solve_ivp(
        _outer_rhs, (np.log(p_top_pa), np.log(p_match_pa)), y0,
        method=solver, rtol=rtol, atol=atol, events=event,
        dense_output=bool(n_sample), args=(*args, budget),
    )
    if not sol.success:
        raise RuntimeError(f"outer integration failed: {sol.message}")

    matched_at_rcb = sol.t_events[0].size > 0
    if matched_at_rcb:
        ln_p_end = float(sol.t_events[0][0])
        r_m, m_m, tau_m = sol.y_events[0][0]
    else:
        ln_p_end = float(sol.t[-1])
        r_m, m_m, tau_m = sol.y[:, -1]
    p_end = float(np.exp(ln_p_end))
    T_m = float(guillot_temperature(tau_m, t_int_k, t_irr, gamma))
    match = MatchState(pressure_pa=p_end, radius_m=r_m, mass_kg=m_m, temperature_k=T_m, tau=tau_m)

    if n_sample:
        ln_grid = np.linspace(np.log(p_top_pa), ln_p_end, n_sample)
        y = sol.sol(ln_grid)
        P, r, m, tau = np.exp(ln_grid), y[0], y[1], y[2]
    else:
        keep = sol.t <= ln_p_end
        P, r, m, tau = np.exp(sol.t[keep]), sol.y[0][keep], sol.y[1][keep], sol.y[2][keep]
    T = np.asarray(guillot_temperature(tau, t_int_k, t_irr, gamma))
    prof = {"P": P, "r": r, "m": m, "T": T, "tau": tau}
    return match, prof, matched_at_rcb


# --------------------------------------------------------------------------
# Interior region
# --------------------------------------------------------------------------
def _interior_rhs(m, y, eos, opacity, luminosity, budget=None):
    if budget is not None:
        budget.spend()
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
    max_rhs_evals: int = 200_000,
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
        args=(eos, opacity, luminosity,
              _EvalBudget(max_rhs_evals, "interior integration")),
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
