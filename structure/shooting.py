"""Shooting on R_p until r(M_c) = R_c (SPEC §5).

For a trial planet radius R_p we integrate the Guillot outer region then the
interior inward in mass; the residual r(M_c) − R_c is driven to zero with a
bracketed root find (Brent). L_int depends on R_p, so it is recomputed each
trial. Relaxation is the documented alternative (future work).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq

from constants import R_EARTH
from eos.base import EOS
from opacity.source import OpacitySource
from structure.integrator import (
    MatchState,
    integrate_interior,
    integrate_outer,
    luminosity_int,
)


@dataclass
class ShootResult:
    radius_p_m: float
    match: MatchState
    outer_profile: dict
    interior_profile: dict
    luminosity: float
    n_iter: int
    residual_m: float


def _residual(
    radius_p_m, eos, opacity, mass_p, mass_core, radius_core, t_int, t_eq, gamma,
    p_top, p_match, rtol, atol, solver,
):
    """r(M_c) − R_c for a trial R_p (returns the residual and the integration state)."""
    match, outer, _ = integrate_outer(
        eos, opacity, radius_p_m, mass_p, t_int, t_eq, gamma, p_top, p_match, rtol, atol, solver
    )
    lum = luminosity_int(radius_p_m, t_int)
    r_core, interior = integrate_interior(
        eos, opacity, match, mass_core, radius_core, lum, rtol, atol, solver
    )
    return r_core - radius_core, match, outer, interior, lum


def shoot_radius(
    eos: EOS,
    opacity: OpacitySource,
    mass_p_kg: float,
    mass_core_kg: float,
    radius_core_m: float,
    t_int_k: float,
    t_eq_k: float,
    gamma: float,
    p_top_pa: float,
    p_match_pa: float,
    rtol: float,
    atol: float,
    solver: str,
    rp_bracket_m: tuple[float, float],
    shoot_rtol: float,
    max_iter: int,
) -> ShootResult:
    """Find R_p such that r(M_c) = R_c to relative tolerance ``shoot_rtol``."""
    lo = max(rp_bracket_m[0], radius_core_m * 1.001)
    hi = rp_bracket_m[1]

    calls = {"n": 0}
    cache: dict[float, tuple] = {}

    def resid(rp):
        calls["n"] += 1
        out = _residual(
            rp, eos, opacity, mass_p_kg, mass_core_kg, radius_core_m, t_int_k, t_eq_k,
            gamma, p_top_pa, p_match_pa, rtol, atol, solver,
        )
        cache[rp] = out
        return out[0]

    # Ensure the bracket straddles a sign change; expand/scan if not.
    f_lo, f_hi = resid(lo), resid(hi)
    if f_lo * f_hi > 0:
        grid = np.linspace(lo, hi, 12)
        vals = [(g, resid(g)) for g in grid]
        bracket = None
        for (a, fa), (b, fb) in zip(vals, vals[1:]):
            if fa * fb <= 0:
                bracket = (a, b)
                break
        if bracket is None:
            raise RuntimeError(
                f"no sign change in residual over R_p ∈ [{lo/R_EARTH:.2f}, {hi/R_EARTH:.2f}] R⊕; "
                f"residuals {[f'{v:.3e}' for _, v in vals]}"
            )
        lo, hi = bracket

    rp_star, r = brentq(
        resid, lo, hi, xtol=shoot_rtol * radius_core_m, rtol=shoot_rtol,
        maxiter=max_iter, full_output=True,
    )
    residual, match, outer, interior, lum = cache[rp_star]
    return ShootResult(
        radius_p_m=rp_star,
        match=match,
        outer_profile=outer,
        interior_profile=interior,
        luminosity=lum,
        n_iter=calls["n"],
        residual_m=residual,
    )
