"""High-level static structure model (SPEC §5).

``build_structure`` assembles a converged interior model for a planet at a
given (T_eq, T_int, opacity source): shoots on R_p, samples the full profile
(outer Guillot region + interior), post-processes derived quantities, extracts
diagnostics (R_p, RCB location(s), photosphere, opacity-clamp fraction,
two-phase fraction, hydrostatic residual), and serializes to npz + json.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from pathlib import Path

import numpy as np

from config import Config
from constants import G, R_EARTH, SIGMA_SB
from eos import AquaEOS, core_mass_si, core_radius_si
from opacity import get_source
from structure.integrator import integrate_interior, integrate_outer, luminosity_int
from structure.shooting import shoot_radius


@dataclass
class StructureResult:
    """Converged model summary + full profile arrays (SI)."""

    name: str
    source: str
    t_eq_k: float
    t_int_k: float
    gamma: float
    mass_core_mearth: float
    f_env: float
    match_pressure_bar: float

    radius_p_m: float
    radius_p_rearth: float
    radius_core_m: float
    mass_p_kg: float
    mass_env_kg: float
    luminosity_w: float

    rcb_pressures_pa: list[float]
    rcb_above_match: bool  # Schwarzschild RCB shallower than P_match (matching caveat)
    photosphere_pressure_pa: float
    photosphere_temperature_k: float

    clamp_fraction: float
    clamp_reliable: bool
    two_phase_fraction: float
    hydrostatic_residual: float
    shoot_iterations: int
    shoot_residual_m: float

    # profile arrays (top -> core), excluded from the json summary
    profile: dict = field(default_factory=dict, repr=False)

    def summary(self) -> dict:
        d = asdict(self)
        d.pop("profile")
        return d


def _post_process(m, r, p, t, eos, opacity, luminosity):
    """Derived quantities along the profile (looped scalar EOS/opacity)."""
    n = m.size
    rho = np.empty(n)
    grad_ad = np.empty(n)
    grad_rad = np.empty(n)
    two_phase = np.zeros(n, dtype=bool)
    for i in range(n):
        st = eos.state(p[i], t[i])
        rho[i] = st.rho
        grad_ad[i] = st.grad_ad
        two_phase[i] = st.two_phase
        kappa = opacity.kappa_R(t[i], st.rho)
        grad_rad[i] = 3.0 * kappa * luminosity * p[i] / (
            16.0 * np.pi * G * m[i] * SIGMA_SB * t[i] ** 4
        )
    is_conv = grad_rad > grad_ad
    return rho, grad_ad, grad_rad, is_conv, two_phase


def _find_rcbs(pressure, is_convective):
    """Pressures where the profile switches convective<->radiative."""
    flips = np.nonzero(np.diff(is_convective.astype(int)) != 0)[0]
    # geometric mean of the bracketing pressures at each flip
    return [float(np.sqrt(pressure[i] * pressure[i + 1])) for i in flips]


def _photosphere(pressure, tau, temperature):
    """(P, T) where τ = 2/3, interpolated in log P; fallback to the top."""
    target = 2.0 / 3.0
    if tau[-1] < target:
        return float(pressure[0]), float(temperature[0])
    i = int(np.searchsorted(tau, target))
    i = min(max(i, 1), tau.size - 1)
    lo, hi = i - 1, i
    w = (target - tau[lo]) / (tau[hi] - tau[lo])
    logp = (1 - w) * np.log(pressure[lo]) + w * np.log(pressure[hi])
    T = (1 - w) * temperature[lo] + w * temperature[hi]
    return float(np.exp(logp)), float(T)


def _hydrostatic_residual(m, r, p):
    """Relative error of ∫(−Gm/4πr⁴)dm vs the actual ΔP over the interior (G3)."""
    rhs = -G * m / (4.0 * np.pi * r**4)
    predicted = np.trapz(rhs, m)
    actual = p[-1] - p[0]
    return abs(predicted - actual) / abs(actual)


def build_structure(
    cfg: Config,
    t_eq_k: float,
    t_int_k: float,
    source_name: str | None = None,
    eos: AquaEOS | None = None,
    n_sample: int = 2000,
) -> StructureResult:
    """Build one converged static structure model."""
    source_name = source_name or cfg.opacity.source
    eos = eos or AquaEOS.from_config(cfg)
    opacity = get_source(source_name, cfg)

    mc_mearth = cfg.core.mass_mearth
    f_env = cfg.envelope.f_env
    mp_mearth = mc_mearth / (1.0 - f_env)
    mass_p = core_mass_si(mp_mearth)
    mass_core = core_mass_si(mc_mearth)
    radius_core = core_radius_si(mc_mearth)

    gamma = cfg.irradiation.gamma
    p_top = cfg.boundary.p_top_bar * 1e5
    p_match = cfg.boundary.match_pressure_bar * 1e5
    num = cfg.numerics

    shot = shoot_radius(
        eos, opacity, mass_p, mass_core, radius_core, t_int_k, t_eq_k, gamma,
        p_top, p_match, num.rtol, num.atol, num.solver,
        (num.rp_bracket_rearth[0] * R_EARTH, num.rp_bracket_rearth[1] * R_EARTH),
        num.shoot_rtol, num.max_shoot_iter,
    )

    # Final sampled profile at the converged radius; measure clamp fraction here.
    opacity.reset_clamp_counter()
    match, outer = integrate_outer(
        eos, opacity, shot.radius_p_m, mass_p, t_int_k, t_eq_k, gamma,
        p_top, p_match, num.rtol, num.atol, num.solver, n_sample=n_sample,
    )
    lum = luminosity_int(shot.radius_p_m, t_int_k)
    _, interior = integrate_interior(
        eos, opacity, match, mass_core, radius_core, lum,
        num.rtol, num.atol, num.solver, n_sample=n_sample,
    )
    clamp_fraction = opacity.clamp_fraction

    # Concatenate top -> core (drop duplicate matching point).
    m = np.concatenate([outer["m"], interior["m"][1:]])
    r = np.concatenate([outer["r"], interior["r"][1:]])
    p = np.concatenate([outer["P"], interior["P"][1:]])
    t = np.concatenate([outer["T"], interior["T"][1:]])
    tau_outer = outer["tau"]

    rho, grad_ad, grad_rad, is_conv, two_phase = _post_process(m, r, p, t, eos, opacity, lum)
    grad = np.gradient(np.log(t), np.log(p))

    rcbs = _find_rcbs(p, is_conv)
    p_phot, t_phot = _photosphere(outer["P"], tau_outer, outer["T"])
    hydro = _hydrostatic_residual(interior["m"], interior["r"], interior["P"])
    threshold = cfg.opacity.clamp_fraction_threshold

    result = StructureResult(
        name=cfg.name,
        source=source_name,
        t_eq_k=t_eq_k,
        t_int_k=t_int_k,
        gamma=gamma,
        mass_core_mearth=mc_mearth,
        f_env=f_env,
        match_pressure_bar=cfg.boundary.match_pressure_bar,
        radius_p_m=shot.radius_p_m,
        radius_p_rearth=shot.radius_p_m / R_EARTH,
        radius_core_m=radius_core,
        mass_p_kg=mass_p,
        mass_env_kg=mass_p - mass_core,
        luminosity_w=lum,
        rcb_pressures_pa=rcbs,
        rcb_above_match=any(rp < p_match for rp in rcbs),
        photosphere_pressure_pa=p_phot,
        photosphere_temperature_k=t_phot,
        clamp_fraction=clamp_fraction,
        clamp_reliable=clamp_fraction <= threshold,
        two_phase_fraction=float(two_phase.mean()),
        hydrostatic_residual=hydro,
        shoot_iterations=shot.n_iter,
        shoot_residual_m=shot.residual_m,
        profile={
            "mass_kg": m, "radius_m": r, "pressure_pa": p, "temperature_k": t,
            "rho_si": rho, "grad": grad, "grad_ad": grad_ad, "grad_rad": grad_rad,
            "is_convective": is_conv, "two_phase": two_phase,
        },
    )
    return result


def save_structure(result: StructureResult, out_dir: Path) -> tuple[Path, Path]:
    """Write profile npz + summary json; return their paths."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    stem = f"{result.name}_{result.source}_Teq{result.t_eq_k:.0f}_Tint{result.t_int_k:.0f}"
    npz_path = out_dir / f"{stem}.npz"
    json_path = out_dir / f"{stem}.json"
    np.savez_compressed(npz_path, **result.profile)
    with json_path.open("w") as fh:
        json.dump(result.summary(), fh, indent=2)
    return npz_path, json_path
