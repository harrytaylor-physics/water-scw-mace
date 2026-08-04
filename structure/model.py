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

# np.trapz was removed in numpy 2.0 (renamed trapezoid); support both.
try:
    from numpy import trapezoid as _trapz
except ImportError:  # numpy < 2.0
    from numpy import trapz as _trapz

from config import Config
from constants import G, R_EARTH, SIGMA_SB
from eos import AquaEOS, core_mass_si, core_radius_si
from opacity import CompositeOpacitySource, get_source
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
    band_cm: tuple[float, float] | None  # opacity band used (A21)
    match_pressure_bar: float  # configured deep cap
    match_pressure_actual_bar: float  # adaptive: min(cap, P_RCB) (A13)
    matched_at_rcb: bool

    radius_p_m: float
    radius_p_rearth: float
    radius_core_m: float
    mass_p_kg: float
    mass_env_kg: float
    luminosity_w: float
    envelope_energy_j: float
    tau_kh_s: float  # Kelvin-Helmholtz time |E|/L_int (SPEC §6, F3)
    interior_entropy_j_kg_k: float  # deep convective adiabat entropy (A23, F3/Phase 2)

    rcb_pressures_pa: list[float]
    rcb_densities_si: list[float]  # ρ at each RCB (A28 scan)
    max_radiative_rho_si: float  # densest radiative point: where opacity still matters
    reaches_fp_grid: bool  # max_radiative_rho ≥ FP grid floor (109 kg/m³) → FP-relevant
    rcb_above_match: bool  # Schwarzschild RCB shallower than P_match (matching caveat)
    photosphere_pressure_pa: float
    photosphere_temperature_k: float

    clamp_fraction: float  # raw fraction of all opacity evals clamped
    clamp_fraction_radiative: float  # clamped fraction in the radiative zone (A26)
    clamp_reliable: bool  # radiative clamp fraction <= threshold
    flat_rho_fraction: float  # dilute-limit ρ-floored fraction (justified; A26)
    two_phase_fraction: float
    hydrostatic_residual: float
    shoot_iterations: int
    shoot_residual_m: float
    composite_seam_ratio: dict | None = None  # FP/line-list κ ratio at seam (A12)

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
    u = np.empty(n)
    entropy = np.empty(n)
    two_phase = np.zeros(n, dtype=bool)
    for i in range(n):
        st = eos.state(p[i], t[i])
        rho[i] = st.rho
        grad_ad[i] = st.grad_ad
        u[i] = st.internal_energy
        entropy[i] = st.entropy
        two_phase[i] = st.two_phase
        kappa = opacity.kappa_R(t[i], st.rho)
        grad_rad[i] = 3.0 * kappa * luminosity * p[i] / (
            16.0 * np.pi * G * m[i] * SIGMA_SB * t[i] ** 4
        )
    is_conv = grad_rad > grad_ad
    return rho, grad_ad, grad_rad, is_conv, two_phase, u, entropy


def _envelope_energy(m, r, u):
    """Total envelope energy E = ∫(u − Gm/r) dm [J] over the envelope.

    Internal + gravitational; the Kelvin–Helmholtz time is |E|/L_int (SPEC §6,
    used for the F3 cooling-time ratio). ``m`` runs top→core (decreasing).
    """
    order = np.argsort(m)
    m_s, r_s, u_s = m[order], r[order], u[order]
    e_int = _trapz(u_s, m_s)
    e_grav = _trapz(-G * m_s / r_s, m_s)
    return float(e_int + e_grav)


def _find_rcbs(pressure, is_convective):
    """Pressures where the profile switches convective<->radiative."""
    flips = np.nonzero(np.diff(is_convective.astype(int)) != 0)[0]
    # geometric mean of the bracketing pressures at each flip
    return [float(np.sqrt(pressure[i] * pressure[i + 1])) for i in flips]


def _rcb_densities(rho, is_convective):
    """ρ at each convective<->radiative flip (geometric mean of neighbours)."""
    flips = np.nonzero(np.diff(is_convective.astype(int)) != 0)[0]
    return [float(np.sqrt(rho[i] * rho[i + 1])) for i in flips]


FP_GRID_FLOOR_SI = 108.796  # FP's least-dense grid point (kg/m³); FP-relevance boundary


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
    predicted = _trapz(rhs, m)
    actual = p[-1] - p[0]
    return abs(predicted - actual) / abs(actual)


def build_structure(
    cfg: Config,
    t_eq_k: float,
    t_int_k: float,
    source_name: str | None = None,
    eos: AquaEOS | None = None,
    n_sample: int = 2000,
    band: tuple[float, float] | None = None,
    rp_guess_rearth: float | None = None,
) -> StructureResult:
    """Build one converged static structure model.

    ``band`` fixes the opacity band; None -> the source's own band (A21).
    Comparison figures pass the common band of all plotted sources so every
    overlaid profile shares one exact band. ``rp_guess_rearth`` warm-starts the
    R_p shooting (cheap for sequences of nearby models).
    """
    source_name = source_name or cfg.opacity.source
    eos = eos or AquaEOS.from_config(cfg)
    opacity = get_source(source_name, cfg, band=band)

    mc_mearth = cfg.core.mass_mearth
    f_env = cfg.envelope.f_env
    mp_mearth = mc_mearth / (1.0 - f_env)
    mass_p = core_mass_si(mp_mearth)
    mass_core = core_mass_si(mc_mearth)
    radius_core = core_radius_si(mc_mearth, cmf=cfg.core.cmf)

    gamma = cfg.irradiation.gamma
    p_top = cfg.boundary.p_top_bar * 1e5
    p_match = cfg.boundary.match_pressure_bar * 1e5
    num = cfg.numerics

    shot = shoot_radius(
        eos, opacity, mass_p, mass_core, radius_core, t_int_k, t_eq_k, gamma,
        p_top, p_match, num.rtol, num.atol, num.solver,
        (num.rp_bracket_rearth[0] * R_EARTH, num.rp_bracket_rearth[1] * R_EARTH),
        num.shoot_rtol, num.max_shoot_iter,
        rp_guess_m=None if rp_guess_rearth is None else rp_guess_rearth * R_EARTH,
    )

    # Final sampled profile at the converged radius; measure clamp fraction here.
    opacity.reset_clamp_counter()
    match, outer, matched_at_rcb = integrate_outer(
        eos, opacity, shot.radius_p_m, mass_p, t_int_k, t_eq_k, gamma,
        p_top, p_match, num.rtol, num.atol, num.solver, n_sample=n_sample,
    )
    lum = luminosity_int(shot.radius_p_m, t_int_k)
    _, interior = integrate_interior(
        eos, opacity, match, mass_core, radius_core, lum,
        num.rtol, num.atol, num.solver, n_sample=n_sample,
    )
    clamp_fraction = opacity.clamp_fraction
    flat_rho_fraction = opacity.flat_rho_fraction

    # Concatenate top -> core (drop duplicate matching point).
    m = np.concatenate([outer["m"], interior["m"][1:]])
    r = np.concatenate([outer["r"], interior["r"][1:]])
    p = np.concatenate([outer["P"], interior["P"][1:]])
    t = np.concatenate([outer["T"], interior["T"][1:]])
    tau_outer = outer["tau"]

    rho, grad_ad, grad_rad, is_conv, two_phase, u, entropy = _post_process(
        m, r, p, t, eos, opacity, lum
    )
    grad = np.gradient(np.log(t), np.log(p))
    energy_j = _envelope_energy(m, r, u)
    tau_kh_s = abs(energy_j) / lum if lum > 0 else float("inf")
    # Interior (deep convective adiabat) specific entropy: the deepest point.
    interior_entropy = float(entropy[-1])

    rcbs = _find_rcbs(p, is_conv)
    rcb_rho = _rcb_densities(rho, is_conv)
    radiative_rho = rho[~is_conv]
    max_rad_rho = float(radiative_rho.max()) if radiative_rho.size else 0.0
    p_phot, t_phot = _photosphere(outer["P"], tau_outer, outer["T"])
    hydro = _hydrostatic_residual(interior["m"], interior["r"], interior["P"])
    threshold = cfg.opacity.clamp_fraction_threshold

    # Reliability is set by clamps in the RADIATIVE zone only: opacity shapes the
    # structure just there; convective clamps (deep interior above the grid T)
    # are physically irrelevant since the profile follows the adiabat (A26).
    radiative = ~is_conv
    if np.any(radiative):
        rad_clamped = opacity.clamped_mask(t[radiative], rho[radiative])
        clamp_fraction_radiative = float(np.mean(rad_clamped))
    else:
        clamp_fraction_radiative = 0.0

    # Composite: log the FP/line-list κ ratio at the seam over the radiative
    # zone (the quantified systematic, A12).
    seam_ratio = None
    if isinstance(opacity, CompositeOpacitySource):
        rad_t = t[~is_conv] if np.any(~is_conv) else t
        temps = np.linspace(float(rad_t.min()), float(rad_t.max()), 5)
        seam_ratio = opacity.seam_ratio(temps)

    result = StructureResult(
        name=cfg.name,
        source=source_name,
        t_eq_k=t_eq_k,
        t_int_k=t_int_k,
        gamma=gamma,
        mass_core_mearth=mc_mearth,
        f_env=f_env,
        band_cm=getattr(opacity, "band_cm", None),
        match_pressure_bar=cfg.boundary.match_pressure_bar,
        match_pressure_actual_bar=match.pressure_pa / 1e5,
        matched_at_rcb=matched_at_rcb,
        radius_p_m=shot.radius_p_m,
        radius_p_rearth=shot.radius_p_m / R_EARTH,
        radius_core_m=radius_core,
        mass_p_kg=mass_p,
        mass_env_kg=mass_p - mass_core,
        luminosity_w=lum,
        envelope_energy_j=energy_j,
        tau_kh_s=tau_kh_s,
        interior_entropy_j_kg_k=interior_entropy,
        rcb_pressures_pa=rcbs,
        rcb_densities_si=rcb_rho,
        max_radiative_rho_si=max_rad_rho,
        reaches_fp_grid=max_rad_rho >= FP_GRID_FLOOR_SI,
        rcb_above_match=any(rp < p_match for rp in rcbs),
        photosphere_pressure_pa=p_phot,
        photosphere_temperature_k=t_phot,
        clamp_fraction=clamp_fraction,
        clamp_fraction_radiative=clamp_fraction_radiative,
        clamp_reliable=clamp_fraction_radiative <= threshold,
        flat_rho_fraction=flat_rho_fraction,
        two_phase_fraction=float(two_phase.mean()),
        hydrostatic_residual=hydro,
        shoot_iterations=shot.n_iter,
        shoot_residual_m=shot.residual_m,
        composite_seam_ratio=seam_ratio,
        profile={
            "mass_kg": m, "radius_m": r, "pressure_pa": p, "temperature_k": t,
            "rho_si": rho, "grad": grad, "grad_ad": grad_ad, "grad_rad": grad_rad,
            "is_convective": is_conv, "two_phase": two_phase, "entropy_j_kg_k": entropy,
        },
    )
    return result


def build_at_entropy(
    cfg: Config,
    s_target_j_kg_k: float,
    t_eq_k: float,
    source_name: str | None = None,
    eos: AquaEOS | None = None,
    t_int_bracket_k: tuple[float, float] = (20.0, 400.0),
    tol_frac: float = 1e-3,
    band: tuple[float, float] | None = None,
    n_sample: int = 800,
) -> StructureResult:
    """Converged model whose interior entropy equals ``s_target_j_kg_k`` (A23).

    Interior entropy increases monotonically with T_int, so we shoot on T_int by
    bisection to match the target adiabat. This is the entropy-parameterized
    sequence used for the fixed-entropy F3 headline and Phase-2 evolution: at a
    fixed interior entropy, a more opaque envelope throttles the flux, giving a
    lower T_int (hence lower L) and a longer cooling time.
    """
    from scipy.optimize import brentq

    eos = eos or AquaEOS.from_config(cfg)

    def entropy_gap(t_int_k: float) -> float:
        res = build_structure(cfg, t_eq_k, float(t_int_k), source_name, eos=eos,
                              n_sample=n_sample, band=band)
        entropy_gap.last = res  # type: ignore[attr-defined]
        return res.interior_entropy_j_kg_k - s_target_j_kg_k

    lo, hi = t_int_bracket_k
    f_lo, f_hi = entropy_gap(lo), entropy_gap(hi)
    if f_lo * f_hi > 0:
        raise ValueError(
            f"interior entropy {s_target_j_kg_k:.1f} J/kg/K not reachable for "
            f"source={source_name}, T_eq={t_eq_k} over T_int∈{t_int_bracket_k}: "
            f"reachable [{f_lo + s_target_j_kg_k:.1f}, {f_hi + s_target_j_kg_k:.1f}]"
        )
    t_int_star = brentq(entropy_gap, lo, hi, rtol=tol_frac, maxiter=60)
    return entropy_gap.last  # type: ignore[attr-defined]


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
