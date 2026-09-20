"""Typed run configuration (SPEC §1, §5, §8).

A run config is a YAML file; :func:`load_config` deep-merges it over the
dataclass defaults, resolves data paths, validates, and returns a frozen
``Config``. All tolerances and physics knobs live here, never hard-coded in
modules (SPEC §5, §10).
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field, fields, is_dataclass
from pathlib import Path
from typing import Any, get_type_hints

import yaml


# --------------------------------------------------------------------------
# Sub-configs
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class Paths:
    """Filesystem locations. Data is never committed (SPEC §1)."""

    data_root: str = "~/subneptune_data"
    # Relative entries are resolved against ``data_root`` (absolute kept as-is).
    # Opacity source subdirectories are governed by ``opacity.registry`` (adding
    # a source is a registry entry, SPEC §3), not by config.
    aqua_table: str = "AQUA/Tables/aqua_eos_rhot_v1_0.dat"
    output_dir: str = "outputs"


@dataclass(frozen=True)
class Core:
    """Inert rocky core (SPEC §4)."""

    mass_mearth: float = 5.0  # M_c in Earth masses
    mr_relation: str = "zeng2016_rock"  # Earth-like rock mass-radius relation
    cmf: float = 0.325  # core (iron) mass fraction; Earth-like default (A9). 0 =
    #                     pure rock (matches Aguichine's no-metallic-core, G5).


@dataclass(frozen=True)
class Envelope:
    """Pure-steam envelope (SPEC §0, §5)."""

    f_env: float = 0.05  # M_env / M_p


@dataclass(frozen=True)
class Irradiation:
    """Outer irradiation boundary (Guillot 2010; SPEC §5)."""

    t_eq_k: tuple[float, ...] = (500.0, 600.0, 700.0, 800.0)  # grid
    # gamma = kappa_visible / kappa_thermal. Fiducial + sensitivity is a knob
    # Ray will ask about; sensitivity scan lives in the structure driver.
    gamma: float = 1.0


@dataclass(frozen=True)
class Interior:
    """Interior thermal state grid (SPEC §5, §8)."""

    t_int_k: tuple[float, ...] = (30.0, 60.0, 100.0, 150.0, 200.0, 300.0)


@dataclass(frozen=True)
class Opacity:
    """Opacity source + table/interpolation policy (SPEC §3)."""

    source: str = "fp"  # registered name: fp | hitemp | pokazatel | composite
    # Common Rosseland/Planck band. Default [100, 4000] cm⁻¹ is the
    # crossover-consistent science band (A25): it excludes the far-band
    # transparency windows (50-100, 4000-4500 cm⁻¹) that dominate the Rosseland
    # harmonic mean and drive a ~2× band systematic, and it reproduces the
    # independent crossover analysis to ~1%. None -> auto-detect the largest
    # band all sources share (A14/A21); [50,4500] is the reported systematic.
    band_cm: tuple[float, float] | None = (100.0, 4000.0)
    # Relative opacity floor before the harmonic Rosseland mean: κ ≥ rel_floor·
    # max(κ) (A25). Guards against zero-transparency windows; inactive for the
    # current data at 1e-6 (window minima are shallower), kept for convention
    # parity and robustness.
    rel_floor: float = 1.0e-6
    # Extrapolation is a gate, not a guess (SPEC §3): clamp to nearest grid
    # value, count clamped calls; a run whose clamped fraction exceeds this
    # threshold is flagged unreliable in its output metadata.
    clamp_fraction_threshold: float = 0.20

    # Composite source (A11/A12): line-list (dilute) blended with FP (dense)
    # across a density window (kg/m³). The window's lower edge must not fall
    # below the FP grid's ρ floor, or the composite hands κ to FP where FP has no
    # data. Seam location/width remain scan parameters (A35: the surviving RCB
    # sits inside the window, so the blend prescription is a stated systematic,
    # scanned by figures/plot_s1_blend_sensitivity.py).
    composite_dilute_source: str = "hitemp"
    composite_dense_source: str = "fp"
    # A54: moved 100-300 -> 10-30 once the extension campaign put FP data down to
    # ρ = 10. Chosen on κ_R seam continuity (see AMENDMENTS): the seam FP/HITEMP
    # ratio falls from 10.5× (100-300) to 3.7× (10-30). The 3-10 window is
    # narrowly better on continuity (3.1×) but its seam at ρ = 5.5 is BELOW the FP
    # grid's floor, so FP is clamped there by construction — disqualified.
    blend_window_kg_m3: tuple[float, float] = (10.0, 30.0)
    # "linear": the existing log-space blend across the window (smoothstep
    # weight, C¹ — unchanged default behaviour, A12/A35). "sharp": a hard switch
    # to FP at the window's lower edge, no blending (the window may be
    # zero-width, e.g. [108.796, 108.796] = switch at the FP grid floor).
    blend_mode: str = "linear"
    # The FP grid's own (T, ρ) floor — the corner below which it has no data, so
    # "in the FP box" and "FP-relevant" are defined against it (A41/A42/A47). Read
    # from config, not hard-coded, because the extension campaign moved it from
    # (1000 K, 108.796) to (650 K, 10); a gate asserts these equal the built FP
    # table's actual minima, so config cannot drift from the data (A54).
    fp_grid_t_floor_k: float = 650.0
    fp_grid_rho_floor_kg_m3: float = 10.0
    # Constant multiplier on the FP table's κ_R and κ_P, applied at table-build
    # time (before any blending) — the opacity-perturbation robustness knob
    # (A37). 1.0 = unperturbed; the manuscript scan is 0.8 / 1.0 / 1.25.
    fp_scale: float = 1.0

    @property
    def composite_blend_rho(self) -> tuple[float, float]:
        """Deprecated alias for :attr:`blend_window_kg_m3` (A35)."""
        return self.blend_window_kg_m3


@dataclass(frozen=True)
class Boundary:
    """Interior<->atmosphere matching (SPEC §5)."""

    match_pressure_bar: float = 10.0  # default matching pressure
    # Sensitivity scan range (1-100 bar), recorded in boundary_sensitivity.md.
    match_pressure_scan_bar: tuple[float, ...] = (1.0, 10.0, 100.0)
    p_top_bar: float = 1.0e-5  # top-of-atmosphere pressure (tau ~ 0)


@dataclass(frozen=True)
class Anchor:
    """Conventions for matching the published literature anchor, gate G5 (A43).

    These are *definition* choices, not physics knobs: how the anchor paper's
    reported quantities map onto ours. Since A46 the defaults are the *corrected*
    mappings established in A43; the pre-A46 conventions remain selectable so the
    audit can still attribute each change one at a time.
    """

    # How Aguichine et al. (2021) T_irr maps to our T_eq. "identity" (the
    # default since A46) passes it through, which is what their Equation (9)
    # actually defines, T_irr = T_eff·√(R⋆/2a) — verified against their own
    # Table 3 (A43). "sqrt2" divides by √2, treating their T_irr as Guillot's
    # irradiation temperature: the original A32 assumption, kept for the audit.
    t_irr_to_t_eq: str = "identity"
    # Pressure level the anchor's reported (transiting) radius refers to. Their
    # atmospheric grid reports the thickness at 0.1 Pa = 1e-6 bar; ours would
    # otherwise stop at boundary.p_top_bar (1e-5 bar), one decade deeper.
    # None -> use boundary.p_top_bar unchanged (the pre-A46 default).
    radius_pressure_bar: float | None = 1.0e-6
    # Their x_core = CMF/(1−WMF) is the core fraction of the *refractory* part;
    # Figure 8 is computed with no metallic core, i.e. x_core = 0 -> cmf = 0.
    cmf: float = 0.0


@dataclass(frozen=True)
class Numerics:
    """Integrator + shooting tolerances (SPEC §5). All in config, not code."""

    solver: str = "LSODA"  # scipy solve_ivp method: LSODA | Radau
    rtol: float = 1e-8
    atol: float = 1e-10
    shoot_rtol: float = 1e-6  # convergence on r(M_c) = R_c(M_c)
    max_shoot_iter: int = 100
    # R_p search bracket, in Earth radii, for the shooting solver.
    rp_bracket_rearth: tuple[float, float] = (1.0, 15.0)


@dataclass(frozen=True)
class Config:
    """Top-level run configuration."""

    name: str = "unnamed"
    paths: Paths = field(default_factory=Paths)
    core: Core = field(default_factory=Core)
    envelope: Envelope = field(default_factory=Envelope)
    irradiation: Irradiation = field(default_factory=Irradiation)
    interior: Interior = field(default_factory=Interior)
    opacity: Opacity = field(default_factory=Opacity)
    boundary: Boundary = field(default_factory=Boundary)
    anchor: Anchor = field(default_factory=Anchor)
    numerics: Numerics = field(default_factory=Numerics)

    # -- resolved absolute paths (convenience) --
    @property
    def data_root_path(self) -> Path:
        return Path(self.paths.data_root).expanduser()

    @property
    def aqua_table_path(self) -> Path:
        return _resolve_under(self.data_root_path, self.paths.aqua_table)

    @property
    def fp_opacity_path(self) -> Path:
        """Where FP spectra + grid.json + crossover_hitemp.json live."""
        return self.data_root_path / "fp_opacity"

    @property
    def output_path(self) -> Path:
        return Path(self.paths.output_dir).expanduser()


# --------------------------------------------------------------------------
# Loading
# --------------------------------------------------------------------------
def _resolve_under(root: Path, entry: str) -> Path:
    """Resolve ``entry`` under ``root`` unless it is already absolute."""
    p = Path(entry).expanduser()
    return p if p.is_absolute() else root / p


def _build(cls: type, data: Any) -> Any:
    """Recursively construct a (possibly nested) dataclass from a dict,
    coercing list values to tuples for tuple-typed fields.

    ``from __future__ import annotations`` makes field types strings, so
    nested-dataclass detection goes through resolved ``get_type_hints``.
    """
    if not is_dataclass(cls):
        return data
    if data is None:
        return cls()
    if not isinstance(data, dict):
        raise TypeError(f"expected mapping for {cls.__name__}, got {type(data).__name__}")

    known = {f.name for f in fields(cls)}
    unknown = set(data) - known
    if unknown:
        raise ValueError(f"unknown config keys for {cls.__name__}: {sorted(unknown)}")
    hints = get_type_hints(cls)

    kwargs: dict[str, Any] = {}
    for key, raw in data.items():
        hint = hints.get(key)
        if isinstance(hint, type) and is_dataclass(hint):
            kwargs[key] = _build(hint, raw)  # nested dataclass
        elif isinstance(raw, list):
            kwargs[key] = tuple(raw)  # YAML lists -> tuples
        else:
            kwargs[key] = raw
    return cls(**kwargs)


def load_config(path: str | Path) -> Config:
    """Load and validate a run config from a YAML file."""
    path = Path(path).expanduser()
    with path.open("r") as fh:
        raw = yaml.safe_load(fh) or {}
    if not isinstance(raw, dict):
        raise TypeError(f"config root must be a mapping, got {type(raw).__name__}")
    cfg = _build(Config, raw)
    _validate(cfg)
    return cfg


def _validate(cfg: Config) -> None:
    """Cheap invariants; physics validation lives in the gates."""
    if not 0.0 < cfg.envelope.f_env < 1.0:
        raise ValueError(f"envelope.f_env must be in (0,1), got {cfg.envelope.f_env}")
    if cfg.core.mass_mearth <= 0:
        raise ValueError("core.mass_mearth must be positive")
    if cfg.opacity.band_cm is not None:
        lo, hi = cfg.opacity.band_cm
        if not lo < hi:
            raise ValueError(f"opacity.band_cm must be increasing, got {cfg.opacity.band_cm}")
    if not 0.0 <= cfg.opacity.rel_floor < 1.0:
        raise ValueError(f"opacity.rel_floor must be in [0,1), got {cfg.opacity.rel_floor}")
    if not 0.0 < cfg.opacity.clamp_fraction_threshold <= 1.0:
        raise ValueError("opacity.clamp_fraction_threshold must be in (0,1]")
    if cfg.opacity.source not in {"fp", "hitemp", "pokazatel", "composite"}:
        raise ValueError(f"unknown opacity.source {cfg.opacity.source!r}")
    if cfg.opacity.blend_mode not in {"linear", "sharp"}:
        raise ValueError(f"opacity.blend_mode must be linear|sharp, got {cfg.opacity.blend_mode!r}")
    lo, hi = cfg.opacity.blend_window_kg_m3
    if not 0.0 < lo <= hi:
        raise ValueError(f"opacity.blend_window_kg_m3 must satisfy 0 < lo <= hi, got {(lo, hi)}")
    if cfg.opacity.blend_mode == "linear" and lo == hi:
        raise ValueError("a zero-width blend window requires opacity.blend_mode='sharp'")
    if cfg.opacity.fp_grid_t_floor_k <= 0 or cfg.opacity.fp_grid_rho_floor_kg_m3 <= 0:
        raise ValueError("opacity.fp_grid_* floors must be positive")
    if cfg.opacity.fp_scale <= 0:
        raise ValueError(f"opacity.fp_scale must be positive, got {cfg.opacity.fp_scale}")
    if cfg.anchor.t_irr_to_t_eq not in {"sqrt2", "identity"}:
        raise ValueError(
            f"anchor.t_irr_to_t_eq must be sqrt2|identity, got {cfg.anchor.t_irr_to_t_eq!r}"
        )
    if cfg.anchor.radius_pressure_bar is not None and cfg.anchor.radius_pressure_bar <= 0:
        raise ValueError("anchor.radius_pressure_bar must be positive or null")
    if not 0.0 <= cfg.anchor.cmf < 1.0:
        raise ValueError(f"anchor.cmf must be in [0,1), got {cfg.anchor.cmf}")
    if cfg.numerics.solver not in {"LSODA", "Radau", "BDF", "RK45"}:
        raise ValueError(f"unsupported numerics.solver {cfg.numerics.solver!r}")
    if not all(t > 0 for t in cfg.irradiation.t_eq_k):
        raise ValueError("irradiation.t_eq_k must be positive")
    if not all(t > 0 for t in cfg.interior.t_int_k):
        raise ValueError("interior.t_int_k must be positive")


def to_dict(cfg: Config) -> dict[str, Any]:
    """Serialize a Config back to a plain dict (for run metadata/JSON)."""
    return asdict(cfg)
