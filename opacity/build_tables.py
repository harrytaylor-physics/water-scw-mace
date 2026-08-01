"""Build gridded opacity tables from per-point spectra (SPEC §2, §3).

For a registered source: load its spectra, compute band-limited Planck and
Rosseland means over the common band (cm²/g), assert the factor-1000
plausibility check (SPEC §2), convert to SI, and wrap in a
:class:`GridOpacitySource`. Sources are cached by (name, band, data_root).
"""

from __future__ import annotations

import numpy as np

from config import Config
from opacity.composite import CompositeOpacitySource
from opacity.loaders import find_spectra, load_spectrum
from opacity.means import planck_mean, rosseland_mean
from opacity.registry import REGISTRY, SourceSpec, get_spec
from opacity.source import GridOpacitySource, OpacitySource
from opacity.table import CM2_PER_G_TO_SI, OpacityTable

# SPEC §2: band-integrated kappa_P must land in ~1e1-1e4 cm²/g. Hard-fail bounds
# are one decade wider each side (avoids edge false-positives while still
# catching a factor-1000 error, which lands kappa_P at ~1e5-1e6 cm²/g).
_PLAUSIBLE_LO_CM2_G = 1.0e0
_PLAUSIBLE_HI_CM2_G = 1.0e5

_cache: dict[tuple[str, tuple[float, float], str], OpacitySource] = {}
_band_cache: dict[str, tuple[float, float]] = {}


def detect_common_band(cfg: Config, sources: list[str] | None = None) -> tuple[float, float]:
    """Largest wavenumber band the given sources share (A14, A21).

    = [max over spectra files of min(wn), min over files of max(wn)], verified
    per file. ``sources`` limits the intersection to those grid sources (band-
    agnostic sources like ``constant`` and undelivered sources are ignored);
    ``None`` uses every registered source with spectra on disk. Per-comparison
    banding (A21): a figure passes exactly the grid sources it plots, so a
    windowed dataset (e.g. POKAZATEL [50,4500]) only constrains comparisons it
    is part of.
    """
    names = sources if sources is not None else list(REGISTRY)
    key = (str(cfg.data_root_path), tuple(sorted(names)))
    if key in _band_cache:
        return _band_cache[key]
    mins: list[float] = []
    maxs: list[float] = []
    for name in names:
        if name not in REGISTRY:  # e.g. "constant" — band-agnostic, skip
            continue
        spec = REGISTRY[name]
        for path in sorted((cfg.data_root_path / spec.subdir).glob(spec.glob)):
            with np.load(path) as d:
                wn = np.asarray(d["wn"], dtype=float)
            mins.append(float(wn.min()))
            maxs.append(float(wn.max()))
    if not mins:
        raise FileNotFoundError(
            f"no opacity spectra found to detect a common band for sources={names}"
        )
    band = (max(mins), min(maxs))
    _band_cache[key] = band
    return band


def resolve_band(cfg: Config, sources: list[str] | None = None) -> tuple[float, float]:
    """Config band if pinned, else the auto-detected band for ``sources`` (A21).

    ``sources`` defaults to all registered sources with spectra. For a single
    source this yields its own full band; for a comparison, pass the plotted
    sources so every curve shares one exact band.
    """
    if cfg.opacity.band_cm is not None:
        return tuple(cfg.opacity.band_cm)  # type: ignore[return-value]
    return detect_common_band(cfg, sources)


def compute_point_means(
    name: str, cfg: Config, band_cm: tuple[float, float]
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Per-point (T, ρ, κ_R, κ_P) means in cm²/g over ``band_cm``.

    Reusable by the table builder and by the band-limitation quantification
    (full-band vs common-band, docs/band_limitation.md).
    """
    spec = get_spec(name)
    spectra_dir = cfg.data_root_path / spec.subdir
    paths = find_spectra(spectra_dir, spec)

    temps, rhos, kappa_r, kappa_p = [], [], [], []
    for path in paths:
        rec = load_spectrum(path, spec, band_cm)
        kappa_p.append(planck_mean(rec.wn_cm, rec.kappa_cm2_g, rec.temperature_k, band_cm))
        kappa_r.append(rosseland_mean(rec.wn_cm, rec.kappa_cm2_g, rec.temperature_k, band_cm))
        temps.append(rec.temperature_k)
        rhos.append(rec.rho_si)
    return (
        np.array(temps),
        np.array(rhos),
        np.array(kappa_r),
        np.array(kappa_p),
    )


def _assert_plausible(name: str, spec: SourceSpec, kappa_p_cm2_g: np.ndarray) -> None:
    """Hard-fail if band-integrated Planck means are implausible (SPEC §2)."""
    lo, hi = float(kappa_p_cm2_g.min()), float(kappa_p_cm2_g.max())
    if lo < _PLAUSIBLE_LO_CM2_G or hi > _PLAUSIBLE_HI_CM2_G:
        raise AssertionError(
            f"opacity source {name!r}: band-integrated kappa_P in "
            f"[{lo:.3g}, {hi:.3g}] cm^2/g is outside the plausible window "
            f"[{_PLAUSIBLE_LO_CM2_G:g}, {_PLAUSIBLE_HI_CM2_G:g}] cm^2/g. "
            f"Applied scale={spec.scale:g} ({spec.scale_reason}). "
            f"A ~1e3 offset indicates the factor-1000 scaling is wrong."
        )


def build_composite(cfg: Config, band: tuple[float, float] | None = None) -> CompositeOpacitySource:
    """Build the composite source: line-list (dilute) blended with FP (dense).

    The two sub-sources share the band common to *them* (not to POKAZATEL or any
    other on-disk dataset), so the composite is not capped by an unrelated
    windowed source (A21).
    """
    op = cfg.opacity
    names = [op.composite_dilute_source, op.composite_dense_source]
    band = band or resolve_band(cfg, names)
    dilute = build_source(op.composite_dilute_source, cfg, band=band)
    dense = build_source(op.composite_dense_source, cfg, band=band)
    return CompositeOpacitySource(
        dilute=dilute, dense=dense, blend_rho_si=tuple(op.composite_blend_rho)
    )


def _default_band(name: str, cfg: Config) -> tuple[float, float]:
    """A source's own band; raise the delivery-specific error if it has no spectra."""
    names = (
        [cfg.opacity.composite_dilute_source, cfg.opacity.composite_dense_source]
        if name == "composite"
        else [name]
    )
    try:
        return resolve_band(cfg, names)
    except FileNotFoundError:
        # Surface the clearer "not yet delivered" message per missing source.
        for n in names:
            find_spectra(cfg.data_root_path / get_spec(n).subdir, get_spec(n))
        raise  # pragma: no cover - find_spectra already raised


def build_source(name: str, cfg: Config, band: tuple[float, float] | None = None) -> OpacitySource:
    """Build a fresh opacity source for ``name`` (no caching).

    ``band`` fixes the Rosseland/Planck band; when ``None`` the source uses its
    *own* band (A21), so a pure-source run is independent of what else is on
    disk. Comparison contexts pass the shared band of the plotted sources.
    """
    if name == "composite":
        return build_composite(cfg, band=band)
    band_cm = band or _default_band(name, cfg)
    temps, rhos, kappa_r, kappa_p = compute_point_means(name, cfg, band_cm)
    _assert_plausible(name, get_spec(name), kappa_p)

    table = OpacityTable(
        temperature_k=temps,
        rho_si=rhos,
        kappa_rosseland_si=kappa_r * CM2_PER_G_TO_SI,
        kappa_planck_si=kappa_p * CM2_PER_G_TO_SI,
        band_cm=band_cm,
    )
    return GridOpacitySource(name, table)


def get_source(name: str, cfg: Config, band: tuple[float, float] | None = None) -> OpacitySource:
    """Build or return a cached opacity source for ``name`` on ``band``.

    ``band`` None -> the source's own band (A21). Cached sources share a clamp
    counter; reset it before a structure integration via
    ``source.reset_clamp_counter()``.
    """
    root = str(cfg.data_root_path)
    # Band-agnostic injected sources (e.g. a constant-κ reference) are keyed with
    # band=None; they are placed in the cache directly and only looked up here.
    if name != "composite" and name not in REGISTRY:
        key = (name, None, root)
        if key in _cache:
            return _cache[key]
        raise KeyError(f"opacity source {name!r} is not registered or cached")
    band_cm = band or _default_band(name, cfg)
    key = (name, band_cm, root)
    if key not in _cache:
        _cache[key] = build_source(name, cfg, band=band_cm)
    return _cache[key]


def clear_cache() -> None:
    _cache.clear()
    _band_cache.clear()
