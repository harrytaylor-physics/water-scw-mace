"""Build gridded opacity tables from per-point spectra (SPEC §2, §3).

For a registered source: load its spectra, compute band-limited Planck and
Rosseland means over the common band (cm²/g), assert the factor-1000
plausibility check (SPEC §2), convert to SI, and wrap in a
:class:`GridOpacitySource`. Sources are cached by (name, band, data_root).
"""

from __future__ import annotations

import numpy as np

from config import Config
from opacity.loaders import find_spectra, load_spectrum
from opacity.means import planck_mean, rosseland_mean
from opacity.registry import SourceSpec, get_spec
from opacity.source import GridOpacitySource
from opacity.table import CM2_PER_G_TO_SI, OpacityTable

# SPEC §2: band-integrated kappa_P must land in ~1e1-1e4 cm²/g. Hard-fail bounds
# are one decade wider each side (avoids edge false-positives while still
# catching a factor-1000 error, which lands kappa_P at ~1e5-1e6 cm²/g).
_PLAUSIBLE_LO_CM2_G = 1.0e0
_PLAUSIBLE_HI_CM2_G = 1.0e5

_cache: dict[tuple[str, tuple[float, float], str], GridOpacitySource] = {}


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


def build_source(name: str, cfg: Config) -> GridOpacitySource:
    """Build a fresh gridded opacity source for ``name`` (no caching)."""
    band_cm = tuple(cfg.opacity.band_cm)  # type: ignore[assignment]
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


def get_source(name: str, cfg: Config) -> GridOpacitySource:
    """Build or return a cached opacity source for ``name``.

    Cached sources share a clamp counter; reset it before a structure
    integration via ``source.reset_clamp_counter()``.
    """
    key = (name, tuple(cfg.opacity.band_cm), str(cfg.data_root_path))
    if key not in _cache:
        _cache[key] = build_source(name, cfg)
    return _cache[key]


def clear_cache() -> None:
    _cache.clear()
