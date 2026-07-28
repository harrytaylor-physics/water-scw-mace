"""Per-point opacity spectrum loading (SPEC §2, §3).

Loads a source's ``T*_opacity.npz`` spectra, applies the source's on-load scale
factor (FP: ×1e-3 for the legacy factor-1000 error), and verifies each file
covers the common band. Returns lightweight records the table builder consumes.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from opacity.registry import SourceSpec


@dataclass
class SpectrumRecord:
    """One (T, ρ) grid point's monochromatic opacity, in cm²/g."""

    temperature_k: float
    rho_si: float  # kg m⁻³
    wn_cm: np.ndarray  # wavenumber, cm⁻¹
    kappa_cm2_g: np.ndarray  # opacity, cm²/g (already scaled)
    path: Path


def load_spectrum(path: Path, spec: SourceSpec, band_cm: tuple[float, float]) -> SpectrumRecord:
    """Load and scale one spectrum npz; assert it spans ``band_cm``."""
    with np.load(path) as d:
        wn = np.asarray(d["wn"], dtype=float)
        kappa = np.asarray(d["kappa"], dtype=float) * spec.scale
        temperature_k = float(d["T"])
        rho_si = float(d["rho_kg_m3"])

    lo, hi = band_cm
    if wn.min() > lo or wn.max() < hi:
        raise ValueError(
            f"{path.name}: wn span [{wn.min():.3g}, {wn.max():.3g}] cm^-1 does not "
            f"cover common band [{lo:g}, {hi:g}] cm^-1"
        )
    if not np.all(np.isfinite(kappa)) or np.any(kappa <= 0):
        raise ValueError(f"{path.name}: opacity must be finite and positive")

    order = np.argsort(wn)  # trapz needs monotone abscissa
    return SpectrumRecord(
        temperature_k=temperature_k,
        rho_si=rho_si,
        wn_cm=wn[order],
        kappa_cm2_g=kappa[order],
        path=path,
    )


def find_spectra(spectra_dir: Path, spec: SourceSpec) -> list[Path]:
    """Return sorted spectra paths, or raise a clear 'not delivered' error."""
    paths = sorted(spectra_dir.glob(spec.glob))
    if not paths:
        raise FileNotFoundError(
            f"opacity source {spec.name!r}: no spectra matching {spec.glob!r} in "
            f"{spectra_dir}. If this is a line-list source, its per-point spectra "
            f"are not yet delivered (see AMENDMENTS A5)."
        )
    return paths
