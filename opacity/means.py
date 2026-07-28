"""Band-limited Rosseland and Planck opacity means (SPEC §3).

Means are computed by integrating monochromatic opacity against the Planck
weight (Planck mean) or the temperature-derivative-of-Planck weight
(Rosseland mean, harmonic) over a wavenumber band. The band is the common band
shared by all sources so the *differential* conclusions between sources are
exact (SPEC §3). Constants come from ``constants.py``.

All spectral integrals are over wavenumber ``wn`` in cm⁻¹. The Planck function
per unit wavenumber is B_wn ∝ wn³ / (exp(x) − 1) with
x = h c wn / (k_B T); constant prefactors cancel in every weighted mean here,
so they are omitted. Opacity units are passed through unchanged (the caller
works in cm²/g; SI conversion happens at table-build time).
"""

from __future__ import annotations

import numpy as np

from constants import C_LIGHT, H_PLANCK, K_B

# Speed of light in cm/s, so that c·wn (wn in cm⁻¹) is a frequency in s⁻¹.
_C_CM_S = C_LIGHT * 100.0


def _x(wn_cm: np.ndarray, temperature_k: float) -> np.ndarray:
    """Dimensionless photon energy x = h c wn / (k_B T)."""
    return H_PLANCK * _C_CM_S * np.asarray(wn_cm, dtype=float) / (K_B * temperature_k)


def planck_weight(wn_cm: np.ndarray, temperature_k: float) -> np.ndarray:
    """Planck weight ∝ B_wn(T) (per unit wavenumber), unnormalized."""
    x = _x(wn_cm, temperature_k)
    return wn_cm**3 / np.expm1(x)


def rosseland_weight(wn_cm: np.ndarray, temperature_k: float) -> np.ndarray:
    """Rosseland weight ∝ dB_wn/dT (per unit wavenumber), unnormalized."""
    x = _x(wn_cm, temperature_k)
    ex = np.exp(x)
    return wn_cm**4 * ex / np.expm1(x) ** 2


def _band_slice(wn_cm: np.ndarray, band_cm: tuple[float, float]) -> np.ndarray:
    lo, hi = band_cm
    mask = (wn_cm >= lo) & (wn_cm <= hi)
    if mask.sum() < 2:
        raise ValueError(
            f"band {band_cm} cm^-1 contains < 2 spectral points "
            f"(wn spans [{wn_cm.min():.3g}, {wn_cm.max():.3g}])"
        )
    return mask


def planck_mean(
    wn_cm: np.ndarray, kappa: np.ndarray, temperature_k: float, band_cm: tuple[float, float]
) -> float:
    """Planck mean opacity ⟨κ⟩_P = ∫κ B dν / ∫B dν over ``band_cm``.

    Units of the result match the units of ``kappa``.
    """
    mask = _band_slice(wn_cm, band_cm)
    wn = wn_cm[mask]
    kap = kappa[mask]
    w = planck_weight(wn, temperature_k)
    return float(np.trapz(kap * w, wn) / np.trapz(w, wn))


def rosseland_mean(
    wn_cm: np.ndarray, kappa: np.ndarray, temperature_k: float, band_cm: tuple[float, float]
) -> float:
    """Rosseland mean opacity via the harmonic weighting
    1/⟨κ⟩_R = ∫(1/κ)(dB/dT) dν / ∫(dB/dT) dν over ``band_cm``.

    Units of the result match the units of ``kappa``. Requires κ > 0 in-band.
    """
    mask = _band_slice(wn_cm, band_cm)
    wn = wn_cm[mask]
    kap = kappa[mask]
    if not np.all(kap > 0):
        raise ValueError("Rosseland mean requires strictly positive in-band opacity")
    w = rosseland_weight(wn, temperature_k)
    num = np.trapz(w, wn)
    den = np.trapz(w / kap, wn)
    return float(num / den)
