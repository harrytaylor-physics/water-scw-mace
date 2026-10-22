"""Per-band absorption integrals, for attributing an FP/line-list ratio (A61).

The dilute cells (ρ = 10-60 kg/m³) show an FP Planck mean 1.32-1.53× bare HITEMP
while their Rosseland ratio falls to 1.08. A band-integrated mean cannot say
*where* that intensity excess sits, and the answer decides which part of the
first-principles calculation it is attributable to. This module splits the band.

The quantity integrated is

    A_band = ∫_band ν (1 − e^{−hcν/k_BT}) κ(ν) dν

i.e. the opacity weighted by the same ν(1 − e^{−hcν/k_BT}) factor that relates a
dipole spectral density to an absorption cross-section. Comparing A_band between
two sources over identical limits is therefore a comparison of integrated band
intensity, not of a radiative mean: no Planck or Rosseland weighting enters, so a
ratio here is not diluted by where the Planck function happens to peak.

Band windows (cm⁻¹) are the standard water assignments:

* ``libration`` 300-1000 — hindered rotation;
* ``bend`` 1400-1800 — the ν₂ fundamental;
* ``window`` 1900-3100 — the transparency window between bend and stretch, where
  what absorption exists is collision-induced;
* ``stretch`` 3200-4000 — the ν₁/ν₃ fundamentals.

How to read the output (this is the interpretation the diagnostic exists to
settle):

* a ratio near-uniform across libration / bend / stretch points at the **dipole
  surface** — an overall intensity scale, affecting all bands alike;
* an excess confined to the **stretch** points at the classical dynamics plus the
  harmonic quantum correction, which act on the high-frequency modes;
* the **window** column is the collision-induced contribution and is expected to
  carry the largest ratio either way, since a line list has almost nothing there.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from constants import C_LIGHT, H_PLANCK, K_B

# np.trapz was removed in numpy 2.0 (renamed trapezoid); support both.
try:
    from numpy import trapezoid as _trapz
except ImportError:  # numpy < 2.0
    from numpy import trapz as _trapz

_C_CM_S = C_LIGHT * 100.0  # so that c·wn (wn in cm⁻¹) is a frequency in s⁻¹

BANDS: dict[str, tuple[float, float]] = {
    "libration": (300.0, 1000.0),
    "bend": (1400.0, 1800.0),
    "window": (1900.0, 3100.0),
    "stretch": (3200.0, 4000.0),
}


@dataclass(frozen=True)
class BandIntegral:
    """One source's band-integrated absorption at one (T, ρ)."""

    name: str
    temperature_k: float
    rho_si: float
    integrals: dict[str, float]  # band -> A_band

    @property
    def stretch_over_bend(self) -> float:
        """A(stretch)/A(bend) — a source's own internal band balance.

        Compared *within* a source this is independent of any overall intensity
        scale, so it separates "everything is brighter" from "the high-frequency
        modes are brighter".
        """
        return self.integrals["stretch"] / self.integrals["bend"]


def absorption_weight(wn_cm: np.ndarray, temperature_k: float) -> np.ndarray:
    """ν(1 − e^{−hcν/k_BT}), cm⁻¹ — the dipole-to-absorption prefactor."""
    wn = np.asarray(wn_cm, dtype=float)
    x = H_PLANCK * _C_CM_S * wn / (K_B * float(temperature_k))
    return wn * -np.expm1(-x)


def band_absorption(
    wn_cm: np.ndarray,
    kappa_cm2_g: np.ndarray,
    temperature_k: float,
    bands: dict[str, tuple[float, float]] | None = None,
) -> dict[str, float]:
    """A_band for each band, integrated on the spectrum's own grid.

    Each source is integrated on the grid it was delivered on (FP at 1 cm⁻¹,
    HITEMP at 0.01 cm⁻¹); the *limits* are identical — the integrand is
    interpolated onto the exact band edges — which is what makes the ratio
    meaningful. Resampling one source onto the other's grid would throw away the
    finer one's structure for no gain, since these are integrals, not point values.
    """
    bands = BANDS if bands is None else bands
    wn = np.asarray(wn_cm, dtype=float)
    kappa = np.asarray(kappa_cm2_g, dtype=float)
    order = np.argsort(wn)
    wn, kappa = wn[order], kappa[order]
    weight = absorption_weight(wn, temperature_k)

    integrand = kappa * weight

    out = {}
    for name, (lo, hi) in bands.items():
        if lo < wn[0] or hi > wn[-1]:
            raise ValueError(
                f"band {name} = [{lo:g}, {hi:g}] cm^-1 is not covered by the spectrum "
                f"(grid spans [{wn.min():.4g}, {wn.max():.4g}])"
            )
        mask = (wn > lo) & (wn < hi)
        if mask.sum() < 2:
            raise ValueError(
                f"band {name} = [{lo:g}, {hi:g}] cm^-1 contains < 2 spectral points "
                f"(grid spans [{wn.min():.4g}, {wn.max():.4g}])"
            )
        # Close the band at exactly [lo, hi] by interpolating the integrand there.
        # Without this the limits fall on whichever grid points happen to lie
        # inside, which is a different pair of limits for each source — FP is on a
        # 1 cm⁻¹ grid, HITEMP on 0.01 — and the ratio stops being like for like.
        edge_vals = np.interp([lo, hi], wn, integrand)
        x = np.concatenate(([lo], wn[mask], [hi]))
        y = np.concatenate(([edge_vals[0]], integrand[mask], [edge_vals[1]]))
        out[name] = float(_trapz(y, x))
    return out
