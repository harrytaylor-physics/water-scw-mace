"""First-principles opacity from a cell-dipole time series (A78-A80).

The production backend (``opacity_pipeline.spectrum_from_dipoles``, preserved
byte-identical in ``ed_data/arc_a78/scw/``) had one defect and one liability:

  * **line 82 applied the harmonic quantum correction twice** (A78). With
    ``Q_h = x/(1-e^-x)``, the standard relation
    ``alpha ~ w (1-e^-x) Itilde_qm`` and ``Itilde_qm = Q_h Itilde_cl`` already
    cancel to ``beta hbar w^2 Itilde_cl``; multiplying by ``Q_h`` again put a
    spurious, strongly temperature-dependent factor into every spectrum.
  * **lines 66 and 79** estimated the spectrum by windowing a *one-sided*
    autocorrelation with a *symmetric* ``np.hanning(L)`` and clipping the result
    at zero. A79 showed the band integrals that came out were a constant 3.52
    times too small, and -- because the clip rectifies noise upward where the
    true cross-section is smallest -- that no scalar can repair the lineshape:
    rescaling reproduces kappa_P to 2% but kappa_R only to 24%.

This module replaces the estimator rather than rescaling it. The normalisation
is unchanged, because A79 showed it was already right:

    sigma(w) = beta w^2 Itilde_2(w) / (6 eps0 c N)

with, stated once and used throughout: SI units; ``M`` the **total cell** dipole
and ``N`` the molecule count; ``C(t) = <M(0).M(t)>`` the **trace**, so the
orientational average of the field projection supplies the explicit 1/3 that
appears here as 1/6 together with the two-sided convention; ``Itilde_2`` the
**two-sided** transform of ``C``; refractive index 1; and ``sigma`` per molecule.
The harmonic correction enters **once**, implicitly, as the cancellation above.

``Itilde_2`` is estimated by the periodogram of the mean-removed series under a
**Hann taper**, normalised by the window's summed square so that Parseval is
preserved exactly (A81). The taper is not cosmetic. An untapered periodogram
treats the record as periodic, so the mismatch between its two ends acts as a
step discontinuity and leaks power as ``omega^-2`` -- which, multiplied by the
``omega^2`` of the absorption formula, is a **flat floor in sigma** right across
the transparency windows. A80 measured that leakage, mistook it for dipole-model
noise, and held propagation on it; A81 identified it and the taper removes it,
dropping the high-frequency floor by three to five orders of magnitude while
leaving every band integral unchanged to 0.09%. The zero clip is dropped: the
estimate is non-negative by construction, so the clip was inert in any case.

Everything downstream of the estimator is kept as delivered -- the 1 cm^-1
output grid, the Gaussian smoothing, the number density and the mass-opacity
conversion -- so a v3/v4 comparison isolates the estimator.
"""

from __future__ import annotations

import numpy as np
from scipy.ndimage import gaussian_filter1d

# CODATA 2018, matching the production backend exactly so that the only
# difference between v3 and v4 is the estimator.
K_B = 1.380649e-23
H_BAR = 1.054571817e-34
C_VAC = 2.99792458e8
C_CM = 2.99792458e10
EPS0 = 8.8541878128e-12
DEBYE_TO_CM = 3.33564e-30
EANG_TO_DEBYE = 1.602176634e-29 / DEBYE_TO_CM
N_AVOGADRO = 6.02214076e23
DEFAULT_MOLAR_MASS_G_MOL = 18.01528

MAX_WN_CM = 6000.0
WN_STEP_OUT = 1.0
ROSS_WN_MIN, ROSS_WN_MAX, ROSS_FLOOR_FRAC = 100.0, 5000.0, 1.0e-6


def two_sided_power_spectrum(m_eang, dt_s, taper=True):
    """Itilde_2(w) of the trace ACF, as a tapered periodogram (A81).

    With ``taper``, a Hann window is applied to the time series and the estimate
    is divided by ``sum(w**2)`` rather than ``n``, which is the normalisation
    that keeps Parseval exact for a tapered record. Without it the record's
    end-to-end mismatch leaks as omega^-2 and shows up as a flat floor in sigma
    (A80's "noise floor"); with it that floor falls by three to five orders of
    magnitude and the band integrals are unchanged.

    ``taper=False`` reproduces the A80 estimator and is kept for the A81
    diagnostic that identifies the leakage.
    """
    m = np.asarray(m_eang, dtype=float) * EANG_TO_DEBYE * DEBYE_TO_CM
    m = m - m.mean(0)
    n = len(m)
    w = np.hanning(n) if taper else np.ones(n)
    ft = np.fft.rfft(m * w[:, None], axis=0)
    itilde = dt_s * np.sum(np.abs(ft) ** 2, axis=1) / np.sum(w ** 2)
    omega = 2.0 * np.pi * np.fft.rfftfreq(n, d=dt_s)
    return omega, itilde


def sigma_to_alpha_kappa(sigma_cm2, density_g_cm3,
                         molar_mass_g_mol=DEFAULT_MOLAR_MASS_G_MOL):
    """Per-molecule cross-section -> absorption coefficient -> mass opacity."""
    n_m3 = (density_g_cm3 * 1.0e6 / molar_mass_g_mol) * N_AVOGADRO
    alpha_cm_inv = (n_m3 * sigma_cm2 * 1.0e-4) / 100.0
    return alpha_cm_inv, alpha_cm_inv / density_g_cm3


def planck_and_rosseland(wn_cm, kappa_cm2_g, t_k):
    """The backend's own means, kept identical so v3 and v4 compare directly."""
    nu = wn_cm * C_CM
    h = 2.0 * np.pi * H_BAR
    bnu = (2.0 * h * nu ** 3 / C_VAC ** 2) / np.maximum(
        np.exp(h * nu / (K_B * t_k)) - 1.0, 1.0e-300)
    x = h * nu / (K_B * t_k)
    ex = np.exp(np.clip(x, None, 700))
    db = (2.0 * h * nu ** 3 / C_VAC ** 2) * (
        ex / np.maximum((ex - 1.0) ** 2, 1.0e-300)) * (x / t_k)
    planck = np.trapezoid(kappa_cm2_g * bnu, nu) / np.trapezoid(bnu, nu)
    m = (wn_cm >= ROSS_WN_MIN) & (wn_cm <= ROSS_WN_MAX)
    k_safe = np.maximum(kappa_cm2_g[m], ROSS_FLOOR_FRAC * np.nanmax(kappa_cm2_g[m]))
    ross = np.trapezoid(db[m], nu[m]) / np.trapezoid(db[m] / k_safe, nu[m])
    return float(planck), float(ross)


def spectrum_from_dipoles(m_eang, t_k, dt_fs, n_molecules, density_g_cm3,
                          smooth_sigma_cm=20.0,
                          molar_mass_g_mol=DEFAULT_MOLAR_MASS_G_MOL):
    """The fixed estimator. `m_eang` is the cell dipole in e.Angstrom."""
    dt_s = dt_fs * 1e-15
    beta = 1.0 / (K_B * t_k)
    omega, itilde = two_sided_power_spectrum(m_eang, dt_s)
    wn = omega / (2.0 * np.pi) / C_CM
    keep = (wn > 0.0) & (wn <= MAX_WN_CM)
    # a single harmonic correction, already folded into the 1/6
    sigma = (beta * omega[keep] ** 2 * itilde[keep]
             / (6.0 * EPS0 * C_VAC * n_molecules)) * 1.0e4
    grid = np.arange(max(wn[keep].min(), 0.0), MAX_WN_CM + WN_STEP_OUT, WN_STEP_OUT)
    raw = np.nan_to_num(np.interp(grid, wn[keep], sigma, left=np.nan, right=np.nan),
                        nan=0.0)
    # no zero clip: a periodogram is non-negative, so it was inert here anyway
    smoothed = gaussian_filter1d(raw, smooth_sigma_cm / WN_STEP_OUT, mode="nearest")
    alpha, kappa = sigma_to_alpha_kappa(smoothed, density_g_cm3, molar_mass_g_mol)
    planck, ross = planck_and_rosseland(grid, kappa, t_k)
    return {"wn": grid, "sigma_raw": raw, "sigma": smoothed, "alpha": alpha,
            "kappa": kappa, "planck": planck, "rosseland": ross,
            "n_frames": int(np.asarray(m_eang).shape[0])}
