"""v5 keeps the sum rule and the leakage floor (A83).

A83 made the v5 spectra live: the A81 estimator (Hann-tapered periodogram,
``sum(w^2)`` normalisation, one harmonic correction, continuous ``omega^2``) at
the A82 smoothing width of 10 cm^-1. Two properties have to survive, and both
were bought at the cost of the whole A78-A81 bug hunt, so they are pinned:

  * the **sum rule** -- each band integral equals ``beta <|Mdot_band|^2>`` on the
    same tapered record, which is what says the normalisation and the single
    harmonic correction are right (A79, A80). Bound 1%; measured 0.09%.
  * the **leakage floor** -- the flat tail that A80 mistook for dipole-model
    noise must stay far below the transparency-window minimum, which is what
    says the taper is doing its job (A81). Bound 100x at every one of the 52
    cells; measured 476x at the worst, (3000 K, 108.796 kg m^-3).

If either fails, the estimator has changed and A79-A83 need re-reading.
"""

from __future__ import annotations

import numpy as np
import pytest

import opacity.dipole_spectrum as ds
from opacity.dipole_spectrum import K_B, two_sided_power_spectrum
from tools.v5_regenerate import DIPOLE_DIR, _series_dir

BANDS_CM = {"rotation": (50.0, 1100.0), "bend": (1000.0, 2400.0),
            "stretch": (2800.0, 4000.0)}
SUM_RULE_TOL = 0.01
FLOOR_MARGIN = 100.0
WINDOW_CM = (1900.0, 3100.0)
PRED = 1e4 / (12.0 * ds.EPS0 * ds.C_VAC * ds.C_CM)
SUM_RULE_CELLS = ((700.0, 30.0), (1000.0, 10.0), (2000.0, 30.0))


@pytest.fixture(scope="module")
def have_bundle(repo_root):
    if not (repo_root / DIPOLE_DIR).is_dir():
        pytest.skip(f"{DIPOLE_DIR} not present; A83's gates need the dipole bundle")
    return repo_root


def test_the_sum_rule_holds_to_one_percent(cfg, have_bundle):
    """int sigma dnu == beta <|Mdot_band|^2>, both on the tapered record."""
    worst = 0.0
    for t_k, rho in SUM_RULE_CELLS:
        with np.load(cfg.fp_opacity_path / f"T{t_k:g}_rho{rho:g}_opacity.npz") as z:
            n_mol, dt_fs = int(z["N"]), float(z["dt_fs"])
        m = np.load(_series_dir(have_bundle, t_k, rho) / "dipoles.npy")
        dt_s = dt_fs * 1e-15
        omega, itilde = two_sided_power_spectrum(m, dt_s)
        wn = omega / (2.0 * np.pi) / ds.C_CM
        sig = ((1.0 / (K_B * t_k)) * omega ** 2 * itilde
               / (6.0 * ds.EPS0 * ds.C_VAC * n_mol) * 1e4)
        mx = np.asarray(m, float) * ds.EANG_TO_DEBYE * ds.DEBYE_TO_CM
        mx = mx - mx.mean(0)
        n = len(mx)
        taper = np.hanning(n)
        for band, (lo, hi) in BANDS_CM.items():
            mask = (wn >= lo) & (wn <= hi)
            total = 0.0
            for k in range(3):
                f = np.fft.rfft(mx[:, k] * taper) * (1j * omega)
                f[~mask] = 0.0
                total += float(np.sum(np.fft.irfft(f, n=n) ** 2) / np.sum(taper ** 2))
            mdot = (1.0 / (K_B * t_k)) * total / n_mol
            got = float(np.trapezoid(sig[mask], wn[mask]))
            worst = max(worst, abs(got / (mdot * PRED) - 1.0))
    assert worst < SUM_RULE_TOL, f"worst sum-rule deviation {100*worst:.2f}%"


def test_the_leakage_floor_stays_far_below_the_window(cfg, have_bundle):
    """At every cell, the flat high-frequency tail is >= 100x below the window."""
    bad = {}
    for path in sorted(cfg.fp_opacity_path.glob("T*_opacity.npz")):
        with np.load(path) as z:
            t_k, rho = float(z["T"]), float(z["rho_kg_m3"])
            n_mol, dt_fs = int(z["N"]), float(z["dt_fs"])
            wn_out, kap = z["wn"].copy(), z["kappa"].copy()
        m = np.load(_series_dir(have_bundle, t_k, rho) / "dipoles.npy")
        omega, itilde = two_sided_power_spectrum(m, dt_fs * 1e-15)
        wn = omega / (2.0 * np.pi) / ds.C_CM
        sig = ((1.0 / (K_B * t_k)) * omega ** 2 * itilde
               / (6.0 * ds.EPS0 * ds.C_VAC * n_mol) * 1e4)
        edges = np.arange(6000.0, 32001.0, 2000.0)
        floor = min(float(np.median(sig[(wn >= lo) & (wn < hi)]))
                    for lo, hi in zip(edges[:-1], edges[1:]))
        _, kap_f = ds.sigma_to_alpha_kappa(np.full_like(wn_out, floor), rho / 1000.0)
        sel = (wn_out >= WINDOW_CM[0]) & (wn_out <= WINDOW_CM[1])
        margin = float(np.min(kap[sel]) / np.median(kap_f[sel]))
        if margin < FLOOR_MARGIN:
            bad[(t_k, rho)] = margin
    assert not bad, f"leakage floor within {FLOOR_MARGIN}x of the window at {bad}"
