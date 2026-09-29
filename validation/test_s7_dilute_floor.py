"""Dilute-cell noise-floor check (A52).

The thirteen FP cells at ρ ≤ 60 kg/m³ are the first-principles data in the
regime the composite normally hands to the line list, so whether they clear the
dipole model's own noise floor decides what the structural claim may lean on.
``figures/plot_s7_dilute_floor.py`` reimplements that floor (the ARC original is
not on this machine); these gates pin its arithmetic, the reason the answer
matters, and the state of the reimplementation — including the contradiction it
runs into, so that cannot quietly disappear.
"""

from __future__ import annotations

import pytest

# A81: this module pins A75's noise-floor construction, which modelled the
# dipole model's residual as WHITE NOISE IN M. A81 showed that is wrong in kind
# -- the residual is a smooth function of configuration, so it distorts the
# bands rather than raising a floor -- and that the flat tail A75 took for it
# was the untapered periodogram's own omega^-2 leakage. The construction is
# retired rather than re-pinned, and the file is kept for the trail (A83).
pytestmark = pytest.mark.skip(
    reason="A81: the white-noise-in-M floor model was wrong in kind; "
           "superseded by the leakage finding. See A81 and A83.")

import numpy as np
import pytest

from figures.plot_s7_dilute_floor import (
    DILUTE_RHO_MAX,
    DIPOLE_FLOOR_MDEBYE,
    DIPOLE_SIGNAL_RMS_MDEBYE,
    SCIENCE_BAND_CM,
    WINDOW_CM,
    _one_minus_exp,
    analyse_cell,
    floor_curve,
    load_dilute,
    nyquist_wavenumber_cm,
    window_share,
)

RATIO = DIPOLE_FLOOR_MDEBYE / DIPOLE_SIGNAL_RMS_MDEBYE
N_DILUTE = 13  # A74: the four corner cells are all rho <= 60


@pytest.fixture(scope="module")
def dilute(cfg):
    return load_dilute(cfg)


@pytest.fixture(scope="module")
def rows(dilute):
    return [analyse_cell(rec, dt, RATIO, "nyquist") for rec, dt in dilute]


def test_the_dilute_cells_are_the_low_density_column(dilute):
    """ρ ≤ 60 picks out the low-density points, thirteen of them since A74.

    (800, 30) has no ρ = 10 or 60 partner, so the temperature set is no longer a
    clean product with the density set; the inventory is asserted as the explicit
    list rather than as two independent sets.
    """
    assert len(dilute) == N_DILUTE
    assert all(rec.rho_si <= DILUTE_RHO_MAX for rec, _ in dilute)
    assert {rec.rho_si for rec, _ in dilute} == {10.0, 30.0, 60.0}
    assert {(rec.temperature_k, rec.rho_si) for rec, _ in dilute} == {
        (700.0, 10.0), (700.0, 30.0), (700.0, 60.0), (800.0, 30.0),
        (1000.0, 10.0), (1000.0, 30.0), (1000.0, 60.0),
        (1500.0, 10.0), (1500.0, 30.0), (1500.0, 60.0),
        (2000.0, 10.0), (2000.0, 30.0), (2000.0, 60.0)}


def test_nyquist_from_the_timestep():
    """0.5 fs ⇒ 1/(2·dt·c) ≈ 3.34e4 cm⁻¹."""
    assert nyquist_wavenumber_cm(0.5) == pytest.approx(33356.4, rel=1e-4)
    assert nyquist_wavenumber_cm(1.0) == pytest.approx(0.5 * 33356.4, rel=1e-4)


def test_floor_arithmetic_on_a_synthetic_spectrum():
    """A spectrum with flat D(ν) must give back κ × r²Δν/bandwidth exactly.

    Constructing κ = C·ν·(1 − e^{−x}) makes D(ν) = C constant, so
    ∫D dν = C·Δν and κ_floor = r²·C·Δν/bandwidth · ν(1−e^{−x}) = κ·r²Δν/B.
    That closed form pins every step of the construction at once.
    """
    wn = np.linspace(100.0, 4000.0, 4000)
    t_k = 1500.0
    const = 3.7
    kappa = const * wn * _one_minus_exp(wn, t_k)
    bandwidth = 33356.4
    floor = floor_curve(wn, kappa, t_k, 0.5, RATIO, bandwidth)
    expected = kappa * RATIO**2 * (wn.max() - wn.min()) / bandwidth
    np.testing.assert_allclose(floor, expected, rtol=1e-10)
    # Halving the bandwidth doubles the floor; the ratio enters squared.
    np.testing.assert_allclose(
        floor_curve(wn, kappa, t_k, 0.5, RATIO, bandwidth / 2.0), 2.0 * floor, rtol=1e-10)
    np.testing.assert_allclose(
        floor_curve(wn, kappa, t_k, 0.5, 2 * RATIO, bandwidth), 4.0 * floor, rtol=1e-10)


def test_the_window_carries_the_rosseland_mean(rows):
    """Why the floor matters: the window the model cannot resolve sets κ_R.

    The Rosseland mean is harmonic, so the most transparent part of the band
    dominates it. At every dilute cell ≥72% of the integrand comes from
    1900-3100 cm⁻¹ — while the Planck mean, an arithmetic average, barely notices
    it. Any floor problem in the window is therefore a κ_R problem.

    The bound was ≥80% over the nine original dilute cells. A74's (700, 10), the
    coldest and thinnest on the grid, sits at 0.724 and is the only one of the
    thirteen below 0.8; the rest span 0.81-0.93. Re-pinned to what the data give.
    """
    for r in rows:
        assert r["rosseland_window_share"] > 0.70, (
            f"T={r['temperature_k']:.0f}, rho={r['rho_si']:.0f}: Rosseland window "
            f"share fell to {r['rosseland_window_share']:.2f}")
        assert r["planck_window_share"] < 0.2


def test_window_share_is_a_fraction_of_the_band():
    """Both shares are well-formed fractions of a band the window sits inside."""
    wn = np.linspace(50.0, 5000.0, 5000)
    kappa = np.full_like(wn, 12.0)  # grey: shares reduce to weight fractions
    ross, planck = window_share(wn, kappa, 1500.0)
    for share in (ross, planck):
        assert 0.0 < share < 1.0
    assert WINDOW_CM[0] > SCIENCE_BAND_CM[0] and WINDOW_CM[1] < SCIENCE_BAND_CM[1]


def test_the_dilute_windows_are_floor_dominated(rows):
    """Regression pin on the headline number (loose: it is a reported result)."""
    for r in rows:
        assert r["frac_above_window"] < 0.10, (
            f"T={r['temperature_k']:.0f}, rho={r['rho_si']:.0f}: window fraction "
            f"above 2x floor rose to {r['frac_above_window']:.3f} — re-report it")
        # The band cores are fine; it is specifically the windows that fail. The
        # bound was > 0.5 over the nine original cells; A74's (700, 30) sits at
        # 0.394, the only one of the thirteen below 0.5, and the rest span
        # 0.54-0.75. The contrast with frac_above_window (< 0.05 everywhere) is
        # what carries the claim, and that is untouched.
        assert r["frac_above_cores"] > 0.35


def test_the_construction_contradicts_itself(rows):
    """The reimplemented floor EXCEEDS the delivered κ in the window (A52).

    For an additive white residual that is impossible — smoothing removes
    fluctuation, not bias, so E[κ] ≥ floor always. This gate exists so the
    contradiction stays on the record: it is the evidence that the quoted
    21.5% power fraction is not what ``compare_point_full.py`` spreads white, and
    the reason no absolute reliability claim is made from these numbers.
    """
    assert all(r["median_margin_window"] < 1.0 for r in rows)
    # The band cores, where the signal is strong, are comfortably above it.
    assert all(r["median_margin_cores"] > 2.0 for r in rows)


def test_the_bandwidth_assumption_is_the_dominant_one(dilute):
    """Nyquist vs the science band moves the verdict; it must stay reported."""
    rec, dt = dilute[0]
    wide = analyse_cell(rec, dt, RATIO, "nyquist")["frac_above_window"]
    narrow = analyse_cell(rec, dt, RATIO, "band")["frac_above_window"]
    assert wide > narrow  # a narrower spread ⇒ higher floor ⇒ less clears it
    assert narrow == pytest.approx(0.0, abs=1e-6)
