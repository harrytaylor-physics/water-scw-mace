"""The A80 pre-propagation gates on the fixed estimator.

A78 found the harmonic quantum correction applied twice; A79 showed the
normalisation was already right and that the residual factor was the spectral
estimator, which no scalar can repair because kappa_R is a lineshape quantity.
`opacity.dipole_spectrum` replaces the estimator (window-free periodogram, one
harmonic correction, no zero clip). A80 gates that replacement before anything
downstream is allowed to read it.

All three gates pass. The third was A80's blocker and A81 removed it: the floor
A80 measured was spectral leakage from an untapered periodogram, not dipole-model
noise, and a Hann taper on the time series drops it by three to five orders of
magnitude while leaving the band integrals exact. The margin is now ~1200x.

Propagation is still not done -- A81 was told not to -- so the registry still
points at v3 and nothing downstream has moved. A gate asserts that.
"""

from __future__ import annotations

import pytest

from tools.v4_regenerate import DIPOLE_DIR, FLOOR_MARGIN_REQUIRED, run


@pytest.fixture(scope="module")
def doc(cfg, repo_root):
    # The 41 MB dipole bundle is raw trajectory product and lives outside git
    # (SPEC §1), like ~/subneptune_data. Skip rather than fail without it.
    if not (repo_root / DIPOLE_DIR).is_dir():
        pytest.skip(f"{DIPOLE_DIR} not present; A80's gates need the dipole bundle")
    return run(cfg, write=False)


def test_every_delivered_cell_matches_a_dipole_series(doc):
    """All 52 resolve, and by reproduction rather than by filename."""
    m = doc["manifest"]
    assert not m["missing"], f"unmatched delivered cells: {m['missing']}"
    assert m["n_matched"] == m["n_delivered"] == 52
    assert m["worst_reproduction_error"] < 1e-4, (
        f"a series was matched only loosely ({m['worst_reproduction_error']:.1e}); "
        f"the manifest is no longer trustworthy")


def test_the_fixed_estimator_reproduces_the_classical_band_strength(doc):
    """Gate (2)(i): the A79 oscillator check, through the v4 estimator."""
    g = doc["gates"]["oscillator"]
    assert g["passes"], (
        f"the estimator gives {g['numeric_over_analytic']:.4f} of the textbook "
        f"classical band strength")


def test_band_integrals_equal_the_dipole_derivative_variance(doc):
    """Gate (2)(ii): int sigma dnu == beta<|Mdot_band|^2> to 2%, at the three
    rho = 30 cells. Mdot is taken spectrally, since central differences carry a
    [sin(w dt)/(w dt)]^2 bias that is a property of the test."""
    g = doc["gates"]["band_integrals"]
    assert g["passes"], f"worst deviation {100*g['worst_deviation']:.2f}% exceeds 2%"


def test_the_dilute_noise_floor_stays_ten_times_below_the_window(doc):
    """Was A80's blocker; A81 identified the cause and it now passes by ~1200x.

    The floor A80 measured was the untapered periodogram's omega^-2 leakage --
    the record's end-to-end mismatch acting as a step -- which the omega^2 of
    the absorption formula turns into a flat floor across the transparency
    windows. A Hann taper on the time series removes it. The floor here is now
    MEASURED from the spectrum's own flat tail rather than assumed from the lag
    test, because A75's white-noise-in-M model of the dipole residual was wrong
    in kind: a smooth function of configuration distorts the bands, it does not
    raise a floor.
    """
    g = doc["gates"]["noise_floor"]
    assert g["worst_margin_at_window_minimum"] >= FLOOR_MARGIN_REQUIRED, (
        f"margin {g['worst_margin_at_window_minimum']:.1f}x")


def test_v5_is_the_live_source(doc):
    """A83 propagated: the registry reads the v5 spectra, regenerated with the
    A81 estimator at the A82 width. Inverted from A80's
    ``test_nothing_has_been_propagated``, which held while the floor gate failed;
    A81 lifted that and A83 switched the registry."""
    from opacity.registry import get_spec

    assert get_spec("fp").subdir == "fp_opacity_v5", (
        f"the live first-principles source is {get_spec('fp').subdir!r}, not v5")
