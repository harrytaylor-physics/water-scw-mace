"""The A80 pre-propagation gates on the fixed estimator.

A78 found the harmonic quantum correction applied twice; A79 showed the
normalisation was already right and that the residual factor was the spectral
estimator, which no scalar can repair because kappa_R is a lineshape quantity.
`opacity.dipole_spectrum` replaces the estimator (window-free periodogram, one
harmonic correction, no zero clip). A80 gates that replacement before anything
downstream is allowed to read it.

Two of the three gates pass and are pinned here. The third fails, and is
recorded as a strict xfail rather than deleted or loosened: the A75 dilute noise
floor is **not** far below the transparency-window opacity once the estimator is
fixed. v3's clip at zero had been rectifying noise upward and the window at the
same time, so v3's window opacity was never a clean measurement and A75's
comparison against it was not well posed. v4 is the first version where the
question can be asked, and the answer is that at (1000 K, 10 kg m^-3) the floor
is 1.8x below the deepest point of the window -- 55% of what is there is
dipole-model noise -- against the 10x the brief requires.

When the floor is reduced (longer trajectories or a better dipole model, both
ARC jobs), this xfail flips to a pass and propagation can proceed. Until then
the registry still points at v3 and nothing downstream has moved.
"""

from __future__ import annotations

import pytest

from tools.v4_regenerate import FLOOR_MARGIN_REQUIRED, run


@pytest.fixture(scope="module")
def doc(cfg):
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


@pytest.mark.xfail(strict=True,
                   reason="A80: with the estimator fixed, the A75 dilute noise "
                          "floor is 1.8x below the window minimum at (1000, 10), "
                          "not the 10x required. v3's windows were never a clean "
                          "measurement, so A75's margin was not well posed. "
                          "Propagation is held until the floor is reduced.")
def test_the_dilute_noise_floor_stays_ten_times_below_the_window(doc):
    g = doc["gates"]["noise_floor"]
    assert g["worst_margin_at_window_minimum"] >= FLOOR_MARGIN_REQUIRED


def test_nothing_has_been_propagated(doc):
    """The registry still points at v3 and no downstream product has moved.
    If this ever fails, someone wired v4 in while the floor gate was failing."""
    from opacity.registry import get_spec

    assert doc["propagated"] is False
    assert get_spec("fp").subdir == "fp_opacity", (
        "the registry points at a v4 directory while the A80 floor gate fails")
