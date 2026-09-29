"""The absolute normalisation of the dipole-ACF cross-section (A79).

A79 derived

    sigma(w) = beta w^2 Itilde_2(w) / (6 eps0 c N)      [two-sided transform]
             = beta w^2 Itilde_1(w) / (3 eps0 c N)      [one-sided]

and concluded that `opacity_pipeline.py` line 82, which writes the one-sided
form, has the **correct** prefactor -- so the factor of ~3.4 A78 left open is
the spectral estimator, not a normalisation.

The whole conclusion turns on the one-/two-sided factor of two, which is exactly
where this class of derivation goes wrong and which no amount of re-reading the
code can settle. So it is pinned numerically rather than by citation: a
thermostatted classical harmonic oscillator with a known dipole derivative and
reduced mass, integrated and compared against the textbook classical band
strength `A = (dmu/dq)^2 / (12 eps0 c^2 m)`.

If anyone later "corrects" the 6 to a 3 (or the 3 to a 6), this fails. The
tolerance is 3%, against a measured 0.6% error from the finite trajectory.
"""

from __future__ import annotations

import pytest

from tools.normalisation_check import synthetic_oscillator_check
# the estimator itself now lives in opacity.dipole_spectrum (A80);
# this check exercises the same expression and stays where A79 put it.

TOLERANCE = 0.03


@pytest.fixture(scope="module")
def check():
    return synthetic_oscillator_check()


def test_the_two_sided_form_divides_by_six(check):
    """beta w^2 Itilde_2 / (6 eps0 c N) reproduces the classical band strength."""
    got = check["two_sided_over_6"]
    assert abs(got - 1.0) < TOLERANCE, (
        f"beta w^2 Itilde_2/(6 eps0 c N) gives {got:.4f} of the textbook "
        f"classical band strength; A79's derivation, and with it the conclusion "
        f"that the production prefactor is correct, depends on this being 1.")


def test_dividing_by_three_would_double_it(check):
    """The one-sided prefactor is only correct with a one-sided transform; used
    with the two-sided one it is 2x too large. This pins the direction of the
    factor, so a sign-of-the-error slip cannot pass."""
    assert abs(check["two_sided_over_3"] - 2.0) < 2 * TOLERANCE
