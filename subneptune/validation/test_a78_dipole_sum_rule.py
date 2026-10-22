"""The A78 source-level sum rule reproduces the delivered band integrals.

A78 read the production formula (`production_campaign/scw/opacity_pipeline.py`,
`spectrum_from_dipoles`) and found it applies the harmonic quantum correction
twice: line 82 multiplies `beta omega^2 Smm` -- which is already the
harmonic-corrected result -- by `Q_harmonic` again. The temperature trend A76
and A77 measured is that spare factor and nothing else.

This gate pins the chain that proves it, end to end:

  1. the formula as read reproduces the delivered `sigma_raw` to float
     precision, so the statement of it is not a paraphrase;
  2. `beta <|Mdot_band|^2>`, computed from `dipoles.npy` with no pipeline in
     between, satisfies the classical sum rule -- so the dipole series is not
     the source of the trend;
  3. the band integral of `sigma` with the spare `Q_harmonic` removed tracks
     `beta <|Mdot_band|^2>` at a temperature-independent ratio, which is the
     "reproduces the delivered band integrals" claim: the delivered integrals
     are recovered from the dipole series once the double correction is undone;
  4. the delivered bend integral still falls steeply, which is the defect
     itself.

Thresholds sit well outside the measured values: reproduction <= 1.7e-7 against
1e-5; |index| <= 0.196 against 0.3; ratio spread <= 5.4% against 10%; the
delivered bend index is -0.741 against -0.5. If (2) or (3) ever fails, the
dipole series or the pipeline has changed and A78 must be rewritten -- do not
loosen these.
"""

from __future__ import annotations

import pytest

import json
from pathlib import Path

from tools.dipole_sum_rule import (A78_MANIFEST, BANDS_CM, a78_cells_missing,
                                   a78_sha256_manifest, run)

REPRODUCTION_BOUND = 1.0e-5
DIPOLE_INDEX_BOUND = 0.3
RATIO_SPREAD_BOUND = 0.10
DELIVERED_BEND_BOUND = -0.5
ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def doc(cfg):
    # The five cells are raw production output, not committed (SPEC §1); like
    # A80's bundle they sit on disk beside the code. Skip rather than fail
    # without them -- but if they are there, they must be the fingerprinted ones.
    missing = a78_cells_missing(ROOT)
    if missing:
        pytest.skip(f"A78's cells are not present ({len(missing)} files missing, "
                    f"e.g. {missing[0]}); see {A78_MANIFEST}")
    expected = json.loads((ROOT / A78_MANIFEST).read_text())
    assert a78_sha256_manifest(ROOT) == expected, (
        f"the A78 cells on disk do not match {A78_MANIFEST}")
    return run(cfg)


def test_the_formula_as_read_reproduces_the_delivered_spectra(doc):
    """If this fails, A78's statement of the production formula is wrong."""
    bad = {c: r["max_rel_error"]["sigma_raw"] for c, r in doc["reproduction"].items()
           if r["max_rel_error"]["sigma_raw"] >= REPRODUCTION_BOUND}
    assert not bad, f"sigma_raw not reproduced to {REPRODUCTION_BOUND}: {bad}"


def test_the_dipole_series_satisfies_the_classical_sum_rule(doc):
    """beta <|Mdot|^2> is temperature-independent, so the defect is downstream
    of the dipoles -- in the formula, not the dynamics or the dipole model."""
    bad = {b: s["index_beta_mdot2"] for b, s in doc["summary"].items()
           if abs(s["index_beta_mdot2"]) >= DIPOLE_INDEX_BOUND}
    assert not bad, (
        f"the dipole series no longer satisfies the sum rule (|index| >= "
        f"{DIPOLE_INDEX_BOUND}): {bad}. A78 attributed the whole trend to the "
        f"spare Q_harmonic on the strength of this; rewrite it, do not relax.")


def test_the_band_integrals_are_recovered_from_the_dipole_series(doc):
    """With the spare Q_harmonic removed, the band integral of sigma is a
    temperature-independent multiple of beta <|Mdot_band|^2>. The constant
    absorbs the pipeline's window; only its constancy is the claim."""
    bad = {b: s["ratio_spread_frac"] for b, s in doc["summary"].items()
           if s["ratio_spread_frac"] >= RATIO_SPREAD_BOUND}
    assert not bad, (
        f"sigma/beta<|Mdot|^2> is no longer temperature-independent "
        f"(spread >= {RATIO_SPREAD_BOUND}): {bad}")


def test_the_delivered_bend_integral_still_carries_the_spare_correction(doc):
    """The defect itself. A silent upstream fix must not pass unnoticed."""
    s = doc["summary"]["bend"]
    assert s["index_sigma_delivered"] < DELIVERED_BEND_BOUND, (
        f"the delivered bend index is now {s['index_sigma_delivered']:.3f}; if "
        f"opacity_pipeline.py has been fixed, retire this gate and rewrite "
        f"A76/A77/A78 rather than loosening the bound.")
    assert set(doc["summary"]) == set(BANDS_CM)
