"""Gates for the per-band absorption diagnostic (A61, Task 22).

The diagnostic decides one manuscript sentence — whether the dilute first-
principles intensity excess is an overall scale or a frequency-dependent tilt —
so both the integrand and the conclusion are pinned here.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from config import load_config
from opacity.bands import BANDS, absorption_weight, band_absorption

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def cfg():
    return load_config(ROOT / "config" / "fiducial.yaml")


@pytest.fixture(scope="module")
def data(cfg):
    import sys

    sys.path.insert(0, str(ROOT))
    from tools.band_diagnostic import rows

    return rows(cfg)


def test_bands_are_disjoint_and_inside_the_science_band():
    """Every band must sit inside [100, 4000], or the comparison is not like for like."""
    spans = sorted(BANDS.values())
    for (lo, hi), (nxt_lo, _) in zip(spans, spans[1:]):
        assert lo < hi <= nxt_lo  # ordered, non-overlapping
    assert min(lo for lo, _ in spans) >= 100.0
    assert max(hi for _, hi in spans) <= 4000.0


def test_absorption_weight_is_the_dipole_to_cross_section_prefactor():
    """ν(1 − e^{−hcν/kT}): → ν when hcν ≫ kT, → ν²hc/kT in the cold limit."""
    from constants import C_LIGHT, H_PLANCK, K_B

    hc_over_k = H_PLANCK * C_LIGHT * 100.0 / K_B  # cm K, CODATA
    wn = np.array([200.0, 1000.0, 3500.0])
    for t_k in (1000.0, 1500.0, 2000.0):
        x = wn * hc_over_k / t_k
        assert np.allclose(absorption_weight(wn, t_k), wn * (1.0 - np.exp(-x)), rtol=1e-12)
    assert absorption_weight(np.array([4000.0]), 50.0)[0] == pytest.approx(4000.0, rel=1e-6)
    small = absorption_weight(np.array([1.0]), 1.0e6)[0]
    assert small == pytest.approx(hc_over_k / 1.0e6, rel=1e-3)


def test_the_two_hc_over_k_constants_are_deliberately_different():
    """Pin the 1.2e-6 gap so nobody "unifies" them and breaks the MT_CKD gate.

    ``opacity/bands.py`` derives hc/k from CODATA via ``constants.py``, which is
    the project rule. ``opacity/mtckd.py`` hard-codes LBLRTM's RADCN2 because it
    must reproduce AER's own binary, and AER's constant is slightly stale.
    """
    from constants import C_LIGHT, H_PLANCK, K_B
    from opacity.mtckd import RADCN2_CM_K

    codata = H_PLANCK * C_LIGHT * 100.0 / K_B
    assert codata == pytest.approx(1.4387769, abs=1e-6)
    assert RADCN2_CM_K == 1.4387752
    assert abs(RADCN2_CM_K / codata - 1.0) == pytest.approx(1.2e-6, rel=0.3)


def test_band_absorption_is_positive_and_additive_over_a_split(cfg):
    """Integrating a band in two halves must equal integrating it whole."""
    import sys

    sys.path.insert(0, str(ROOT))
    from tools.band_diagnostic import _load

    wn, kappa = _load(cfg, "fp", 1000.0, 10.0)
    whole = band_absorption(wn, kappa, 1000.0, {"b": (1400.0, 1800.0)})["b"]
    halves = band_absorption(wn, kappa, 1000.0,
                             {"lo": (1400.0, 1600.0), "hi": (1600.0, 1800.0)})
    assert whole > 0
    assert whole == pytest.approx(halves["lo"] + halves["hi"], rel=1e-9)


def test_six_dilute_points_are_all_present(data):
    assert len(data) == 6
    assert {(r["T"], r["rho"]) for r in data} == {
        (t, rho) for t in (1000.0, 1500.0, 2000.0) for rho in (10.0, 30.0)
    }


def test_the_ratio_rises_monotonically_with_band_centre(data):
    """The finding: not an overall intensity scale but a frequency-dependent tilt.

    At every one of the six points the FP/line-list band ratio increases from the
    libration band through bend and window to the stretch, and crosses unity
    between bend and window. A uniform dipole-surface scale would give a flat row.
    """
    order = ("libration", "bend", "window", "stretch")
    for r in data:
        vals = [r["ratio"][b] for b in order]
        assert vals == sorted(vals), f"not monotone at T={r['T']}, rho={r['rho']}: {vals}"
        assert vals[0] < 1.0 < vals[-1]  # deficit at low frequency, excess at high


def test_fp_is_markedly_more_stretch_weighted_than_the_line_list(data):
    """The scale-free version of the same statement, source by source."""
    for r in data:
        assert r["stretch_over_bend_fp"] / r["stretch_over_bend_line"] > 3.0


def test_the_verdict_is_the_tilt_not_a_uniform_scale(data):
    import sys

    sys.path.insert(0, str(ROOT))
    from tools.band_diagnostic import _verdict

    verdict, monotone = _verdict(data)
    assert monotone
    assert "DIPOLE SURFACE" not in verdict
    assert "CLASSICAL DYNAMICS" in verdict and "MONOTONICALLY" in verdict
