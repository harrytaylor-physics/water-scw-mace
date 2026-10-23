"""Opacity gates (SPEC §3, §7).

Covers:
  * G4 opacity-swap null test (opacity level): identical tables -> identical κ
    to machine precision. The structure-level G4 (identical structure) lives in
    the structure gate.
  * factor-1000 plausibility (SPEC §2): scaled FP means are plausible; a raw
    (unscaled) table trips the hard-fail assertion.
  * band-limited means: known grid-point values, positivity.
  * interpolation exactness at grid vertices.
  * extrapolation gate: clamp counting and reported fraction.
  * registry behaviour: unknown source errors; undelivered line-list sources
    raise a clear error.
"""

from __future__ import annotations

import numpy as np
import pytest

from opacity import build_source, compute_point_means, get_spec, resolve_band
from opacity.build_tables import _assert_plausible, detect_common_band
from opacity.loaders import find_spectra, load_spectrum
from opacity.means import planck_mean, rosseland_mean
from opacity.table import CM2_PER_G_TO_SI

# Fixed reference band for pinned numeric assertions (independent of auto-detect).
_REF_BAND = (50.0, 4500.0)


# --------------------------------------------------------------------------
# G4 — opacity-swap null test (opacity level)
# --------------------------------------------------------------------------
def test_g4_identical_tables_identical_opacity(cfg):
    a = build_source("fp", cfg)
    b = build_source("fp", cfg)  # independent build, identical data

    rng = np.random.default_rng(0)
    # Query points spanning the grid interior and outside (clamped) region.
    temperature = rng.uniform(900.0, 3100.0, size=500)
    rho = 10.0 ** rng.uniform(1.5, 3.5, size=500)

    np.testing.assert_array_equal(a.kappa_R(temperature, rho), b.kappa_R(temperature, rho))
    np.testing.assert_array_equal(a.kappa_P(temperature, rho), b.kappa_P(temperature, rho))


# --------------------------------------------------------------------------
# Factor-1000 plausibility (SPEC §2)
# --------------------------------------------------------------------------
def test_scaled_fp_planck_means_are_plausible(cfg):
    _, _, _, kappa_p = compute_point_means("fp", cfg, resolve_band(cfg))
    assert kappa_p.min() >= 1e1 * 0.5  # comfortably inside the ~1e1-1e4 window
    assert kappa_p.max() <= 1e4 * 1.5
    # build_source runs the same assertion internally; must not raise.
    build_source("fp", cfg)


def test_raw_unscaled_means_trip_the_assertion(cfg):
    """The plausibility assertion still catches a factor-1000 error (SPEC §2).

    Since A50 the FP files are in correct units and the registry scale is 1.0,
    so there is no scale left to undo: the legacy error is re-introduced
    explicitly instead. The guard itself is unchanged — it is what would catch a
    legacy ×1000 file that somehow reached the table builder.
    """
    _, _, _, kappa_p = compute_point_means("fp", cfg, resolve_band(cfg))
    assert get_spec("fp").scale == 1.0  # A50: the correction lives in the data
    legacy = kappa_p * 1.0e3
    with pytest.raises(AssertionError, match="factor-1000"):
        _assert_plausible("fp", get_spec("fp"), legacy)


# --------------------------------------------------------------------------
# Band-limited means
# --------------------------------------------------------------------------
def test_known_grid_point_means(cfg):
    band = _REF_BAND  # pinned band for the reference values below
    spec = get_spec("fp")
    path = next(p for p in find_spectra(cfg.spectra_dir(spec.subdir), spec)
                if p.name.startswith("T1500_rho300"))
    rec = load_spectrum(path, spec, band)
    kp = planck_mean(rec.wn_cm, rec.kappa_cm2_g, rec.temperature_k, band)
    kr = rosseland_mean(rec.wn_cm, rec.kappa_cm2_g, rec.temperature_k, band)
    # Reference values computed independently at (T=1500, rho=300).
    assert kp == pytest.approx(565.6, rel=1e-3)  # A83: v5 (A81 estimator, sd = 10 cm^-1)
    assert kr == pytest.approx(83.737, rel=1e-3)  # A83: v5
    assert kr < kp  # Rosseland (harmonic) below Planck here


def test_means_positive_and_finite(cfg):
    _, _, kr, kp = compute_point_means("fp", cfg, resolve_band(cfg))
    assert np.all(np.isfinite(kr)) and np.all(kr > 0)
    assert np.all(np.isfinite(kp)) and np.all(kp > 0)


# --------------------------------------------------------------------------
# Interpolation exactness at grid vertices
# --------------------------------------------------------------------------
def test_interpolation_recovers_grid_values(cfg):
    src = build_source("fp", cfg)
    t = src.table.temperature_k
    rho = src.table.rho_si
    kr = src.kappa_R(t, rho)
    kp = src.kappa_P(t, rho)
    np.testing.assert_allclose(kr, src.table.kappa_rosseland_si, rtol=1e-9)
    np.testing.assert_allclose(kp, src.table.kappa_planck_si, rtol=1e-9)


# --------------------------------------------------------------------------
# Extrapolation gate (SPEC §3): clamp counting
# --------------------------------------------------------------------------
def test_clamp_counting(cfg):
    src = build_source("fp", cfg)
    src.reset_clamp_counter()
    # Three points inside the grid hull -> no clamps.
    src.kappa_R(np.array([1500.0, 2000.0, 2500.0]), np.array([300.0, 500.0, 800.0]))
    assert src.clamp_fraction == 0.0

    src.reset_clamp_counter()
    # Two points far outside the grid -> both clamped.
    src.kappa_R(np.array([200.0, 5000.0]), np.array([1.0, 1e5]))
    assert src.table.clamp.clamped_calls == 2
    assert src.clamp_fraction == 1.0


def test_clamp_uses_nearest_grid_value(cfg):
    src = build_source("fp", cfg)
    # Very high T and rho -> nearest vertex is the hottest, densest point.
    i = np.argmax(src.table.rho_si)
    kr = src.kappa_R(1e4, 1e5)
    assert kr == pytest.approx(src.table.kappa_rosseland_si[i], rel=1e-9)


# --------------------------------------------------------------------------
# Registry
# --------------------------------------------------------------------------
def test_unknown_source_raises():
    with pytest.raises(KeyError, match="unknown opacity source"):
        get_spec("bogus")


def test_line_list_sources_build(cfg):
    """POKAZATEL (27 points) and HITEMP (90) both build.

    HITEMP was 76 = 11 T × 7 ρ − 1 (T = 500, ρ = 0.1 was never produced) until A60
    added the T = 2250 and 2750 K rows at all seven densities, making 13 × 7 − 1.
    """
    pk = build_source("pokazatel", cfg)
    assert pk.n_points == 27
    ht = build_source("hitemp", cfg)
    assert ht.n_points == 90  # 13 T × 7 ρ − 1 (T=500, ρ=0.1 missing)
    t_grid = set(np.unique(ht.table.temperature_k).tolist())
    assert {2250.0, 2750.0} <= t_grid
    for src in (pk, ht):
        kr = src.kappa_R(1500.0, 300.0)  # a shared grid point
        assert kr > 0 and kr < 1e5


# --------------------------------------------------------------------------
# Band coverage guard
# --------------------------------------------------------------------------
def test_band_must_be_covered(cfg):
    spec = get_spec("fp")
    path = find_spectra(cfg.spectra_dir(spec.subdir), spec)[0]
    with pytest.raises(ValueError, match="does not cover"):
        load_spectrum(path, spec, (50.0, 7000.0))  # 7000 > FP max wn (~6000)


def test_si_conversion_factor():
    assert CM2_PER_G_TO_SI == 0.1  # 1 cm^2/g = 0.1 m^2/kg


# --------------------------------------------------------------------------
# Common-band auto-detection (A14)
# --------------------------------------------------------------------------
def test_detect_common_band(cfg):
    lo, hi = detect_common_band(cfg)
    # No args = all available. POKAZATEL is windowed to [50, 4500], narrower
    # than FP's [1.34, 6000], so FP∩POKAZATEL is [50, 4499.95] (A16).
    assert lo == pytest.approx(50.0, abs=0.1)
    assert hi == pytest.approx(4499.95, abs=1.0)
    assert lo < hi


def test_per_comparison_bands(cfg):
    """detect_common_band resolves per set of sources (A21): a windowed dataset
    only constrains comparisons it is part of. Independent of the config pin."""
    fp_only = detect_common_band(cfg, ["fp"])
    both_pk = detect_common_band(cfg, ["fp", "pokazatel"])
    both_ht = detect_common_band(cfg, ["fp", "hitemp"])

    assert fp_only[0] == pytest.approx(1.336, abs=0.05)  # FP's own full band
    assert fp_only[1] == pytest.approx(6000.11, abs=1.0)
    assert both_pk == pytest.approx((50.0, 4499.95), abs=1.0)  # POKAZATEL caps it
    # HITEMP does NOT cap FP at 4500; FP∩HITEMP = [50, 6000].
    #
    # The lower edge was 10 until A60 added the T = 2250 and 2750 K rows. Those
    # were computed from 50 cm⁻¹ because the 0-50 cm⁻¹ HITEMP line file is not in
    # the local cache and fetching it needs a HITRAN login, so the grid's
    # max-of-minima moved 10 -> 50. This is a band-detection change only: the
    # science band is pinned at [100, 4000] (A25), and with truncation = 50 a line
    # below 50 cm⁻¹ cannot reach past 100, so no mean moves. Re-running the two
    # rows with `--wavenum-min 10` where that file exists restores 10.
    assert both_ht[0] == pytest.approx(50.0, abs=1.0)
    assert both_ht[1] == pytest.approx(6000.0, abs=2.0)


def test_resolve_band_prefers_explicit(cfg):
    from dataclasses import replace

    # Config pins the science band [100,4000] (A25) -> resolve_band returns it.
    assert resolve_band(cfg, ["fp"]) == (100.0, 4000.0)
    # None -> auto-detect per source set.
    auto = replace(cfg, opacity=replace(cfg.opacity, band_cm=None))
    assert resolve_band(auto, ["fp"]) == detect_common_band(auto, ["fp"])
