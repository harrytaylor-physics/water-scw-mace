"""FP units fixed at source, and the merged 48-point grid (A50).

``tools/fix_fp_units.py`` moved the legacy factor-1000 correction out of the
loader and into the data, and merged the 21 extension spectra. These gates pin
the three things that operation had to preserve:

(a) the κ the loader returns at (1000 K, 129.744 kg/m³) is bit-for-bit the value
    it returned before the rewrite — the numbers below were captured from the
    live pipeline immediately before ``fix_fp_units.py`` ran;
(b) the old and new spectra are mutually consistent where they nearly coincide
    in (T, ρ);
(c) all 48 files load with finite, positive κ and carry the units marker.

Plus the two invariants that keep the fix from silently regressing: the loader
REFUSES an unmarked file, and the rewritten files still satisfy α = n·σ (the
A50 evidence that σ/σ_raw were already correct and must not be scaled).
"""

from __future__ import annotations

import numpy as np
import pytest

from opacity.build_tables import compute_point_means
from opacity.loaders import load_spectrum
from opacity.registry import get_spec
from tools.fix_fp_units import (
    UNITS_MARKER,
    alpha_over_n_sigma,
    audit,
)

BAND_CM = (100.0, 4000.0)
N_FP_POINTS = 48  # 27 legacy + 21 extension

# (a) Captured from the pipeline immediately BEFORE the rewrite, on BAND_CM with
# rel_floor 1e-6. Once the data is fixed there is no live "before" left, so the
# baseline lives here.
BASELINE_POINT = (1000.0, 129.744)
BASELINE_KAPPA_R_CM2_G = 3.5409291039e01
BASELINE_KAPPA_P_CM2_G = 3.8343749591e02

# The new grid's domain (A50): the box the extension campaign opened up.
GRID_T_RANGE_K = (650.0, 3000.0)
GRID_RHO_RANGE_SI = (10.0, 2228.71)


@pytest.fixture(scope="module")
def means(cfg):
    """Per-point FP means on the science band: (T, ρ, κ_R, κ_P) in cm²/g."""
    return compute_point_means("fp", cfg, BAND_CM)


def _index_of(temps, rhos, point):
    i = int(np.argmin(np.abs(temps - point[0]) + np.abs(rhos - point[1])))
    assert temps[i] == pytest.approx(point[0])
    assert rhos[i] == pytest.approx(point[1])
    return i


def test_the_grid_is_the_merged_48_points(means):
    """27 legacy + 21 extension, spanning the enlarged domain."""
    temps, rhos, _, _ = means
    assert temps.size == N_FP_POINTS
    assert (temps.min(), temps.max()) == pytest.approx(GRID_T_RANGE_K)
    assert (rhos.min(), rhos.max()) == pytest.approx(GRID_RHO_RANGE_SI)


def test_a_loader_value_unchanged_by_the_rewrite(means):
    """Gate (a): the fix moved the ×1e-3 from the loader into the data, nothing else."""
    temps, rhos, kappa_r, kappa_p = means
    i = _index_of(temps, rhos, BASELINE_POINT)
    assert kappa_r[i] == pytest.approx(BASELINE_KAPPA_R_CM2_G, rel=1e-9)
    assert kappa_p[i] == pytest.approx(BASELINE_KAPPA_P_CM2_G, rel=1e-9)


def test_b_old_and_new_spectra_are_consistent(means):
    """Gate (b): κ_P(1000, 129.744) / κ_P(1000, 100) lands in [0.7, 1.3].

    A legacy point and an extension point 30% apart in density, at the same
    temperature: if either carried the wrong unit convention this ratio would be
    ~1e-3 or ~1e3, not order unity.
    """
    temps, rhos, kappa_r, kappa_p = means
    i = _index_of(temps, rhos, (1000.0, 129.744))  # legacy
    j = _index_of(temps, rhos, (1000.0, 100.0))  # extension
    assert 0.7 <= kappa_p[i] / kappa_p[j] <= 1.3
    assert 0.7 <= kappa_r[i] / kappa_r[j] <= 1.3


def test_c_every_spectrum_loads_finite_positive_and_marked(cfg):
    """Gate (c): all 48 files load, κ finite and > 0, marker present."""
    rows = audit(cfg.data_root_path / get_spec("fp").subdir)
    assert len(rows) == N_FP_POINTS
    for r in rows:
        assert r["marked"], f"{r['name']} has no {UNITS_MARKER!r} marker"
        assert r["n_nan"] == 0, f"{r['name']} has non-finite kappa"
        assert r["n_nonpos"] == 0, f"{r['name']} has non-positive kappa"
        assert r["kappa_lo"] > 0.0
    # Every file must still span the science band.
    assert max(r["wn_lo"] for r in rows) <= BAND_CM[0]
    assert min(r["wn_hi"] for r in rows) >= BAND_CM[1]


def test_the_loader_refuses_an_unmarked_spectrum(cfg, tmp_path):
    """A legacy ×1000 file can never be silently re-ingested."""
    spec = get_spec("fp")
    src = sorted((cfg.data_root_path / spec.subdir).glob(spec.glob))[0]
    with np.load(src) as d:
        payload = {k: d[k] for k in d.files if k != UNITS_MARKER}
    unmarked = tmp_path / src.name
    np.savez(unmarked, **payload)

    with pytest.raises(ValueError, match=UNITS_MARKER):
        load_spectrum(unmarked, spec, BAND_CM)
    # ... and the marked original still loads fine through the same path.
    assert load_spectrum(src, spec, BAND_CM).kappa_cm2_g.min() > 0.0


def test_no_on_load_scale_survives():
    """The correction lives in the data now, not in the registry (A50)."""
    assert get_spec("fp").scale == 1.0
    assert get_spec("fp").units_marker == UNITS_MARKER


def test_sigma_was_not_scaled(cfg):
    """α = n·σ in every rewritten file — the A50 decision, pinned.

    The factor-1000 entered at the α step, so σ/σ_raw were already correct
    cm²/molecule. Scaling them would have broken this identity and corrupted the
    cross-sections the dilute noise-floor check (A52) reads.
    """
    spec = get_spec("fp")
    paths = sorted((cfg.data_root_path / spec.subdir).glob(spec.glob))
    assert len(paths) == N_FP_POINTS
    for path in paths:
        with np.load(path) as d:
            assert alpha_over_n_sigma(d) == pytest.approx(1.0, rel=1e-6), (
                f"{path.name}: alpha/(n*sigma) != 1 — sigma and alpha no longer "
                f"share one unit convention"
            )
