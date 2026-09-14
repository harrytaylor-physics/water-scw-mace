"""Extended domain-scan gates (A47).

The A47 scan opens the planet space in every direction that could push a
line-list radiative zone into the FP data box — T_eq 300–1300 K, M_c 5–15 M⊕,
wmf 0.05–0.50. Three things have to hold for its table to mean anything:

  * the mass convention is the one stated (wmf ≡ f_env ≡ M_env/M_p, and
    M_p = M_c/(1−wmf)) — otherwise the water-world cells are a different planet
    from the one reported;
  * a cell that does not converge, or whose profile leaves the AQUA table, is
    *skipped and logged* rather than silently reported — AQUA clips instead of
    raising (A2), so nothing else would notice;
  * the CSV round-trips, since the scan costs ~35 min and is re-read, not re-run.
"""

from __future__ import annotations

import csv

import numpy as np
import pytest

import figures.plot_s4_domain_scan as s4
from eos import AquaEOS
from structure.diagnostics import aqua_coverage


@pytest.fixture(scope="module")
def eos(cfg):
    return AquaEOS.from_config(cfg)


def test_extended_grid_is_the_amended_one():
    """A47's grid, pinned: 6 x 2 x 3 x 3 = 108 cells."""
    assert s4.EXT_T_EQ_GRID == (300.0, 500.0, 700.0, 900.0, 1100.0, 1300.0)
    assert s4.EXT_T_INT_GRID == (30.0, 50.0)
    assert s4.EXT_M_CORE_GRID == (5.0, 10.0, 15.0)
    assert s4.EXT_WMF_GRID == (0.05, 0.20, 0.50)
    n = (len(s4.EXT_T_EQ_GRID) * len(s4.EXT_T_INT_GRID) * len(s4.EXT_M_CORE_GRID)
         * len(s4.EXT_WMF_GRID))
    assert n == 108
    # The A42 default grid is untouched by the generalization.
    assert s4.T_EQ_GRID == (500.0, 600.0, 700.0, 800.0)
    assert s4.T_INT_GRID == (20.0, 30.0, 40.0, 50.0, 75.0)
    assert s4.WMF_GRID == (0.05, 0.10)


def test_aqua_coverage_flags_what_the_eos_would_silently_clip(eos):
    """AquaEOS clips out-of-table queries; the scan must not be fooled (A47)."""
    t_lo, t_hi = eos.T_range_k
    inside = aqua_coverage(eos, np.array([1e5, 1e9]), np.array([500.0, 2000.0]))
    assert inside["inside"] and inside["summary"] == "inside the table"

    too_hot = aqua_coverage(eos, np.array([1e7]), np.array([10.0 * t_hi]))
    assert not too_hot["inside"] and too_hot["n_above_t"] == 1
    too_cold = aqua_coverage(eos, np.array([1e7]), np.array([0.5 * t_lo]))
    assert not too_cold["inside"] and too_cold["n_below_t"] == 1
    too_thin = aqua_coverage(eos, np.array([1e-9]), np.array([500.0]))
    assert not too_thin["inside"] and too_thin["n_below_p"] == 1
    # And the EOS really would have clipped rather than raised.
    assert np.isfinite(eos.state(1e-9, 500.0).rho)


def test_a_failed_cell_is_skipped_and_logged(cfg, monkeypatch, capsys):
    """A cell whose structure does not converge is recorded, not dropped."""
    def boom(*a, **kw):
        raise RuntimeError("shooting did not bracket R_p")

    monkeypatch.setattr(s4, "build_structure", boom)
    rows = s4.run_scan(cfg, eos=object(), t_eq_grid=(700.0,), t_int_grid=(30.0,),
                       wmf_grid=(0.05,), m_core_grid=(5.0,), n_sample=50)
    assert len(rows) == 1
    row = rows[0]
    assert row["skipped"] is True
    assert "did not converge" in row["skip_reason"]
    assert "hitemp" in row["skip_reason"]
    assert "SKIP" in capsys.readouterr().out
    # Blank source columns keep the CSV rectangular.
    for key in ("hitemp_radius_rearth", "composite_clamp_fraction_all"):
        assert row[key] == ""


def test_a_cell_off_the_aqua_table_is_skipped(cfg, monkeypatch):
    """The EOS clips silently, so an off-table profile must be caught here."""
    class FakeResult:
        profile = {"pressure_pa": np.array([1e5, 1e6]),
                   "temperature_k": np.array([1e9, 1e9])}
        radius_p_rearth = 2.0

    monkeypatch.setattr(s4, "build_structure", lambda *a, **kw: FakeResult())
    monkeypatch.setattr(s4, "aqua_coverage",
                        lambda *a, **kw: {"inside": False, "summary": "T too high"})
    rows = s4.run_scan(cfg, eos=object(), t_eq_grid=(700.0,), t_int_grid=(30.0,),
                       wmf_grid=(0.05,), m_core_grid=(5.0,), n_sample=50, verbose=False)
    assert rows[0]["skipped"] is True
    assert "outside the AQUA table" in rows[0]["skip_reason"]


def test_mass_convention_is_wmf_equals_f_env(cfg, monkeypatch):
    """A47: f_env is already M_env/M_p, so wmf ≡ f_env and M_p = M_c/(1−wmf)."""
    seen = []

    class FakeResult:
        profile = {"pressure_pa": np.array([1e5, 1e6]),
                   "temperature_k": np.array([500.0, 600.0]),
                   "rho_si": np.array([1.0, 2.0]), "mass_kg": np.array([2.0, 1.0]),
                   "is_convective": np.array([True, True])}
        radius_p_rearth = 2.0
        rcb_densities_si: list = []
        max_radiative_rho_si = 0.0
        clamp_fraction_all = 0.0
        clamp_fraction_radiative = 0.0

    def fake_build(cfg_c, t_eq, t_int, name, **kw):
        seen.append((cfg_c.core.mass_mearth, cfg_c.envelope.f_env))
        return FakeResult()

    monkeypatch.setattr(s4, "build_structure", fake_build)
    monkeypatch.setattr(s4, "aqua_coverage", lambda *a, **kw: {"inside": True})
    rows = s4.run_scan(cfg, eos=object(), t_eq_grid=(700.0,), t_int_grid=(30.0,),
                       wmf_grid=(0.50,), m_core_grid=(15.0,), n_sample=50, verbose=False)
    assert seen and all(mc == 15.0 and f_env == 0.50 for mc, f_env in seen)
    row = rows[0]
    assert row["wmf"] == row["f_env"] == 0.50
    assert float(row["mass_p_mearth"]) == pytest.approx(15.0 / 0.5)  # 30 M⊕


def test_csv_round_trips(tmp_path, cfg, monkeypatch):
    """The scan is re-read from CSV, not re-run (~35 min), so parsing is a gate."""
    class FakeResult:
        profile = {"pressure_pa": np.array([1e5]), "temperature_k": np.array([500.0]),
                   "rho_si": np.array([1.0]), "mass_kg": np.array([2.0, 1.0]),
                   "is_convective": np.array([True, True])}
        radius_p_rearth = 2.0
        rcb_densities_si: list = []
        max_radiative_rho_si = 0.0
        clamp_fraction_all = 0.25
        clamp_fraction_radiative = 0.0

    monkeypatch.setattr(s4, "build_structure", lambda *a, **kw: FakeResult())
    monkeypatch.setattr(s4, "aqua_coverage", lambda *a, **kw: {"inside": True})
    rows = s4.run_scan(cfg, eos=object(), t_eq_grid=(700.0, 900.0), t_int_grid=(30.0,),
                       wmf_grid=(0.05,), m_core_grid=(5.0,), n_sample=50, verbose=False)
    path = tmp_path / "scan.csv"
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    back = s4.load_rows(path)
    assert len(back) == 2
    for orig, r in zip(rows, back):
        assert r["t_eq_k"] == orig["t_eq_k"]
        assert r["skipped"] is False
        assert r["hitemp_n_in_box_layers"] == 0
        assert float(r["composite_clamp_fraction_all"]) == pytest.approx(0.25)

    grid = s4._grid_of(back, {"mass_core_mearth": 5.0, "wmf": 0.05},
                       "composite_clamp_fraction_all", (700.0, 900.0), (30.0,))
    assert grid.shape == (1, 2)
    assert np.allclose(grid, 0.25)


def test_rankings_ignore_skipped_cells():
    """Skipped rows must not win a top-ten place with a blank value."""
    live = {"skipped": False, "hitemp_in_box_radiative_mass_fraction": "0.002",
            "mass_core_mearth": 5.0, "wmf": 0.05, "t_eq_k": 700.0, "t_int_k": 30.0}
    dead = {**live, "skipped": True, "hitemp_in_box_radiative_mass_fraction": ""}
    ranked = s4._rank([dead, live], "hitemp_in_box_radiative_mass_fraction", 10)
    assert ranked == [live]
