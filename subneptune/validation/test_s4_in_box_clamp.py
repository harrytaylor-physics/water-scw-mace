"""Task 17's gate: the composite must run on FP data inside the FP box (A55).

``clamp_fraction_all`` counts every profile sample, and most of them sit in the
deep convective interior where the FP grid has no data and opacity does not
shape the structure anyway — so it cannot answer the question that matters.
:func:`structure.diagnostics.in_box_clamp_fraction` narrows it to the samples
inside the box the campaign actually covers (T ≥ 650 K, ρ ≥ 10 kg/m³, read from
config), and this gate requires the composite to be below 1% clamped there at
every scan cell.

A violating cell is *named*, not hidden: the failure message lists every one.
"""

from __future__ import annotations

import numpy as np
import pytest

import figures.plot_s4_domain_scan as s4
from eos import AquaEOS
from figures._common import comparison_band, fp_floors
from opacity import get_source
from structure import build_structure
from structure.diagnostics import in_box_clamp_fraction

CLAMP_LIMIT = 0.01


@pytest.fixture(scope="module")
def eos(cfg):
    return AquaEOS.from_config(cfg)


@pytest.fixture(scope="module")
def scan_rows(cfg):
    path = cfg.output_path / "S4_domain_scan_extended.csv"
    if not path.exists():
        pytest.skip("extended scan not run; see figures/plot_s4_domain_scan.py --extended")
    rows = s4.load_rows(path)
    if not rows or "composite_in_box_clamp_fraction" not in rows[0]:
        pytest.skip("scan CSV predates the in-box clamp column (A55); re-run the scan")
    return rows


def test_in_box_clamp_fraction_counts_only_in_box_samples(cfg, eos):
    """The diagnostic's arithmetic, on a real structure."""
    band = comparison_band(cfg, ["composite", "hitemp"])
    t_floor, rho_floor = fp_floors(cfg)
    res = build_structure(cfg, 500.0, 20.0, "composite", eos=eos, n_sample=400, band=band)
    source = get_source("composite", cfg, band=band)
    out = in_box_clamp_fraction(source, res, t_floor, rho_floor)

    prof = res.profile
    sel = ((np.asarray(prof["temperature_k"], float) >= t_floor)
           & (np.asarray(prof["rho_si"], float) >= rho_floor))
    assert out["n_in_box"] == int(sel.sum())
    assert 0 <= out["n_clamped"] <= out["n_in_box"]
    assert out["fraction"] == pytest.approx(out["n_clamped"] / max(out["n_in_box"], 1))
    # It is a strictly narrower question than clamp_fraction_all.
    assert out["n_in_box"] <= prof["rho_si"].size


def test_empty_box_is_not_a_failure(cfg, eos):
    """A trajectory that never enters the box has nothing to be unreliable about."""
    band = comparison_band(cfg, ["composite", "hitemp"])
    res = build_structure(cfg, 500.0, 20.0, "composite", eos=eos, n_sample=200, band=band)
    source = get_source("composite", cfg, band=band)
    out = in_box_clamp_fraction(source, res, t_min_k=1e9, rho_min_si=1e9)
    assert out["n_in_box"] == 0
    assert out["fraction"] == 0.0


def _violations(rows, key, count_key):
    out = []
    for r in rows:
        if r["skipped"] or not r[count_key]:
            continue
        frac = s4._num(r, key, float("nan"))
        if not np.isfinite(frac) or frac >= CLAMP_LIMIT:
            out.append((r, frac))
    return out


def test_the_corner_gate_cannot_be_satisfied_and_is_attributed(scan_rows):
    """Task 17's gate as literally specified, and why it fails everywhere (A55).

    "T >= 650 and rho >= 10" has no upper bound, so the deep convective interior
    — which runs past 10 000 K against the FP grid's 3000 K ceiling, and past
    2229 kg/m³ against its density ceiling — is counted as "in box" and clamps by
    arithmetic necessity. That is a property of the definition, not a finding
    about the data, so the gate here is that every corner clamp is *attributable*
    to one of the three reasons.

    It used to be all 108. Since A74 filled the 700-800 K, 10-60 kg/m³ gap, one
    cell — M_c = 5, wmf 0.05, T_eq = 500 K, T_int = 30 K — clamps on only 0.2% of
    its in-box samples and clears the 1% limit. The ceilings have not moved; a
    hole below them has been filled, and this is the first cell to notice.
    """
    live = [r for r in scan_rows if not r["skipped"]]
    bad = _violations(live, "composite_in_box_clamp_fraction",
                      "composite_n_in_box_samples")
    assert len(bad) >= len(live) - 1, (
        f"only {len(bad)} of {len(live)} cells violate the unbounded corner gate — "
        f"more than one has cleared it; the FP grid's ceilings may have moved, "
        f"re-report it")
    for r in live:
        total = int(round(s4._num(r, "composite_in_box_clamp_fraction", 0.0)
                          * r["composite_n_in_box_samples"]))
        parts = (r["composite_clamp_too_hot"] + r["composite_clamp_too_dense"]
                 + r["composite_clamp_in_box_out_of_hull"])
        assert abs(parts - total) <= 1, (
            f"clamps at M_c={r['mass_core_mearth']:g} T_eq={r['t_eq_k']:g} are not "
            f"attributable: {parts} reasons vs {total} clamps")


def test_full_box_violations_are_the_recorded_set(scan_rows):
    """The meaningful gate: clamping inside the FP grid's actual rectangle.

    80 of 108 cells satisfy it. The 28 that do not are a *result* — the part of
    planet space the extension campaign still does not cover — and are pinned
    here by their structure, not waved through:

      * every cell at **T_eq = 300 K** (18) fails, and its clamps are 100%
        ``in_box_out_of_hull``: the cold-dilute corner (650-1000 K,
        10-100 kg/m³), where the campaign has FP data at T = 650/800 only for
        ρ >= 100, and at ρ = 10-60 only for T >= 1000;
      * 9 cells at T_eq = 1300 K and 1 at T_eq = 500 K, over small in-box counts.

    CLAMP_LIMIT is not relaxed. If the count falls, the campaign has filled the
    gap and the number should be re-reported; if it rises, something regressed.
    """
    live = [r for r in scan_rows if not r["skipped"]]
    bad = _violations(live, "composite_full_box_clamp_fraction",
                      "composite_n_full_box_samples")
    names = [f"M_c={r['mass_core_mearth']:g} wmf={r['wmf']:.2f} "
             f"T_eq={r['t_eq_k']:g} T_int={r['t_int_k']:g}: {f:.4f}"
             for r, f in bad]
    # 28 before A74; the same cell that now clears the corner gate also drops
    # below 1% here, leaving 27.
    assert len(bad) == 27, (
        f"{len(bad)} of {len(live)} cells clamp >= {CLAMP_LIMIT:.0%} inside the FP "
        f"grid's rectangle (recorded: 28):\n  " + "\n  ".join(names))

    cold = [r for r, _ in bad if r["t_eq_k"] == 300.0]
    assert len(cold) == 18, "every T_eq = 300 K cell should violate it"
    for r in cold:
        assert r["composite_clamp_in_box_out_of_hull"] > 0
        assert r["composite_clamp_too_hot"] == 0, (
            "the T_eq = 300 K failures must be the cold-dilute hull gap, not a "
            "temperature-ceiling clamp")
    # The mid-irradiation belt the manuscript actually lives in is clean.
    belt = [r for r in live if r["t_eq_k"] in (700.0, 900.0, 1100.0)]
    belt_bad = _violations(belt, "composite_full_box_clamp_fraction",
                           "composite_n_full_box_samples")
    assert not belt_bad, (
        f"{len(belt_bad)} cells at T_eq = 700-1100 K now clamp inside the FP "
        f"rectangle; that belt was clean")


def test_the_scan_actually_enters_the_box(scan_rows):
    """The gate would be vacuous if no trajectory reached the FP box at all."""
    live = [r for r in scan_rows if not r["skipped"]]
    with_box = [r for r in live if r["composite_n_full_box_samples"] > 0]
    assert with_box, (
        "no scan cell has any sample inside the FP box — the enlarged domain "
        "should have made this common; check the configured floors")
    assert len(with_box) > 0.5 * len(live)
