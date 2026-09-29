"""Band-strength sum-rule gates (A77).

The integrated band strength of a classical trajectory is temperature-
independent: with `alpha ~ beta w^2 I_cl`, `int alpha dnu ~ beta <|mu_dot|^2>
= beta (kT/m) sum (dmu/dQ)^2`, and the line list's net (absorption minus
stimulated emission) band strength is T-independent by the same rule. Fitting
`int kappa dnu ~ T^p` over a window wide enough to hold a whole band must
therefore give p = 0.

Two gates, one per side of the comparison:

  * the line list must satisfy it -- this is the control, and if it ever fails
    the measurement itself is broken, not the first-principles spectra;
  * the first-principles vibrational bands must keep failing it, because A76/A77
    rest on that failure and a silent upstream fix would invalidate them
    without anything else noticing.

The second gate is pinned to the **delivered** kappa with no quantum correction
applied, i.e. the quantity the .npz files actually contain. A77 established that
the production ACF -> cross-section formula is not on this machine, and the two
candidates for it differ by Q_harmonic; gating the delivered array keeps the gate
independent of that open question. A row only counts if its temperatures span at least a factor of two, because a
log-log index over a narrower range is not a measurement of a trend. Two rows
fail that test and both are three hot points: rho = 100 (650-1000 K, 1.54x) and
rho = 800 (1500-2500 K, 1.67x). The second is the reason the rule is needed --
dropping its coldest point moves its index from -0.31 to **+0.07**, so it
carries no information about the trend either way. The six qualifying rows are
stable: p in [-0.825, -0.628] with a median of -0.715, and no single
leave-one-out fit rises above -0.535. Thresholds sit outside that with margin.
"""

from __future__ import annotations

import numpy as np
import pytest

from tools.sum_rule_trace import _index, _integral, _rows

BEND_CM = (1000.0, 2400.0)
MIN_T_SPAN = 2.0
LINE_LIST_BOUND = 0.1
FP_BEND_BOUND = -0.40
FP_BEND_MEDIAN_BOUND = -0.55


def _spans_enough(rows) -> bool:
    ts = [t for t, _ in rows]
    return max(ts) / min(ts) >= MIN_T_SPAN


def _bend_index(cfg, source, rows):
    ts = [t for t, _ in rows]
    vals = []
    for _, path in rows:
        with np.load(path) as d:
            vals.append(_integral(d["wn"], d["kappa"], *BEND_CM))
    return _index(ts, vals)


def test_the_line_list_obeys_the_band_strength_sum_rule(cfg):
    """HITEMP is the control: |p| < 0.1 in every density row with >= 3 points."""
    rows = {k: v for k, v in _rows(cfg, "hitemp").items() if _spans_enough(v)}
    assert len(rows) >= 6, f"expected the full HITEMP row set, got {len(rows)}"
    bad = {}
    for rho, r in rows.items():
        p = _bend_index(cfg, "hitemp", r)
        assert p is not None, f"rho = {rho}: no index"
        if abs(p) >= LINE_LIST_BOUND:
            bad[rho] = p
    assert not bad, (
        f"the line list must satisfy the sum rule it is the control for; "
        f"|p| >= {LINE_LIST_BOUND} at {bad}")


def test_the_first_principles_vibrational_bands_still_violate_it(cfg):
    """A76/A77's finding, pinned: every first-principles row fails the rule in
    the bend, and fails it by a wide margin in the median."""
    rows = {k: v for k, v in _rows(cfg, "fp").items() if _spans_enough(v)}
    assert len(rows) >= 6, (
        f"expected >= 6 first-principles rows spanning {MIN_T_SPAN}x in "
        f"temperature, got {len(rows)}")
    indices = {}
    for rho, r in rows.items():
        p = _bend_index(cfg, "fp", r)
        assert p is not None, f"rho = {rho}: no index"
        indices[rho] = p
    ok = {rho: p for rho, p in indices.items() if p >= FP_BEND_BOUND}
    assert not ok, (
        f"A76/A77 rest on the first-principles bend violating the sum rule; "
        f"p >= {FP_BEND_BOUND} at {ok}. If this is an upstream fix, retire the "
        f"gate and rewrite A76/A77 -- do not loosen the bound.")
    median = float(np.median(list(indices.values())))
    assert median < FP_BEND_MEDIAN_BOUND, (
        f"median bend index {median:.3f} is no longer below "
        f"{FP_BEND_MEDIAN_BOUND}; see above")
