"""Gates for the density-scaling test, the wing cut-off, and the release tables.

Three of these check a *claim* rather than a number: that the line-list control
behaves as isolated lines must (A70), that the cut-off scan reproduces the
delivered spectra where it overlaps them (A71), and that the released table is
the model's own table rather than a second fit to it (A73). A number that moved
because the physics moved is a finding; one of these breaking means the argument
has lost its footing.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from config import load_config

ROOT = Path(__file__).resolve().parents[1]
RELEASE = ROOT / "data" / "release"


@pytest.fixture(scope="module")
def cfg():
    return load_config(ROOT / "config" / "fiducial.yaml")


def _out(cfg, name):
    p = cfg.output_path / name
    if not p.exists():
        pytest.skip(f"outputs/{name} not built on this machine")
    return json.loads(p.read_text())


# -- the rules, independent of whether the outputs exist -------------------
def test_lambda_is_a_density_and_nothing_else():
    """Lambda = 1 at 307 kg/m^3, recomputed here from the definition.

    tools/rho_scaling.py folds the constants into one expression; this rebuilds
    it the long way, so a slip in either is visible.
    """
    from constants import MMW_H2O, N_A
    from tools.rho_scaling import COLLISION_DIAMETER_M as d

    m = MMW_H2O / N_A
    for lam, rho in ((1.0, 306.7), (0.1, 3067.3)):
        n = rho / m
        mfp = 1.0 / (np.sqrt(2.0) * n * np.pi * d ** 2)
        assert mfp / d == pytest.approx(lam, rel=1e-3)


def test_the_index_is_two_for_exactly_pairwise_absorption():
    from tools.rho_scaling import _index

    # alpha = C rho^2  =>  q = alpha/rho^2 = C, flat, index 2
    assert _index(10.0, 3.0, 1000.0, 3.0) == pytest.approx(2.0)
    # alpha = C rho    =>  q ~ 1/rho, index 1
    assert _index(10.0, 1.0 / 10.0, 1000.0, 1.0 / 1000.0) == pytest.approx(1.0)


def test_mesa_log_r_uses_cgs_density():
    """log R = log10(rho[g/cm^3] / T6^3); the unit is the whole trap."""
    from tools.release_tables import log_r

    # 1 g/cm^3 at 10^6 K is R = 1 by construction
    assert float(log_r(1000.0, 1.0e6)) == pytest.approx(0.0)
    # and a steam envelope lands far outside any MESA table
    assert float(log_r(10.0, 1000.0)) == pytest.approx(7.0)


# -- the outputs ----------------------------------------------------------
def test_line_list_scaling_is_the_dilute_limit(cfg):
    """Bare HITEMP must be alpha ~ rho^1: isolated lines, no pair term.

    This is the control for A70. If it drifts from 1 the scaling machinery is
    measuring something other than what it claims to.
    """
    doc = _out(cfg, "rho_scaling.json")
    fit = doc["fitted_slope_dense"]["per_source"]["hitemp"]
    for key in ("window_1000_1200", "rosseland"):
        slopes = [v[key]["slope"] for v in fit.values()]
        assert min(slopes) > 0.9 and max(slopes) < 1.2, (
            f"bare HITEMP {key} index {min(slopes):.2f}-{max(slopes):.2f}, "
            f"not the rho^1 of isolated lines")


def test_continuum_scaling_is_pairwise_where_it_dominates(cfg):
    """HITEMP + MT_CKD must approach rho^2 at the cold end, rho^1 at the hot."""
    doc = _out(cfg, "rho_scaling.json")
    fit = doc["fitted_slope_dense"]["per_source"]["hitemp_mtckd"]
    cold = fit["500"]["rosseland"]["slope"]
    hot = fit["3000"]["rosseland"]["slope"]
    assert cold > 1.8, f"MT_CKD at 500 K has index {cold:.2f}, not pair-wise"
    assert hot < 1.2, f"MT_CKD at 3000 K has index {hot:.2f}, not line-dominated"


# The two diagnostics A70's claims rest on. The one-pair sigma has to be the
# conservative estimate for these, or the per-pair breakdown would be claiming
# more significance than the data support. It is NOT conservative for
# window_1000_1200 (0.118 against a pooled 0.141) — recorded, and harmless there
# because no breakdown density is quoted for that diagnostic.
_PAIR_MUST_BE_CONSERVATIVE = ("window_2550", "rosseland")


def test_the_pair_scatter_estimate_is_the_conservative_one(cfg):
    doc = _out(cfg, "rho_scaling.json")
    for key in _PAIR_MUST_BE_CONSERVATIVE:
        pair = doc["scatter"]["per_diagnostic"][key]["sigma_ln_alpha"]
        pooled = doc["scatter_pooled"]["per_diagnostic"][key]["sigma_ln_alpha"]
        assert pair >= pooled, (
            f"{key}: the pair estimate {pair:.3f} is below the pooled {pooled:.3f}, "
            f"so the breakdown densities quoted for it overstate their significance")


def test_where_the_pair_estimate_is_not_conservative(cfg):
    """Pinned so the exception stays visible rather than becoming folklore."""
    doc = _out(cfg, "rho_scaling.json")
    lenient = [k for k, v in doc["scatter"]["per_diagnostic"].items()
               if v["sigma_ln_alpha"]
               < doc["scatter_pooled"]["per_diagnostic"][k]["sigma_ln_alpha"]]
    assert lenient == ["window_1000_1200"], (
        f"the pair estimate is anti-conservative for {lenient}, not the single "
        f"diagnostic A70 records")


def test_wing_cutoff_reproduces_the_delivered_spectra(cfg):
    """The 50 cm^-1 column is the delivered convention and must match it."""
    doc = _out(cfg, "wing_cutoff.json")
    checked = 0
    for r in doc["results"]:
        if "reproduces_delivered_to" not in r:
            continue
        checked += 1
        # Exact (0.0) at four of the six points and 1e-9 at a fifth. The sixth,
        # (3000 K, 2228.7), differs by 1.0e-3 — and that is the reduced spectral
        # range, not the truncation: recomputed over the delivered
        # [50, 11000] cm^-1 the same column reproduces the delivered value
        # exactly. RADIS builds its lineshape-interpolation grid from the set of
        # loaded lines, so dropping 4600-11000 cm^-1 perturbs every line's
        # profile slightly, and most at the state whose widths span the widest
        # range. Measured, not assumed; see A71.
        assert r["reproduces_delivered_to"] < 2e-3, (
            f"T={r['T']:g}, rho={r['rho']:g}: the 50 cm^-1 column differs from the "
            f"delivered kappa_R by {r['reproduces_delivered_to']:.2e}")
    assert checked >= 1, "no point could be checked against the delivered spectra"


def test_wing_cutoff_raises_kappa_r(cfg):
    """A wider cut-off moves kappa_R up, at every point — the A71 direction."""
    doc = _out(cfg, "wing_cutoff.json")
    truncs = sorted(float(t) for t in doc["truncations_cm"])
    for r in doc["results"]:
        k = [r["by_truncation"][f"{t:g}"]["kappa_R_cm2_g"] for t in truncs
             if f"{t:g}" in r["by_truncation"]]
        assert all(b >= a for a, b in zip(k, k[1:])), (
            f"T={r['T']:g}, rho={r['rho']:g}: kappa_R is not monotone in the "
            f"cut-off: {k}")


def test_aqua_correction_leaves_density_and_gradient_alone(cfg):
    """The stated scope of the entropy fix, checked rather than taken on trust."""
    doc = _out(cfg, "aqua_entropy_check.json")
    if "points" not in doc:
        pytest.skip(doc.get("blocked", "no points"))
    for r in doc["points"]:
        assert abs(r["rho_check_frac"]) < 1e-3
        assert r["grad_ad_median_abs_frac_change"] < 5e-3


# -- the release ----------------------------------------------------------
def _release(name):
    p = RELEASE / name
    if not p.exists():
        pytest.skip(f"{p.relative_to(ROOT)} not built "
                    f"(python tools/release_tables.py config/fiducial.yaml)")
    return p


@pytest.mark.parametrize("stem", ("subneptune_kappa_R_mesa", "subneptune_kappa_P_mesa"))
def test_release_table_parses_as_mesa_expects(stem):
    """Header fields and row widths, per kap/private/load_kap.f90."""
    from tools.release_tables import MESA_FORM, MESA_VERSION, read_mesa

    d = read_mesa(_release(f"{stem}.txt"))
    assert d["form"] == MESA_FORM
    assert d["version"] >= MESA_VERSION
    assert d["log_kappa"].shape == (d["log_rs"].size, d["log_ts"].size)
    assert np.all(np.diff(d["log_rs"]) > 0) and np.all(np.diff(d["log_ts"]) > 0)
    assert d["log_rs"][0] == pytest.approx(d["header"][5])
    assert d["log_ts"][-1] == pytest.approx(d["header"][9])


def test_release_table_is_the_models_own_table(cfg):
    """Every filled cell must equal what the structure model would read there.

    A release that disagrees with the paper's own model is worse than no release.
    """
    from opacity.build_tables import get_source
    from opacity.table import CM2_PER_G_TO_SI
    from tools.release_tables import BAND_CM, NO_DATA, read_mesa

    d = read_mesa(_release("subneptune_kappa_R_mesa.txt"))
    src = get_source("fp", cfg, band=BAND_CM)
    worst, n = 0.0, 0
    for i, lr in enumerate(d["log_rs"]):
        for j, lt in enumerate(d["log_ts"]):
            v = d["log_kappa"][i, j]
            if v == NO_DATA:
                continue
            t_k = 10.0 ** lt
            rho = 10.0 ** lr * (t_k / 1.0e6) ** 3 * 1000.0
            model = np.log10(src.kappa_R(t_k, rho) / CM2_PER_G_TO_SI)
            worst = max(worst, abs(model - v))
            n += 1
    assert n > 50, f"only {n} filled cells to check"
    assert worst < 1e-4, f"released table differs from the model by {worst:.2e} dex"


def test_release_never_extrapolates():
    """Cells outside the computed hull carry the sentinel, not a number."""
    from tools.release_tables import NO_DATA, read_mesa

    d = read_mesa(_release("subneptune_kappa_R_mesa.txt"))
    k = d["log_kappa"]
    assert (k == NO_DATA).any(), "no cell is marked as outside the hull"
    assert (k != NO_DATA).any(), "every cell is marked as outside the hull"
    # the filled region must be contiguous along each temperature row
    for j in range(k.shape[1]):
        filled = np.nonzero(k[:, j] != NO_DATA)[0]
        if filled.size:
            assert filled[-1] - filled[0] + 1 == filled.size, (
                f"row {j} has a hole: a convex hull cannot produce one")


def test_release_fit_matches_its_own_residual(cfg):
    """The published coefficients must reproduce the published RMS."""
    from opacity.build_tables import compute_point_means
    from tools.release_tables import BAND_CM, DENSE_RHO_MIN

    doc = json.loads(_release("subneptune_kappa_fit.json").read_text())
    t_pt, rho_pt, kr, _ = compute_point_means("fp", cfg, BAND_CM)
    for name, mask in (("dense", rho_pt >= DENSE_RHO_MIN),
                       ("dilute", rho_pt < DENSE_RHO_MIN)):
        f = doc["fit"][name]
        if "blocked" in f:
            continue
        pred = (f["kappa_0_cm2_g"]
                * (rho_pt[mask] / f["rho_0_kg_m3"]) ** f["a_rho_exponent"]
                * (t_pt[mask] / f["T_0_K"]) ** f["b_T_exponent"])
        rms = float(np.sqrt((np.log10(pred / kr[mask]) ** 2).mean()))
        assert rms == pytest.approx(f["rms_log10_residual"], rel=1e-6)
        assert int(mask.sum()) == f["n_points"]
