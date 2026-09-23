"""Gates for the MT_CKD self-continuum and the ``hitemp_mtckd`` source (A58).

The central gate is not a self-consistency check: it compares this module's
coefficients against AER's own compiled ``cntnm`` v4.3 binary, whose printed
``WATER.COEF`` columns are the reference continuum coefficients with and without
the radiation term. Those reference outputs are checked in under
``validation/data/mtckd/`` so the gate runs anywhere, with
``tools/extract_mtckd_ref.py`` and ``validation/data/mtckd/README.md`` recording
how they were produced.

At 296 K the comparison is close to tautological — that is the state the
reference table is defined at — so it is run at 650, 1000 and 2000 K as well,
which exercises the (T_ref/T)^n(ν) path end to end. The binary prints 5
significant figures in single precision, so ~1e-4 is the floor of what the
comparison can resolve; the brief's gate is 1%.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from config import load_config
from constants import MMW_H2O, N_A
from opacity import get_source
from opacity.build_tables import compute_point_means
from opacity.mtckd import (
    load_reference,
    radiation_term,
    self_continuum_coefficient,
    self_continuum_kappa,
)
from opacity.registry import get_spec

BAND = (100.0, 4000.0)
REF_DIR = Path(__file__).parent / "data" / "mtckd"
# AER's cntnm prints 5 significant figures from a single-precision calculation.
CNTNM_TOL = 1.0e-3
# The brief's gate.
GATE_TOL = 1.0e-2


@pytest.fixture(scope="module")
def cfg():
    return load_config(Path(__file__).resolve().parents[1] / "config" / "fiducial.yaml")


@pytest.fixture(scope="module")
def ref(cfg):
    return load_reference(cfg.mtckd_ref_path)


def _read_water_coef(path: Path) -> np.ndarray:
    """(wn, self_no_rad, foreign_no_rad, self_with_rad, foreign_with_rad) rows.

    The file is a long provenance header followed by the coefficient block; the
    block is the final run of strictly-increasing 5-column numeric lines.
    """
    rows = []
    for line in path.read_text(errors="ignore").splitlines():
        fields = line.split()
        if len(fields) != 5:
            continue
        try:
            rows.append([float(x) for x in fields])
        except ValueError:
            continue
    arr = np.array(rows)
    i = len(arr) - 1
    while i > 0 and arr[i - 1, 0] < arr[i, 0]:
        i -= 1
    return arr[i:]


def _cntnm_cases():
    return sorted(REF_DIR.glob("WATER.COEF.T*"))


def test_the_reference_table_is_the_aer_distributed_one(ref):
    """The npz is MT_CKD 4.3, on AER's own 10 cm⁻¹ grid, at 1013 mb / 296 K."""
    assert "4.3" in ref.title
    assert ref.ref_temp_k == pytest.approx(296.0)
    assert ref.ref_press_mb == pytest.approx(1013.0)
    d_wn = np.diff(ref.wavenumber_cm)
    assert np.allclose(d_wn, 10.0)
    assert ref.wavenumber_cm.min() <= 0.0 and ref.wavenumber_cm.max() >= 20000.0
    # Reference density ~0.7418 kg/m³ — the density every envelope cell is
    # compared against, and the reason the extrapolation is as large as it is.
    assert ref.ref_rho_si == pytest.approx(0.7418, rel=1e-3)
    assert len(ref.source_sha256) == 64


def test_radiation_term_matches_the_analytic_tanh_away_from_the_branches():
    """R(ν,T) = ν·tanh(hcν/2k_BT) — RADFN's branches are an implementation detail."""
    from opacity.mtckd import RADCN2_CM_K

    wn = np.linspace(100.0, 4000.0, 4001)
    for t_k in (296.0, 650.0, 1000.0, 3000.0):
        x = wn * RADCN2_CM_K / t_k
        inner = (x > 0.02) & (x < 9.0)  # away from the branch switches
        analytic = wn * np.tanh(0.5 * x)
        assert np.allclose(radiation_term(wn, t_k)[inner], analytic[inner], rtol=1e-12)


@pytest.mark.parametrize("path", _cntnm_cases(), ids=lambda p: p.name.split(".")[-1])
def test_coefficients_reproduce_aers_compiled_cntnm_binary(ref, path):
    """The brief's gate: <1% against MT_CKD's own output, at 296-2000 K.

    Compared on AER's own 10 cm⁻¹ nodes, where their XINT and ours are both the
    identity, so this tests the physics (reference table, T dependence, radiation
    term) rather than an interpolation.
    """
    t_k = float(path.name.split("T")[-1])
    arr = _read_water_coef(path)
    wn, self_no_rad, self_with_rad = arr[:, 0], arr[:, 1], arr[:, 3]
    in_band = (wn >= BAND[0]) & (wn <= BAND[1]) & (self_with_rad > 0)
    assert in_band.sum() > 300

    coeff_ref = ref.self_absco_ref * (ref.ref_temp_k / t_k) ** ref.self_texp
    ours_no_rad = np.interp(wn[in_band], ref.wavenumber_cm, coeff_ref)
    ours_with_rad = self_continuum_coefficient(wn[in_band], t_k, ref)

    dev_no_rad = np.abs(ours_no_rad / self_no_rad[in_band] - 1.0).max()
    dev_with_rad = np.abs(ours_with_rad / self_with_rad[in_band] - 1.0).max()
    assert dev_no_rad < GATE_TOL and dev_with_rad < GATE_TOL
    # ...and in fact to the binary's print precision, which is the real result.
    assert dev_no_rad < CNTNM_TOL and dev_with_rad < CNTNM_TOL


def test_the_gate_state_is_the_briefs_296k_1bar_equivalent_density(ref):
    """At 296 K and ρ ≈ 0.7418 kg/m³ the density scaling is unity by definition.

    That is the brief's stated gate point; it is recorded here explicitly so the
    (much stronger) hot-temperature gates above are not mistaken for it.
    """
    wn = np.array([500.0, 1000.0, 2500.0, 3500.0])
    rho_ref = ref.ref_rho_si
    n_ref_cm3 = rho_ref * 1e-3 * N_A / (MMW_H2O * 1.0e3)
    assert n_ref_cm3 == pytest.approx(ref.ref_number_density_cm3, rel=1e-12)
    # κ_cont·ρ/n = C·R exactly when n = n_ref.
    kappa = self_continuum_kappa(wn, 296.0, rho_ref, ref)
    implied = kappa * (rho_ref * 1e-3) / n_ref_cm3
    assert np.allclose(implied, self_continuum_coefficient(wn, 296.0, ref), rtol=1e-12)


def test_continuum_mass_opacity_is_linear_in_density(ref):
    """α_cont ∝ ρ² (Mlawer Eq. 3) ⇒ κ_cont ∝ ρ. The whole reason it matters."""
    wn = np.linspace(200.0, 3800.0, 601)
    base = self_continuum_kappa(wn, 650.0, 100.0, ref)
    for factor in (2.0, 5.0, 12.0):
        assert np.allclose(self_continuum_kappa(wn, 650.0, 100.0 * factor, ref),
                           base * factor, rtol=1e-12)


def test_continuum_falls_steeply_with_temperature(ref):
    """(296/T)^n with n = 1.4-10.3 ⇒ the continuum is a cold-cell effect."""
    wn = np.array([1000.0])
    vals = [self_continuum_kappa(wn, t, 300.0, ref)[0] for t in (650.0, 1000.0, 2000.0)]
    assert vals[0] > vals[1] > vals[2]
    assert vals[0] / vals[2] > 100.0


def test_xint_choice_does_not_move_the_band_limited_means(ref):
    """AER's XINT is used for fidelity; quantify what the choice is worth.

    Pointwise the two schemes differ by up to ~11% at the sharpest features of
    the 10 cm⁻¹ reference grid — which is why XINT is reproduced rather than
    substituted — but the band-limited means the science path uses are insensitive
    to it at the sub-percent level. Both statements are asserted, because only the
    second licenses reading the result as physics.
    """
    from opacity.means import planck_mean, rosseland_mean

    wn = np.linspace(100.0, 4000.0, 39001)
    t_k = 650.0
    coeff_ref = ref.self_absco_ref * (ref.ref_temp_k / t_k) ** ref.self_texp
    rad = radiation_term(wn, t_k)
    xint = self_continuum_coefficient(wn, t_k, ref)
    linear = np.interp(wn, ref.wavenumber_cm, coeff_ref) * rad
    pointwise = np.abs(xint / linear - 1.0).max()
    assert 0.01 < pointwise < 0.15  # real, and bounded

    band = (100.0, 4000.0)
    for mean in (planck_mean, rosseland_mean):
        a = mean(wn, xint, t_k, band)
        b = mean(wn, linear, t_k, band)
        assert abs(a / b - 1.0) < 0.01


# -- the density scaling, end to end, against the same binary --------------
# cntnm takes a pressure and converts it with its own ideal-gas Loschmidt
# constant, so the state is specified by the pressure reproducing the cell's
# number density. See validation/data/mtckd/README.md.
_OPTDPT_CASES = {
    "CNTNM.OPTDPT.T1000_rho129.744": (1000.0, 129.744),
    "CNTNM.OPTDPT.T650_rho300": (650.0, 300.0),
}
# cntnm converts its input pressure to a column with its own hard-coded constants:
#   WTOT = ALOSMT * (P_mb/1013) * (273/T) * L
# Against CODATA (and 273.15 K) that column is low by ~3.1e-4, which is the whole
# of the residual below. The *density ratio* term cancels exactly — cntnm forms it
# from pressure as (P/1013)(296/T), and n/n_ref is the same expression — so this
# single factor, not the ρ² scaling, is what the comparison is left with.
CNTNM_ALOSMT = 2.68675e19  # molecules cm⁻³
CNTNM_COLUMN_BIAS = CNTNM_ALOSMT * 1e6 * 273.0 * 1.380649e-23 / (1013.0 * 100.0)


def _read_optdpt(path: Path) -> np.ndarray:
    rows = []
    for line in path.read_text(errors="ignore").splitlines():
        fields = line.split()
        if len(fields) != 2:
            continue
        try:
            rows.append([float(x) for x in fields])
        except ValueError:
            continue
    arr = np.array(rows)
    i = len(arr) - 1
    while i > 0 and arr[i - 1, 0] < arr[i, 0]:
        i -= 1
    return arr[i:]


def test_cntnms_own_column_constants_are_low_by_five_parts_in_ten_thousand():
    """Pin the known offset rather than absorb it into the tolerance below."""
    assert CNTNM_COLUMN_BIAS == pytest.approx(0.99969, abs=1e-5)


@pytest.mark.parametrize("name,state", sorted(_OPTDPT_CASES.items()),
                         ids=lambda v: v if isinstance(v, str) else "")
def test_density_scaling_reproduces_cntnm_at_envelope_densities(ref, name, state):
    """The gate the WATER.COEF files cannot give: alpha ∝ rho² at 175x and 404x rho_ref.

    With VMR = 1 and a 1 cm path, cntnm's printed optical depth is the self
    continuum alone (W_dry = 0, so no foreign and no O2/N2/CO2), i.e. exactly
    kappa_cont·rho·L. This exercises the number-density ratio, the temperature
    power law, the radiation term and the units in one comparison, at densities
    two to three orders of magnitude above anything MT_CKD was fitted to.
    """
    t_k, rho_si = state
    arr = _read_optdpt(REF_DIR / name)
    wn, tau = arr[:, 0], arr[:, 1]
    in_band = (wn >= BAND[0]) & (wn <= BAND[1]) & (tau > 0)
    assert in_band.sum() > 300

    assert rho_si / ref.ref_rho_si > 100.0  # the extrapolation this is testing
    ours = self_continuum_kappa(wn[in_band], t_k, rho_si, ref) * (rho_si * 1e-3)
    dev = ours / tau[in_band] - 1.0
    assert np.abs(dev).max() < GATE_TOL
    assert np.abs(dev).max() < CNTNM_TOL
    # The residual is a systematic offset, not scatter: it is positive at ~97% of
    # wavenumbers and its median is cntnm's own column-constant bias (3.1e-4); the
    # spread about it is the binary's single precision. If the density scaling were
    # wrong this would be a power of rho/rho_ref — 175x or 404x — not 3e-4.
    assert (dev > 0).mean() > 0.95
    assert np.median(dev) == pytest.approx(1.0 / CNTNM_COLUMN_BIAS - 1.0, abs=1e-4)


# -- the registered source ------------------------------------------------
def test_hitemp_mtckd_is_a_registry_entry_over_the_same_spectra():
    """Adding the source was a registry entry, not a new data delivery (SPEC §3)."""
    bare, cont = get_spec("hitemp"), get_spec("hitemp_mtckd")
    assert cont.subdir == bare.subdir and cont.glob == bare.glob
    assert cont.continuum == "mtckd_self" and bare.continuum is None
    assert cont.dilute_flat == bare.dilute_flat


def test_continuum_source_is_never_less_opaque_than_the_bare_line_list(cfg):
    """An added, positive continuum can only raise both means, at every point."""
    _, _, kr_bare, kp_bare = compute_point_means("hitemp", cfg, BAND)
    _, _, kr_cont, kp_cont = compute_point_means("hitemp_mtckd", cfg, BAND)
    assert np.all(kr_cont >= kr_bare * (1.0 - 1e-12))
    assert np.all(kp_cont >= kp_bare * (1.0 - 1e-12))


def test_the_continuum_dominates_the_cold_dense_corner_and_not_the_hot_one(cfg):
    """The result A58 rests on: the share is a strong function of (T, ρ)."""
    bare = get_source("hitemp", cfg, band=BAND)
    cont = get_source("hitemp_mtckd", cfg, band=BAND)
    share = lambda t, r: 1.0 - bare.kappa_R(t, r) / cont.kappa_R(t, r)  # noqa: E731
    assert share(650.0, 500.0) > 0.5  # cold + dense: continuum-dominated
    assert share(3000.0, 1000.0) < 0.15  # hot: the T power law has killed it
    assert share(1000.0, 10.0) < 0.05  # dilute: quadratic in ρ, so negligible


def test_the_line_continuum_decomposition_is_not_well_posed_on_most_of_the_grid(cfg):
    """MT_CKD is a residual to a 25 cm⁻¹-truncated line calculation. Check where
    that premise survives — this is a caveat on the result, so it is pinned.

    The self-broadened Lorentz half-width at the cell's own (ideal-gas) pressure
    is compared with LBLRTM's line-by-line bound. Only the thinnest cells are
    inside it; everywhere else the lines are broader than the truncation the
    continuum is defined against.
    """
    from opacity.crossover import crossover_v3_rows
    from opacity.mtckd import LBLRTM_LINE_BOUND_CM, lorentz_hwhm_cm

    assert LBLRTM_LINE_BOUND_CM == 25.0
    # sanity: the estimator is monotone in P and falls with T
    assert lorentz_hwhm_cm(650.0, 1e7) > lorentz_hwhm_cm(650.0, 1e6)
    assert lorentz_hwhm_cm(3000.0, 1e7) < lorentz_hwhm_cm(650.0, 1e7)

    rows = crossover_v3_rows(cfg, BAND)
    well_posed = [r for r in rows if r["decomposition_well_posed"]]
    assert 0 < len(well_posed) <= 5  # a handful of dilute cells, no more
    assert all(r["is_dilute"] for r in well_posed)
    # and the cold-dense corner, where the continuum matters most, is not among them
    cold = [r for r in rows if r["group"] == "cold_dense"]
    assert not any(r["decomposition_well_posed"] for r in cold)
    assert min(r["lorentz_hwhm_cm"] for r in cold) > LBLRTM_LINE_BOUND_CM


def test_line_strength_is_conserved_across_broadening_pressure(cfg):
    """RADIS renormalises truncated profiles — the truncation redistributes, not discards.

    This is the measurement behind A59's correction. RADIS cuts each line at
    ``truncation`` (50 cm⁻¹ here) and then divides by the truncated profile's own
    area, so the integrated line strength survives. If it discarded the wings
    instead, the integral of kappa would fall with pressure roughly as
    (2/pi)*arctan(truncation/HWHM) — a factor ~4 across this range, not 1e-4.

    Consequence, and the reason this is gated rather than left as a code reading:
    the high-pressure spectra are NOT lower bounds on absorption. They carry a
    line-shape systematic whose sign is undetermined, wherever the half-width
    exceeds the truncation.
    """
    import numpy as np
    from opacity.loaders import load_spectrum
    from opacity.registry import get_spec

    try:
        from numpy import trapezoid as trapz
    except ImportError:  # numpy < 2
        from numpy import trapz

    spec = get_spec("hitemp")
    spectra = cfg.data_root_path / spec.subdir
    integrals, pressures = [], []
    for rho in (0.1, 1.0, 10.0, 100.0, 1000.0):
        path = spectra / f"T1000_rho{rho:g}.npz"
        rec = load_spectrum(path, spec, BAND)
        integrals.append(trapz(rec.kappa_cm2_g, rec.wn_cm))
        with np.load(path) as d:
            pressures.append(float(d["p_bar"]))

    integrals = np.asarray(integrals)
    assert max(pressures) / min(pressures) > 1e3  # the range must be wide to bite
    assert np.ptp(integrals) / integrals.mean() < 1e-3


def test_a_spectrum_refuses_to_load_without_a_reference_table(cfg):
    """No silent fallback to bare line absorption if the table is missing."""
    from opacity.loaders import find_spectra, load_spectrum

    spec = get_spec("hitemp_mtckd")
    path = find_spectra(cfg.data_root_path / spec.subdir, spec)[0]
    with pytest.raises(ValueError, match="continuum"):
        load_spectrum(path, spec, BAND, mtckd_ref_path=None)
    with pytest.raises(FileNotFoundError, match="extract_mtckd_ref"):
        load_spectrum(path, spec, BAND, mtckd_ref_path=Path("/nonexistent/mtckd.npz"))
