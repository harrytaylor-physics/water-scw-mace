"""MT_CKD water-vapour **self**-continuum, added to a line-list spectrum (A58).

Why. The manuscript's Fig. 2b ratios are first-principles against *bare* HITEMP:
line absorption only, no continuum. Every planetary radiative-transfer model adds
an empirical continuum on top of its line list, so a reviewer will ask for the
comparison against line list + continuum as well. This module supplies that
second reference.

What is implemented. MT_CKD 4.3's self-continuum, exactly as AER's own
``src/mt_ckd_h2o_module.f90`` computes it:

    C(ν, T) = C_ref(ν) · (T_ref / T)^n(ν)                          (Mlawer 2023 Eq. 2)
    α_cont(ν) = n_H2O · C(ν, T) · (n_H2O / n_ref) · R(ν, T)        (Eqs. 1, 3)
    κ_cont(ν) = α_cont(ν) / ρ                                      [cm² g⁻¹]

with C_ref and the exponents n(ν) read verbatim from AER's distributed reference
table (``tools/extract_mtckd_ref.py``), T_ref = 296 K, and n_ref the number
density of the reference state (1013 mb, 296 K). R(ν, T) = ν·tanh(hcν / 2k_BT) is
LBLRTM's RADFN radiation term (Eq. 4), reproduced branch for branch.

α_cont is **quadratic in density** (one factor of n_H2O is the absorber column,
the other the collision partner), so the mass opacity κ_cont is *linear* in ρ —
unlike line absorption, whose mass opacity is ρ-independent at fixed broadening.
That is the whole reason the continuum matters at envelope densities.

**Foreign continuum is identically zero here and that is not a modelling
choice.** AER's foreign term carries a factor (1 − x_H2O); the envelope is pure
steam, x_H2O = 1. Self-only is forced by the composition.

**The continuum is extrapolated, hard, in both variables.** MT_CKD is fitted to
terrestrial-atmosphere measurements: T ≈ 250-350 K and densities at or below the
1013 mb / 296 K reference (ρ_ref = 0.7418 kg m⁻³). It is used here at
T = 500-3000 K and ρ = 0.1-1000 kg m⁻³, i.e. up to ~1350× the reference density
and ~10× the reference temperature, with a power law in T fitted over a ~100 K
span carried over a ~2700 K one. Nothing about that is endorsed by AER. It is
done because it is what a planetary model does, and the point of the exercise is
to show the reader what that costs.

Gate (``validation/test_mtckd_continuum.py``): the coefficients this module
computes reproduce AER's own compiled ``cntnm`` v4.3 binary — whose printed
``WATER.COEF`` columns are the reference coefficients with and without the
radiation term — to better than 1e-4 over 100-4000 cm⁻¹ at 296, 650, 1000 and
2000 K. The 296 K row alone would be near-tautological; the hot rows exercise the
(T_ref/T)^n path end to end.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import numpy as np

from constants import K_B, MMW_H2O, N_A

# hc/k_B in cm K — LBLRTM's RADCN2, the constant its RADFN is written against.
RADCN2_CM_K = 1.4387752

# Default location of the extracted reference table, relative to ``data_root``.
MTCKD_REF_RELPATH = "mtckd/mtckd_v4.3_self_ref.npz"


@dataclass(frozen=True)
class MTCKDReference:
    """AER's MT_CKD self-continuum reference table, as distributed."""

    wavenumber_cm: np.ndarray  # reference grid, 10 cm⁻¹ spacing
    self_absco_ref: np.ndarray  # cm² molecule⁻¹, at the reference density
    self_texp: np.ndarray  # dimensionless exponent n(ν)
    ref_press_mb: float
    ref_temp_k: float
    title: str
    source_sha256: str

    @property
    def ref_number_density_cm3(self) -> float:
        """n_ref: number density of the reference state, molecules cm⁻³.

        AER's Fortran scales by ``rho_rat = (p/p_ref)·(T_ref/T)``, which *is*
        n/n_ref for an ideal gas. Steam at envelope densities is not ideal, so we
        form the number-density ratio directly (the physically meaningful one,
        and what Mlawer 2023 Eq. 3 means) and never touch pressure.
        """
        return self.ref_press_mb * 100.0 / (K_B * self.ref_temp_k) * 1.0e-6

    @property
    def ref_rho_si(self) -> float:
        """ρ_ref: mass density of the reference state, kg m⁻³ (≈ 0.7418)."""
        return self.ref_number_density_cm3 * 1.0e6 * MMW_H2O / N_A


@lru_cache(maxsize=4)
def load_reference(path: str | Path) -> MTCKDReference:
    """Load the extracted MT_CKD reference npz (cached by path)."""
    path = Path(path).expanduser()
    if not path.exists():
        raise FileNotFoundError(
            f"MT_CKD self-continuum reference table not found at {path}. Generate it "
            f"once from AER's distributed netCDF with\n"
            f"    python tools/extract_mtckd_ref.py --out {path}\n"
            f"(see AMENDMENTS A58)."
        )
    with np.load(path, allow_pickle=False) as d:
        return MTCKDReference(
            wavenumber_cm=np.asarray(d["wavenumber_cm"], dtype=float),
            self_absco_ref=np.asarray(d["self_absco_ref"], dtype=float),
            self_texp=np.asarray(d["self_texp"], dtype=float),
            ref_press_mb=float(d["ref_press_mb"]),
            ref_temp_k=float(d["ref_temp_k"]),
            title=str(d["title"]),
            source_sha256=str(d["source_sha256"]),
        )


def radiation_term(wn_cm: np.ndarray, temperature_k: float) -> np.ndarray:
    """LBLRTM RADFN: R(ν, T) = ν·tanh(hcν / 2k_BT), cm⁻¹ (Mlawer 2023 Eq. 4).

    Reproduces AER's three branches exactly, including the small-argument form
    (which is where the ``0.5`` the Fortran comments flag as "required" lives).
    """
    wn = np.asarray(wn_cm, dtype=float)
    x = wn * RADCN2_CM_K / float(temperature_k)  # = hcν / k_B T
    with np.errstate(over="ignore"):
        mid = wn * (1.0 - np.exp(-x)) / (1.0 + np.exp(-x))
    return np.where(x <= 0.01, 0.5 * x * wn, np.where(x <= 10.0, mid, wn))


def _xint(
    x_nodes: np.ndarray, y_nodes: np.ndarray, x_out: np.ndarray
) -> np.ndarray:
    """AER's XINT 4-point interpolation off a uniform grid, vectorized.

    MT_CKD interpolates its 10 cm⁻¹ reference coefficients onto the output mesh
    with ``myxint``; reproducing that scheme rather than substituting a linear
    interpolation keeps this module's spectra identical to the model's wherever
    the two are compared. Points outside the inner nodes fall back to linear.
    """
    dx = float(x_nodes[1] - x_nodes[0])
    n = x_nodes.size
    x_out = np.asarray(x_out, dtype=float)

    # j = index of the node at or below x_out, restricted so j-1 and j+2 exist.
    j = np.floor((x_out - x_nodes[0]) / dx).astype(np.int64)
    inner = (j >= 1) & (j <= n - 3)
    jc = np.clip(j, 1, n - 3)
    p = (x_out - (x_nodes[0] + dx * jc)) / dx

    c = (3.0 - 2.0 * p) * p * p
    b = 0.5 * p * (1.0 - p)
    b1 = b * (1.0 - p)
    b2 = b * p
    cubic = (
        -y_nodes[jc - 1] * b1
        + y_nodes[jc] * (1.0 - c + b2)
        + y_nodes[jc + 1] * (c + b1)
        - y_nodes[jc + 2] * b2
    )
    return np.where(inner, cubic, np.interp(x_out, x_nodes, y_nodes))


def self_continuum_coefficient(
    wn_cm: np.ndarray, temperature_k: float, ref: MTCKDReference
) -> np.ndarray:
    """C(ν,T)·R(ν,T): the reference-density coefficient, cm² molecule⁻¹.

    This is exactly the "with radiation field / self" column AER's ``cntnm``
    binary writes to ``WATER.COEF``, and is what the gate compares against. It
    carries no density scaling — that is applied by :func:`self_continuum_kappa`.
    """
    coeff_ref = ref.self_absco_ref * (ref.ref_temp_k / float(temperature_k)) ** ref.self_texp
    return _xint(ref.wavenumber_cm, coeff_ref, wn_cm) * radiation_term(wn_cm, temperature_k)


def self_continuum_kappa(
    wn_cm: np.ndarray, temperature_k: float, rho_si: float, ref: MTCKDReference
) -> np.ndarray:
    """MT_CKD self-continuum **mass opacity** κ_cont(ν), cm² g⁻¹.

    α_cont = n·C·(n/n_ref)·R is quadratic in number density, so κ = α/ρ is linear
    in ρ. Returns a non-negative array (AER's reference coefficients are positive
    across 100-4000 cm⁻¹; the clip only guards the table edges).
    """
    rho_g_cm3 = float(rho_si) * 1.0e-3  # kg m⁻³ -> g cm⁻³
    n_cm3 = rho_g_cm3 * N_A / (MMW_H2O * 1.0e3)  # molecules cm⁻³
    coeff = self_continuum_coefficient(wn_cm, temperature_k, ref)
    alpha_cm = n_cm3 * coeff * (n_cm3 / ref.ref_number_density_cm3)  # cm⁻¹
    return np.maximum(alpha_cm / rho_g_cm3, 0.0)


# Registered continuum models, by the name a ``SourceSpec.continuum`` carries.
CONTINUUM_MODELS = {"mtckd_self": self_continuum_kappa}


def apply_continuum(
    name: str,
    wn_cm: np.ndarray,
    kappa_cm2_g: np.ndarray,
    temperature_k: float,
    rho_si: float,
    ref_path: str | Path,
) -> np.ndarray:
    """Add continuum model ``name`` to a loaded line-list spectrum (cm² g⁻¹)."""
    try:
        model = CONTINUUM_MODELS[name]
    except KeyError:
        raise KeyError(
            f"unknown continuum model {name!r}; registered: {sorted(CONTINUUM_MODELS)}"
        ) from None
    ref = load_reference(ref_path)
    return kappa_cm2_g + model(wn_cm, temperature_k, rho_si, ref)


# --------------------------------------------------------------------------
# Is the line/continuum decomposition even well posed here? (A58)
# --------------------------------------------------------------------------
# MT_CKD is not "extra absorption" in the abstract: it is defined as the residual
# beyond a line-by-line calculation done to LBLRTM's convention — every line
# truncated at LBLRTM_LINE_BOUND_CM from centre, with the Lorentz value at that
# bound ("the plinth") subtracted. That decomposition only means anything while
# the Lorentz half-width is small compared with the bound. These constants let a
# caller check whether it is, at the cell in question.
#
# gamma_self and its exponent are intensity-weighted means over 2.8e6 HITEMP H2O
# lines sampled across 150-4150 cm^-1 (the `selbrd` and `Tdpair` columns of the
# RADIS HITEMP2010 cache). Intensity weighting is the right one: strong lines set
# the opacity. The unweighted medians are 0.172 cm^-1/atm and 0.41, so this is the
# pessimistic end of a factor ~2.6 spread, and the estimate is an order-of-
# magnitude diagnostic, not a line-shape calculation.
LBLRTM_LINE_BOUND_CM = 25.0
HITEMP_GAMMA_SELF_CM_ATM = 0.455  # at 296 K, intensity-weighted
HITEMP_GAMMA_TEXP = 0.529
ATM_PA = 101325.0


def lorentz_hwhm_cm(temperature_k: float, pressure_pa: float) -> float:
    """Representative self-broadened Lorentz half-width, cm⁻¹.

    γ(T, P) = γ_self(296 K) · (P/1 atm) · (296/T)^n, the standard HITRAN form.
    """
    return float(
        HITEMP_GAMMA_SELF_CM_ATM
        * (float(pressure_pa) / ATM_PA)
        * (296.0 / float(temperature_k)) ** HITEMP_GAMMA_TEXP
    )


def decomposition_is_well_posed(temperature_k: float, pressure_pa: float) -> bool:
    """True where γ < LBLRTM's 25 cm⁻¹ line-by-line bound.

    Where it is False, "line list + MT_CKD" is not a decomposition of anything:
    the lines are broader than the truncation the continuum is defined against,
    and the continuum's own fit says nothing about that regime.
    """
    return lorentz_hwhm_cm(temperature_k, pressure_pa) < LBLRTM_LINE_BOUND_CM
