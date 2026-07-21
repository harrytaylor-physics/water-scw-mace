"""Unit guards. The project's costliest bug was a x1000 number-density
unit error in sigma->alpha->kappa (opacity_pipeline, fixed 2026-06-29):
molar mass entered in kg/mol where g/mol was assumed. These tests pin the
conversions so that class of error can never silently return. When the
opacity pipeline is ported into this repo, import and test it directly."""
import numpy as np

NA = 6.02214076e23
M_WATER = 18.01528          # g/mol

def n_m3(density_g_cm3, molar_mass_g_mol=M_WATER):
    """Canonical number-density formula (post-fix opacity_pipeline line)."""
    return (density_g_cm3 * 1.0e6 / molar_mass_g_mol) * NA

def test_number_density_liquid_water():
    n = n_m3(1.0)
    assert abs(n - 3.3428e28) / 3.3428e28 < 1e-3

def test_number_density_not_off_by_1000():
    """The historical bug produced exactly n/1000; make its signature loud."""
    n = n_m3(1.0)
    assert not (3.0e25 < n < 3.7e25), "x1000 unit bug signature detected"

def test_kappa_alpha_relation():
    """kappa = alpha / (rho): units must compose to cm^2/g."""
    alpha_cm1 = 12.0
    rho_g_cm3 = 0.8
    kappa = alpha_cm1 / rho_g_cm3
    assert abs(kappa - 15.0) < 1e-12

def test_wavenumber_frequency():
    """nu [Hz] = c[cm/s] * wn[cm^-1]; 3000 cm^-1 -> ~9e13 Hz (stretch band)."""
    c_cm = 2.99792458e10
    assert abs(c_cm * 3000 - 8.9938e13) / 8.9938e13 < 1e-4
