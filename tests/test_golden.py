"""Golden regression: pin model outputs to committed reference values.
First run writes tests/golden/*.json (COMMIT these); later runs compare.
Catches silent drift: swapped model files, LFS corruption, mace-torch
version changes altering evaluation."""
import numpy as np
from conftest import golden_check

def test_potential_golden(calc, water_dimer):
    a = water_dimer; a.calc = calc
    vals = {"energy_eV": float(a.get_potential_energy()),
            "fmax_eV_A": float(np.abs(a.get_forces()).max())}
    ref, first = golden_check("potential_dimer", vals, tol=1e-6)
    if first:
        import pytest; pytest.skip("golden pinned on first run — commit tests/golden/")
    assert abs(vals["energy_eV"] - ref["energy_eV"]) < 1e-6
    assert abs(vals["fmax_eV_A"] - ref["fmax_eV_A"]) < 1e-5

def test_dipole_golden(dipole_calc, water_dimer):
    from test_dipole import _mu
    mu = _mu(dipole_calc, water_dimer)
    vals = {"mu_norm": float(np.linalg.norm(mu))}
    ref, first = golden_check("dipole_dimer", vals, tol=1e-6)
    if first:
        import pytest; pytest.skip("golden pinned on first run — commit tests/golden/")
    assert abs(vals["mu_norm"] - ref["mu_norm"]) < 1e-5
