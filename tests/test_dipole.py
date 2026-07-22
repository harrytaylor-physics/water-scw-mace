import numpy as np
from test_invariances import _rotation

def _mu(calc, atoms):
    atoms = atoms.copy()
    calc.calculate(atoms)            # DipoleMACE has no energy head
    for key in ("dipole", "dipole_moment"):
        if key in calc.results:
            return np.asarray(calc.results[key]).reshape(3)
    raise KeyError(f"no dipole in results: {list(calc.results)}")

def test_dipole_shape_and_finite(dipole_calc, water_dimer):
    mu = _mu(dipole_calc, water_dimer)
    assert mu.shape == (3,)
    assert np.isfinite(mu).all()
    assert 0.01 < np.linalg.norm(mu) < 50.0   # e*A-scale sanity for 2 waters

def test_dipole_rotational_equivariance(dipole_calc, water_dimer):
    """The physics test: rotate the structure, the dipole must co-rotate.
    mu(R x) = R mu(x). Catches corrupted/mis-exported dipole models."""
    mu1 = _mu(dipole_calc, water_dimer)
    R = _rotation(11)
    b = water_dimer.copy()
    b.positions = b.positions @ R.T
    b.cell = np.array(b.cell) @ R.T
    mu2 = _mu(dipole_calc, b)
    assert np.abs(mu2 - R @ mu1).max() < 1e-4

def test_dipole_translation_invariance(dipole_calc, water_dimer):
    mu1 = _mu(dipole_calc, water_dimer)
    b = water_dimer.copy(); b.positions += [1.3, 0.7, -2.1]; b.wrap()
    mu2 = _mu(dipole_calc, b)
    assert np.abs(mu2 - mu1).max() < 1e-4
