import numpy as np

def _rotation(seed=3):
    rng = np.random.default_rng(seed)
    q = rng.normal(size=4); q /= np.linalg.norm(q)
    w, x, y, z = q
    return np.array([[1-2*(y*y+z*z), 2*(x*y-w*z), 2*(x*z+w*y)],
                     [2*(x*y+w*z), 1-2*(x*x+z*z), 2*(y*z-w*x)],
                     [2*(x*z-w*y), 2*(y*z+w*x), 1-2*(x*x+y*y)]])

def test_translation_invariance(calc, water_dimer):
    a = water_dimer; a.calc = calc
    e1 = a.get_potential_energy()
    a.positions += np.array([1.7, -0.4, 2.2]); a.wrap()
    assert abs(a.get_potential_energy() - e1) < 1e-6

def test_rotation_invariance(calc, water_dimer):
    """Rotate positions AND cell: energy identical, forces co-rotate."""
    a = water_dimer.copy(); a.calc = calc
    e1, f1 = a.get_potential_energy(), a.get_forces()
    R = _rotation()
    b = water_dimer.copy(); b.calc = calc
    b.positions = a.positions @ R.T
    b.cell = np.array(a.cell) @ R.T
    e2, f2 = b.get_potential_energy(), b.get_forces()
    assert abs(e2 - e1) < 1e-6
    assert np.abs(f2 - f1 @ R.T).max() < 1e-5

def test_permutation_invariance(calc, water_dimer):
    a = water_dimer.copy(); a.calc = calc
    e1 = a.get_potential_energy()
    order = [3, 4, 5, 0, 1, 2]          # swap the two molecules
    b = water_dimer[order]; b.calc = calc
    assert abs(b.get_potential_energy() - e1) < 1e-6

def test_supercell_extensivity(calc, single_water_small_box):
    """E(2x2x2 replica) = 8 * E(cell), with L/2 < r_max in the small cell —
    exercises multi-image neighbour handling as in the dense production
    boxes (L ~ 9.5 A < 2*r_max = 12 A)."""
    a = single_water_small_box; a.calc = calc
    e1 = a.get_potential_energy()
    b = a.repeat((2, 2, 2)); b.calc = calc
    e8 = b.get_potential_energy()
    assert abs(e8 - 8 * e1) < 1e-4
