import numpy as np
import pytest

def test_energy_force_consistency(calc, water_dimer):
    """F = -dE/dr is an identity in MACE (autograd); finite difference must agree."""
    a = water_dimer; a.calc = calc
    f0 = a.get_forces()[0, 2]
    h = 1e-4
    p = a.positions.copy()
    a.positions[0, 2] = p[0, 2] + h; ep = a.get_potential_energy()
    a.positions[0, 2] = p[0, 2] - h; em = a.get_potential_energy()
    a.positions[:] = p
    assert abs(-(ep - em) / (2 * h) - f0) < 5e-3

def test_dimer_binds_and_no_fusion(calc, water_dimer):
    a = water_dimer; a.calc = calc
    e_close = a.get_potential_energy()
    b = water_dimer.copy(); b.calc = calc
    b.positions[3:, 2] += 4.0
    e_far = b.get_potential_energy()
    assert e_close < e_far
    assert np.abs(a.get_forces()).max() < 20.0

def test_close_approach_repulsive_wall(calc, water_dimer):
    """Compressed O-H at 0.70 A (the deliberately-admitted proton-transfer
    regime): forces must be finite and energy must rise steeply."""
    a = water_dimer.copy(); a.calc = calc
    e_eq = a.get_potential_energy()
    b = water_dimer.copy(); b.calc = calc
    vec = b.positions[1] - b.positions[0]
    b.positions[1] = b.positions[0] + vec / np.linalg.norm(vec) * 0.70
    e_comp = b.get_potential_energy()
    f = b.get_forces()
    assert np.isfinite(f).all() and np.abs(f).max() < 1e3
    assert e_comp > e_eq + 1.0     # strongly uphill, no spurious well

def test_determinism(calc, water_dimer):
    a = water_dimer; a.calc = calc
    assert a.get_potential_energy() == a.get_potential_energy()

@pytest.mark.slow
def test_structure_regression(calc):
    """3 ps NVT at 673 K / 1000 kg/m3: nearest O-O in the liquid range."""
    import math
    from ase import units, Atoms
    from ase.md.langevin import Langevin
    from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
    rng = np.random.default_rng(1)
    n_mol, rho, T = 32, 1000.0, 673.0
    M, NA = 18.01528, 6.02214076e23
    L = (n_mol * M / NA / 1000 / rho) ** (1/3) * 1e10
    r, ang = 0.9572, math.radians(104.52)
    base = np.array([[0,0,0],[r*math.sin(ang/2),0,r*math.cos(ang/2)],
                     [-r*math.sin(ang/2),0,r*math.cos(ang/2)]])
    ns = math.ceil(n_mol ** (1/3)); sp = L/ns
    sym, pos, n = [], [], 0
    for i in range(ns):
        for j in range(ns):
            for k in range(ns):
                if n >= n_mol: continue
                c = np.array([(i+.5)*sp,(j+.5)*sp,(k+.5)*sp])
                for el, v in zip("OHH", base):
                    sym.append(el); pos.append(c+v+rng.normal(0,.05,3))
                n += 1
    a = Atoms(sym, positions=pos, cell=[L,L,L], pbc=True); a.calc = calc
    MaxwellBoltzmannDistribution(a, temperature_K=T)
    dyn = Langevin(a, timestep=0.5*units.fs, temperature_K=T, friction=0.02)
    frames = []
    dyn.attach(lambda: frames.append(a.copy()), interval=200)
    dyn.run(6000)
    goo = []
    for fr in frames[10:]:
        dm = fr.get_all_distances(mic=True)
        iO = [i for i,s in enumerate(fr.get_chemical_symbols()) if s=="O"]
        d = dm[np.ix_(iO,iO)][np.triu_indices(len(iO),1)]
        goo.append(d.min())
    assert 2.2 < np.median(goo) < 3.4
