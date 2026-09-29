# wbox.py -- shared helpers: build an N-water periodic box, total pressure, block SE.
# Imported by iso_pressure.py and finite_size.py so the box packer and pressure
# definition are identical everywhere (no copy-paste drift).
import numpy as np
from ase import Atoms
from scipy.spatial import cKDTree

Mw = 18.015e-3            # kg/mol
NA = 6.02214076e23        # 1/mol
EVA3_to_GPa = 160.21766   # 1 eV/A^3 in GPa

def box_length(N_mol, rho_gcc):
    """Cubic box edge (Angstrom) for N_mol waters at density rho (g/cc)."""
    return (N_mol * Mw / ((rho_gcc * 1000.0) * NA)) ** (1.0 / 3.0) * 1e10

def _wmol():
    d = 0.9572; a = np.deg2rad(104.52)
    return np.array([[0, 0, 0], [d, 0, 0], [d * np.cos(a), d * np.sin(a), 0]])

def _rot(rng):
    A = rng.standard_normal((3, 3)); Q, R = np.linalg.qr(A)
    Q = Q @ np.diag(np.sign(np.diag(R) + 1e-12))
    if np.linalg.det(Q) < 0: Q[:, 0] = -Q[:, 0]
    return Q

def build_box(N_mol, rho_gcc, seed=0, dmin=1.2, max_repair=200):
    """
    Place N_mol rigid waters on a shuffled cubic grid with random orientations,
    then repair any inter-molecular clash (< dmin A) by re-randomising the
    offending molecules. Works for N from 64 to many thousands (KD-tree clash
    check is O(N log N), not O(N^2)). Density is exact; the box is cubic.
    """
    L = box_length(N_mol, rho_gcc)
    n = 1
    while n ** 3 < N_mol: n += 1
    sp = L / n
    rng = np.random.default_rng(seed)
    sites = np.array([[(i + .5) * sp, (j + .5) * sp, (k + .5) * sp]
                      for i in range(n) for j in range(n) for k in range(n)])
    rng.shuffle(sites); sites = sites[:N_mol]

    mols = np.stack([(_wmol() @ _rot(rng).T) + sites[m] for m in range(N_mol)])  # (N_mol,3,3)

    def clashes(M):
        P = (M.reshape(-1, 3)) % L
        molid = np.repeat(np.arange(len(M)), 3)
        tree = cKDTree(P, boxsize=L)
        bad = set()
        for i, j in tree.query_pairs(dmin):
            if molid[i] != molid[j]:
                bad.add(molid[i]); bad.add(molid[j])
        return bad

    for _ in range(max_repair):
        bad = clashes(mols)
        if not bad: break
        for m in bad:
            mols[m] = (_wmol() @ _rot(rng).T) + sites[m]

    pos = mols.reshape(-1, 3)
    sym = ["O", "H", "H"] * N_mol
    at = Atoms(symbols=sym, positions=pos, cell=[L, L, L], pbc=True)
    return at, L

def total_pressure_GPa(atoms):
    """Thermodynamic pressure = virial (potential) + ideal-gas (kinetic), in GPa."""
    s = atoms.get_stress(voigt=True)                              # virial only, eV/A^3
    P_vir = -(s[0] + s[1] + s[2]) / 3.0
    P_kin = (2.0 / 3.0) * atoms.get_kinetic_energy() / atoms.get_volume()
    return (P_vir + P_kin) * EVA3_to_GPa

def block_se(Ps, nb=10):
    """Standard error of the mean from nb block averages (decorrelated estimate)."""
    blk = np.array_split(np.asarray(Ps), nb)
    bm = np.array([b.mean() for b in blk])
    return bm.std(ddof=1) / np.sqrt(nb)
