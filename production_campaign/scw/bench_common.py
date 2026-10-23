# bench_common.py -- shared helpers for run_traj / analyze_traj / eos / finite-size scripts.
import numpy as np

EVA3_to_GPa = 160.21766          # 1 eV/A^3 in GPa
Mw = 18.015e-3                    # kg/mol (water)
NA = 6.02214076e23

def box_length_A(nmol, rho_gcc):
    mass = nmol * Mw / NA
    V = mass / (rho_gcc * 1000.0)
    return V ** (1.0/3.0) * 1e10

def build_box(nmol, rho_gcc, seed=0):
    """nmol waters on a jittered grid with random orientations; equilibration melts it."""
    from ase import Atoms
    L = box_length_A(nmol, rho_gcc)
    n = int(np.ceil(nmol ** (1.0/3.0)))
    spacing = L / n
    rng = np.random.default_rng(seed)
    def rot():
        A = rng.standard_normal((3,3)); Q,R = np.linalg.qr(A)
        Q = Q @ np.diag(np.sign(np.diag(R) + 1e-12))
        if np.linalg.det(Q) < 0: Q[:,0] = -Q[:,0]
        return Q
    d = 0.9572; ang = np.deg2rad(104.52)
    wat = np.array([[0,0,0],[d,0,0],[d*np.cos(ang), d*np.sin(ang), 0]])
    centers = np.array([[(i+.5)*spacing,(j+.5)*spacing,(k+.5)*spacing]
                        for i in range(n) for j in range(n) for k in range(n)])[:nmol]
    pos = []; sym = []
    for c in centers:
        m = wat @ rot().T; m = m - m[0] + c
        pos.extend(m); sym += ["O","H","H"]
    return Atoms(symbols=sym, positions=np.array(pos), cell=[L,L,L], pbc=True), L

def make_calc(which, device="cuda", dtype="float32"):
    from mace.calculators import MACECalculator
    MODEL = {"v4":"water_scw_v4.model", "foundation":"mace-mpa-0-medium.model"}[which]
    return MACECalculator(model_paths=MODEL, device=device, default_dtype=dtype)

def pressure_GPa(atoms):
    """Total pressure = virial + ideal-gas (kinetic), in GPa."""
    s = atoms.get_stress(voigt=True)
    P_vir = -(s[0]+s[1]+s[2]) / 3.0
    P_kin = (2.0/3.0) * atoms.get_kinetic_energy() / atoms.get_volume()
    return (P_vir + P_kin) * EVA3_to_GPa

def block_se(a, nb=10):
    a = np.asarray(a, float)
    if len(a) < nb:
        return float(a.std() / np.sqrt(max(len(a),1)))
    bm = np.array([b.mean() for b in np.array_split(a, nb)])
    return float(bm.std(ddof=1) / np.sqrt(nb))

def rdf_OO(framesO, L, rmax=6.0, nbins=300, chunk=256):
    """O-O RDF from (nO,3) frames, minimum image, memory-safe for thousands of atoms."""
    rmax = min(rmax, L/2.0)
    edges = np.linspace(0.0, rmax, nbins+1); rc = 0.5*(edges[1:]+edges[:-1])
    hist = np.zeros(nbins)
    nfr = len(framesO); nO = framesO[0].shape[0]
    for o in framesO:
        for s in range(0, nO, chunk):
            a = o[s:s+chunk]
            dvec = a[:,None,:] - o[None,:,:]
            dvec -= L*np.round(dvec/L)
            r = np.sqrt((dvec*dvec).sum(-1)).ravel()
            r = r[(r > 1e-6) & (r < rmax)]
            hist += np.histogram(r, bins=edges)[0]
    numdens = nO / L**3
    shell = 4*np.pi*rc**2*(edges[1]-edges[0])
    return rc, hist / (nfr * nO * numdens * shell)

def iapws_P_GPa(T, rho_gcc):
    try:
        from iapws import IAPWS95
        return IAPWS95(T=T, rho=rho_gcc*1000.0).P * 1e-3
    except Exception:
        return None

def aqua_P_GPa(T, rho_gcc, path="/home/shug7909/scw/aqua_eos_rhot_v1_0.dat"):
    """Nearest-node AQUA pressure (GPa); returns (P, rho_node, T_node) or (None, None, None)."""
    try:
        A = np.loadtxt(path, skiprows=24, usecols=(0,1,2))   # rho[kg/m3] T[K] P[Pa]
        rho_kg = rho_gcc*1000.0
        d = (np.log10(A[:,0])-np.log10(rho_kg))**2 + (np.log10(A[:,1])-np.log10(T))**2
        i = d.argmin()
        return A[i,2]*1e-9, A[i,0], A[i,1]
    except Exception:
        return None, None, None
