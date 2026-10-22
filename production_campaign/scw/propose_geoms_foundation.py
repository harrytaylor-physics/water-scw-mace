#!/usr/bin/env python3
"""Phase 1 (bulk): propose geometries with a MACE foundation model -- NO DFT dynamics.

For each `bulk_finetune` cell in phase0_grid.csv, build a box at (T, rho), run short
NVT MD with a foundation MACE potential, and write decorrelated frames to a per-cell
pool.xyz. These pools are the geometries that Phase 2 single-point-labels at revPBE0-D3.

Requires (on an ARC htc GPU node): mace-torch, ase, torch+cuda.
    pip install mace-torch ase
"""
import argparse, csv, os, numpy as np
from ase import Atoms
from ase.io import write
from ase.md.langevin import Langevin
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase import units

NA = 6.02214076e23; M = 18.01528e-3  # kg/mol
rOH = 0.9572; half = np.deg2rad(104.52)/2
TEMPL = np.array([[0,0,0],
                  [ rOH*np.sin(half),0, rOH*np.cos(half)],
                  [-rOH*np.sin(half),0, rOH*np.cos(half)]])

def box_len(nmol, rho):                      # Angstrom
    return (nmol*M/NA/rho*1e30)**(1/3)

def rand_rot():
    q=np.random.normal(size=4); q/=np.linalg.norm(q); w,x,y,z=q
    return np.array([[1-2*(y*y+z*z),2*(x*y-z*w),2*(x*z+y*w)],
                     [2*(x*y+z*w),1-2*(x*x+z*z),2*(y*z-x*w)],
                     [2*(x*z-y*w),2*(y*z+x*w),1-2*(x*x+y*y)]])

def build_box(nmol, L, jitter=0.3, **_):
    """Grid-initialise nmol waters in an L^3 box (robust at any density; MD relaxes it)."""
    m = int(np.ceil(nmol ** (1/3)))
    s = L / m
    sites = [(i, j, k) for i in range(m) for j in range(m) for k in range(m)][:nmol]
    sym, pos = [], []
    for (i, j, k) in sites:
        c = (np.array([i, j, k]) + 0.5) * s + (np.random.rand(3) - 0.5) * jitter
        mol = (rand_rot() @ TEMPL.T).T + c
        sym += ["O", "H", "H"]; pos += list(mol)
    return Atoms("".join(sym), positions=np.array(pos), cell=[L, L, L], pbc=True)

def make_calc(model):
    from mace.calculators import mace_mp
    # model="medium" = MACE-MP-0; for high-pressure stability prefer MACE-MPA-0
    # (pass its downloaded path, e.g. model="/path/to/mace-mpa-0-medium.model").
    return mace_mp(model=model, dispersion=False, default_dtype="float32", device="cuda")

def run_cell(T, rho, nmol, equil_ps, sample_ps, dump_fs, calc, outdir):
    L=box_len(nmol, rho)
    min_oo = 2.4 if rho < 600 else 2.1
    at=build_box(nmol, L, min_oo=min_oo); at.calc=calc
    MaxwellBoltzmannDistribution(at, temperature_K=T)
    dt=0.5*units.fs
    dyn=Langevin(at, dt, temperature_K=T, friction=0.01/units.fs)
    dyn.run(int(equil_ps*1000/0.5))                       # discard equilibration
    outfile=os.path.join(outdir, f"pool_T{int(T)}_rho{rho:.3g}.xyz")
    every=int(dump_fs/0.5)
    n_dump=int(sample_ps*1000/0.5)
    def dump(a=at):
        write(outfile, a, append=True)
    dyn.attach(dump, interval=every)
    dyn.run(n_dump)
    return outfile, L

def auto_nmol(rho, L0, nmin, nmax):
    """Pick molecule count so the cubic box is ~L0 Angstrom, clamped to [nmin,nmax]."""
    n = rho * (L0 * 1e-10) ** 3 * NA / M
    return int(min(nmax, max(nmin, round(n))))

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--grid", default="phase0_grid.csv")
    ap.add_argument("--regime", default="bulk_finetune")
    ap.add_argument("--nmol", type=int, default=0,        # 0 = auto (variable N per density)
                    help="fixed molecule count, or 0 to size each cell to ~L0")
    ap.add_argument("--L0", type=float, default=14.0)     # target box edge (A)
    ap.add_argument("--nmol_min", type=int, default=16)
    ap.add_argument("--nmol_max", type=int, default=64)
    ap.add_argument("--max_L", type=float, default=28.0,  # above this -> too dilute, leave to clusters
                    help="skip a cell if even nmol_min gives a box larger than this")
    ap.add_argument("--equil_ps", type=float, default=5.0)
    ap.add_argument("--sample_ps", type=float, default=15.0)
    ap.add_argument("--dump_fs", type=float, default=25.0)   # decorrelated frames
    ap.add_argument("--model", default="medium")
    ap.add_argument("--out", default="pools")
    a=ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    calc=make_calc(a.model)
    rows=[r for r in csv.DictReader(open(a.grid)) if r["regime"]==a.regime]
    print(f"{len(rows)} {a.regime} cells; nmol={'auto' if a.nmol==0 else a.nmol}")
    for r in rows:
        T=float(r["T_K"]); rho=float(r["rho_kg_m3"])
        nmol = a.nmol if a.nmol > 0 else auto_nmol(rho, a.L0, a.nmol_min, a.nmol_max)
        L = box_len(nmol, rho)
        if L > a.max_L:
            print(f"  SKIP T={T:.0f} rho={rho:.3g}: even {nmol} mol gives L={L:.1f}A (dilute -> clusters)")
            continue
        f,L=run_cell(T, rho, nmol, a.equil_ps, a.sample_ps, a.dump_fs, calc, a.out)
        print(f"  T={T:.0f} rho={rho:.3g} nmol={nmol} L={L:.2f}A -> {f}")

if __name__=="__main__":
    main()
