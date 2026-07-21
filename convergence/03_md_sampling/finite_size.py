#!/usr/bin/env python3
"""Finite-size check: MACE NVT at several box sizes; mean pressure
(virial + kinetic) and g_OO first peak vs N.
Usage: finite_size.py model.model --nmols 32 64 128 256 --temp 1500 --rho 1200 --ps 5"""
import argparse, math
import numpy as np
from ase import units, Atoms
from ase.md.bussi import Bussi
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from mace.calculators import MACECalculator

M, NA = 18.01528, 6.02214076e23

def water_box(n_mol, rho, seed=0):
    L = (n_mol * M / NA / 1000.0 / rho) ** (1/3) * 1e10
    rng = np.random.default_rng(seed)
    r, ang = 0.9572, math.radians(104.52)
    base = np.array([[0,0,0],[r*math.sin(ang/2),0,r*math.cos(ang/2)],
                     [-r*math.sin(ang/2),0,r*math.cos(ang/2)]])
    ns = math.ceil(n_mol**(1/3)); sp = L/ns
    sym, pos, n = [], [], 0
    for i in range(ns):
        for j in range(ns):
            for k in range(ns):
                if n >= n_mol: continue
                c = np.array([(i+.5)*sp,(j+.5)*sp,(k+.5)*sp])
                q = rng.normal(size=4); q /= np.linalg.norm(q)
                w,x,y,z = q
                R = np.array([[1-2*(y*y+z*z),2*(x*y-w*z),2*(x*z+w*y)],
                              [2*(x*y+w*z),1-2*(x*x+z*z),2*(y*z-w*x)],
                              [2*(x*z-w*y),2*(y*z+w*x),1-2*(x*x+y*y)]])
                for el, v in zip("OHH", base):
                    sym.append(el); pos.append(c + R@v)
                n += 1
    return Atoms(sym, positions=pos, cell=[L,L,L], pbc=True)

def goo_peak(frames, L):
    hist = np.zeros(120); edges = np.linspace(0.5, min(6.0, L/2), 121)
    for a in frames:
        iO = [i for i,s in enumerate(a.get_chemical_symbols()) if s=="O"]
        dm = a.get_all_distances(mic=True)[np.ix_(iO,iO)]
        hist += np.histogram(dm[np.triu_indices(len(iO),1)], bins=edges)[0]
    r = 0.5*(edges[1:]+edges[:-1])
    shell = 4*np.pi*r**2*np.diff(edges)
    return float(r[np.argmax(hist/shell)])

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("model")
    ap.add_argument("--nmols", nargs="+", type=int, default=[32, 64, 128])
    ap.add_argument("--temp", type=float, default=1500.0)
    ap.add_argument("--rho", type=float, default=1200.0)
    ap.add_argument("--ps", type=float, default=5.0)
    ap.add_argument("--device", default="cuda")
    args = ap.parse_args()
    calc = MACECalculator(model_paths=args.model, device=args.device,
                          default_dtype="float32")
    print(f"{'N_mol':>6} {'P (GPa)':>16} {'g_OO peak (A)':>14}")
    for nm in args.nmols:
        atoms = water_box(nm, args.rho); atoms.calc = calc
        MaxwellBoltzmannDistribution(atoms, temperature_K=args.temp)
        dyn = Bussi(atoms, timestep=0.5*units.fs, temperature_K=args.temp,
                    taut=100*units.fs)
        ps_, frames = [], []
        def grab():
            s = atoms.get_stress(voigt=True); V = atoms.get_volume()
            ps_.append((-(s[0]+s[1]+s[2])/3 + len(atoms)*units.kB*args.temp/V)/units.GPa)
            frames.append(atoms.copy())
        dyn.attach(grab, interval=100)
        dyn.run(int(args.ps*1000/0.5))
        burn = len(ps_)//3
        p = np.array(ps_[burn:])
        print(f"{nm:6d} {p.mean():10.3f} +/- {p.std()/len(p)**.5:5.3f} "
              f"{goo_peak(frames[burn:], atoms.cell.lengths()[0]):14.3f}")

if __name__ == "__main__":
    main()
