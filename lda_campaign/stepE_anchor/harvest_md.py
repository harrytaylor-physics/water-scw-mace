#!/usr/bin/env python3
"""Step D-lite harvest: v1 MD at the dense-hot branch, 15 frames each -
v1 generates its own failure configurations (soft O-O wall)."""
import argparse, math, pathlib
import numpy as np
from ase import units, Atoms
from ase.io import write, read
from ase.md.bussi import Bussi
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from mace.calculators import MACECalculator

POINTS = [(2500, 1520.77), (2750, 1841.02), (3000, 1043.59),
          (3000, 1343.91), (3000, 1730.66), (3000, 2228.71)]
M, NA = 18.01528, 6.02214076e23

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", required=True)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--n-mol", type=int, default=64)
    ap.add_argument("--settle-ps", type=float, default=1.0)
    ap.add_argument("--run-ps", type=float, default=10.0)
    ap.add_argument("--n-frames", type=int, default=15)
    ap.add_argument("--subdir", default="harvest")
    args = ap.parse_args()
    calc = MACECalculator(model_paths=args.model, device="cuda",
                          default_dtype="float32")
    work = pathlib.Path(args.workdir)
    for T, rho in POINTS:
        tag = f"T{T}_rho{rho:g}"
        d = work / args.subdir / tag
        d.mkdir(parents=True, exist_ok=True)
        # start from the selected pool geometry for this point (exists from step A)
        src = sorted((work / "selected").glob(f"{tag}*/sel_0000.xyz")) or \
              sorted((work / "selected" / tag).glob("sel_0000.xyz")) or \
              sorted((work / "propose" / tag).glob("traj.xyz"))
        if src:
            atoms = read(src[0], index=-1)
        else:
            raise SystemExit(f"{tag}: no starting geometry found under selected/ or propose/")
        atoms.calc = calc
        MaxwellBoltzmannDistribution(atoms, temperature_K=T)
        dyn = Bussi(atoms, timestep=0.5*units.fs, temperature_K=T,
                    taut=100*units.fs)
        dyn.run(int(args.settle_ps*2000))
        grabbed = []
        interval = max(1, int(args.run_ps*2000/args.n_frames))
        for k in range(args.n_frames):
            dyn.run(interval)
            grabbed.append(atoms.copy())
            write(d / f"anchor_{k:02d}.xyz", atoms)
        Ts = [a.get_kinetic_energy()/(1.5*len(a)*units.kB) for a in grabbed]
        print(f"{tag}: {args.n_frames} frames, T = {np.mean(Ts):.0f} K (target {T})")
    print("anchor MD done")

if __name__ == "__main__":
    main()
