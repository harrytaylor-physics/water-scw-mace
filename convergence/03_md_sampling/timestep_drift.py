#!/usr/bin/env python3
"""NVE energy-drift vs timestep for a MACE potential.
Usage: timestep_drift.py model.model init.xyz --dts 0.25 0.5 1.0 --ps 5 --temp 1500"""
import argparse
import numpy as np
from ase import units
from ase.io import read
from ase.md.verlet import VelocityVerlet
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from mace.calculators import MACECalculator

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("model"); ap.add_argument("init")
    ap.add_argument("--dts", nargs="+", type=float, default=[0.25, 0.5, 1.0])
    ap.add_argument("--ps", type=float, default=5.0)
    ap.add_argument("--temp", type=float, default=1500.0)
    ap.add_argument("--device", default="cuda")
    args = ap.parse_args()
    calc = MACECalculator(model_paths=args.model, device=args.device,
                          default_dtype="float64")   # fp64: measure the integrator, not the dtype
    base = read(args.init)
    print(f"{'dt (fs)':>8} {'drift (meV/atom/ps)':>22}   (gate: |drift| <= 2)")
    for dt in args.dts:
        atoms = base.copy(); atoms.calc = calc
        MaxwellBoltzmannDistribution(atoms, temperature_K=args.temp)
        dyn = VelocityVerlet(atoms, timestep=dt * units.fs)
        es, ts = [], []
        def grab():
            es.append(atoms.get_total_energy()); ts.append(dyn.nsteps * dt / 1000)
        dyn.attach(grab, interval=20)
        dyn.run(int(args.ps * 1000 / dt))
        drift = np.polyfit(ts, np.array(es) * 1000 / len(atoms), 1)[0]
        print(f"{dt:8.2f} {drift:22.4f}   {'PASS' if abs(drift) <= 2 else 'FAIL'}")

if __name__ == "__main__":
    main()
