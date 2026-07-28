#!/usr/bin/env python3
"""Overlay water_scw_v4 on the +-D3 dimer curves: did MACE learn the tail?"""
import argparse, pathlib
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from analyze_dispersion import parse, first_out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("scan"); ap.add_argument("--model", required=True)
    ap.add_argument("--device", default="cuda")
    args = ap.parse_args()
    from ase import Atoms
    from mace.calculators import MACECalculator
    calc = MACECalculator(model_paths=args.model, device=args.device,
                          default_dtype="float64")
    root = pathlib.Path(args.scan)
    rs, em, e1, e0 = [], [], [], []
    for d in sorted(root.glob("dimer/r*_d3")):
        E1, _, c1 = parse(first_out(d))
        E0, _, c0 = parse(first_out(str(d).replace("_d3", "_nod3")))
        if not (c1 and c0): continue
        from analyze_dispersion import read_coords
        lines = read_coords(d)
        sym = [l.split()[0] for l in lines]
        pos = np.array([l.split()[1:4] for l in lines], float)
        a = Atoms(sym, positions=pos, cell=[15,15,15], pbc=True); a.calc = calc
        rs.append(float(d.name[1:-3])); em.append(a.get_potential_energy())
        e1.append(E1); e0.append(E0)
    rs = np.array(rs)
    em, e1, e0 = (np.array(x) - np.array(x)[-1] for x in (em, e1, e0))
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(rs, e1*1000, "k-", label="revPBE0-D3 (training level)")
    ax.plot(rs, e0*1000, "k--", label="revPBE0 (no D3)")
    ax.plot(rs, em*1000, "k.", ms=7, label="water_scw_v4 (MACE)")
    ax.set_xlabel("O-O separation (A)"); ax.set_ylabel("binding energy (meV)")
    ax.legend(frameon=False)
    fig.tight_layout(); fig.savefig(root / "mace_dispersion.png", dpi=170)
    tail = rs > 4.0
    r1 = np.sqrt(np.mean((em[tail]-e1[tail])**2))*1000
    r0 = np.sqrt(np.mean((em[tail]-e0[tail])**2))*1000
    print(f"tail (r>4 A) RMS: MACE vs +D3 {r1:.1f} meV | vs noD3 {r0:.1f} meV")
    print("verdict:", "MACE learned the dispersion tail"
          if r1 < r0 else "tail closer to no-D3 - NOT fully learned (log it)")

if __name__ == "__main__":
    main()
