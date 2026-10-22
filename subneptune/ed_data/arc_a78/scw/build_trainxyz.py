#!/usr/bin/env python3
"""Phase 2 assembly: CP2K revPBE0-D3 single-point outputs -> MACE train.xyz.

Walks labels/*/f*/output.out, parses the energy, atomic forces and Berry-phase
dipole, pairs them with the geometry (coord.xyz) and cell (ABC in input.inp),
converts to MACE units, and writes one extended-XYZ training file.

Units: CP2K gives energy in Hartree and forces in Hartree/Bohr; we write eV and
eV/Ang. Dipole is printed in Debye; we write e*Ang. (The Berry-phase dipole is
defined modulo the cell quantum -- fine to store, handle the gauge at ACF time.)

Safe to run while the array is still going: incomplete outputs are skipped, so
re-run when squeue is empty to pick up the rest.

    python build_trainxyz.py                 # -> train.xyz
Then fine-tune:  mace_run_train --foundation_model .../mace-mpa-0-medium.model \
                   --train_file train.xyz --energy_key energy --forces_key forces ...
"""
import argparse, glob, os, re, numpy as np
from ase.io import read, write
from ase.calculators.singlepoint import SinglePointCalculator

HA   = 27.211386245988          # eV / Hartree
BOHR = 0.529177210903           # Ang / Bohr
FCONV = HA / BOHR               # Hartree/Bohr -> eV/Ang
DEBYE = 0.20819434              # e*Ang / Debye

FLT     = r"[-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?"
ABC_PAT = re.compile(r"ABC\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)")
EN_PAT  = re.compile(r"ENERGY\| Total FORCE_EVAL \( QS \) energy.*?:\s*(-?[\d.]+)")
DIP_PAT = re.compile(rf"Dipole moment \[Debye\].*?X=\s*({FLT})\s+Y=\s*({FLT})\s+Z=\s*({FLT})", re.S)

def parse_forces(text, natoms):
    idx = text.find("ATOMIC FORCES")
    if idx < 0:
        return None
    forces = []
    for ln in text[idx:].splitlines()[1:]:
        if "SUM OF ATOMIC FORCES" in ln:
            break
        p = ln.split()
        if len(p) >= 6 and p[0].isdigit():
            forces.append([float(p[3]), float(p[4]), float(p[5])])
    return np.array(forces) if len(forces) == natoms else None

def parse_one(outdir):
    out = os.path.join(outdir, "output.out")
    cx  = os.path.join(outdir, "coord.xyz")
    inp = os.path.join(outdir, "input.inp")
    if not (os.path.exists(out) and os.path.exists(cx)):
        return None
    text = open(out).read()
    if "PROGRAM ENDED" not in text:                 # still running / crashed
        return None
    em = EN_PAT.search(text)
    if not em:
        return None
    atoms = read(cx)
    abc = ABC_PAT.search(open(inp).read())
    if abc:
        atoms.set_cell([float(abc.group(i)) for i in (1, 2, 3)])
        atoms.set_pbc(True)
    farr = parse_forces(text, len(atoms))
    if farr is None:
        return None
    atoms.calc = SinglePointCalculator(atoms,
                                       energy=float(em.group(1)) * HA,
                                       forces=farr * FCONV)
    dm = DIP_PAT.search(text)
    if dm:
        atoms.info["dipole"] = [float(dm.group(i)) * DEBYE for i in (1, 2, 3)]
    atoms.info["config_type"] = os.path.basename(os.path.dirname(outdir))
    return atoms

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--labels", default="labels")
    ap.add_argument("--out", default="train.xyz")
    a = ap.parse_args()

    dirs = sorted(glob.glob(os.path.join(a.labels, "*", "f*")))
    good, skipped, nodip = [], 0, 0
    for d in dirs:
        at = parse_one(d)
        if at is None:
            skipped += 1
            continue
        if "dipole" not in at.info:
            nodip += 1
        good.append(at)
    if not good:
        raise SystemExit("no completed single points parsed yet -- is the array still running?")
    write(a.out, good, format="extxyz")
    epa = [at.get_potential_energy() / len(at) for at in good]
    print(f"parsed {len(good)} configs -> {a.out}")
    print(f"  skipped {skipped} (incomplete/failed), {nodip} had no dipole")
    print(f"  energy/atom: {min(epa):.3f} .. {max(epa):.3f} eV   ({len(good[0])} atoms/frame)")

if __name__ == "__main__":
    main()
