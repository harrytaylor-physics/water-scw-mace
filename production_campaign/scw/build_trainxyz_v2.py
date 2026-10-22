#!/usr/bin/env python3
"""All labels -> one train file, with sigma->0 energy extrapolation for smeared outputs.

E(sigma->0) = F - 0.5*E_entropic  (Gillan extrapolation; E_entropic = -T*S, as printed
by CP2K "Electronic entropic energy"). OT outputs have no entropic line -> energy used
as-is. Configs with |F|max > --fmax are dropped (the v3 detonator lesson)."""
import argparse, glob, os, re, numpy as np
from ase.io import read, write
from ase.calculators.singlepoint import SinglePointCalculator
HA=27.211386245988; BOHR=0.529177210903; FCONV=HA/BOHR; DEBYE=0.20819434
FLT=r"[-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?"
ABC=re.compile(r"ABC\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)")
EN=re.compile(r"ENERGY\| Total FORCE_EVAL \( QS \) energy.*?:\s*(-?[\d.]+)")
ENT=re.compile(rf"Electronic entropic energy:\s*({FLT})")
DIP=re.compile(rf"Dipole moment \[Debye\].*?X=\s*({FLT})\s+Y=\s*({FLT})\s+Z=\s*({FLT})",re.S)
def parse_forces(text,n):
    i=text.find("ATOMIC FORCES")
    if i<0: return None
    F=[]
    for ln in text[i:].splitlines()[1:]:
        if "SUM OF ATOMIC FORCES" in ln: break
        p=ln.split()
        if len(p)>=6 and p[0].isdigit(): F.append([float(p[3]),float(p[4]),float(p[5])])
    return np.array(F) if len(F)==n else None
def parse_one(d):
    out,cx,inp=(os.path.join(d,f) for f in("output.out","coord.xyz","input.inp"))
    if not(os.path.exists(out) and os.path.exists(cx)): return None,None
    t=open(out).read()
    if "PROGRAM ENDED" not in t: return None,None
    if "SCF run NOT converged" in t: return None,None
    em=list(EN.finditer(t))
    if not em: return None,None
    F_ha=float(em[-1].group(1))
    ents=list(ENT.finditer(t))
    corr=0.0; smeared=False
    if ents:
        smeared=True; corr=-0.5*float(ents[-1].group(1))   # E0 = F - 0.5*ent
    at=read(cx)
    m=ABC.search(open(inp).read())
    if m: at.set_cell([float(m.group(i)) for i in(1,2,3)]); at.set_pbc(True)
    Fa=parse_forces(t,len(at))
    if Fa is None: return None,None
    at.calc=SinglePointCalculator(at,energy=(F_ha+corr)*HA,forces=Fa*FCONV)
    dm=DIP.search(t)
    if dm: at.info["dipole"]=[float(dm.group(i))*DEBYE for i in(1,2,3)]
    at.info["config_type"]=os.path.basename(os.path.dirname(d))
    return at,smeared
ap=argparse.ArgumentParser()
ap.add_argument("--labels",nargs="+",default=["labels","labels_rcb","labels_aimd"])
ap.add_argument("--out",default="train_v4.xyz"); ap.add_argument("--fmax",type=float,default=50.0)
a=ap.parse_args()
good=[]; nsm=0; skip=0; hot=0; dropped=0
for root in a.labels:
    for d in sorted(glob.glob(os.path.join(root,"*","f*"))):
        at,sm=parse_one(d)
        if at is None: skip+=1; continue
        if np.abs(at.get_forces()).max()>a.fmax: dropped+=1; continue
        good.append(at); nsm+=bool(sm)
write(a.out,good,format="extxyz")
epa=[at.get_potential_energy()/len(at) for at in good]
print(f"{len(good)} configs -> {a.out}  ({nsm} smeared/corrected, {len(good)-nsm} OT)")
print(f"skipped {skip} incomplete, dropped {dropped} with |F|>{a.fmax}")
print(f"E/atom {min(epa):.3f}..{max(epa):.3f} eV")
