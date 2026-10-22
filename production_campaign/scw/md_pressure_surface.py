#!/usr/bin/env python3
"""
md_pressure_surface.py -- MACE-MD virial pressure vs AQUA EoS, as an independent check
that the force field reproduces the equation of state.

Re-evaluates the stored production trajectories ($DATA/scw_prod/*_05fs/prod.xyz) with the
water_scw_v4 force field (NO new MD), averages the virial pressure per cell, and compares
to AQUA P(T,rho). Writes md_pressure.csv + md_pressure_parity.png.

GPU batch job (MACE eval). Stress convention matches validate_v2.py:
    P = -(sxx+syy+szz)/3 * EV_A3_TO_PA   [Pa],  stress from get_stress(voigt=True) in eV/A^3.
"""
import os, glob, re, csv, numpy as np
from ase.io import iread
from aqua_eos import AquaEOS

DATA     = os.environ["DATA"]
MODEL    = "water_scw_v4.model"
AQUA_DAT = "aqua_eos_rhot_v1_0.dat"
STRIDE   = 50         # sample every 50th frame
MAXFRAMES= 200        # cap frames per cell (pressure converges fast; bounds runtime)
EV_A3_TO_PA = 1.602176634e11
NA=6.02214076e23; M=18.01528e-3   # kg/mol

def main():
    from mace.calculators import MACECalculator
    calc = MACECalculator(model_paths=MODEL, device="cuda", default_dtype="float64")
    eos  = AquaEOS(AQUA_DAT)

    rows=[]
    for d in sorted(glob.glob(os.path.join(DATA,"scw_prod","*_05fs"))):
        xyz=os.path.join(d,"prod.xyz")
        if not os.path.exists(xyz): continue
        m=re.search(r"T([0-9.]+)_rho([0-9.]+)", os.path.basename(d))
        if not m: continue
        T=float(m.group(1)); rho=float(m.group(2))
        Ps=[]; n=0
        for i,at in enumerate(iread(xyz)):
            if i % STRIDE: continue
            if at.cell.rank < 3 or at.get_volume() <= 0:        # rebuild cell from rho if absent
                L=(len(at)/3.0*M/(rho*NA))**(1/3)*1e10
                at.set_cell([L,L,L]); at.set_pbc(True)
            at.calc=calc
            s=at.get_stress(voigt=True)                          # eV/A^3
            Ps.append(-(s[0]+s[1]+s[2])/3.0*EV_A3_TO_PA)         # Pa
            n+=1
            if n>=MAXFRAMES: break
        if not Ps:
            print(f"  {os.path.basename(d):26s}  (no frames)"); continue
        Ps=np.array(Ps)
        P_md=Ps.mean()/1e9; P_err=Ps.std()/np.sqrt(len(Ps))/1e9   # GPa
        P_aqua=float(eos.P(T,rho))/1e9                            # GPa
        rows.append((os.path.basename(d),T,rho,P_md,P_err,P_aqua))
        print(f"  T={T:6.0f} rho={rho:8.1f}  P_MD={P_md:8.3f} +/- {P_err:5.3f} GPa   "
              f"P_AQUA={P_aqua:8.3f} GPa   ratio={P_md/P_aqua:5.2f}  (n={len(Ps)})", flush=True)

    with open("md_pressure.csv","w",newline="") as f:
        w=csv.writer(f); w.writerow(["cell","T_K","rho_kg_m3","P_MD_GPa","P_MD_err_GPa","P_AQUA_GPa"])
        w.writerows(rows)

    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    Pm=np.array([r[3] for r in rows]); Pe=np.array([r[4] for r in rows]); Pa=np.array([r[5] for r in rows])
    fig,ax=plt.subplots(figsize=(6,6))
    ax.errorbar(Pa,Pm,yerr=Pe,fmt="o",ms=5,capsize=2)
    lim=[min(Pa.min(),Pm.min())*0.6, max(Pa.max(),Pm.max())*1.4]
    ax.plot(lim,lim,"k--",lw=1,label="1:1")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel("AQUA EoS pressure / GPa"); ax.set_ylabel("MACE-MD virial pressure / GPa")
    ax.set_title("MACE force field vs AQUA equation of state"); ax.legend()
    plt.tight_layout(); plt.savefig("md_pressure_parity.png",dpi=150)
    print("wrote md_pressure.csv + md_pressure_parity.png")

if __name__=="__main__":
    main()
