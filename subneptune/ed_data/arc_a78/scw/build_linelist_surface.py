"""
build_linelist_surface.py -- HITEMP & ExoMol Rosseland surfaces at the MACE state points,
both lines-only and lines+MT_CKD, then the breakdown-ratio figures. RUN AS A BATCH JOB.

Per state point (T, rho):
    lines    : HITEMP / ExoMol discrete lines at fixed reference broadening (T-dependent)
    +MT_CKD  : add the MT_CKD self continuum evaluated at the cell's actual (T, rho)
Outputs: ratio_*.png (lines only) and ratio_*_mtckd.png (lines + empirical continuum),
plus linelist_surface.npz with all Rosseland arrays.
"""
import os, numpy as np
from aqua_eos import AquaEOS
from opacity_means import rosseland_mean
from make_ratio_figures import make_all
from hitemp_opacity import HitempCalculator
from exomolop_opacity import exomolop_kappa
from mtckd_continuum import MTCKD

AQUA_DAT    = "aqua_eos_rhot_v1_0.dat"
EXOMOLOP_H5 = "POKAZATEL.TauREx.h5"
MTCKD_NPZ   = "mtckd_ref.npz"
MACE_NPZ    = "/home/shug7909/scw/opacity_surface.npz"
NUMIN, NUMAX, WSTEP = 10.0, 11000.0, 0.01
EXO_CLIP = 200.0

def _get(d,*names):
    for n in names:
        if n in d.files: return np.asarray(d[n],float)
    raise KeyError(f"none of {names} in {d.files}")

def main():
    eos=AquaEOS(AQUA_DAT); mt=MTCKD(MTCKD_NPZ)
    H=HitempCalculator(NUMIN,NUMAX,WSTEP)   # load HITEMP databank once
    d=np.load(MACE_NPZ,allow_pickle=True); print("MACE npz keys:",d.files)
    T=_get(d,"T","Tgrid","temperature"); rho=_get(d,"rho","rho_kg_m3","density")
    R_mace=_get(d,"rosseland","Rosseland","rosseland_cm2_g","R_rosseland")

    R_hL=np.empty_like(R_mace); R_xL=np.empty_like(R_mace)   # lines only
    R_hM=np.empty_like(R_mace); R_xM=np.empty_like(R_mace)   # lines + MT_CKD
    for i,(Ti,ri) in enumerate(zip(T,rho)):
        P_atm=float(eos.P(Ti,ri))/101325.0
        wn_h,kh=H.kappa(Ti)
        kc_h=mt.kappa_self(wn_h,Ti,ri)
        R_hL[i]=rosseland_mean(wn_h,kh,Ti); R_hM[i]=rosseland_mean(wn_h,kh+kc_h,Ti)
        wn_x,kx=exomolop_kappa(EXOMOLOP_H5,Ti); m=wn_x>=EXO_CLIP
        kc_x=mt.kappa_self(wn_x[m],Ti,ri)
        R_xL[i]=rosseland_mean(wn_x[m],kx[m],Ti); R_xM[i]=rosseland_mean(wn_x[m],kx[m]+kc_x,Ti)
        print(f"  T={Ti:6.0f} rho={ri:8.1f} P={P_atm:9.1f}atm | "
              f"MACE={R_mace[i]:.2e}  HIT(L/M)={R_hL[i]:.2e}/{R_hM[i]:.2e}  "
              f"EXO(L/M)={R_xL[i]:.2e}/{R_xM[i]:.2e}",flush=True)

    np.savez("linelist_surface.npz",T=T,rho=rho,R_mace=R_mace,
             R_hitemp=R_hL,R_exomol=R_xL,R_hitemp_mtckd=R_hM,R_exomol_mtckd=R_xM)
    make_all(T,rho,R_mace,R_hL,R_xL,outdir=".",suffix="")        # lines only
    make_all(T,rho,R_mace,R_hM,R_xM,outdir=".",suffix="_mtckd")  # lines + continuum
    print("done -> linelist_surface.npz, ratio_*.png (lines), ratio_*_mtckd.png (lines+MT_CKD)")

if __name__=="__main__":
    main()
