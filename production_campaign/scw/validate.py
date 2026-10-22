"""validate.py -- check each data module independently on ONE state point. Run: python validate.py
Skips (does not crash on) anything not yet installed/downloaded, and reports status per module."""
import numpy as np
T0, RHO0 = 3000.0, 108.8

print("== AQUA ==")
try:
    from aqua_eos import AquaEOS
    eos=AquaEOS("aqua_eos_rhot_v1_0.dat")   # or AQUA/Tables/aqua_eos_rhot_v1_0.dat
    P_atm=float(eos.P(T0,RHO0))/101325.0
    print(f"  OK   P({T0:.0f}K,{RHO0} kg/m3) = {P_atm:.1f} atm")
except Exception as e:
    P_atm=None; print("  SKIP/FAIL:", type(e).__name__, e)

print("== HITEMP (hapi) ==")
try:
    from hitemp_opacity import register_hitemp, hitemp_kappa
    name=register_hitemp("H2O_HITEMP.par")
    wn,kh=hitemp_kappa(name,T0)
    print(f"  OK   HITEMP kappa max = {kh.max():.3e} cm2/g  on {wn[0]:.0f}-{wn[-1]:.0f} cm-1")
except Exception as e:
    wn=None; print("  SKIP/FAIL:", type(e).__name__, e)

print("== ExoMolOP ==")
try:
    from exomolop_opacity import exomolop_kappa
    wx,kx=exomolop_kappa("POKAZATEL.TauREx.h5",T0)
    print(f"  OK   ExoMol kappa max = {kx.max():.3e} cm2/g")
except Exception as e:
    print("  SKIP/FAIL:", type(e).__name__, e)

print("== MT_CKD (optional, second pass) ==")
print("  AER repo is Fortran + netCDF; no pip package. Wiring deferred until lines-only works.")
