import os, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from opacity_pipeline import spectrum_from_dipoles, c_cm_s, h, kB, _trapz

df = pd.read_csv("h2o32_3000K_rho01_final_dipoles.csv")
Ma = df[['Mx_D','My_D','Mz_D']].to_numpy()
ra = spectrum_from_dipoles(Ma, T_K=2998.07, dt_fs=0.5, n_molecules=32, density_g_cm3=0.1,
                           units="debye", smooth_sigma_cm=20.0)
M = np.load(os.path.join(os.environ['DATA'],'scw_prod','T3000_rho108.796_05fs','dipoles.npy'))
rv = spectrum_from_dipoles(M, T_K=3000, dt_fs=0.5, n_molecules=64, density_g_cm3=0.108796,
                           units="eAng", smooth_sigma_cm=20.0)

wn = ra['wn_cm-1']; sa = ra['sigma_cm2_per_molecule']
sv = np.interp(wn, rv['wn_cm-1'], rv['sigma_cm2_per_molecule'])
nu = wn*c_cm_s; T=3000.0; B = (nu**3)/np.maximum(np.exp(h*nu/(kB*T))-1, 1e-300)
print(f"{'band':16s}{'AIMD sig-int':>14s}{'v3 sig-int':>14s}{'Pwt ratio':>11s}")
for lo,hi,nm in [(0,900,'librational'),(900,2000,'bend'),(2000,4200,'stretch'),(4200,6000,'hi-freq tail')]:
    m=(wn>=lo)&(wn<hi)
    ia=_trapz(sa[m],wn[m]); iv=_trapz(sv[m],wn[m])
    pa=_trapz(sa[m]*B[m],nu[m]); pv=_trapz(sv[m]*B[m],nu[m])
    print(f"{nm:16s}{ia:14.3e}{iv:14.3e}{pv/max(pa,1e-300):11.2f}")
plt.figure(figsize=(10,6))
plt.semilogy(wn, np.maximum(sa,1e-30), lw=2.2, label="AIMD Berry")
plt.semilogy(wn, np.maximum(sv,1e-30), lw=2.2, label="v3 MACE 0.5fs")
plt.xlim(0,6000); plt.ylim(1e-23,None); plt.xlabel("cm-1"); plt.ylabel("sigma cm2/molecule")
plt.legend(); plt.title("3000K/0.1: shape comparison (nb=1, 0.5fs)"); plt.tight_layout()
plt.savefig(os.path.join(os.environ['DATA'],'scw_prod','T3000_rho108.796_05fs','band_compare.png'), dpi=140)
print("wrote band_compare.png")
