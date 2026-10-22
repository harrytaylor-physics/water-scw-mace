import sys, os, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from opacity_pipeline import spectrum_from_dipoles
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

AIMD_CSV = sys.argv[1]
NB = int(sys.argv[2]) if len(sys.argv) > 2 else 8
CELL = os.path.join(os.environ["DATA"], "scw_prod", "T3000_rho108.796")

df = pd.read_csv(AIMD_CSV)
Ma = np.column_stack([df["Mx_D"].to_numpy(), df["My_D"].to_numpy(), df["Mz_D"].to_numpy()])
ra = spectrum_from_dipoles(Ma, T_K=2998.07, dt_fs=0.5, n_molecules=32,
                           density_g_cm3=0.1, units="debye", smooth_sigma_cm=20.0)
print("=" * 64)
print("SELF-CHECK (your AIMD dipoles through my pipeline):")
print(f"   Planck    {ra['planck_mean_cm2_g']:.4e}  target 9.4332e+04  ratio {ra['planck_mean_cm2_g']/9.4332e4:.3f}")
print(f"   Rosseland {ra['rosseland_mean_cm2_g']:.4e}  target 2.9860e+04  ratio {ra['rosseland_mean_cm2_g']/2.9860e4:.3f}")
print("   (both ~1.0 => pipeline faithful; the v3 gap is physical)")
print("=" * 64)

dip = np.load(os.path.join(CELL, "dipoles.npy"))
print("v3 MACE-MD block-averaging sweep (1 fs dipoles):")
for nb in (1, 4, 8, 16):
    rv = spectrum_from_dipoles(dip, T_K=3000, dt_fs=1.0, n_molecules=64,
                               density_g_cm3=0.108796, units="eAng", smooth_sigma_cm=20.0, n_blocks=nb)
    print(f"   nb={nb:2d}  Planck {rv['planck_mean_cm2_g']:.3e} ({rv['planck_mean_cm2_g']/9.4332e4:.2f}x)"
          f"   Rosseland {rv['rosseland_mean_cm2_g']:.3e} ({rv['rosseland_mean_cm2_g']/2.9860e4:.2f}x)")
print("=" * 64)

rv = spectrum_from_dipoles(dip, T_K=3000, dt_fs=1.0, n_molecules=64,
                           density_g_cm3=0.108796, units="eAng", smooth_sigma_cm=20.0, n_blocks=NB)
plt.figure(figsize=(10, 6))
plt.semilogy(ra["wn_cm-1"], np.maximum(ra["sigma_cm2_per_molecule"], 1e-30), lw=2.4, label="AIMD Berry")
plt.semilogy(rv["wn_cm-1"], np.maximum(rv["sigma_cm2_per_molecule"], 1e-30), lw=2.4, label=f"MACE-MD v3 (nb={NB})")
plt.xlim(0, 6000); plt.ylim(1e-23, None)
plt.xlabel("wavenumber / cm$^{-1}$"); plt.ylabel("sigma / cm$^2$ molecule$^{-1}$")
plt.legend(); plt.title("3000 K, 0.1 g/cc: v3 vs AIMD (identical pipeline)")
plt.tight_layout(); plt.savefig(os.path.join(CELL, "compare_3000.png"), dpi=140)
print("wrote", os.path.join(CELL, "compare_3000.png"))
