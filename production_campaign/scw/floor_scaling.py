import os, numpy as np
from opacity_pipeline import spectrum_from_dipoles
dip = np.load(os.path.join(os.environ["DATA"], "scw_prod", "T3000_rho108.796", "dipoles.npy"))
print("coarsening the existing 1 fs data (noise floor should grow with dt):")
for stride in (1, 2, 3):
    M = dip[::stride]; dtf = 1.0 * stride
    r = spectrum_from_dipoles(M, T_K=3000, dt_fs=dtf, n_molecules=64,
                              density_g_cm3=0.108796, units="eAng", smooth_sigma_cm=20.0, n_blocks=1)
    fl_hi = np.median(r["sigma_cm2_per_molecule"][(r["wn_cm-1"] > 4500) & (r["wn_cm-1"] < 6000)])
    print(f"  dt={dtf:.1f}fs n={len(M):5d}  Planck {r['planck_mean_cm2_g']/9.4332e4:.2f}x"
          f"  Rosseland {r['rosseland_mean_cm2_g']/2.9860e4:.2f}x  hi-freq floor {fl_hi:.2e}")
