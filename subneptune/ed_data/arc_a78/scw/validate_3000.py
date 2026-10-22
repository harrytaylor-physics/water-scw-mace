"""
validate_3000.py  --  compare MACE-MD v3 opacity against the AIMD reference at 3000 K / 0.1 g/cc.

Usage:
    python validate_3000.py [aimd_database_compare.csv]

Always prints the v3 Planck/Rosseland means vs the AIMD targets (94,332 / 29,860 cm2/g).
If the AIMD spectral CSV path is given (h2o32_3000K_rho01_final_database_compare.csv),
also writes an overlay of sigma(omega): v3 vs AIMD vs POKAZATEL vs HITEMP.
"""
import sys, os, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

CELL = os.path.join(os.environ["DATA"], "scw_prod", "T3000_rho108.796")
AIMD_PLANCK, AIMD_ROSS = 9.4332e4, 2.9860e4   # from the means CSV

v3 = np.load(os.path.join(CELL, "opacity_v3.npz"))
P, R = float(v3["planck"]), float(v3["rosseland"])
print("="*60)
print(f"  v3   Planck = {P:.4e}   Rosseland = {R:.4e}  cm2/g")
print(f"  AIMD Planck = {AIMD_PLANCK:.4e}   Rosseland = {AIMD_ROSS:.4e}  cm2/g")
print(f"  ratio v3/AIMD:  Planck {P/AIMD_PLANCK:.2f}x   Rosseland {R/AIMD_ROSS:.2f}x")
print("="*60)

if len(sys.argv) > 1 and os.path.exists(sys.argv[1]):
    import pandas as pd
    df = pd.read_csv(sys.argv[1])
    wn_a = df["wn_cm-1"].to_numpy()
    sig_a = df["aimd_sigma_cm2_per_molecule"].to_numpy()
    wn_v, sig_v = v3["wn"], v3["sigma"]
    plt.figure(figsize=(10, 6))
    plt.semilogy(wn_a, np.maximum(sig_a, 1e-30), lw=2.4, label="AIMD Berry (revPBE-D3)")
    plt.semilogy(wn_v, np.maximum(sig_v, 1e-30), lw=2.4, label="MACE-MD v3 dipole")
    for col, lab in [("exocross_doppler_sigma_cm2_per_molecule", "POKAZATEL Doppler"),
                     ("hitemp_sigma_cm2_per_molecule", "HITEMP")]:
        if col in df.columns:
            plt.semilogy(wn_a, np.maximum(df[col].to_numpy(), 1e-30), lw=1.2, alpha=0.7, label=lab)
    plt.xlim(0, 6000); plt.ylim(1e-24, None)
    plt.xlabel("wavenumber / cm$^{-1}$"); plt.ylabel("$\\sigma$ / cm$^2$ molecule$^{-1}$")
    plt.legend(); plt.title("3000 K, 0.1 g/cc:  MACE-MD v3 vs AIMD reference")
    plt.tight_layout()
    outp = os.path.join(CELL, "validate_3000.png")
    plt.savefig(outp, dpi=140); print("wrote", outp)
else:
    print("(no spectral CSV given -> means-only check; pass the AIMD database_compare.csv path for the band overlay)")
