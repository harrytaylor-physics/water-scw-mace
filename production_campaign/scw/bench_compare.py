# bench_compare.py  --  python bench_compare.py <tag> <T_K> <rho_g_cc>
import sys, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

tag=sys.argv[1]; T=float(sys.argv[2]); rho=float(sys.argv[3])      # rho in g/cc
v4=np.load(f"bench/{tag}_v4/bench.npz"); fo=np.load(f"bench/{tag}_foundation/bench.npz")

# NOTE: bench.npz now stores the TOTAL pressure (virial + ideal-gas), since bench_md.py /
# bench_pressure.py were fixed to add the kinetic term. So we read P_GPa directly here.
# Re-run the benchmarks with the fixed scripts before trusting these numbers -- any old
# bench.npz from before the fix holds virial-only (negative) pressures.

# --- reference pressure: IAPWS-95 (<=1 GPa anchor) and AQUA (carries the dense data) ---
P_iapws=None
try:
    from iapws import IAPWS95
    P_iapws=IAPWS95(T=T, rho=rho*1000).P*1e-3            # MPa -> GPa
except Exception as e:
    print("IAPWS-95 unavailable (pip install iapws on login node):", e)

P_aqua=None
try:
    A=np.loadtxt("/home/shug7909/scw/aqua_eos_rhot_v1_0.dat", skiprows=24, usecols=(0,1,2))  # rho[kg/m3], T[K], P[Pa]
    rho_kg=rho*1000.0
    d=(np.log10(A[:,0])-np.log10(rho_kg))**2 + (np.log10(A[:,1])-np.log10(T))**2  # nearest node in log space
    P_aqua=A[d.argmin(),2]*1e-9                                            # Pa -> GPa
    print(f"  AQUA node: rho={A[d.argmin(),0]:.1f} kg/m3  T={A[d.argmin(),1]:.1f} K")
except Exception as e:
    print("AQUA lookup skipped:", e)

# --- ambient experimental peak positions, shown only as a GUIDE (not a target at 673/1000 K) ---
AMB_RDF1, AMB_RDF2 = 2.80, 4.5                            # ambient O-O 1st/2nd shell (A)
EXP_BANDS={"H-bond":180,"libration":600,"bend":1650,"stretch(liq)":3400,"stretch(gas)":3657}

fig,ax=plt.subplots(1,2,figsize=(13,5))

# RDF panel
ax[0].plot(v4["rc"],v4["gOO"],lw=1.6,label="v4 (fine-tuned)")
ax[0].plot(fo["rc"],fo["gOO"],lw=1.6,label="foundation")
ax[0].axvline(AMB_RDF1,ls=":",c="grey",alpha=.6)
ax[0].axvline(AMB_RDF2,ls=":",c="grey",alpha=.6)
ax[0].text(AMB_RDF2,ax[0].get_ylim()[1]*0.92,"ambient 2nd shell\n(expect collapsed here)",
           fontsize=7,color="grey",ha="center")
ax[0].set_xlabel("r (A)"); ax[0].set_ylabel("g_OO(r)")
ax[0].set_title(f"{tag}  {T:.0f} K, {rho:.2f} g/cc : O-O RDF\n(dotted = ambient peak positions, guide only)")
ax[0].legend()

# VDOS panel
ax[1].plot(v4["wn"],v4["dos"]/v4["dos"].max(),lw=1.4,label="v4")
ax[1].plot(fo["wn"],fo["dos"]/fo["dos"].max(),lw=1.4,label="foundation")
for n,x in EXP_BANDS.items():
    ax[1].axvline(x,ls=":",c="grey",alpha=.5); ax[1].text(x,1.02,n,rotation=90,fontsize=7)
ax[1].set_xlim(0,4200); ax[1].set_xlabel("cm-1"); ax[1].set_ylabel("VDOS (norm)")
ax[1].set_title("vibrational density of states\n(stretch should blue-shift toward gas with T)")
ax[1].legend()
plt.tight_layout(); plt.savefig(f"bench/{tag}_compare.png",dpi=140)

# --- pressure verdict (P_GPa is already the total pressure) ---
v4_P, v4_std = float(v4['P_GPa']), float(v4['P_std'])
fo_P, fo_std = float(fo['P_GPa']), float(fo['P_std'])
print(f"\n{tag}  T={T:.0f} K  rho={rho:.2f} g/cc  pressure (GPa, total = virial + ideal-gas):")
print(f"  v4         = {v4_P:+.3f} +/- {v4_std:.3f}")
print(f"  foundation = {fo_P:+.3f} +/- {fo_std:.3f}")
if P_iapws is not None: print(f"  IAPWS-95   = {P_iapws:+.3f}   (experimental standard; trust <~1 GPa)")
if P_aqua  is not None: print(f"  AQUA       = {P_aqua:+.3f}   (carries QMD data at the dense end)")
if v4_P < 0:
    print("  !! v4 pressure is negative at this density -- likely an OLD virial-only bench.npz;")
    print("     re-run bench_md.py / bench_pressure.py (now include the kinetic term).")
print(f"  -> wrote bench/{tag}_compare.png")
