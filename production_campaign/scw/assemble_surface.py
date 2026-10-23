"""assemble_surface.py -- gather every cell's opacity_v3.npz into the opacity surface."""
import os, glob, csv, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

files = sorted(glob.glob(os.path.join(os.environ["DATA"], "scw_prod", "*_05fs", "opacity_v3.npz")))
if not files:
    raise SystemExit("no opacity_v3.npz found in *_05fs cells")

cells = []
for f in files:
    d = np.load(f)
    cells.append(dict(cell=os.path.basename(os.path.dirname(f)),
                      T=float(d["T"]), rho=float(d["rho_kg_m3"]),
                      P=float(d["planck"]), R=float(d["rosseland"]),
                      wn=d["wn"], kappa=d["kappa"]))
cells.sort(key=lambda c: (c["T"], c["rho"]))

with open("opacity_surface_means.csv", "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["cell","T_K","rho_kg_m3","Planck_cm2_g","Rosseland_cm2_g"])
    for c in cells: w.writerow([c["cell"], f"{c['T']:.0f}", f"{c['rho']:.3f}", f"{c['P']:.6e}", f"{c['R']:.6e}"])
print(f"wrote opacity_surface_means.csv  ({len(cells)} cells)")
for c in cells:
    print(f"  {c['cell']:24s} T={c['T']:.0f} rho={c['rho']:8.2f}  P={c['P']:.3e}  R={c['R']:.3e}")

plt.figure(figsize=(10,6)); cmap=plt.cm.viridis; Ts=sorted({c["T"] for c in cells})
for c in cells:
    col=cmap((c["T"]-min(Ts))/max(max(Ts)-min(Ts),1))
    plt.semilogy(c["wn"], np.maximum(c["kappa"],1e-2), lw=1.1, color=col, label=f"{c['T']:.0f}K {c['rho']:.0f}")
plt.xlim(0,6000); plt.xlabel("wavenumber / cm$^{-1}$"); plt.ylabel("$\\kappa$ / cm$^2$ g$^{-1}$")
plt.legend(fontsize=6, ncol=2); plt.title("MACE-MD continuum opacity (revPBE0-D3, v3 dipole) -- condensed cells")
plt.tight_layout(); plt.savefig("opacity_surface_kappa.png", dpi=140); print("wrote opacity_surface_kappa.png")

fig, ax = plt.subplots(1,2, figsize=(12,5))
for c in cells:
    col=cmap((c["T"]-min(Ts))/max(max(Ts)-min(Ts),1))
    ax[0].scatter(c["rho"], c["P"], color=col, s=40); ax[1].scatter(c["rho"], c["R"], color=col, s=40)
for a,ttl in zip(ax,["Planck mean","Rosseland mean"]):
    a.set_xlabel("density / kg m$^{-3}$"); a.set_ylabel("$\\kappa$ / cm$^2$ g$^{-1}$")
    a.set_title(ttl); a.set_xscale("log"); a.set_yscale("log")
sm=plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(min(Ts),max(Ts))); fig.colorbar(sm, ax=ax, label="T / K")
plt.savefig("opacity_surface_means.png", dpi=140); print("wrote opacity_surface_means.png")

wn0=cells[0]["wn"]; kap=np.vstack([np.interp(wn0,c["wn"],c["kappa"]) for c in cells])
np.savez("opacity_surface.npz", wn=wn0, kappa=kap,
         T=np.array([c["T"] for c in cells]), rho=np.array([c["rho"] for c in cells]),
         planck=np.array([c["P"] for c in cells]), rosseland=np.array([c["R"] for c in cells]),
         cells=np.array([c["cell"] for c in cells]))
print("wrote opacity_surface.npz")
