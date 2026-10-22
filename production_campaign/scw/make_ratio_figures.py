"""
make_ratio_figures.py -- three (T, rho) ratio maps with the line-list breakdown contour.

ratio maps (colour = log10 ratio, diverging about ratio=1; bold contour = ratio 1):
   1) MACE / ExoMol     -> where the first-principles continuum overtakes POKAZATEL(+MT_CKD)
   2) MACE / HITEMP     -> ... overtakes HITEMP(+MT_CKD)
   3) ExoMol / HITEMP   -> where the two line lists disagree (completeness floor)

The ratio=1 contour on (1) and (2) is the breakdown line: inside it (red) the
water-vapour continuum dominates and the line list can no longer represent the opacity.

Call:
   make_all(T, rho, R_mace, R_hitemp, R_exomol, outdir=".")
all inputs are 1-D arrays over the SAME state points (Rosseland means, cm^2/g).
"""
import os, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, matplotlib.tri as mtri
from matplotlib.colors import TwoSlopeNorm

def _panel(ax, x, y, ratio, title):
    z=np.log10(ratio); tri=mtri.Triangulation(x,y)
    vmax=max(0.3,np.nanmax(np.abs(z))); norm=TwoSlopeNorm(vcenter=0,vmin=-vmax,vmax=vmax)
    cf=ax.tricontourf(tri,z,levels=np.linspace(-vmax,vmax,21),cmap="RdBu_r",norm=norm,extend="both")
    # breakdown line (ratio=1) + guides at x2 and x10
    for lv,style in [(0.0,dict(colors="k",linewidths=2.2)),
                     (np.log10(2),dict(colors="k",linewidths=0.8,linestyles="--")),
                     (1.0,dict(colors="k",linewidths=0.8,linestyles=":"))]:
        try: ax.tricontour(tri,z,levels=[lv],**style)
        except Exception: pass
    ax.scatter(x,y,c="k",s=8,zorder=4)
    ax.set_title(title); ax.set_xlabel("density / kg m$^{-3}$"); ax.set_ylabel("temperature / K")
    xt=[2,2.3,2.6,3,3.3]; ax.set_xticks(xt); ax.set_xticklabels([f"{10**t:.0f}" for t in xt])
    return cf

def make_all(T, rho, R_mace, R_hitemp, R_exomol, outdir=".", suffix=""):
    T=np.asarray(T,float); x=np.log10(np.asarray(rho,float))
    specs=[("MACE_over_ExoMol", R_mace/R_exomol, "MACE / ExoMol  (POKAZATEL+MT_CKD)"),
           ("MACE_over_HITEMP", R_mace/R_hitemp, "MACE / HITEMP  (+MT_CKD)"),
           ("ExoMol_over_HITEMP", R_exomol/R_hitemp, "ExoMol / HITEMP")]
    for name,ratio,title in specs:
        fig,ax=plt.subplots(figsize=(7.2,5.6))
        cf=_panel(ax,x,T,np.asarray(ratio,float),title)
        cb=fig.colorbar(cf,ax=ax); cb.set_label("log$_{10}$ ratio   (bold = breakdown, ratio=1)")
        plt.tight_layout(); plt.savefig(os.path.join(outdir,f"ratio_{name}{suffix}.png"),dpi=150); plt.close()
        print("wrote", f"ratio_{name}{suffix}.png")
