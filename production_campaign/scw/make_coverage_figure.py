# make_coverage_figure.py -- coverage of the ACTUAL training set vs production MD visits.
#   python make_coverage_figure.py --train train_v4.xyz
# Standalone replacement for the coverage panel in make_missing_figures.py: instead of
# reconstructing the selection from the pools, it reads every configuration in the real
# training file (bulk + RCB + hot/AIMD, all ~1300), computes the same two descriptors,
# and overlays them on the production trajectories. Writes figs/coverage_training_vs_production.png
import sys, os, re, numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from ase.io import read

ANCHORS=[("scw673",673.0,0.66),("hot1000",1000.0,0.80)]
RCUT=3.4
os.makedirs("figs",exist_ok=True)

def descriptors_O(O,L):
    d=O[:,None,:]-O[None,:,:]; d-=L*np.round(d/L)
    r=np.sqrt((d*d).sum(-1)); np.fill_diagonal(r,np.inf)
    return float(r.min(1).mean()), float((r<RCUT).sum(1).mean())

def temp_of(at):
    for k in ("config_type","regime","source"):
        if k in at.info:
            m=re.search(r'[Tt](\d{3,4})',str(at.info[k]))
            if m: return float(m.group(1))
    return np.nan

i=sys.argv.index("--train"); train_path=sys.argv[i+1]
ats=read(train_path,index=":")
pts=[]; skipped=0
for at in ats:
    c=np.array(at.cell)
    if not np.any(c): skipped+=1; continue
    L=float(np.mean(np.diag(c)))
    isO=np.array(at.get_chemical_symbols())=="O"
    if isO.sum()<2: skipped+=1; continue
    d1,d2=descriptors_O(at.get_positions()[isO],L)
    pts.append((d1,d2,temp_of(at)))
pts=np.array(pts)
print(f"training set: {len(pts)} configs plotted ({skipped} skipped: no cell / <2 O)")
print(f"  with T tag: {np.isfinite(pts[:,2]).sum()}, untagged: {(~np.isfinite(pts[:,2])).sum()}")

fig,ax=plt.subplots(1,len(ANCHORS),figsize=(6.8*len(ANCHORS),5.6),squeeze=False); ax=ax[0]
for j,(tag,T,rho) in enumerate(ANCHORS):
    f=f"traj/{tag}/v4_n64.frames.npz"
    if os.path.exists(f):
        D=np.load(f); L=float(D["L"])
        prod=np.array([descriptors_O(O,L) for O in D["framesO"]])
        ax[j].scatter(prod[:,0],prod[:,1],s=8,c="0.75",zorder=1,
                      label=f"production MD visits ({len(prod)} frames)")
        # nearest-training-point distance for each production frame (normalised axes)
        tr=pts[:,:2]; span=prod.max(0)-prod.min(0)+1e-12
        dd=np.sqrt((((prod[:,None,:]-tr[None,:,:])/span)**2).sum(-1)).min(1)
        print(f"  {tag}: nearest-training-point distance (axis-normalised) "
              f"median={np.median(dd):.3f}, 95th pct={np.percentile(dd,95):.3f}")
    else:
        ax[j].text(.5,.5,"no production trajectory",ha="center",transform=ax[j].transAxes)
    tagged=np.isfinite(pts[:,2])
    sc=ax[j].scatter(pts[tagged,0],pts[tagged,1],s=16,c=pts[tagged,2],cmap="viridis",
                     edgecolors="k",linewidths=0.25,zorder=2,label="training configs (colour = T)")
    if (~tagged).any():
        ax[j].scatter(pts[~tagged,0],pts[~tagged,1],s=16,c="k",marker="x",zorder=2,
                      label="training configs (no T tag)")
    ax[j].set_xlabel("mean nearest O-O distance (A)")
    ax[j].set_ylabel("mean first-shell O-O coordination")
    ax[j].set_title(f"{tag} ({T:.0f} K, {rho:.2f} g/cc)")
    ax[j].legend(fontsize=8)
    plt.colorbar(sc,ax=ax[j],label="training cell T (K)")
plt.tight_layout()
fn="figs/coverage_training_vs_production.png"; plt.savefig(fn,dpi=140); plt.close(fig)
print(f"wrote {fn}")
