# make_missing_figures.py -- the three remaining evidence figures:
#   (1) label distributions: per-atom energies + force magnitudes of the TRAINING SET itself
#   (2) snapshot gallery: rendered atomic configurations from representative pool cells
#   (3) coverage scatter: training snapshots overlaid on the configurations production MD visits
# Usage (from ~/scw):
#   python make_missing_figures.py --train <training_set.xyz>      # (1) needs the labelled set
#   python make_missing_figures.py                                  # (2)+(3) only
# Analysis only: no GPU, no MD. Outputs to figs/. Missing inputs are skipped with a reason.
import sys, os, re, glob, numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from ase.io import read


os.makedirs("figs", exist_ok=True)
MADE=[]; SKIP=[]
def done(fn): MADE.append(fn); print(f"  wrote {fn}")
def skip(w,y): SKIP.append((w,y)); print(f"  SKIP {w}: {y}")

ANCHORS=[("scw673",673.0,0.66),("hot1000",1000.0,0.80)]
FORCE_CUT=50.0     # eV/A hard filter used at training time
RCUT=3.4           # first-shell O-O cutoff (A)

# ---------------- (1) label distributions ----------------
def get_energy(at):
    for k in ("REF_energy","energy","dft_energy","DFT_energy"):
        if k in at.info: return float(at.info[k])
    try: return float(at.get_potential_energy())
    except Exception: return None

def get_forces(at):
    for k in ("REF_forces","forces","dft_forces","DFT_forces"):
        if k in at.arrays: return np.asarray(at.arrays[k])
    try: return np.asarray(at.get_forces())
    except Exception: return None

def get_regime(at,i):
    name=None
    for k in ("config_type","regime","source"):
        if k in at.info: name=str(at.info[k]); break
    if name is None: return "all"
    m=re.search(r'[Tt](\d+)',name)
    if not m: return "other"
    T=int(m.group(1))
    if T<=2250: return "bulk pools (650-2250 K)"
    if T<=3000: return "RCB fill (2500-3000 K)"
    return "hot/AIMD (3000-6000 K)"

def fig_labels(train_path):
    if not train_path: return skip("label distributions","no --train <file> given")
    if not os.path.exists(train_path): return skip("label distributions",f"{train_path} not found")
    ats=read(train_path,index=":")
    E={}; F={}
    nE=nF=0
    for i,at in enumerate(ats):
        reg=get_regime(at,i)
        e=get_energy(at); f=get_forces(at)
        if e is not None:
            E.setdefault(reg,[]).append(e/len(at)); nE+=1
        if f is not None:
            F.setdefault(reg,[]).append(np.linalg.norm(f,axis=1)); nF+=1
    if nE==0 and nF==0:
        return skip("label distributions","no energy/forces found (keys tried: REF_energy/energy/..., REF_forces/forces/...)")
    regs=sorted(set(list(E)+list(F)))
    fig,ax=plt.subplots(1,2,figsize=(13,5))
    for reg in regs:
        if reg in E and len(E[reg])>1:
            ax[0].hist(E[reg],bins=40,histtype="step",lw=1.5,label=f"{reg} ({len(E[reg])})")
    ax[0].set_xlabel("energy per atom (eV)"); ax[0].set_ylabel("configurations")
    ax[0].set_title(f"labelled training set: energies ({nE} configs)")
    if len(regs)>1 or "all" not in regs: ax[0].legend(fontsize=8)
    allF=[]
    for reg in regs:
        if reg in F:
            fmag=np.concatenate(F[reg]); allF.append(fmag)
            ax[1].hist(fmag,bins=np.logspace(-2,np.log10(max(FORCE_CUT*1.5,fmag.max()+1)),60),
                       histtype="step",lw=1.5,label=f"{reg}")
    ax[1].axvline(FORCE_CUT,ls="--",c="red",label=f"{FORCE_CUT:.0f} eV/A cut")
    ax[1].set_xscale("log"); ax[1].set_yscale("log")
    ax[1].set_xlabel("|force| per atom (eV/A)"); ax[1].set_ylabel("atoms")
    ax[1].set_title("labelled training set: force magnitudes")
    ax[1].legend(fontsize=8)
    over=(np.concatenate(allF)>FORCE_CUT).sum() if allF else 0
    plt.tight_layout(); fn="figs/label_distributions.png"; plt.savefig(fn,dpi=140); plt.close(fig); done(fn)
    print(f"    atoms above {FORCE_CUT} eV/A cut in the final set: {over} (should be 0)")

# ---------------- (2) snapshot gallery ----------------
def parse_Trho(stem):
    mT=re.search(r'[Tt](\d+(?:\.\d+)?)',stem); mr=re.search(r'rho([0-9.]+(?:e[+-]?\d+)?)',stem)
    return (float(mT.group(1)) if mT else np.nan, float(mr.group(1)) if mr else np.nan)

def fig_gallery(max_cells=6):
    from ase.visualize.plot import plot_atoms
    xyzs=[]
    for d in ("pools","pools_rcb"):
        xyzs+=sorted(glob.glob(os.path.join(d,"pool*.xyz")))
    if not xyzs: return skip("snapshot gallery","no pool*.xyz under pools/ or pools_rcb/")
    cells=[]
    for x in xyzs:
        stem=os.path.splitext(os.path.basename(x))[0]; T,rho=parse_Trho(stem)
        cells.append((T,rho,x,stem))
    cells.sort()
    # pick spread: coolest-dilute, coolest-dense, mid, hottest-dilute, hottest-dense...
    idx=np.unique(np.linspace(0,len(cells)-1,min(max_cells,len(cells))).astype(int))
    pick=[cells[i] for i in idx]
    ncol=3; nrow=int(np.ceil(len(pick)/ncol))
    fig,axs=plt.subplots(nrow,ncol,figsize=(4.6*ncol,4.2*nrow),squeeze=False)
    for i,(T,rho,x,stem) in enumerate(pick):
        at=read(x,index=-1); at.wrap()            # final frame, wrapped into cell
        a=axs[i//ncol][i%ncol]
        plot_atoms(at,a,radii=0.45,rotation=("12x,8y,0z"))
        a.set_title(f"T={T:.0f} K, rho={rho:.0f} kg/m3",fontsize=10)
        a.set_xticks([]); a.set_yticks([])
    for j in range(len(pick),nrow*ncol): axs[j//ncol][j%ncol].axis("off")
    fig.suptitle("training configurations: final frame of representative cells")
    plt.tight_layout(rect=[0,0,1,0.95])
    fn="figs/snapshot_gallery.png"; plt.savefig(fn,dpi=130); plt.close(fig); done(fn)

# ---------------- (3) coverage: training snapshots vs production visits ----------------
def descriptors_O(O,L):
    """per-frame (mean nearest O-O distance, mean first-shell coordination)"""
    d=O[:,None,:]-O[None,:,:]; d-=L*np.round(d/L)
    r=np.sqrt((d*d).sum(-1)); np.fill_diagonal(r,np.inf)
    return float(r.min(1).mean()), float((r<RCUT).sum(1).mean())

def pool_selected_descriptors(per_cell=12):
    pts=[]
    for d in ("pools","pools_rcb"):
        for x in sorted(glob.glob(os.path.join(d,"pool*.xyz"))):
            frames=read(x,index=":")
            c=np.array(frames[0].cell); L=float(np.mean(np.diag(c))) if np.any(c) else None
            if L is None: continue
            isO=np.array(frames[0].get_chemical_symbols())=="O"
            sel=np.linspace(0,len(frames)-1,per_cell).astype(int)
            stem=os.path.splitext(os.path.basename(x))[0]; T,_=parse_Trho(stem)
            for s in sel:
                O=frames[s].get_positions()[isO]
                pts.append((*descriptors_O(O,L),T))
    return np.array(pts) if pts else None

def fig_coverage(per_cell=12):
    train=pool_selected_descriptors(per_cell)
    if train is None: return skip("coverage scatter","no pool cells found")
    fig,ax=plt.subplots(1,len(ANCHORS),figsize=(6.5*len(ANCHORS),5.5),squeeze=False); ax=ax[0]
    any_prod=False
    for j,(tag,T,rho) in enumerate(ANCHORS):
        f=f"traj/{tag}/v4_n64.frames.npz"
        if os.path.exists(f):
            D=np.load(f); L=float(D["L"])
            prod=np.array([descriptors_O(O,L) for O in D["framesO"]])
            ax[j].scatter(prod[:,0],prod[:,1],s=8,c="0.7",label=f"production MD visits ({len(prod)} frames)")
            any_prod=True
        else:
            ax[j].text(.5,.5,"no production trajectory",ha="center",transform=ax[j].transAxes)
        sc=ax[j].scatter(train[:,0],train[:,1],s=18,c=train[:,2],cmap="viridis",
                         edgecolors="k",linewidths=0.3,label="training snapshots (colour = T)")
        ax[j].set_xlabel("mean nearest O-O distance (A)")
        ax[j].set_ylabel("mean first-shell O-O coordination")
        ax[j].set_title(f"{tag} ({T:.0f} K, {rho:.2f} g/cc)")
        ax[j].legend(fontsize=8)
        plt.colorbar(sc,ax=ax[j],label="training cell T (K)")
    if not any_prod: skip("coverage scatter","no traj/<tag>/v4_n64.frames.npz for any anchor")
    plt.tight_layout(); fn="figs/coverage_training_vs_production.png"
    plt.savefig(fn,dpi=140); plt.close(fig); done(fn)

# ---------------- run ----------------
train_path=None
if "--train" in sys.argv:
    i=sys.argv.index("--train"); train_path=sys.argv[i+1]
per_cell=12
if "--per_cell" in sys.argv:
    i=sys.argv.index("--per_cell"); per_cell=int(sys.argv[i+1])

print("== (1) label distributions =="); fig_labels(train_path)
print("== (2) snapshot gallery ==");    fig_gallery()
print("== (3) coverage scatter ==");    fig_coverage(per_cell)

print(f"\nMANIFEST: {len(MADE)} figures")
for m in MADE: print(f"  + {m}")
for w,y in SKIP: print(f"  - skipped {w}: {y}")
