# make_training_figures.py -- the two model-quality figures:
#   (A) training curve from the MACE log  -> figs/training_curve.png   (overfitting check)
#   (B) force parity: model vs DFT labels -> figs/force_parity.png    (force accuracy)
# Usage (from ~/scw):
#   python make_training_figures.py --log logs/<v4_log>                      # (A) login node, instant
#   python make_training_figures.py --model water_scw_v4.model --train train_v4.xyz   # (B) needs GPU
# Both can be given at once. Outputs to figs/.
import sys, os, re, numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

os.makedirs("figs", exist_ok=True)

def parse_args():
    a=sys.argv[1:]; d={}
    for k in ("--log","--model","--train"):
        if k in a:
            i=a.index(k); d[k[2:]]=a[i+1]
    d["swa_start"]=int(a[a.index("--swa_start")+1]) if "--swa_start" in a else 120
    return d

# ---------------- (A) training curve from MACE log ----------------
def fig_training_curve(logpath, swa_start):
    txt=open(logpath, errors="ignore").read()
    # MACE logs per-epoch validation lines in a few formats across versions; try them all.
    pats=[
        # e.g. "Epoch 12: loss=0.0123, RMSE_E_per_atom=5.3 meV, RMSE_F=42.1 meV / A"
        r"Epoch\s+(\d+):\s*loss=\s*([\deE+.-]+).*?RMSE_E(?:_per_atom)?=\s*([\deE+.-]+)\s*meV.*?RMSE_F=\s*([\deE+.-]+)",
        # e.g. "Epoch 12: head: Default, loss=..., RMSE_E_per_atom=..., RMSE_F=..."
        r"Epoch\s+(\d+).*?loss=\s*([\deE+.-]+).*?RMSE_E(?:_per_atom)?=\s*([\deE+.-]+).*?RMSE_F=\s*([\deE+.-]+)",
    ]
    rows=[]
    for p in pats:
        rows=[(int(e),float(l),float(re_),float(rf)) for e,l,re_,rf in re.findall(p,txt)]
        if rows: break
    if not rows:
        print("  SKIP training curve: no per-epoch lines matched.")
        print("  Paste 3 example 'Epoch ...' lines from the log and I will adapt the parser.")
        ex=[l for l in txt.splitlines() if re.search(r"Epoch\s+\d",l)][:3]
        for l in ex: print("   sample:",l[:160])
        return
    rows.sort()
    ep=np.array([r[0] for r in rows]); loss=np.array([r[1] for r in rows])
    rmseE=np.array([r[2] for r in rows]); rmseF=np.array([r[3] for r in rows])
    # keep last entry per epoch (logs sometimes repeat on restart)
    _,idx=np.unique(ep[::-1],return_index=True); idx=len(ep)-1-idx
    ep,loss,rmseE,rmseF=ep[idx],loss[idx],rmseE[idx],rmseF[idx]

    fig,ax=plt.subplots(1,2,figsize=(13,5))
    ax[0].semilogy(ep,loss,lw=1.4)
    ax[0].axvline(swa_start,ls=":",c="grey"); ax[0].text(swa_start,ax[0].get_ylim()[1],"SWA",fontsize=8,va="top")
    ax[0].set_xlabel("epoch"); ax[0].set_ylabel("training loss (log)"); ax[0].set_title("training loss")
    ax[1].plot(ep,rmseF,lw=1.4,label="validation force RMSE (meV/A)")
    ax2=ax[1].twinx(); ax2.plot(ep,rmseE,lw=1.2,c="C1",label="validation energy RMSE (meV/atom)")
    ax[1].axvline(swa_start,ls=":",c="grey")
    ax[1].set_xlabel("epoch"); ax[1].set_ylabel("force RMSE (meV/A)")
    ax2.set_ylabel("energy RMSE (meV/atom)")
    ax[1].set_title("validation error vs epoch (10% held-out split)\nflat/decreasing = no overfitting")
    h1,l1=ax[1].get_legend_handles_labels(); h2,l2=ax2.get_legend_handles_labels()
    ax[1].legend(h1+h2,l1+l2,fontsize=8)
    plt.tight_layout(); fn="figs/training_curve.png"; plt.savefig(fn,dpi=140); plt.close(fig)
    print(f"  wrote {fn}   (final: loss={loss[-1]:.4g}, F-RMSE={rmseF[-1]:.1f} meV/A, E-RMSE={rmseE[-1]:.2f} meV/atom)")

# ---------------- (B) force parity: model vs DFT labels ----------------
def fig_force_parity(model_path, train_path):
    from ase.io import read
    from mace.calculators import MACECalculator
    ats=read(train_path,index=":")
    calc=MACECalculator(model_paths=model_path, device="cuda", default_dtype="float32")
    ref=[]; pred=[]; reg=[]
    for i,at in enumerate(ats):
        F=None
        for k in ("forces","REF_forces","dft_forces"):
            if k in at.arrays: F=np.asarray(at.arrays[k]); break
        if F is None and at.calc is not None:
            try: F=np.asarray(at.get_forces())
            except Exception: F=None
        if F is None: continue
        at2=at.copy(); at2.calc=calc
        Fp=at2.get_forces()
        ref.append(F.ravel()); pred.append(Fp.ravel())
        name=str(at.info.get("config_type","all"))
        m=re.search(r"[Tt](\d{3,4})",name)
        T=int(m.group(1)) if m else None
        reg.append("bulk pools (650-2250 K)" if T and T<=2250 else
                   "RCB fill (2500-3000 K)" if T and T<=3000 else
                   "hot/AIMD (3000-6000 K)" if T else "other")
        if (i+1)%100==0: print(f"    {i+1}/{len(ats)} configs evaluated", flush=True)
    if not ref: print("  SKIP force parity: no force labels found"); return
    regs=sorted(set(reg))
    fig,ax=plt.subplots(figsize=(7.5,7))
    lims=[0,0]
    print("  per-regime force RMSE (meV/A):")
    for rname in regs:
        R=np.concatenate([r for r,g in zip(ref,reg) if g==rname])
        P=np.concatenate([p for p,g in zip(pred,reg) if g==rname])
        rmse=1000*np.sqrt(np.mean((P-R)**2))
        n=min(len(R),20000)                       # subsample for plotting only
        j=np.random.default_rng(0).choice(len(R),n,replace=False)
        ax.plot(R[j],P[j],".",ms=2,alpha=0.35,label=f"{rname} (RMSE {rmse:.0f} meV/A)")
        lims=[min(lims[0],R.min()),max(lims[1],R.max())]
        print(f"    {rname:>14}: {rmse:8.1f}   ({len(R)//3} atoms)")
    ax.plot(lims,lims,"k-",lw=0.8)
    ax.set_xlabel("DFT force component (eV/A)"); ax.set_ylabel("model force component (eV/A)")
    ax.set_title("force parity: water_scw_v4 vs revPBE0-D3 labels")
    ax.legend(fontsize=8,markerscale=4)
    plt.tight_layout(); fn="figs/force_parity.png"; plt.savefig(fn,dpi=140); plt.close(fig)
    print(f"  wrote {fn}")

d=parse_args()
if "log" in d:
    print("== (A) training curve =="); fig_training_curve(d["log"], d["swa_start"])
if "model" in d and "train" in d:
    print("== (B) force parity =="); fig_force_parity(d["model"], d["train"])
if "log" not in d and "model" not in d:
    print(__doc__ or "usage: --log <macelog> and/or --model <model> --train <xyz>")
