# make_all_figures.py -- ONE script that generates ALL validation figures into figs/.
#   python make_all_figures.py            (run from ~/scw; analysis only, no GPU, no MD)
# Reads only data that already exists:
#   bench/<tag>_<model>/bench.npz   -> anchor-point RDF + VDOS vs experiment
#   traj/<tag>/<model>_n*.frames.npz-> anchor pressures, finite-size/long-range overlays
#   eos/<tag>/*.npz                 -> isotherm P(rho) vs IAPWS/AQUA + bulk modulus
#   pools/, pools_rcb/  *.xyz       -> physicality maps + per-cell audits + RDF grid
# Missing inputs are skipped with a message; a manifest of made/skipped prints at the end.
import os, re, glob, numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import bench_common as bc

ANCHORS=[("scw673",673.0,0.66),("hot1000",1000.0,0.80)]
POOL_DIRS=[("pools",12),("pools_rcb",12)]          # (dir, per_cell selection to mark)
os.makedirs("figs",exist_ok=True)
MADE=[]; SKIP=[]

AMB_RDF1,AMB_RDF2=2.80,4.5
EXP_BANDS={"H-bond":180,"libration":600,"bend":1650,"stretch(liq)":3400,"stretch(gas)":3657}

def done(fn): MADE.append(fn); print(f"  wrote {fn}")
def skip(what,why): SKIP.append((what,why)); print(f"  SKIP {what}: {why}")

# ---------- 1) anchor points: RDF + VDOS vs experiment ----------
def fig_bench(tag,T,rho):
    try:
        v4=np.load(f"bench/{tag}_v4/bench.npz"); fo=np.load(f"bench/{tag}_foundation/bench.npz")
    except Exception as e: return skip(f"bench {tag}",f"missing bench.npz ({e})")
    fig,ax=plt.subplots(1,2,figsize=(13,5))
    ax[0].plot(v4["rc"],v4["gOO"],lw=1.6,label="v4 (fine-tuned)")
    ax[0].plot(fo["rc"],fo["gOO"],lw=1.6,label="foundation")
    ax[0].axvline(AMB_RDF1,ls=":",c="grey",alpha=.6); ax[0].axvline(AMB_RDF2,ls=":",c="grey",alpha=.6)
    ax[0].set_xlabel("r (A)"); ax[0].set_ylabel("g_OO(r)")
    ax[0].set_title(f"{tag} {T:.0f} K {rho:.2f} g/cc: O-O RDF (dotted = ambient guides)")
    ax[0].legend()
    ax[1].plot(v4["wn"],v4["dos"]/v4["dos"].max(),lw=1.3,label="v4")
    ax[1].plot(fo["wn"],fo["dos"]/fo["dos"].max(),lw=1.3,label="foundation")
    for n,x in EXP_BANDS.items():
        ax[1].axvline(x,ls=":",c="grey",alpha=.5); ax[1].text(x,1.02,n,rotation=90,fontsize=7)
    ax[1].set_xlim(0,4200); ax[1].set_xlabel("cm-1"); ax[1].set_ylabel("VDOS (norm)")
    ax[1].set_title("vibrations vs experimental band positions"); ax[1].legend()
    plt.tight_layout(); fn=f"figs/bench_{tag}.png"; plt.savefig(fn,dpi=140); plt.close(fig); done(fn)

# ---------- 2) anchor pressures vs experiment (the missing figure) ----------
def fig_anchor_pressure():
    fig,ax=plt.subplots(figsize=(8,5)); got=False
    for i,(tag,T,rho) in enumerate(ANCHORS):
        vals={}
        for w in ("v4","foundation"):
            f=f"traj/{tag}/{w}_n64.frames.npz"
            if os.path.exists(f):
                D=np.load(f); vals[w]=(float(D["Ps"].mean()), bc.block_se(D["Ps"],10))
        if not vals: skip(f"anchor pressure {tag}","no traj/*_n64.frames.npz"); continue
        got=True
        Pi=bc.iapws_P_GPa(T,rho); Pa,_,_=bc.aqua_P_GPa(T,rho)
        x0=i*2.0
        if "v4" in vals: ax.errorbar([x0-0.25],[vals["v4"][0]],yerr=[vals["v4"][1]],fmt="o",c="C0",capsize=4,label="v4" if i==0 else None)
        if "foundation" in vals: ax.errorbar([x0+0.25],[vals["foundation"][0]],yerr=[vals["foundation"][1]],fmt="s",c="C1",capsize=4,label="foundation" if i==0 else None)
        if Pi is not None: ax.hlines(Pi,x0-0.55,x0+0.55,colors="k",linestyles="--",label="IAPWS-95" if i==0 else None)
        if Pa is not None: ax.hlines(Pa,x0-0.55,x0+0.55,colors="grey",linestyles=":",label="AQUA" if i==0 else None)
    if not got: plt.close(fig); return
    ax.axhline(0,c="grey",lw=.5)
    ax.set_xticks([0,2.0]); ax.set_xticklabels([f"{t} ({T:.0f} K, {r:.2f} g/cc)" for t,T,r in ANCHORS],fontsize=9)
    ax.set_ylabel("P (GPa, total)"); ax.set_title("anchor-point pressure vs experimental references")
    ax.legend(fontsize=9); plt.tight_layout()
    fn="figs/pressure_anchor_points.png"; plt.savefig(fn,dpi=140); plt.close(fig); done(fn)

# ---------- 3) isotherms: P(rho) + bulk modulus vs references ----------
def fig_isotherm(tag,T):
    files=sorted(glob.glob(f"eos/{tag}/*_rho*.npz"))
    if not files: return skip(f"isotherm {tag}","no eos/<tag>/*.npz")
    data={"v4":[],"foundation":[]}
    for f in files:
        d=np.load(f); data[str(d["which"])].append((float(d["rho"]),float(d["P"]),float(d["SE"]),bool(d["converged"])))
    for w in data: data[w].sort()
    rhos=sorted({r for w in data for (r,_,_,_) in data[w]})
    rr=np.linspace(min(rhos),max(rhos),60)
    iap=[bc.iapws_P_GPa(T,r) for r in rr]; aq=[bc.aqua_P_GPa(T,r)[0] for r in rr]
    fig,ax=plt.subplots(1,2,figsize=(13,5))
    ax[0].plot(rr,[p if p is not None else np.nan for p in iap],"k-",lw=1.2,label="IAPWS-95")
    ax[0].plot(rr,[p if p is not None else np.nan for p in aq],"--",c="grey",lw=1.2,label="AQUA")
    for w,mk in [("v4","o"),("foundation","s")]:
        if not data[w]: continue
        R=np.array([x[0] for x in data[w]]); P=np.array([x[1] for x in data[w]]); E=np.array([x[2] for x in data[w]])
        C=np.array([x[3] for x in data[w]])
        ax[0].errorbar(R,P,yerr=E,fmt=mk+"-",capsize=3,label=w)
        if (~C).any(): ax[0].scatter(R[~C],P[~C],s=130,facecolors="none",edgecolors="red")
    ax[0].axhline(0,c="grey",lw=.5); ax[0].set_xlabel("rho (g/cc)"); ax[0].set_ylabel("P (GPa)")
    ax[0].set_title(f"{tag}: {T:.0f} K isotherm (red ring = not converged)"); ax[0].legend(fontsize=8)
    def K(R,P): R=np.asarray(R); return R*np.gradient(np.asarray(P),R)
    for w,mk in [("v4","o"),("foundation","s")]:
        if data[w]:
            R=[x[0] for x in data[w]]; P=[x[1] for x in data[w]]
            ax[1].plot(R,K(R,P),mk+"-",label=w)
    ax[1].plot(rr,K(rr,np.array([p if p is not None else np.nan for p in aq])),"--",c="grey",label="AQUA")
    ax[1].set_xlabel("rho (g/cc)"); ax[1].set_ylabel("K_T (GPa)")
    ax[1].set_title("bulk modulus (slope; offset-insensitive)"); ax[1].legend(fontsize=8)
    plt.tight_layout(); fn=f"figs/isotherm_{tag}.png"; plt.savefig(fn,dpi=140); plt.close(fig); done(fn)

# ---------- 4) finite-size / long-range overlays ----------
def fig_overlay(tag,T,rho):
    files=sorted(glob.glob(f"traj/{tag}/*_n*.frames.npz"))
    if not files: return skip(f"overlay {tag}","no traj/<tag>/*.frames.npz")
    res=[]
    for f in files:
        D=np.load(f); framesO=D["framesO"]
        if len(framesO)==0: continue
        L=float(D["L"]); numdens=framesO.shape[1]/L**3
        rc,g=bc.rdf_OO(list(framesO),L,rmax=L/2,nbins=400)
        dr=rc[1]-rc[0]; h=g-1.0; win=0.5*(1+np.cos(np.pi*rc/(L/2)))
        kk=np.linspace(2*np.pi/L,8.0,300)
        Sk=1.0+4*np.pi*numdens*np.array([np.sum(h*win*rc*np.sin(k*rc)/k)*dr for k in kk])
        res.append(dict(which=str(D["which"]),n=int(D["nmol"]),rc=rc,g=g,kk=kk,Sk=Sk,
                        P=float(D["Ps"].mean()),SE=bc.block_se(D["Ps"],10)))
    if not res: return skip(f"overlay {tag}","frames files empty")
    res.sort(key=lambda r:(r["which"],r["n"]))
    fig,ax=plt.subplots(1,3,figsize=(16,5))
    for r in res:
        ax[0].plot(r["rc"],r["g"],lw=1.1,label=f"{r['which']} N={r['n']}")
        ax[1].plot(r["kk"],r["Sk"],lw=1.1,label=f"{r['which']} N={r['n']}")
    ax[0].axhline(1,c="grey",lw=.5); ax[0].set_xlabel("r (A)"); ax[0].set_ylabel("g_OO(r)")
    ax[0].set_title("RDF out to each box's L/2"); ax[1].set_xlim(0,4)
    ax[1].set_xlabel("k (1/A)"); ax[1].set_ylabel("S(k)"); ax[1].set_title("structure factor to lowest k")
    for w,mk in [("v4","o"),("foundation","s")]:
        rr=[r for r in res if r["which"]==w]
        if rr: ax[2].errorbar([r["n"] for r in rr],[r["P"] for r in rr],yerr=[r["SE"] for r in rr],
                              fmt=mk+"-",capsize=3,label=w)
    Pref=bc.iapws_P_GPa(T,rho)
    if Pref is not None: ax[2].axhline(Pref,ls="--",c="grey",label="IAPWS-95")
    ax[2].set_xscale("log"); ax[2].set_xlabel("N molecules (log)"); ax[2].set_ylabel("P (GPa)")
    ax[2].set_title("pressure vs system size")
    for x in ax: x.legend(fontsize=7)
    plt.tight_layout(); fn=f"figs/overlay_{tag}.png"; plt.savefig(fn,dpi=140); plt.close(fig); done(fn)

# ---------- 5) pool physicality map + per-cell audits + RDF grid ----------
def parse_Trho(stem):
    mT=re.search(r'[Tt](\d+(?:\.\d+)?)',stem); mr=re.search(r'rho([0-9.]+(?:e[+-]?\d+)?)',stem)
    return (float(mT.group(1)) if mT else np.nan, float(mr.group(1)) if mr else np.nan)

def audit_pools(pdir,per_cell,grid_res):
    from ase.io import read
    xyzs=sorted(glob.glob(os.path.join(pdir,"pool*.xyz")))
    if not xyzs: return skip(f"audit {pdir}","no pool*.xyz")
    mapr=[]
    for xyz in xyzs:
        stem=os.path.splitext(os.path.basename(xyz))[0]
        frames=read(xyz,index=":"); c=np.array(frames[0].cell)
        L=float(np.mean(np.diag(c))) if np.any(c) else None
        sym=np.array(frames[0].get_chemical_symbols()); isO=sym=="O"
        P=np.array([f.get_positions() for f in frames]); O=P[:,isO,:]; nfr=len(P)
        def mind(fr):
            d=fr[:,None,:]-fr[None,:,:]; d-=L*np.round(d/L)
            r=np.sqrt((d*d).sum(-1)); np.fill_diagonal(r,np.inf); return r.min()
        md=np.array([mind(fr) for fr in P])
        def cnum(o):
            d=o[:,None,:]-o[None,:,:]; d-=L*np.round(d/L)
            r=np.sqrt((d*d).sum(-1)); np.fill_diagonal(r,np.inf); return (r<3.4).sum(1).mean()
        cn=np.array([cnum(o) for o in O])
        sel=np.linspace(0,nfr-1,per_cell).astype(int)
        thirds=[range(0,nfr//3),range(nfr//3,2*nfr//3),range(2*nfr//3,nfr)]
        rdfs=[bc.rdf_OO([O[i] for i in t][::max(1,len(t)//40)],L,rmax=min(6.0,L/2),nbins=150) for t in thirds]
        U=O.copy()
        for k in range(1,nfr):
            dd=O[k]-O[k-1]; dd-=L*np.round(dd/L); U[k]=U[k-1]+dd
        msd=((U-U[0])**2).sum(-1).mean(1)
        fig,ax=plt.subplots(2,2,figsize=(12,8)); fig.suptitle(stem)
        ax[0,0].plot(md,lw=0.8); ax[0,0].axhline(0.96,ls=":",c="green",label="O-H bond")
        ax[0,0].axhline(0.70,ls="--",c="red",label="danger")
        ax[0,0].scatter(sel,md[sel],s=25,c="k",zorder=5,label="selected frames")
        ax[0,0].set_xlabel("frame"); ax[0,0].set_ylabel("min dist (A)")
        ax[0,0].set_title(f"physicality (min={md.min():.2f} A)"); ax[0,0].legend(fontsize=7)
        for (rc,g),lab in zip(rdfs,["first third","middle third","last third"]):
            ax[0,1].plot(rc,g,lw=1.2,label=lab)
        ax[0,1].axhline(1,c="grey",lw=.5); ax[0,1].set_xlabel("r (A)"); ax[0,1].set_ylabel("g_OO(r)")
        ax[0,1].set_title("structure by time slice (overlap = stationary)"); ax[0,1].legend(fontsize=7)
        ax[1,0].plot(msd,lw=1.0); ax[1,0].scatter(sel,msd[sel],s=25,c="k",zorder=5)
        ax[1,0].set_xlabel("frame"); ax[1,0].set_ylabel("O MSD (A^2)")
        ax[1,0].set_title(f"exploration (reaches {msd[-1]:.0f} A^2)")
        ax[1,1].hist(cn,bins=30,color="0.7",label="all frames")
        for s_ in sel: ax[1,1].axvline(cn[s_],c="k",lw=1.0,alpha=0.8)
        cover=(cn[sel].max()-cn[sel].min())/max(cn.max()-cn.min(),1e-9)
        ax[1,1].set_xlabel("first-shell O-O coordination"); ax[1,1].set_ylabel("frames")
        ax[1,1].set_title(f"selected frames representative? span covered: {100*cover:.0f}%")
        plt.tight_layout(rect=[0,0,1,0.96])
        fn=f"figs/audit_{stem}.png"; plt.savefig(fn,dpi=130); plt.close(fig); done(fn)
        T,rho=parse_Trho(stem)
        mapr.append(dict(stem=stem,T=T,rho=rho,dmin=md.min(),rc=rdfs[2][0],g=rdfs[2][1]))
        grid_res.append(mapr[-1])
    # physicality map for this directory
    fig,ax=plt.subplots(figsize=(7,5.5))
    rr=np.array([r["rho"] for r in mapr]); TT=np.array([r["T"] for r in mapr]); dm=np.array([r["dmin"] for r in mapr])
    sc=ax.scatter(rr,TT,c=np.clip(dm,0.4,1.0),cmap="RdYlGn",vmin=0.5,vmax=1.0,s=160,edgecolors="k",linewidths=0.6)
    for r in mapr:
        if r["dmin"]<0.70:
            ax.scatter([r["rho"]],[r["T"]],s=320,facecolors="none",edgecolors="red",linewidths=2)
    ax.set_xscale("log"); ax.set_xlabel("density (kg/m^3, log)"); ax.set_ylabel("temperature (K)")
    ax.set_title(f"{pdir}: physicality map (colour = min interatomic distance)")
    plt.colorbar(sc,ax=ax,label="min dist (A)"); plt.tight_layout()
    fn=f"figs/pool_map_{pdir}.png"; plt.savefig(fn,dpi=140); plt.close(fig); done(fn)

def fig_rdf_grid(grid_res):
    if not grid_res: return skip("rdf grid","no pool cells audited")
    res=sorted(grid_res,key=lambda r:(r["T"],r["rho"]))
    n=len(res); ncol=min(6,n); nrow=int(np.ceil(n/ncol))
    fig,axs=plt.subplots(nrow,ncol,figsize=(2.8*ncol,2.3*nrow),squeeze=False)
    for i,r in enumerate(res):
        a=axs[i//ncol][i%ncol]; a.plot(r["rc"],r["g"],lw=1.0); a.axhline(1,c="grey",lw=.4)
        a.set_title(f"T={r['T']:.0f} rho={r['rho']:.0f}",fontsize=8); a.set_xlim(0,6); a.tick_params(labelsize=6)
    for j in range(n,nrow*ncol): axs[j//ncol][j%ncol].axis("off")
    fig.suptitle("g_OO(r) across all training cells")
    plt.tight_layout(rect=[0,0,1,0.95]); fn="figs/audit_rdf_grid.png"; plt.savefig(fn,dpi=130); plt.close(fig); done(fn)

# ---------------- run everything ----------------
print("== anchor-point benchmarks (RDF + VDOS) ==")
for tag,T,rho in ANCHORS: fig_bench(tag,T,rho)
print("== anchor-point pressures vs experiment ==")
fig_anchor_pressure()
print("== isotherms ==")
for tag,T,rho in ANCHORS: fig_isotherm(tag,T)
print("== finite-size / long-range overlays ==")
for tag,T,rho in ANCHORS: fig_overlay(tag,T,rho)
print("== training-pool audits ==")
grid_res=[]
for pdir,pc in POOL_DIRS: audit_pools(pdir,pc,grid_res)
fig_rdf_grid(grid_res)

print(f"\nMANIFEST: {len(MADE)} figures written to figs/")
for m in MADE: print(f"  + {m}")
if SKIP:
    print("skipped:")
    for w,y in SKIP: print(f"  - {w}: {y}")
