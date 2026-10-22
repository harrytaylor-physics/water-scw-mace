# analyze_traj.py -- ALL analyses from saved trajectories; never reruns MD, needs no GPU.
#   python analyze_traj.py <tag>                    # analyse every trajectory under traj/<tag>/
#   python analyze_traj.py <tag> --compare          # + finite-size / long-range overlay figure
# Reads traj/<tag>/<which>_n<nmol>.frames.npz written by run_traj.py and produces:
#   * total pressure with block-SE (per trajectory)
#   * O-O RDF out to L/2, and S(k) down to k=2pi/L (per trajectory)
#   * with --compare: one figure overlaying g(r), S(k), and P across box sizes/models
import sys, os, glob, numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import bench_common as bc

tag=sys.argv[1]; do_cmp = "--compare" in sys.argv
files=sorted(glob.glob(f"traj/{tag}/*_n*.frames.npz"))
if not files: sys.exit(f"no trajectories under traj/{tag}/ -- run run_traj.py first")

res=[]
for f in files:
    D=np.load(f)
    which=str(D["which"]); nmol=int(D["nmol"]); L=float(D["L"]); T=float(D["T"]); rho=float(D["rho"])
    framesO=D["framesO"]; Ps=D["Ps"]
    if len(framesO)==0: print(f"{f}: no frames yet (run still in progress?) -- skipping"); continue
    P=float(Ps.mean()); SE=bc.block_se(Ps,10) if len(Ps)>=10 else float("nan")
    rmax=L/2.0
    rc,g=bc.rdf_OO(list(framesO), L, rmax=rmax, nbins=400)
    numdens=framesO.shape[1]/L**3
    dr=rc[1]-rc[0]; h=g-1.0; win=0.5*(1+np.cos(np.pi*rc/rmax))
    kmin=2*np.pi/L; kk=np.linspace(kmin,8.0,300)
    Sk=1.0+4*np.pi*numdens*np.array([np.sum(h*win*rc*np.sin(k*rc)/k)*dr for k in kk])
    tail=np.where(np.abs(h)>0.02)[0]; xi=rc[tail[-1]] if len(tail) else 0.0
    res.append(dict(which=which,nmol=nmol,L=L,T=T,rho=rho,P=P,SE=SE,
                    rc=rc,g=g,kk=kk,Sk=Sk,xi=xi,nfr=len(framesO)))
    print(f"{which:>10} N={nmol:5d} L={L:5.1f}A [{len(framesO)} frames]: "
          f"P={P:+.3f}+/-{SE:.3f} GPa  RDF reach xi={xi:.1f}A (L/2={rmax:.1f})  "
          f"S(kmin={kmin:.2f})={Sk[0]:.2f}")

if do_cmp and len(res)>1:
    fig,ax=plt.subplots(1,3,figsize=(16,5))
    order=sorted(res, key=lambda r:(r["which"],r["nmol"]))
    for r in order:
        lbl=f"{r['which']} N={r['nmol']}"
        ax[0].plot(r["rc"], r["g"], lw=1.2, label=lbl)
        ax[1].plot(r["kk"], r["Sk"], lw=1.2, label=lbl)
    ax[0].axhline(1,c="grey",lw=.5); ax[0].set_xlabel("r (A)"); ax[0].set_ylabel("g_OO(r)")
    ax[0].set_title("RDF out to each box's L/2\n(curves agreeing + flat tails = nothing missed)")
    ax[1].set_xlabel("k (1/A)"); ax[1].set_ylabel("S(k)")
    ax[1].set_title("structure factor to each box's lowest k\n(no low-k upturn = no long-range correlation)")
    ax[1].set_xlim(0, 4)
    for w,mk in [("v4","o"),("foundation","s")]:
        rr=[r for r in order if r["which"]==w]
        if rr: ax[2].errorbar([r["nmol"] for r in rr],[r["P"] for r in rr],
                              yerr=[r["SE"] for r in rr],fmt=mk+"-",capsize=3,label=w)
    Pref=bc.iapws_P_GPa(res[0]["T"], res[0]["rho"])
    if Pref is not None: ax[2].axhline(Pref, ls="--", c="grey", label="IAPWS-95")
    ax[2].set_xscale("log"); ax[2].set_xlabel("N molecules (log)"); ax[2].set_ylabel("P (GPa)")
    ax[2].set_title("pressure vs system size")
    for x in ax: x.legend(fontsize=7)
    plt.tight_layout(); os.makedirs("lr",exist_ok=True)
    plt.savefig(f"lr/{tag}_overlay.png", dpi=140)
    print(f"wrote lr/{tag}_overlay.png")
