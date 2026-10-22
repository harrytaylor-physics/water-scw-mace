# pool_overview.py -- one-glance physicality map across all pool cells.
#   python pool_overview.py <dir> [--glob 'pool*.xyz']
# Builds a (T, rho) grid map coloured by minimum interatomic distance (the blow-up / bad-
# geometry detector, robust at every density) plus per-cell min-distance traces, and prints
# a worst-first table. No energies needed -- this is about geometry, not thermodynamics.
import sys, os, re, numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from ase.io import read
from analyze_aimd import discover, L_from_rho_path, iact   # reuse the tested helpers

OH = 0.96          # O-H bond length (Angstrom): healthy min-distance sits near here
DANGER = 0.70      # below this => atoms crashing together / diverged geometry

def parse_Trho(stem):
    mT=re.search(r'[Tt](\d+(?:\.\d+)?)', stem)
    mr=re.search(r'rho([0-9.]+(?:e[+-]?\d+)?)', stem)
    return (float(mT.group(1)) if mT else np.nan,
            float(mr.group(1)) if mr else np.nan)

def get_L(frames, xyz, nat):
    c=np.array(frames[0].cell)
    if np.any(c): return float(np.mean(np.diag(c)))
    return L_from_rho_path(xyz, nat//3)

def min_dist_trace(P, L):
    """Per-frame global minimum interatomic distance (minimum image)."""
    out=np.empty(len(P))
    for k,fr in enumerate(P):
        d=fr[:,None,:]-fr[None,:,:]
        if L: d-=L*np.round(d/L)
        r=np.sqrt((d*d).sum(-1)); np.fill_diagonal(r, np.inf)
        out[k]=r.min()
    return out

def coord_Neff(O, L, rcut=3.4):
    """Robust diversity proxy: decorrelation of the mean O-O coordination number,
    which fluctuates on the physical collision timescale in both gas and liquid."""
    c=np.empty(len(O))
    for k,o in enumerate(O):
        d=o[:,None,:]-o[None,:,:]
        if L: d-=L*np.round(d/L)
        r=np.sqrt((d*d).sum(-1)); np.fill_diagonal(r, np.inf)
        c[k]=(r<rcut).sum()/len(o)
    return len(c)/iact(c)

def main():
    args=sys.argv[1:]; pattern=None
    if "--glob" in args:
        i=args.index("--glob"); pattern=args[i+1]; del args[i:i+2]
    if not args: sys.exit("usage: python pool_overview.py <dir> [--glob 'pool*.xyz']")
    root=args[0]
    jobs=discover(root, pattern)
    if not jobs: sys.exit(f"no pool/xyz files found under {root}")

    rows=[]; traces=[]
    for xyz,_ in jobs:
        stem=os.path.splitext(os.path.basename(xyz))[0]
        frames=read(xyz, index=":"); nat=len(frames[0])
        sym=np.array(frames[0].get_chemical_symbols()); isO=sym=="O"
        P=np.array([f.get_positions() for f in frames])
        finite = np.isfinite(P).all()
        L=get_L(frames, xyz, nat)
        T,rho=parse_Trho(stem)
        md=min_dist_trace(P, L) if finite else np.array([np.nan])
        dmin=float(np.nanmin(md)); frac_bad=float(np.mean(md<DANGER))
        O=P[:,isO,:]
        U=O.copy()
        if L:
            for k in range(1,len(O)):
                dd=O[k]-O[k-1]; dd-=L*np.round(dd/L); U[k]=U[k-1]+dd
        msd=float(((U-U[0])**2).sum(-1).mean(1)[-1])
        neff=coord_Neff(O, L)
        healthy = finite and dmin>=DANGER
        rows.append(dict(stem=stem,T=T,rho=rho,nat=nat,L=L,dmin=dmin,
                         frac_bad=frac_bad,msd=msd,neff=neff,finite=finite,healthy=healthy))
        traces.append((stem, np.linspace(0,1,len(md)), md, healthy))

    rows.sort(key=lambda r: (r["finite"], r["dmin"]))   # worst first

    # ---- figure: (T,rho) health map  +  min-distance traces ----
    fig,ax=plt.subplots(1,2,figsize=(14,6))
    rr=np.array([r["rho"] for r in rows]); TT=np.array([r["T"] for r in rows])
    dm=np.array([r["dmin"] for r in rows])
    sc=ax[0].scatter(rr, TT, c=np.clip(dm,0.4,1.0), cmap="RdYlGn", vmin=0.5, vmax=1.0,
                     s=160, edgecolors="k", linewidths=0.6, zorder=3)
    for r in rows:
        if not r["healthy"]:
            ax[0].scatter([r["rho"]],[r["T"]], s=320, facecolors="none",
                          edgecolors="red", linewidths=2.0, zorder=4)
            ax[0].annotate(f"{r['dmin']:.2f}", (r["rho"],r["T"]),
                           textcoords="offset points", xytext=(7,6), fontsize=7, color="red")
    ax[0].set_xscale("log"); ax[0].set_xlabel("density (kg/m^3, log)"); ax[0].set_ylabel("temperature (K)")
    ax[0].set_title("pool physicality map\ncolour = min interatomic distance (A); red ring = blow-up")
    cb=plt.colorbar(sc, ax=ax[0]); cb.set_label("min interatomic distance (A)")

    for stem,x,md,healthy in traces:
        ax[1].plot(x, md, lw=(1.6 if not healthy else 0.7),
                   color=("red" if not healthy else "0.6"),
                   alpha=(0.9 if not healthy else 0.5),
                   label=(stem if not healthy else None))
    ax[1].axhline(OH, ls=":", c="green", label=f"O-H bond ~{OH:.2f} A")
    ax[1].axhline(DANGER, ls="--", c="red", label=f"danger <{DANGER:.2f} A")
    ax[1].set_xlabel("fraction of trajectory"); ax[1].set_ylabel("min interatomic distance (A)")
    ax[1].set_ylim(0, max(OH*1.4, np.nanmax([t[2].max() for t in traces if np.isfinite(t[2]).any()])))
    ax[1].set_title("per-cell min-distance traces (red = ever dips into danger)")
    if any(not h for *_,h in traces): ax[1].legend(fontsize=7, loc="upper right")

    plt.tight_layout(); plt.savefig("pool_overview.png", dpi=140)
    print("wrote pool_overview.png\n")

    # ---- worst-first table ----
    print(f"{'cell':<26}{'T':>6}{'rho':>9}{'N':>5}{'minDist':>9}{'%<0.7':>7}{'MSD':>8}{'coordNeff':>10}  flag")
    for r in rows:
        flag = "BLOW-UP" if not r["finite"] else ("crash" if r["dmin"]<DANGER else "ok")
        print(f"{r['stem']:<26}{r['T']:>6.0f}{r['rho']:>9.1f}{r['nat']//3:>5}"
              f"{r['dmin']:>9.2f}{100*r['frac_bad']:>6.1f}%{r['msd']:>8.0f}{r['neff']:>10.0f}  {flag}")

if __name__=="__main__":
    main()
