# analyze_aimd.py -- quick health check on the AIMD trajectories that fed the ML training.
#   python analyze_aimd.py <path> [L_angstrom]
# <path> is a CP2K run directory (auto-finds *-pos-1.xyz and *-1.ener) or a single xyz file.
# Optional L (cubic box edge, Angstrom) enables minimum-image structural descriptors when the
# xyz carries no Lattice header.  No GPU needed -- runs on the login node.
#
# It answers: did the run converge / stay stable (energy + temperature + conserved-qty drift),
# do the atoms actually move (O-atom MSD), and are the sampled frames diverse or near-identical
# (per-frame structural descriptor -> PCA scatter + a redundancy ratio).
import sys, os, glob, re, numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from ase.io import read

HARTREE_eV = 27.211386245988

def parse_ener(path):
    rows=[]
    for ln in open(path):
        ln=ln.strip()
        if not ln or ln.startswith("#"): continue
        try: rows.append([float(x) for x in ln.split()])
        except ValueError: pass
    if not rows: return None
    return np.array(rows)   # cols: step time_fs kin_au temp_K pot_au cons_au used_s

def descriptor(o, L, edges):
    """rotation/translation-invariant per-frame fingerprint: histogram of O-O distances."""
    d=o[:,None,:]-o[None,:,:]
    if L: d-=L*np.round(d/L)
    r=np.sqrt((d*d).sum(-1))[np.triu_indices(len(o),1)]
    h,_=np.histogram(r, bins=edges, density=True); return h

def L_from_rho_path(xyz, nmol):
    """Cubic box edge (Angstrom) inferred from a 'rho<NUMBER>' tag in the path (kg/m^3),
    e.g. .../cell_T3000_rho1044/...  ->  rho = 1044 kg/m^3. Returns None if no tag found."""
    m=re.search(r'rho([0-9]+(?:\.[0-9]+)?(?:e[+-]?[0-9]+)?)', xyz)
    if not m: return None
    rho_kg=float(m.group(1)); Mw=18.015e-3; NA=6.02214076e23
    V=nmol*Mw/(rho_kg*NA)            # m^3
    return V**(1.0/3.0)*1e10         # Angstrom

def iact(x):
    """Integrated autocorrelation time (in frame units) of a 1-D series, summed to the
    first non-positive lag. tau=1 => independent frames; tau>>1 => strongly correlated."""
    x=np.asarray(x,float); x=x-x.mean(); n=len(x); v=x.dot(x)/n
    if n<4 or v<=0: return 1.0
    ac=np.correlate(x,x,mode="full")[n-1:]/(v*np.arange(n,0,-1))   # normalized, lag>=0
    tau=1.0
    for k in range(1,n):
        if ac[k]<=0: break
        tau+=2.0*ac[k]
    return float(tau)

def analyze_one(xyz, ener=None, L=None):
    stem=make_stem(xyz)
    frames=read(xyz, index=":")
    nfr=len(frames); nat=len(frames[0])
    sym=np.array(frames[0].get_chemical_symbols()); isO=sym=="O"
    # box edge, in priority order: explicit CLI arg > xyz Lattice header > rho in path name
    L_src = "arg" if L is not None else "none"
    if L is None:
        c=np.array(frames[0].cell)
        if np.any(c): L=float(np.mean(np.diag(c))); L_src="xyz-cell"
    if L is None:
        Lp=L_from_rho_path(xyz, nat//3)              # nat//3 water molecules
        if Lp is not None: L=Lp; L_src="rho-name"
    P=np.array([f.get_positions() for f in frames])     # (nfr,nat,3)
    O=P[:,isO,:]                                         # (nfr,nO,3)

    fig,ax=plt.subplots(2,2,figsize=(12,8)); fig.suptitle(stem)

    # (A) temperature + (B) energies from .ener
    drift_meV=None; Tmean=Tstd=None
    if ener is not None:
        E=parse_ener(ener)
        if E is not None and E.shape[1]>=6:
            t=E[:,1]/1000.0                              # ps
            ax[0,0].plot(t, E[:,3]); ax[0,0].set_xlabel("time (ps)"); ax[0,0].set_ylabel("T (K)")
            half=len(E)//2; Tmean=E[half:,3].mean(); Tstd=E[half:,3].std()
            ax[0,0].set_title(f"temperature  (2nd half: {Tmean:.0f}+/-{Tstd:.0f} K)")
            pot=(E[:,4]-E[0,4])*HARTREE_eV*1000.0/nat    # meV/atom rel.
            cons=(E[:,5]-E[0,5])*HARTREE_eV*1000.0/nat
            ax[0,1].plot(t, pot, label="potential")
            ax[0,1].plot(t, cons, label="conserved qty")
            if len(t)>2:
                slope=np.polyfit(t, cons, 1)[0]          # meV/atom per ps
                drift_meV=slope
                ax[0,1].set_title(f"energy drift  (cons: {slope:+.3f} meV/atom/ps)")
            ax[0,1].set_xlabel("time (ps)"); ax[0,1].set_ylabel("rel. energy (meV/atom)")
            ax[0,1].legend(fontsize=8)
    else:
        ax[0,0].text(.5,.5,"no .ener file",ha="center"); ax[0,1].text(.5,.5,"no .ener file",ha="center")

    # (C) O-atom MSD (unwrapped if box known)
    if L:
        U=np.empty_like(O); U[0]=O[0]
        for k in range(1,nfr):
            dd=O[k]-O[k-1]; dd-=L*np.round(dd/L); U[k]=U[k-1]+dd
    else:
        U=O
    msd=((U-U[0])**2).sum(-1).mean(1)                    # (nfr,)
    fr=np.arange(nfr)
    ax[1,0].plot(fr, msd); ax[1,0].set_xlabel("frame"); ax[1,0].set_ylabel("O MSD (A^2)")
    moved=msd[-1]
    ax[1,0].set_title(f"O MSD: reaches {moved:.2f} A^2 "
                      + ("(atoms diffuse)" if moved>1.0 else "(atoms barely move)"))

    # (D) configurational diversity: descriptor PCA scatter + redundancy ratio
    sub=np.linspace(0,nfr-1,min(nfr,500)).astype(int)    # cap for speed
    rmax=(L/2 if L else 6.0); edges=np.linspace(0,rmax,60)
    D=np.array([descriptor(O[i], L, edges) for i in sub])
    Dc=D-D.mean(0)
    try:
        _,_,Vt=np.linalg.svd(Dc, full_matrices=False); proj=Dc@Vt[:2].T
        sc=ax[1,1].scatter(proj[:,0], proj[:,1], c=sub, cmap="viridis", s=12)
        plt.colorbar(sc, ax=ax[1,1], label="frame")
    except np.linalg.LinAlgError:
        Vt=None; proj=None; ax[1,1].text(.5,.5,"SVD failed",ha="center")
    # scale-free redundancy: effective number of INDEPENDENT frames from the descriptor
    # autocorrelation along the slowest (PC1) mode.  N_eff = n_analysed / tau.
    n_an=len(sub)
    pc1 = proj[:,0] if proj is not None else D[:,0]
    tau = iact(pc1)
    N_eff = n_an/tau
    frac  = N_eff/n_an                                   # fraction of frames that are independent
    ax[1,1].set_xlabel("PC1"); ax[1,1].set_ylabel("PC2")
    ax[1,1].set_title(f"frame diversity  (tau={tau:.1f} fr, N_eff={N_eff:.0f}/{n_an}, {100*frac:.0f}% indep.)")

    plt.tight_layout(rect=[0,0,1,0.97])
    fn=f"aimd_{stem}.png"; plt.savefig(fn, dpi=130); plt.close(fig)

    verdict=("redundant (frames strongly correlated)" if frac<0.05 else
             "well-sampled / diverse" if frac>0.30 else "moderately correlated")
    motion =("atoms barely translate" if moved<0.5 else "atoms diffuse")
    print(f"{stem}: nfr={nfr} nat={nat} L={('%.2f'%L) if L else 'n/a'}({L_src})  "
          f"T={'%.0f+/-%.0f'%(Tmean,Tstd) if Tmean else 'n/a'}  "
          f"cons-drift={('%.3f meV/atom/ps'%drift_meV) if drift_meV is not None else 'n/a'}  "
          f"MSD_end={moved:.2f}A^2 ({motion})  N_eff={N_eff:.0f}/{n_an} -> {verdict}  -> {fn}")
    return dict(stem=stem, nfr=nfr, Tmean=Tmean, Tstd=Tstd, drift=drift_meV,
                msd=moved, tau=tau, N_eff=N_eff, frac=frac, verdict=verdict)

def make_stem(xyz):
    """Choose a unique, informative label: prefer whichever of the filename or parent
    directory carries the cell identity (a 'rho...' or 'T####' tag). Handles both AIMD
    cells (all named aimd-pos-1.xyz in distinct dirs) and pools (pool_T*_rho*.xyz in one dir)."""
    parent=os.path.basename(os.path.normpath(os.path.dirname(xyz)))
    filestem=os.path.splitext(os.path.basename(xyz))[0]
    info=lambda s: ("rho" in s.lower()) or bool(re.search(r"t\d{3,4}", s.lower()))
    if info(filestem): return filestem      # e.g. pool_T1000_rho130
    if info(parent):   return parent        # e.g. cell_T3000_rho1044
    return f"{parent}_{filestem}" if parent else filestem

def find_ener(xyz):
    """Match a CP2K .ener to a trajectory: same project base, else any *-1.ener in the dir."""
    base=xyz.split("-pos-1.xyz")[0]
    en=base+"-1.ener"
    if os.path.exists(en): return en
    d=os.path.dirname(xyz) or "."
    cands=sorted(glob.glob(os.path.join(d,"*-1.ener")) or glob.glob(os.path.join(d,"*.ener")))
    return cands[0] if cands else None

def discover(path, pattern=None):
    """Find trajectory/pool xyz files under <path> (recursively if a directory).
    Tries CP2K trajectories first, then loose '*pos*', then foundation 'pool*' files.
    A custom `pattern` (from --glob) overrides the chain."""
    if os.path.isfile(path):
        return [(path, find_ener(path))]
    if pattern:
        xyzs=sorted(glob.glob(os.path.join(path,"**",pattern), recursive=True))
    else:
        xyzs=[]
        for pat in ("*-pos-1.xyz", "*pos*.xyz", "pool*.xyz"):
            xyzs=sorted(glob.glob(os.path.join(path,"**",pat), recursive=True))
            if xyzs: break
    return [(x, find_ener(x)) for x in xyzs]

def quick_count(xyz):
    """Cheap frame count without loading the whole trajectory with ASE."""
    try:
        with open(xyz) as f:
            nat=int(f.readline().split()[0])
        with open(xyz) as f:
            lines=sum(1 for _ in f)
        return nat, lines//(nat+2)
    except Exception:
        return None, None

def main():
    args=[a for a in sys.argv[1:]]
    do_list = "--list" in args; args=[a for a in args if a!="--list"]
    mx=200; pattern=None
    if "--max" in args:
        i=args.index("--max"); mx=int(args[i+1]); del args[i:i+2]
    if "--glob" in args:
        i=args.index("--glob"); pattern=args[i+1]; del args[i:i+2]
    if not args:
        sys.exit("usage: python analyze_aimd.py <dir-or-xyz> [L_angstrom] [--list] [--max N] [--glob PAT]\n"
                 "  --list : just show which files were found (no analysis)\n"
                 "  --glob : custom filename pattern, e.g. --glob 'pool*.xyz'\n"
                 "  searches recursively; finds *-pos-1.xyz, then *pos*.xyz, then pool*.xyz")
    path=args[0]; L=float(args[1]) if len(args)>1 else None

    jobs=discover(path, pattern)
    if not jobs:
        sys.exit(f"no *-pos-1.xyz (or *pos*.xyz) trajectories found under '{path}'.\n"
                 f"Try a wider root, e.g.  python analyze_aimd.py . --list")

    print(f"found {len(jobs)} trajectorie(s) under '{path}':")
    for xyz,en in jobs:
        nat,nfr=quick_count(xyz)
        print(f"  {xyz}   [{nfr if nfr is not None else '?'} frames, "
              f"{nat if nat is not None else '?'} atoms]   ener={'yes' if en else 'no'}")
    if do_list:
        print("\n(--list only: re-run without --list to analyse)")
        return

    if len(jobs)>mx:
        print(f"\n(analysing first {mx} of {len(jobs)}; raise with --max)")
        jobs=jobs[:mx]
    print()
    summ=[analyze_one(x,e,L) for x,e in jobs]
    if len(summ)>1:
        print("\nsummary (watch for: large cons-drift = unstable; tiny MSD + low % indep = redundant sampling):")
        for s in sorted(summ, key=lambda s:s['frac']):
            print(f"  {s['stem']:<28} MSD={s['msd']:6.2f}A^2  "
                  f"N_eff={s['N_eff']:5.0f} ({100*s['frac']:3.0f}% indep)  {s['verdict']}")

if __name__=="__main__":
    main()
