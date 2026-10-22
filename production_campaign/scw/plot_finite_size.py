# plot_finite_size.py -- gather fs/*.npz and show whether properties converge with N.
# usage: python plot_finite_size.py <tag>
import sys, glob, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

tag=sys.argv[1]
def load(which):
    rows=[]
    for f in glob.glob(f"fs/{tag}_{which}_N*.npz"):
        z=np.load(f); rows.append((int(z["N"]),float(z["P_mean"]),float(z["P_se"]),
                                   z["rc"],z["gOO"],float(z["L"]),int(z["natoms"])))
    return sorted(rows, key=lambda r:r[0])

fig,ax=plt.subplots(1,2,figsize=(13,5))
for which,c in (("v4","C0"),("foundation","C1")):
    rows=load(which)
    if not rows: continue
    N =np.array([r[0] for r in rows]); P=np.array([r[1] for r in rows]); SE=np.array([r[2] for r in rows])
    ax[0].errorbar(N,P,yerr=SE,fmt="o-",c=c,capsize=3,label=which)
    # RDF overlay (first peak position/height printed)
    print(f"\n{which} finite-size:")
    for n,pm,se,rc,g,L,na in rows:
        ipk=np.argmax(g); ax[1].plot(rc,g,lw=1.2,label=f"{which} N={n} ({na} at.)" if which=="v4" else None)
        print(f"   N={n:5d} ({na:6d} atoms) L={L:6.2f}A  P={pm:+.3f}+/-{se:.3f} GPa  "
              f"1st peak: g={g[ipk]:.3f} at r={rc[ipk]:.2f}A")
    # convergence relative to largest box
    Pinf=P[-1]
    print(f"   -> P drift vs largest box: " + ", ".join(f"N={n}:{pp-Pinf:+.3f}" for n,pp in zip(N,P)))

ax[0].set_xscale("log"); ax[0].set_xlabel("N molecules (log)"); ax[0].set_ylabel("P (GPa)")
ax[0].set_title(f"{tag}: pressure vs system size (should plateau)"); ax[0].legend(fontsize=8)
ax[1].set_xlabel("r (A)"); ax[1].set_ylabel("g_OO(r)")
ax[1].set_title("O-O RDF vs system size (tail -> 1; larger box reaches larger r)")
ax[1].axhline(1,c="grey",lw=.5,ls=":"); ax[1].legend(fontsize=7)
plt.tight_layout(); plt.savefig(f"fs/{tag}_finite_size.png",dpi=140)
print(f"\n-> wrote fs/{tag}_finite_size.png")
