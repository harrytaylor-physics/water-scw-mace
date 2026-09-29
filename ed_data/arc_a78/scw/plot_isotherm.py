# plot_isotherm.py -- gather iso/*.npz for one isotherm and plot P(rho) vs AQUA/IAPWS,
# with the bulk modulus K_T = rho dP/drho from a spline through the model points.
# usage: python plot_isotherm.py <isotag> <T_K>
import sys, glob, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline, griddata

isotag=sys.argv[1]; T=float(sys.argv[2])
AQUA="/home/shug7909/scw/aqua_eos_rhot_v1_0.dat"

def load(which):
    d={}
    for f in glob.glob(f"iso/{isotag}_{which}_rho*.npz"):
        z=np.load(f); d[float(z["rho"])]=(float(z["P_mean"]),float(z["P_se"]),float(z["last25"]))
    r=np.array(sorted(d)); 
    if len(r)==0: return r,r,r,r
    P =np.array([d[x][0] for x in r]); SE=np.array([d[x][1] for x in r]); L25=np.array([d[x][2] for x in r])
    return r,P,SE,L25

# --- AQUA isotherm by log-log interpolation onto a fine rho grid ---
def aqua_curve(rgrid):
    try:
        A=np.loadtxt(AQUA, skiprows=24, usecols=(0,1,2))           # rho[kg/m3],T[K],P[Pa]
        pts=np.column_stack([np.log10(A[:,0]), np.log10(A[:,1])])
        val=A[:,2]*1e-9                                            # GPa
        q=np.column_stack([np.log10(rgrid*1000), np.full_like(rgrid,np.log10(T))])
        return griddata(pts,val,q,method="linear")
    except Exception as e:
        print("AQUA curve skipped:",e); return np.full_like(rgrid,np.nan)

def iapws_pts(rgrid):
    try:
        from iapws import IAPWS95
        out=[]
        for r in rgrid:
            try: out.append(IAPWS95(T=T,rho=r*1000).P*1e-3)
            except Exception: out.append(np.nan)
        return np.array(out)
    except Exception as e:
        print("IAPWS skipped:",e); return np.full_like(rgrid,np.nan)

fig,ax=plt.subplots(1,2,figsize=(13,5))
rgrid=np.linspace(0.05,1.4,200)
ax[0].plot(rgrid,aqua_curve(rgrid),"-",c="k",lw=1.5,label="AQUA")
ip=iapws_pts(rgrid); 
if np.isfinite(ip).any(): ax[0].plot(rgrid,ip,"--",c="grey",lw=1.2,label="IAPWS-95")

colors={"v4":"C0","foundation":"C1"}
for which in ("v4","foundation"):
    r,P,SE,L25=load(which)
    if len(r)==0: continue
    ax[0].errorbar(r,P,yerr=SE,fmt="o-",c=colors[which],capsize=3,label=which)
    drift=np.abs(L25-P)>2*SE
    if drift.any(): ax[0].plot(r[drift],P[drift],"x",c="red",ms=11,mew=2,label=f"{which}: not converged")
    # bulk modulus K_T = rho dP/drho via smoothing spline
    if len(r)>=4:
        spl=UnivariateSpline(r,P,k=3,s=len(r)*np.mean(SE)**2)
        K=r*spl.derivative()(r)
        ax[1].plot(r,K,"o-",c=colors[which],label=f"{which}")
        print(f"\n{which} bulk modulus K_T = rho dP/drho (GPa):")
        for x,k in zip(r,K): print(f"   rho={x:.3f}  P={spl(x):+.3f}  K_T={k:+.3f}")

ax[0].axhline(0,c="k",lw=.5); ax[0].set_xlabel("rho (g/cc)"); ax[0].set_ylabel("P (GPa)")
ax[0].set_title(f"{isotag}: {T:.0f} K isotherm  (P = virial + ideal-gas)"); ax[0].legend(fontsize=8)
ax[1].set_xlabel("rho (g/cc)"); ax[1].set_ylabel("K_T (GPa)")
ax[1].set_title("bulk modulus  K_T = rho (dP/drho)_T"); ax[1].legend(fontsize=8)
plt.tight_layout(); plt.savefig(f"iso/{isotag}_isotherm.png",dpi=140)
print(f"\n-> wrote iso/{isotag}_isotherm.png")
