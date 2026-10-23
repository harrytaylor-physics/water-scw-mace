# finite_size.py -- pressure + O-O RDF at fixed (T, rho) for a given system size.
# usage: python finite_size.py <v4|foundation> <T_K> <rho_g_cc> <N_mol> <tag>
# writes fs/<tag>_<which>_N<Nmol>.npz  with N, P_mean, P_se, rc, gOO, L
import sys, os, numpy as np
from ase import units
from ase.md.langevin import Langevin
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution, Stationary, ZeroRotation
from mace.calculators import MACECalculator
from scipy.spatial import cKDTree
import wbox

which=sys.argv[1]; T=float(sys.argv[2]); rho=float(sys.argv[3]); Nmol=int(sys.argv[4]); tag=sys.argv[5]
MODEL={"v4":"water_scw_v4.model","foundation":"mace-mpa-0-medium.model"}[which]
at,L=wbox.build_box(Nmol, rho, seed=0)
at.calc=MACECalculator(model_paths=MODEL, device="cuda", default_dtype="float32")

# cost ~ N * steps; pressure & RDF converge fast, so shorten the long runs for big boxes
if   Nmol<=128:  n_eq,n_pr = 80000,120000     # 40+60 ps (matches iso)
elif Nmol<=1000: n_eq,n_pr = 20000, 40000     # 10+20 ps
else:            n_eq,n_pr =  8000, 16000     #  4+ 8 ps  (big boxes self-average over atoms)

MaxwellBoltzmannDistribution(at,temperature_K=T); Stationary(at); ZeroRotation(at)
Langevin(at,0.1*units.fs,temperature_K=T,friction=0.1).run(2000)
Langevin(at,0.5*units.fs,temperature_K=T,friction=0.02).run(n_eq)
dyn=Langevin(at,0.5*units.fs,temperature_K=T,friction=0.01)

# online O-O RDF accumulation (don't store positions -> memory-safe for big N)
isO=np.array([s=="O" for s in at.get_chemical_symbols()])
nO=int(isO.sum()); rmax=min(L/2.0, 8.0); nbins=400
edges=np.linspace(0,rmax,nbins+1); rc=0.5*(edges[1:]+edges[:-1])
shell=4*np.pi*rc**2*(edges[1]-edges[0]); numdens=nO/L**3
hist=np.zeros(nbins); nfr=[0]; Ps=[]

def grab_P(a=at): Ps.append(wbox.total_pressure_GPa(a))
def grab_RDF(a=at):
    o=(a.get_positions()[isO]) % L
    tree=cKDTree(o, boxsize=L)
    pairs=tree.query_pairs(rmax, output_type='ndarray')
    if len(pairs):
        d=o[pairs[:,0]]-o[pairs[:,1]]; d-=L*np.round(d/L)
        r=np.sqrt((d*d).sum(1))
        hist[:]+=np.histogram(r,bins=edges)[0]
    nfr[0]+=1

dyn.attach(grab_P,  interval=10)    # 5 fs
dyn.attach(grab_RDF, interval=100)  # 50 fs
dyn.run(n_pr)

Ps=np.array(Ps)
gOO=hist/(nfr[0]*0.5*nO*numdens*shell)
P_mean=float(Ps.mean()); P_se=float(wbox.block_se(Ps))
os.makedirs("fs",exist_ok=True)
np.savez(f"fs/{tag}_{which}_N{Nmol}.npz",
         N=Nmol, natoms=len(at), P_mean=P_mean, P_se=P_se, P_std=float(Ps.std()),
         rc=rc, gOO=gOO, L=L, T=T, rho=rho, which=which)
print(f"[fs] {tag} {which} N={Nmol} ({len(at)} atoms) rho={rho:.2f} T={T:.0f}: "
      f"P={P_mean:+.3f} SE={P_se:.3f}  L={L:.2f}A  rmax={rmax:.2f}A")
