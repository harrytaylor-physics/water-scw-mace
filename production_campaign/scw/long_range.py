# long_range.py -- probe for long-range correlations that a small (64-mol) box cannot hold.
#   python long_range.py <v4|foundation> <T_K> <rho_g_cc> <nmol> <tag>
# Runs a large-box MD and computes (1) the O-O RDF out to the box half-length (far beyond
# the ~7 A a 64-mol box allows) and (2) the structure factor S(k) into the low-k regime.
# Flat g(r)->1 well before L/2 and a flat, finite S(k) at low k  ==>  no long-range structure
# the small box was missing.  A low-k upturn in S(k) or persistent g(r) oscillations = there is.
import sys, os, numpy as np
from ase import units
from ase.md.langevin import Langevin
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution, Stationary, ZeroRotation
import bench_common as bc

which=sys.argv[1]; T=float(sys.argv[2]); rho=float(sys.argv[3]); nmol=int(sys.argv[4]); tag=sys.argv[5]
out=f"lr/{tag}"; os.makedirs(out, exist_ok=True)

# run length shrinks with N (large boxes self-average), same schedule as finite_size.py
if   nmol<=64:  equil,prod=80000,120000
elif nmol<=640: equil,prod=20000,24000
else:           equil,prod=8000,8000

at,L=bc.build_box(nmol, rho, seed=0); at.calc=bc.make_calc(which)
isO=np.array([s=="O" for s in at.get_chemical_symbols()])
MaxwellBoltzmannDistribution(at,temperature_K=T); Stationary(at); ZeroRotation(at)
Langevin(at,0.1*units.fs,temperature_K=T,friction=0.1).run(2000)
Langevin(at,0.5*units.fs,temperature_K=T,friction=0.02).run(equil)
dyn=Langevin(at,0.5*units.fs,temperature_K=T,friction=0.01)
framesO=[]
dyn.attach(lambda a=at: framesO.append(a.get_positions()[isO].copy()), interval=100)  # every 50 fs
dyn.run(prod)

nO=isO.sum(); numdens=nO/L**3
rmax=L/2.0                                    # the whole point: reach out to L/2, not 6 A
rc,gOO=bc.rdf_OO(framesO, L, rmax=rmax, nbins=400)

# ---- structure factor S(k) from the RDF (radial Fourier transform) ----
# S(k) = 1 + 4*pi*rho * integral_0^inf [g(r)-1] r sin(kr)/k dr   (isotropic liquid)
# window g(r)-1 with a smooth cutoff at rmax to suppress truncation ripple.
kmin=2*np.pi/L                                # smallest wavevector the box can represent
kk=np.linspace(kmin, 8.0, 300)                # 1/A ; low-k end is the long-range probe
dr=rc[1]-rc[0]; h=gOO-1.0
win=0.5*(1+np.cos(np.pi*rc/rmax))             # Hann taper -> 0 at rmax
Sk=1.0 + 4*np.pi*numdens*np.array([np.sum(h*win*rc*np.sin(k*rc)/k)*dr for k in kk])

np.savez(f"{out}/{which}_n{nmol}.npz", which=which,nmol=nmol,rho=rho,T=T,L=L,
         rc=rc,gOO=gOO,kk=kk,Sk=Sk,numdens=numdens)

# correlation length: largest r where |g(r)-1| still exceeds a small threshold
tail=np.where(np.abs(h)>0.02)[0]
xi = rc[tail[-1]] if len(tail) else 0.0
print(f"{tag} {which} N={nmol} L={L:.1f}A  RDF reach={rmax:.1f}A")
print(f"  g(r) last exceeds |g-1|>0.02 at r = {xi:.1f} A   (box half-length = {rmax:.1f} A)")
print(f"  -> {'structure extends near/beyond L/2: box may be truncating it' if xi>0.8*rmax else 'g(r) flat well inside L/2: no long-range structure missed'}")
print(f"  S(k) at lowest k ({kmin:.2f} 1/A) = {Sk[0]:.3f}   (a low-k upturn would signal long-range correlation)")
print(f"  -> wrote {out}/{which}_n{nmol}.npz")
