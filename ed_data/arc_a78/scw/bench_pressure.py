# bench_pressure.py  --  python bench_pressure.py <v4|foundation> <T_K> <rho_g_cc> <tag>
import sys, numpy as np
from ase import Atoms, units
from ase.md.langevin import Langevin
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution, Stationary, ZeroRotation
from mace.calculators import MACECalculator

which=sys.argv[1]; T=float(sys.argv[2]); rho=float(sys.argv[3]); tag=sys.argv[4]
MODEL={"v4":"water_scw_v4.model","foundation":"mace-mpa-0-medium.model"}[which]
N=64; Mw=18.015e-3; NA=6.02214076e23
L=(N*Mw/((rho*1000)*NA))**(1/3)*1e10

# --- build box (same packer) ---
n=4; spacing=L/n; rng=np.random.default_rng(0)
def rot():
    A=rng.standard_normal((3,3)); Q,R=np.linalg.qr(A); Q=Q@np.diag(np.sign(np.diag(R)+1e-12))
    if np.linalg.det(Q)<0: Q[:,0]=-Q[:,0]
    return Q
def w(): d=0.9572; a=np.deg2rad(104.52); return np.array([[0,0,0],[d,0,0],[d*np.cos(a),d*np.sin(a),0]])
def ok(m,p):
    if not p: return True
    P=np.array(p)
    for x in m:
        dd=P-x; dd-=L*np.round(dd/L)
        if (np.einsum('ij,ij->i',dd,dd)<1.2**2).any(): return False
    return True
pos=[]; sym=[]
for i in range(n):
  for j in range(n):
    for k in range(n):
      c=np.array([(i+.5)*spacing,(j+.5)*spacing,(k+.5)*spacing])
      for _ in range(400):
          m=w()@rot().T; m=m-m[0]+c
          if ok(m,pos): break
      pos.extend(m); sym+=["O","H","H"]
at=Atoms(symbols=sym,positions=pos,cell=[L,L,L],pbc=True)
at.calc=MACECalculator(model_paths=MODEL,device="cuda",default_dtype="float32")

EVA3_to_GPa=160.21766
def P_now(a=at):
    # total pressure = virial (potential) + ideal-gas (kinetic) contributions
    s=a.get_stress(voigt=True)                              # virial only (eV/A^3)
    P_vir=-(s[0]+s[1]+s[2])/3.0                             # eV/A^3
    P_kin=(2.0/3.0)*a.get_kinetic_energy()/a.get_volume()  # ideal-gas term (eV/A^3)
    return (P_vir+P_kin)*EVA3_to_GPa

# --- LONG equilibration: 40 ps, melt the grid lattice ---
MaxwellBoltzmannDistribution(at,temperature_K=T); Stationary(at); ZeroRotation(at)
Langevin(at,0.1*units.fs,temperature_K=T,friction=0.1).run(2000)
Langevin(at,0.5*units.fs,temperature_K=T,friction=0.02).run(80000)     # 40 ps NVT equil
print(f"after equil: instantaneous P = {P_now():.3f} GPa")

# --- production NVT, sample pressure, keep running average ---
dyn=Langevin(at,0.5*units.fs,temperature_K=T,friction=0.01)
Ps=[]
def grab(a=at): Ps.append(P_now(a))
dyn.attach(grab, interval=10)              # every 5 fs
dyn.run(120000)                            # 60 ps
Ps=np.array(Ps)

# block-averaged mean + error, and a convergence trace
run=np.cumsum(Ps)/np.arange(1,len(Ps)+1)
nb=10; blk=np.array_split(Ps,nb); bm=np.array([b.mean() for b in blk])
print(f"\n{tag} {which}  T={T:.0f} rho={rho:.2f}")
print(f"  <P> = {Ps.mean():.3f} GPa (total)   block-SE = {bm.std()/np.sqrt(nb):.3f}")
print(f"  running avg last 25%: {run[int(0.75*len(run)):].mean():.3f} GPa")
print(f"  (sign check: ASE stress convention -> positive P means compression)")
np.save(f"bench/P_{tag}_{which}.npy", Ps)
