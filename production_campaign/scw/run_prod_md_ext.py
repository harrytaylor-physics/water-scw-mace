import sys, numpy as np, os
from ase import Atoms, units
from ase.io import write
from ase.md.langevin import Langevin
from ase.md.verlet import VelocityVerlet
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution, Stationary, ZeroRotation
from mace.calculators import MACECalculator

T=float(sys.argv[1]); rho=float(sys.argv[2])
N=int(sys.argv[3]) if len(sys.argv)>3 else 64; Mw=18.015e-3; NA=6.02214076e23
L=(N*Mw/(rho*NA))**(1/3)*1e10
if T>3500.0:
    pass  # dilute skip disabled for extension campaign
if L>28.0:
    pass  # dilute skip disabled for extension campaign
tag=f"T{T:.0f}_rho{rho:g}_05fs"
out=os.path.join(os.environ["DATA"],"scw_prod",tag)
os.makedirs(out, exist_ok=True)

n=int(round(N**(1/3))); assert n**3==N, f"N must be a perfect cube, got {N}"; spacing=L/n; rng=np.random.default_rng(0)
def rand_rot():
    A=rng.standard_normal((3,3)); Q,R=np.linalg.qr(A); Q=Q@np.diag(np.sign(np.diag(R)+1e-12))
    if np.linalg.det(Q)<0: Q[:,0]=-Q[:,0]
    return Q
def water():
    d=0.9572; a=np.deg2rad(104.52); return np.array([[0,0,0],[d,0,0],[d*np.cos(a),d*np.sin(a),0]])
def minok(w,pos):
    if not pos: return True
    P=np.array(pos)
    for x in w:
        dd=P-x; dd-=L*np.round(dd/L)
        if (np.einsum('ij,ij->i',dd,dd)<1.3**2).any(): return False
    return True
pos=[]; sym=[]
for i in range(n):
  for j in range(n):
    for k in range(n):
      c=np.array([(i+.5)*spacing,(j+.5)*spacing,(k+.5)*spacing])
      for _ in range(300):
          w=water()@rand_rot().T; w=w-w[0]+c
          if minok(w,pos): break
      pos.extend(w); sym+=["O","H","H"]
at=Atoms(symbols=sym,positions=pos,cell=[L,L,L],pbc=True)
# FAILFAST: abort in seconds if the cell is not what was requested
assert len(at)==3*N, f"placed {len(at)} atoms, expected {3*N}"
rho_real=(len(at)/3*18.015e-3/6.02214076e23)/(L*1e-10)**3
assert abs(rho_real-rho)/rho<0.01, f"realised density {rho_real:.1f} != requested {rho:.1f}"
print(f"FAILFAST OK: N={len(at)//3} mol, L={L:.2f} A, rho={rho_real:.2f} kg/m3", flush=True)
at.calc=MACECalculator(model_paths="water_scw_v4.model",device="cuda",default_dtype="float32",enable_cueq=True)

MaxwellBoltzmannDistribution(at,temperature_K=T); Stationary(at); ZeroRotation(at)
Langevin(at,0.1*units.fs,temperature_K=T,friction=0.1).run(2000)
Langevin(at,0.5*units.fs,temperature_K=T,friction=0.02).run(20000)     # 10 ps NVT equil
write(f"{out}/prod.xyz",at)
dyn=VelocityVerlet(at,0.5*units.fs)
log=open(f"{out}/thermo.log","w")
def thermo(a=at):
    ek=a.get_kinetic_energy()/len(a); Ti=ek/(1.5*units.kB)
    log.write(f"{dyn.nsteps*0.5:.1f} {Ti:.1f} {a.get_potential_energy()/len(a):.6f}\n"); log.flush()
dyn.attach(thermo, interval=2000)
dyn.attach(lambda a=at: write(f"{out}/prod.xyz",a,append=True), interval=1)   # every 0.5 fs
dyn.run(60000)                                                         # 30 ps NVE
print(f"done {tag}  L={L:.3f} A  frames={60000//1+1}")
