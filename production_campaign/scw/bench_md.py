# bench_md.py  --  python bench_md.py <v4|foundation> <T_K> <rho_g_cc> <tag>
import sys, os, numpy as np
from ase import Atoms, units
from ase.md.langevin import Langevin
from ase.md.verlet import VelocityVerlet
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution, Stationary, ZeroRotation
from mace.calculators import MACECalculator

which=sys.argv[1]; T=float(sys.argv[2]); rho=float(sys.argv[3]); tag=sys.argv[4]
MODEL={"v4":"water_scw_v4.model","foundation":"mace-mpa-0-medium.model"}[which]
N=64; Mw=18.015e-3; NA=6.02214076e23
L=(N*Mw/((rho*1000)*NA))**(1/3)*1e10                      # rho in g/cc -> kg/m3 -> box A
out=f"bench/{tag}_{which}"; os.makedirs(out,exist_ok=True)
EVA3_to_GPa=160.21766

# --- build 64-water box (same packer as production) ---
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

# --- equilibrate (NVT) then production (NVE) ---
MaxwellBoltzmannDistribution(at,temperature_K=T); Stationary(at); ZeroRotation(at)
Langevin(at,0.1*units.fs,temperature_K=T,friction=0.1).run(2000)
Langevin(at,0.5*units.fs,temperature_K=T,friction=0.02).run(20000)     # 10 ps NVT
dyn=VelocityVerlet(at,0.5*units.fs)
V=[]; POS=[]; PRESS=[]
M=at.get_masses()
def grab(a=at):
    V.append(a.get_velocities().copy()); POS.append(a.get_positions().copy())
    s=a.get_stress(voigt=True)                                          # virial only (eV/A^3)
    P_vir=-(s[0]+s[1]+s[2])/3.0                                         # eV/A^3
    P_kin=(2.0/3.0)*a.get_kinetic_energy()/a.get_volume()              # ideal-gas term (eV/A^3)
    PRESS.append((P_vir+P_kin)*EVA3_to_GPa)                             # TOTAL pressure -> GPa
dyn.attach(grab, interval=2)                                            # 1 fs
dyn.run(40000)                                                          # 20 ps NVE
V=np.array(V); POS=np.array(POS); PRESS=np.array(PRESS)

# --- O-O RDF ---
isO=np.array([s=="O" for s in sym]); rmax=min(L/2,6.0); nb=300
edges=np.linspace(0,rmax,nb+1); rc=0.5*(edges[1:]+edges[:-1]); hist=np.zeros(nb)
for P in POS[::5]:
    o=P[isO]; d=o[:,None,:]-o[None,:,:]; d-=L*np.round(d/L)
    r=np.sqrt((d**2).sum(-1)); r=r[np.triu_indices(len(o),1)]
    hist+=np.histogram(r,bins=edges)[0]
nfr=len(POS[::5]); nO=isO.sum(); numdens=nO/L**3
shell=4*np.pi*rc**2*(edges[1]-edges[0])
gOO=hist/(nfr*0.5*nO*numdens*shell)

# --- VDOS from mass-weighted velocity ACF ---
dt=1e-15                                                                # 1 fs sampling
v=V-V.mean(0)
vacf=np.zeros(len(v))
for i in range(len(M)):
    for k in range(3):
        x=v[:,i,k]; f=np.fft.rfft(x,n=2*len(x)); a=np.fft.irfft(f*np.conj(f))[:len(x)].real
        vacf+=M[i]*a
vacf/=(len(v)-np.arange(len(v)))
win=np.hanning(2*len(vacf))[len(vacf):]; nfft=2**16
dos=np.abs(np.fft.rfft(vacf*win,n=nfft))
wn=np.fft.rfftfreq(nfft,dt)/2.99792458e10

np.savez(f"{out}/bench.npz", rc=rc, gOO=gOO, wn=wn, dos=dos,
         P_GPa=PRESS.mean(), P_std=PRESS.std(), T=T, rho=rho, which=which)
print(f"{tag} {which}: <P>={PRESS.mean():.3f}+/-{PRESS.std():.3f} GPa (total)  L={L:.2f}A")
