import numpy as np, os
from ase import Atoms, units
from ase.io import write
from ase.md.langevin import Langevin
from ase.md.verlet import VelocityVerlet
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution, Stationary, ZeroRotation
from mace.calculators import MACECalculator

# (T K, rho kg/m3): dilute->dense spread inside the v2 force-field domain.
# 1500/10 ~ your 1500K reference point; 3000/100 gets added once v4 passes.
CELLS=[(1500,10.0),(800,61.0),(1000,130.0),(2250,176.0),(2250,1256.0)]
NSEG=25; SEGLEN=8; GAP_FS=200          # 25 segments x 8 frames @1 fs -> 200 frames/cell

calc=MACECalculator(model_paths="water_scw_v2.model",device="cuda",default_dtype="float32")
Mw=18.015e-3; NA=6.02214076e23; rng=np.random.default_rng(1)
def rand_rot():
    A=rng.standard_normal((3,3)); Q,R=np.linalg.qr(A); Q=Q@np.diag(np.sign(np.diag(R)+1e-12))
    if np.linalg.det(Q)<0: Q[:,0]=-Q[:,0]
    return Q
def water():
    d=0.9572; a=np.deg2rad(104.52); return np.array([[0,0,0],[d,0,0],[d*np.cos(a),d*np.sin(a),0]])
os.makedirs("snippets",exist_ok=True)
for T,rho in CELLS:
    N=64; L=(N*Mw/(rho*NA))**(1/3)*1e10
    n=4; sp=L/n; pos=[]; sym=[]
    def minok(w):
        if not pos: return True
        P=np.array(pos)
        for x in w:
            dd=P-x; dd-=L*np.round(dd/L)
            if (np.einsum('ij,ij->i',dd,dd)<1.3**2).any(): return False
        return True
    for i in range(n):
      for j in range(n):
        for k in range(n):
          c=np.array([(i+.5)*sp,(j+.5)*sp,(k+.5)*sp])
          for _ in range(300):
              w=water()@rand_rot().T; w=w-w[0]+c
              if minok(w): break
          pos.extend(w); sym+=["O","H","H"]
    at=Atoms(symbols=sym,positions=pos,cell=[L,L,L],pbc=True); at.calc=calc
    MaxwellBoltzmannDistribution(at,temperature_K=T); Stationary(at); ZeroRotation(at)
    Langevin(at,0.1*units.fs,temperature_K=T,friction=0.1).run(2000)
    Langevin(at,0.5*units.fs,temperature_K=T,friction=0.02).run(6000)      # 3 ps equil
    out=f"snippets/snip_T{T}_rho{rho:g}.xyz"
    dyn=VelocityVerlet(at,0.5*units.fs)
    for s in range(NSEG):
        for f in range(SEGLEN):
            write(out,at,append=True)
            dyn.run(2)                                  # 1 fs between frames in a segment
        dyn.run(int(GAP_FS/0.5))                        # decorrelate between segments
    print(f"{out}: {NSEG*SEGLEN} frames  L={L:.2f} A", flush=True)
print("snippets done")
