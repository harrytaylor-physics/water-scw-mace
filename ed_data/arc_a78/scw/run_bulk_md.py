import numpy as np, csv, os
from ase import Atoms, units
from ase.io import write
from ase.md.langevin import Langevin
from ase.md.verlet import VelocityVerlet
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution, Stationary, ZeroRotation
from mace.calculators import MACECalculator

rows=[r for r in csv.DictReader(open("phase0_grid.csv")) if r["regime"]=="bulk_finetune"]
cell=max(rows, key=lambda r: float(r["rho_kg_m3"]))      # densest bulk cell: 2250 K / 1256
T=float(cell["T_K"]); rho=float(cell["rho_kg_m3"])
N=64; Mw=18.015e-3; NA=6.02214076e23
L=(N*Mw/(rho*NA))**(1/3)*1e10
n=4; spacing=L/n
rng=np.random.default_rng(0)
def rand_rot():
    A=rng.standard_normal((3,3)); Q,R=np.linalg.qr(A); Q=Q@np.diag(np.sign(np.diag(R)+1e-12))
    if np.linalg.det(Q)<0: Q[:,0]=-Q[:,0]
    return Q
def water():
    d=0.9572; a=np.deg2rad(104.52)
    return np.array([[0,0,0],[d,0,0],[d*np.cos(a),d*np.sin(a),0]])
def minok(w,pos,L):
    if not pos: return True
    P=np.array(pos)
    for a in w:
        d=P-a; d-=L*np.round(d/L)
        if (np.einsum('ij,ij->i',d,d)<1.3**2).any(): return False
    return True
pos=[]; sym=[]
for i in range(n):
  for j in range(n):
    for k in range(n):
      c=np.array([(i+.5)*spacing,(j+.5)*spacing,(k+.5)*spacing])
      for _ in range(200):
          w=water()@rand_rot().T; w=w-w[0]+c
          if minok(w,pos,L): break
      pos.extend(w); sym+=["O","H","H"]
atoms=Atoms(symbols=sym,positions=pos,cell=[L,L,L],pbc=True)
atoms.calc=MACECalculator(model_paths="water_scw_v2.model",device="cuda",default_dtype="float32")
os.makedirs("md/bulk_dense",exist_ok=True)
MaxwellBoltzmannDistribution(atoms,temperature_K=T); Stationary(atoms); ZeroRotation(atoms)
Langevin(atoms,0.1*units.fs,temperature_K=T,friction=0.1).run(2000)     # 0.2 ps soft start
Langevin(atoms,0.5*units.fs,temperature_K=T,friction=0.02).run(20000)   # 10 ps equilibrate
write("md/bulk_dense/prod.xyz",atoms)
prod=VelocityVerlet(atoms,0.5*units.fs)
prod.attach(lambda a=atoms: write("md/bulk_dense/prod.xyz",a,append=True), interval=4)  # 2 fs
prod.run(60000)                                                          # 30 ps
print("done  T=%g rho=%g L=%.3f A"%(T,rho,L))
