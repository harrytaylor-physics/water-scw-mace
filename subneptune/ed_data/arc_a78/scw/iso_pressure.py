# iso_pressure.py -- one converged pressure point for an isotherm sweep.
# usage: python iso_pressure.py <v4|foundation> <T_K> <rho_g_cc> <isotag>
# writes iso/<isotag>_<which>_rho<rho>.npz  with P_mean, P_se, last25 (convergence check)
import sys, os, numpy as np
from ase import units
from ase.md.langevin import Langevin
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution, Stationary, ZeroRotation
from mace.calculators import MACECalculator
import wbox

which=sys.argv[1]; T=float(sys.argv[2]); rho=float(sys.argv[3]); isotag=sys.argv[4]
MODEL={"v4":"water_scw_v4.model","foundation":"mace-mpa-0-medium.model"}[which]
N=64
at,L=wbox.build_box(N, rho, seed=0)
at.calc=MACECalculator(model_paths=MODEL, device="cuda", default_dtype="float32")

# 40 ps NVT equilibration (melt the grid), then 60 ps NVT pressure sampling
MaxwellBoltzmannDistribution(at,temperature_K=T); Stationary(at); ZeroRotation(at)
Langevin(at,0.1*units.fs,temperature_K=T,friction=0.1).run(2000)
Langevin(at,0.5*units.fs,temperature_K=T,friction=0.02).run(80000)     # 40 ps equil
dyn=Langevin(at,0.5*units.fs,temperature_K=T,friction=0.01)
Ps=[]
dyn.attach(lambda a=at: Ps.append(wbox.total_pressure_GPa(a)), interval=10)   # every 5 fs
dyn.run(120000)                                                              # 60 ps
Ps=np.array(Ps)

run=np.cumsum(Ps)/np.arange(1,len(Ps)+1)
P_mean=float(Ps.mean()); P_se=float(wbox.block_se(Ps)); last25=float(run[int(0.75*len(run)):].mean())
os.makedirs("iso",exist_ok=True)
np.savez(f"iso/{isotag}_{which}_rho{rho:.3f}.npz",
         P_mean=P_mean, P_se=P_se, P_std=float(Ps.std()), last25=last25,
         T=T, rho=rho, which=which, L=L)
conv = "OK" if abs(last25-P_mean)<2*P_se else "CHECK(drifting)"
print(f"[iso] {isotag} {which} rho={rho:.3f} T={T:.0f}: P={P_mean:+.3f} SE={P_se:.3f} "
      f"last25={last25:+.3f} -> {conv}  L={L:.2f}A")
