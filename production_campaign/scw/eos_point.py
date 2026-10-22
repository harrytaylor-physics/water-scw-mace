# eos_point.py -- one EoS pressure point for an isotherm sweep.
#   python eos_point.py <v4|foundation> <T_K> <rho_g_cc> <tag>
import sys, os, numpy as np
from ase import units
from ase.md.langevin import Langevin
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution, Stationary, ZeroRotation
import bench_common as bc

which=sys.argv[1]; T=float(sys.argv[2]); rho=float(sys.argv[3]); tag=sys.argv[4]
out=f"eos/{tag}"; os.makedirs(out, exist_ok=True)

at,L = bc.build_box(64, rho, seed=0)
at.calc = bc.make_calc(which)

MaxwellBoltzmannDistribution(at, temperature_K=T); Stationary(at); ZeroRotation(at)
Langevin(at, 0.1*units.fs, temperature_K=T, friction=0.1).run(2000)
Langevin(at, 0.5*units.fs, temperature_K=T, friction=0.02).run(80000)     # 40 ps equil
dyn = Langevin(at, 0.5*units.fs, temperature_K=T, friction=0.01)
Ps=[]
dyn.attach(lambda a=at: Ps.append(bc.pressure_GPa(a)), interval=10)        # every 5 fs
dyn.run(60000)                                                            # 60 ps
Ps=np.array(Ps)

run = np.cumsum(Ps)/np.arange(1, len(Ps)+1)
P   = float(Ps.mean()); SE = bc.block_se(Ps, 10)
last25 = float(run[int(0.75*len(run)):].mean())
drift  = abs(last25 - P)
converged = bool(drift < 2*SE)

np.savez(f"{out}/{which}_rho{rho:.3f}.npz",
         which=which, rho=rho, T=T, L=L, P=P, SE=SE,
         last25=last25, drift=drift, converged=converged, Ps=Ps)
print(f"{tag} {which} rho={rho:.3f} T={T:.0f}: P={P:+.3f} +/- {SE:.3f} GPa  "
      f"last25%={last25:+.3f} drift={drift:.3f} converged={converged}")
