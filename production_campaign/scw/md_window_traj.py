import numpy as np
from ase.io import read, write
from ase import units
from ase.md.langevin import Langevin
from ase.md.verlet import VelocityVerlet
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution, Stationary
from mace.calculators import MACECalculator
T = 1000.0; NFRAME = 501; STRIDE = 4          # 4 x 0.5 fs = 2 fs between frames
at = read("snippets/snip_T1000_rho130.xyz", index=-1)
assert len(at) == 192, f"expected 192 atoms, got {len(at)}"
at.calc = MACECalculator(model_paths="water_scw_v4.model", device="cuda", default_dtype="float32")
MaxwellBoltzmannDistribution(at, temperature_K=T, rng=np.random.default_rng(1)); Stationary(at)
Langevin(at, 0.5*units.fs, temperature_K=T, friction=0.02, rng=np.random.default_rng(2)).run(4000)  # 2 ps
dyn = VelocityVerlet(at, 0.5*units.fs); out = "window_traj_T1000_rho130.xyz"
for i in range(NFRAME):
    if i > 0: dyn.run(STRIDE)
    at.info["step"] = i*STRIDE; at.info["time_fs"] = i*STRIDE*0.5
    write(out, at, append=(i > 0))
    if i % 100 == 0: print(f"frame {i}  t={i*2} fs  T={at.get_temperature():.0f} K", flush=True)
print("done")
