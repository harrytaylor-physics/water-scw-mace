# 03 — MD sampling convergence

Three independent convergence paths:

**traj_convergence.py** — is the trajectory long enough? Running averages
of thermo-observables with block-averaged error bars; equilibration-time
estimate. Criterion: running mean stable within its block error over the
final third. Input: any whitespace thermo file (step time T E ...).

**timestep_drift.py** — is dt small enough? Short NVE runs at several dt,
energy-conservation drift per atom. Criterion: <=2 meV/atom/ps at your
production dt. Runs MACE MD directly (GPU advised).

**finite_size.py** — is the box big enough? MACE MD at several N, compares
mean pressure and g_OO. Criterion: pressure within block error and RDF
first peak within ~2% across sizes.

**Reference results:** production 60 ps @ 0.5 fs, 64
molecules. Finite-size flat (pressure/RDF/S(k)) over a 32x volume range,
64 -> 2048 molecules; 6400 OOM'd on 44 GB GPUs. NVE drift gate 2
meV/atom/ps (GF.1). Timestep 0.5 fs vs coarser cross-checked in the
_05fs production redo. Pressure includes the ideal-gas kinetic term; the scripts here include it.
