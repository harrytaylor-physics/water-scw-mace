# convergence/

Runnable convergence tests for every stage of the pipeline — and the
project's convergence diary. Each subfolder has: a README stating what to
converge, the criterion, and the reference results; and core scripts that run the check on different systems.

Scripts are cluster-agnostic: they generate inputs / analyse outputs /
run analysis locally, and never submit jobs. How you run the generated
inputs (SLURM, PBS, laptop) can be customised. Dependencies: numpy, ase,
matplotlib; CP2K for 01; mace-torch + a GPU helps for 02/03.

| Stage | Converge | Criterion | Project reference |
|---|---|---|---|
| 01 electronic structure | plane-wave cutoff, rel cutoff, basis | <=1 meV/atom, <=10 meV/A vs high-cutoff ref | 600/60 Ry, DZVP-MOLOPT-SR-GTH |
| 02 training data | dataset size, seed variance | validation force-RMSE plateau | 1193 cfgs (v4), 1080 (lda_v1) |
| 03 MD sampling | traj length, box size, timestep | block-averaged means stable; NVE drift <=2 meV/atom/ps | 60 ps, 64 mol (flat to 32x vol), 0.5 fs |
| 04 spectra | ACF window, traj length, smoothing | band positions/intensities stable | maxlag 8000, nfft 2^17, 20 cm^-1 |
| 05 post-processing | Rosseland window/floor, KK cutoff | means stable <2%; KK drift <5% | 100-5000 cm^-1, 1e-6 floor |

Diary entries (what the project actually did) live in each README under
"Reference results" — including the honest ones, e.g. the LDA campaign
inherited 600/60 from the hybrid benchmarks rather than re-converging
(logged decision), and 01 is the tool to verify that inheritance.
