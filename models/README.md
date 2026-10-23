# models/

Trained model files (git-lfs). See `architecture.svg` for how a MACE
potential works and where the repo's validations attach to it.

![MACE architecture](architecture.svg)

## Inventory & provenance

| File | What | Labels (functional) | Training set | Status |
|---|---|---|---|---|
| water_scw_v4.model | production potential | revPBE0-D3 (CP2K) | train_v4.xyz, 1193 cfgs | validated; apply <= 3000 K. Stage-one model (epoch 109), sha256 dc4a8977...; held-out 57.8 meV/A, 20.8 meV/atom ([mlip_heldout_errors.json](mlip_heldout_errors.json)) |
| water_scw_v4_compiled.model | TorchScript export of v4 | — | — | for LAMMPS / no-Python inference |
| mlip_heldout_errors.json | v4 error tables | — | — | both tables from the training log (stage one = production; stage two SWA, not used), and the A94 reproduction |
| water_dipole_v2.model | dipole model | MLWF revPBE0-D3 | dipole label set | production; behind opacity surface |
| mace-mpa-0-medium.model | foundation model (upstream) | MPA-0 | — | fine-tune parent; never used raw for production |
| water_lda_v1.model | LDA potential | LDA/PADE | train_lda_v1.xyz, 1080 cfgs | in training (functional-sensitivity campaign) |

Provenance rule: every model added here gets a row with its label level and
training set, and the training command lives in the campaign scripts.
Never mix label levels within one training set.

## Production model provenance

`water_scw_v4.model`, sha256
`dc4a8977d44e58cf959901b39a594abd3daf69fb09e5622bbe7e8d3e74b7f06c`.
- **Origin:** written by MACE's stage-one pass at 2026-06-11 22:08:48, from
  epoch 109 ([`logs/train_v4_7924966.out`](../production_campaign/arc_scripts/logs/train_v4_7924966.out),
  lines 121-122). The internal archive name is `water_scw_v4`.
- **Identity:** identical on 95/95 tensors to `checkpoints/water_scw_v4_run-1.model`,
  and on 64/95 to `water_scw_v4_stagetwo.model` (sha256
  `380e20a8595d456ead18794e9410a13843271637fc68adaba92bfe67e9cf92f8`).
- **Loaded by:** `run_prod_md.py:44`, `run_prod_md_05fs.py:44` and
  `run_prod_md_ext.py:49` (in [`production_campaign/scw/`](../production_campaign/scw/)).
- **On ARC:** the only files of that name (`~/scw`, `~/water-scw-mace/models`)
  are byte-identical.

`water_scw_v4_stagetwo.model` was not used for production.

## Quickstart (ASE)

All group workflows use the models through ASE. Three common patterns:

**1. Load and evaluate a structure**
```python
from mace.calculators import MACECalculator
from ase.io import read
calc = MACECalculator(model_paths="water_scw_v4.model", device="cuda")  # or "cpu"
atoms = read("my_config.xyz")          # needs cell + pbc set
atoms.calc = calc
print(atoms.get_potential_energy())    # eV
print(atoms.get_forces())              # eV/A
```

**2. Run MD**
```python
from ase import units
from ase.md.bussi import Bussi          # CSVR thermostat (ASE >= 3.23)
from ase.md.verlet import VelocityVerlet
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution

MaxwellBoltzmannDistribution(atoms, temperature_K=1500)
equil = Bussi(atoms, timestep=0.5*units.fs, temperature_K=1500, taut=100*units.fs)
equil.run(10000)                        # 5 ps CSVR equilibration
prod = VelocityVerlet(atoms, timestep=0.5*units.fs)
prod.run(120000)                        # 60 ps NVE production
```
Use `default_dtype="float32"` in the calculator for ~2x faster MD; float64
for parity/benchmark evaluation. Pressure must include the ideal-gas
kinetic term: P = -tr(stress)/3 + N k_B T / V.

**3. An EOS point** — see `benchmarks/run_benchmarks.py` (the `eos_point`
function is a complete worked example).

You can also run the compiled model through LAMMPS (pair_style mace) if
you want — not documented here.

## Why ASE and not LAMMPS?

(1) The whole pipeline is Python — 
selection, dipole evaluation, spectra, analysis all share one process and
one numpy world; ASE keeps MD in the same stack with no export step.
(2) The native torch model is always exactly what training produced;
the compiled/LAMMPS route adds a conversion that can silently lag the
canonical weights. (3) At our system sizes (64-2048 molecules) the GPU
model evaluation dominates the cost and ASE's Python MD loop overhead is
negligible - LAMMPS's strengths (domain decomposition, multi-GPU scaling,
huge N) only pay off far above our sizes. If you need million-atom or
multi-node MD, use the compiled model in LAMMPS.

## Conditions of use

T <= 3000 K and densities within the trained envelope (~100-2230 kg/m^3).
No ices, no reactive >3000 K regime (at a certain planetary depth the envelope is sufficiently optical thick and the opacity falls out of the Schwarzschild's radiative transfer equation) . Full guidance: ../TESTING.md.
