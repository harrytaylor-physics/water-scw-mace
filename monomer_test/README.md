# Isolated-molecule dipole test (A64)

Does `water_dipole_v3` return the right dipole **derivative** for a gas-phase
water molecule?

The dilute-limit band diagnostic (A61) found the first-principles spectrum
carrying **2.4–3.9×** the line list's stretch-band intensity, **0.35–0.50×** its
rotational-band intensity, and a stretch/bend balance **3.3–3.6×** the line
list's — the same at 1,000, 1,500 and 2,000 K. A ratio that is mode-dependent but
temperature-independent is what a wrong dipole derivative looks like, and
`water_dipole_v3` was trained almost entirely on dense configurations, so it may
return a liquid-like dμ/dr when a bond stretches in the gas. Band intensity goes
as |dμ/dQ|², so a stretch derivative too large by ~1.8× would produce the
observed excess on its own.

The isolated molecule settles it: 13 geometries, the model's dipole against
revPBE0/MLWF at the *same* geometries, and against experiment (μ = 1.855 D).

## Files

| file | runs on | what it does |
|---|---|---|
| `make_geoms.py` | anywhere | writes `g00..g12.xyz` (3 atoms, 20 Å box) + `geoms.json` |
| `eval_model.py` | **compute node** | `water_dipole_v3` dipoles → `model_dipoles.json` (Debye) |
| `make_cp2k_inputs.py` | login node | `g*/dip.inp` + `g*/coord.xyz` from your MLWF template |
| `harvest_dft.py` | login node | CP2K Wannier outputs → `dft_dipoles.json` (Debye) |
| `compare.py` | login node | derivatives, ratios, verdict → `compare.json` |

`g00` is the equilibrium (r_OH = 0.9572 Å, HOH = 104.52°). `g01`–`g06` are ±
pairs for symmetric stretch (±0.01 Å, both bonds), asymmetric stretch (±0.01 Å,
opposed) and bend (±1.0°). `g07`–`g12` repeat those three pairs about a stretched
reference, **r_OH = 1.00 Å** — the length the bonds sample at 1,000 K, so an
error that is anharmonic (right at equilibrium, wrong where the spectrum lives)
still shows. Orientation is identical across all thirteen, so the dipole
*vectors* compare directly, not only their magnitudes.

## Running it on ARC, in order

**1. Model dipoles — compute node.** `torch` does not import on a login node.

```bash
cat > monomer_eval.sh <<'EOF'
#!/bin/bash
#SBATCH --job-name=monomer-dipole
#SBATCH --partition=short
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --export=NONE
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
module purge
module load Anaconda3/2025.06-1
cd "$SLURM_SUBMIT_DIR"
python eval_model.py --model ~/scw/water_dipole_v3.model --device cpu
EOF
sbatch --export=NONE monomer_eval.sh          # bash -l is implied by --export=NONE + module
```

Thirteen 3-atom forward passes; CPU is ample. Add `--gres=gpu:1` and
`--device cuda` only if the checkpoint refuses to load on CPU.

**2. CP2K inputs — login node.**

```bash
python make_cp2k_inputs.py --template ~/scw/mlwf_test/mlwf.inp
```

The template is read and edited in exactly one place, the cell, with the same
substitution `make_dipole_labels_v3.py` uses (`re.sub(r"ABC\s+[\d. ]+", ...)`)
set to `ABC 20.0000 20.0000 20.0000`. Nothing else is parsed or assumed — the DFT
here has to be the same calculation that produced the training labels, or a
model/DFT disagreement could just be a settings difference.

**3. Wannier runs — sbatch array of 13.**

```bash
cat > monomer_dft.sh <<'EOF'
#!/bin/bash
#SBATCH --job-name=monomer-mlwf
#SBATCH --partition=short
#SBATCH --time=00:20:00
#SBATCH --ntasks=8
#SBATCH --array=0-12
#SBATCH --export=NONE
#SBATCH --output=%x-%A_%a.out
#SBATCH --error=%x-%A_%a.err
module purge
module load CP2K/2023.1-foss-2023a
cd "$SLURM_SUBMIT_DIR/$(printf 'g%02d' "$SLURM_ARRAY_TASK_ID")"
srun cp2k.psmp -i dip.inp -o dip.out
EOF
sbatch monomer_dft.sh
```

> This wrapper follows the resources the brief specifies (CP2K/2023.1-foss-2023a,
> 8 tasks, 20 min). It is **not** copied from `htc_dipole_label3.sh`, which is on
> ARC and was not readable from here — check the partition, account and `srun`
> invocation against it before submitting, and prefer its versions of those.

**4. Harvest and compare — login node, numpy only.**

```bash
python harvest_dft.py --builder ~/scw/build_dipole_xyz_v3.py
python compare.py
```

`harvest_dft.py` takes `dip()` from the label builder by reading the file and
executing it **only up to the line `out=[]`** — never importing it. Everything
after that marker is the harvest loop that rewrites the training set; truncating
there takes the definitions and leaves the side effects behind.

## Reading the result

`compare.py` prints μ at equilibrium (model / DFT / experiment), dμ/dQ by central
differences for all three modes at both references with |model|/|DFT| per mode,
and the squared-magnitude stretch/bend ratio per source — the intensity proxy the
dilute-limit result is about, with the bend derivative converted from D/deg to
D/Å using the arc length r·dθ so the ratio is dimensionless.

- **stretch-derivative ratio near 1** → the dipole model is *not* the cause of the
  stretch excess; look at the dynamics or the line list.
- **1.5–2 or more** → it is: a liquid-like dμ/dr in the gas, squared into
  intensity.

## What has been checked from here

- All 13 geometries build, and `config_from_atoms` → `AtomicData.from_config` →
  `DataLoader` → `batch.to_dict()` runs for every one of them
  (`python eval_model.py --dry-run`, which needs no checkpoint). Only the forward
  pass and the `"dipole"` output key are untested locally, because
  `water_dipole_v3.model` is on ARC.
- `compare.py`'s arithmetic was verified against synthetic inputs carrying a
  known 1.8× stretch bias: it recovered 1.800 on every stretch derivative, 1.000
  on bend, and an intensity-proxy ratio of exactly 3.240 = 1.8².
