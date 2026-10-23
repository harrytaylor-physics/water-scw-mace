#!/bin/bash
#SBATCH --job-name=subnep-broadening
#SBATCH --partition=short
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#
# ARC hand-off for Task 21 (A59): HITEMP broadening-convention sensitivity.
#
#   sbatch tools/broadening_sensitivity.sbatch
#
# Only needed if the Mac cannot run it — RADIS 0.16.3 and the HITEMP cache are
# both present locally, so the default route is
#
#   python tools/broadening_sensitivity.py config/fiducial.yaml
#
# Two things this wrapper gets that the Mac does not:
#
#   * the 0-50 cm^-1 HITEMP file, which the local ~/.radisdb cache is missing.
#     With it, pass --wavenum-min 10 to reproduce the delivered spectra's line
#     loading range exactly. Without it the science band [100, 4000] is still
#     unaffected: with truncation = 50 cm^-1 a line below 50 cannot reach past
#     100. The local run therefore uses 50 and says so.
#   * memory headroom. Six points x three conventions at wstep = 0.01 cm^-1 over
#     [10, 11000] is ~1.1e6 points per spectrum plus the full H2O line set.
#
# NOTE on the RADIS API: ~/scw/hitemp_opacity.py and generate_linelist_tables.py
# set `sf.pressure` between calls. RADIS >= 0.16 made SpectrumFactory slotted and
# that raises AttributeError; the pressure is a per-call argument to eq_spectrum
# now. tools/broadening_sensitivity.py uses the new form. If ARC's RADIS is older
# than 0.16 the new form still works — `pressure=` has been an eq_spectrum kwarg
# throughout — so no branch is needed.

set -euo pipefail

module purge
module load Anaconda3/2025.06-1

REPO="${REPO:-$HOME/subneptune}"
CONFIG="${CONFIG:-config/fiducial.yaml}"

# RADIS reads its line cache from $DATA/radisdb (the convention hitemp_opacity.py
# and generate_linelist_tables.py already use on this cluster).
export DATA="${DATA:-$HOME/scw}"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"

cd "$REPO"

echo "host=$(hostname)  date=$(date -Is)"
echo "python=$(which python)  DATA=$DATA"
python -c "import radis; print('radis', radis.__version__)"

# Broadening pressures first: cheap, and it fails fast if AQUA is not readable.
python tools/broadening_sensitivity.py "$CONFIG" --dry-run

python tools/broadening_sensitivity.py "$CONFIG" --wavenum-min 10

echo "done $(date -Is)"
echo "copy back with:"
echo "  scp -r <user>@htc-login.arc.ox.ac.uk:$REPO/outputs/broadening_sensitivity.json \\"
echo "         <user>@htc-login.arc.ox.ac.uk:$REPO/outputs/broadening \\"
echo "         ~/subneptune/outputs/"
echo "then: python figures/plot_ed1_broadening.py $CONFIG"
