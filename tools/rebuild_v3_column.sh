#!/bin/bash
# A85/A87: rebuild the SAME producers against the v3 spectra, fully isolated, so
# the "v3" column of a v3->v5 comparison is a genuine v3 value.
#
# Three things have to be redirected, not one. cfg.output_path covers the scan
# and evolution CSVs and the ledger; but the figure and Extended Data ledgers go
# to figs/ and the release products to data/release/, both hard-coded repo paths
# until A87 made them redirectable. A86's v3 column was wrong precisely because
# figs/ was not isolated: it carried v5 figure numbers, so those groups compared
# v5 against v5 and reported "unchanged".
set -e
cd "$(dirname "$0")/.."
CONF=/tmp/scratch-501/fiducial_v3.yaml
CAND=/tmp/scratch-501/candidate_v3.yaml
export SUBNEPTUNE_FIGS_DIR="$PWD/figs_v3"
export SUBNEPTUNE_RELEASE_DIR="$PWD/data_release_v3"
rm -rf outputs_v3 figs_v3 data_release_v3
cp -R outputs outputs_v3          # seed the auxiliaries producers read but never write
mkdir -p figs_v3 data_release_v3
python - <<'PY'
import pathlib, re
for src,dst in (("config/fiducial.yaml","/tmp/scratch-501/fiducial_v3.yaml"),
                ("config/candidate.yaml","/tmp/scratch-501/candidate_v3.yaml")):
    s=pathlib.Path(src).read_text()
    s=re.sub(r'^(\s*)output_dir:.*$', r'\1output_dir: outputs_v3', s, flags=re.M)
    pathlib.Path(dst).write_text(s)
PY
cp config/schema.py /tmp/scratch-501/schema.bak
cp opacity/registry.py /tmp/scratch-501/registry.bak
trap 'cp /tmp/scratch-501/schema.bak config/schema.py; cp /tmp/scratch-501/registry.bak opacity/registry.py; echo "sources restored"' EXIT
sed -i '' 's/FP_OPACITY_SUBDIR = "fp_opacity_v5"/FP_OPACITY_SUBDIR = "fp_opacity"/' config/schema.py
grep -q 'FP_OPACITY_SUBDIR = "fp_opacity"' config/schema.py || { echo "flip failed"; exit 1; }
echo "=== v3 producers -> outputs_v3/, figs_v3/, data_release_v3/ ==="
python figures/plot_s4_domain_scan.py $CONF --extended >/dev/null
python tools/corner_boundaries.py $CONF >/dev/null
python figures/plot_s8_rcb_flip.py $CONF >/dev/null
python figures/plot_s9_case_evolution.py $CONF >/dev/null
python figures/plot_s6_candidate_evolution.py $CAND >/dev/null   # candidate config
python tools/corner_cells.py $CONF --crossover >/dev/null
python tools/band_diagnostic.py $CONF >/dev/null
python tools/text_checks.py $CONF >/dev/null
python tools/rho_scaling.py $CONF >/dev/null
python tools/wing_cutoff.py $CONF --summarise-only >/dev/null
python tools/release_tables.py $CONF >/dev/null
for f in figures/plot_*.py figures/fig1.py figures/fig2.py figures/fig3.py figures/fig4.py; do
  case "$(basename $f)" in
    plot_s6_candidate_evolution.py|plot_s4_domain_scan.py) continue ;;
  esac
  python $f $CONF >/dev/null 2>&1 || echo "    FAIL $(basename $f)"
done
# ledger BEFORE ed_index, which merges ed.* into it
# A87: the ED *figures* populate ed_numbers.json and must be built here. They
# were not, in either A86 pass: v5 inherited 341 ed.* entries from an earlier run
# and v3 started empty, so the Extended Data group was never actually compared.
python figures/ed_items.py $CONF >/dev/null 2>&1 || echo "    FAIL ed_items"
python figures/ed_tables.py $CONF >/dev/null 2>&1 || echo "    FAIL ed_tables"
python tools/manuscript_numbers.py $CONF >/dev/null
python figures/ed_index.py $CONF >/dev/null 2>&1 || echo "    FAIL ed_index"
cp outputs_v3/manuscript_numbers.json outputs/manuscript_numbers_v3_rebuilt.json
echo "=== v3 column written to outputs/manuscript_numbers_v3_rebuilt.json ==="
