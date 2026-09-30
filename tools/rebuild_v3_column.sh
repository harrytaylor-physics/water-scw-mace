#!/bin/bash
# A85: rebuild the SAME producers against the v3 spectra, into their own output
# directory, so the "v3" column of any v3->v5 comparison is a genuine v3 value
# rather than whatever file happened to be on disk. A83's v3 column was not:
# outputs/S6_candidate_evolution.csv had not been rebuilt since commit 4834c22.
#
# Option A: flip FP_OPACITY_SUBDIR and point output_dir at outputs_v3/, run, flip
# back. Every producer honours cfg.output_path, and the only hard-coded
# "outputs/" strings in the tree are provenance labels, not write paths.
set -e
cd "$(dirname "$0")/.."
SCHEMA=config/schema.py
CONF=/tmp/scratch-501/fiducial_v3.yaml
CAND=/tmp/scratch-501/candidate_v3.yaml
mkdir -p outputs_v3
python - <<'PY'
import pathlib, re
for src,dst in (("config/fiducial.yaml","/tmp/scratch-501/fiducial_v3.yaml"),
                ("config/candidate.yaml","/tmp/scratch-501/candidate_v3.yaml")):
    s=pathlib.Path(src).read_text()
    if re.search(r'^\s*output_dir:', s, re.M):
        s=re.sub(r'^(\s*)output_dir:.*$', r'\1output_dir: outputs_v3', s, flags=re.M)
    else:
        s=re.sub(r'^(paths:\s*)$', r'\1\n  output_dir: outputs_v3', s, flags=re.M)
    pathlib.Path(dst).write_text(s)
print("v3 configs written")
PY
cp $SCHEMA /tmp/scratch-501/schema.bak
trap 'cp /tmp/scratch-501/schema.bak '"$SCHEMA"'; echo "schema restored"' EXIT
sed -i '' 's/FP_OPACITY_SUBDIR = "fp_opacity_v5"/FP_OPACITY_SUBDIR = "fp_opacity"/' $SCHEMA
grep -q 'FP_OPACITY_SUBDIR = "fp_opacity"' $SCHEMA || { echo "flip failed"; exit 1; }
echo "=== rebuilding producers against v3 spectra into outputs_v3/ ==="
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
for f in outputs/broadening_sensitivity.json outputs/aqua_entropy_check.json \
         outputs/aguichine2025_compare.json outputs/crossover_v3_before_a60.json \
         outputs/crossover_v3_after_a60.json; do
  [ -f "$f" ] && cp "$f" outputs_v3/   # FP-independent; copied so the ledger runs
done
python tools/manuscript_numbers.py $CONF >/dev/null
cp outputs_v3/manuscript_numbers.json outputs/manuscript_numbers_v3_rebuilt.json
echo "=== v3 column written to outputs/manuscript_numbers_v3_rebuilt.json ==="
