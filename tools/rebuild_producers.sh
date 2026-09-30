#!/bin/bash
# A85: rebuild every producer the manuscript ledger reads, in dependency order.
# The ledger READS stored files; it does not rebuild them, which is how a stale
# CSV survived the v5 switch. Run this before the ledger, always.
#   bash tools/rebuild_producers.sh [outputs_dir_override]
set -e
cd "$(dirname "$0")/.."
C=config/fiducial.yaml
P=python
echo "--- 1 extended domain scan (108 boundaries; needs --extended, ~5 min)"
$P figures/plot_s4_domain_scan.py $C --extended >/dev/null
echo "--- 2 corner boundaries"; $P tools/corner_boundaries.py $C >/dev/null
echo "--- 3 S8 RCB flip";      $P figures/plot_s8_rcb_flip.py $C >/dev/null
echo "--- 4 S9 case evolution"; $P figures/plot_s9_case_evolution.py $C >/dev/null
# NOTE: S6 takes config/candidate.yaml, not fiducial. Its first irradiation point
# is T_eq = 700 K; running it on fiducial (500 K) puts the entropy-matched solve
# outside the default bracket and was the source of A83's spurious 6.07 (A84).
echo "--- 5 S6 candidate evolution (candidate.yaml, NOT fiducial)"
$P figures/plot_s6_candidate_evolution.py config/candidate.yaml >/dev/null
echo "--- 6 crossover";        $P tools/corner_cells.py $C --crossover >/dev/null
echo "--- 7 band diagnostic";  $P tools/band_diagnostic.py $C >/dev/null
echo "--- 8 text checks";      $P tools/text_checks.py $C >/dev/null
echo "--- 9 rho^2 scaling";    $P tools/rho_scaling.py $C >/dev/null
echo "--- 10 wing cut-off";    $P tools/wing_cutoff.py $C --summarise-only >/dev/null
echo "--- 11 release tables";  $P tools/release_tables.py $C >/dev/null
echo "--- 12 figures"
for f in figures/plot_*.py figures/fig1.py figures/fig2.py figures/fig3.py figures/fig4.py; do
  case "$(basename $f)" in
    plot_s6_candidate_evolution.py) continue ;;          # done above, other config
    plot_s4_domain_scan.py) continue ;;                  # done above, --extended
  esac
  $P $f $C >/dev/null 2>&1 || echo "    FAIL $(basename $f)"
done
# ORDER MATTERS HERE. ed_index.py MERGES its numbers into
# outputs/manuscript_numbers.json under ed.* keys, so the ledger has to be built
# BEFORE it, not after -- running the ledger last overwrites the merge and
# silently drops 341 ed.* entries. ed_tables feeds ed_index, so it comes first.
echo "--- 13 ED tables";  $P figures/ed_tables.py $C >/dev/null 2>&1 || echo "    FAIL ed_tables"
echo "--- 14 ledger";     $P tools/manuscript_numbers.py $C >/dev/null
echo "--- 15 ED merge";   $P figures/ed_index.py $C >/dev/null 2>&1 || echo "    FAIL ed_index"
echo "--- rebuild complete"
