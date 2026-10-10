#!/usr/bin/env bash
# Purpose: Re-run every analysis stage that consumes the processed h5ad objects (quantification and
#   preprocessing excluded), writing tables to results/tables_rerun and objects to data/processed_rerun.
# Input: data/processed_rerun/*.h5ad (a copy of data/processed), data/counts, data/gse164241, frozen/, ref/.
# Output: results/tables_rerun/*.csv, logs/rerun_<stage>.log, logs/rerun_status.tsv.
set -uo pipefail
W="${APSC_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")/.." && pwd)}"
cd "$W"
export APSC_TABDIR="$W/results/tables_rerun"
export APSC_OUTDIR="$W/data/processed_rerun"
export APSC_PRIVATE="$W/docs/private/tables_rerun"
export NUMBA_CACHE_DIR="$W/.numba_cache" OMP_NUM_THREADS=5 NUMBA_NUM_THREADS=5 APSC_REQUANT_DEAD=1
mkdir -p "$APSC_TABDIR" "$APSC_OUTDIR" "$APSC_PRIVATE" logs
STATUS=logs/rerun_status.tsv
: > $STATUS
STAGES="03b_fix_compartment 03c_recheck_discovery 04_discovery_scores_and_screen 04b_state_labels_discovery 05_validation_and_projection 05e_rho_prev_cellset 09_geo_vs_kb 10_receptors_by_state 11_coactivity 12_coupling_by_compartment 13_gse164241 14_robustness 14c_robustness_integration 17_receptor_power 18_stromal_markerdefined 20b_coverage_after_fix 25_validation_states_within 28_manuscript_numbers 32_coactivity_ruleset 33_ms_figure_inputs 33_extract_figure_sources 34_lolo_rulebased 37_diag_fig2c"
for s in $STAGES; do
  t0=$(date +%s)
  printf "[%s] %s start\n" "$(date -u +%FT%TZ)" "$s"
  if python -u "scripts/$s.py" > "logs/rerun_$s.log" 2>&1; then
    printf "%s\tOK\t%s\n" "$s" "$(( $(date +%s) - t0 ))" >> $STATUS
    printf "[%s] %s OK (%ss)\n" "$(date -u +%FT%TZ)" "$s" "$(( $(date +%s) - t0 ))"
  else
    printf "%s\tFAILED\t%s\n" "$s" "$(( $(date +%s) - t0 ))" >> $STATUS
    printf "[%s] %s FAILED (%ss) — continuing\n" "$(date -u +%FT%TZ)" "$s" "$(( $(date +%s) - t0 ))"
  fi
done
printf "[%s] driver done: %s OK / %s FAILED\n" "$(date -u +%FT%TZ)" \
  "$(grep -c OK $STATUS || true)" "$(grep -c FAILED $STATUS || true)"
