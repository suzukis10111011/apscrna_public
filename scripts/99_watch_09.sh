#!/usr/bin/env bash
# Purpose: Monitor completion of stage 08 (kb requantification for 5 discovery samples) every 10 minutes; when
#   all verify.json files exist, auto-run stage 09 and dependent pipelines
# Input: logs/08.lock, logs/08_requant.log (read-only), data/counts_discovery/Sample*/verify.json
# Output: logs/08.lock (no writes), downstream outputs from stages 03/04/09/09b
# 99_watch_09: Monitor every 10 minutes for the completion of 08 (re-quantification of kb for the 5 discovery samples), and once aligned, automatically execute step 09.
# Never write to 08's lock (logs/08.lock) or log (logs/08_requant.log) (read only).
#
# Monitoring target is data/counts_discovery/Sample*/verify.json (self-judgment 16: the instructed data/counts/ is
# for the verification cohort CAP1-3; 08's output destination is counts_discovery).
#
# Steps after completion (each stage skips if output already exists, so restarting does not cause duplicate execution):
#   1) APSC_SOURCE=kb    03 -> 04   (kb all barcodes + same QC)
#   2) APSC_SOURCE=kbgeo 03 -> 04   (restricted to GEO barcodes)
#   3) 09  comparison table
#   4) 09b REPORT_day3.md / docs/dev_log_ja.md / 16_methods_numbers.csv
set -uo pipefail
W="${APSC_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")/.." && pwd)}"
PY=${APSC_PYTHON:-python}
CNT=$W/data/counts_discovery
L=$W/logs/99_watch_09.log
LOCK=$W/logs/09.lock
MAXWAIT=$((72*3600))
cd "$W"
export NUMBA_CACHE_DIR=$W/.numba_cache
export OMP_NUM_THREADS=5 OPENBLAS_NUM_THREADS=5 MKL_NUM_THREADS=5 NUMBA_NUM_THREADS=5
say(){ echo "[$(date -u +%FT%TZ)] $*" >> $L; }
beat(){ date -u +%FT%TZ > $W/logs/watch09_heartbeat.txt; }

mkdir "$LOCK" 2>/dev/null || { say "別の watcher が稼働中（$LOCK）— 起動しない"; exit 0; }
trap 'rmdir "$LOCK" 2>/dev/null' EXIT

ready(){ ls $CNT/Sample*/verify.json 2>/dev/null | wc -l | tr -d ' '; }
say "=== watcher 起動（10分周期、最大 $((MAXWAIT/3600)) 時間）現在 $(ready)/5 ==="
beat

waited=0
while :; do
  n=$(ready)
  if [ "$n" -ge 5 ]; then say "5検体そろった（$n/5）— 09 工程に進む"; break; fi
  if grep -q "=== 08 DONE" $W/logs/08_requant.log 2>/dev/null; then
    say "08 が完了行を出したが verify.json は $n/5 — kb 版は不完全なので 09 は実行しない"
    say "（03 は5検体を前提にしているため。不足分を 08 で埋めてから watcher を再起動すること）"
    beat; exit 1
  fi
  beat
  sleep 600; waited=$((waited+600))
  if [ $waited -ge $MAXWAIT ]; then say "$((MAXWAIT/3600)) 時間待っても揃わず — 打ち切り（$n/5）"; beat; exit 1; fi
  if [ $((waited % 3600)) -eq 0 ]; then say "待機中 $((waited/3600))h — $n/5、空き $(df -k $W | tail -1 | awk '{print int($4/1048576)}')G"; fi
done

run(){   # run <label> <marker file to skip if already done> <log> <command...>
  local lab=$1 mark=$2 log=$3; shift 3
  if [ -s "$mark" ]; then say "$lab: 既に出力あり（$(basename $mark)）— 飛ばす"; return 0; fi
  say "$lab: 開始"
  if "$@" > "$log" 2>&1; then say "$lab: 完了"; beat; return 0
  else say "$lab: 失敗（$(basename $log) の末尾: $(tail -2 "$log" | tr '\n' ' ' | cut -c1-200)）"; beat; return 1; fi
}

TAB=$W/results/tables; OUT=$W/data/processed
APSC_SOURCE=kb    run "03 kb"    "$OUT/all_cells_discovery_kb.h5ad"    "$W/logs/03_kb.log"    env APSC_SOURCE=kb    $PY -u scripts/03_preprocess_discovery.py
APSC_SOURCE=kb    run "04 kb"    "$OUT/stromal_discovery_kb.h5ad"      "$W/logs/04_kb.log"    env APSC_SOURCE=kb    $PY -u scripts/04_discovery_scores_and_screen.py
APSC_SOURCE=kbgeo run "03 kbgeo" "$OUT/all_cells_discovery_kbgeo.h5ad" "$W/logs/03_kbgeo.log" env APSC_SOURCE=kbgeo $PY -u scripts/03_preprocess_discovery.py
APSC_SOURCE=kbgeo run "04 kbgeo" "$OUT/stromal_discovery_kbgeo.h5ad"   "$W/logs/04_kbgeo.log" env APSC_SOURCE=kbgeo $PY -u scripts/04_discovery_scores_and_screen.py

# 09 and 09b are run every time (to incorporate the increased input)
say "09: 開始"; $PY -u scripts/09_geo_vs_kb.py    > $W/logs/09.log  2>&1 && say "09: 完了"  || say "09: 失敗"
beat
say "09b: 開始"; $PY -u scripts/09b_report_day3.py > $W/logs/09b.log 2>&1 && say "09b: 完了" || say "09b: 失敗"
beat

git -C $W add -A >/dev/null 2>&1
git -C $W commit -qm "09: kb requant vs GEO matrix comparison (auto, 99_watch_09)" >/dev/null 2>&1
say "commit: $(git -C $W log --oneline | head -1)"
WS="${APSC_WORKSPACE:-}"        # optional: where to copy deliverables; skipped when unset
if [ -d "$WS" ]; then
  [ -n "$WS" ] && cp $W/REPORT_day3.md $W/docs/dev_log_ja.md "$WS/" 2>/dev/null
  for f in 09_geo_vs_kb 09_tf_agreement_geo_vs_kb 09_stromal_jaccard 16_methods_numbers; do
    [ -n "$WS" ] && cp $TAB/$f.csv "$WS/" 2>/dev/null
  done
  say "workspace に成果物を配置"
fi
say "=== watcher DONE（kb $(ls $TAB/04_*_kb.csv 2>/dev/null | wc -l | tr -d ' ') 表 / kbgeo $(ls $TAB/04_*_kbgeo.csv 2>/dev/null | wc -l | tr -d ' ') 表）==="
beat
