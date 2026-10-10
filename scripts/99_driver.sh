#!/usr/bin/env bash
# Purpose: Entry point for the unattended driver; execs the current version (99_driver_v2.sh).
# Input: none (the superseded v1 body is kept below as a comment, unexecuted).
# Output: none of its own; see logs/99_driver.log written by the version it execs.
# (The old version is saved as 99_driver_v1.sh.bak. The following is the restarted version)
W="${APSC_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")/.." && pwd)}"
exec bash "$W/scripts/99_driver_v2.sh"
# ---- Old version starts here (not executed) ----
# 99: Self-running driver. Executes the remainder in sequence even after the agent's turn has ended.
#  Dependencies: 03 -> 04 -> 04b, 01c -> 05. If 05 has the outputs of 04/04b, it also performs (3)② and state projection.
#  Log failed steps and proceed to the next independent step (principle per user instructions).
set -uo pipefail
W="${APSC_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")/.." && pwd)}"
cd $W
export NUMBA_CACHE_DIR=$W/.numba_cache; mkdir -p $NUMBA_CACHE_DIR
L=$W/logs/99_driver.log
say(){ echo "[$(date -u +%FT%TZ)] $*" >> $L; echo "[$(date -u +%FT%TZ)] $*" >> $W/STATUS_APPEND.log; }
beat(){ date -u +%FT%TZ > $W/logs/driver_heartbeat.txt; }   # For liveness check (updated every 60 seconds)
beat; say "=== driver 起動 (pid $$) ==="
run(){ # $1=script $2=tag
  local s=$1 tag=$2
  say "$tag: start"
  if python -u "$s" > "$W/logs/${tag}.log" 2>&1; then
    say "$tag: DONE"; git add -A >/dev/null 2>&1; git commit -qm "$tag done" >/dev/null 2>&1; return 0
  else
    say "$tag: FAILED (末尾: $(tail -3 "$W/logs/${tag}.log" | tr '\n' ' ' | cut -c1-300))"; return 1
  fi
}

# --- Discovery series: wait for completion of 03 (up to 3 hours) ---
w=0
while [ ! -s data/processed/all_cells_discovery.h5ad ]; do
  beat; sleep 60; w=$((w+60))
  [ $w -ge 10800 ] && { say "03 待機タイムアウト — 探索系を中止"; break; }
done
if [ -s data/processed/all_cells_discovery.h5ad ]; then
  say "03: 出力確認"
  if [ -s results/tables/04_tf_screen_discovery_candidates.csv ]; then say "04: 既存出力あり — skip"; else
    run scripts/04_discovery_scores_and_screen.py 04 || say "04b: スキップ（依存先 04 失敗）"; fi
  if [ -s data/processed/stromal_discovery_labelled.h5ad ]; then say "04b: 既存出力あり — skip"; else
    [ -s results/tables/04_tf_screen_discovery_candidates.csv ] && run scripts/04b_state_labels_discovery.py 04b; fi
else
  say "04/04b: スキップ（依存先 03 失敗）"
fi

# --- Validation series: wait for completion of 3 runs of 01c (up to 20 hours) ---
w=0
while [ "$(ls data/counts/*/verify.json 2>/dev/null | wc -l | tr -d ' ')" -lt 3 ]; do
  beat; sleep 300; w=$((w+300))
  if [ $w -ge 72000 ]; then
    say "01c 待機タイムアウト（$(ls data/counts/*/verify.json 2>/dev/null | wc -l | tr -d ' ')/3 完了）— 揃った分で 05 を実行"
    break
  fi
done
n=$(ls data/counts/*/verify.json 2>/dev/null | wc -l | tr -d ' ')
say "検証カウント完了数: $n/3"
if [ "$n" -ge 1 ]; then
  run scripts/05_validation_and_projection.py 05
else
  say "05: スキップ（依存先 01c 失敗、カウント0件）"
fi

# --- Comparison with previous CAP1 ---
if [ -s data/counts/CAP1/verify.json ]; then
  python - <<'PY' 2>&1 | tee -a $L
import json
v = json.load(open("$W/data/counts/CAP1/verify.json"))
prev = dict(barcodes=1953956, genes=35689, nonzero=18856281, n_processed=435676976, p_pseudoaligned=69.0)
print("=== CAP1 前回 vs 今回 ===")
for k, p in prev.items():
    c = v.get(k)
    same = (c == p) if isinstance(p, int) else (abs((c or 0) - p) < 0.5)
    print(f"  {k:16s} 前回 {p} / 今回 {c} -> {'一致' if same else '不一致'}")
print("今回 total_umi:", v.get("total_umi"), "(前回値は記録なし)")
json.dump(dict(previous=prev, current=v), open("$W/results/tables/01_CAP1_prev_vs_new.json","w"), indent=1)
PY
  git add -A >/dev/null 2>&1; git commit -qm "CAP1 comparison" >/dev/null 2>&1
fi

# 06_report.py handles figures for the application document, so it is not included in the public repository (docs/private/scripts/).
  if [ -f scripts/06_report.py ]; then run scripts/06_report.py 06 || say "06: レポート生成失敗"; \
  else say "06: 非公開のため省略（docs/private/scripts/06_report.py）"; fi
say "=== 99 driver 終了 ==="
