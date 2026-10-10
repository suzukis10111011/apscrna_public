#!/usr/bin/env bash
# Purpose: Wait for 02b_quant.sh to finish, then run steps 05 through 06 unattended.
# Input: none detected
# Output: none detected
# 99 v3: Wait for 02b_quant.sh to complete, then run 05 -> 06.
# v2 bug: `kill -0 $(cat quant_running.pid)` returned "absent" even though 02b was alive,
#   immediately exiting the wait loop, skipping 05 with quantification at 0/3, and terminating.
#   Since process queries are unstable in this sandbox, **do not use pid at all**.
# Completion is judged by any of: (a) all three verify.json files present, (b) "02b DONE" in the 02b log, (c) progress stalls for 3 hours.
set -uo pipefail
W="${APSC_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")/.." && pwd)}"
cd $W
export NUMBA_CACHE_DIR=$W/.numba_cache; mkdir -p $NUMBA_CACHE_DIR
PY=${APSC_PYTHON:-python}
L=$W/logs/99_driver.log; QL=$W/logs/02b_quant.log; P=$W/logs/quant_progress.txt
STALL=10800
say(){ echo "[$(date -u +%FT%TZ)] [v3] $*" >> $L; echo "[$(date -u +%FT%TZ)] $*" >> $W/STATUS_APPEND.log; }
beat(){ date -u +%FT%TZ > $W/logs/driver_heartbeat.txt; }
ncnt(){ ls $W/data/counts/*/verify.json 2>/dev/null | wc -l | tr -d ' '; }
run(){ local s=$1 tag=$2
  say "$tag: start"
  if $PY -u "$s" > "$W/logs/${tag}.log" 2>&1; then
    say "$tag: DONE"; git add -A >/dev/null 2>&1; git commit -qm "$tag done" >/dev/null 2>&1; return 0
  fi
  say "$tag: FAILED ($(tail -3 "$W/logs/${tag}.log" | tr '\n' ' ' | cut -c1-300))"; return 1
}
beat; say "=== driver v3 起動 ==="

while :; do
  beat
  [ "$(ncnt)" -ge 3 ] && { say "定量 3/3 完了"; break; }
  if grep -q "02b DONE" "$QL" 2>/dev/null; then say "02b の終了を検出（完了 $(ncnt)/3）"; break; fi
  if [ -f "$P" ]; then
    age=$(( $(date +%s) - $(stat -f%m "$P") ))
    [ $age -ge $STALL ] && { say "進捗が $((age/3600)) 時間停止 — 待機打ち切り（完了 $(ncnt)/3）"; break; }
  fi
  sleep 300
done
say "定量完了数: $(ncnt)/3"

if [ "$(ncnt)" -ge 1 ]; then
  if [ -s $W/data/processed/stromal_validation.h5ad ]; then say "05: 既存出力あり — skip"
  else run scripts/05_validation_and_projection.py 05; fi
else
  say "05: スキップ（定量0件）"
fi

if [ -s $W/data/counts/CAP3/verify.json ]; then
  # The previous kb record uses the internal label "CAP1", but its actual entity is run SRR32119705.
  # Per ref/run_donor_map_v2.csv, SRR32119705 = CAP3, so cross-check against CAP3.
  $PY - <<'PY' >> $L 2>&1
import json, sys
sys.path.insert(0, "$W/scripts")
from _utils import load_run_donor_map
_m, r2d, _d2r = load_run_donor_map()
donor = r2d["SRR32119705"]
v = json.load(open(f"$W/data/counts/{donor}/verify.json"))
prev = dict(barcodes=1953956, genes=35689, nonzero=18856281, n_processed=435676976, p_pseudoaligned=69.0)
out = {k: dict(previous=p, current=v.get(k),
               same=(v.get(k) == p) if isinstance(p, int) else (abs((v.get(k) or 0) - p) < 0.5))
       for k, p in prev.items()}
rec = dict(run="SRR32119705", donor_now=donor, previous_internal_label="CAP1",
           note="前回パイプラインの内部ラベルは誤った対応表由来。突き合わせはラン番号で行う。",
           previous=prev, current=v, comparison=out)
print("SRR32119705 前回 vs 今回:", json.dumps(out, ensure_ascii=False))
json.dump(rec, open("$W/results/tables/02_SRR32119705_prev_vs_new.json","w"), indent=1)
PY
fi

# 06_report.py handles figures for the application document, so it is not included in the public repository (docs/private/scripts/).
  if [ -f scripts/06_report.py ]; then run scripts/06_report.py 06 || say "06: レポート生成失敗"; \
  else say "06: 非公開のため省略（docs/private/scripts/06_report.py）"; fi
beat; say "=== driver v3 終了（定量 $(ncnt)/3）==="
