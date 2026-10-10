#!/usr/bin/env bash
# Purpose: Orchestrate quantization and verification pipeline with progress monitoring and idempotent restarts
# Input: none (reads logs/quant_progress.txt for state; executes 02_quant.sh → 05_* → 06_*)
# Output: verify.json, quantization logs in logs/, stage-specific artifacts from 02/05/06
# 99 v2: Restart-version driver. Auto-runs 02_quant.sh -> 05 -> 06.
# Changes from v1:
#   - python specified with full path (in v1, `python` resolved to ./.venv/apsc/bin/python and failed)
#   - Waiting is not a "fixed 20-hour timeout" but "abort if no progress for 3 hours"
#     (progress = update of logs/quant_progress.txt. 02_quant.sh touches it at each stage)
#   - Idempotent: stages with verify.json / output already present are skipped
set -uo pipefail
W="${APSC_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")/.." && pwd)}"
cd $W
export NUMBA_CACHE_DIR=$W/.numba_cache; mkdir -p $NUMBA_CACHE_DIR
PY=${APSC_PYTHON:-python}
L=$W/logs/99_driver.log
STALL=10800                                   # 3 hours
say(){ echo "[$(date -u +%FT%TZ)] $*" >> $L; echo "[$(date -u +%FT%TZ)] $*" >> $W/STATUS_APPEND.log; }
beat(){ date -u +%FT%TZ > $W/logs/driver_heartbeat.txt; }
run(){ local s=$1 tag=$2
  say "$tag: start"
  if $PY -u "$s" > "$W/logs/${tag}.log" 2>&1; then
    say "$tag: DONE"; git add -A >/dev/null 2>&1; git commit -qm "$tag done" >/dev/null 2>&1; return 0
  fi
  say "$tag: FAILED ($(tail -3 "$W/logs/${tag}.log" | tr '\n' ' ' | cut -c1-300))"; return 1
}
beat; say "=== driver v2 起動 (pid $$) ==="

# ---- Quantification (run 02_quant.sh in the background and monitor progress) ----
ncnt(){ ls $W/data/counts/*/verify.json 2>/dev/null | wc -l | tr -d ' '; }
if [ "$(ncnt)" -lt 3 ]; then
  if [ ! -f $W/logs/quant_running.pid ] || ! kill -0 "$(cat $W/logs/quant_running.pid 2>/dev/null)" 2>/dev/null; then
    say "02_quant: 起動"
    nohup bash $W/scripts/02b_quant.sh > $W/logs/02_nohup.log 2>&1 &
    echo $! > $W/logs/quant_running.pid
  else
    say "02_quant: 既に稼働中 (pid $(cat $W/logs/quant_running.pid))"
  fi
  P=$W/logs/quant_progress.txt
  while [ "$(ncnt)" -lt 3 ]; do
    beat; sleep 300
    if [ -f "$P" ]; then
      age=$(( $(date +%s) - $(stat -f%m "$P") ))
      if [ $age -ge $STALL ]; then say "02_quant: 進捗が $((age/3600)) 時間無い — 待機を打ち切る（完了 $(ncnt)/3）"; break; fi
    fi
    if ! kill -0 "$(cat $W/logs/quant_running.pid 2>/dev/null)" 2>/dev/null; then
      say "02_quant: プロセス終了を検出（完了 $(ncnt)/3）"; break
    fi
  done
fi
say "定量完了数: $(ncnt)/3"

# ---- Verification side ----
if [ "$(ncnt)" -ge 1 ]; then
  if [ -s $W/data/processed/stromal_validation.h5ad ]; then say "05: 既存出力あり — skip"
  else run scripts/05_validation_and_projection.py 05; fi
else
  say "05: スキップ（定量0件）"
fi

# ---- Comparison with previous CAP1 ----
if [ -s $W/data/counts/CAP1/verify.json ]; then
  $PY - <<'PY' >> $L 2>&1
import json
v = json.load(open("$W/data/counts/CAP1/verify.json"))
prev = dict(barcodes=1953956, genes=35689, nonzero=18856281, n_processed=435676976, p_pseudoaligned=69.0)
out = {}
for k, p in prev.items():
    c = v.get(k)
    out[k] = dict(previous=p, current=c, same=(c == p) if isinstance(p, int) else (abs((c or 0) - p) < 0.5))
print("CAP1 前回 vs 今回:", json.dumps(out, ensure_ascii=False))
json.dump(dict(previous=prev, current=v, comparison=out),
          open("$W/results/tables/02_CAP1_prev_vs_new.json", "w"), indent=1)
PY
fi

# 06_report.py handles figures for the application document, so it is not included in the public repository (docs/private/scripts/).
  if [ -f scripts/06_report.py ]; then run scripts/06_report.py 06 || say "06: レポート生成失敗"; \
  else say "06: 非公開のため省略（docs/private/scripts/06_report.py）"; fi
beat; say "=== 99 driver v2 終了（定量 $(ncnt)/3）==="
