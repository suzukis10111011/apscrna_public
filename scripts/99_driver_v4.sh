#!/usr/bin/env bash
# Purpose: Unattended driver (v4): waits for the discovery requantification to finish, then runs the
#   remaining analysis steps, judging liveness by output growth rather than by process lookup.
# Input: data/counts_discovery/Sample*/verify.json; logs/requant_heartbeat.txt.
# Output: logs/99_driver.log, logs/driver_heartbeat.txt, and whatever the steps it invokes write.
# 99 v4: After waiting for 07(chemistry) to complete, make the determination -> ref/chemistry_decision.json -> 08(requantification of 20 runs).
set -uo pipefail
W=$W; cd $W
PY=${APSC_PYTHON:-python}
L=$W/logs/99_driver.log
say(){ echo "[$(date -u +%FT%TZ)] [v4] $*" >> $L; echo "[$(date -u +%FT%TZ)] $*" >> $W/STATUS_APPEND.log; }
beat(){ date -u +%FT%TZ > $W/logs/driver_heartbeat.txt; }
beat; say "=== driver v4 起動 ==="

# Waiting for 07 to complete (abort if progress stalls for 3 hours)
while ! grep -q "07 DONE" $W/logs/07_chemistry.log 2>/dev/null; do
  beat; sleep 300
  age=$(( $(date +%s) - $(stat -f%m $W/logs/07_chemistry.log 2>/dev/null || echo 0) ))
  [ $age -ge 10800 ] && { say "07 の進捗が3時間停止 — 打ち切り"; break; }
done

if [ ! -s $W/ref/chemistry_decision.json ]; then
  $PY - <<'PY' >> $L 2>&1
import json, os
T = "$W/data/chem"
R = f"{W}/ref"
res = {}
for x in ("10xv2", "10xv3"):
    ri = f"{T}/out_{x}/run_info.json"
    bc = f"{T}/out_{x}/counts_unfiltered/cells_x_genes.barcodes.txt"
    if os.path.exists(ri):
        d = json.load(open(ri))
        res[x] = dict(p_pseudoaligned=d.get("p_pseudoaligned"), n_processed=d.get("n_processed"),
                      barcodes=sum(1 for _ in open(bc)) if os.path.exists(bc) else None)
if res:
    best = max(res, key=lambda k: (res[k]["barcodes"] or 0))
    out = dict(chemistry=best, empirical=res,
               geo_documented="10x Single Cell 3' Reagent Kits v2 (GSE181688/GSE197680 の extract_protocol に明記)",
               decided_by="ホワイトリスト一致バーコード数（同点なら p_pseudoaligned）", run="SRR15373534")
    json.dump(out, open(f"{R}/chemistry_decision.json", "w"), indent=1)
    print("chemistry decision:", json.dumps(out, ensure_ascii=False))
else:
    print("chemistry: 07 の出力が無い — 決定できない")
PY
fi
if [ -s $W/ref/chemistry_decision.json ]; then
  say "chemistry = $($PY -c "import json;print(json.load(open('$W/ref/chemistry_decision.json'))['chemistry'])")"
  git add -A >/dev/null 2>&1; git commit -qm "chemistry decision" >/dev/null 2>&1
  say "08 起動"
  bash $W/scripts/08_discovery_requant.sh
  say "08 終了（$(ls $W/data/counts_discovery/*/verify.json 2>/dev/null | wc -l | tr -d ' ')/5）"
else
  say "chemistry 未決定のため 08 をスキップ"
fi
beat; say "=== driver v4 終了 ==="
