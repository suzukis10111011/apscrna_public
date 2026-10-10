"""Purpose: Guard common analysis operations with thread limits and process heartbeat checks
Input: none (monitors logs/requant_heartbeat.txt and output directory growth)
Output: none (enforces constraints, no persistent output)

Common guard on the analysis side: thread limits and liveness checks for requantification.

Lesson from 2026-09-23: the old version only watched updates to logs/requant_heartbeat.txt, but since
08 only beats at stage boundaries, during a 65-minute .sra fetch or a multi-hour fastq-dump it appears
"stopped." In fact we misjudged a live 08 as dead this way (ps is unreliable in this sandbox).
So even if the heartbeat is stale, if the requantification output directory is growing, we treat it as alive.
"""
import os, time
for v in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","NUMEXPR_NUM_THREADS","NUMBA_NUM_THREADS"):
    os.environ[v] = "4"
os.environ.setdefault("NUMBA_CACHE_DIR", os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), ".numba_cache"))
W = os.environ.get("APSC_ROOT") or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WATCH = [f"{W}/data/dsra", f"{W}/data/dfq", f"{W}/data/counts_discovery"]


def _bytes():
    t = 0
    for d in WATCH:
        for root, _, fs in os.walk(d):
            for f in fs:
                try: t += os.path.getsize(os.path.join(root, f))
                except OSError: pass
    return t


def requant_complete(n_expect=5):
    """08 が完走したか。verify.json が n_expect 個あり、08 のログに DONE 行があれば完了。"""
    import glob, re
    n = len(glob.glob(f"{W}/data/counts_discovery/Sample*/verify.json"))
    lg = f"{W}/logs/08_requant.log"
    done = os.path.exists(lg) and re.search(r"08 DONE counted=(\d+)/", open(lg, errors="replace").read()) is not None
    return n >= n_expect and done


def requant_state(limit_min=10, grow_s=45, quiet=False):
    """'complete'（08 完走）/ 'absent'（未起動）/ 'beating' / 'growing'（heartbeat は古いが
    出力が増えている）/ 'stalled'。

    2026-09-26 の教訓: 完了と停止を区別していなかったため、5検体を完走した後の解析が
    「停止」として中断された。完了判定を先に置く。"""
    if requant_complete():
        return "complete"
    p = f"{W}/logs/requant_heartbeat.txt"
    if not os.path.exists(p) and not any(os.path.exists(d) for d in WATCH):
        return "absent"
    if os.path.exists(p) and (time.time() - os.path.getmtime(p)) / 60 <= limit_min:
        return "beating"
    a = _bytes(); time.sleep(grow_s); b = _bytes()
    if not quiet:
        print(f"[guard] heartbeat は古いが出力を {grow_s}s 監視: {(b-a)/1e6:+.1f} MB")
    return "growing" if b > a else "stalled"


def require_heartbeat(step):
    s = requant_state()
    if s in ("complete", "absent", "beating", "growing"):
        print(f"[guard] 再定量 = {s} — {step} を続行（4スレッド制限）")
        return
    raise SystemExit(f"[guard] 再定量が停止（heartbeat 古く出力も増えない）— {step} を中断（ドライバは再起動しない）")
