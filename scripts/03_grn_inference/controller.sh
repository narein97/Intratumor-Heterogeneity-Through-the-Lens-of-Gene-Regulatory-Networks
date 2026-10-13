#!/usr/bin/env bash
# Loops submit_all.sh and waits for the queue to drain until every patient is gathered.
# Usage: DATA_ROOT=clones_out RESULTS_ROOT=grn_results nohup bash controller.sh > controller.log 2>&1 &

set -uo pipefail

export DATA_ROOT="${DATA_ROOT:-/scratch/${USER}/clones_out}"
export RESULTS_ROOT="${RESULTS_ROOT:-/scratch/${USER}/grn_results}"
export BATCH="${BATCH:-12}"
export INFLIGHT="${INFLIGHT:-6}"
export INNER="${INNER:-15}"
INTERVAL="${INTERVAL:-300}"   # seconds between queue checks (5 min)

# jobs of THIS pipeline still queued/running; prints -1 if scheduler unreachable
myjobs() {
    local out rc
    out="$(squeue --me -h -o '%j' 2>/dev/null)"; rc=$?
    [ "$rc" -ne 0 ] && { echo -1; return; }
    printf '%s\n' "$out" | grep -cE '^(co_grn|co_gather)$'
}
hms() { printf '%dh%02dm' $(( $1/3600 )) $(( ($1%3600)/60 )); }

echo "[controller] start $(date)  BATCH=$BATCH INFLIGHT=$INFLIGHT INNER=$INNER INTERVAL=${INTERVAL}s"
wave=0
while true; do
    out="$(bash submit_all.sh)"
    echo "$out"
    grep -q "ALL DONE" <<<"$out" && { echo "[controller] ALL DONE at $(date)"; break; }
    wave=$((wave + 1))
    wave_start=$(date +%s)
    echo "[controller] wave $wave launched at $(date); polling every ${INTERVAL}s until queue drains ..."

    # poll until this pipeline's queue is empty -- NEVER give up, just keep reporting
    while true; do
        sleep "$INTERVAL"
        n="$(myjobs)"
        elapsed=$(( $(date +%s) - wave_start ))
        g="$(find "$RESULTS_ROOT" -name .gathered 2>/dev/null | wc -l)"
        if [ "$n" -lt 0 ]; then
            echo "[controller] wave $wave: squeue unreachable (scheduler blip); $(hms $elapsed) in, gathered=$g -- retrying"
        elif [ "$n" -eq 0 ]; then
            echo "[controller] wave $wave DRAINED after $(hms $elapsed); gathered=$g/257"
            break
        else
            echo "[controller] wave $wave: $n job(s) still in queue; $(hms $elapsed) elapsed, gathered=$g/257"
        fi
    done
done
echo "[controller] exiting $(date)"
