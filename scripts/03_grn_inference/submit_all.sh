#!/usr/bin/env bash
# Submits one 201-task CellOracle array (plus a gather job) per patient not yet done.
# Usage: DATA_ROOT=clones_out RESULTS_ROOT=grn_results bash submit_all.sh

DATA_ROOT="${DATA_ROOT:-/scratch/${USER}/clones_out}"   # augmented patient h5ads
RESULTS_ROOT="${RESULTS_ROOT:-/scratch/${USER}/grn_results}"
BATCH="${BATCH:-12}"        # patients submitted per invocation (keep *201 < 5000)
INFLIGHT="${INFLIGHT:-6}"   # patients running concurrently (rolling dependency)
INNER="${INNER:-15}"        # replicate-tasks running at once within a patient (%)
# ---------------------------------------------------------------------------
set -uo pipefail
mkdir -p logs "$RESULTS_ROOT"

mapfile -t ALL < <(find "$DATA_ROOT" -name "*_clones.h5ad" | sort)
[ "${#ALL[@]}" -eq 0 ] && { echo "no *_clones.h5ad under $DATA_ROOT"; exit 1; }

# select patients not yet gathered
todo=()
for h in "${ALL[@]}"; do
  rel="${h#$DATA_ROOT/}"; res="$RESULTS_ROOT/${rel%_clones.h5ad}"
  [ -f "$res/.gathered" ] || todo+=("$h")
done
echo "patients: total=${#ALL[@]} remaining=${#todo[@]}  (launching up to $BATCH)"
[ "${#todo[@]}" -eq 0 ] && { echo "ALL DONE."; exit 0; }

declare -a gids
i=0
for h in "${todo[@]:0:$BATCH}"; do
  rel="${h#$DATA_ROOT/}"; res="$RESULTS_ROOT/${rel%_clones.h5ad}"; mkdir -p "$res"
  dep=""
  if [ "$i" -ge "$INFLIGHT" ]; then
    dep="--dependency=afterany:${gids[$((i-INFLIGHT))]}"   # gate on a patient INFLIGHT ago
  fi
  aid=$(sbatch $dep --parsable --array=0-200%"$INNER" celloracle_array.slurm "$h" "$res")
  gid=$(sbatch --parsable --dependency=afterany:"$aid" gather.slurm "$h" "$res")
  gids[$i]="$gid"
  echo "  [$((i+1))] $(basename "$h")  array=$aid gather=$gid ${dep:+(waits on ${gids[$((i-INFLIGHT))]})}"
  i=$((i+1))
done

remaining=$(( ${#todo[@]} - i ))
echo "submitted $i patient(s); $remaining still remaining -- re-run submit_all.sh to continue."
