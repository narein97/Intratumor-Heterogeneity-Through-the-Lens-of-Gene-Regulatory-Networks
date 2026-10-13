#!/usr/bin/env bash
# Runs compute_bootstrap_support.py over all patients as a SLURM array.
# Usage: DATA_ROOT=clones_out bash submit_bootstrap_support.sh

set -uo pipefail

DATA_ROOT="${DATA_ROOT:-/scratch/${USER}/clones_out}"
FILE_LIST="${FILE_LIST:-patients_to_process.txt}"
CONCURRENCY="${CONCURRENCY:-50}"   # concurrent array tasks -- generous given
                                   # each task is 1 CPU / 4GB / ~30s, well
                                   # under the 1536-core / 2304GB QOS caps
                                   # described in GRN_pipeline_summary.md

mkdir -p logs

echo "Discovering patients under $DATA_ROOT (excluding Ovarian)..."
python compute_bootstrap_support.py --root "$DATA_ROOT" --write-file-list "$FILE_LIST"

N=$(wc -l < "$FILE_LIST")
if [ "$N" -eq 0 ]; then
    echo "No patients found under $DATA_ROOT -- nothing to submit."
    exit 0
fi

LAST=$((N - 1))
echo "Submitting array 0-$LAST%${CONCURRENCY} ($N patients, up to $CONCURRENCY concurrent)..."
sbatch --array="0-${LAST}%${CONCURRENCY}" bootstrap_support_array.slurm "$FILE_LIST"
