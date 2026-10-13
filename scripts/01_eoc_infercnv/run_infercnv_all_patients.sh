#!/usr/bin/env bash
# =============================================================================
# run_infercnv_all_patients.sh
#
# Runs run_infercnv_zhang2022.R sequentially across every patient found in
# split_patients/, skipping any patient that's already completed.
#
# Sequential by design, not parallel: inferCNV is already configured to use
# num_threads=10 internally (hardcoded in run_infercnv_zhang2022.R to match
# the paper's own parameters) -- running several inferCNV instances
# concurrently would oversubscribe CPU rather than speed things up, unless
# you have an unusually large machine. If you do, and want to parallelize
# across patients, reduce num_threads in run_infercnv_zhang2022.R first and
# then add job-slot parallelism -- not done here to keep this simple.
#
# Resume: a patient is considered complete if
#   infercnv_out/<patient>/subclone_assignments_filtered.tsv
# already exists (the final artifact run_infercnv_zhang2022.R produces on
# a successful run).
#
# One patient failing does not stop the batch -- it's logged and the script
# moves on to the next patient, with a final pass/fail summary at the end.
#
# Usage:
#   bash run_infercnv_all_patients.sh
# =============================================================================

set -uo pipefail   # NOTE: deliberately no -e -- a failing Rscript call for
                    # one patient must not kill the whole batch; each call
                    # is checked explicitly instead (see run_one_patient).

SPLIT_DIR="split_patients"
OUT_DIR="infercnv_out"
GENE_ORDER_FILE="gencode_v25_gene_order.txt"
STROMAL_LABELS="Stromal"
SEED=1
LOG_DIR="logs/infercnv"

log() { echo "[$(date '+%H:%M:%S')] $*"; }

mkdir -p "$LOG_DIR"

[[ -d "$SPLIT_DIR" ]] || { log "[ERROR] $SPLIT_DIR not found -- run split_gse165897_by_patient.py first"; exit 1; }
[[ -f "$GENE_ORDER_FILE" ]] || { log "[ERROR] $GENE_ORDER_FILE not found"; exit 1; }

# Discover patients from whatever's actually been split, rather than
# re-parsing the combined annotation file -- this naturally only includes
# patients we actually have per-patient input files for.
mapfile -t PATIENTS < <(
    find "$SPLIT_DIR" -maxdepth 1 -name "*_counts.tsv" -exec basename {} \; \
        | sed 's/_counts\.tsv$//' | sort
)

if [[ ${#PATIENTS[@]} -eq 0 ]]; then
    log "[ERROR] No *_counts.tsv files found in $SPLIT_DIR"
    exit 1
fi

log "Found ${#PATIENTS[@]} patient(s): ${PATIENTS[*]}"

N_OK=0; N_SKIP=0; N_FAIL=0
FAILED_PATIENTS=()

run_one_patient() {
    local PATIENT="$1"
    local COMPLETE_MARKER="${OUT_DIR}/${PATIENT}/subclone_assignments_filtered.tsv"
    local LOG="${LOG_DIR}/${PATIENT}.log"

    if [[ -f "$COMPLETE_MARKER" ]]; then
        log "  [SKIP] $PATIENT -- already completed ($COMPLETE_MARKER exists)"
        N_SKIP=$(( N_SKIP + 1 ))
        return 0
    fi

    log "  [RUN]  $PATIENT -- log: $LOG"
    Rscript run_infercnv_zhang2022.R \
        --counts        "${SPLIT_DIR}/${PATIENT}_counts.tsv" \
        --annotations   "${SPLIT_DIR}/${PATIENT}_annotations.tsv" \
        --gene_order_file "$GENE_ORDER_FILE" \
        --stromal_labels "$STROMAL_LABELS" \
        --out_dir       "${OUT_DIR}/${PATIENT}" \
        --seed          "$SEED" \
        > "$LOG" 2>&1

    if [[ $? -eq 0 && -f "$COMPLETE_MARKER" ]]; then
        log "  [OK]   $PATIENT"
        N_OK=$(( N_OK + 1 ))
    else
        log "  [FAIL] $PATIENT -- see $LOG"
        N_FAIL=$(( N_FAIL + 1 ))
        FAILED_PATIENTS+=("$PATIENT")
    fi
}

for PATIENT in "${PATIENTS[@]}"; do
    run_one_patient "$PATIENT"
done

log ""
log "════════════════════════════════════════"
log "  Completed : $N_OK"
log "  Skipped   : $N_SKIP"
log "  Failed    : $N_FAIL"
[[ $N_FAIL -gt 0 ]] && log "  Failed patients: ${FAILED_PATIENTS[*]}"
log "════════════════════════════════════════"

[[ $N_FAIL -eq 0 ]] && exit 0 || exit 1
