#!/usr/bin/env bash
# Runs run_infercnv_zhang2022.R on every patient in split_patients/, several at a time (largest first, resumable).
# Usage: MAX_PARALLEL=10 bash run_infercnv_all_patients_parallel.sh

set -uo pipefail   # deliberately no -e

SPLIT_DIR="split_patients"
OUT_DIR="infercnv_out"
GENE_ORDER_FILE="gencode_v25_gene_order.txt"
STROMAL_LABELS="Stromal"
MALIGNANT_LABELS="EOC"
SEED=1
LOG_DIR="logs/infercnv"
MAX_PARALLEL="${MAX_PARALLEL:-10}"   # 100 cores / num_threads=10 per job

log() { echo "[$(date '+%H:%M:%S')] $*"; }

mkdir -p "$LOG_DIR"

[[ -d "$SPLIT_DIR" ]] || { log "[ERROR] $SPLIT_DIR not found -- run split_gse165897_by_patient.py first"; exit 1; }
[[ -f "$GENE_ORDER_FILE" ]] || { log "[ERROR] $GENE_ORDER_FILE not found"; exit 1; }

log "MAX_PARALLEL=$MAX_PARALLEL (each inferCNV run uses num_threads=10 internally --"
log "  so $MAX_PARALLEL concurrent jobs uses up to $((MAX_PARALLEL * 10)) cores)"

# ---- Discover patients and their cell counts (header field count - 1),
#      sort largest-first ------------------------------------------------
declare -A PATIENT_SIZE
for f in "$SPLIT_DIR"/*_counts.tsv; do
    [[ -e "$f" ]] || continue
    PATIENT=$(basename "$f" | sed 's/_counts\.tsv$//')
    N_FIELDS=$(head -1 "$f" | awk -F'\t' '{print NF}')
    PATIENT_SIZE["$PATIENT"]=$(( N_FIELDS - 1 ))  # -1 for the gene-name column
done

if [[ ${#PATIENT_SIZE[@]} -eq 0 ]]; then
    log "[ERROR] No *_counts.tsv files found in $SPLIT_DIR"
    exit 1
fi

mapfile -t PATIENTS < <(
    for p in "${!PATIENT_SIZE[@]}"; do
        echo "${PATIENT_SIZE[$p]} $p"
    done | sort -rn -k1,1 | awk '{print $2}'
)

log "Found ${#PATIENTS[@]} patient(s), largest-first: ${PATIENTS[*]}"

# ---- Job-slot-limited launch ---------------------------------------------
run_one_patient() {
    local PATIENT="$1"
    local COMPLETE_MARKER="${OUT_DIR}/${PATIENT}/subclone_assignments_filtered.tsv"
    local LOG="${LOG_DIR}/${PATIENT}.log"

    Rscript run_infercnv_zhang2022.R \
        --counts        "${SPLIT_DIR}/${PATIENT}_counts.tsv" \
        --annotations   "${SPLIT_DIR}/${PATIENT}_annotations.tsv" \
        --gene_order_file "$GENE_ORDER_FILE" \
        --stromal_labels "$STROMAL_LABELS" \
        --malignant_labels "$MALIGNANT_LABELS" \
        --out_dir       "${OUT_DIR}/${PATIENT}" \
        --seed          "$SEED" \
        > "$LOG" 2>&1
}

N_LAUNCHED=0
N_SKIPPED=0
declare -a TO_RUN=()

for PATIENT in "${PATIENTS[@]}"; do
    COMPLETE_MARKER="${OUT_DIR}/${PATIENT}/subclone_assignments_filtered.tsv"
    if [[ -f "$COMPLETE_MARKER" ]]; then
        log "  [SKIP] $PATIENT -- already completed"
        N_SKIPPED=$(( N_SKIPPED + 1 ))
    else
        TO_RUN+=("$PATIENT")
    fi
done

log "${#TO_RUN[@]} patient(s) to run, $N_SKIPPED already complete"

N_ACTIVE=0
for PATIENT in "${TO_RUN[@]}"; do
    if [[ $N_ACTIVE -ge $MAX_PARALLEL ]]; then
        wait -n   # free up a slot as soon as ANY running job finishes
        N_ACTIVE=$(( N_ACTIVE - 1 ))
    fi
    log "  [LAUNCH] $PATIENT (${PATIENT_SIZE[$PATIENT]} cells) -- log: ${LOG_DIR}/${PATIENT}.log  [$((N_ACTIVE+1))/${MAX_PARALLEL} slots]"
    run_one_patient "$PATIENT" &
    N_ACTIVE=$(( N_ACTIVE + 1 ))
    N_LAUNCHED=$(( N_LAUNCHED + 1 ))
done

log "All ${#TO_RUN[@]} job(s) launched -- waiting for the remaining ones to finish..."
wait   # final, unconditional wait for everything still running

N_OK=0; N_FAIL=0
FAILED_PATIENTS=()
for PATIENT in "${TO_RUN[@]}"; do
    if [[ -f "${OUT_DIR}/${PATIENT}/subclone_assignments_filtered.tsv" ]]; then
        N_OK=$(( N_OK + 1 ))
    else
        N_FAIL=$(( N_FAIL + 1 ))
        FAILED_PATIENTS+=("$PATIENT")
    fi
done

log ""
log "════════════════════════════════════════"
log "  Completed this run : $N_OK"
log "  Already done before: $N_SKIPPED"
log "  Failed             : $N_FAIL"
[[ $N_FAIL -gt 0 ]] && log "  Failed patients: ${FAILED_PATIENTS[*]} -- check logs in $LOG_DIR/"
log "════════════════════════════════════════"

[[ $N_FAIL -eq 0 ]] && exit 0 || exit 1
