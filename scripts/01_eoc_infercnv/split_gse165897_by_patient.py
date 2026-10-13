#!/usr/bin/env python3
"""Splits the GSE165897 count matrix and cell annotations into per-patient files.

Usage:
    python split_gse165897_by_patient.py --counts GSE165897_UMIcounts_HGSOC.tsv --annotations GSE165897_cellInfo_HGSOC.tsv --all --out_dir split_patients
"""

import argparse
import os
import pandas as pd
from datetime import datetime


def log(msg):
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}")


parser = argparse.ArgumentParser()
parser.add_argument("--counts", required=True, help="GSE165897_UMIcounts_HGSOC.tsv (all patients combined)")
parser.add_argument("--annotations", required=True, help="GSE165897_cellInfo_HGSOC.tsv (all patients combined)")
parser.add_argument("--patient_id", default=None, help="Single patient_id to extract (e.g. EOC372)")
parser.add_argument("--all", action="store_true", help="Extract every patient found in the annotation file")
parser.add_argument("--out_dir", required=True)
args = parser.parse_args()

if not args.patient_id and not args.all:
    raise SystemExit("Specify --patient_id <ID> or --all")

os.makedirs(args.out_dir, exist_ok=True)

# =============================================================================
# 1. Load the full annotation file (small -- ~52k rows, fine to load whole)
# =============================================================================
log(f"Loading annotations from {args.annotations}")
annot = pd.read_csv(args.annotations, sep="\t")
required_cols = {"cell", "patient_id", "cell_type"}
missing = required_cols - set(annot.columns)
if missing:
    raise SystemExit(f"Annotation file is missing expected column(s): {missing}")

target_patients = [args.patient_id] if args.patient_id else sorted(annot["patient_id"].unique())
log(f"Target patient(s): {len(target_patients)}")

log(f"Reading header of {args.counts} to map cell barcode -> column index...")
with open(args.counts) as f:
    header = f.readline().rstrip("\n").split("\t")

# First column is the (unnamed) gene-name index column -- always column 0,
# regardless of what it's actually labeled (blank in the real file).
gene_col_idx = 0
barcode_to_colidx = {barcode: i for i, barcode in enumerate(header) if i != gene_col_idx}
log(f"Count matrix has {len(barcode_to_colidx)} total cells across all patients")

for patient in target_patients:
    patient_annot = annot[annot["patient_id"] == patient]
    patient_cells = set(patient_annot["cell"])

    matched_cols = [barcode_to_colidx[c] for c in patient_cells if c in barcode_to_colidx]
    unmatched = patient_cells - set(barcode_to_colidx.keys())
    if unmatched:
        log(f"  [WARN] {patient}: {len(unmatched)} cells in annotation not found "
            f"in count matrix header (dropped): {list(unmatched)[:3]}...")

    if not matched_cols:
        log(f"  [SKIP] {patient}: no matching cells found in count matrix")
        continue

    usecols = [gene_col_idx] + sorted(matched_cols)
    log(f"  {patient}: {len(matched_cols)} cells -- reading matching columns from count matrix...")

    counts_subset = pd.read_csv(args.counts, sep="\t", usecols=usecols)
    # Rename the blank/unnamed gene-index column to something explicit
    counts_subset = counts_subset.rename(columns={counts_subset.columns[0]: "gene"})

    out_counts = os.path.join(args.out_dir, f"{patient}_counts.tsv")
    counts_subset.to_csv(out_counts, sep="\t", index=False)

    out_annot = os.path.join(args.out_dir, f"{patient}_annotations.tsv")
    patient_annot[["cell", "cell_type"]].to_csv(out_annot, sep="\t", index=False, header=False)

    log(f"    -> {out_counts}  ({counts_subset.shape[0]} genes x {counts_subset.shape[1]-1} cells)")
    log(f"    -> {out_annot}")

log("Done.")
