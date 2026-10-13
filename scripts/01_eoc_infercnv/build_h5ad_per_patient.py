#!/usr/bin/env python3
"""Builds one h5ad per EOC patient (malignant + stromal cells, inferCNV copy-number profiles in .obsm['cnp']).

Usage:
    python build_h5ad_per_patient.py --split_dir split_patients --full_annotations GSE165897_cellInfo_HGSOC.tsv --infercnv_dir infercnv_out --all --out_dir h5ad_out
"""

import argparse
import os
import glob
import pandas as pd
import numpy as np
import anndata as ad
from scipy import sparse
from datetime import datetime


def log(msg):
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}")


parser = argparse.ArgumentParser()
parser.add_argument("--counts", default=None, help="<patient>_counts.tsv (single-patient mode)")
parser.add_argument("--split_dir", default=None, help="Directory of split per-patient counts files (--all mode)")
parser.add_argument("--full_annotations", required=True, help="GSE165897_cellInfo_HGSOC.tsv (ALL patients, all columns)")
parser.add_argument("--patient_id", default=None)
parser.add_argument("--all", action="store_true")
parser.add_argument("--infercnv_dir", required=True)
parser.add_argument("--cnv_filename", default="infercnv.observations.txt",
                   help="filename (within each patient's infercnv_dir) to use as the "
                        "CN profile source. REASONED CHOICE, not independently confirmed "
                        "for every inferCNV version -- see add_cnp()'s docstring below.")
parser.add_argument("--out_dir", required=True)
args = parser.parse_args()

if not args.all and not (args.counts and args.patient_id):
    raise SystemExit("Specify (--counts <file> --patient_id <ID>) for one patient, or --split_dir ... --all")

os.makedirs(args.out_dir, exist_ok=True)

log(f"Loading full annotation file {args.full_annotations}...")
full_annot = pd.read_csv(args.full_annotations, sep="\t")
required_cols = {"cell", "patient_id", "cell_type"}
missing = required_cols - set(full_annot.columns)
if missing:
    raise SystemExit(f"Full annotation file missing expected column(s): {missing}")

if args.all:
    if not args.split_dir:
        raise SystemExit("--all requires --split_dir")
    counts_files = sorted(glob.glob(os.path.join(args.split_dir, "*_counts.tsv")))
    jobs = [(f, os.path.basename(f).replace("_counts.tsv", "")) for f in counts_files]
else:
    jobs = [(args.counts, args.patient_id)]

log(f"Target patient(s): {len(jobs)}")


def load_cnv_observations(path):
    """
    Load inferCNV's observations output, auto-orienting (genes-as-rows,
    cells-as-columns is inferCNV's standard convention; transposed to
    cells x genes here to match .obsm's orientation).

    MOVED here from update_h5ad_cnp_bootstrap_shuffle.py: CNP must exist
    BEFORE clone_sweep_zhang2022.py runs (it's an input to Leiden
    clustering), whereas bootstrap/shuffle can only meaningfully be built
    AFTER the final clone assignment is known -- these two used to be
    bundled into one later step, which was backwards for anyone re-deriving
    clones via clone_sweep_zhang2022.py rather than trusting inferCNV's
    raw subclone output as final.
    """
    df = pd.read_csv(path, sep=r"\s+", index_col=0)
    return df.T


def add_cnp(adata, infercnv_dir, cnv_filename):
    """
    REASONED CHOICE, not independently confirmed: given denoise=TRUE and
    HMM=TRUE were both used, and 3CA's own CNA_matrix (which clone_sweep.py
    consumes) is a continuous signal, this defaults to
    infercnv.observations.txt. On newer inferCNV versions this file may not
    exist at all (confirmed absent in at least one real run directory under
    1.18.x) -- if so, this logs a warning and skips CNP for that patient
    rather than silently producing an empty/wrong result. Override via
    --cnv_filename if you've identified the right file for your version.
    """
    cnv_path = os.path.join(infercnv_dir, cnv_filename)
    if not os.path.exists(cnv_path):
        log(f"  [SKIP cnp] {cnv_path} not found")
        return False

    cnv = load_cnv_observations(cnv_path)  # cells x genes
    clone = adata.obs["clone"].astype(str).to_numpy()
    malignant_cells = adata.obs_names[clone != "N"]

    shared_genes = sorted(set(adata.var_names) & set(cnv.columns))
    if not shared_genes:
        log(f"  [SKIP cnp] no genes shared between {cnv_path} and adata.var_names")
        return False

    shared_malignant = [c for c in malignant_cells if c in cnv.index]
    missing_malignant = set(malignant_cells) - set(cnv.index)
    if missing_malignant:
        log(f"  [WARN cnp] {len(missing_malignant)} malignant cells have no entry "
            f"in {cnv_filename} -- left as NaN like non-malignant cells: "
            f"{list(missing_malignant)[:3]}")

    cnp_df = pd.DataFrame(np.nan, index=adata.obs_names, columns=shared_genes, dtype="float32")
    cnp_df.loc[shared_malignant, :] = cnv.loc[shared_malignant, shared_genes].to_numpy().astype(np.float32)

    adata.obsm["cnp"] = cnp_df.to_numpy()
    adata.uns["cnp_genes"] = list(shared_genes)
    log(f"  [ok cnp] {len(shared_genes)} genes, {len(shared_malignant)}/{len(malignant_cells)} "
        f"malignant cells with real values")
    return True


def build_one_patient(counts_path, patient):
    subclone_path = os.path.join(args.infercnv_dir, patient, "subclone_assignments_filtered.tsv")

    if not os.path.exists(counts_path):
        log(f"  [SKIP] {patient} -- counts file not found: {counts_path}")
        return

    log(f"=== {patient} ===")

    # ---- 1. Expression matrix: genes (rows) x cells (columns) -> transpose
    log("  Loading counts matrix...")
    counts = pd.read_csv(counts_path, sep="\t", index_col=0)
    expr = counts.T
    log(f"    {counts.shape[0]} genes x {counts.shape[1]} cells (as stored)")

    # ---- 2. Filter the FULL annotation file to this patient's cells,
    #         using an exact match on patient_id (never substring-matching
    #         the barcode -- e.g. "EOC3" is a literal substring of "EOC372")
    patient_annot = full_annot[full_annot["patient_id"] == patient].set_index("cell")

    missing_from_annot = set(expr.index) - set(patient_annot.index)
    if missing_from_annot:
        log(f"  [WARN] {len(missing_from_annot)} cells in the count matrix have "
            f"no entry in the full annotation file for patient_id=={patient} "
            f"-- these cells will be dropped: {list(missing_from_annot)[:3]}...")

    # ---- 3. Load inferCNV subclone assignment (may be entirely absent if
    #         inferCNV hasn't been run for this patient yet) ------------------
    if os.path.exists(subclone_path):
        subclone = pd.read_csv(subclone_path, sep="\t").set_index("cell")
        subclone_map = subclone["subclone"]
    else:
        log(f"  [WARN] No subclone file at {subclone_path} -- all EOC cells "
            f"will have a missing clone value")
        subclone_map = pd.Series(dtype=object)

    # ---- 4. Build the clone column, explicitly gated on cell_type --------
    def assign_clone(cell, cell_type):
        if cell_type == "Stromal":
            return "N"
        elif cell_type == "EOC":
            return subclone_map.get(cell, np.nan)
        else:
            return None  # Immune -- irrelevant, these rows get dropped next anyway

    patient_annot = patient_annot.copy()
    patient_annot["clone"] = [
        assign_clone(cell, ct) for cell, ct in zip(patient_annot.index, patient_annot["cell_type"])
    ]

    # ---- 5. Drop Immune cells, cells missing from the annotation file, AND
    #         EOC cells with no inferCNV subclone assignment (rather than
    #         keeping them with clone=NaN) -----------------------------------
    keep_mask = (patient_annot["cell_type"] != "Immune") & (patient_annot["clone"].notna())
    keep_cells = patient_annot.index[keep_mask]
    keep_cells = [c for c in keep_cells if c in expr.index]

    n_immune_dropped = (patient_annot["cell_type"] == "Immune").sum()
    n_eoc_no_clone = ((patient_annot["cell_type"] == "EOC") & (patient_annot["clone"].isna())).sum()
    n_stromal = (patient_annot["cell_type"] == "Stromal").sum()
    n_eoc = (patient_annot["cell_type"] == "EOC").sum()

    log(f"  Cell type breakdown before filtering: EOC={n_eoc}, Stromal={n_stromal}, "
        f"Immune={n_immune_dropped} (dropped)")
    if n_eoc_no_clone > 0:
        log(f"  [NOTE] {n_eoc_no_clone} EOC cells have no inferCNV subclone "
            f"assignment and are DROPPED (not kept with a missing clone value)")

    expr_filtered = expr.loc[keep_cells]
    obs = patient_annot.loc[keep_cells]

    log(f"  Final cell count after removing Immune: {len(keep_cells)}")

    # ---- 6. Build the AnnData object ---------------------------------------
    X = sparse.csr_matrix(expr_filtered.values.astype(np.float32))
    var = pd.DataFrame(index=expr_filtered.columns.tolist())

    adata = ad.AnnData(X=X, obs=obs.reset_index().rename(columns={"index": "cell_barcode"}), var=var)
    adata.obs_names = obs.index.tolist()  # keep the real cell barcode as the obs index
    adata.var_names_make_unique()

    # ---- 7. CN profile -- must be present before clone_sweep_zhang2022.py
    #         can run, since it's an input to Leiden clustering there ------
    add_cnp(adata, os.path.join(args.infercnv_dir, patient), args.cnv_filename)

    out_path = os.path.join(args.out_dir, f"{patient}.h5ad")
    adata.write_h5ad(out_path)
    log(f"  -> {out_path}  ({adata.n_obs} cells x {adata.n_vars} genes)")


for counts_path, patient in jobs:
    build_one_patient(counts_path, patient)

log("Done.")
