#!/usr/bin/env python3
"""Keeps significant (p < 0.05) CellOracle edges and counts in how many of the 100 bootstraps each is recovered.

Usage:
    python compute_bootstrap_support.py --root clones_out
"""

import argparse
import glob
import os
import logging

import numpy as np
import pandas as pd
import anndata as ad

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)-7s %(message)s",
                    datefmt="%H:%M:%S")
log = logging.getLogger("bootstrap_support")


def output_path(source_path, suffix="_with_support"):
    base, ext = os.path.splitext(source_path)
    return f"{base}{suffix}{ext}"


def decode_edgelist(arr, columns, genes, clone_categories):
    df = pd.DataFrame(arr, columns=columns)
    df["clone"] = clone_categories[df["clone"].astype(int)]
    df["source"] = genes[df["source"].astype(int)]
    df["target"] = genes[df["target"].astype(int)]
    return df


def compute_support_for_patient(source_path, p_thresh=0.05, force=False):
    out_path = output_path(source_path)

    if not force and os.path.exists(out_path):
        log.info("  [skip] %s -- output already exists (%s)", source_path, out_path)
        return None

    # Read-only: this variable is never passed to any write call, anywhere
    # in this function, under any code path.
    adata = ad.read_h5ad(source_path)

    required = ["celloracle_original", "celloracle_bootstrap", "celloracle_columns",
                "celloracle_genes", "clone_categories"]
    missing = [k for k in required if k not in adata.uns]
    if missing:
        log.warning("  [skip] %s -- missing .uns key(s): %s", source_path, missing)
        return None

    columns = list(adata.uns["celloracle_columns"])
    genes = np.asarray(adata.uns["celloracle_genes"]).astype(str)
    clone_categories = np.asarray(adata.uns["clone_categories"]).astype(str)

    orig_df = decode_edgelist(adata.uns["celloracle_original"], columns, genes, clone_categories)
    orig_sig = orig_df[orig_df["p"] < p_thresh].copy()

    boot_edge_sets = []
    for rep_key, rep_arr in adata.uns["celloracle_bootstrap"].items():
        rep_df = decode_edgelist(rep_arr, columns, genes, clone_categories)
        rep_sig = rep_df[rep_df["p"] < p_thresh]
        boot_edge_sets.append(set(zip(rep_sig["clone"], rep_sig["source"], rep_sig["target"])))

    n_boot = len(boot_edge_sets)

    def support(row):
        key = (row["clone"], row["source"], row["target"])
        return sum(1 for s in boot_edge_sets if key in s)

    orig_sig["bootstrap_support"] = orig_sig.apply(support, axis=1)

    log.info("  %s: %d/%d edges significant (p<%.2f), %d bootstrap replicates found",
             source_path, len(orig_sig), len(orig_df), p_thresh, n_boot)

    adata.uns["celloracle_original_filtered"] = orig_sig.reset_index(drop=True)
    tmp_path = out_path + ".tmp.h5ad"
    adata.write_h5ad(tmp_path)
    os.replace(tmp_path, out_path)
    log.info("  -> %s (source untouched)", out_path)
    return orig_sig


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--h5ad", default=None)
    p.add_argument("--root", default=None, help="walk this dir for *_clones.h5ad")
    p.add_argument("--p-thresh", type=float, default=0.05)
    p.add_argument("--exclude-cancer-type", nargs="*", default=["Ovarian"])
    p.add_argument("--force", action="store_true",
                   help="recompute even if the _with_support output already exists")
    p.add_argument("--write-file-list", default=None)
    p.add_argument("--file-list", default=None)
    p.add_argument("--array-index", type=int, default=None)
    args = p.parse_args()

    def discover():
        files = sorted(glob.glob(f"{args.root}/**/*_clones.h5ad", recursive=True))
        files = [f for f in files if "_with_support" not in f]  # never re-glob our own output
        if args.exclude_cancer_type:
            before = len(files)
            files = [f for f in files if not any(f"/{ct}/" in f for ct in args.exclude_cancer_type)]
            log.info("Found %d patient h5ad(s), excluded %d matching %s -> %d remaining",
                     before, before - len(files), args.exclude_cancer_type, len(files))
        return files

    if args.write_file_list:
        if not args.root:
            p.error("--write-file-list requires --root")
        files = discover()
        with open(args.write_file_list, "w") as f:
            f.write("\n".join(files) + "\n")
        log.info("Wrote %d path(s) -> %s (submit an array job with --array=0-%d)",
                 len(files), args.write_file_list, len(files) - 1)
        return

    if args.file_list is not None and args.array_index is not None:
        with open(args.file_list) as f:
            all_files = [line.strip() for line in f if line.strip()]
        if args.array_index >= len(all_files):
            log.warning("array-index %d >= %d files -- nothing to do", args.array_index, len(all_files))
            return
        files = [all_files[args.array_index]]
    elif args.h5ad:
        files = [args.h5ad]
    elif args.root:
        files = discover()
    else:
        p.error("Provide --h5ad <file>, --root <dir>, --write-file-list, or --file-list + --array-index")
        return

    for f in files:
        compute_support_for_patient(f, args.p_thresh, args.force)


if __name__ == "__main__":
    main()
