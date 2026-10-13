#!/usr/bin/env python3
"""Folds one patient's per-replicate CellOracle GRNs into its h5ad (.uns['celloracle_*']).

Usage:
    python gather_grn.py --h5ad P_clones.h5ad --results-dir DIR
"""

import os
import argparse
import logging
import numpy as np

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s %(levelname)-7s %(message)s",
                    datefmt="%H:%M:%S")
log = logging.getLogger("co_gather")


def unit_path(results_dir, kind, rep):
    name = "original" if kind == "original" else "%s_%03d" % (kind, rep)
    return os.path.join(results_dir, kind, name + ".npz")


def load_units(results_dir, n_boot, n_shuffle):
    """Return dict name -> loaded npz (or None if missing), and value_cols."""
    units = {"original": (unit_path(results_dir, "original", -1))}
    for i in range(n_boot):
        units["bootstrap_%03d" % i] = unit_path(results_dir, "bootstrap", i)
    for i in range(n_shuffle):
        units["shuffle_%03d" % i] = unit_path(results_dir, "shuffle", i)
    loaded, value_cols = {}, None
    for name, path in units.items():
        if os.path.exists(path):
            try:
                z = np.load(path, allow_pickle=False)
                loaded[name] = {k: z[k] for k in z.files}
                if value_cols is None:
                    value_cols = [str(c) for c in z["value_cols"]]
            except Exception as e:
                log.warning("corrupt unit %s: %s", path, e)
                loaded[name] = None
        else:
            loaded[name] = None
    return loaded, value_cols


def encode(unit, gene_code, clone_code):
    """Re-code a unit's edgelist to global gene/clone codes -> float32 2-D array."""
    genes = unit["genes"].astype(str)
    clones = unit["clones"].astype(str)
    e = unit["edges"]                       # [local_clone, local_src, local_tgt]
    g2 = np.array([gene_code[g] for g in genes])
    c2 = np.array([clone_code.get(c, -1) for c in clones])
    clone_col = c2[e[:, 0]].astype(np.float32)
    src_col = g2[e[:, 1]].astype(np.float32)
    tgt_col = g2[e[:, 2]].astype(np.float32)
    return np.column_stack([clone_col, src_col, tgt_col,
                            unit["values"].astype(np.float32)])


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--h5ad", required=True)
    p.add_argument("--results-dir", required=True)
    p.add_argument("--n-boot", type=int, default=100)
    p.add_argument("--n-shuffle", type=int, default=100)
    p.add_argument("--force", action="store_true")
    args = p.parse_args()

    import anndata as ad
    adata = ad.read_h5ad(args.h5ad)
    if "celloracle_original" in adata.uns and not args.force:
        log.info("[skip] already gathered: %s", args.h5ad)
        return

    loaded, value_cols = load_units(args.results_dir, args.n_boot, args.n_shuffle)
    present = {k: v for k, v in loaded.items() if v is not None}
    missing = [k for k, v in loaded.items() if v is None]
    if not present:
        log.warning("no unit checkpoints found under %s -- nothing to gather",
                    args.results_dir)
        return
    if missing:
        log.warning("%d/%d units missing (folding the %d present): %s ...",
                    len(missing), len(loaded), len(present), missing[:5])

    # unified gene vocabulary across all present units
    all_genes = sorted(set().union(*[set(u["genes"].astype(str)) for u in present.values()]))
    gene_code = {g: i for i, g in enumerate(all_genes)}
    clone_code = {str(c): i for i, c in
                  enumerate(np.asarray(adata.uns["clone_categories"]).astype(str))}

    adata.uns["celloracle_genes"] = np.asarray(all_genes, dtype=object)
    adata.uns["celloracle_columns"] = np.asarray(
        ["clone", "source", "target"] + (value_cols or []), dtype=object)

    if loaded["original"] is not None:
        adata.uns["celloracle_original"] = encode(loaded["original"], gene_code, clone_code)
    boot = {str(i): encode(loaded["bootstrap_%03d" % i], gene_code, clone_code)
            for i in range(args.n_boot) if loaded["bootstrap_%03d" % i] is not None}
    shuf = {str(i): encode(loaded["shuffle_%03d" % i], gene_code, clone_code)
            for i in range(args.n_shuffle) if loaded["shuffle_%03d" % i] is not None}
    adata.uns["celloracle_bootstrap"] = boot
    adata.uns["celloracle_shuffle"] = shuf
    if missing:
        adata.uns["celloracle_missing"] = np.asarray(missing, dtype=object)

    tmp = args.h5ad + ".tmp.h5ad"
    adata.write_h5ad(tmp)
    os.replace(tmp, args.h5ad)              # atomic swap; never leave a half-written h5ad
    log.info("[ok] gathered %s: original=%s boot=%d shuffle=%d genes=%d%s",
             os.path.basename(args.h5ad),
             "celloracle_original" in adata.uns, len(boot), len(shuf),
             len(all_genes), (" MISSING=%d" % len(missing)) if missing else "")


if __name__ == "__main__":
    main()
