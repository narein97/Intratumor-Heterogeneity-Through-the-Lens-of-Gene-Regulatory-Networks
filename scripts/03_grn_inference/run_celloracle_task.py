#!/usr/bin/env python3
"""Runs CellOracle on one unit (original, one bootstrap or one shuffle replicate) of one patient.

Usage:
    python run_celloracle_task.py --h5ad P_clones.h5ad --results-dir DIR --kind {original|bootstrap|shuffle} [--rep i]
"""

import os
import argparse
import logging
import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s %(levelname)-7s %(message)s",
                    datefmt="%H:%M:%S")
log = logging.getLogger("co_task")


def unit_path(results_dir, kind, rep):
    sub = kind if kind == "original" else os.path.join(kind, "%s_%03d" % (kind, rep))
    name = "original" if kind == "original" else "%s_%03d" % (kind, rep)
    return os.path.join(results_dir, kind, name + ".npz")


def reconstruct_unit(adata, kind, rep):
    """Return (X_raw_counts ndarray, labels str-array) for the requested unit."""
    cats = np.asarray(adata.uns["clone_categories"]).astype(str)
    clone = adata.obs["clone"].astype(str).to_numpy()
    X = adata.X
    if kind == "original":
        return X, clone
    if kind == "bootstrap":
        idx = np.asarray(adata.uns["bootstrap"]["indices"][rep])
        return X[idx], clone[idx]
    if kind == "shuffle":
        codes = np.asarray(adata.uns["shuffle"]["labels"][rep])
        return X, cats[codes]
    raise ValueError("unknown kind %s" % kind)


def run_grn(X, labels, var_names, base_GRN):
    """Reference genGRN pipeline (default params), run on ALL clones. -> links_dict.

    var_names MUST be the real gene symbols: CellOracle matches the base GRN to
    gene names *inside* get_links, so the names have to be set before any step.
    """
    import scanpy as sc
    import scipy.sparse as sp
    import celloracle as co

    Xd = X.toarray() if sp.issparse(X) else np.asarray(X)
    df = pd.DataFrame(Xd, columns=[str(g) for g in var_names]).astype("float64")
    adata = sc.AnnData(df)
    adata.var_names = [str(g) for g in var_names]
    adata.var_names_make_unique()
    adata.obs["clone"] = pd.Categorical(labels.astype(str))

    sc.pp.filter_genes(adata, min_counts=1)
    sc.pp.normalize_per_cell(adata, key_n_counts="n_counts_all")
    fr = sc.pp.filter_genes_dispersion(adata.X, flavor="cell_ranger",
                                       n_top_genes=min(adata.n_vars, 2000), log=False)
    adata = adata[:, fr.gene_subset]
    sc.pp.normalize_per_cell(adata)
    adata.raw = adata
    adata.layers["raw_count"] = adata.raw.X.copy()
    sc.pp.log1p(adata)
    sc.pp.scale(adata)
    sc.tl.pca(adata, svd_solver="arpack")

    oracle = co.Oracle()
    adata.X = adata.layers["raw_count"].copy()
    oracle.import_anndata_as_raw_count(adata=adata, cluster_column_name="clone",
                                       embedding_name="X_pca")
    oracle.import_TF_data(TF_info_matrix=base_GRN)
    oracle.perform_PCA()
    n_cell = oracle.adata.shape[0]
    k = max(1, int(0.025 * n_cell))
    try:
        n_comps = np.where(np.diff(np.diff(
            np.cumsum(oracle.pca.explained_variance_ratio_)) > 0.002))[0][0]
        n_comps = int(min(max(n_comps, 2), 50))
    except Exception:
        n_comps = min(50, oracle.adata.obsm["X_pca"].shape[1])
        log.warning("n_comps auto-selection failed; using %d", n_comps)
    oracle.knn_imputation(n_pca_dims=n_comps, k=k, balanced=True,
                          b_sight=k * 8, b_maxl=k * 4, n_jobs=4)
    links = oracle.get_links(cluster_name_for_GRN_unit="clone", alpha=10,
                             verbose_level=0)
    return links.links_dict


def links_to_npz(links_dict, var_names, out_path):
    """Pack per-clone link tables into one compact edgelist .npz."""
    clones = list(links_dict.keys())
    frames = []
    for ci, cl in enumerate(clones):
        df = links_dict[cl].copy()
        df.insert(0, "_clone", ci)
        frames.append(df)
    alldf = pd.concat(frames, ignore_index=True)
    value_cols = [c for c in alldf.columns
                  if c not in ("_clone", "source", "target")]
    # gene vocabulary local to this unit
    genes = pd.unique(pd.concat([alldf["source"], alldf["target"]]).astype(str))
    gcode = {g: i for i, g in enumerate(genes)}
    edges = np.column_stack([
        alldf["_clone"].to_numpy(np.int32),
        alldf["source"].astype(str).map(gcode).to_numpy(np.int32),
        alldf["target"].astype(str).map(gcode).to_numpy(np.int32)])
    values = alldf[value_cols].to_numpy(np.float32)
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    tmp = out_path + ".tmp.npz"
    np.savez_compressed(tmp,
                        genes=np.asarray(genes, dtype="U"),
                        clones=np.asarray([str(c) for c in clones], dtype="U"),
                        edges=edges, values=values,
                        value_cols=np.asarray(value_cols, dtype="U"))
    os.replace(tmp, out_path)               # atomic: partial writes never look complete


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--h5ad", required=True)
    p.add_argument("--results-dir", required=True)
    p.add_argument("--kind", required=True, choices=["original", "bootstrap", "shuffle"])
    p.add_argument("--rep", type=int, default=-1)
    p.add_argument("--force", action="store_true")
    args = p.parse_args()

    out = unit_path(args.results_dir, args.kind, args.rep)
    if os.path.exists(out) and not args.force:
        log.info("[skip] already done: %s", out)
        return

    import anndata as ad
    adata = ad.read_h5ad(args.h5ad)
    X, labels = reconstruct_unit(adata, args.kind, args.rep)

    import celloracle as co
    base_GRN = co.data.load_human_promoter_base_GRN()

    var = np.asarray(adata.var_names).astype(str)
    log.info("[run] %s kind=%s rep=%d cells=%d clones=%s",
             os.path.basename(args.h5ad), args.kind, args.rep,
             X.shape[0], sorted(set(labels)))
    links_dict = run_grn(X, labels, var, base_GRN)
    links_to_npz(links_dict, var, out)
    log.info("[ok]  wrote %s", out)


if __name__ == "__main__":
    main()
