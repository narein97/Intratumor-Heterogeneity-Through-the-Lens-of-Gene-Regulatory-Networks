#!/usr/bin/env python3
"""Builds motif-pruned SCENIC regulons from the CellOracle GRNs (and optionally exports CellOracle edge lists as TSVs).

Usage:
    python run_scenic_from_celloracle.py --h5ad clones_out/Breast/Wu_et_al_2021/CID4495_clones.h5ad --out scenic_out --motif-file MOTIFS.tbl --db-glob '/path/to/hg38/*feather' [--aucell] [--export-celloracle]
"""
import os
import re
import glob
import argparse
import logging
import pickle as pkl
import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)-7s %(message)s", datefmt="%H:%M:%S")
log = logging.getLogger("scenic_from_co")

MIN_CELLS = 20   # clones with fewer cells are skipped (as in the original script)


def decode_edgelist(arr, columns, genes, clone_categories):
    df = pd.DataFrame(np.asarray(arr), columns=columns)
    df["clone"] = clone_categories[df["clone"].astype(int)]
    df["source"] = genes[df["source"].astype(int)]
    df["target"] = genes[df["target"].astype(int)]
    return df


def sample_names_from_path(h5ad_path):
    base = os.path.basename(h5ad_path)
    sample = re.sub(r"(_with_support)?\.h5ad$", "", base)
    sample = re.sub(r"_clones$", "", sample)
    dataset = os.path.basename(os.path.dirname(os.path.abspath(h5ad_path)))
    return dataset, sample


def expression_frame(adata_sub):
    """cells x genes DataFrame; log-normalise when X holds raw integer counts."""
    import scanpy as sc
    import scipy.sparse as sp
    a = adata_sub.copy()
    X = a.X
    vals = X.data if sp.issparse(X) else np.asarray(X)
    if vals.size and np.allclose(vals, np.round(vals)):
        sc.pp.normalize_total(a, target_sum=1e4)
        sc.pp.log1p(a)
    X = a.X.toarray() if sp.issparse(a.X) else np.asarray(a.X)
    return pd.DataFrame(X, index=a.obs_names, columns=a.var_names)


def export_celloracle(adata, out_dir, columns, genes, cats):
    os.makedirs(out_dir, exist_ok=True)
    whole = decode_edgelist(adata.uns["celloracle_original"], columns, genes, cats)
    boots = adata.uns.get("celloracle_bootstrap", {})
    boot_dfs = {int(k): decode_edgelist(v, columns, genes, cats) for k, v in boots.items()}
    keep = [c for c in ["source", "target", "coef_mean", "coef_abs", "p", "-logp"] if c in whole.columns]
    for clone, d in whole.groupby("clone"):
        d[keep].to_csv(os.path.join(out_dir, f"CL{clone}-whole-celloracle-net.tsv"), sep="\t", index=False)
    for k in sorted(boot_dfs):
        for clone, d in boot_dfs[k].groupby("clone"):
            d[keep].to_csv(os.path.join(out_dir, f"CL{clone}-{k + 1}-celloracle-net.tsv"), sep="\t", index=False)
    log.info("exported CellOracle edge lists (whole + %d bootstrap) to %s", len(boot_dfs), out_dir)


def run_patient(h5ad_path, out_root, motif_file, db_glob, p_thresh=None, aucell_flag=False,
                export=False, dataset=None, sample=None, num_workers=10):
    import anndata as ad
    ds, sm = sample_names_from_path(h5ad_path)
    dataset, sample = dataset or ds, sample or sm
    out_dir = os.path.join(out_root, dataset, sample)
    os.makedirs(out_dir, exist_ok=True)

    adata = ad.read_h5ad(h5ad_path)
    need = ["celloracle_original", "celloracle_columns", "celloracle_genes", "clone_categories"]
    miss = [k for k in need if k not in adata.uns]
    if miss:
        log.warning("[skip] %s: missing .uns key(s) %s (run the GRN-inference step first)", h5ad_path, miss)
        return
    columns = [str(c) for c in adata.uns["celloracle_columns"]]
    genes = np.asarray(adata.uns["celloracle_genes"]).astype(str)
    cats = np.asarray(adata.uns["clone_categories"]).astype(str)

    if export:
        export_celloracle(adata, out_dir, columns, genes, cats)

    from pyscenic.utils import modules_from_adjacencies
    from pyscenic.prune import prune2df, df2regulons
    from pyscenic.aucell import aucell
    from ctxcore.rnkdb import FeatherRankingDatabase as RankingDatabase
    from dask.diagnostics import ProgressBar

    edges = decode_edgelist(adata.uns["celloracle_original"], columns, genes, cats)
    dbs = [RankingDatabase(fname=f, name=os.path.splitext(os.path.basename(f))[0]) for f in glob.glob(db_glob)]
    if not dbs:
        raise SystemExit("no ranking databases matched --db-glob %r" % db_glob)

    clone_col = adata.obs["clone"].astype(str)
    for clone in sorted(clone_col.unique()):
        net_path = os.path.join(out_dir, f"CL{clone}-whole-scenic-net.tsv")
        if os.path.exists(net_path):
            log.info("CL%s: output exists, skipping", clone); continue
        sub = adata[(clone_col == clone).to_numpy()]
        if sub.n_obs < MIN_CELLS:
            log.info("CL%s: only %d cells, skipping", clone, sub.n_obs); continue

        adj = edges[edges["clone"] == clone].copy()
        if p_thresh is not None:
            adj = adj[adj["p"] < p_thresh]
        # same renaming as the original script: source->TF, -logp->importance
        adj = adj.rename(columns={"source": "TF", "-logp": "importance"})[["TF", "target", "importance"]]
        data = expression_frame(sub)
        adj = adj[adj["TF"].isin(data.columns) & adj["target"].isin(data.columns)]
        if adj.empty:
            log.warning("CL%s: no usable edges, skipping", clone); continue
        tf_names = sorted(set(adj["TF"]) | set(adj["target"]))
        ex_matrix = data.loc[:, tf_names]

        modules = list(modules_from_adjacencies(adj, ex_matrix))
        with ProgressBar():
            df = prune2df(dbs, modules, motif_file, client_or_address="custom_multiprocessing", num_workers=num_workers)
        if df is None or len(df) == 0:
            log.warning("CL%s: pruning returned no enriched motifs, skipping", clone); continue
        regulons = df2regulons(df)

        df.to_csv(os.path.join(out_dir, f"CL{clone}-whole-scenic-motifs.csv"), sep="\t")
        with open(os.path.join(out_dir, f"CL{clone}-whole-scenic-regulons.p"), "wb") as fh:
            pkl.dump(regulons, fh)
        rows = [{"Gene1": reg.name.split("(")[0], "Gene2": tgt, "weight": w}
                for reg in regulons for tgt, w in reg.gene2weight.items()]
        pd.DataFrame(rows).to_csv(net_path, sep="\t", index=False)     # written last: also the "done" marker
        if aucell_flag:
            aucell(ex_matrix, regulons, num_workers=4).to_csv(
                os.path.join(out_dir, f"CL{clone}-whole-scenic-AUCell_scores.tsv"), sep="\t")
        log.info("CL%s: %d regulons, %d TF-target pairs", clone, len(regulons), len(rows))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--h5ad", help="one patient h5ad containing .uns['celloracle_*'] (step 3/4 output)")
    g.add_argument("--root", help="directory searched recursively for *_clones.h5ad (ignores *_with_support.h5ad)")
    ap.add_argument("--out", required=True, help="output root; files go to <out>/<Dataset>/<Sample>/")
    ap.add_argument("--motif-file", help="cisTarget motif annotation table (.tbl)")
    ap.add_argument("--db-glob", help="glob for the cisTarget ranking databases, e.g. '/path/hg38/*feather'")
    ap.add_argument("--p-thresh", type=float, default=None, help="keep only CellOracle edges with p < this (default: all edges, as in the original)")
    ap.add_argument("--aucell", action="store_true", help="also write SCENIC AUCell scores")
    ap.add_argument("--export-celloracle", action="store_true", help="also write CL<clone>-{whole|i}-celloracle-net.tsv files for downstream analysis")
    ap.add_argument("--dataset"); ap.add_argument("--sample")
    ap.add_argument("--num-workers", type=int, default=10)
    a = ap.parse_args()
    if not (a.motif_file and a.db_glob):
        ap.error("--motif-file and --db-glob are required")
    files = [a.h5ad] if a.h5ad else sorted(p for p in glob.glob(os.path.join(a.root, "**", "*_clones.h5ad"), recursive=True))
    for f in files:
        log.info("Processing %s", f)
        run_patient(f, a.out, a.motif_file, a.db_glob, a.p_thresh, a.aucell, a.export_celloracle,
                    a.dataset, a.sample, a.num_workers)


if __name__ == "__main__":
    main()
