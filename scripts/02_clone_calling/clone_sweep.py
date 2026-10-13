#!/usr/bin/env python3
"""Calls tumour clones per patient from 3CA copy-number profiles (Leiden over-clustering, then permutation-test merging).

Usage:
    python clone_sweep.py --data-root data --gene-order hg38_gencode_v27.txt --materialize --out clones_out    # add --inspect first to check the data layout
"""

import os
import glob
import json
import argparse
import logging
import warnings

import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s %(levelname)-7s %(message)s",
                    datefmt="%H:%M:%S")
log = logging.getLogger("clone_sweep")

# Heavy imports (scanpy/matplotlib) are deferred into main() so --help is cheap.

# ---------------------------------------------------------------------------
# Configuration / heuristics
# ---------------------------------------------------------------------------

# Candidate column names in the 3CA Cells.csv (matched case-insensitively).
SAMPLE_COLS   = ["sample", "samples", "sample_name", "sampleid",
                 "patient", "patient_id", "donor", "orig.ident"]
CELLTYPE_COLS = ["cell_type", "celltype", "cell.type", "cell_type_label",
                 "annotation", "cell_type_original", "majority_voting"]
CELLNAME_COLS = ["cell_name", "cell", "cellname", "barcode", "name",
                 "cell_id", "unnamed: 0"]

DEFAULT_MALIGNANT = ["Malignant"]

# Approximate hg38 centromere midpoints (Mb) for p/q arm assignment.
# Stable, well-known values; replace with exact coordinates if you need them.
HG38_CENTROMERE_MB = {
    "chr1": 123.4, "chr2": 93.9, "chr3": 90.9, "chr4": 50.0, "chr5": 48.8,
    "chr6": 59.8, "chr7": 60.1, "chr8": 45.2, "chr9": 43.0, "chr10": 39.8,
    "chr11": 53.4, "chr12": 35.5, "chr13": 17.7, "chr14": 17.2, "chr15": 19.0,
    "chr16": 36.8, "chr17": 25.1, "chr18": 18.5, "chr19": 26.2, "chr20": 28.1,
    "chr21": 12.0, "chr22": 15.0, "chrX": 61.0, "chrY": 10.4,
}
CHROM_ORDER = ["chr%d" % i for i in range(1, 23)] + ["chrX", "chrY"]


# ---------------------------------------------------------------------------
# Dataset discovery & file location
# ---------------------------------------------------------------------------

def find_datasets(data_root, cancer_subset=None, dataset_subset=None):
    """Yield (cancer_type, dataset_name, dataset_dir) for every dataset folder."""
    out = []
    if not os.path.isdir(data_root):
        log.error("data root does not exist: %s", data_root)
        return out
    for cancer in sorted(os.listdir(data_root)):
        cdir = os.path.join(data_root, cancer)
        if not os.path.isdir(cdir):
            continue
        if cancer_subset and cancer not in cancer_subset:
            continue
        for ds in sorted(os.listdir(cdir)):
            ddir = os.path.join(cdir, ds)
            if not os.path.isdir(ddir):
                continue
            if dataset_subset and ds not in dataset_subset:
                continue
            out.append((cancer, ds, ddir))
    return out


def locate_files(ddir):
    """
    Inspect one dataset folder and return what we can use.
    reason is None when usable, else a human-readable skip reason.
    """
    info = {"exp_dir": None, "mtx": None, "genes": None, "cells": None,
            "cna": None, "nested_exp": False, "nested_cna": False, "reason": None}

    exp_dirs = [d for d in glob.glob(os.path.join(ddir, "Data_*"))
                if os.path.isdir(d)]
    if not exp_dirs:
        info["reason"] = "no extracted Data_* directory (only a tarball?)"
        return info
    exp_dir = exp_dirs[0]
    info["exp_dir"] = exp_dir

    subdirs = [d for d in glob.glob(os.path.join(exp_dir, "*")) if os.path.isdir(d)]
    mtxs = (glob.glob(os.path.join(exp_dir, "Exp_data_*.mtx")) +
            glob.glob(os.path.join(exp_dir, "exp_data_*.mtx")))
    if subdirs and not mtxs:
        info["nested_exp"] = True
        info["reason"] = ("nested expression (%s) -- multi-platform/batch, skip in v1"
                          % ", ".join(os.path.basename(s) for s in subdirs))
        return info
    if len(mtxs) > 1:
        info["nested_exp"] = True
        info["reason"] = "multiple Exp_data_*.mtx in one folder -- skip in v1"
        return info
    if not mtxs:
        info["reason"] = "no Exp_data_*.mtx found"
        return info
    info["mtx"] = mtxs[0]

    genes = glob.glob(os.path.join(exp_dir, "[Gg]enes*.txt"))
    cells = glob.glob(os.path.join(exp_dir, "[Cc]ells*.csv"))
    if not genes or not cells:
        info["reason"] = "missing Genes.txt or Cells.csv next to the matrix"
        return info
    info["genes"] = sorted(genes)[0]
    info["cells"] = sorted(cells)[0]

    cna_single = (glob.glob(os.path.join(ddir, "CNA_matrix_*.csv.gz")) +
                  glob.glob(os.path.join(ddir, "CNA_matrix_*.csv")))
    cna_dirs = [d for d in glob.glob(os.path.join(ddir, "CNA_matrix_*"))
                if os.path.isdir(d)]
    if cna_single:
        info["cna"] = sorted(cna_single)[0]
    elif cna_dirs:
        info["nested_cna"] = True
        info["reason"] = "CNA matrix split into groups/platforms -- skip in v1"
        return info
    else:
        info["reason"] = "no CNA matrix (PNG-only dataset)"
        return info

    return info


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------

def _pick_col(columns, candidates):
    low = {str(c).lower(): c for c in columns}
    for cand in candidates:
        if cand in low:
            return low[cand]
    return None


def load_expression(info, sc):
    """Build a cells x genes AnnData from the 3CA flat files."""
    adata = sc.read_mtx(info["mtx"])

    genes = pd.read_csv(info["genes"], header=None)[0].astype(str).tolist()
    cells_df = pd.read_csv(info["cells"])

    n_genes, n_cells = len(genes), len(cells_df)
    if adata.shape == (n_genes, n_cells):
        adata = adata.T
    elif adata.shape == (n_cells, n_genes):
        pass
    else:
        raise ValueError("matrix %s matches neither cells(%d)xgenes(%d)"
                         % (adata.shape, n_cells, n_genes))

    idx_col = _pick_col(cells_df.columns, CELLNAME_COLS) or cells_df.columns[0]
    cells_df = cells_df.set_index(idx_col)
    adata.obs = cells_df.copy()
    adata.obs_names = cells_df.index.astype(str)
    adata.var_names = pd.Index(genes).astype(str)
    adata.var_names_make_unique()

    sample_col   = _pick_col(adata.obs.columns, SAMPLE_COLS)
    celltype_col = _pick_col(adata.obs.columns, CELLTYPE_COLS)
    return adata, sample_col, celltype_col


def load_cna(cna_path, obs_names):
    """Load CNA matrix oriented to cells (rows) x genes (cols)."""
    df = pd.read_csv(cna_path, index_col=0)
    df.index = df.index.astype(str)
    df.columns = df.columns.astype(str)
    cells = set(map(str, obs_names))
    as_is = len(cells.intersection(df.index))
    transp = len(cells.intersection(df.columns))
    if transp > as_is:
        df = df.T
    return df


# ---------------------------------------------------------------------------
# Genome ordering & arm-level aggregation
# ---------------------------------------------------------------------------

def load_gene_order(path):
    """Parse the gencode-style gene order table (gene, chrom, start, ...)."""
    go = pd.read_csv(path, names=range(5), sep="\t")
    go = go.rename(columns={0: "gene", 1: "chrom", 2: "start"})
    go["gene"] = go["gene"].astype(str)
    go["chrom"] = go["chrom"].astype(str)
    go["start"] = pd.to_numeric(go["start"], errors="coerce")
    go = go.dropna(subset=["start"]).drop_duplicates("gene").set_index("gene")
    return go[["chrom", "start"]]


def genome_ordered_genes(var_names, gene_order):
    """var_names ordered by chromosome then start position (for the heatmap)."""
    common = [g for g in var_names if g in gene_order.index]
    sub = gene_order.loc[common].copy()
    sub["chrom"] = pd.Categorical(sub["chrom"], categories=CHROM_ORDER, ordered=True)
    sub = sub.dropna(subset=["chrom"]).sort_values(["chrom", "start"])
    return sub.index.tolist()


def gene_to_bin(var_names, gene_order, level="arm", bin_mb=10.0,
                centromeres=HG38_CENTROMERE_MB):
    """Map each gene to a genomic bin. Returns (gene->bin dict, ordered bins)."""
    common = [g for g in var_names if g in gene_order.index]
    sub = gene_order.loc[common]
    mapping = {}
    for gene, row in sub.iterrows():
        chrom, start = row["chrom"], row["start"]
        if chrom not in CHROM_ORDER:
            continue
        if level == "chromosome":
            mapping[gene] = chrom
        elif level == "bin":
            mapping[gene] = "%s:%d" % (chrom, int(start // (bin_mb * 1e6)))
        else:  # arm
            cen = centromeres.get(chrom, np.nan) * 1e6
            arm = "p" if (not np.isnan(cen) and start < cen) else "q"
            mapping[gene] = chrom + arm

    def sort_key(b):
        chrom = b.split(":")[0].rstrip("pq")
        ci = CHROM_ORDER.index(chrom) if chrom in CHROM_ORDER else 99
        sub_key = 0
        if level == "arm":
            sub_key = 0 if b.endswith("p") else 1
        elif level == "bin":
            sub_key = int(b.split(":")[1])
        return (ci, sub_key)

    bins = sorted(set(mapping.values()), key=sort_key)
    return mapping, bins


def arm_profiles(cna_cells_by_gene, gene_bin_map, bins):
    """Collapse cells x genes CNA to cells x bins (mean per bin)."""
    keep = [g for g in cna_cells_by_gene.columns if g in gene_bin_map]
    sub = cna_cells_by_gene[keep]
    binlabels = pd.Index([gene_bin_map[g] for g in keep], name="bin")
    prof = sub.T.groupby(binlabels).mean().T
    prof = prof.reindex(columns=[b for b in bins if b in prof.columns])
    return prof


# ---------------------------------------------------------------------------
# Scores
# ---------------------------------------------------------------------------

def distinctness(arm_prof, labels):
    """
    Pairwise distance (1 - pearson r) between clone arm-level centroids.
    Returns (min_distance, mean_distance); min = closest/worst-separated pair.

    Robust to flat centroids: if a clone's profile has zero variance,
    correlation is undefined, but a flat clone IS distinct from a CNA-bearing
    one, so that pair is assigned the maximal distance (1.0) rather than NaN.
    """
    df = pd.DataFrame(arm_prof.values, index=np.asarray(labels),
                      columns=arm_prof.columns)
    cent = df.groupby(level=0).mean()
    if cent.shape[0] < 2:
        return np.nan, np.nan

    M = np.nan_to_num(cent.values, nan=0.0)
    mean = M.mean(axis=1, keepdims=True)
    std = M.std(axis=1, keepdims=True)
    with np.errstate(invalid="ignore", divide="ignore"):
        Z = (M - mean) / std
    n, b = M.shape
    pairs = []
    for i in range(n):
        for j in range(i + 1, n):
            if std[i, 0] == 0 or std[j, 0] == 0:
                pairs.append(1.0)            # flat vs non-flat -> distinct
            else:
                r = float(np.clip(np.dot(Z[i], Z[j]) / b, -1.0, 1.0))
                pairs.append(1.0 - r)
    pairs = np.asarray(pairs)
    return float(np.min(pairs)), float(np.mean(pairs))


def stability_ari(X_pca, full_labels, resolution, sc, anndata,
                  reps=5, frac=0.8, n_neighbors=20, seed=0):
    """Mean ARI between the full partition and partitions on subsampled cells."""
    from sklearn.metrics import adjusted_rand_score
    rng = np.random.default_rng(seed)
    n = X_pca.shape[0]
    full = np.asarray(full_labels)
    aris = []
    for r in range(reps):
        idx = rng.choice(n, size=max(2, int(frac * n)), replace=False)
        sub = anndata.AnnData(X=np.zeros((len(idx), 1), dtype=np.float32))
        sub.obsm["X_pca"] = X_pca[idx]
        try:
            sc.pp.neighbors(sub, use_rep="X_pca",
                            n_neighbors=min(n_neighbors, len(idx) - 1))
            sc.tl.leiden(sub, resolution=resolution, key_added="c",
                         random_state=seed + r)
        except Exception:
            continue
        aris.append(adjusted_rand_score(full[idx], sub.obs["c"].values))
    return float(np.mean(aris)) if aris else np.nan


# ---------------------------------------------------------------------------
# Plotting (montages across resolutions)
# ---------------------------------------------------------------------------

def clone_colors(labels, plt):
    uniq = sorted(pd.unique(labels),
                  key=lambda x: int(x) if str(x).isdigit() else str(x))
    pal = plt.get_cmap("tab10")
    return {c: pal(i % 10) for i, c in enumerate(uniq)}


def heatmap_panel(ax, cna_ordered, labels, plt, vlim=5.0):
    """Genome-ordered CNA heatmap (cells sorted by clone) with a clone sidebar."""
    order = np.argsort([int(l) if str(l).isdigit() else str(l) for l in labels])
    mat = cna_ordered[order]
    lab = np.asarray(labels)[order]
    ax.imshow(mat, aspect="auto", cmap="bwr", vmin=-vlim, vmax=vlim,
              interpolation="nearest")
    ax.set_xticks([]); ax.set_yticks([])
    cmap = clone_colors(lab, plt)
    side = ax.inset_axes([-0.045, 0, 0.03, 1])
    strip = np.array([[cmap[l][:3] for l in lab]]).transpose(1, 0, 2)
    side.imshow(strip, aspect="auto"); side.set_xticks([]); side.set_yticks([])


def make_montages(sample_dir, panels, plt):
    """Build a UMAP grid and a genome-ordered-heatmap grid across resolutions."""
    n = len(panels)
    if n == 0:
        return
    ncol = min(5, n)
    nrow = int(np.ceil(n / ncol))

    fig, axes = plt.subplots(nrow, ncol, figsize=(3.2 * ncol, 3.0 * nrow),
                             squeeze=False)
    for k, p in enumerate(panels):
        ax = axes[k // ncol][k % ncol]
        cmap = clone_colors(p["labels"], plt)
        ax.scatter(p["umap_xy"][:, 0], p["umap_xy"][:, 1], s=3,
                   c=[cmap[l] for l in p["labels"]], linewidths=0)
        ax.set_title("res=%.3g  %s" % (p["res"], p["scores"]), fontsize=8)
        ax.set_xticks([]); ax.set_yticks([])
    for k in range(n, nrow * ncol):
        axes[k // ncol][k % ncol].axis("off")
    fig.tight_layout()
    fig.savefig(os.path.join(sample_dir, "umap_sweep.png"), dpi=130)
    plt.close(fig)

    fig, axes = plt.subplots(nrow, ncol, figsize=(3.6 * ncol, 3.0 * nrow),
                             squeeze=False)
    for k, p in enumerate(panels):
        ax = axes[k // ncol][k % ncol]
        heatmap_panel(ax, p["cna_ordered"], p["labels"], plt)
        ax.set_title("res=%.3g  %s" % (p["res"], p["scores"]), fontsize=8)
    for k in range(n, nrow * ncol):
        axes[k // ncol][k % ncol].axis("off")
    fig.tight_layout()
    fig.savefig(os.path.join(sample_dir, "cna_heatmap_sweep.png"), dpi=130)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Per-sample sweep
# ---------------------------------------------------------------------------

def sweep_sample(adata_s, cna_s, gene_order, args, sc, anndata, plt, sample_dir):
    """Resolution sweep for one sample. Returns a list of per-resolution rows."""
    adata_s = adata_s.copy()
    cna_raw = cna_s.loc[adata_s.obs_names, adata_s.var_names]   # raw CNA, cells x genes
    adata_s.X = cna_raw.to_numpy()
    sc.pp.scale(adata_s, max_value=10)
    n_comps = max(2, min(args.n_pcs, adata_s.n_obs - 1, adata_s.n_vars - 1))
    sc.tl.pca(adata_s, n_comps=n_comps)
    sc.pp.neighbors(adata_s, use_rep="X_pca",
                    n_neighbors=min(args.n_neighbors, adata_s.n_obs - 1))
    sc.tl.umap(adata_s, random_state=args.seed)
    X_pca = adata_s.obsm["X_pca"]
    umap_xy = adata_s.obsm["X_umap"]

    gbin, bins = gene_to_bin(adata_s.var_names, gene_order,
                             level=args.aggregate, bin_mb=args.bin_mb)
    # Distinctness on RAW CNA (pre-scale): per-gene z-scoring forces two-clone
    # centroids to correlation -1 (distinct == 2.0 for every k=2 split), so the
    # unscaled copy is used here to make the score informative at k=2.
    arm = arm_profiles(cna_raw, gbin, bins)
    # Heatmap uses RAW CNA (pre-scale) so shared and clone-specific copy-number
    # events show as actual gains/losses; white is centred on the baseline.
    ordered_genes = genome_ordered_genes(adata_s.var_names, gene_order)
    cna_ordered_full = cna_raw[ordered_genes].to_numpy()

    rows, panels = [], []
    for res in args.resolutions:
        sc.tl.leiden(adata_s, resolution=res, key_added="clone",
                     random_state=args.seed)
        labels = adata_s.obs["clone"].astype(str).values
        k = len(np.unique(labels))
        sizes = pd.Series(labels).value_counts().sort_index().to_dict()

        sil = d_min = d_mean = stab = np.nan
        if k >= 2:
            from sklearn.metrics import silhouette_score
            try:
                sil = float(silhouette_score(X_pca, labels))
            except Exception:
                pass
            d_min, d_mean = distinctness(arm, labels)
            if not args.skip_stability:
                stab = stability_ari(X_pca, labels, res, sc, anndata,
                                     reps=args.stability_reps,
                                     frac=args.stability_frac,
                                     n_neighbors=args.n_neighbors,
                                     seed=args.seed)

        rows.append({"resolution": res, "n_clones": k,
                     "min_clone_size": min(sizes.values()) if sizes else 0,
                     "clone_sizes": json.dumps(sizes),
                     "silhouette": sil, "distinct_min": d_min,
                     "distinct_mean": d_mean, "stability_ari": stab})

        if not args.no_plots:
            def _f(v):
                return ("%.2f" % v) if v == v else "na"
            panels.append({"res": res, "umap_xy": umap_xy, "labels": labels,
                           "cna_ordered": cna_ordered_full,
                           "scores": "k=%d sil=%s dmin=%s st=%s" % (
                               k, _f(sil), _f(d_min), _f(stab))})

    if not args.no_plots:
        try:
            make_montages(sample_dir, panels, plt)
        except Exception as e:
            log.warning("    montage failed: %s", e)
    return rows


# ---------------------------------------------------------------------------
# Prune mode: over-cluster once, then significance-merge (CHOIR-style)
# ---------------------------------------------------------------------------

def make_gaussian_test(n_null, max_cells, seed0):
    """
    Return test(Xa, Xb) -> (observed_separation, p_value) implementing the sc-SHC
    principle: the null hypothesis is that A and B are ONE Gaussian blob, not two
    clones. We compare the observed A-vs-B separation (Mahalanobis distance between
    centroids, in the pooled-covariance metric) to the separation obtained by
    splitting a single fitted Gaussian along its top principal component.

    This fixes over-splitting: over-clustering a homogeneous clone produces sub-
    groups that a permutation/classifier test would call "significant" (the split
    was chosen to be separable), but whose separation does NOT exceed what
    splitting one Gaussian yields -- so here they merge. Large p => merge.
    Operates in the CNA-PCA space where the clusters were actually formed.
    """
    def test(Xa, Xb, seed=0):
        rng = np.random.default_rng(seed0 + seed)
        Xa = np.nan_to_num(np.asarray(Xa, float))
        Xb = np.nan_to_num(np.asarray(Xb, float))

        def sub(X):
            if len(X) > max_cells:
                return X[rng.choice(len(X), max_cells, replace=False)]
            return X
        Xa, Xb = sub(Xa), sub(Xb)
        if len(Xa) < 3 or len(Xb) < 3:          # too few cells -> merge
            return 0.0, 1.0

        X = np.vstack([Xa, Xb])
        n, d = X.shape
        cov = np.cov(X, rowvar=False) + np.eye(d) * 1e-6
        Cinv = np.linalg.pinv(cov)

        def sep(A, B):                           # Mahalanobis distance of centroids
            diff = A.mean(0) - B.mean(0)
            return float(np.sqrt(max(diff @ Cinv @ diff, 0.0)))

        obs = sep(Xa, Xb)
        mu = X.mean(0)
        null = np.empty(n_null)
        for i in range(n_null):
            sim = rng.multivariate_normal(mu, cov, size=n)
            simc = sim - sim.mean(0)
            # split the simulated blob along its top principal component (best axis)
            try:
                _, _, vt = np.linalg.svd(simc, full_matrices=False)
                proj = simc @ vt[0]
            except Exception:
                proj = simc[:, 0]
            gA, gB = sim[proj <= 0], sim[proj > 0]
            null[i] = sep(gA, gB) if len(gA) > 1 and len(gB) > 1 else 0.0
        pval = (1 + np.sum(null >= obs)) / (1 + n_null)
        return obs, float(pval)

    return test


def prune_clusters(X, labels, test, alpha, seed, log_fn=None):
    """
    Agglomerative significance merge in feature space X (cells x dims, aligned to
    labels). Repeatedly test the closest remaining pair of clusters (by centroid
    distance); merge if their split is NOT more separated than a single Gaussian
    blob (Bonferroni-corrected over current pairs), else stop -- if even the
    closest pair is a real split, all pairs are. Collapses to one clone when no
    split beats the Gaussian null. Returns relabelled 0..k-1 assignments.
    """
    labels = np.asarray(labels).astype(object)
    X = np.nan_to_num(np.asarray(X, float))

    while True:
        uniq = list(pd.unique(labels))
        if len(uniq) < 2:
            break
        cents = {c: X[labels == c].mean(0) for c in uniq}
        pairs = [(uniq[i], uniq[j])
                 for i in range(len(uniq)) for j in range(i + 1, len(uniq))]
        pairs.sort(key=lambda p: np.linalg.norm(cents[p[0]] - cents[p[1]]))
        a, b = pairs[0]                          # closest (least separated) pair
        obs, p = test(X[labels == a], X[labels == b], seed=seed)
        p_adj = min(1.0, p * len(pairs))         # Bonferroni over current pairs

        if p_adj > alpha:
            labels[labels == b] = a              # not beyond noise -> merge
            if log_fn:
                log_fn("merge %s+%s (sep=%.2f p_adj=%.3f)" % (a, b, obs, p_adj))
        else:
            if log_fn:
                log_fn("closest pair %s,%s distinct (sep=%.2f p_adj=%.3f) -> stop"
                       % (a, b, obs, p_adj))
            break

    remap = {c: str(i) for i, c in enumerate(pd.unique(labels))}
    return np.array([remap[c] for c in labels])


def _prune_core(adata_s, cna_s, gene_order, args, sc, anndata):
    """Embed in CNA space, over-cluster once, significance-merge. Returns a dict
    with labels (aligned to adata_s.obs_names order) and plotting aux."""
    adata_s = adata_s.copy()
    cna_raw = cna_s.loc[adata_s.obs_names, adata_s.var_names]
    adata_s.X = cna_raw.to_numpy()
    sc.pp.scale(adata_s, max_value=10)
    n_comps = max(2, min(args.n_pcs, adata_s.n_obs - 1, adata_s.n_vars - 1))
    sc.tl.pca(adata_s, n_comps=n_comps)
    sc.pp.neighbors(adata_s, use_rep="X_pca",
                    n_neighbors=min(args.n_neighbors, adata_s.n_obs - 1))
    sc.tl.umap(adata_s, random_state=args.seed)

    gbin, bins = gene_to_bin(adata_s.var_names, gene_order,
                             level=args.aggregate, bin_mb=args.bin_mb)
    arm = arm_profiles(cna_raw, gbin, bins).fillna(0.0)
    ordered_genes = genome_ordered_genes(adata_s.var_names, gene_order)
    cna_ordered_full = cna_raw[ordered_genes].to_numpy()

    sc.tl.leiden(adata_s, resolution=args.prune_resolution,
                 key_added="clone_init", random_state=args.seed)
    init_labels = adata_s.obs["clone_init"].astype(str).values
    n_init = len(np.unique(init_labels))

    test = make_gaussian_test(args.prune_perm, args.prune_max_cells, args.seed)
    logs = []
    final_labels = prune_clusters(adata_s.obsm["X_pca"], init_labels, test,
                                  args.prune_alpha, args.seed, log_fn=logs.append)
    for m in logs:
        log.info("       %s", m)
    log.info("       pruned %d -> %d clone(s)", n_init, len(np.unique(final_labels)))

    return {"labels": final_labels, "n_init": n_init,
            "obs_names": list(adata_s.obs_names),
            "X_pca": adata_s.obsm["X_pca"], "umap_xy": adata_s.obsm["X_umap"],
            "arm": arm, "cna_ordered_full": cna_ordered_full}


def prune_sample(adata_s, cna_s, gene_order, args, sc, anndata, plt, sample_dir):
    """Over-cluster at a high resolution, then significance-merge to a clone count."""
    core = _prune_core(adata_s, cna_s, gene_order, args, sc, anndata)
    final_labels = core["labels"]
    k = len(np.unique(final_labels))
    sizes = pd.Series(final_labels).value_counts().sort_index().to_dict()

    sil = d_min = d_mean = np.nan
    if k >= 2:
        from sklearn.metrics import silhouette_score
        try:
            sil = float(silhouette_score(core["X_pca"], final_labels))
        except Exception:
            pass
        d_min, d_mean = distinctness(core["arm"], final_labels)

    row = {"resolution": args.prune_resolution, "n_init_clusters": core["n_init"],
           "n_clones": k, "min_clone_size": min(sizes.values()) if sizes else 0,
           "clone_sizes": json.dumps(sizes), "silhouette": sil,
           "distinct_min": d_min, "distinct_mean": d_mean, "stability_ari": np.nan}

    if not args.no_plots:
        def _f(v):
            return ("%.2f" % v) if v == v else "na"
        panel = {"res": args.prune_resolution, "umap_xy": core["umap_xy"],
                 "labels": final_labels, "cna_ordered": core["cna_ordered_full"],
                 "scores": "k=%d (from %d) sil=%s dmin=%s" % (
                     k, core["n_init"], _f(sil), _f(d_min))}
        try:
            make_montages(sample_dir, [panel], plt)
        except Exception as e:
            log.warning("    montage failed: %s", e)

    return [row]


# ---------------------------------------------------------------------------
# Materialize mode: write one clone-labelled .h5ad per dataset (uses prune)
# ---------------------------------------------------------------------------

def _sanitize_for_h5ad(adata):
    """
    h5py cannot serialise object columns that mix strings with NaN/None (common in
    3CA Cells.csv metadata). Rebuild obs & var as fresh DataFrames where every
    object/categorical column is converted to plain strings ("NA" for missing).
    Numeric columns pass through unchanged. Rebuilding (rather than overwriting
    columns in place) avoids Copy-on-Write / dtype-inference quirks across pandas
    versions.
    """
    import pandas as _pd
    import numpy as _np
    from pandas.api import types as _pdt

    def clean_df(df):
        cols = {}
        for col in df.columns:
            s = df[col]
            if _pdt.is_numeric_dtype(s) or _pdt.is_bool_dtype(s):
                cols[col] = s.to_numpy()            # numeric/bool pass through
            else:                                    # object / category / StringDtype
                isna = _pd.isna(s).to_numpy()
                vals = s.astype(object).to_numpy()
                cols[col] = _np.array(
                    ["NA" if m else str(v) for v, m in zip(vals, isna)],
                    dtype=object)
        return _pd.DataFrame(cols, index=df.index.astype(str))

    adata.obs = clean_df(adata.obs)
    adata.var = clean_df(adata.var)


def materialize_mode(args, sc, anndata, plt, gene_order):
    """
    For each dataset, prune each patient's malignant cells into clones (prune
    setting), then write ONE .h5ad per sample containing only the cells that
    belong to clones with >= --min-clone-cells cells. Each output .h5ad has:
        .X               raw RNA counts (all genes), for the kept-clone cells
        .obs['clone']    integer clone labels 0,1,2,... (relabelled contiguously)
        .obsm['cnp']     raw CNA profile for those cells (cells x CNA-genes)
        .uns['cnp_genes'] the CNA gene names (columns of .obsm['cnp'])
    Cells in clones below the threshold (and non-malignant cells) are removed.
    Also writes a global manifest of kept clones for building job arrays.
    """
    manifest = []
    for cancer, ds, ddir in find_datasets(args.data_root, args.cancer, args.dataset):
        info = locate_files(ddir)
        if info["reason"] is not None:
            log.info("[skip] %s/%s -- %s", cancer, ds, info["reason"])
            continue
        try:
            adata_full, scol, ccol = load_expression(info, sc)   # raw RNA, all cells/genes
        except Exception as e:
            log.warning("[skip] %s/%s -- load error: %s", cancer, ds, e)
            continue
        if ccol is None or scol is None:
            log.warning("[skip] %s/%s -- no cell_type/sample column", cancer, ds)
            continue
        if not adata_full.obs[ccol].isin(args.malignant_label).any():
            log.warning("[skip] %s/%s -- no malignant cells", cancer, ds)
            continue

        cna = load_cna(info["cna"], adata_full.obs_names)
        shared_cells = list(set(adata_full.obs_names) & set(cna.index))
        shared_genes = list(set(adata_full.var_names) & set(cna.columns))
        if len(shared_cells) < args.min_malignant or len(shared_genes) < 100:
            log.info("[skip] %s/%s -- tiny overlap (cells=%d genes=%d)",
                     cancer, ds, len(shared_cells), len(shared_genes))
            continue
        cna = cna.loc[shared_cells, shared_genes]
        shared_set = set(shared_cells)
        log.info("[run]  %s/%s  cells=%d genes=%d", cancer, ds,
                 adata_full.n_obs, adata_full.n_vars)

        for sample in sorted(adata_full.obs[scol].astype(str).unique()):
            mask = ((adata_full.obs[scol].astype(str) == sample) &
                    (adata_full.obs[ccol].isin(args.malignant_label)) &
                    (adata_full.obs_names.isin(shared_set)))
            n_mal = int(mask.sum())
            if n_mal < args.min_malignant:
                log.info("    sample %s  malignant=%d < %d -> skip",
                         sample, n_mal, args.min_malignant)
                continue
            log.info("    sample %s  malignant=%d", sample, n_mal)

            # cluster in CNA space (shared genes); labels align to these cells
            adata_clust = adata_full[mask.values, shared_genes].copy()
            cna_s = cna.loc[adata_clust.obs_names]
            core = _prune_core(adata_clust, cna_s, gene_order, args, sc, anndata)
            labels = pd.Series(core["labels"], index=core["obs_names"])

            # drop MALIGNANT clones below the cell threshold (their cells removed)
            sizes = labels.value_counts()
            kept = [c for c, n in sizes.items() if n >= args.min_clone_cells]
            if len(kept) < 2:
                log.info("       %d malignant clone(s) >= %d cells -> skip "
                         "(need >= 2 clones)", len(kept), args.min_clone_cells)
                continue

            # non-malignant cells of THIS sample -> kept as clone 'N'
            nonmal = adata_full.obs_names[
                (adata_full.obs[scol].astype(str) == sample) &
                (~adata_full.obs[ccol].isin(args.malignant_label))]
            if len(nonmal) < args.min_normal_cells:
                log.info("       %d non-malignant cells < %d -> skip",
                         len(nonmal), args.min_normal_cells)
                continue

            keep_mal = labels[labels.isin(kept)].index

            # relabel kept clones contiguously to 0,1,2,... (largest -> 0)
            order = sizes.loc[kept].sort_values(ascending=False).index.tolist()
            remap = {old: str(i) for i, old in enumerate(order)}
            new_labels = labels.loc[keep_mal].map(remap)     # index = kept malignant cells

            # output = kept malignant clones (0,1,2,...) + all non-malignant ('N')
            out_cells = list(keep_mal) + list(nonmal)
            adata_out = adata_full[out_cells].copy()           # raw RNA, all genes

            clone_series = pd.Series("N", index=adata_out.obs_names, dtype=object)
            clone_series.loc[new_labels.index] = new_labels.values
            adata_out.obs["clone"] = pd.Categorical(clone_series.values)

            # raw CNA in obsm: real values for malignant cells, NaN for non-malignant
            cnp_df = pd.DataFrame(np.nan, index=adata_out.obs_names,
                                  columns=shared_genes, dtype="float32")
            cnp_df.loc[new_labels.index, :] = cna.loc[new_labels.index, shared_genes].values
            adata_out.obsm["cnp"] = cnp_df.to_numpy()
            adata_out.uns["cnp_genes"] = list(shared_genes)

            outdir = os.path.join(args.out, cancer, ds)
            os.makedirs(outdir, exist_ok=True)
            h5 = os.path.join(outdir, "%s_clones.h5ad" % sample)
            _sanitize_for_h5ad(adata_out)
            adata_out.write_h5ad(h5)
            counts = adata_out.obs["clone"].value_counts().sort_index().to_dict()
            log.info("       -> %s  (%d clones, %d malignant + %d normal='N': %s)",
                     h5, len(kept), len(keep_mal), len(nonmal), counts)

            # optional heatmap/UMAP montage of the malignant clones only (+/-5 scale)
            if not args.no_plots:
                try:
                    montage_dir = os.path.join(outdir, sample)
                    os.makedirs(montage_dir, exist_ok=True)
                    ordered_genes = genome_ordered_genes(shared_genes, gene_order)
                    cna_ord = cna.loc[new_labels.index, ordered_genes].to_numpy()
                    pos = {c: i for i, c in enumerate(core["obs_names"])}
                    umap_idx = [pos[c] for c in new_labels.index]
                    panel = {"res": args.prune_resolution,
                             "umap_xy": core["umap_xy"][umap_idx],
                             "labels": new_labels.values, "cna_ordered": cna_ord,
                             "scores": "k=%d nmal=%d" % (len(kept), len(keep_mal))}
                    make_montages(montage_dir, [panel], plt)
                except Exception as e:
                    log.warning("       montage failed: %s", e)

            for clone, n in counts.items():
                if clone == "N":
                    continue
                manifest.append({"cancer": cancer, "dataset": ds, "sample": sample,
                                 "clone": int(clone), "n_cells": int(n),
                                 "h5ad": os.path.abspath(h5)})

    os.makedirs(args.out, exist_ok=True)
    man = os.path.join(args.out, "clones_manifest.csv")
    pd.DataFrame(manifest).to_csv(man, index=False)
    log.info("\nWrote %d sample file(s), %d clones total. Manifest: %s",
             len({(m['cancer'], m['dataset'], m['sample']) for m in manifest}),
             len(manifest), man)


# ---------------------------------------------------------------------------
# Inspect mode
# ---------------------------------------------------------------------------

def inspect_mode(args, sc):
    log.info("INSPECT mode -- no clustering will run.\n")
    rows = []
    for cancer, ds, ddir in find_datasets(args.data_root, args.cancer, args.dataset):
        info = locate_files(ddir)
        rec = {"cancer": cancer, "dataset": ds,
               "usable": info["reason"] is None, "reason": info["reason"] or ""}
        if info["reason"] is None:
            try:
                adata, scol, ccol = load_expression(info, sc)
                cna = load_cna(info["cna"], adata.obs_names)
                shared_cells = len(set(adata.obs_names) & set(cna.index))
                shared_genes = len(set(adata.var_names) & set(cna.columns))
                mal = int(adata.obs[ccol].isin(args.malignant_label).sum()) if ccol else 0
                ctypes = (adata.obs[ccol].value_counts().head(8).to_dict()
                          if ccol else {})
                rec.update({"sample_col": scol, "celltype_col": ccol,
                            "n_cells": adata.n_obs, "n_genes": adata.n_vars,
                            "shared_cells": shared_cells,
                            "shared_genes": shared_genes, "malignant_cells": mal})
                log.info("[OK]   %-14s %-26s sample=%s celltype=%s cells=%d genes=%d "
                         "shared(cells=%d genes=%d) malignant=%d",
                         cancer, ds, scol, ccol, adata.n_obs, adata.n_vars,
                         shared_cells, shared_genes, mal)
                log.info("       cell types: %s", ctypes)
                if scol is None or ccol is None:
                    log.warning("       !! sample/celltype column NOT detected -- "
                                "inspect Cells.csv and extend *_COLS lists")
            except Exception as e:
                rec["reason"], rec["usable"] = "load error: %s" % e, False
                log.warning("[ERR]  %-14s %-26s %s", cancer, ds, e)
        else:
            log.info("[skip] %-14s %-26s %s", cancer, ds, info["reason"])
        rows.append(rec)
    os.makedirs(args.out, exist_ok=True)
    out = os.path.join(args.out, "inspect_report.csv")
    pd.DataFrame(rows).to_csv(out, index=False)
    log.info("\nWrote %s", out)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def build_argparser():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--data-root", default="data")
    p.add_argument("--out", default="clone_sweep_out")
    p.add_argument("--gene-order", default=None,
                   help="gencode-style table (gene, chrom, start). Required for the sweep.")
    p.add_argument("--cancer", nargs="+", default=None,
                   help="restrict to these cancer-type folders")
    p.add_argument("--dataset", nargs="+", default=None,
                   help="restrict to these dataset folder names (for parallel jobs)")
    p.add_argument("--malignant-label", nargs="+", default=DEFAULT_MALIGNANT)
    p.add_argument("--resolutions", nargs="+", type=float,
                   default=[0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0])
    p.add_argument("--min-malignant", type=int, default=100,
                   help="min malignant cells in a sample to cluster it")
    p.add_argument("--n-pcs", type=int, default=20)
    p.add_argument("--n-neighbors", type=int, default=20)
    p.add_argument("--aggregate", choices=["arm", "chromosome", "bin"],
                   default="arm")
    p.add_argument("--bin-mb", type=float, default=10.0)
    p.add_argument("--stability-reps", type=int, default=5)
    p.add_argument("--stability-frac", type=float, default=0.8)
    p.add_argument("--skip-stability", action="store_true")
    p.add_argument("--no-plots", action="store_true")
    p.add_argument("--seed", type=int, default=0)
    # clone-number selection
    p.add_argument("--select", choices=["sweep", "prune"], default="sweep",
                   help="sweep = per-resolution scores (default); "
                        "prune = over-cluster once, then significance-merge to a clone count")
    p.add_argument("--prune-resolution", type=float, default=1.0,
                   help="[prune] Leiden resolution to over-cluster at before merging")
    p.add_argument("--prune-alpha", type=float, default=0.05,
                   help="[prune] significance for keeping a split (Bonferroni). "
                        "RAISE it for MORE clones (less merging), lower for fewer.")
    p.add_argument("--prune-perm", type=int, default=200,
                   help="[prune] Gaussian-null simulations per merge test")
    p.add_argument("--prune-max-cells", type=int, default=400,
                   help="[prune] max cells per cluster used in the test (for speed)")
    # materialize (write clone-labelled h5ad per dataset for CellOracle)
    p.add_argument("--materialize", action="store_true",
                   help="prune each patient into clones and write one "
                        "<sample>_clones.h5ad per patient (raw RNA in .X, raw CNA "
                        "in .obsm['cnp'], integer obs['clone']); cells in clones "
                        "below --min-clone-cells are dropped")
    p.add_argument("--min-clone-cells", type=int, default=50,
                   help="[materialize] keep only clones with >= this many cells")
    p.add_argument("--min-normal-cells", type=int, default=50,
                   help="[materialize] skip a patient with fewer than this many "
                        "non-malignant (N) cells")
    p.add_argument("--inspect", action="store_true",
                   help="report what the loader finds per dataset; no clustering")
    return p


def main():
    args = build_argparser().parse_args()
    warnings.filterwarnings("ignore")

    import scanpy as sc
    import anndata
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    sc.settings.verbosity = 0
    np.random.seed(args.seed)
    os.makedirs(args.out, exist_ok=True)

    if args.inspect:
        inspect_mode(args, sc)
        return

    if not args.gene_order:
        log.error("--gene-order is required for the sweep (genome ordering / arm bins).")
        return
    gene_order = load_gene_order(args.gene_order)

    if args.materialize:
        materialize_mode(args, sc, anndata, plt, gene_order)
        return

    summary = []
    for cancer, ds, ddir in find_datasets(args.data_root, args.cancer, args.dataset):
        info = locate_files(ddir)
        if info["reason"] is not None:
            log.info("[skip] %s/%s -- %s", cancer, ds, info["reason"])
            summary.append({"cancer": cancer, "dataset": ds, "sample": "",
                            "status": "skip", "reason": info["reason"]})
            continue

        try:
            adata, scol, ccol = load_expression(info, sc)
        except Exception as e:
            log.warning("[skip] %s/%s -- load error: %s", cancer, ds, e)
            summary.append({"cancer": cancer, "dataset": ds, "sample": "",
                            "status": "error", "reason": str(e)})
            continue

        if scol is None or ccol is None:
            log.warning("[skip] %s/%s -- sample/celltype column not detected", cancer, ds)
            summary.append({"cancer": cancer, "dataset": ds, "sample": "",
                            "status": "skip", "reason": "no sample/celltype column"})
            continue

        cna = load_cna(info["cna"], adata.obs_names)
        shared_cells = list(set(adata.obs_names) & set(cna.index))
        shared_genes = list(set(adata.var_names) & set(cna.columns))
        if len(shared_cells) < args.min_malignant or len(shared_genes) < 100:
            reason = "tiny overlap (cells=%d genes=%d)" % (
                len(shared_cells), len(shared_genes))
            log.info("[skip] %s/%s -- %s", cancer, ds, reason)
            summary.append({"cancer": cancer, "dataset": ds, "sample": "",
                            "status": "skip", "reason": reason})
            continue

        adata = adata[shared_cells, shared_genes].copy()
        cna = cna.loc[shared_cells, shared_genes]
        log.info("[run]  %s/%s  cells=%d genes=%d", cancer, ds,
                 adata.n_obs, adata.n_vars)

        for sample in sorted(adata.obs[scol].astype(str).unique()):
            mask = ((adata.obs[scol].astype(str) == sample) &
                    (adata.obs[ccol].isin(args.malignant_label)))
            n_mal = int(mask.sum())
            if n_mal < args.min_malignant:
                summary.append({"cancer": cancer, "dataset": ds, "sample": sample,
                                "status": "skip",
                                "reason": "malignant=%d < %d" % (n_mal, args.min_malignant)})
                continue

            adata_s = adata[mask.values].copy()
            cna_s = cna.loc[adata_s.obs_names]
            sample_dir = os.path.join(args.out, cancer, ds, sample)
            os.makedirs(sample_dir, exist_ok=True)
            log.info("    sample %s  malignant=%d", sample, n_mal)

            try:
                if args.select == "prune":
                    rows = prune_sample(adata_s, cna_s, gene_order, args,
                                        sc, anndata, plt, sample_dir)
                else:
                    rows = sweep_sample(adata_s, cna_s, gene_order, args,
                                        sc, anndata, plt, sample_dir)
            except Exception as e:
                log.warning("    %s failed for %s: %s", args.select, sample, e)
                summary.append({"cancer": cancer, "dataset": ds, "sample": sample,
                                "status": "error", "reason": str(e)})
                continue

            sdf = pd.DataFrame(rows)
            sdf.insert(0, "sample", sample)
            sdf.insert(0, "dataset", ds)
            sdf.insert(0, "cancer", cancer)
            sdf.to_csv(os.path.join(sample_dir, "scores.csv"), index=False)
            for r in rows:
                summary.append({"cancer": cancer, "dataset": ds, "sample": sample,
                                "status": "ok", **r})

    pd.DataFrame(summary).to_csv(os.path.join(args.out, "sweep_summary.csv"),
                                 index=False)
    log.info("\nDone. Per-sample outputs under %s/<cancer>/<dataset>/<sample>/", args.out)
    log.info("Global table: %s", os.path.join(args.out, "sweep_summary.csv"))


if __name__ == "__main__":
    main()
