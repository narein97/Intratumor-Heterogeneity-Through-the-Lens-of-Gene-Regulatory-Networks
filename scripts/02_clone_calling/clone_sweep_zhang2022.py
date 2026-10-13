#!/usr/bin/env python3
"""Same clone-calling procedure for the inferCNV-derived EOC h5ads.

Usage:
    python clone_sweep_zhang2022.py --h5ad_dir h5ad_out --gene-order gencode_v25_gene_order.txt --materialize --out clones_out
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
log = logging.getLogger("clone_sweep_zhang2022")

# =============================================================================
# Dataset-agnostic helpers.
# =============================================================================

HG38_CENTROMERE_MB = {
    "chr1": 123.4, "chr2": 93.9, "chr3": 90.9, "chr4": 50.0, "chr5": 48.8,
    "chr6": 59.8, "chr7": 60.1, "chr8": 45.2, "chr9": 43.0, "chr10": 39.8,
    "chr11": 53.4, "chr12": 35.5, "chr13": 17.7, "chr14": 17.2, "chr15": 19.0,
    "chr16": 36.8, "chr17": 25.1, "chr18": 18.5, "chr19": 26.2, "chr20": 28.1,
    "chr21": 12.0, "chr22": 15.0, "chrX": 61.0, "chrY": 10.4,
}
CHROM_ORDER = ["chr%d" % i for i in range(1, 23)] + ["chrX", "chrY"]


def load_gene_order(path):
    go = pd.read_csv(path, names=range(5), sep="\t")
    go = go.rename(columns={0: "gene", 1: "chrom", 2: "start"})
    go["gene"] = go["gene"].astype(str)
    go["chrom"] = go["chrom"].astype(str)
    go["start"] = pd.to_numeric(go["start"], errors="coerce")
    go = go.dropna(subset=["start"]).drop_duplicates("gene").set_index("gene")
    return go[["chrom", "start"]]


def genome_ordered_genes(var_names, gene_order):
    common = [g for g in var_names if g in gene_order.index]
    sub = gene_order.loc[common].copy()
    sub["chrom"] = pd.Categorical(sub["chrom"], categories=CHROM_ORDER, ordered=True)
    sub = sub.dropna(subset=["chrom"]).sort_values(["chrom", "start"])
    return sub.index.tolist()


def gene_to_bin(var_names, gene_order, level="arm", bin_mb=10.0,
                centromeres=HG38_CENTROMERE_MB):
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
        else:
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
    keep = [g for g in cna_cells_by_gene.columns if g in gene_bin_map]
    sub = cna_cells_by_gene[keep]
    binlabels = pd.Index([gene_bin_map[g] for g in keep], name="bin")
    prof = sub.T.groupby(binlabels).mean().T
    prof = prof.reindex(columns=[b for b in bins if b in prof.columns])
    return prof


def distinctness(arm_prof, labels):
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
                pairs.append(1.0)
            else:
                r = float(np.clip(np.dot(Z[i], Z[j]) / b, -1.0, 1.0))
                pairs.append(1.0 - r)
    pairs = np.asarray(pairs)
    return float(np.min(pairs)), float(np.mean(pairs))


def stability_ari(X_pca, full_labels, resolution, sc, anndata,
                  reps=5, frac=0.8, n_neighbors=20, seed=0):
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


def clone_colors(labels, plt):
    uniq = sorted(pd.unique(labels),
                  key=lambda x: int(x) if str(x).isdigit() else str(x))
    pal = plt.get_cmap("tab10")
    return {c: pal(i % 10) for i, c in enumerate(uniq)}


def heatmap_panel(ax, cna_ordered, labels, plt, vcenter=1.0):
    """Genome-ordered CNA heatmap (cells sorted by clone) + clone sidebar.

    vmin/vmax are the ACTUAL min/max of this patient's data (not a fixed
    +/-5 assumption, which was tuned for 3CA's own CNA_matrix scale and
    produced an all-one-color heatmap on inferCNV's raw output, which is
    centered around 1.0 -- multiplicative scale, 1=neutral, >1=gain,
    <1=loss -- not around 0). TwoSlopeNorm is used specifically because a
    plain Normalize(vmin, vmax) always puts the colormap's center at the
    arithmetic midpoint of vmin/vmax, which won't land on 1.0 unless the
    data happens to be symmetric around it -- TwoSlopeNorm supports an
    explicit, independent pivot point instead.
    """
    import matplotlib.colors as mcolors

    order = np.argsort([int(l) if str(l).isdigit() else str(l) for l in labels])
    mat = cna_ordered[order]
    lab = np.asarray(labels)[order]

    finite = mat[np.isfinite(mat)]
    vmin = float(np.min(finite)) if len(finite) else vcenter - 1
    vmax = float(np.max(finite)) if len(finite) else vcenter + 1

    # TwoSlopeNorm requires vmin < vcenter < vmax strictly -- defend against
    # a patient whose data doesn't actually straddle 1.0 (e.g. very low
    # variance / mostly-flat CNV signal), rather than letting it crash.
    if vmin >= vcenter or vmax <= vcenter:
        ax.text(0.5, 0.5, "no variation around center=%.2f\n(min=%.3f max=%.3f)"
               % (vcenter, vmin, vmax), ha="center", va="center", fontsize=7,
               transform=ax.transAxes)
        ax.set_xticks([]); ax.set_yticks([])
        return

    norm = mcolors.TwoSlopeNorm(vmin=vmin, vcenter=vcenter, vmax=vmax)
    ax.imshow(mat, aspect="auto", cmap="bwr", norm=norm, interpolation="nearest")
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


def make_gaussian_test(n_null, max_cells, seed0):
    def test(Xa, Xb, seed=0):
        rng = np.random.default_rng(seed0 + seed)
        Xa = np.nan_to_num(np.asarray(Xa, float))
        Xb = np.nan_to_num(np.asarray(Xb, float))

        def sub(X):
            if len(X) > max_cells:
                return X[rng.choice(len(X), max_cells, replace=False)]
            return X
        Xa, Xb = sub(Xa), sub(Xb)
        if len(Xa) < 3 or len(Xb) < 3:
            return 0.0, 1.0

        X = np.vstack([Xa, Xb])
        n, d = X.shape
        cov = np.cov(X, rowvar=False) + np.eye(d) * 1e-6
        Cinv = np.linalg.pinv(cov)

        def sep(A, B):
            diff = A.mean(0) - B.mean(0)
            return float(np.sqrt(max(diff @ Cinv @ diff, 0.0)))

        obs = sep(Xa, Xb)
        mu = X.mean(0)
        null = np.empty(n_null)
        for i in range(n_null):
            sim = rng.multivariate_normal(mu, cov, size=n)
            simc = sim - sim.mean(0)
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
        a, b = pairs[0]
        obs, p = test(X[labels == a], X[labels == b], seed=seed)
        p_adj = min(1.0, p * len(pairs))
        if p_adj > alpha:
            labels[labels == b] = a
            if log_fn:
                log_fn("merge %s+%s (sep=%.2f p_adj=%.3f)" % (a, b, obs, p_adj))
        else:
            if log_fn:
                log_fn("closest pair %s,%s distinct (sep=%.2f p_adj=%.3f) -> stop"
                       % (a, b, obs, p_adj))
            break
    remap = {c: str(i) for i, c in enumerate(pd.unique(labels))}
    return np.array([remap[c] for c in labels])



def find_patient_h5ads(h5ad_dir, patient_subset=None):
    files = sorted(glob.glob(os.path.join(h5ad_dir, "*.h5ad")))
    if patient_subset:
        wanted = set(patient_subset)
        files = [f for f in files
                if os.path.basename(f).replace(".h5ad", "") in wanted]
    return files


def load_patient_malignant_cnp(h5ad_path, ad):
    """
    Load one patient's h5ad and extract the malignant-cell CNA profile,
    matching the shape/orientation clone_sweep.py's sweep_sample/_prune_core
    expect (a cells x genes DataFrame, plus a minimal AnnData "shell" with
    that same matrix in .X for scanpy's scale/PCA/neighbors/UMAP/leiden).

    Malignant = obs['clone'] != 'N' (the convention build_h5ad_per_patient.py
    already established), NOT a separate cell_type/malignant-label lookup --
    this h5ad's clone column already encodes exactly that distinction.

    Returns None (with a logged reason) if the h5ad is missing what's needed,
    rather than raising, so a batch run can skip one bad patient without
    stopping.
    """
    adata = ad.read_h5ad(h5ad_path)

    if "clone" not in adata.obs.columns:
        log.warning("  [skip] %s -- no obs['clone']", os.path.basename(h5ad_path))
        return None
    if "cnp" not in adata.obsm or "cnp_genes" not in adata.uns:
        log.warning("  [skip] %s -- no obsm['cnp']/uns['cnp_genes'] "
                    "(run update_h5ad_cnp_bootstrap_shuffle.py first)",
                    os.path.basename(h5ad_path))
        return None

    malignant_mask = (adata.obs["clone"].astype(str) != "N").to_numpy()
    n_malignant = int(malignant_mask.sum())
    if n_malignant == 0:
        log.warning("  [skip] %s -- 0 malignant cells", os.path.basename(h5ad_path))
        return None

    cnp_genes = list(adata.uns["cnp_genes"])
    cnp_all = np.asarray(adata.obsm["cnp"])
    cnp_malignant = cnp_all[malignant_mask]

    if np.isnan(cnp_malignant).all():
        log.warning("  [skip] %s -- obsm['cnp'] is all-NaN for malignant cells "
                    "(the CN profile step likely didn't find a matching "
                    "inferCNV output file for this patient -- check "
                    "update_h5ad_cnp_bootstrap_shuffle.py's logs/--cnv_filename)",
                    os.path.basename(h5ad_path))
        return None

    malignant_cells = adata.obs_names[malignant_mask]
    cna_df = pd.DataFrame(cnp_malignant, index=malignant_cells, columns=cnp_genes)

    finite_vals = cna_df.to_numpy()
    finite_vals = finite_vals[np.isfinite(finite_vals)]
    if len(finite_vals) > 0:
        log.info("    CNP value range for this patient: min=%.3f max=%.3f mean=%.3f "
                "(heatmap color scale assumes roughly -5..+5, centered at 0 -- "
                "if these numbers don't look like that, that's very likely the "
                "all-red cause, not the plotting code itself)",
                finite_vals.min(), finite_vals.max(), finite_vals.mean())

    return {
        "adata_full": adata,             # full patient object (for obs metadata, non-malignant cells)
        "malignant_mask": malignant_mask,
        "malignant_cells": malignant_cells,
        "cna_df": cna_df,                # cells x genes, malignant only
        "n_malignant": n_malignant,
    }


# =============================================================================
# NEW: treatment-phase proportions per clone -- the reporting requirement
# =============================================================================

def treatment_phase_proportions(obs, clone_labels, clone_index, phase_col="treatment_phase"):
    """
    obs: full patient obs DataFrame (indexed by cell barcode)
    clone_labels: array of clone labels, aligned to clone_index
    clone_index: cell barcodes clone_labels correspond to (e.g. malignant cells only)

    Returns a tidy DataFrame: clone, treatment_phase, n_cells, proportion
    -- proportion is within-clone (each clone's own rows sum to 1.0).
    """
    if phase_col not in obs.columns:
        log.warning("    [treatment-phase] obs has no '%s' column -- skipping", phase_col)
        return pd.DataFrame(columns=["clone", phase_col, "n_cells", "proportion"])

    df = pd.DataFrame({
        "clone": np.asarray(clone_labels).astype(str),
        phase_col: obs.loc[clone_index, phase_col].astype(str).to_numpy(),
    })
    counts = df.groupby(["clone", phase_col]).size().rename("n_cells").reset_index()
    totals = counts.groupby("clone")["n_cells"].transform("sum")
    counts["proportion"] = counts["n_cells"] / totals
    return counts.sort_values(["clone", phase_col]).reset_index(drop=True)


def print_treatment_phase_table(prop_df, patient):
    if prop_df.empty:
        return
    log.info("    Treatment-phase composition per clone (%s):", patient)
    for clone, grp in prop_df.groupby("clone"):
        parts = ", ".join(
            "%s=%.0f%%(n=%d)" % (row.treatment_phase, row.proportion * 100, row.n_cells)
            for row in grp.itertuples()
        )
        log.info("      clone %s: %s", clone, parts)


# =============================================================================
# Sweep mode (per-resolution scores, no h5ad written)
# =============================================================================

def sweep_one_patient(loaded, gene_order, args, sc, anndata, plt, sample_dir):
    cna_raw = loaded["cna_df"]
    adata_s = anndata.AnnData(X=cna_raw.to_numpy().astype(np.float32))
    adata_s.obs_names = cna_raw.index.astype(str)
    adata_s.var_names = cna_raw.columns.astype(str)

    sc.pp.scale(adata_s, max_value=10)
    n_comps = max(2, min(args.n_pcs, adata_s.n_obs - 1, adata_s.n_vars - 1))
    sc.tl.pca(adata_s, n_comps=n_comps)
    sc.pp.neighbors(adata_s, use_rep="X_pca",
                    n_neighbors=min(args.n_neighbors, adata_s.n_obs - 1))
    sc.tl.umap(adata_s, random_state=args.seed)
    X_pca = adata_s.obsm["X_pca"]
    umap_xy = adata_s.obsm["X_umap"]

    gbin, bins = gene_to_bin(cna_raw.columns, gene_order,
                             level=args.aggregate, bin_mb=args.bin_mb)
    arm = arm_profiles(cna_raw, gbin, bins)
    ordered_genes = genome_ordered_genes(cna_raw.columns, gene_order)
    cna_ordered_full = cna_raw[ordered_genes].to_numpy()

    rows, panels = [], []
    for res in args.resolutions:
        sc.tl.leiden(adata_s, resolution=res, key_added="clone", random_state=args.seed)
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
                                     reps=args.stability_reps, frac=args.stability_frac,
                                     n_neighbors=args.n_neighbors, seed=args.seed)

        rows.append({"resolution": res, "n_clones": k,
                     "min_clone_size": min(sizes.values()) if sizes else 0,
                     "clone_sizes": json.dumps(sizes), "silhouette": sil,
                     "distinct_min": d_min, "distinct_mean": d_mean, "stability_ari": stab})

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


# =============================================================================
# Prune / materialize mode
# =============================================================================

def _prune_core_one_patient(loaded, gene_order, args, sc, anndata):
    cna_raw = loaded["cna_df"]
    adata_s = anndata.AnnData(X=cna_raw.to_numpy().astype(np.float32))
    adata_s.obs_names = cna_raw.index.astype(str)
    adata_s.var_names = cna_raw.columns.astype(str)

    sc.pp.scale(adata_s, max_value=10)
    n_comps = max(2, min(args.n_pcs, adata_s.n_obs - 1, adata_s.n_vars - 1))
    sc.tl.pca(adata_s, n_comps=n_comps)
    sc.pp.neighbors(adata_s, use_rep="X_pca",
                    n_neighbors=min(args.n_neighbors, adata_s.n_obs - 1))
    sc.tl.umap(adata_s, random_state=args.seed)

    gbin, bins = gene_to_bin(cna_raw.columns, gene_order,
                             level=args.aggregate, bin_mb=args.bin_mb)
    arm = arm_profiles(cna_raw, gbin, bins).fillna(0.0)
    ordered_genes = genome_ordered_genes(cna_raw.columns, gene_order)
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


def materialize_one_patient(patient, h5ad_path, gene_order, args, sc, anndata, plt, ad):
    loaded = load_patient_malignant_cnp(h5ad_path, ad)
    if loaded is None:
        return None
    if loaded["n_malignant"] < args.min_malignant:
        log.info("  [skip] %s -- malignant=%d < %d", patient, loaded["n_malignant"], args.min_malignant)
        return None

    log.info("=== %s ===  malignant=%d", patient, loaded["n_malignant"])
    core = _prune_core_one_patient(loaded, gene_order, args, sc, anndata)
    labels = pd.Series(core["labels"], index=core["obs_names"])

    sizes = labels.value_counts()
    kept = [c for c, n in sizes.items() if n >= args.min_clone_cells]
    if len(kept) < 1:
        log.info("  [skip] %s -- 0 clone(s) >= %d cells", patient, args.min_clone_cells)
        return None

    # Relabel kept clones contiguously (largest -> "0")
    order = sizes.loc[kept].sort_values(ascending=False).index.tolist()
    remap = {old: str(i) for i, old in enumerate(order)}
    keep_mal = labels[labels.isin(kept)].index
    new_labels = labels.loc[keep_mal].map(remap)

    adata_full = loaded["adata_full"]
    clone_leiden = pd.Series(np.nan, index=adata_full.obs_names, dtype=object)
    clone_leiden.loc[adata_full.obs_names[~loaded["malignant_mask"]]] = "N"
    clone_leiden.loc[new_labels.index] = new_labels.values
    adata_full.obs["clone_leiden"] = clone_leiden.values

    n_excluded_small = int(adata_full.obs["clone_leiden"].isna().sum())
    if n_excluded_small > 0:
        log.info("  Dropping %d malignant cell(s) in clones < %d cells "
                "(not kept with clone_leiden=NaN)", n_excluded_small, args.min_clone_cells)
        adata_full = adata_full[adata_full.obs["clone_leiden"].notna()].copy()

    # Treatment-phase composition per new clone (malignant clones only --
    # 'N' stromal cells aren't a "clone" in the sense this is meant to report on)
    prop_df = treatment_phase_proportions(adata_full.obs, new_labels.values, new_labels.index)
    print_treatment_phase_table(prop_df, patient)

    out_path = os.path.join(args.out, "%s.h5ad" % patient)
    os.makedirs(args.out, exist_ok=True)
    adata_full.write_h5ad(out_path)
    log.info("  -> %s  (%d clone(s): %s, %d cells total after drop)", out_path, len(kept),
             dict(sizes.loc[kept]), adata_full.n_obs)

    if not args.no_plots:
        try:
            sample_dir = os.path.join(args.out, "plots", patient)
            os.makedirs(sample_dir, exist_ok=True)
            ordered_genes = genome_ordered_genes(loaded["cna_df"].columns, gene_order)
            cna_ord = loaded["cna_df"].loc[new_labels.index, ordered_genes].to_numpy()
            pos = {c: i for i, c in enumerate(core["obs_names"])}
            umap_idx = [pos[c] for c in new_labels.index]
            panel = {"res": args.prune_resolution, "umap_xy": core["umap_xy"][umap_idx],
                     "labels": new_labels.values, "cna_ordered": cna_ord,
                     "scores": "k=%d (from %d)" % (len(kept), core["n_init"])}
            make_montages(sample_dir, [panel], plt)
        except Exception as e:
            log.warning("  montage failed: %s", e)

    manifest_rows = []
    for clone, n in sizes.loc[kept].items():
        new_id = remap[clone]
        clone_prop = prop_df[prop_df["clone"] == new_id]
        manifest_rows.append({
            "patient": patient, "clone": new_id, "n_cells": int(n),
            "treatment_phase_breakdown": json.dumps(
                dict(zip(clone_prop["treatment_phase"], clone_prop["proportion"].round(3)))),
            "h5ad": os.path.abspath(out_path),
        })
    return manifest_rows


# =============================================================================
# main
# =============================================================================

def build_argparser():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--h5ad_dir", required=True,
                   help="directory of per-patient h5ads (build_h5ad_per_patient.py output, "
                        "already updated with obsm['cnp']/uns['cnp_genes'] via "
                        "update_h5ad_cnp_bootstrap_shuffle.py)")
    p.add_argument("--patients", nargs="+", default=None,
                   help="restrict to these patient IDs (filenames without .h5ad)")
    p.add_argument("--out", default="clone_sweep_zhang2022_out")
    p.add_argument("--gene-order", default=None,
                   help="gencode-style table (gene, chrom, start) -- e.g. "
                        "gencode_v25_gene_order.txt. Required unless --inspect.")
    p.add_argument("--resolutions", nargs="+", type=float,
                   default=[0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0])
    p.add_argument("--min-malignant", type=int, default=50,
                   help="min malignant cells in a patient to cluster it")
    p.add_argument("--n-pcs", type=int, default=20)
    p.add_argument("--n-neighbors", type=int, default=20)
    p.add_argument("--aggregate", choices=["arm", "chromosome", "bin"], default="arm")
    p.add_argument("--bin-mb", type=float, default=10.0)
    p.add_argument("--stability-reps", type=int, default=5)
    p.add_argument("--stability-frac", type=float, default=0.8)
    p.add_argument("--skip-stability", action="store_true")
    p.add_argument("--no-plots", action="store_true")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--select", choices=["sweep", "prune"], default="sweep")
    p.add_argument("--prune-resolution", type=float, default=1.0)
    p.add_argument("--prune-alpha", type=float, default=0.05,
                   help="RAISE for MORE clones (less merging), lower for fewer.")
    p.add_argument("--prune-perm", type=int, default=200)
    p.add_argument("--prune-max-cells", type=int, default=400)
    p.add_argument("--materialize", action="store_true",
                   help="prune each patient into clones, write clone_leiden into "
                        "the h5ad (alongside the existing inferCNV-derived obs['clone']), "
                        "and report treatment-phase composition per clone")
    p.add_argument("--min-clone-cells", type=int, default=50)
    p.add_argument("--inspect", action="store_true",
                   help="report what's found per patient h5ad; no clustering")
    return p


def inspect_mode(args, ad):
    files = find_patient_h5ads(args.h5ad_dir, args.patients)
    log.info("Found %d h5ad(s) in %s", len(files), args.h5ad_dir)
    for f in files:
        patient = os.path.basename(f).replace(".h5ad", "")
        try:
            adata = ad.read_h5ad(f)
        except Exception as e:
            log.info("  %s -- READ ERROR: %s", patient, e)
            continue
        has_clone = "clone" in adata.obs.columns
        has_cnp = "cnp" in adata.obsm and "cnp_genes" in adata.uns
        n_mal = int((adata.obs["clone"].astype(str) != "N").sum()) if has_clone else -1
        has_phase = "treatment_phase" in adata.obs.columns
        log.info("  %-12s cells=%-6d clone_col=%-5s cnp=%-5s malignant=%-6d treatment_phase_col=%s",
                 patient, adata.n_obs, has_clone, has_cnp, n_mal, has_phase)


def main():
    args = build_argparser().parse_args()
    warnings.filterwarnings("ignore")

    import anndata as ad
    import scanpy as sc
    import anndata
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    sc.settings.verbosity = 0
    np.random.seed(args.seed)
    os.makedirs(args.out, exist_ok=True)

    if args.inspect:
        inspect_mode(args, ad)
        return

    if not args.gene_order:
        log.error("--gene-order is required (genome ordering / arm bins).")
        return
    gene_order = load_gene_order(args.gene_order)

    files = find_patient_h5ads(args.h5ad_dir, args.patients)
    log.info("Found %d h5ad(s) to process", len(files))

    if args.materialize:
        manifest = []
        for f in files:
            patient = os.path.basename(f).replace(".h5ad", "")
            try:
                rows = materialize_one_patient(patient, f, gene_order, args, sc, anndata, plt, ad)
            except Exception as e:
                log.warning("[error] %s -- %s", patient, e)
                continue
            if rows:
                manifest.extend(rows)
        man_path = os.path.join(args.out, "clones_manifest.csv")
        pd.DataFrame(manifest).to_csv(man_path, index=False)
        log.info("\nWrote %d patient h5ad(s), %d clone(s) total. Manifest: %s",
                 len({r["patient"] for r in manifest}), len(manifest), man_path)
        return

    summary = []
    for f in files:
        patient = os.path.basename(f).replace(".h5ad", "")
        loaded = load_patient_malignant_cnp(f, ad)
        if loaded is None:
            summary.append({"patient": patient, "status": "skip"})
            continue
        if loaded["n_malignant"] < args.min_malignant:
            log.info("[skip] %s -- malignant=%d < %d", patient, loaded["n_malignant"], args.min_malignant)
            summary.append({"patient": patient, "status": "skip",
                            "reason": "malignant=%d < %d" % (loaded["n_malignant"], args.min_malignant)})
            continue

        sample_dir = os.path.join(args.out, patient)
        os.makedirs(sample_dir, exist_ok=True)
        log.info("=== %s ===  malignant=%d", patient, loaded["n_malignant"])
        try:
            if args.select == "prune":
                core = _prune_core_one_patient(loaded, gene_order, args, sc, anndata)
                rows = [{"n_init_clusters": core["n_init"],
                        "n_clones": len(np.unique(core["labels"]))}]
            else:
                rows = sweep_one_patient(loaded, gene_order, args, sc, anndata, plt, sample_dir)
        except Exception as e:
            log.warning("  %s failed for %s: %s", args.select, patient, e)
            summary.append({"patient": patient, "status": "error", "reason": str(e)})
            continue

        sdf = pd.DataFrame(rows)
        sdf.insert(0, "patient", patient)
        sdf.to_csv(os.path.join(sample_dir, "scores.csv"), index=False)
        for r in rows:
            summary.append({"patient": patient, "status": "ok", **r})

    pd.DataFrame(summary).to_csv(os.path.join(args.out, "sweep_summary.csv"), index=False)
    log.info("\nDone. Per-patient outputs under %s/<patient>/", args.out)
    log.info("Global table: %s", os.path.join(args.out, "sweep_summary.csv"))


if __name__ == "__main__":
    main()
