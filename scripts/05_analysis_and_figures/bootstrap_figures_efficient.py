#!/usr/bin/env python3
"""Bootstrap support, clone-vs-normal Jaccard, intra/inter Hamming and clone-size Wasserstein figures over the cohort.

Usage:
    python bootstrap_figures_efficient.py --root clones_out --out bootstrap_figs_out
"""

import argparse
import glob
import itertools
import json
import os
import logging
import multiprocessing as mp
import textwrap

import numpy as np
import pandas as pd
import anndata as ad
import networkx as nx
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from scipy.stats import mannwhitneyu, wasserstein_distance, binomtest, wilcoxon, spearmanr

from intra_vs_inter_efficient import (
    edge_sets_to_matrix,
    pairwise_hamming,
    jaccard,
    cancer_type_from_path,
)

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)-7s %(message)s",
                    datefmt="%H:%M:%S")
log = logging.getLogger("bootstrap_figures")

PLOT_DPI = 300


def _wrap_title(text, width=38):
    """Wrap long plot titles (which interpolate a cancer-type name and can
    run past a narrow figure's width) onto multiple lines instead of
    letting matplotlib clip them at the figure edge."""
    return "\n".join(textwrap.wrap(text, width=width))


# =============================================================================
# Data loading -- decode ONCE per patient, share across all 4 figures
# =============================================================================

def decode_and_filter(arr, columns, genes, clone_categories, p_thresh):
    df = pd.DataFrame(arr, columns=columns)
    df["clone"] = clone_categories[df["clone"].astype(int)]
    df["source"] = genes[df["source"].astype(int)]
    df["target"] = genes[df["target"].astype(int)]
    return df[df["p"] < p_thresh]


def to_per_clone_edge_sets(df, all_clones):
    out = {c: set() for c in all_clones}
    for clone, grp in df.groupby("clone"):
        if clone in out:
            out[clone] = set(zip(grp["source"], grp["target"]))
    return out


def to_per_clone_weighted_edges(df, all_clones):
    """{clone: {(source,target): coef_mean}} -- same p<0.05-filtered edges
    as to_per_clone_edge_sets, but keeping the coef_mean weight (needed only
    by Figure 4's Wasserstein divergence)."""
    out = {c: {} for c in all_clones}
    for clone, grp in df.groupby("clone"):
        if clone in out:
            out[clone] = dict(zip(zip(grp["source"], grp["target"]), grp["coef_mean"]))
    return out


def load_patient_full(h5ad_path, p_thresh):
    """
    One h5ad read, decoded once, shared by all 4 figures. Returns None only
    if the file lacks what EVERY figure needs; otherwise returns a dict
    whose 'has_bootstrap' / 'has_clone_sizes' flags tell each figure's
    per-patient function whether ITS extra requirement was met.
    """
    adata = ad.read_h5ad(h5ad_path, backed="r")
    core_required = ["celloracle_original", "celloracle_columns", "celloracle_genes",
                     "clone_categories"]
    if any(k not in adata.uns for k in core_required):
        return None

    columns = list(adata.uns["celloracle_columns"])
    genes = np.asarray(adata.uns["celloracle_genes"]).astype(str)
    clone_categories = np.asarray(adata.uns["clone_categories"]).astype(str)
    all_clones = list(clone_categories)

    original_df = decode_and_filter(adata.uns["celloracle_original"], columns, genes,
                                    clone_categories, p_thresh)
    original_sets = to_per_clone_edge_sets(original_df, all_clones)
    original_weighted = to_per_clone_weighted_edges(original_df, all_clones)

    has_bootstrap = "celloracle_bootstrap" in adata.uns
    bootstrap_sets = None
    if has_bootstrap:
        rep_keys_sorted = sorted(adata.uns["celloracle_bootstrap"].keys(), key=int)
        bootstrap_sets = {c: [] for c in all_clones}
        for rk in rep_keys_sorted:
            df = decode_and_filter(adata.uns["celloracle_bootstrap"][rk], columns, genes,
                                   clone_categories, p_thresh)
            per_clone = to_per_clone_edge_sets(df, all_clones)
            for c in all_clones:
                bootstrap_sets[c].append(per_clone[c])

    clone_counts = {}
    has_clone_sizes = "clone" in adata.obs
    if has_clone_sizes:
        clone_counts = adata.obs["clone"].astype(str).value_counts().to_dict()

    return {
        "clones": all_clones,
        "original_sets": original_sets,
        "original_weighted": original_weighted,
        "has_bootstrap": has_bootstrap,
        "bootstrap_sets": bootstrap_sets,
        "has_clone_sizes": has_clone_sizes,
        "clone_counts": clone_counts,
    }


# =============================================================================
# Figure 1: bootstrap support, unique vs shared edges
# =============================================================================

def compute_support_unique_shared(original_sets, bootstrap_sets, all_clones):
    """Bootstrap support of unique vs shared edges per non-normal clone."""
    unique_support, shared_support = [], []
    for clone in all_clones:
        if clone == "N":
            continue
        other_edges = set()
        for c2 in all_clones:
            if c2 != clone:
                other_edges |= original_sets[c2]
        reps = bootstrap_sets[clone]
        n_boot = len(reps)
        for e in original_sets[clone]:
            support = sum(1 for rep in range(n_boot) if e in reps[rep])
            is_unique = e not in other_edges
            (unique_support if is_unique else shared_support).append(support)
    return unique_support, shared_support


def compute_support_patient_stats(unique_support, shared_support):
    """cells 26/29: per-patient Wasserstein distance and Mann-Whitney U
    p-value (alternative='less': unique stochastically < shared) between
    the unique- and shared-edge support distributions. None if either side
    is empty (degenerate patient -- e.g. every edge is unique, or none are)."""
    if not unique_support or not shared_support:
        return None, None
    wdist = wasserstein_distance(unique_support, shared_support)
    _, pval = mannwhitneyu(unique_support, shared_support, alternative="less")
    return wdist, pval


def plot_support_boxplot(unique_vals, shared_vals, title, out_path):
    rows = ([{"support": v, "group": "unique"} for v in unique_vals] +
            [{"support": v, "group": "shared"} for v in shared_vals])
    df = pd.DataFrame(rows)
    pval = np.nan
    if unique_vals and shared_vals:
        _, pval = mannwhitneyu(unique_vals, shared_vals, alternative="less")

    plt.figure(figsize=(6, 6.5))
    ax = sns.boxplot(data=df, x="group", y="support", palette="tab10", showfliers=False)
    ax.tick_params(left=True)
    subtitle = f"MWU (unique<shared) p={pval:.2e}" if not np.isnan(pval) else "MWU: n/a"
    ax.set_title(f"{_wrap_title(title)}\n{subtitle}")
    ax.set_xlabel("")
    ax.set_ylabel("Bootstrap support (/100)")
    plt.tight_layout()
    plt.savefig(out_path, dpi=PLOT_DPI)
    plt.close()
    log.info("    %s: n_unique=%d n_shared=%d MWU p=%s", title, len(unique_vals), len(shared_vals),
             f"{pval:.3e}" if not np.isnan(pval) else "n/a")


def plot_per_patient_boxplot_by_cancer(values_by_cancer, ylabel, title, out_path, log_scale=False):
    """Generic 'one box per cancer type' plot, reused for Figure 1's
    per-patient Wasserstein-distance and MWU-p-value distributions, and for
    Figure 4's clone-size-divergence boxplot. Same overlapping-label fix as
    intra_vs_inter_efficient.py's stddev boxplot: rotated labels + a figure
    width that scales with the number of categories."""
    cancers = [c for c, vals in values_by_cancer.items() if len(vals) > 0]
    if not cancers:
        log.warning("    [skip plot] %s -- no non-empty categories", title)
        return
    n_cancer = max(len(cancers), 1)
    fig_w = max(7, 0.55 * n_cancer)
    plt.figure(figsize=(fig_w, 6))
    sns.boxplot(data=[values_by_cancer[c] for c in cancers], palette="tab10")
    ax = plt.gca()
    ax.set_xticks(range(len(cancers)))
    ax.set_xticklabels(cancers, rotation=45, ha="right", rotation_mode="anchor")
    if log_scale:
        ax.set_yscale("log")
    plt.ylabel(ylabel)
    plt.xlabel("Cancer type")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(out_path, dpi=PLOT_DPI)
    plt.close()


def plot_grouped_boxplot_by_cancer(values_a_by_cancer, values_b_by_cancer, label_a, label_b,
                                   ylabel, title, out_path, palette=None):
    """One combined boxplot: cancer type on the x-axis, groups side by side."""
    cancers = sorted(set(values_a_by_cancer) | set(values_b_by_cancer))
    cancers = [c for c in cancers if values_a_by_cancer.get(c) or values_b_by_cancer.get(c)]
    if not cancers:
        log.warning("    [skip plot] %s -- no data", title)
        return
    rows = []
    for c in cancers:
        for v in values_a_by_cancer.get(c, []):
            rows.append({"cancer": c, "value": v, "group": label_a})
        for v in values_b_by_cancer.get(c, []):
            rows.append({"cancer": c, "value": v, "group": label_b})
    df = pd.DataFrame(rows)

    n_cancer = max(len(cancers), 1)
    fig_w = max(8, 0.9 * n_cancer)
    plt.figure(figsize=(fig_w, 6.5))
    ax = sns.boxplot(data=df, x="cancer", y="value", hue="group", order=cancers,
                     palette=palette or "tab10", showfliers=False)
    ax.set_xticks(range(len(cancers)))
    ax.set_xticklabels(cancers, rotation=45, ha="right", rotation_mode="anchor")
    plt.ylabel(ylabel)
    plt.xlabel("Cancer type")
    plt.title(_wrap_title(title, width=60))
    plt.legend(title="", bbox_to_anchor=(1.02, 1), loc="upper left", borderaxespad=0.)
    plt.tight_layout()
    plt.savefig(out_path, dpi=PLOT_DPI)
    plt.close()
    log.info("    %s: %d cancer type(s) plotted", title, len(cancers))


# =============================================================================
# Figure 2: clone-vs-clone vs normal-vs-clone Jaccard
# =============================================================================

def compute_clone_normal_jaccard(original_sets, all_clones):
    """Mean pairwise Jaccard over non-normal clone pairs and over normal-vs-clone pairs."""
    clone_groups = [c for c in all_clones if c != "N"]

    clone_jaccards = []
    for g1, g2 in itertools.combinations(clone_groups, 2):
        j = jaccard(original_sets[g1], original_sets[g2])
        if not np.isnan(j):
            clone_jaccards.append(j)
    avg_clone = float(np.mean(clone_jaccards)) if clone_jaccards else None

    avg_norm = None
    if "N" in all_clones:
        norm_jaccards = []
        for g in clone_groups:
            j = jaccard(original_sets[g], original_sets["N"])
            if not np.isnan(j):
                norm_jaccards.append(j)
        avg_norm = float(np.mean(norm_jaccards)) if norm_jaccards else None

    return avg_clone, avg_norm


def plot_clone_normal_jaccard_boxplot(clone_vals, norm_vals, title, out_path):
    rows = ([{"value": v, "group": "clone v clone"} for v in clone_vals] +
            [{"value": v, "group": "normal v clone"} for v in norm_vals])
    df = pd.DataFrame(rows)
    pval = np.nan
    if clone_vals and norm_vals:
        _, pval = mannwhitneyu(clone_vals, norm_vals, alternative="greater")

    plt.figure(figsize=(6, 6.5))
    ax = sns.boxplot(data=df, x="group", y="value", palette=sns.color_palette("tab10")[4:6])
    ax.tick_params(left=True)
    subtitle = f"MWU (clone>normal) p={pval:.2e}" if not np.isnan(pval) else "MWU: n/a"
    ax.set_title(f"{_wrap_title(title)}\n{subtitle}")
    ax.set_xlabel("")
    ax.set_ylabel("Pairwise Jaccard similarity")
    plt.tight_layout()
    plt.savefig(out_path, dpi=PLOT_DPI)
    plt.close()
    log.info("    %s: n_clone=%d n_normal=%d MWU p=%s", title, len(clone_vals), len(norm_vals),
             f"{pval:.3e}" if not np.isnan(pval) else "n/a")


# =============================================================================
# Figure 3: intra vs inter category Hamming distance (original GRN only)
# =============================================================================

def compute_intra_inter_original(original_sets, all_clones):
    """Intra- vs inter-clone Hamming distances from the original (non-bootstrap) GRN."""
    clone_groups = [c for c in all_clones if c != "N"]
    all_edges = set()
    for c in all_clones:
        all_edges |= original_sets[c]
    if not all_edges:
        return [], []
    edge_index = {e: i for i, e in enumerate(sorted(all_edges))}
    n_edges = len(edge_index)

    vectors = {c: edge_sets_to_matrix([original_sets[c]], edge_index) for c in all_clones}

    intra = []
    for c1, c2 in itertools.combinations(clone_groups, 2):
        h = pairwise_hamming(vectors[c1], vectors[c2], n_edges)
        intra.append(float(h[0, 0]))

    inter = []
    if "N" in all_clones:
        for c in clone_groups:
            h = pairwise_hamming(vectors[c], vectors["N"], n_edges)
            inter.append(float(h[0, 0]))

    return intra, inter


def plot_split_violin_with_stars(intra_dict, inter_dict, title, out_path):
    """Split violin of two groups per patient with per-patient median stars."""
    star_colors = ["#AA5042", "#D5A021", "#1098F7", "#061A40", "#631A86",
                   "#226F54", "#CE4865", "#C6AD94", "#E3D7FF"]

    rows = []
    for patient, vals in intra_dict.items():
        for v in vals:
            rows.append({"distance": v, "type": "intra", "pos": 0})
    for patient, vals in inter_dict.items():
        for v in vals:
            rows.append({"distance": v, "type": "inter", "pos": 0})
    df = pd.DataFrame(rows)
    if df.empty:
        log.warning("    [skip plot] %s -- no data", title)
        return None

    group1 = df.loc[df["type"] == "intra", "distance"]
    group2 = df.loc[df["type"] == "inter", "distance"]
    pval = None
    if len(group1) > 0 and len(group2) > 0:
        _, pval = mannwhitneyu(group1, group2)

    plt.figure(figsize=(8, 7))
    vp = sns.violinplot(x="pos", y="distance", hue="type", data=df, split=True,
                        inner="quart", fill=False,
                        palette={"intra": "#9372B2", "inter": "#845B53"}, legend=False)
    vp.tick_params(left=True)

    for i, patient in enumerate(intra_dict.keys()):
        color = star_colors[i % len(star_colors)]
        temp_rows = []
        if intra_dict.get(patient):
            temp_rows.append({"pos": 0, "distance": float(np.mean(intra_dict[patient])), "type": "intra"})
        if inter_dict.get(patient):
            temp_rows.append({"pos": 0, "distance": float(np.mean(inter_dict[patient])), "type": "inter"})
        if not temp_rows:
            continue
        temp_df = pd.DataFrame(temp_rows)
        star_ax = sns.boxplot(x="pos", y="distance", hue="type", data=temp_df,
                              showfliers=False, showbox=False, width=0.05,
                              linecolor=color, linewidth=2, legend=False,
                              medianprops=dict(label="median_star"))
        for line in star_ax.get_lines():
            if line.get_label() == "median_star":
                line.set_linestyle("")
                line.set_marker("*")
                line.set_markersize(9)

    legend_handles = [mpatches.Patch(color="#9372B2", label="Intra (clone-clone)"),
                      mpatches.Patch(color="#845B53", label="Inter (clone-Normal)")]
    plt.legend(handles=legend_handles, title=f"n={len(intra_dict)} patients\n(stars = per-patient medians)",
              bbox_to_anchor=(1.05, 1), loc="upper left", borderaxespad=0.)

    subtitle = f"MWU p={pval:.2e}" if pval is not None else "MWU: n/a"
    plt.title(f"{_wrap_title(title, width=45)}\n{subtitle}")
    plt.ylabel("Hamming distance")
    plt.xlabel("")
    plt.tight_layout()
    plt.savefig(out_path, dpi=PLOT_DPI)
    plt.close()
    log.info("    %s: n_patients=%d MWU p=%s", title, len(intra_dict),
             f"{pval:.3e}" if pval is not None else "n/a")
    return float(pval) if pval is not None else None


# =============================================================================
# Figure 4: clone size vs GRN divergence (Wasserstein)
# =============================================================================

def compute_wasserstein(G1, G2, weight, s):
    """Wasserstein distance between the edge-weight distributions of source s in two GRNs."""
    t1 = {v for _, v in G1.out_edges(s)}
    t2 = {v for _, v in G2.out_edges(s)}
    targets = list(t1 | t2)

    v1 = []
    v2 = []

    for t in targets:
        w1 = G1[s][t][weight] if G1.has_edge(s, t) else 0
        w2 = G2[s][t][weight] if G2.has_edge(s, t) else 0

        v1.append(w1)
        v2.append(w2)

    v1 = np.array(v1, dtype=float)
    v2 = np.array(v2, dtype=float)

    if v1.sum() > 0:
        v1 = v1 / np.sum(np.abs(v1))
    if v2.sum() > 0:
        v2 = v2 / np.sum(np.abs(v2))

    return np.sum(np.abs(v1 - v2))


def clone_wasserstein(G1, G2, weight="coef_mean"):
    """Mean per-source Wasserstein distance between two clone GRNs."""
    src1 = {u for u, _ in G1.edges()}
    src2 = {u for u, _ in G2.edges()}

    common_sources = src1 & src2
    if len(common_sources) == 0:
        return 0

    total = 0
    for s in common_sources:
        s_dist = compute_wasserstein(G1, G2, weight, s)
        total += s_dist

    return total / len(common_sources)


def build_weighted_graph(weighted_dict):
    G = nx.DiGraph()
    for (s, t), w in weighted_dict.items():
        G.add_edge(s, t, coef_mean=w)
    return G


def compute_clone_size_divergence(original_weighted, clone_counts, all_clones):
    """cells 51/53: pick the largest and smallest non-N clone by cell
    count, compute each one's Wasserstein divergence from the Normal GRN.
    Returns None if there's no 'N' clone, fewer than 2 non-N clones with a
    nonzero cell count, or the largest/smallest are the same clone."""
    if "N" not in all_clones or "N" not in clone_counts:
        return None
    clone_groups = [c for c in all_clones if c != "N"]
    counts = {c: clone_counts.get(c, 0) for c in clone_groups}
    counts = {c: v for c, v in counts.items() if v > 0}
    if len(counts) < 2:
        return None
    clone_large = max(counts, key=counts.get)
    clone_small = min(counts, key=counts.get)
    if clone_large == clone_small:
        return None

    G_large = build_weighted_graph(original_weighted[clone_large])
    G_small = build_weighted_graph(original_weighted[clone_small])
    G_N = build_weighted_graph(original_weighted["N"])

    W_large = clone_wasserstein(G_large, G_N)
    W_small = clone_wasserstein(G_small, G_N)

    return {"D": W_large - W_small,
            "size_large": counts[clone_large], "W_large": W_large,
            "size_small": counts[clone_small], "W_small": W_small,
            "size_normal": clone_counts.get("N", 0)}


def plot_large_small_boxplot(large_vals, small_vals, title, out_path):
    if not large_vals and not small_vals:
        log.warning("    [skip plot] %s -- no data", title)
        return
    rows = ([{"value": v, "group": "large clone"} for v in large_vals] +
            [{"value": v, "group": "small clone"} for v in small_vals])
    df = pd.DataFrame(rows)
    plt.figure(figsize=(6, 6.5))
    ax = sns.boxplot(data=df, x="group", y="value", palette=sns.color_palette("husl", 8)[1:3])
    ax.tick_params(left=True)
    ax.set_title(_wrap_title(title))
    ax.set_xlabel("")
    ax.set_ylabel("Wasserstein distance to Normal")
    plt.tight_layout()
    plt.savefig(out_path, dpi=PLOT_DPI)
    plt.close()


def plot_clone_size_bars(patient_sizes, cancer, out_dir):
    """Clone-size bar charts (stacked proportion and absolute count)."""
    if not patient_sizes:
        log.warning("    [skip plot] clone size bars (%s) -- no data", cancer)
        return
    df = pd.DataFrame(patient_sizes).T
    df = df[["normal", "large", "small"]]
    colors = ["#4C72B0", "#DD8452", "#55A868"]

    n_patient = len(df)
    fig_w = max(7, 0.5 * n_patient)

    prop_df = df.div(df.sum(axis=1), axis=0)
    plt.figure(figsize=(fig_w, 6))
    ax = prop_df.plot.bar(stacked=True, ax=plt.gca(), color=colors)
    ax.legend(title="Clone size class", loc="upper left", bbox_to_anchor=(1.0, 1.0))
    plt.ylabel("Proportion of cells")
    plt.xlabel("Patient")
    plt.xticks(rotation=45, ha="right", rotation_mode="anchor")
    plt.title(_wrap_title(f"Clone-size proportion of cells per patient ({cancer})", width=50))
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, f"clone_size_proportion_{cancer}.png"), dpi=PLOT_DPI)
    plt.close()

    plt.figure(figsize=(fig_w, 6))
    ax = df.plot.bar(ax=plt.gca(), color=colors)
    ax.locator_params(axis="y", nbins=20)
    ax.legend(title="Clone size class", loc="upper left", bbox_to_anchor=(1.0, 1.0))
    plt.ylabel("Cell count")
    plt.xlabel("Patient")
    plt.xticks(rotation=45, ha="right", rotation_mode="anchor")
    plt.title(_wrap_title(f"Clone-size cell count per patient ({cancer})", width=50))
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, f"clone_size_count_{cancer}.png"), dpi=PLOT_DPI)
    plt.close()
    log.info("    Clone-size bars (%s): %d patient(s)", cancer, n_patient)


# =============================================================================
# Per-patient worker (runs in a separate process)
# =============================================================================

def process_patient(task):
    """All per-patient math for all 4 figures, done once per h5ad read.
    Only load+compute happens here; every plt.* call and every cohort-level
    statistical test stays in the parent process (matplotlib figures can't
    cross a process boundary, and pool.map's order-preserving return keeps
    the parent's aggregation identical to a serial run -- see main())."""
    f, root, p_thresh = task
    patient = os.path.basename(f).replace("_clones.h5ad", "")
    cancer = cancer_type_from_path(f, root)

    loaded = load_patient_full(f, p_thresh)
    if loaded is None:
        return {"file": f, "patient": patient, "cancer": cancer, "skipped": True}

    all_clones = loaded["clones"]
    result = {"file": f, "patient": patient, "cancer": cancer, "skipped": False}

    # Figure 1
    if loaded["has_bootstrap"]:
        unique_supp, shared_supp = compute_support_unique_shared(
            loaded["original_sets"], loaded["bootstrap_sets"], all_clones)
        wdist, pval = compute_support_patient_stats(unique_supp, shared_supp)
        result["fig1"] = {"unique_support": unique_supp, "shared_support": shared_supp,
                          "wasserstein": wdist, "mwu_pvalue": pval}
    else:
        result["fig1"] = None

    # Figure 2
    avg_clone_jacc, avg_norm_jacc = compute_clone_normal_jaccard(loaded["original_sets"], all_clones)
    result["fig2"] = {"clone_jacc": avg_clone_jacc, "norm_jacc": avg_norm_jacc}

    # Figure 3
    intra, inter = compute_intra_inter_original(loaded["original_sets"], all_clones)
    result["fig3"] = {"intra": intra, "inter": inter}

    # Figure 4
    if loaded["has_clone_sizes"]:
        div = compute_clone_size_divergence(loaded["original_weighted"], loaded["clone_counts"], all_clones)
        result["fig4"] = div
    else:
        result["fig4"] = None

    return result


# =============================================================================
# main
# =============================================================================

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--root", required=True)
    p.add_argument("--out", default="bootstrap_figs_out")
    p.add_argument("--p-thresh", type=float, default=0.05)
    p.add_argument("--exclude-cancer-type", nargs="*", default=["Temporal"],
                   help="Cancer-type directory name(s) to skip (path segment, e.g. Temporal).")
    p.add_argument("--n-workers", type=int, default=max(1, (os.cpu_count() or 2) - 1),
                   help="Local worker processes (1 = sequential).")
    args = p.parse_args()
    os.makedirs(args.out, exist_ok=True)
    for sub in ["fig1_bootstrap_support", "fig2_clone_normal_jaccard",
               "fig3_intra_inter_hamming", "fig4_clone_size_divergence"]:
        os.makedirs(os.path.join(args.out, sub), exist_ok=True)

    files = sorted(glob.glob(f"{args.root}/**/*_clones.h5ad", recursive=True))
    files = [f for f in files if "_with_support" not in f]
    if args.exclude_cancer_type:
        files = [f for f in files if not any(f"/{ct}/" in f for ct in args.exclude_cancer_type)]
    log.info("Processing %d patient(s) with %d worker process(es)", len(files), args.n_workers)

    tasks = [(f, args.root, args.p_thresh) for f in files]
    if args.n_workers > 1 and len(tasks) > 1:
        with mp.Pool(processes=min(args.n_workers, len(tasks))) as pool:
            # order-preserving map -- see process_patient()'s docstring and
            # intra_vs_inter_efficient.py's main() for why this matters.
            results = pool.map(process_patient, tasks)
    else:
        results = [process_patient(t) for t in tasks]

    # ---- Aggregate, in the same fixed order as `files`, by cancer type ----
    fig1_unique_by_cancer, fig1_shared_by_cancer = {}, {}
    fig1_wdist_by_cancer, fig1_pval_by_cancer = {}, {}
    fig2_clone_by_cancer, fig2_norm_by_cancer = {}, {}
    fig3_intra_dict_by_cancer, fig3_inter_dict_by_cancer = {}, {}  # {cancer: {patient_key: [values]}}
    fig4_large_by_cancer, fig4_small_by_cancer = {}, {}
    fig4_D_all, fig4_sizes_all, fig4_Ws_all = [], [], []
    fig4_patient_sizes_by_cancer = {}  # {cancer: {patient: {'normal','large','small'}}}

    n_ok = 0
    for res in results:
        if res["skipped"]:
            log.warning("  [skip] %s -- missing required .uns keys", res["file"])
            continue
        n_ok += 1
        patient, cancer = res["patient"], res["cancer"]

        if res["fig1"] is not None:
            fig1_unique_by_cancer.setdefault(cancer, []).extend(res["fig1"]["unique_support"])
            fig1_shared_by_cancer.setdefault(cancer, []).extend(res["fig1"]["shared_support"])
            if res["fig1"]["wasserstein"] is not None:
                fig1_wdist_by_cancer.setdefault(cancer, []).append(res["fig1"]["wasserstein"])
                fig1_pval_by_cancer.setdefault(cancer, []).append(res["fig1"]["mwu_pvalue"])

        f2 = res["fig2"]
        if f2["clone_jacc"] is not None:
            fig2_clone_by_cancer.setdefault(cancer, []).append(f2["clone_jacc"])
        if f2["norm_jacc"] is not None:
            fig2_norm_by_cancer.setdefault(cancer, []).append(f2["norm_jacc"])

        f3 = res["fig3"]
        if f3["intra"] or f3["inter"]:
            # key carries the cancer type so merging every cancer's dict
            # together for the "all patients" panel can never collide two
            # different cancer types' same-named patient IDs.
            key = f"{cancer}::{patient}"
            fig3_intra_dict_by_cancer.setdefault(cancer, {})[key] = f3["intra"]
            fig3_inter_dict_by_cancer.setdefault(cancer, {})[key] = f3["inter"]

        f4 = res["fig4"]
        if f4 is not None:
            fig4_large_by_cancer.setdefault(cancer, []).append(f4["W_large"])
            fig4_small_by_cancer.setdefault(cancer, []).append(f4["W_small"])
            fig4_D_all.append(f4["D"])
            fig4_sizes_all.extend([f4["size_large"], f4["size_small"]])
            fig4_Ws_all.extend([f4["W_large"], f4["W_small"]])
            fig4_patient_sizes_by_cancer.setdefault(cancer, {})[patient] = {
                "normal": f4["size_normal"], "large": f4["size_large"], "small": f4["size_small"]}

        log.info("  %s (%s): fig1=%s fig2=%s fig3=%d/%d fig4=%s", patient, cancer,
                 "ok" if res["fig1"] else "skip",
                 "ok" if (f2["clone_jacc"] is not None or f2["norm_jacc"] is not None) else "skip",
                 len(f3["intra"]), len(f3["inter"]),
                 "ok" if f4 is not None else "skip")

    log.info("%d/%d patient(s) contributed to at least one figure", n_ok, len(files))

    stats_summary = {}

    # ---- Figure 1: bootstrap support, unique vs shared ----
    d1 = os.path.join(args.out, "fig1_bootstrap_support")
    all_unique = [v for vs in fig1_unique_by_cancer.values() for v in vs]
    all_shared = [v for vs in fig1_shared_by_cancer.values() for v in vs]
    plot_support_boxplot(all_unique, all_shared, "Bootstrap support: unique vs shared edges (all patients)",
                         os.path.join(d1, "support_all_patients.png"))
    for cancer in fig1_unique_by_cancer:
        plot_support_boxplot(fig1_unique_by_cancer[cancer], fig1_shared_by_cancer.get(cancer, []),
                             f"Bootstrap support: unique vs shared edges ({cancer})",
                             os.path.join(d1, f"support_{cancer}.png"))
    plot_per_patient_boxplot_by_cancer(
        fig1_wdist_by_cancer, "Wasserstein distance",
        "Per-patient distance between unique/shared support distributions",
        os.path.join(d1, "wasserstein_by_cancer_type.png"))
    plot_per_patient_boxplot_by_cancer(
        fig1_pval_by_cancer, "MWU p-value (unique<shared)",
        "Per-patient significance of unique vs shared support difference",
        os.path.join(d1, "mwu_pvalue_by_cancer_type.png"), log_scale=True)
    plot_grouped_boxplot_by_cancer(
        fig1_unique_by_cancer, fig1_shared_by_cancer, "unique", "shared",
        "Bootstrap support (/100)", "Bootstrap support: unique vs shared edges, by cancer type",
        os.path.join(d1, "support_by_cancer_type.png"), palette="tab10")
    if all_unique and all_shared:
        _, p_all = mannwhitneyu(all_unique, all_shared, alternative="less")
        stats_summary["fig1_pooled_mwu_p_unique_lt_shared"] = float(p_all)

    # ---- Figure 2: clone-vs-clone vs normal-vs-clone Jaccard ----
    d2 = os.path.join(args.out, "fig2_clone_normal_jaccard")
    all_clone_j = [v for vs in fig2_clone_by_cancer.values() for v in vs]
    all_norm_j = [v for vs in fig2_norm_by_cancer.values() for v in vs]
    plot_clone_normal_jaccard_boxplot(all_clone_j, all_norm_j,
                                      "Pairwise Jaccard: clone-clone vs normal-clone (all patients)",
                                      os.path.join(d2, "jaccard_all_patients.png"))
    for cancer in set(list(fig2_clone_by_cancer) + list(fig2_norm_by_cancer)):
        plot_clone_normal_jaccard_boxplot(fig2_clone_by_cancer.get(cancer, []),
                                          fig2_norm_by_cancer.get(cancer, []),
                                          f"Pairwise Jaccard: clone-clone vs normal-clone ({cancer})",
                                          os.path.join(d2, f"jaccard_{cancer}.png"))
    plot_grouped_boxplot_by_cancer(
        fig2_clone_by_cancer, fig2_norm_by_cancer, "clone v clone", "normal v clone",
        "Pairwise Jaccard similarity", "Clone-clone vs normal-clone Jaccard, by cancer type",
        os.path.join(d2, "jaccard_by_cancer_type.png"), palette=sns.color_palette("tab10")[4:6])
    if all_clone_j and all_norm_j:
        _, p_all = mannwhitneyu(all_clone_j, all_norm_j, alternative="greater")
        stats_summary["fig2_pooled_mwu_p_clone_gt_normal"] = float(p_all)

    # ---- Figure 3: intra vs inter Hamming distance (original GRN) ----
    d3 = os.path.join(args.out, "fig3_intra_inter_hamming")
    all_intra_dict = {k: v for pd_ in fig3_intra_dict_by_cancer.values() for k, v in pd_.items()}
    all_inter_dict = {k: v for pd_ in fig3_inter_dict_by_cancer.values() for k, v in pd_.items()}
    if all_intra_dict:
        p_all = plot_split_violin_with_stars(
            all_intra_dict, all_inter_dict,
            "Intra (clone-clone) vs Inter (clone-Normal) Hamming, original GRN (all patients)",
            os.path.join(d3, "hamming_all_patients.png"))
        if p_all is not None:
            stats_summary["fig3_pooled_mwu_p"] = p_all
    for cancer in fig3_intra_dict_by_cancer:
        plot_split_violin_with_stars(
            fig3_intra_dict_by_cancer[cancer], fig3_inter_dict_by_cancer.get(cancer, {}),
            f"Intra (clone-clone) vs Inter (clone-Normal) Hamming, original GRN ({cancer})",
            os.path.join(d3, f"hamming_{cancer}.png"))

    # ---- Figure 4: clone size vs GRN divergence ----
    d4 = os.path.join(args.out, "fig4_clone_size_divergence")
    all_large = [v for vs in fig4_large_by_cancer.values() for v in vs]
    all_small = [v for vs in fig4_small_by_cancer.values() for v in vs]
    plot_large_small_boxplot(all_large, all_small,
                             "Wasserstein divergence from Normal: large vs small clone (all patients)",
                             os.path.join(d4, "divergence_all_patients.png"))
    for cancer in set(list(fig4_large_by_cancer) + list(fig4_small_by_cancer)):
        plot_large_small_boxplot(fig4_large_by_cancer.get(cancer, []), fig4_small_by_cancer.get(cancer, []),
                                 f"Wasserstein divergence from Normal: large vs small clone ({cancer})",
                                 os.path.join(d4, f"divergence_{cancer}.png"))
    for cancer, patient_sizes in fig4_patient_sizes_by_cancer.items():
        plot_clone_size_bars(patient_sizes, cancer, d4)

    n_D = len(fig4_D_all)
    if n_D >= 1:
        T = sum(1 for d in fig4_D_all if d > 0)
        binom_result = binomtest(T, n=n_D, p=0.5, alternative="greater")
        stats_summary["fig4_n_patients"] = n_D
        stats_summary["fig4_binom_p_D_gt_0"] = float(binom_result.pvalue)
        if n_D >= 2 and any(d != 0 for d in fig4_D_all):
            wilcoxon_result = wilcoxon(fig4_D_all, alternative="greater")
            stats_summary["fig4_wilcoxon_p_D_gt_0"] = float(wilcoxon_result.pvalue)
        if len(fig4_sizes_all) >= 3:
            rho, spearman_p = spearmanr(fig4_sizes_all, fig4_Ws_all)
            stats_summary["fig4_spearman_rho_size_vs_W"] = float(rho)
            stats_summary["fig4_spearman_p_size_vs_W"] = float(spearman_p)
        log.info("    Figure 4 (all patients, n=%d): binom p(D>0)=%.3e", n_D, binom_result.pvalue)

    stats_path = os.path.join(args.out, "stats_summary.json")
    with open(stats_path, "w") as fh:
        json.dump(stats_summary, fh, indent=2)
    log.info("Wrote pooled statistical tests to %s", stats_path)

    log.info("Done. Outputs under %s/", args.out)


if __name__ == "__main__":
    main()
