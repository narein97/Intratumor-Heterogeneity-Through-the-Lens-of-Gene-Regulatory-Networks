#!/usr/bin/env python3
"""Intra- vs inter-clone bootstrap Hamming distances and the label-shuffle null over the cohort.

Usage:
    python intra_vs_inter_efficient.py --root clones_out --out intra_vs_inter_out
"""

import argparse
import glob
import os
import logging
import multiprocessing as mp

import numpy as np
import pandas as pd
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu

PLOT_DPI = 300
NORMAL_CLONE_LABEL = "N"  # same convention used throughout this pipeline
                          # (bootstrap_figures_efficient.py, the EOC/treatment-
                          # phase scripts, etc.)

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)-7s %(message)s",
                    datefmt="%H:%M:%S")
log = logging.getLogger("intra_vs_inter")


# =============================================================================
# Data loading -- decode to per-clone, per-replicate edge SETS
# =============================================================================

def decode_and_filter(arr, columns, genes, clone_categories, p_thresh):
    df = pd.DataFrame(arr, columns=columns)
    df["clone"] = clone_categories[df["clone"].astype(int)]
    df["source"] = genes[df["source"].astype(int)]
    df["target"] = genes[df["target"].astype(int)]
    return df[df["p"] < p_thresh]


def to_per_clone_edge_sets(df, all_clones):
    """{clone: set of (source,target)}. Every clone in all_clones gets an
    entry, EMPTY if it has no significant edges in this replicate -- never
    omitted, to preserve index alignment across replicates."""
    out = {c: set() for c in all_clones}
    for clone, grp in df.groupby("clone"):
        if clone in out:
            out[clone] = set(zip(grp["source"], grp["target"]))
    return out


def load_patient_grn_sets(h5ad_path, p_thresh=0.05):
    """
    Returns None if required keys are missing, else:
        {"original": {clone: edge_set},
         "bootstrap": {clone: [edge_set_rep0, ..., edge_set_repN]},  # index-aligned across clones
         "shuffle":   {clone: [edge_set_rep0, ..., edge_set_repN]}}  # index-aligned across clones
    """
    adata = ad.read_h5ad(h5ad_path, backed="r")
    required = ["celloracle_original", "celloracle_bootstrap", "celloracle_shuffle",
                "celloracle_columns", "celloracle_genes", "clone_categories"]
    if any(k not in adata.uns for k in required):
        return None

    columns = list(adata.uns["celloracle_columns"])
    genes = np.asarray(adata.uns["celloracle_genes"]).astype(str)
    clone_categories = np.asarray(adata.uns["clone_categories"]).astype(str)
    all_clones = list(clone_categories)

    original_df = decode_and_filter(adata.uns["celloracle_original"], columns, genes, clone_categories, p_thresh)
    original = to_per_clone_edge_sets(original_df, all_clones)

    def load_replicates(uns_key):
        rep_keys_sorted = sorted(adata.uns[uns_key].keys(), key=int)
        per_clone_lists = {c: [] for c in all_clones}
        for rk in rep_keys_sorted:
            df = decode_and_filter(adata.uns[uns_key][rk], columns, genes, clone_categories, p_thresh)
            per_clone = to_per_clone_edge_sets(df, all_clones)
            for c in all_clones:
                per_clone_lists[c].append(per_clone[c])  # always appends, even if empty set
        return per_clone_lists

    bootstrap = load_replicates("celloracle_bootstrap")
    shuffle = load_replicates("celloracle_shuffle")

    return {"original": original, "bootstrap": bootstrap, "shuffle": shuffle, "clones": all_clones}


# =============================================================================
# Analysis 1: intra vs inter bootstrap Hamming distance
# =============================================================================

def edge_sets_to_matrix(list_of_sets, edge_index):
    """(n_replicates, n_edges) boolean matrix, using a SHARED edge_index
    (dict edge->column) so different clones' matrices share column meaning."""
    M = np.zeros((len(list_of_sets), len(edge_index)), dtype=np.int32)
    for i, s in enumerate(list_of_sets):
        for e in s:
            M[i, edge_index[e]] = 1
    return M


def pairwise_hamming(M1, M2, n_edges):
    """All (i,j) Hamming proportions between rows of M1 and rows of M2,
    via |u|+|v|-2*dot(u,v))/n_edges -- verified identical to
    scipy.spatial.distance.hamming on the equivalent binary vectors."""
    sizes1 = M1.sum(axis=1)
    sizes2 = M2.sum(axis=1)
    dot = M1 @ M2.T
    return (sizes1[:, None] + sizes2[None, :] - 2 * dot) / n_edges


def compute_intra_inter_for_patient(bootstrap_dict):
    """Intra- vs inter-clone bootstrap Hamming distances for one patient (normal excluded)."""
    clones = sorted(c for c in bootstrap_dict.keys() if c != NORMAL_CLONE_LABEL)
    if len(clones) <= 1:
        return [], []

    all_edges = set()
    for c in clones:
        for s in bootstrap_dict[c]:
            all_edges |= s
    if not all_edges:
        return [], []
    edge_index = {e: i for i, e in enumerate(sorted(all_edges))}
    n_edges = len(edge_index)

    matrices = {c: edge_sets_to_matrix(bootstrap_dict[c], edge_index) for c in clones}

    intra_all = []
    for c in clones:
        M = matrices[c]
        h = pairwise_hamming(M, M, n_edges)  # includes diagonal (i==j) -- matches original
        intra_all.extend(h.flatten().tolist())

    inter_all = []
    for i in range(len(clones)):
        for j in range(i + 1, len(clones)):
            h = pairwise_hamming(matrices[clones[i]], matrices[clones[j]], n_edges)
            inter_all.extend(h.flatten().tolist())

    return intra_all, inter_all


# =============================================================================
# Analysis 2: shuffle-based Jaccard null model (Version A)
# =============================================================================

def jaccard(s1, s2):
    union = s1 | s2
    if not union:
        return np.nan  # both empty -- undefined, matches a 0/0 in the original
    return len(s1 & s2) / len(union)


def compute_shuffle_jaccard_for_patient(loaded):
    """Jaccard of the clone GRNs vs the shuffled-label null, in standard deviations, for one patient."""
    clones = sorted(c for c in loaded["clones"] if c != NORMAL_CLONE_LABEL)
    results = {}
    for i in range(len(clones)):
        for j in range(i + 1, len(clones)):
            c1, c2 = clones[i], clones[j]
            original_j = jaccard(loaded["original"][c1], loaded["original"][c2])
            n = min(len(loaded["shuffle"][c1]), len(loaded["shuffle"][c2]))
            shuffle_j = [jaccard(loaded["shuffle"][c1][k], loaded["shuffle"][c2][k]) for k in range(n)]
            results[f"{c1}-{c2}"] = {"original": original_j, "shuffle": shuffle_j}
    return results



def cancer_type_from_path(path, root):
    rel = os.path.relpath(path, root)
    return rel.split(os.sep)[0]


def plot_split_violin(intra_all, inter_all, title, out_path):
    rows = ([{"distance": v, "type": "intra"} for v in intra_all] +
            [{"distance": v, "type": "inter"} for v in inter_all])
    df = pd.DataFrame(rows)
    df["x"] = "all"

    stat, pval = mannwhitneyu(intra_all, inter_all)
    log.info("    %s: Mann-Whitney U=%.1f p=%.3e (n_intra=%d, n_inter=%d)",
             title, stat, pval, len(intra_all), len(inter_all))

    plt.figure(figsize=(5, 6))
    ax = sns.violinplot(x="x", y="distance", hue="type", data=df, split=True,
                        inner="quartile", palette={"intra": "#029E73", "inter": "#D55E00"})
    ax.set_title(f"{title}\nMann-Whitney p={pval:.2e}")
    ax.set_xlabel("")
    plt.tight_layout()
    plt.savefig(out_path, dpi=PLOT_DPI)
    plt.close()


def plot_shuffle_results(all_shuffle_results, title, out_dir):
    """all_shuffle_results: {patient_label: {pair_key: {"original":val, "shuffle":[...]}}}"""
    os.makedirs(out_dir, exist_ok=True)
    palette = ["#2f4f4f", "#8b4513", "#228b22", "#000080", "#ff0000", "#ffd700",
              "#7fff00", "#00ffff", "#ff00ff", "#1e90ff", "#ffe4b5", "#ff69b4"]

    stddevs = []
    for patient, pair_results in all_shuffle_results.items():
        plt.figure(figsize=(7, 5))
        keys = list(pair_results.keys())
        for idx, (pair, data) in enumerate(pair_results.items()):
            shuffle_vals = [v for v in data["shuffle"] if not np.isnan(v)]
            if not shuffle_vals or np.isnan(data["original"]):
                continue
            color = palette[idx % len(palette)]
            sns.histplot(shuffle_vals, stat="proportion", color=color, alpha=0.6)
            plt.axvline(data["original"], linestyle="--", linewidth=2, color=color)
            std = np.std(shuffle_vals)
            if std > 0:
                stddevs.append(abs(np.mean(shuffle_vals) - data["original"]) / std)
        plt.title(str(patient))
        plt.xlabel("Jaccard Similarity")
        plt.ylabel("Frequency")
        plt.tight_layout()
        plt.savefig(os.path.join(out_dir, f"{patient.replace('/', '_')}_shuffle_hist.png"), dpi=PLOT_DPI)
        plt.close()

    return stddevs


def process_patient(task):
    """Per-patient worker unit, run in a separate process by main()'s Pool.

    Does ONLY the per-patient math (loading + intra/inter + shuffle Jaccard);
    all cross-patient aggregation and every plt.* call stay in the parent
    process, since matplotlib figures cannot cross a process boundary and
    aggregation order matters for reproducibility (see main()).

    Must be a module-level function (not a closure/nested def) so it can be
    pickled and sent to worker processes.
    """
    f, root, p_thresh = task
    patient = os.path.basename(f).replace("_clones.h5ad", "")
    cancer = cancer_type_from_path(f, root)
    loaded = load_patient_grn_sets(f, p_thresh)
    if loaded is None:
        return {"file": f, "patient": patient, "cancer": cancer, "skipped": True}

    intra, inter = compute_intra_inter_for_patient(loaded["bootstrap"])
    shuffle_results = compute_shuffle_jaccard_for_patient(loaded)
    return {"file": f, "patient": patient, "cancer": cancer, "skipped": False,
            "intra": intra, "inter": inter, "shuffle_results": shuffle_results}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--root", required=True)
    p.add_argument("--out", default="intra_vs_inter_out")
    p.add_argument("--p-thresh", type=float, default=0.05)
    p.add_argument("--exclude-cancer-type", nargs="*", default=["Temporal"],
                   help="Cancer-type directory name(s) to skip (path segment, e.g. Temporal).")
    p.add_argument("--n-workers", type=int, default=max(1, (os.cpu_count() or 2) - 1),
                   help="Local worker processes (1 = sequential).")
    args = p.parse_args()
    os.makedirs(args.out, exist_ok=True)

    files = sorted(glob.glob(f"{args.root}/**/*_clones.h5ad", recursive=True))
    files = [f for f in files if "_with_support" not in f]
    if args.exclude_cancer_type:
        files = [f for f in files if not any(f"/{ct}/" in f for ct in args.exclude_cancer_type)]
    log.info("Processing %d patient(s) with %d worker process(es)", len(files), args.n_workers)

    intra_by_cancer, inter_by_cancer = {}, {}
    shuffle_by_cancer = {}

    tasks = [(f, args.root, args.p_thresh) for f in files]
    if args.n_workers > 1 and len(tasks) > 1:
        with mp.Pool(processes=min(args.n_workers, len(tasks))) as pool:
            results = pool.map(process_patient, tasks)
    else:
        results = [process_patient(t) for t in tasks]

    for res in results:
        if res["skipped"]:
            log.warning("  [skip] %s -- missing required .uns keys", res["file"])
            continue
        patient, cancer = res["patient"], res["cancer"]
        intra, inter = res["intra"], res["inter"]
        shuffle_results = res["shuffle_results"]

        intra_by_cancer.setdefault(cancer, []).extend(intra)
        inter_by_cancer.setdefault(cancer, []).extend(inter)
        shuffle_by_cancer.setdefault(cancer, {})[patient] = shuffle_results

        log.info("  %s (%s): %d intra, %d inter values; %d clone-pair(s) for shuffle",
                 patient, cancer, len(intra), len(inter), len(shuffle_results))

    # ---- Plot 1: split violin -- all-patients, and per cancer type ----
    all_intra = [v for vs in intra_by_cancer.values() for v in vs]
    all_inter = [v for vs in inter_by_cancer.values() for v in vs]
    plot_split_violin(all_intra, all_inter, "Intra vs Inter Hamming Distance (all patients)",
                      os.path.join(args.out, "violin_all_patients.png"))
    for cancer in intra_by_cancer:
        plot_split_violin(intra_by_cancer[cancer], inter_by_cancer[cancer],
                          f"Intra vs Inter Hamming Distance ({cancer})",
                          os.path.join(args.out, f"violin_{cancer}.png"))

    # ---- Plot 2: shuffle histograms + stddev boxplot -- all-patients, and per cancer type ----
    all_shuffle_flat = {f"{c}/{p}": r for c, pr in shuffle_by_cancer.items() for p, r in pr.items()}
    stddevs_all = plot_shuffle_results(all_shuffle_flat, "all patients",
                                       os.path.join(args.out, "shuffle_all_patients"))

    stddevs_by_cancer = {}
    for cancer, patients in shuffle_by_cancer.items():
        stddevs_by_cancer[cancer] = plot_shuffle_results(patients, cancer,
                                                         os.path.join(args.out, f"shuffle_{cancer}"))

    n_cancer = max(len(stddevs_by_cancer), 1)
    # Widen the figure with the number of cancer types so rotated labels have
    # room to sit below their own box instead of colliding with neighbors.
    fig_w = max(7, 0.55 * n_cancer)
    plt.figure(figsize=(fig_w, 6))
    sns.boxplot(data=[v for v in stddevs_by_cancer.values()], palette="tab10")
    ax = plt.gca()
    ax.set_xticks(range(len(stddevs_by_cancer)))
    ax.set_xticklabels(list(stddevs_by_cancer.keys()), rotation=45, ha="right",
                       rotation_mode="anchor")
    plt.ylabel("No. of standard deviations")
    plt.xlabel("Cancer type")
    plt.title("Deviation of original Jaccard from shuffle null (by cancer type)")
    plt.tight_layout()
    plt.savefig(os.path.join(args.out, "stddev_boxplot_by_cancer_type.png"), dpi=PLOT_DPI)
    plt.close()

    log.info("Done. Outputs under %s/", args.out)


if __name__ == "__main__":
    main()
