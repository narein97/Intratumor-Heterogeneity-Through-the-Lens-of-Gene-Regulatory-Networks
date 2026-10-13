#!/usr/bin/env python3
"""Clone-vs-clone and normal-vs-clone Jaccard boxplots, pooled and per cancer type.

Usage:
    python plot_jaccard_clone_vs_normal.py --root clones_out --out plots/
"""

import argparse
import glob
import os
import logging

import numpy as np
import pandas as pd
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)-7s %(message)s",
                    datefmt="%H:%M:%S")
log = logging.getLogger("jaccard_plot")


def cancer_type_from_path(h5ad_path, root):
    """clones_out/<Cancer>/<Dataset>/<Patient>_clones.h5ad -> <Cancer>"""
    rel = os.path.relpath(h5ad_path, root)
    return rel.split(os.sep)[0]


def per_patient_jaccard(h5ad_path):
    """Returns (avg clone-vs-clone, avg normal-vs-clone) Jaccard for one patient, or None if there are fewer than 2 malignant clones."""
    adata = ad.read_h5ad(h5ad_path)
    if "celloracle_original_filtered" not in adata.uns:
        log.warning("  [skip] %s -- run compute_bootstrap_support.py first "
                    "(no celloracle_original_filtered)", h5ad_path)
        return None

    df = adata.uns["celloracle_original_filtered"]
    clones = sorted(df["clone"].unique(), key=lambda c: (c == "N", c))
    malignant = [c for c in clones if c != "N"]
    if len(malignant) < 2:
        return None

    edge_sets = {c: set(zip(df.loc[df["clone"] == c, "source"],
                            df.loc[df["clone"] == c, "target"]))
                for c in clones}

    def jaccard(a, b):
        union = edge_sets[a] | edge_sets[b]
        return len(edge_sets[a] & edge_sets[b]) / len(union) if union else np.nan

    clone_jaccards = [jaccard(a, b) for i, a in enumerate(malignant)
                      for b in malignant[i + 1:]]
    normal_jaccards = ([jaccard(c, "N") for c in malignant]
                       if "N" in edge_sets else [])

    avg_clone = float(np.nanmean(clone_jaccards)) if clone_jaccards else np.nan
    avg_normal = float(np.nanmean(normal_jaccards)) if normal_jaccards else np.nan
    return avg_clone, avg_normal


def collect_all_patients(root):
    files = sorted(glob.glob(f"{root}/**/*_clones.h5ad", recursive=True))
    log.info("Found %d patient h5ad(s) under %s", len(files), root)

    rows = []
    for f in files:
        result = per_patient_jaccard(f)
        if result is None:
            continue
        avg_clone, avg_normal = result
        cancer = cancer_type_from_path(f, root)
        patient = os.path.basename(f).replace("_clones.h5ad", "")
        rows.append({"patient": patient, "cancer_type": cancer,
                    "clone v clone": avg_clone, "normal v clone": avg_normal})

    df = pd.DataFrame(rows)
    log.info("Computed Jaccard values for %d/%d patient(s) (rest skipped: "
            "<2 malignant clones or missing filtered edgelist)", len(df), len(files))
    return df


def plot_all_patients(df, out_dir):
    plot_df = df.melt(id_vars=["patient", "cancer_type"],
                      value_vars=["clone v clone", "normal v clone"],
                      var_name="Group", value_name="Jaccard similarity")
    plt.figure(figsize=(5, 6))
    ax = sns.boxplot(data=plot_df, x="Group", y="Jaccard similarity", palette="tab10")
    ax.tick_params(left=True)
    plt.title(f"Pairwise Jaccard similarity, clone vs. normal\n(all {df['patient'].nunique()} patients)")
    plt.tight_layout()
    out_path = os.path.join(out_dir, "jaccard_clone_vs_normal_ALL.png")
    plt.savefig(out_path, dpi=200)
    plt.close()
    log.info("-> %s", out_path)


def plot_by_cancer_type(df, out_dir):
    plot_df = df.melt(id_vars=["patient", "cancer_type"],
                      value_vars=["clone v clone", "normal v clone"],
                      var_name="Group", value_name="Jaccard similarity")
    cancer_types = sorted(df["cancer_type"].unique())
    n = len(cancer_types)
    ncol = min(5, n)
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.2 * ncol, 4 * nrow), squeeze=False)
    for i, ct in enumerate(cancer_types):
        ax = axes[i // ncol][i % ncol]
        sub = plot_df[plot_df["cancer_type"] == ct]
        sns.boxplot(data=sub, x="Group", y="Jaccard similarity", ax=ax, palette="tab10")
        ax.set_title(f"{ct} (n={df[df['cancer_type']==ct]['patient'].nunique()})", fontsize=10)
        ax.tick_params(axis="x", rotation=20)
    for i in range(n, nrow * ncol):
        axes[i // ncol][i % ncol].axis("off")
    fig.suptitle("Pairwise Jaccard similarity, clone vs. normal, by cancer type")
    fig.tight_layout()
    out_path = os.path.join(out_dir, "jaccard_clone_vs_normal_BY_CANCER_TYPE.png")
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    log.info("-> %s", out_path)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--root", required=True)
    p.add_argument("--out", default=".")
    args = p.parse_args()
    os.makedirs(args.out, exist_ok=True)

    df = collect_all_patients(args.root)
    if df.empty:
        log.error("No usable patients found -- did you run compute_bootstrap_support.py first?")
        return

    plot_all_patients(df, args.out)
    plot_by_cancer_type(df, args.out)


if __name__ == "__main__":
    main()
