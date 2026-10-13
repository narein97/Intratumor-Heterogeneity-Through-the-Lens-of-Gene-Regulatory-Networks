#!/usr/bin/env python3
"""Pre- vs post-chemotherapy clone-pair GRN divergence and gene out-degree shifts in the EOC cohort.

Usage:
    python treatment_phase_comparison_zhang2022.py --root clones_out/Ovarian --out treatment_phase_out --pre-label treatment-naive --post-label post-NACT --min-abs-cells 50
"""

import argparse
import glob
import itertools
import json
import os
import logging
import multiprocessing as mp

import numpy as np
import pandas as pd
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from scipy.stats import mannwhitneyu

from bootstrap_figures_efficient import (
    decode_and_filter,
    to_per_clone_edge_sets,
    to_per_clone_weighted_edges,
    compute_wasserstein,
    clone_wasserstein,
    build_weighted_graph,
    PLOT_DPI,
    _wrap_title,
)
from intra_vs_inter_efficient import edge_sets_to_matrix, pairwise_hamming, jaccard

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)-7s %(message)s",
                    datefmt="%H:%M:%S")
log = logging.getLogger("treatment_phase_comparison")


# =============================================================================
# Clone -> phase assignment
# =============================================================================

def assign_clone_phase(obs, clone_col, phase_col, pre_label, post_label, min_abs_cells):
    """
    Returns (phase_of_clone, purity_rows):
      phase_of_clone: {clone_label(str): 'pre'|'post'|'mixed'|'insufficient_cells'}
      purity_rows: one dict per non-N clone, for the audit manifest -- still
      reports n_pre/n_post/proportions for transparency even though the
      decision below is driven entirely by the absolute min_abs_cells
      floor, not by proportion (see the module docstring's point 1 for the
      full rule and rationale, confirmed with the user).

    n_other counts cells whose phase_col value is neither pre_label nor
    post_label (a typo, a third timepoint, missing data, ...).

    Returns ({}, []) if clone_col or phase_col isn't present in this
    patient's obs at all (the caller treats that as a skip, not a crash).
    """
    if clone_col not in obs.columns or phase_col not in obs.columns:
        return {}, []

    raw_phase = obs[phase_col].astype(str)
    df = pd.DataFrame({
        "clone": obs[clone_col].astype(str).to_numpy(),
        "phase": raw_phase.to_numpy(),
    })

    is_pre = df["phase"].str.contains(pre_label, case=False, na=False, regex=False)
    is_post = df["phase"].str.contains(post_label, case=False, na=False, regex=False)
    both = is_pre & is_post
    if both.any():
        log.warning("    %d cell(s) matched BOTH pre-label %r and post-label %r in obs['%s'] "
                    "-- counted as neither (ambiguous, not silently assigned)",
                    int(both.sum()), pre_label, post_label, phase_col)
    df["is_pre"] = is_pre & ~both
    df["is_post"] = is_post & ~both

    phase_of_clone = {}
    rows = []
    for clone, grp in df.groupby("clone"):
        if clone == "N":
            continue
        n_total = len(grp)
        n_pre = int(grp["is_pre"].sum())
        n_post = int(grp["is_post"].sum())
        n_other = n_total - n_pre - n_post
        prop_pre = n_pre / n_total if n_total else 0.0
        prop_post = n_post / n_total if n_total else 0.0

        pre_ok = n_pre >= min_abs_cells
        post_ok = n_post >= min_abs_cells
        if pre_ok and post_ok:
            assigned = "mixed"
        elif pre_ok:
            assigned = "pre"
        elif post_ok:
            assigned = "post"
        else:
            assigned = "insufficient_cells"

        phase_of_clone[clone] = assigned
        rows.append({"clone": clone, "n_cells": n_total, "n_pre": n_pre, "n_post": n_post,
                     "n_other": n_other, "proportion_pre": round(prop_pre, 4),
                     "proportion_post": round(prop_post, 4), "assigned_phase": assigned})

    if (df["is_pre"].sum() == 0) and (df["is_post"].sum() == 0):
        examples = df["phase"].unique()[:5].tolist()
        log.warning("    obs['%s'] matched neither pre-label %r nor post-label %r for ANY "
                    "cell in this patient -- example raw values seen: %s -- check --phase-col/"
                    "--pre-label/--post-label against these", phase_col, pre_label, post_label,
                    examples)

    return phase_of_clone, rows


# =============================================================================
# Data loading -- one h5ad read, decoded once
# =============================================================================

def load_patient_for_phase_comparison(h5ad_path, p_thresh, clone_col, phase_col,
                                      pre_label, post_label, min_abs_cells):
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

    phase_of_clone, purity_rows = assign_clone_phase(
        adata.obs, clone_col, phase_col, pre_label, post_label, min_abs_cells)

    return {
        "clones": all_clones,
        "original_sets": original_sets,
        "original_weighted": original_weighted,
        "phase_of_clone": phase_of_clone,
        "purity_rows": purity_rows,
    }


# =============================================================================
# Within-phase vs cross-phase clone-pair metrics
# =============================================================================

def compute_phase_pair_metrics(original_sets, original_weighted, phase_of_clone, all_clones):
    """
    For every pair of non-N clones that both received a 'pre' or 'post'
    assignment ('mixed' and 'insufficient_cells' clones are dropped
    entirely -- assigning them to either phase would be misleading, which
    is the whole point of the min-abs-cells floor), classify the pair as
    'within' (both pre, or both post) or 'cross' (one pre, one post), and
    compute Hamming distance, Jaccard similarity, and Wasserstein
    divergence for that pair.
    """
    usable = [c for c in all_clones if c != "N" and phase_of_clone.get(c) in ("pre", "post")]
    result = {"hamming": {"within": [], "cross": []},
              "jaccard": {"within": [], "cross": []},
              "wasserstein": {"within": [], "cross": []}}
    if len(usable) < 2:
        return result

    all_edges = set()
    for c in usable:
        all_edges |= original_sets[c]
    edge_index = {e: i for i, e in enumerate(sorted(all_edges))}
    n_edges = len(edge_index)
    vectors = {c: edge_sets_to_matrix([original_sets[c]], edge_index) for c in usable}

    for c1, c2 in itertools.combinations(usable, 2):
        group = "within" if phase_of_clone[c1] == phase_of_clone[c2] else "cross"

        if n_edges > 0:
            h = pairwise_hamming(vectors[c1], vectors[c2], n_edges)
            result["hamming"][group].append(float(h[0, 0]))

        j = jaccard(original_sets[c1], original_sets[c2])
        if not np.isnan(j):
            result["jaccard"][group].append(j)

        G1 = build_weighted_graph(original_weighted[c1])
        G2 = build_weighted_graph(original_weighted[c2])
        result["wasserstein"][group].append(float(clone_wasserstein(G1, G2)))

    return result


# =============================================================================
# Gene-level regulatory activity shift (degree centrality), pre vs post
# =============================================================================

def build_phase_network(original_weighted, clones_in_phase):
    """Combine one or more clones' weighted edge dicts ({(source,target):
    coef_mean}) into a single directed graph representing that phase. An
    edge present in more than one contributing clone gets its coef_mean
    AVERAGED across those clones (not summed), so a phase with more
    surviving clones doesn't mechanically inflate edge WEIGHT relative to
    a phase with fewer. Edge PRESENCE, however, is a plain union across
    the contributing clones -- and a union's size grows with the number
    of sets unioned, so a phase combining more clones will generally end
    up with a bigger network (more distinct edges) than one combining
    fewer, purely as a set-theoretic effect, independent of any real
    biological difference. That means raw degree (a count of edges
    touching a gene) IS confounded by how many clones fed each phase's
    network -- see compute_gene_degree_shift's normalized *_share_pct
    columns, which divide by that phase's total edge count specifically
    to correct for this. Only used for the gene-centrality-shift analysis
    below -- the pairwise divergence metrics above stay strictly clone vs
    clone, never pooled like this."""
    combined_sum, combined_n = {}, {}
    for c in clones_in_phase:
        for e, w in original_weighted[c].items():
            combined_sum[e] = combined_sum.get(e, 0.0) + w
            combined_n[e] = combined_n.get(e, 0) + 1
    averaged = {e: combined_sum[e] / combined_n[e] for e in combined_sum}
    return build_weighted_graph(averaged)


def compute_gene_degree_shift(G_pre, G_post):
    """Per-gene out-degree ("master regulator" activity -- how many genes
    it regulates) and in-degree (how many regulators control it) in the
    combined pre- and post-treatment phase networks, plus the pre->post
    delta. A gene present in only one phase's network gets degree 0 in
    the other (not dropped) -- a gene that's a hub pre-treatment and
    completely disconnected post-treatment (or vice versa) is exactly the
    kind of finding this is meant to surface, not something to filter
    out before the user sees it.

    Raw degree alone is confounded by network size (see
    build_phase_network's docstring: a phase combining more clones tends
    to have a bigger unioned network regardless of biology), so every
    gene also gets its degree expressed as a SHARE of that phase
    network's total edge count (out_degree_pre_share_pct = 100 *
    out_degree_pre / total_edges_pre, etc.) and the delta of that share
    (delta_out_degree_share_pct = post_share - pre_share). This delta is
    well-defined even at a zero baseline, unlike a ratio, which is what
    makes it the right statistic for genes that are entirely absent from
    one phase's network -- those are additionally flagged via
    out_degree_emergence / in_degree_emergence: 'gained_denovo' (0 in
    pre, >0 in post), 'lost_denovo' (>0 in pre, 0 in post), or 'none'.
    total_edges_pre/post are also returned per-row (same value repeated
    on every row) purely for downstream auditing -- so a consumer of the
    per-patient CSV can see, without re-deriving it, how big a network
    each share_pct was computed against."""
    total_pre = G_pre.number_of_edges()
    total_post = G_post.number_of_edges()
    all_genes = sorted(set(G_pre.nodes()) | set(G_post.nodes()))
    rows = []
    for g in all_genes:
        out_pre = G_pre.out_degree(g) if G_pre.has_node(g) else 0
        in_pre = G_pre.in_degree(g) if G_pre.has_node(g) else 0
        out_post = G_post.out_degree(g) if G_post.has_node(g) else 0
        in_post = G_post.in_degree(g) if G_post.has_node(g) else 0

        out_pre_share = 100.0 * out_pre / total_pre if total_pre else 0.0
        out_post_share = 100.0 * out_post / total_post if total_post else 0.0
        in_pre_share = 100.0 * in_pre / total_pre if total_pre else 0.0
        in_post_share = 100.0 * in_post / total_post if total_post else 0.0

        def _emergence(pre_val, post_val):
            if pre_val == 0 and post_val > 0:
                return "gained_denovo"
            if pre_val > 0 and post_val == 0:
                return "lost_denovo"
            return "none"

        rows.append({
            "gene": g,
            "out_degree_pre": out_pre, "out_degree_post": out_post,
            "delta_out_degree": out_post - out_pre,
            "out_degree_pre_share_pct": round(out_pre_share, 6),
            "out_degree_post_share_pct": round(out_post_share, 6),
            "delta_out_degree_share_pct": round(out_post_share - out_pre_share, 6),
            "out_degree_emergence": _emergence(out_pre, out_post),
            "in_degree_pre": in_pre, "in_degree_post": in_post,
            "delta_in_degree": in_post - in_pre,
            "in_degree_pre_share_pct": round(in_pre_share, 6),
            "in_degree_post_share_pct": round(in_post_share, 6),
            "delta_in_degree_share_pct": round(in_post_share - in_pre_share, 6),
            "in_degree_emergence": _emergence(in_pre, in_post),
            "total_edges_pre": total_pre, "total_edges_post": total_post,
        })
    return rows


# =============================================================================
# Plotting
# =============================================================================

def plot_within_cross_boxplot(within_vals, cross_vals, ylabel, title, out_path, palette=None):
    if not within_vals and not cross_vals:
        log.warning("    [skip plot] %s -- no data (0 within-phase and 0 cross-phase pairs "
                    "across all patients -- check clone_phase_purity_manifest.csv: likely "
                    "every patient has <2 usable, non-mixed clones)", title)
        return
    rows = ([{"value": v, "group": "within-phase"} for v in within_vals] +
            [{"value": v, "group": "cross-phase"} for v in cross_vals])
    df = pd.DataFrame(rows)
    pval = np.nan
    if within_vals and cross_vals:
        _, pval = mannwhitneyu(within_vals, cross_vals)

    n_groups = df["group"].nunique()
    full_palette = list(palette or sns.color_palette("tab10")[:2])

    plt.figure(figsize=(6, 6.5))
    ax = sns.boxplot(data=df, x="group", y="value", hue="group", legend=False,
                     palette=full_palette[:n_groups], showfliers=False)
    ax.tick_params(left=True)
    subtitle = f"MWU p={pval:.2e}" if not np.isnan(pval) else "MWU: n/a"
    ax.set_title(f"{_wrap_title(title)}\n{subtitle}")
    ax.set_xlabel("")
    ax.set_ylabel(ylabel)
    plt.tight_layout()
    plt.savefig(out_path, dpi=PLOT_DPI)
    plt.close()
    log.info("    %s: n_within=%d n_cross=%d MWU p=%s", title, len(within_vals), len(cross_vals),
             f"{pval:.3e}" if not np.isnan(pval) else "n/a")


def plot_split_violin_phase(within_dict, cross_dict, title, out_path):
    """
    Same star-overlay design as bootstrap_figures_efficient.py's
    plot_split_violin_with_stars (one median-star marker per patient over
    a split violin), reimplemented here under its own labels
    ('within-phase' / 'cross-phase', purple/brown) rather than imported --
    reusing a same-shaped plot function with the wrong axis/legend labels
    was exactly the bug caught earlier in this pipeline, not repeating
    that here.

    within_dict / cross_dict: {patient: [values]}.
    Returns the pooled Mann-Whitney U p-value, or None if there wasn't
    enough data to compute one.
    """
    star_colors = ["#AA5042", "#D5A021", "#1098F7", "#061A40", "#631A86",
                   "#226F54", "#CE4865", "#C6AD94", "#E3D7FF", "#7A9E7E", "#B85C38"]

    rows = []
    for patient, vals in within_dict.items():
        for v in vals:
            rows.append({"distance": v, "type": "within", "pos": 0})
    for patient, vals in cross_dict.items():
        for v in vals:
            rows.append({"distance": v, "type": "cross", "pos": 0})
    df = pd.DataFrame(rows)
    if df.empty:
        log.warning("    [skip plot] %s -- no data", title)
        return None

    group1 = df.loc[df["type"] == "within", "distance"]
    group2 = df.loc[df["type"] == "cross", "distance"]
    pval = None
    if len(group1) > 0 and len(group2) > 0:
        _, pval = mannwhitneyu(group1, group2)

    plt.figure(figsize=(8, 7))
    vp = sns.violinplot(x="pos", y="distance", hue="type", data=df, split=True,
                        inner="quart", fill=False,
                        palette={"within": "#9372B2", "cross": "#845B53"}, legend=False)
    vp.tick_params(left=True)

    all_patients = sorted(set(within_dict) | set(cross_dict))
    for i, patient in enumerate(all_patients):
        color = star_colors[i % len(star_colors)]
        temp_rows = []
        if within_dict.get(patient):
            temp_rows.append({"pos": 0, "distance": float(np.mean(within_dict[patient])), "type": "within"})
        if cross_dict.get(patient):
            temp_rows.append({"pos": 0, "distance": float(np.mean(cross_dict[patient])), "type": "cross"})
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

    legend_handles = [mpatches.Patch(color="#9372B2", label="Within-phase (pre-pre or post-post)"),
                      mpatches.Patch(color="#845B53", label="Cross-phase (pre-post)")]
    plt.legend(handles=legend_handles, title=f"n={len(all_patients)} patients\n(stars = per-patient medians)",
              bbox_to_anchor=(1.05, 1), loc="upper left", borderaxespad=0.)

    subtitle = f"MWU p={pval:.2e}" if pval is not None else "MWU: n/a"
    plt.title(f"{_wrap_title(title, width=45)}\n{subtitle}")
    plt.ylabel("Hamming distance")
    plt.xlabel("")
    plt.tight_layout()
    plt.savefig(out_path, dpi=PLOT_DPI)
    plt.close()
    log.info("    %s: n_patients=%d MWU p=%s", title, len(all_patients),
             f"{pval:.3e}" if pval is not None else "n/a")
    return float(pval) if pval is not None else None


def plot_gene_degree_shift(cohort_df, out_path, top_n=20, min_patients=2):
    """Horizontal diverging bar chart of the genes with the largest, most
    consistent regulatory out-degree shift between phases, restricted to
    genes seen (i.e. present in at least one phase's network) in
    >=min_patients patients -- a shift seen in only one patient is still
    in the per-patient CSVs, just not reliable enough to headline here.

    Ranked by mean_delta_out_degree_share_pct (each gene's out-degree as
    a SHARE of that phase network's total edges, post minus pre), not raw
    mean_delta_out_degree. Raw degree is confounded by network size: a
    phase that happens to combine more clones into its union network
    (build_phase_network) will show inflated degree for most genes at
    once, independent of real biology -- see that function's and
    compute_gene_degree_shift's docstrings. The share-of-network delta is
    also well-defined for genes absent from one phase (a 0% baseline),
    so it doubles as the right ranking for de novo emergence genes,
    without needing a separate ratio-based special case."""
    sub = cohort_df[cohort_df["n_patients"] >= min_patients].copy()
    if sub.empty:
        log.warning("    [skip plot] gene degree shift -- no gene appears in >=%d patient(s) "
                    "(cohort too small/sparse for a recurring-hub-gene view; see the full "
                    "per-gene, per-patient CSVs under gene_degree_shift/ instead)", min_patients)
        return
    sub = sub.reindex(sub["mean_delta_out_degree_share_pct"].abs().sort_values(ascending=False).index)
    sub = sub.head(top_n).iloc[::-1]  # smallest magnitude at the bottom of a horizontal chart
    colors = ["#C0392B" if v > 0 else "#2E5EAA" for v in sub["mean_delta_out_degree_share_pct"]]

    fig_h = max(4, 0.35 * len(sub))
    plt.figure(figsize=(8, fig_h))
    plt.barh(sub["gene"], sub["mean_delta_out_degree_share_pct"], color=colors)
    plt.axvline(0, color="black", linewidth=0.8)
    plt.xlabel("Mean delta out-degree, as % share of that phase's total network edges "
               "(post - pre), across patients")
    plt.title(_wrap_title(
        f"Genes with the most consistent regulatory out-degree shift, post- vs "
        f"pre-treatment, size-normalized (n>={min_patients} patient(s) each)", width=55))
    legend_handles = [mpatches.Patch(color="#C0392B", label="Higher out-degree share POST-treatment"),
                      mpatches.Patch(color="#2E5EAA", label="Higher out-degree share PRE-treatment")]
    plt.legend(handles=legend_handles, loc="lower right", fontsize=8)
    plt.tight_layout()
    plt.savefig(out_path, dpi=PLOT_DPI)
    plt.close()
    log.info("    Gene degree-shift plot: %d gene(s) shown (of %d with n_patients>=%d), "
             "ranked by normalized share delta", len(sub),
             len(cohort_df[cohort_df["n_patients"] >= min_patients]), min_patients)


# =============================================================================
# Per-patient worker (runs in a separate process)
# =============================================================================

def process_patient(task):
    """Module-level so it's picklable for multiprocessing.Pool. Does only
    load+compute; every plt.* call and every cohort-level statistical test
    stays in the parent process (see main())."""
    (f, p_thresh, clone_col, phase_col, pre_label, post_label, min_abs_cells) = task
    patient = os.path.basename(f).replace("_clones.h5ad", "")

    loaded = load_patient_for_phase_comparison(
        f, p_thresh, clone_col, phase_col, pre_label, post_label, min_abs_cells)
    if loaded is None:
        return {"file": f, "patient": patient, "skipped": True,
                "reason": "missing celloracle_original/_columns/_genes/clone_categories"}
    if not loaded["phase_of_clone"]:
        return {"file": f, "patient": patient, "skipped": True,
                "reason": f"missing obs['{clone_col}'] or obs['{phase_col}']"}

    metrics = compute_phase_pair_metrics(
        loaded["original_sets"], loaded["original_weighted"],
        loaded["phase_of_clone"], loaded["clones"])

    n_pre = sum(1 for v in loaded["phase_of_clone"].values() if v == "pre")
    n_post = sum(1 for v in loaded["phase_of_clone"].values() if v == "post")
    n_mixed = sum(1 for v in loaded["phase_of_clone"].values() if v == "mixed")
    n_insufficient = sum(1 for v in loaded["phase_of_clone"].values() if v == "insufficient_cells")

    gene_shift_rows = None
    pre_clones = [c for c, v in loaded["phase_of_clone"].items() if v == "pre"]
    post_clones = [c for c, v in loaded["phase_of_clone"].items() if v == "post"]
    if pre_clones and post_clones:
        G_pre = build_phase_network(loaded["original_weighted"], pre_clones)
        G_post = build_phase_network(loaded["original_weighted"], post_clones)
        gene_shift_rows = compute_gene_degree_shift(G_pre, G_post)

    return {"file": f, "patient": patient, "skipped": False, "metrics": metrics,
            "purity_rows": loaded["purity_rows"], "n_pre_clones": n_pre,
            "n_post_clones": n_post, "n_mixed_clones": n_mixed,
            "n_insufficient_clones": n_insufficient,
            "gene_shift_rows": gene_shift_rows}


# =============================================================================
# main
# =============================================================================

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--root", required=True,
                   help="directory of the 11 Zhang2022 patients' *_clones.h5ad files")
    p.add_argument("--out", default="treatment_phase_out")
    p.add_argument("--p-thresh", type=float, default=0.05)
    p.add_argument("--clone-col", default="clone",
                   help="obs column identifying each cell's clone -- must match whichever "
                        "column was used when the CellOracle GRNs stored in .uns were "
                        "inferred (e.g. 'clone' for the inferCNV-original clones, "
                        "'clone_leiden' for the Leiden-refined ones -- see "
                        "add_bootstrap_shuffle_zhang2022.py's own --clone-col)")
    p.add_argument("--phase-col", default="treatment_phase",
                   help="obs column to read the treatment phase from. Confirmed against the "
                        "real data: obs['treatment_phase'] already exists and holds clean "
                        "values 'treatment-naive' / 'post-NACT' (obs['sample'] also carries "
                        "the phase, e.g. 'EOC1005_primary_Peritoneum', but treatment_phase is "
                        "the intended, already-clean column). --pre-label/--post-label are "
                        "matched as a case-insensitive SUBSTRING of this column (not exact "
                        "equality) purely for robustness to whitespace/casing -- there's no "
                        "risk of 'treatment-naive' and 'post-NACT' colliding as substrings.")
    p.add_argument("--pre-label", default="treatment-naive",
                   help="substring of obs[--phase-col] marking treatment-naive cells")
    p.add_argument("--post-label", default="post-NACT",
                   help="substring of obs[--phase-col] marking post-NACT cells")
    p.add_argument("--min-abs-cells", type=int, default=50,
                   help="absolute cell-count floor used to assign each clone to a phase: "
                        "if both phases reach this many cells the clone is 'mixed' "
                        "(excluded); if only one phase reaches it the clone is assigned "
                        "to that phase regardless of purity; if neither phase reaches it "
                        "the clone is 'insufficient_cells' (excluded, too small to trust "
                        "a GRN comparison from). See module docstring for full rationale.")
    p.add_argument("--gene-shift-min-patients", type=int, default=2,
                   help="the top-gene-shift plot only shows genes seen (present in at least "
                        "one phase's network) in >= this many patients; the full per-gene "
                        "table (gene_degree_shift/cohort_gene_degree_shift_summary.csv) "
                        "always includes every gene regardless of this cutoff")
    p.add_argument("--gene-shift-top-n", type=int, default=20,
                   help="number of genes shown in the top-gene-shift bar chart")
    p.add_argument("--n-workers", type=int, default=max(1, (os.cpu_count() or 2) - 1),
                   help="Multiprocessing workers (no SLURM -- plain local "
                        "multiprocessing.Pool). Set to 1 to force sequential.")
    args = p.parse_args()
    os.makedirs(args.out, exist_ok=True)

    files = sorted(glob.glob(f"{args.root}/**/*_clones.h5ad", recursive=True))
    files = [f for f in files if "_with_support" not in f]
    log.info("Processing %d patient(s) with %d worker process(es)", len(files), args.n_workers)
    log.info("clone-col=%s phase-col=%s pre-label=%s post-label=%s min-abs-cells=%d",
             args.clone_col, args.phase_col, args.pre_label, args.post_label, args.min_abs_cells)

    tasks = [(f, args.p_thresh, args.clone_col, args.phase_col,
             args.pre_label, args.post_label, args.min_abs_cells) for f in files]
    if args.n_workers > 1 and len(tasks) > 1:
        with mp.Pool(processes=min(args.n_workers, len(tasks))) as pool:
            # order-preserving map -- see process_patient()'s docstring and
            # intra_vs_inter_efficient.py's main() for why this matters
            # (byte-identical parallel-vs-serial output).
            results = pool.map(process_patient, tasks)
    else:
        results = [process_patient(t) for t in tasks]

    hamming_within_dict, hamming_cross_dict = {}, {}
    all_jacc_within, all_jacc_cross = [], []
    all_wass_within, all_wass_cross = [], []
    purity_manifest_rows = []
    cohort_gene_out_deltas, cohort_gene_in_deltas = {}, {}
    cohort_gene_out_share_deltas, cohort_gene_in_share_deltas = {}, {}
    denovo_emergence_events = []
    patients_missing_a_phase = []
    n_ok = 0
    n_gene_shift = 0

    gene_shift_dir = os.path.join(args.out, "gene_degree_shift")
    os.makedirs(gene_shift_dir, exist_ok=True)

    for res in results:
        patient = res["patient"]
        if res["skipped"]:
            log.warning("  [skip] %s -- %s", patient, res["reason"])
            continue
        n_ok += 1
        for row in res["purity_rows"]:
            purity_manifest_rows.append({"patient": patient, **row})

        m = res["metrics"]
        if m["hamming"]["within"] or m["hamming"]["cross"]:
            hamming_within_dict[patient] = m["hamming"]["within"]
            hamming_cross_dict[patient] = m["hamming"]["cross"]
        all_jacc_within.extend(m["jaccard"]["within"])
        all_jacc_cross.extend(m["jaccard"]["cross"])
        all_wass_within.extend(m["wasserstein"]["within"])
        all_wass_cross.extend(m["wasserstein"]["cross"])

        if res["gene_shift_rows"] is not None:
            n_gene_shift += 1
            pd.DataFrame(res["gene_shift_rows"]).to_csv(
                os.path.join(gene_shift_dir, f"{patient}_gene_degree_shift.csv"), index=False)
            for row in res["gene_shift_rows"]:
                cohort_gene_out_deltas.setdefault(row["gene"], []).append(row["delta_out_degree"])
                cohort_gene_in_deltas.setdefault(row["gene"], []).append(row["delta_in_degree"])
                cohort_gene_out_share_deltas.setdefault(row["gene"], []).append(
                    row["delta_out_degree_share_pct"])
                cohort_gene_in_share_deltas.setdefault(row["gene"], []).append(
                    row["delta_in_degree_share_pct"])
                if row["out_degree_emergence"] != "none":
                    denovo_emergence_events.append({
                        "patient": patient, "gene": row["gene"], "degree_type": "out",
                        "direction": row["out_degree_emergence"],
                        "degree_pre": row["out_degree_pre"], "degree_post": row["out_degree_post"],
                        "share_pct_pre": row["out_degree_pre_share_pct"],
                        "share_pct_post": row["out_degree_post_share_pct"]})
                if row["in_degree_emergence"] != "none":
                    denovo_emergence_events.append({
                        "patient": patient, "gene": row["gene"], "degree_type": "in",
                        "direction": row["in_degree_emergence"],
                        "degree_pre": row["in_degree_pre"], "degree_post": row["in_degree_post"],
                        "share_pct_pre": row["in_degree_pre_share_pct"],
                        "share_pct_post": row["in_degree_post_share_pct"]})
        else:
            patients_missing_a_phase.append(patient)

        log.info("  %s: pre-clones=%d post-clones=%d mixed(excluded)=%d "
                 "insufficient_cells(excluded)=%d -- "
                 "%d within-phase pair(s), %d cross-phase pair(s), gene-shift=%s",
                 patient, res["n_pre_clones"], res["n_post_clones"], res["n_mixed_clones"],
                 res["n_insufficient_clones"],
                 len(m["jaccard"]["within"]), len(m["jaccard"]["cross"]),
                 "ok" if res["gene_shift_rows"] is not None else "skip (missing a phase)")

    log.info("%d/%d patient(s) contributed to at least one comparison", n_ok, len(files))
    log.info("Gene centrality shift computed for %d/%d patient(s) (needs >=1 pre-clone AND "
             ">=1 post-clone surviving the min-abs-cells floor); missing one phase entirely: %s",
             n_gene_shift, n_ok, ", ".join(patients_missing_a_phase) if patients_missing_a_phase
             else "none")

    manifest_path = os.path.join(args.out, "clone_phase_purity_manifest.csv")
    pd.DataFrame(purity_manifest_rows).to_csv(manifest_path, index=False)
    n_pre_total = sum(1 for r in purity_manifest_rows if r["assigned_phase"] == "pre")
    n_post_total = sum(1 for r in purity_manifest_rows if r["assigned_phase"] == "post")
    n_mixed_total = sum(1 for r in purity_manifest_rows if r["assigned_phase"] == "mixed")
    n_insufficient_total = sum(1 for r in purity_manifest_rows
                                if r["assigned_phase"] == "insufficient_cells")
    log.info("Clone-phase assignment across all patients: pre=%d post=%d mixed(excluded)=%d "
             "insufficient_cells(excluded)=%d (min-abs-cells=%d)",
             n_pre_total, n_post_total, n_mixed_total, n_insufficient_total, args.min_abs_cells)
    log.info("Wrote clone phase-purity manifest to %s", manifest_path)

    stats_summary = {}

    # ---- Jaccard ----
    plot_within_cross_boxplot(
        all_jacc_within, all_jacc_cross, "Pairwise Jaccard similarity",
        "Original-GRN Jaccard: within-phase vs cross-phase clone pairs (all patients)",
        os.path.join(args.out, "jaccard_within_vs_cross.png"))
    if all_jacc_within and all_jacc_cross:
        _, p_all = mannwhitneyu(all_jacc_within, all_jacc_cross)
        stats_summary["jaccard_pooled_mwu_p"] = float(p_all)

    # ---- Wasserstein ----
    plot_within_cross_boxplot(
        all_wass_within, all_wass_cross, "Wasserstein divergence",
        "Original-GRN Wasserstein divergence: within-phase vs cross-phase clone pairs (all patients)",
        os.path.join(args.out, "wasserstein_within_vs_cross.png"),
        palette=sns.color_palette("husl", 8)[1:3])
    if all_wass_within and all_wass_cross:
        _, p_all = mannwhitneyu(all_wass_within, all_wass_cross)
        stats_summary["wasserstein_pooled_mwu_p"] = float(p_all)

    # ---- Hamming ----
    if hamming_within_dict:
        p_all = plot_split_violin_phase(
            hamming_within_dict, hamming_cross_dict,
            "Original-GRN Hamming distance: within-phase vs cross-phase clone pairs (all patients)",
            os.path.join(args.out, "hamming_within_vs_cross.png"))
        if p_all is not None:
            stats_summary["hamming_pooled_mwu_p"] = p_all

    cohort_gene_rows = []
    for gene, out_deltas in cohort_gene_out_deltas.items():
        in_deltas = cohort_gene_in_deltas[gene]
        out_share_deltas = cohort_gene_out_share_deltas[gene]
        in_share_deltas = cohort_gene_in_share_deltas[gene]
        cohort_gene_rows.append({
            "gene": gene, "n_patients": len(out_deltas),
            "mean_delta_out_degree_share_pct": round(float(np.mean(out_share_deltas)), 4),
            "n_patients_gained_out_degree_share": sum(1 for d in out_share_deltas if d > 0),
            "n_patients_lost_out_degree_share": sum(1 for d in out_share_deltas if d < 0),
            "mean_delta_in_degree_share_pct": round(float(np.mean(in_share_deltas)), 4),
            "n_patients_gained_in_degree_share": sum(1 for d in in_share_deltas if d > 0),
            "n_patients_lost_in_degree_share": sum(1 for d in in_share_deltas if d < 0),
            # raw (network-size-confounded) columns, kept for reference only:
            "mean_delta_out_degree": round(float(np.mean(out_deltas)), 4),
            "n_patients_gained_out_degree": sum(1 for d in out_deltas if d > 0),
            "n_patients_lost_out_degree": sum(1 for d in out_deltas if d < 0),
            "mean_delta_in_degree": round(float(np.mean(in_deltas)), 4),
            "n_patients_gained_in_degree": sum(1 for d in in_deltas if d > 0),
            "n_patients_lost_in_degree": sum(1 for d in in_deltas if d < 0),
        })
    if cohort_gene_rows:
        cohort_gene_df = pd.DataFrame(cohort_gene_rows)
        cohort_gene_df = cohort_gene_df.reindex(
            cohort_gene_df["mean_delta_out_degree_share_pct"].abs().sort_values(ascending=False).index)
        cohort_summary_path = os.path.join(gene_shift_dir, "cohort_gene_degree_shift_summary.csv")
        cohort_gene_df.to_csv(cohort_summary_path, index=False)
        log.info("Wrote cohort gene degree-shift summary (%d gene(s), ranked by normalized "
                 "share delta) to %s", len(cohort_gene_df), cohort_summary_path)
        plot_gene_degree_shift(
            cohort_gene_df, os.path.join(gene_shift_dir, "top_gene_degree_shift.png"),
            top_n=args.gene_shift_top_n, min_patients=args.gene_shift_min_patients)

        if denovo_emergence_events:
            emergence_path = os.path.join(gene_shift_dir, "denovo_emergence_events.csv")
            pd.DataFrame(denovo_emergence_events).to_csv(emergence_path, index=False)
            n_gained = sum(1 for e in denovo_emergence_events if e["direction"] == "gained_denovo")
            n_lost = sum(1 for e in denovo_emergence_events if e["direction"] == "lost_denovo")
            log.info("Wrote %d de novo emergence event(s) (%d gained_denovo, %d lost_denovo -- "
                     "a gene going from 0 edges in one phase to >0 in the other, or vice versa; "
                     "each row is a single patient/gene event, NOT pooled across patients -- see "
                     "module docstring) to %s", len(denovo_emergence_events), n_gained, n_lost,
                     emergence_path)
        else:
            log.info("No de novo emergence events (no gene went from 0 edges in one phase to "
                     ">0 in the other, in either direction, in any patient)")
    else:
        log.warning("No patient had both a pre-clone and a post-clone -- gene degree-shift "
                    "analysis produced nothing (see clone_phase_purity_manifest.csv)")

    stats_path = os.path.join(args.out, "stats_summary.json")
    with open(stats_path, "w") as fh:
        json.dump(stats_summary, fh, indent=2)
    log.info("Wrote pooled statistical tests to %s", stats_path)

    log.info("Done. Outputs under %s/", args.out)


if __name__ == "__main__":
    main()
