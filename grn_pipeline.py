"""
GRN Pipeline
============
Orchestrates four stages for one or more datasets and methods:

  1. Bootstrap data generation   (create_bootstrap_data)
  2. Shuffle data generation      (create_shuffle_data)
  3. GRN inference                (run_grn_parallel)
  4. Bootstrap support scoring    (compute_bootstrap_support)

Usage
-----
    python grn_pipeline.py --config pipeline_config.yaml [--stages STAGE ...]

To run only specific stages:
    python grn_pipeline.py --config pipeline_config.yaml --stages bootstrap grn_inference
"""

import os
import sys
import pickle
import random
import argparse
import importlib.util
import multiprocessing as mp

import numpy as np
import pandas as pd
import scanpy as sc
import networkx as nx
import yaml
from tqdm import tqdm

mp.set_start_method("spawn", force=True)


# ─────────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────────

def parse_args():
    parser = argparse.ArgumentParser(description="End-to-end GRN inference pipeline.")
    parser.add_argument(
        "--config", type=str, required=True,
        help="Path to the YAML pipeline configuration file.",
    )
    parser.add_argument(
        "--stages", type=str, nargs="*",
        choices=["bootstrap", "shuffle", "grn_inference", "bootstrap_support"],
        default=None,
        help=(
            "Stages to run. Overrides the 'run' flags in the config. "
            "If omitted, the config 'run' flags are respected."
        ),
    )
    return parser.parse_args()


# ─────────────────────────────────────────────────────────────────────────────
# Config loading
# ─────────────────────────────────────────────────────────────────────────────

def load_config(path: str) -> dict:
    with open(path) as f:
        return yaml.safe_load(f)


def should_run(stage_name: str, config: dict, override_stages: list | None) -> bool:
    """Return True if this stage should execute."""
    if override_stages is not None:
        return stage_name in override_stages
    return config["stages"][stage_name].get("run", True)


# ─────────────────────────────────────────────────────────────────────────────
# Stage 1 — Bootstrap data generation
# ─────────────────────────────────────────────────────────────────────────────

def run_bootstrap_stage(datasets: dict, data_dir: str, cfg: dict):
    """Generate bootstrap replicates for every sample in every dataset."""
    n_bootstraps   = cfg.get("n_bootstraps", 100)
    cell_proportion = cfg.get("cell_proportion", 0.9)

    print("\n" + "=" * 60)
    print("STAGE 1 — Bootstrap data generation")
    print("=" * 60)

    for dataset, samples in datasets.items():
        print(f"\nDataset: {dataset}")
        boot_root = os.path.join(data_dir, dataset, f"bootstrap_{dataset}")
        os.makedirs(boot_root, exist_ok=True)

        for sample, h5_path in samples.items():
            print(f"  Sample: {sample}")

            if not os.path.exists(h5_path):
                print(f"  h5ad not found: {h5_path}. Skipping.")
                continue

            adata = sc.read_h5ad(h5_path)

            if "inferred_clone" not in adata.obs:
                print(f"  'inferred_clone' missing. Skipping.")
                continue

            counts       = adata.obs["inferred_clone"].value_counts()
            valid_clones = counts[counts > 0].index.tolist()

            if "N" not in valid_clones:
                print(f"  No clone 'N' found. Skipping.")
                continue

            print(f"  Valid clones: {valid_clones}")
            sample_out = os.path.join(boot_root, sample)
            os.makedirs(sample_out, exist_ok=True)

            for i in tqdm(range(1, n_bootstraps + 1), desc=f"  Bootstrapping {sample}"):
                selection = []
                for clone in valid_clones:
                    idx    = adata.obs[adata.obs["inferred_clone"] == clone].index
                    n      = len(idx)
                    sample_n = max(1, int(cell_proportion * n))
                    chosen = np.random.choice(idx, size=sample_n, replace=False)
                    selection += chosen.tolist()

                adata[selection].copy().write_h5ad(
                    os.path.join(sample_out, f"{sample}_bootstrap_{i}.h5ad")
                )

    print("\nBootstrap stage complete.")


# ─────────────────────────────────────────────────────────────────────────────
# Stage 2 — Shuffle data generation
# ─────────────────────────────────────────────────────────────────────────────

def run_shuffle_stage(datasets: dict, data_dir: str, cfg: dict):
    """Generate shuffle replicates for every sample in every dataset."""
    n_shuffles  = cfg.get("n_shuffles", 100)
    random_seed = cfg.get("random_seed", None)

    print("\n" + "=" * 60)
    print("STAGE 2 — Shuffle data generation")
    print("=" * 60)

    for dataset, samples in datasets.items():
        print(f"\nDataset: {dataset}")
        shuffle_root = os.path.join(data_dir, dataset, f"shuffle_{dataset}")
        os.makedirs(shuffle_root, exist_ok=True)

        for sample in samples:
            print(f"  Sample: {sample}")

            # Shuffles are derived from the whole-sample bootstrap file
            h5_path = os.path.join(
                data_dir, dataset, f"bootstrap_{dataset}", sample, f"{sample}_whole.h5ad"
            )
            if not os.path.exists(h5_path):
                print(f"  Whole h5ad not found: {h5_path}. Skipping.")
                continue

            sample_out = os.path.join(shuffle_root, sample)
            os.makedirs(sample_out, exist_ok=True)

            for i in tqdm(range(1, n_shuffles + 1), desc=f"  Shuffling {sample}"):
                adata = sc.read_h5ad(h5_path)

                if "inferred_clone" not in adata.obs:
                    print(f"  'inferred_clone' missing. Skipping.")
                    break

                rng     = np.random.default_rng(random_seed + i if random_seed is not None else None)
                mask    = adata.obs["inferred_clone"] != "N"
                counts  = adata.obs.loc[mask, "inferred_clone"].value_counts()
                indices = adata.obs.loc[mask].index.to_numpy()
                rng.shuffle(indices)

                start = 0
                for label, count in counts.items():
                    adata.obs.loc[indices[start : start + count], "inferred_clone"] = label
                    start += count

                adata.write_h5ad(os.path.join(sample_out, f"{sample}_shuffle_{i}.h5ad"))

    print("\nShuffle stage complete.")


# ─────────────────────────────────────────────────────────────────────────────
# Stage 3 — GRN inference
# ─────────────────────────────────────────────────────────────────────────────

# ── Shared preprocessing ──────────────────────────────────────────────────────

def _preprocess(adata):
    sc.pp.filter_genes(adata, min_counts=1)
    sc.pp.normalize_per_cell(adata, key_n_counts="n_counts_all")
    filter_result = sc.pp.filter_genes_dispersion(
        adata.X, flavor="cell_ranger",
        n_top_genes=min(len(adata.var_names), 2000), log=False,
    )
    adata = adata[:, filter_result.gene_subset]
    sc.pp.normalize_per_cell(adata)
    adata.raw = adata
    adata.layers["raw_count"] = adata.raw.X.copy()
    sc.pp.log1p(adata)
    sc.pp.scale(adata)
    sc.tl.pca(adata, svd_solver="arpack")
    return adata


def _filter_hvg(expr, n_top_genes=8000):
    if expr.shape[1] <= n_top_genes:
        return expr
    ad = sc.AnnData(expr)
    sc.pp.filter_cells(ad, min_genes=200)
    sc.pp.filter_genes(ad, min_cells=3)
    if len(ad.var_names) <= n_top_genes:
        return expr[list(ad.var_names)]
    flavor_kwargs = {"span": 0.5} if len(ad.obs_names) <= 20 else {}
    sc.pp.highly_variable_genes(ad, n_top_genes=n_top_genes, flavor="seurat_v3", **flavor_kwargs)
    return expr[ad.var.index[ad.var["highly_variable"]].tolist()]


# ── Built-in worker functions ─────────────────────────────────────────────────

def _celloracle_worker(args):
    import celloracle as co
    h5ad_path, output_dir = args
    sample = os.path.basename(h5ad_path).replace(".h5ad", "")
    print(f"[CellOracle START] {sample}")

    adata = sc.read_h5ad(h5ad_path)
    if "inferred_clone" not in adata.obs:
        print(f"  'inferred_clone' missing in {sample}. Skipping.")
        return
    adata.obs["clone"] = adata.obs["inferred_clone"].astype(str)

    if all(
        os.path.exists(os.path.join(output_dir, f"CL{k}-{sample}-celloracle-net.tsv"))
        for k in adata.obs["clone"].unique()
    ):
        print(f"  All outputs exist for {sample}. Skipping.")
        return

    base_GRN = co.data.load_human_promoter_base_GRN()
    adata    = _preprocess(adata)

    oracle = co.Oracle()
    adata.X = adata.layers["raw_count"].copy()
    oracle.import_anndata_as_raw_count(adata=adata, cluster_column_name="clone", embedding_name="X_pca")
    oracle.import_TF_data(TF_info_matrix=base_GRN)
    oracle.perform_PCA()

    n_cell = oracle.adata.shape[0]
    k      = max(10, int(0.025 * n_cell))
    second_diff = np.diff(np.diff(np.cumsum(oracle.pca.explained_variance_ratio_)))
    candidates  = np.where(second_diff > 0.002)[0]
    n_comps     = int(candidates[0]) if len(candidates) > 0 else 50

    oracle.knn_imputation(n_pca_dims=n_comps, k=k, balanced=True, b_sight=k * 8, b_maxl=k * 4, n_jobs=4)
    net = oracle.get_links(cluster_name_for_GRN_unit="clone", alpha=10, verbose_level=0)

    for key, df in net.links_dict.items():
        df.to_csv(os.path.join(output_dir, f"CL{key}-{sample}-celloracle-net.tsv"), sep="\t")
    print(f"[CellOracle DONE] {sample}")


def _grnboost_worker(h5ad_path, output_dir):
    from arboreto.algo import grnboost2
    sample = os.path.basename(h5ad_path).replace(".h5ad", "")
    print(f"[GRNBoost START] {sample}")

    adata = sc.read_h5ad(h5ad_path)
    if "inferred_clone" not in adata.obs:
        print(f"  'inferred_clone' missing in {sample}. Skipping.")
        return

    for clone in adata.obs["inferred_clone"].unique():
        out_file = os.path.join(output_dir, f"CL{clone}-{sample}-grnboost-net.tsv")
        if os.path.exists(out_file):
            continue
        adata_sub = adata[adata.obs["inferred_clone"] == clone].copy()
        expr = pd.DataFrame(
            adata_sub.X.toarray() if hasattr(adata_sub.X, "toarray") else adata_sub.X,
            index=adata_sub.obs_names, columns=adata_sub.var_names,
        )
        expr = _filter_hvg(expr)
        if expr.shape[0] < 10:
            continue
        genes   = expr.columns.tolist()
        network = grnboost2(expression_data=expr, tf_names=genes)
        network.to_csv(out_file, sep="\t", index=False)
    print(f"[GRNBoost DONE] {sample}")


def _load_custom_worker(script_path: str):
    """
    Dynamically load a custom GRN script and return its worker callable.

    The script must define a top-level function named `run_grn_for_file`
    with the signature:

        def run_grn_for_file(h5ad_path: str, output_dir: str) -> None

    It should write one TSV per clone following the naming convention:
        CL{clone}-{sample}-{method}-net.tsv
    """
    spec   = importlib.util.spec_from_file_location("custom_grn_module", script_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    if not hasattr(module, "run_grn_for_file"):
        raise AttributeError(
            f"Custom script '{script_path}' must define a function "
            "`run_grn_for_file(h5ad_path, output_dir)`."
        )
    return module.run_grn_for_file


def _dispatch_grn_inference(method_cfg: dict, h5ad_files: list, output_dir: str, n_jobs: int):
    """Route to the correct inference backend for one method."""
    name   = method_cfg["name"]
    script = method_cfg.get("script", None)

    h5ad_shuffled = list(h5ad_files)
    random.shuffle(h5ad_shuffled)

    if script:
        # Custom method — run sequentially via the user-supplied worker
        worker = _load_custom_worker(script)
        for f in tqdm(h5ad_shuffled, desc=f"  [{name}] custom inference"):
            worker(f, output_dir)

    elif name == "celloracle":
        from multiprocessing import Pool
        work_items = [(f, output_dir) for f in h5ad_shuffled]
        with Pool(processes=n_jobs) as pool:
            list(tqdm(pool.imap_unordered(_celloracle_worker, work_items), total=len(work_items)))

    elif name == "grnboost":
        from joblib import Parallel, delayed
        Parallel(n_jobs=n_jobs, backend="loky")(
            delayed(_grnboost_worker)(f, output_dir) for f in h5ad_shuffled
        )

    else:
        raise ValueError(
            f"Unknown method '{name}'. Either use 'celloracle', 'grnboost', 'combined', "
            "or provide a 'script' path in the config."
        )


def run_grn_stage(datasets: dict, data_dir: str, grn_dir: str, methods: list, cfg: dict):
    """Run GRN inference for every dataset × method combination."""
    n_jobs = cfg.get("n_jobs", 4)

    print("\n" + "=" * 60)
    print("STAGE 3 — GRN inference")
    print("=" * 60)

    for dataset, samples in datasets.items():
        boot_root = os.path.join(data_dir, dataset, f"bootstrap_{dataset}")

        for method_cfg in methods:
            name = method_cfg["name"]
            if name == "combined":
                continue  # combined is a support-scoring step, not inference

            output_dir = os.path.join(grn_dir, dataset, name)
            os.makedirs(output_dir, exist_ok=True)

            # Collect all bootstrap h5ad files across all samples for this dataset
            h5ad_files = []
            for sample in samples:
                sample_dir = os.path.join(boot_root, sample)
                if not os.path.isdir(sample_dir):
                    print(f"  Bootstrap dir missing for {sample}. Skipping.")
                    continue
                h5ad_files += [
                    os.path.join(sample_dir, f)
                    for f in os.listdir(sample_dir)
                    if f.endswith(".h5ad")
                ]

            if not h5ad_files:
                print(f"  No h5ad files found for dataset '{dataset}'. Skipping.")
                continue

            print(f"\n  Dataset: {dataset} | Method: {name} | Files: {len(h5ad_files)}")
            _dispatch_grn_inference(method_cfg, h5ad_files, output_dir, n_jobs)

    print("\nGRN inference stage complete.")


# ─────────────────────────────────────────────────────────────────────────────
# Stage 4 — Bootstrap support scoring
# ─────────────────────────────────────────────────────────────────────────────

def _load_grn_as_graph(network: pd.DataFrame, method: str) -> nx.DiGraph:
    if method == "celloracle":
        edges = network.loc[network["p"] < 0.05]
        return nx.from_pandas_edgelist(
            edges, source="source", target="target",
            edge_attr="coef_mean", create_using=nx.DiGraph,
        )
    if method == "grnboost":
        thresh = np.histogram(network["importance"])[1][1]
        edges  = network.loc[network["importance"] >= thresh]
        return nx.from_pandas_edgelist(
            edges, source="TF", target="target",
            edge_attr="importance", create_using=nx.DiGraph,
        )
    raise ValueError(f"Cannot auto-parse GRN for method '{method}'. "
                     "Custom methods should supply their own graph loader.")


def _flag_unique_edges(G_df: dict, all_grns: dict, clones: list) -> dict:
    clone_unique = {}
    for clone in clones:
        unique = set(all_grns[clone].edges())
        for other in clones:
            if other != clone:
                unique -= set(all_grns[other].edges())
        clone_unique[clone] = unique

    for clone in clones:
        G_df[clone]["edge"]   = list(zip(G_df[clone]["source"], G_df[clone]["target"]))
        G_df[clone]["unique"] = G_df[clone]["edge"].isin(clone_unique[clone])
    return G_df


def _compute_support_for_method(method: str, grn_dir: str, clones: list, rate: float, n_bootstraps: int) -> dict:
    """Load bootstrap + whole GRNs, accumulate support counts, return per-clone DataFrames."""
    support_key = f"support_{rate}"

    boot_grns = {clone: [] for clone in clones}
    for clone in clones:
        for i in range(1, n_bootstraps + 1):
            path    = os.path.join(grn_dir, f"CL{clone}-{i}-{method}-net.tsv")
            network = pd.read_csv(path, sep="\t")
            boot_grns[clone].append(_load_grn_as_graph(network, method))

    all_grns = {}
    for clone in clones:
        path    = os.path.join(grn_dir, f"CL{clone}-whole-{method}-net.tsv")
        network = pd.read_csv(path, sep="\t")
        all_grns[clone] = _load_grn_as_graph(network, method)

    for clone in clones:
        for G_boot in boot_grns[clone]:
            for u, v in G_boot.edges():
                if all_grns[clone].has_edge(u, v):
                    data = all_grns[clone].edges[u, v]
                    data[support_key] = data.get(support_key, 0) + 1

    G_df = {}
    for clone in clones:
        rows = {"source": [], "target": [], support_key: []}
        for u, v, data in all_grns[clone].edges(data=True):
            rows["source"].append(u)
            rows["target"].append(v)
            rows[support_key].append(data.get(support_key, 0))
        G_df[clone] = pd.DataFrame(rows)

    return _flag_unique_edges(G_df, all_grns, clones)


def _combine_methods(G_df_a: dict, G_df_b: dict, clones: list, rate: float) -> dict:
    """Intersect two methods' results and average their support scores."""
    support_key = f"support_{rate}"
    all_grns_combined = {}
    G_df_avg = {}

    for clone in clones:
        merged = pd.merge(G_df_a[clone], G_df_b[clone], how="inner", on=["source", "target"], suffixes=("_a", "_b"))
        merged[support_key] = (merged[f"{support_key}_a"] + merged[f"{support_key}_b"]) / 2
        drop = [c for c in merged.columns if c.endswith(("_a", "_b")) or c in ("edge_x", "edge_y", "unique_x", "unique_y")]
        G_df_avg[clone] = merged.drop(columns=drop, errors="ignore")

        G_a = nx.from_pandas_edgelist(G_df_a[clone], source="source", target="target", create_using=nx.DiGraph)
        G_b = nx.from_pandas_edgelist(G_df_b[clone], source="source", target="target", create_using=nx.DiGraph)
        all_grns_combined[clone] = nx.intersection(G_a, G_b)

    return _flag_unique_edges(G_df_avg, all_grns_combined, clones)


def run_support_stage(datasets: dict, grn_dir: str, support_dir: str, methods: list, cfg: dict):
    """Compute bootstrap support scores for every dataset × method combination."""
    clones       = cfg.get("clones", [])
    rate         = cfg.get("rate", 0.9)
    n_bootstraps = cfg.get("n_bootstraps", 100)

    if not clones:
        raise ValueError("bootstrap_support.clones must be set in the config.")

    print("\n" + "=" * 60)
    print("STAGE 4 — Bootstrap support scoring")
    print("=" * 60)

    # Separate out 'combined' so we can handle it after per-method results are ready
    single_methods  = [m for m in methods if m["name"] != "combined"]
    run_combined    = any(m["name"] == "combined" for m in methods)

    for dataset in datasets:
        print(f"\n  Dataset: {dataset}")
        dataset_support_dir = os.path.join(support_dir, dataset)
        os.makedirs(dataset_support_dir, exist_ok=True)

        per_method_results = {}

        for method_cfg in single_methods:
            name    = method_cfg["name"]
            method_grn_dir = os.path.join(grn_dir, dataset, name)
            print(f"    Method: {name}")

            G_df = _compute_support_for_method(
                method=name,
                grn_dir=method_grn_dir,
                clones=clones,
                rate=rate,
                n_bootstraps=n_bootstraps,
            )
            per_method_results[name] = G_df

            out_path = os.path.join(dataset_support_dir, f"G_df_{name}.pkl")
            with open(out_path, "wb") as f:
                pickle.dump(G_df, f)
            print(f"    Saved: {out_path}")

        if run_combined:
            method_names = [m["name"] for m in single_methods]
            if len(method_names) < 2:
                print("    WARNING: 'combined' requires at least 2 non-combined methods. Skipping.")
            else:
                # Combine the first two available methods (typically celloracle + grnboost)
                name_a, name_b = method_names[0], method_names[1]
                print(f"    Combining: {name_a} ∩ {name_b}")
                G_df_combined = _combine_methods(
                    per_method_results[name_a],
                    per_method_results[name_b],
                    clones=clones,
                    rate=rate,
                )
                out_path = os.path.join(dataset_support_dir, "G_df_combined.pkl")
                with open(out_path, "wb") as f:
                    pickle.dump(G_df_combined, f)
                print(f"    Saved: {out_path}")

    print("\nBootstrap support stage complete.")


# ─────────────────────────────────────────────────────────────────────────────
# Pipeline entry point
# ─────────────────────────────────────────────────────────────────────────────

def main():
    args   = parse_args()
    config = load_config(args.config)

    # Unpack config
    datasets    = {name: info["samples"] for name, info in config["datasets"].items()}
    data_dir    = config["output"]["data_dir"]
    grn_dir     = config["output"]["grn_dir"]
    support_dir = config["output"]["support_dir"]
    methods     = config.get("methods", [{"name": "celloracle"}, {"name": "grnboost"}])

    for d in [data_dir, grn_dir, support_dir]:
        os.makedirs(d, exist_ok=True)

    override = args.stages  # None means "use config run flags"

    if should_run("bootstrap", config, override):
        run_bootstrap_stage(datasets, data_dir, config["stages"]["bootstrap"])

    if should_run("shuffle", config, override):
        run_shuffle_stage(datasets, data_dir, config["stages"]["shuffle"])

    if should_run("grn_inference", config, override):
        run_grn_stage(datasets, data_dir, grn_dir, methods, config["stages"]["grn_inference"])

    if should_run("bootstrap_support", config, override):
        run_support_stage(datasets, grn_dir, support_dir, methods, config["stages"]["bootstrap_support"])

    print("\n" + "=" * 60)
    print("Pipeline complete.")
    print("=" * 60)


if __name__ == "__main__":
    main()
