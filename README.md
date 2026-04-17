# Intratumor-Heterogeneity-Through-the-Lens-of-Gene-Regulatory-Networks

An end-to-end pipeline for gene regulatory network (GRN) inference from single-cell RNA-seq data with clonal annotations. The pipeline bootstraps and shuffles input data, runs GRN inference in parallel, and computes per-edge bootstrap support scores across one or more datasets and methods.

---

## Overview

The pipeline runs four sequential stages:

```
Input .h5ad files
       │
       ▼
┌─────────────────────┐
│  1. Bootstrap data  │  Sample 90% of cells per clone × 100 replicates
└─────────────────────┘
       │
       ▼
┌─────────────────────┐
│  2. Shuffle data    │  Randomly permute clone labels to generate null replicates
└─────────────────────┘
       │
       ▼
┌─────────────────────┐
│  3. GRN inference   │  Run CellOracle, GRNBoost2, or a custom method in parallel
└─────────────────────┘
       │
       ▼
┌─────────────────────┐
│  4. Bootstrap       │  Score each edge in the whole-sample GRN by how often
│     support         │  it appears across bootstrap replicates
└─────────────────────┘
       │
       ▼
Per-clone DataFrames saved as .pkl files
```

Each stage can be run independently or toggled on/off in the config file.

---

## File Structure

```
grn_pipeline.py        # Main pipeline runner
pipeline_config.yaml   # User-facing configuration (edit this)
```

The individual module scripts (`create_bootstrap_data.py`, `create_shuffle_data.py`, `run_grn_parallel.py`, `compute_bootstrap_support.py`) are not required at runtime — their logic is fully integrated into `grn_pipeline.py`.

---

## Requirements

```
scanpy
numpy
pandas
networkx
pyyaml
tqdm
joblib
```

Method-specific dependencies (only needed if the corresponding method is enabled):

| Method | Package |
|---|---|
| `celloracle` | `celloracle` |
| `grnboost` | `arboreto` |

---

## Quick Start

**1. Edit the config** to point to your data and set your parameters:

```yaml
datasets:
  My_Dataset:
    samples:
      sample_A: /path/to/sample_A_clone_correct.h5ad
      sample_B: /path/to/sample_B_clone_correct.h5ad

output:
  data_dir:    /path/to/output/data
  grn_dir:     /path/to/output/grns
  support_dir: /path/to/output/support

methods:
  - name: celloracle
  - name: grnboost

stages:
  bootstrap_support:
    clones: [N, C1, C2]   # must match inferred_clone values in your .h5ad files
```

**2. Run the full pipeline:**

```bash
python grn_pipeline.py --config pipeline_config.yaml
```

**3. Run specific stages only:**

```bash
python grn_pipeline.py --config pipeline_config.yaml --stages bootstrap grn_inference
```

---

## Configuration Reference

### `datasets`

Defines one or more datasets to process. Each dataset key is a free-form label used to name output subdirectories. `samples` maps sample names to their source `.h5ad` file paths.

The `.h5ad` files must have an `inferred_clone` column in `adata.obs`.

```yaml
datasets:
  Dataset_Name:
    samples:
      sample_name: /path/to/file_clone_correct.h5ad
```

### `output`

All output directories. Subdirectories per dataset and method are created automatically.

| Key | Contents |
|---|---|
| `data_dir` | Bootstrap and shuffle `.h5ad` replicates |
| `grn_dir` | Inferred GRN `.tsv` files |
| `support_dir` | Bootstrap support score `.pkl` files |

### `methods`

List of GRN inference methods to run. Built-in options are `celloracle`, `grnboost`, and `combined`. See the [Custom Methods](#custom-methods) section to add your own.

```yaml
methods:
  - name: celloracle
  - name: grnboost
  - name: combined    # computes the intersection of celloracle and grnboost results
```

`combined` is only relevant to Stage 4 (bootstrap support). It intersects the edge sets of the two single methods and averages their support scores. It is silently ignored during Stage 3 (GRN inference).

### `stages`

Each stage has a `run` boolean flag and its own parameters.

#### `bootstrap`

| Parameter | Default | Description |
|---|---|---|
| `run` | `true` | Whether to run this stage |
| `n_bootstraps` | `100` | Number of bootstrap replicates per sample |
| `cell_proportion` | `0.9` | Fraction of cells sampled per clone per replicate |

#### `shuffle`

| Parameter | Default | Description |
|---|---|---|
| `run` | `true` | Whether to run this stage |
| `n_shuffles` | `100` | Number of shuffle replicates per sample |
| `random_seed` | `42` | RNG seed; set to `null` for non-deterministic results |

#### `grn_inference`

| Parameter | Default | Description |
|---|---|---|
| `run` | `true` | Whether to run this stage |
| `n_jobs` | `4` | Number of parallel worker processes |

#### `bootstrap_support`

| Parameter | Default | Description |
|---|---|---|
| `run` | `true` | Whether to run this stage |
| `rate` | `0.9` | Bootstrap sampling rate label (should match `bootstrap.cell_proportion`) |
| `n_bootstraps` | `100` | Number of replicates to load (should match `bootstrap.n_bootstraps`) |
| `clones` | — | **Required.** List of clone labels present in `inferred_clone` |

---

## Output Structure

```
data_dir/
└── {dataset}/
    ├── bootstrap_{dataset}/
    │   └── {sample}/
    │       ├── {sample}_bootstrap_1.h5ad
    │       ├── {sample}_bootstrap_2.h5ad
    │       └── ...
    └── shuffle_{dataset}/
        └── {sample}/
            ├── {sample}_shuffle_1.h5ad
            └── ...

grn_dir/
└── {dataset}/
    └── {method}/
        ├── CL{clone}-{sample}_bootstrap_1-{method}-net.tsv
        ├── CL{clone}-{sample}_bootstrap_2-{method}-net.tsv
        └── ...

support_dir/
└── {dataset}/
    ├── G_df_celloracle.pkl
    ├── G_df_grnboost.pkl
    └── G_df_combined.pkl      # only if combined is in methods
```

### Support score pickles

Each `.pkl` file is a `dict[clone -> pd.DataFrame]` with the following columns:

| Column | Description |
|---|---|
| `source` | Source gene (TF) |
| `target` | Target gene |
| `support_{rate}` | Number of bootstrap replicates in which this edge was present (0–100) |
| `edge` | `(source, target)` tuple |
| `unique` | `True` if this edge does not appear in any other clone's GRN |

---

## Custom Methods

To use a GRN inference method other than CellOracle or GRNBoost2, provide a Python script that defines the following function:

```python
def run_grn_for_file(h5ad_path: str, output_dir: str) -> None:
    ...
```

The function receives the path to a single `.h5ad` bootstrap replicate and the directory where output TSV files should be written. Output files must follow the naming convention:

```
CL{clone}-{sample}-{method}-net.tsv
```

Then register it in the config:

```yaml
methods:
  - name: my_method
    script: /path/to/my_grn_method.py
```

The pipeline will import the script dynamically and call `run_grn_for_file` for each bootstrap file. The method name (`my_method`) is used as a directory label and as the `{method}` token in output filenames.

> **Note:** Bootstrap support scoring (Stage 4) uses built-in edge-filtering logic for `celloracle` and `grnboost`. For a custom method, you will need to extend `_load_grn_as_graph()` in `grn_pipeline.py` to handle your TSV format, or ensure your TSV follows the CellOracle (`source`, `target`, `p`, `coef_mean`) or GRNBoost2 (`TF`, `target`, `importance`) column schema.

---

## Input Requirements

Each `.h5ad` file must satisfy the following:

- `adata.obs` contains an `inferred_clone` column with string clone labels.
- At least one clone must be labelled `"N"` (the normal/reference clone).
- Expression data is stored in `adata.X` (raw counts are expected prior to normalisation; the pipeline normalises internally for CellOracle).

---

## Resuming a Partial Run

Each stage checks whether its outputs already exist before recomputing:

- **Bootstrap / Shuffle:** Output `.h5ad` files are written per replicate; already-written files are not regenerated.
- **GRN inference:** Per-clone TSV outputs are checked before running inference on a file; existing files are skipped.
- **Bootstrap support:** Recomputed from scratch per run. Re-run only this stage with `--stages bootstrap_support` if needed.

To re-run only specific stages:

```bash
# Re-run GRN inference and support scoring only
python grn_pipeline.py --config pipeline_config.yaml --stages grn_inference bootstrap_support
```
