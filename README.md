# Intratumor Heterogeneity Through the Lens of Gene Regulatory Networks — analysis code

Code to reproduce the analyses and figures of the manuscript. Per-clone gene regulatory networks (GRNs) are inferred with CellOracle from
single-cell RNA-seq, with bootstrap and label-shuffle controls. The pipeline takes any dataset from the Curated Cancer Cell Atlas (3CA) that has per-cell copy-number profiles: copy it into the data tree and steps 2-5 produce its clone-labelled h5ads, GRNs, bootstrap support and figures.
The manuscript uses 253 patients across 12 cancer types (list in `data_manifests/cohort_253_patients.csv`), plus two longitudinal cohorts: a stage IV breast cancer patient (pre-/post-treatment) and the ovarian cancer cohort of Zhang et al. 2022 (GSE165897, pre-/post-neoadjuvant chemotherapy).

## Pipeline

```
 raw scRNA-seq + CNA profiles
        │
  [1] EOC only: split GSE165897 → inferCNV → per-patient h5ad
        │
  [2] clone calling   → <Patient>_clones.h5ad   (obs['clone'] = 0,1,2,… or 'N' for normal)
        │
  [3] bootstrap/shuffle definitions → CellOracle (201 units per patient, SLURM) → gather  → .uns['celloracle_*']
        │                                   └→ [3b] SCENIC regulons on the CellOracle GRNs
  [4] bootstrap support (p < 0.05 edges)   → <Patient>_clones_with_support.h5ad
        │
  [5] analysis + figure scripts → PNGs → composite figures
  [6] regulons / pathways / meta-program networks / survival, and the stage IV longitudinal analysis (notebooks)
```

## Layout

| Path | Contents |
|---|---|
| `scripts/01_eoc_infercnv/` | GSE165897 splitting, inferCNV (parameters of Zhang et al. 2022), per-patient h5ad builder |
| `scripts/02_clone_calling/` | `clone_sweep.py` (3CA), `clone_sweep_zhang2022.py` (EOC), `make_pooled_tumor_h5ad.py`, helper scripts |
| `scripts/03_grn_inference/` | bootstrap/shuffle definitions, CellOracle worker and SLURM scripts, gather step, SCENIC |
| `scripts/04_bootstrap_support/` | per-edge bootstrap support |
| `scripts/05_analysis_and_figures/` | cohort-wide figures and statistics, treatment-phase analysis, composite figures |
| `notebooks/pathway_regulon_survival.ipynb` | regulons, pathways, meta-program networks, centrality, TCGA survival (Figs 7–8) |
| `notebooks/longitudinal_stage4_breast.ipynb` | stage IV breast cancer: pre- vs post-treatment GRN divergence, neighbourhood similarity, centrality (Fig 5) |
| `data_manifests/cohort_253_patients.csv` | the 253 patients (`cancer, dataset, sample, n_clones, tag`) |
| `data_manifests/distinctness_classified.csv` | per-sample result of the visual distinctness pre-filter on the 3CA cohort |
| `environment.yml`, `r-requirements.md` | software environments |

Each script starts with a one-line description and its run command.

## Software

```bash
conda env create -f environment.yml && conda activate celloracle_env
```
EOC inferCNV additionally needs R with inferCNV **1.4.0** (`r-requirements.md`). The cluster scripts (`*.slurm`, `controller.sh`, `submit_*.sh`) were written for a SLURM cluster; edit the account, partition and environment lines.

## Data needed

| Dataset | What | Source |
|---|---|---|
| 3CA scRNA-seq + per-cell CNAs (any 3CA dataset; 253 patients, 12 cancer types in the manuscript) | `Exp_data_*.mtx`, gene and cell files, per-cell CNA matrix, laid out as `data/<Cancer>/<Dataset>/`, as `clone_sweep.py` expects (`--inspect` checks it) | Curated Cancer Cell Atlas |
| GSE165897 (Zhang et al. 2022) | `GSE165897_UMIcounts_HGSOC.tsv`, `GSE165897_cellInfo_HGSOC.tsv` | GEO |
| GENCODE v25 gene order file | for EOC inferCNV and `clone_sweep_zhang2022.py` | gencodegenes.org |
| hg38 gene order file (`hg38_gencode_v27.txt`) | for `clone_sweep.py` | GENCODE v27 |
| TCGA-BRCA / TCGA-GBM expression and clinical data | fetched live in the survival section of the notebook | UCSC Xena |
| cisTarget motif table and hg38 ranking databases | for the SCENIC step | aertslab |

No data are shipped (`.gitignore` excludes `*.h5ad`, `*.npz`).

## Running it

Run each command from the directory containing the script.

### 1. EOC only — inferCNV and per-patient h5ads
```bash
cd scripts/01_eoc_infercnv
python split_gse165897_by_patient.py --counts GSE165897_UMIcounts_HGSOC.tsv --annotations GSE165897_cellInfo_HGSOC.tsv --all --out_dir split_patients
MAX_PARALLEL=10 bash run_infercnv_all_patients_parallel.sh          # wraps run_infercnv_zhang2022.R; resumable
python build_h5ad_per_patient.py --split_dir split_patients --full_annotations GSE165897_cellInfo_HGSOC.tsv --infercnv_dir infercnv_out --all --out_dir h5ad_out
# if your inferCNV version does not write infercnv.observations.txt:  Rscript extract_cnv_from_infercnv_obj.R infercnv_out/<patient>
```

### 2. Clone calling
Clones are called with one automated procedure for all patients: scale → PCA(20) → neighbours(20) → Leiden at resolution 1.0 (over-clustering) → iterative merging of the closest cluster pair
unless a permutation test (200 permutations, Bonferroni-adjusted α = 0.05) shows the pair is distinct; CNAs aggregated to chromosome arms; clones with < 50 cells dropped; patients with < 2 clones skipped.
```bash
cd scripts/02_clone_calling
python clone_sweep.py --inspect
python clone_sweep.py --data-root data --gene-order hg38_gencode_v27.txt --materialize --out clones_out     # clones_out/<Cancer>/<Dataset>/<Sample>_clones.h5ad
```
To add a dataset, put its 3CA folder at `data/<Cancer>/<Dataset>/` and rerun this command (restrict to it with `--cancer <Cancer> --dataset <Dataset>`); steps 3-5 then pick up every `*_clones.h5ad` under `clones_out`.
EOC (same procedure on the inferCNV-derived profiles):
```bash
python clone_sweep_zhang2022.py --h5ad_dir h5ad_out --gene-order gencode_v25_gene_order.txt --materialize --out clones_out
cd clones_out && python <repo>/scripts/02_clone_calling/make_clone_copies.py      # writes <Patient>_clones.h5ad with the clone-sweep clones
mkdir -p Ovarian/Zhang2022 && mv *_clones.h5ad Ovarian/Zhang2022/
```

### 3. CellOracle GRNs with bootstrap and shuffle controls
```bash
cd scripts/03_grn_inference
python add_bootstrap_shuffle.py --root clones_out --inplace
python add_bootstrap_shuffle_zhang2022.py --root clones_out/Ovarian --inplace
DATA_ROOT=clones_out RESULTS_ROOT=grn_results nohup bash controller.sh >> controller.log 2>&1 &
python check_grn_completeness.py --root clones_out
```
Each patient gets 100 bootstrap resamples (90 % of each clone's cells) and 100 label shuffles, run as 201 CellOracle units on SLURM and gathered into the patient h5ad
(`.uns['celloracle_original' | '_bootstrap' | '_shuffle']`; edge-list columns `[clone, source, target, coef_mean, coef_abs, p, -logp]`, coded against `clone_categories` / `celloracle_genes`).
A single unit can be run directly: `python run_celloracle_task.py --h5ad P_clones.h5ad --results-dir DIR --kind {original|bootstrap|shuffle} [--rep i]`.

**3b. SCENIC regulons on the CellOracle GRNs** (used for the meta-program networks of Figs 7–8):
```bash
python run_scenic_from_celloracle.py --h5ad clones_out/Breast/Wu_et_al_2021/CID4495_clones.h5ad --out scenic_out \
    --motif-file motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl --db-glob '/path/to/hg38/*feather' --export-celloracle
```
Writes `scenic_out/<Dataset>/<Sample>/CL<clone>-whole-scenic-net.tsv` (`Gene1, Gene2, weight`) and, with `--export-celloracle`, `CL<clone>-{whole|1..100}-celloracle-net.tsv`.

### 4. Bootstrap support
```bash
cd scripts/04_bootstrap_support
python compute_bootstrap_support.py --root clones_out          # or: bash submit_bootstrap_support.sh  (SLURM array)
python validate_bootstrap_results.py --root clones_out
```
Writes `*_with_support.h5ad` next to each patient file; the source files are never modified. The ovarian cohort is analysed separately (`--root clones_out/Ovarian --exclude-cancer-type`).

### 5. Analyses and figures
```bash
cd scripts/05_analysis_and_figures
python bootstrap_figures_efficient.py --root clones_out --out bootstrap_figs_out
python intra_vs_inter_efficient.py    --root clones_out --out intra_vs_inter_out
python plot_jaccard_clone_vs_normal.py --root clones_out --out plots/
python treatment_phase_comparison_zhang2022.py --root clones_out/Ovarian --out treatment_phase_out --pre-label treatment-naive --post-label post-NACT --min-abs-cells 50
python make_composite_figures.py --bootstrap-figs bootstrap_figs_out --intra-inter intra_vs_inter_out --treatment-phase treatment_phase_out --out figures
python make_supplementary_composites.py --bootstrap-figs bootstrap_figs_out --out supp_figs
```
`make_composite_figures.py` writes `Fig2.png`, `Fig3.png` and `Fig6.png`.

### 6. Regulons, pathways, meta-program networks, centrality and survival
`notebooks/pathway_regulon_survival.ipynb` (Figs 7–8, Supplementary Figs S5, S7–S9). Set `DATA_ROOT` and `GENE_ORDER_FILE` in the first code cell (or `ITH_DATA_ROOT` / `ITH_GENE_ORDER`), then run it once for TNBC and once for glioblastoma
by editing `SAMPLE` / `specimen` / `datasets`. It uses the SCENIC and CellOracle edge lists from step 3b and the gene-set tables listed in the notebook's first cell, and expects the h5ad conventions listed there (`obs['inferred_clone']`, `uns['cnp']`).

### 7. Stage IV breast cancer longitudinal analysis
Run the stage IV samples (`pre`, `post1`, `post2`) as a dataset through steps 2-4 into their own root, then create the pooled-tumour copies (all malignant cells as one clone `T`, normal `N`) and run steps 3-4 on them too:
```bash
python scripts/02_clone_calling/make_pooled_tumor_h5ad.py --in-dir stage4_out/Breast/Stage4 --out-dir stage4_out/Breast/Stage4_pooled
```
`notebooks/longitudinal_stage4_breast.ipynb` (Fig 5, Supplementary Fig S3) reads `*_clones_with_support.h5ad` from `stage4_out/Breast/Stage4` and `stage4_out/Breast/Stage4_pooled` (set `ROOT`, `CLONE_DIR`, `POOLED_DIR` in its first cell, or `ITH_STAGE4_ROOT`).
Keep this root separate from `clones_out` so the cohort-wide scripts do not pick it up.

## Figure map

| Manuscript figure | Produced by |
|---|---|
| Figs 1, 4 (TNBC/CRC clone-vs-normal, intra/inter-clone) | `bootstrap_figures_efficient.py`, `intra_vs_inter_efficient.py`, `plot_jaccard_clone_vs_normal.py` |
| Fig 2 (extended cohort: clone vs normal) | `bootstrap_figures_efficient.py` → `make_composite_figures.py` (`Fig2.png`) |
| Fig 3 (extended cohort: intra/inter, shuffle, clone-size divergence) | `intra_vs_inter_efficient.py` + `bootstrap_figures_efficient.py` → `make_composite_figures.py` (`Fig3.png`) |
| Fig 5 (stage IV breast, longitudinal) | `notebooks/longitudinal_stage4_breast.ipynb` |
| Fig 6 (pre- vs post-NACT, EOC) | `treatment_phase_comparison_zhang2022.py` → `make_composite_figures.py` (`Fig6.png`) |
| Figs 7–8, Supplementary Figs S5, S7–S9 | `notebooks/pathway_regulon_survival.ipynb` |
| Supplementary Figs S11–S18 | `bootstrap_figures_efficient.py`, `intra_vs_inter_efficient.py` → `make_supplementary_composites.py` |
| Supplementary Tables S1–S3 | `data_manifests/cohort_253_patients.csv` |


## Citing / data availability

