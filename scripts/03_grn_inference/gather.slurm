#!/bin/bash
# Folds one patient's replicate GRNs into its h5ad.
# Usage: sbatch --dependency=afterany:<arrayjobid> gather.slurm <patient_clones.h5ad> <results_dir>
#SBATCH --job-name=co_gather
#SBATCH --account=commons
#SBATCH --partition=commons
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=logs/%x_%j.out

set -uo pipefail

H5AD="$1"; RESULTS="$2"

# activate the env by explicit path (see celloracle_array.slurm for why)
module load Mamba/23.11.0-0
export PATH="$HOME/.conda/envs/celloracle_env/bin:$PATH"
export PYTHONNOUSERSITE=1
export LD_LIBRARY_PATH="$HOME/.conda/envs/celloracle_env/lib:${LD_LIBRARY_PATH:-}"

# per-task private cache so concurrent tasks never share genomepy/diskcache SQLite
export XDG_CACHE_HOME="/tmp/gpcache_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID:-0}"
mkdir -p "$XDG_CACHE_HOME"

python gather_grn.py --h5ad "$H5AD" --results-dir "$RESULTS"

# sentinel so the launcher can skip this patient on re-runs (fast, no h5ad load)
if python -c "import anndata,sys; a=anndata.read_h5ad('$H5AD'); sys.exit(0 if 'celloracle_original' in a.uns else 1)"; then
    touch "$RESULTS/.gathered"
    echo "gathered OK -> touched $RESULTS/.gathered"
fi
