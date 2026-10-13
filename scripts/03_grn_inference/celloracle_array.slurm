#!/bin/bash
# One array task = one CellOracle replicate of one patient (0 original, 1-100 bootstrap, 101-200 shuffle).
# Usage: sbatch --array=0-200%15 celloracle_array.slurm <patient_clones.h5ad> <results_dir>
#SBATCH --job-name=co_grn
#SBATCH --account=commons
#SBATCH --partition=commons
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=06:00:00          # one replicate; well inside the 1-day QOS
#SBATCH --output=logs/%x_%A_%a.out

set -uo pipefail

H5AD="$1"; RESULTS="$2"
t="${SLURM_ARRAY_TASK_ID:?run via sbatch --array}"
if   [ "$t" -eq 0   ]; then KIND=original;  REP=-1
elif [ "$t" -le 100 ]; then KIND=bootstrap; REP=$((t-1))
else                        KIND=shuffle;   REP=$((t-101)); fi

# activate the env by explicit path (do NOT rely on `conda activate`, and do NOT
# `module load GCCcore` -- both caused the pkg_resources / xgboost libstdc++ breakage)
module load Mamba/23.11.0-0
export PATH="$HOME/.conda/envs/celloracle_env/bin:$PATH"
export PYTHONNOUSERSITE=1                                   # ignore ~/.local leak
export LD_LIBRARY_PATH="$HOME/.conda/envs/celloracle_env/lib:${LD_LIBRARY_PATH:-}"

# per-task private cache so concurrent tasks never share genomepy/diskcache SQLite
export XDG_CACHE_HOME="/tmp/gpcache_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID:-0}"
mkdir -p "$XDG_CACHE_HOME"

# let CellOracle/sklearn use exactly the cores we requested
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export NUMEXPR_NUM_THREADS=$SLURM_CPUS_PER_TASK

python run_celloracle_task.py \
    --h5ad "$H5AD" --results-dir "$RESULTS" --kind "$KIND" --rep "$REP"
