#!/bin/bash
# One array task = bootstrap support for one patient (array range set by submit_bootstrap_support.sh).
# Usage: sbatch --array=0-N bootstrap_support_array.slurm patients_to_process.txt
#SBATCH --job-name=boot_support
#SBATCH --output=logs/boot_support_%A_%a.out
#SBATCH --error=logs/boot_support_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00

module load Mamba/23.11.0-0
export PATH="$HOME/.conda/envs/celloracle_env/bin:$PATH"
export PYTHONNOUSERSITE=1
export LD_LIBRARY_PATH="$HOME/.conda/envs/celloracle_env/lib:${LD_LIBRARY_PATH:-}"

FILE_LIST="$1"

python compute_bootstrap_support.py \
    --file-list "$FILE_LIST" \
    --array-index "$SLURM_ARRAY_TASK_ID"
