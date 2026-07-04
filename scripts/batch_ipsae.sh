#!/bin/bash
#SBATCH --job-name=ipsae_array
#SBATCH --partition=cpu              # adjust
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --array=0-1149             # set XXX = number_of_dirs-1
#SBATCH -o logs/ipsae/ipsae_array_%A_%a.out
#SBATCH -e logs/ipsae/ipsae_array_%A_%a.err

set -euo pipefail

mkdir -p ./logs/ipsae/
# -----------------------------
# HARD-SET CONFIGURATION
# -----------------------------
# Base directory containing complex subdirectories
BASE_DIR="/data1/choderaj/deepa2/rotation_project/alphapulldown/STK11/structures_output"

# Path to ipsae.py
IPSAE_SCRIPT="/data1/choderaj/deepa2/scripts/ipsae.py"

# AlphaFold2 recommended cutoffs
PAE_CUTOFF=15
DIST_CUTOFF=15

# Env setup if needed
# module load ...
source activate AlphaPulldown
# -----------------------------

echo "SLURM_ARRAY_TASK_ID = ${SLURM_ARRAY_TASK_ID}"
echo "BASE_DIR = ${BASE_DIR}"

# Build list of complex directories
mapfile -t DIRS < <(ls -d "${BASE_DIR}"/*/)

NUM_DIRS=${#DIRS[@]}
if [[ "${SLURM_ARRAY_TASK_ID}" -ge "${NUM_DIRS}" ]]; then
    echo "Task ID ${SLURM_ARRAY_TASK_ID} >= number of dirs ${NUM_DIRS}, exiting."
    exit 0
fi

DIR="${DIRS[$SLURM_ARRAY_TASK_ID]}"
echo "Processing directory: ${DIR}"

cd "${DIR}"

# For each result_model_i .pkl, call ipsae.py with the matching unrelaxed_model_i .pdb
for PKL in result_model_*_multimer_v3_pred_*.pkl; do
    [[ -e "$PKL" ]] || continue

    FNAME=$(basename "$PKL")
    IDX="${FNAME#result_model_}"
    IDX="${IDX%%_multimer_v3_pred_*}"

    PDB="unrelaxed_model_${IDX}_multimer_v3_pred_0.pdb"
    if [[ ! -f "$PDB" ]]; then
        echo "  Missing PDB for model ${IDX}, skipping."
        continue
    fi

    echo "  Running ipSAE on model ${IDX}: ${PKL} + ${PDB}"
    python "${IPSAE_SCRIPT}" "${PKL}" "${PDB}" "${PAE_CUTOFF}" "${DIST_CUTOFF}"
    # No parsing, no redirection: ipsae.py writes its own outputs.
done

echo "Done with directory: ${DIR}"

