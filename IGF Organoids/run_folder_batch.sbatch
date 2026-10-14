#!/bin/bash
#SBATCH --job-name=IMS_Mosaic_Fields
#SBATCH --account=parent0
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=1-00:00:00
#SBATCH --output=logs/ims_%A_%a.out
#SBATCH --error=logs/ims_%A_%a.err

set -euo pipefail
: "${SLURM_ARRAY_TASK_ID:?Please submit using sbatch --array=0 or --array=0-N%2}"
: "${1:?Usage: sbatch --array=0-N run_folder_batch.sbatch <absolute_manifest.json> <absolute_run_folder_batch.py>}"
: "${2:?Pass the absolute run_folder_batch.py path as the second argument}"

module load singularity/4.4.1
export SINGULARITY_BINDPATH="/gpfs:/gpfs,/scratch:/scratch"

source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate ims-mosaic

MANIFEST="$1"
BATCH_RUNNER="$2"

if [[ "$MANIFEST" != /* || "$BATCH_RUNNER" != /* ]]; then
    echo "Manifest and Python runner paths must be absolute" >&2
    exit 2
fi
if [[ ! -f "$MANIFEST" || ! -f "$BATCH_RUNNER" ]]; then
    echo "Manifest or Python runner does not exist: $MANIFEST ; $BATCH_RUNNER" >&2
    exit 2
fi

echo "IMS field array task: ${SLURM_ARRAY_TASK_ID}"
echo "Node: $(hostname)"
echo "Manifest: ${MANIFEST}"
echo "Batch runner: ${BATCH_RUNNER}"
echo "Start: $(date)"
python "$BATCH_RUNNER" run --manifest "$MANIFEST" --index "$SLURM_ARRAY_TASK_ID" --resume
echo "Finish: $(date)"
