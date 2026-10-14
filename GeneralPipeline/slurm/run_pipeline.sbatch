#!/bin/bash
#SBATCH --job-name=MicroscopyPipeline
#SBATCH --account=parent0
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=1-00:00:00
#SBATCH --output=logs/general_%A_%a.out
#SBATCH --error=logs/general_%A_%a.err

set -euo pipefail
: "${SLURM_ARRAY_TASK_ID:?Submit with an Slurm --array index}"
: "${1:?Pass the absolute GeneralPipeline manifest path}"
: "${2:?Pass the absolute GeneralPipeline/run_pipeline.py path}"

module load singularity/4.4.1
export SINGULARITY_BINDPATH="/gpfs:/gpfs,/scratch:/scratch"

source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate ims-mosaic

MANIFEST="$1"
PIPELINE_RUNNER="$2"
if [[ "$MANIFEST" != /* || "$PIPELINE_RUNNER" != /* ]]; then
    echo "Manifest and pipeline runner must be absolute paths" >&2
    exit 2
fi
if [[ ! -f "$MANIFEST" || ! -f "$PIPELINE_RUNNER" ]]; then
    echo "Manifest or pipeline runner does not exist: $MANIFEST ; $PIPELINE_RUNNER" >&2
    exit 2
fi

echo "General microscopy array task: ${SLURM_ARRAY_TASK_ID}"
echo "Node: $(hostname)"
echo "Manifest: ${MANIFEST}"
echo "Start: $(date)"
python "$PIPELINE_RUNNER" resume --manifest "$MANIFEST" --index "$SLURM_ARRAY_TASK_ID"
echo "Finish: $(date)"
