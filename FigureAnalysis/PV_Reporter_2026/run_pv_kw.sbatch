#!/bin/bash
#SBATCH --job-name=PV_Kruskal_Wallis
#SBATCH --account=parent0
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=pv_kw_%j.out
#SBATCH --error=pv_kw_%j.err

set -euo pipefail

source ~/miniconda3/etc/profile.d/conda.sh
conda activate ims-mosaic

cd ~/pv_panelA_kw_pkg

echo "Starting PV statistical analysis"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Start time: $(date)"

python pv_panelA_kw_and_BC_descriptives.py \
  --batch-dir ~/neuron_ims_batch50_maskfix \
  --strict

echo "Analysis completed successfully"
echo "End time: $(date)"
