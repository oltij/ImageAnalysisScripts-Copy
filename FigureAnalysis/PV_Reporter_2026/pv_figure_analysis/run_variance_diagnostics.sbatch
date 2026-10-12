#!/bin/bash -l
#SBATCH --job-name=PV_Variance_QC
#SBATCH --account=parent0
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=pv_variance_%j.out
#SBATCH --error=pv_variance_%j.err
set -euo pipefail
conda activate ims-mosaic
export OPENBLAS_NUM_THREADS=2
export OMP_NUM_THREADS=2
cd "$HOME/pv_figure_analysis"
python pv_variance_diagnostics.py --batch-dir "$HOME/neuron_ims_batch50_maskfix" --strict
