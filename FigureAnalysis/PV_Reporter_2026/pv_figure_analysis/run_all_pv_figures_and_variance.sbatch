#!/bin/bash -l
#SBATCH --job-name=PV_Figures_Variance
#SBATCH --account=parent0
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=01:20:00
#SBATCH --output=pv_all_figures_%j.out
#SBATCH --error=pv_all_figures_%j.err
set -euo pipefail
conda activate ims-mosaic
export OPENBLAS_NUM_THREADS=2
export OMP_NUM_THREADS=2
cd "$HOME/pv_figure_analysis"
python -u run_all_pv_figures_and_variance.py \
  --batch-dir "$HOME/neuron_ims_batch50_maskfix" \
  --include-inference
