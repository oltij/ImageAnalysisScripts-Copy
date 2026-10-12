#!/bin/bash -l
#SBATCH --job-name=PV_publication_figures
#SBATCH --account=parent0
#SBATCH --partition=standard
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=00:45:00
#SBATCH --output=pv_publication_%j.out
#SBATCH --error=pv_publication_%j.err
set -euo pipefail
cd "$HOME/pv_figure_analysis"
conda activate ims-mosaic
export OPENBLAS_NUM_THREADS=2
export OMP_NUM_THREADS=2
python -u pv_figures.py --batch-dir "$HOME/neuron_ims_batch50_maskfix" --strict --redraw-only --no-normality
