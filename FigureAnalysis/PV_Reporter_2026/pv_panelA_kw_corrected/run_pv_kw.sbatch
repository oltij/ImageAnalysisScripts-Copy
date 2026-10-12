#!/bin/bash
#SBATCH --job-name=PV_KW_Final
#SBATCH --account=parent0
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=01:00:00
#SBATCH --output=pv_kw_final_%j.out
#SBATCH --error=pv_kw_final_%j.err
set -euo pipefail
source ~/miniconda3/etc/profile.d/conda.sh
conda activate ims-mosaic
cd ~/pv_panelA_kw_corrected
python -u pv_panelA_kw_and_BC_descriptives.py --batch-dir ~/neuron_ims_batch50_maskfix --strict
