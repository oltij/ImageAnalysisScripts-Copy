# Publication figures V1/V2 + variance figures and statistical diagnostics

This archive is a **full** replacement package for `~/pv_figure_analysis`, not a patch. It includes the original complete `pv_figures.py` unchanged from the earlier BC-only/legends package, the full statistical diagnostics script with **new explicit sample-variance visualizations**, existing verifier/Slurm scripts, and a new **one-command controller**.

## Outputs and conditions

**Scientific measurements:** same preexisting 35-acquisition `per_acquisition_metrics.csv`; no segmentation, colocalization, stitching, measurement recomputation or CSV changes. Inputs resolved from `~/neuron_ims_batch50_maskfix/prepared.json` to the actual scratch directory. Strict validation of every acquisition's index, line, batch, and enhancer.

**Publication figures (12 figure sets, 36 files):** Version 1 (experiment-colored circles) and Version 2 (experiment-colored dots, shape reflects H1 circle / CL32 square / CL23 triangle). Both retain no marker outlines, black mean bars for B/C, 1.5×IQR Tukey whiskers with individual points, matching red experiment palette, editable SVG text, PDF, and 600-DPI PNG. Figures: A PV area, A supplementary density, B specificity/coverage, C specificity/coverage, combined ABC, and BC-only. Standalone A panels have experiment legends; the BC-only panel has experiment legend in V1 and both experiment and cell-line legends in V2.

**Variance/statistical figures (15 figure sets, 45 files):** all six metrics (A PV area; A supplementary cell density; BiVe3 specificity and coverage; BiVe4 specificity and coverage) are analyzed separately:

* `all_metrics_distributions` — original raw distributions (2×3 grid).
* `sample_sd_by_cell_line` and `sample_sd_by_experiment` — SD bars annotated with sample variance.
* `additive_model_residual_qq` and `residuals_vs_fitted` — experiment-adjusted residuals.
* **NEW** `sample_variance_by_cell_line` — direct sample variance (s²) bars for all six outcomes.
* **NEW** `sample_variance_by_experiment` — direct variance bars using experiment red shades.
* **NEW** `sample_variance_cell_line_by_experiment_heatmaps` — 3×3 matrices for each metric; no information is fabricated for empty cells or n<2.
* **NEW** `absolute_deviations_by_cell_line` — individual absolute deviations from each cell-line median, plus mean absolute deviation bar (descriptive component of Brown–Forsythe test).
* **NEW** `A_PV_area_variance_detail`, `A_supp_PV_density_variance_detail`, `B_BiVe3_specificity_variance_detail`, `B_BiVe3_coverage_variance_detail`, `C_BiVe4_specificity_variance_detail`, `C_BiVe4_coverage_variance_detail` — each outcome's own 3-part variance comparison (by cell line / by batch / cell line × batch).

**Tables (7 CSVs):** `descriptive_variances.csv`, `variance_matrix_cell_line_experiment.csv`, `variance_homogeneity_tests.csv`, `shapiro_normality.csv`, `additive_model_residual_diagnostics.csv`, `per_acquisition_residuals.csv`, `validated_sample_assignments.csv`. The `README_results.md` summarizes outcomes. Optional (`--include-inference`) extra exploratory Panel A omnibus tests and HC3 pairwise comparisons.

Each sample variance uses `ddof=1` with squared units. **If n<2, s² is undefined and reported as `N/A`, not zero.** Within each heatmap, dark shading indicates relatively larger variance *for that outcome*, not a shared scale across metrics. These are descriptive, not group variance estimates with uncertainty intervals. Brown–Forsythe and Fligner tests have limited power in small groups and cannot prove equal variance; residual model assumptions must be considered before choosing inferential tests.

**Biological caution:** Images/acquisitions are NOT automatically independent organoids. If multiple images belong to the same organoid, aggregate or fit a hierarchical model before treating inference results as biological replicates.

## Install

On Mac Terminal after downloading ZIP:

```bash
scp ~/Downloads/pv_figures_all_with_variance_visualizations.zip \
  oltij@greatlakes-xfer.arc-ts.umich.edu:~/
```

On Great Lakes:

```bash
conda activate ims-mosaic
cp ~/pv_figure_analysis/pv_figures.py \
  ~/pv_figure_analysis/pv_figures.py.backup_$(date +%Y%m%d_%H%M%S)
cd ~
unzip -oq pv_figures_all_with_variance_visualizations.zip
cd ~/pv_figure_analysis
```

## ONE COMMAND: generate everything

```bash
conda activate ims-mosaic
cd ~/pv_figure_analysis
python run_all_pv_figures_and_variance.py \
  --batch-dir ~/neuron_ims_batch50_maskfix \
  --include-inference
```

Optional: `--stats-only` to avoid regenerating publication figures, `--publication-only` to generate only publication figures, omit `--include-inference` to generate descriptive statistics, variance/normality diagnostics and charts but no formal Panel A exploratory p-value table.

Or submit all plots to Slurm:

```bash
sbatch run_all_pv_figures_and_variance.sbatch
```

## Data locations

The figures are automatically placed in:

```text
/gpfs/accounts/parent_root/parent0/oltij/Neuron_Rebuttal_IMS_Batch50_MaskFix_20261009/_PV_selected_figures/
```

The 45 statistical diagnostic figure files and tables are placed in its `variance_diagnostics/` subfolder. The controller invokes both verifiers after processing. The source CSV is read-only.

## Download results to your Mac

```bash
mkdir -p "$HOME/Desktop/PV_All_Publication_and_Variance"
rsync -avh --progress \
  --include='*/' \
  --include='*_publication_v1.svg' --include='*_publication_v1.pdf' --include='*_publication_v1.png' \
  --include='*_publication_v2.svg' --include='*_publication_v2.pdf' --include='*_publication_v2.png' \
  --include='variance_diagnostics/*.csv' \
  --include='variance_diagnostics/*.md' \
  --include='variance_diagnostics/*.svg' \
  --include='variance_diagnostics/*.pdf' \
  --include='variance_diagnostics/*.png' \
  --exclude='*' \
  'oltij@greatlakes-xfer.arc-ts.umich.edu:/gpfs/accounts/parent_root/parent0/oltij/Neuron_Rebuttal_IMS_Batch50_MaskFix_20261009/_PV_selected_figures/' \
  "$HOME/Desktop/PV_All_Publication_and_Variance/"
open "$HOME/Desktop/PV_All_Publication_and_Variance"
```

You can also download all statistical diagnostics separately via `rsync` (see `README_VARIANCE.md`).
