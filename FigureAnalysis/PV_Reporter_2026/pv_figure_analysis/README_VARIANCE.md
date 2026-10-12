# PV organoid Panel A / B / C — variance and assumption audit

A **complete, read-only** exploratory statistics workflow for the **same 35 selected acquisitions** used by `pv_figures.py`. This package includes the latest full plotting script `pv_figures.py` (unchanged from the BC-only plus legends package) **and** an extended `pv_variance_diagnostics.py`. It does **not** rerun microscopy segmentation, registration, matching, or calculate new biological measurements. It does **not** overwrite `per_acquisition_metrics.csv` or existing publication figures.

## What is assessed?

- A: PV-positive projected area (%) for all 35 acquisitions, pooled across BiVe3/BiVe4 as in Panel A.
- A supplementary: nucleus-associated PV-positive cells per mm².
- B: BiVe3 specificity (%) and coverage (%) separately (17 acquisitions).
- C: BiVe4 specificity (%) and coverage (%) separately (18 acquisitions).

**Sample variance** uses the unbiased degrees-of-freedom convention `ddof=1`. For each outcome the report includes sample size and missingness, mean/median, variance, standard deviation (SD), SEM, IQR, coefficient of variation (CV when meaningful), and min/max, grouped by **cell line**, **experiment ID**, **cell line × experiment ID**, and for A also **enhancer** and **cell line × enhancer**. Differences in variance are tested using Brown–Forsythe (Levene with median center) and Fligner–Killeen across cell lines and, separately, experimental batches. Shapiro–Wilk results are included by group, but small n limits their interpretation. Separate additive models (`metric ~ cell_line + experiment ID`) produce residual SD/variance, Shapiro diagnostics, Q–Q plots and residual-versus-fitted plots.

For each metric, the experiment colors match the manuscript: **17** `#F7A1A3`, **24** `#DF4C50`, **RA02** `#801C31`; marker shapes encode H1 circle / CL32 square / CL23 triangle. The statistical/QC figures do not replace your manuscript figures.

### Optional exploratory Panel A inference

Use `--include-inference` to also output ordinary one-way ANOVA, **Welch's ANOVA**, Kruskal–Wallis, experiment-adjusted additive linear-model joint F tests (conventional and HC3 heteroskedasticity-robust covariance) and three HC3 pairwise contrasts, Holm-adjusted **within each model**. A second adjusted-model sensitivity check includes enhancer (`cell_line + experiment + enhancer`). The tests are **alternatives/sensitivity checks, not a prespecified final analysis**, and their differing p-values should not be used to select a preferred significant test after viewing the outcome. Results assume independent biological experimental units. The robust F reference is approximate, particularly for small samples.

**Crucial:** If multiple acquisitions came from sections/fields of the same organoid, acquisition-level tests **pseudoreplicate**. Confirm organoid identity before interpreting any p-values; aggregate organoid-level values or use a hierarchical model as needed. Batch adjustment cannot fully fix confounding caused by missing cell-line×experiment combinations. Neither a nonsignificant Levene/Brown–Forsythe result nor a nonsignificant Shapiro result proves an assumption true.

## Install from Mac to Great Lakes

For installation of the newest **all-figures** package, follow `README_ALL.md` instead. The commands below describe running this statistics script independently.

On **Mac Terminal**, after downloading the ZIP:

```bash
scp ~/Downloads/pv_figures_all_with_variance_visualizations.zip \
  oltij@greatlakes-xfer.arc-ts.umich.edu:~/
```

On **Great Lakes**:

```bash
conda activate ims-mosaic
cp ~/pv_figure_analysis/pv_figures.py ~/pv_figure_analysis/pv_figures.py.backup_$(date +%Y%m%d_%H%M%S)
cd ~
unzip -oq pv_figures_all_with_variance_visualizations.zip
cd ~/pv_figure_analysis
```

No extra package installation should be required for `ims-mosaic`: the script uses `numpy`, `pandas`, `scipy`, `matplotlib` already used by your figure workflow.

## Generate the variance / assumption report (no image reanalysis)

```bash
conda activate ims-mosaic
cd ~/pv_figure_analysis
python pv_variance_diagnostics.py \
  --batch-dir ~/neuron_ims_batch50_maskfix \
  --strict
```

This reads the actual scratch location from `~/neuron_ims_batch50_maskfix/prepared.json` and then reads the previously extracted `per_acquisition_metrics.csv`. It validates block ID, cell line, enhancer, and the exact 35-image selection, failing instead of silently mismatching samples.

**Optional** if you want exploratory Panel A test results in addition to the variance/assumption report:

```bash
python pv_variance_diagnostics.py \
  --batch-dir ~/neuron_ims_batch50_maskfix \
  --strict --include-inference
```

Optionally override the CSV location:

```bash
python pv_variance_diagnostics.py --csv /full/path/per_acquisition_metrics.csv --strict
```

To submit through Slurm instead, from the same folder:

```bash
sbatch run_variance_diagnostics.sbatch
```

### Outputs

All results go into the NEW subdirectory:

`/gpfs/accounts/parent_root/parent0/oltij/Neuron_Rebuttal_IMS_Batch50_MaskFix_20261009/_PV_selected_figures/variance_diagnostics/`

- `README_results.md` — readable overall summary, variance/SD by cell line, Brown–Forsythe and Fligner results, residual diagnostics and caveats
- `variance_matrix_cell_line_experiment.csv` — all 54 combinations across six measurements, three cell lines, and three experiments, with missing or n=1 variances marked undefined
- `descriptive_variances.csv` — group-level **sample variance**, SD, IQR, CV, SEM, counts, etc.
- `variance_homogeneity_tests.csv` — Brown–Forsythe/Levene, Fligner and largest/smallest sample variance ratios
- `shapiro_normality.csv` — within-group and pooled normality screening
- `additive_model_residual_diagnostics.csv` — batch-adjusted residual variances, residual Shapiro and model DF
- `per_acquisition_residuals.csv` — fitted values and residuals for every acquisition and metric
- `validated_sample_assignments.csv` — checked block/line/enhancer metadata
- **Fifteen three-format statistical/QC figure sets:** the five existing sets (`all_metrics_distributions`, `sample_sd_by_cell_line`, `sample_sd_by_experiment`, `additive_model_residual_qq`, `residuals_vs_fitted`), four new variance overviews (`sample_variance_by_cell_line`, `sample_variance_by_experiment`, `sample_variance_cell_line_by_experiment_heatmaps`, `absolute_deviations_by_cell_line`), and six new per-measurement `*_variance_detail` sets; each SVG, PDF and 600-DPI PNG (45 files)
- If `--include-inference`: `A_exploratory_omnibus_tests.csv`, `A_exploratory_HC3_pairwise.csv`

### Verify

```bash
python verify_variance_outputs.py \
  --results /gpfs/accounts/parent_root/parent0/oltij/Neuron_Rebuttal_IMS_Batch50_MaskFix_20261009/_PV_selected_figures/variance_diagnostics
```

The verifier checks CSV columns, 6 metrics, 35 validated acquisitions and that each of the 45 figure files exists, including editable SVG text.

## Copy diagnostics back to your Mac

Run in your **Mac Terminal**:

```bash
mkdir -p "$HOME/Desktop/PV_Statistical_Diagnostics"
rsync -avh --progress \
  'oltij@greatlakes-xfer.arc-ts.umich.edu:/gpfs/accounts/parent_root/parent0/oltij/Neuron_Rebuttal_IMS_Batch50_MaskFix_20261009/_PV_selected_figures/variance_diagnostics/' \
  "$HOME/Desktop/PV_Statistical_Diagnostics/"
open "$HOME/Desktop/PV_Statistical_Diagnostics"
```

## Manuscript figures remain separately available

To regenerate the two current publication versions, as before:

```bash
python pv_figures.py \
  --batch-dir ~/neuron_ims_batch50_maskfix \
  --strict --redraw-only --no-normality
```

This package does not change that script's calculations or figure formatting.
