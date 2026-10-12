# PV Panel A Kruskal–Wallis — fixed for actual Great Lakes selected metrics

This version uses the **real** column names shown in `per_acquisition_metrics.csv` and fixes the previous error. It retains only Panel A hypothesis testing and Panel B/C descriptive statistics. Existing B/C figures are not changed.

## Endpoints

- **A density:** `pv_nuclear_positive_density_mm2` directly from selected metrics, plotted as cells/mm².
- **A nuclear-associated PV area:** counts unique `X,Y` pixels of full PV segmentation objects in `tables/filtered_object_a_pixels` from `source_pv_nuclear_h5`. Divides by `organoid_area_pixels`, multiplies by 100. **This excludes PV objects that were NOT Hoechst-associated**; it is NOT PV+ nuclei fraction. The PV-positive objects' full area is retained, *not only the pixels overlapping a nucleus*.
- **B/C:** `specificity_percent` and `coverage_percent`, subset by `enhancer` BiVe3 or BiVe4; descriptive statistics only.

## Statistical method

Kruskal–Wallis tests compare three cell lines separately for both Panel A endpoints. Pairwise tests are two-sided Mann–Whitney U tests with Holm correction **within each endpoint**. These nonparametric tests do NOT adjust for experimental batch, and the design is unbalanced (not every cell line present in every experiment). The CSV measurements are acquisition-level; please verify one acquisition equals one biologically independent organoid before inference.

Plots contain boxplots (median/IQR, 1.5×IQR whiskers), all individual data points colored by experiment, with cell-line shape markers, mean ± SD overlaid, overall KW H and p, plus three pairwise Holm-adjusted p brackets. Plot formats: editable SVG and PDF, 600-DPI PNG.

## Install from Mac

```bash
scp ~/Downloads/pv_panelA_kw_corrected.zip oltij@greatlakes-xfer.arc-ts.umich.edu:~/
```

## Install and submit from Great Lakes

```bash
cd ~
mkdir -p ~/pv_panelA_kw_corrected
unzip -oq ~/pv_panelA_kw_corrected.zip -d ~/pv_panelA_kw_corrected
cd ~/pv_panelA_kw_corrected
sbatch run_pv_kw.sbatch
```

Monitor:

```bash
squeue -u "$USER"
cat ~/pv_panelA_kw_corrected/pv_kw_final_*.out
cat ~/pv_panelA_kw_corrected/pv_kw_final_*.err
```

## Output on Turbo/scratch

`<stage_root>/_PV_selected_figures/panelA_kw/`

- `A_density_KruskalWallis.svg/.pdf/.png`
- `A_nuclear_associated_area_KruskalWallis.svg/.pdf/.png`
- `Panel_A_KruskalWallis_omnibus.csv`
- `Panel_A_pairwise_MannWhitney_Holm.csv`
- `Panel_A_descriptive_mean_SD.csv`
- `BC_descriptive_mean_SD.csv`
- `per_acquisition_metrics_with_nuclear_associated_PV_area.csv`
- `panelA_kw_results_and_BC_descriptives.md`

## Download from Mac

```bash
mkdir -p "$HOME/Desktop/PV_KW_Final"
rsync -avh --progress \
  'oltij@greatlakes-xfer.arc-ts.umich.edu:/gpfs/accounts/parent_root/parent0/oltij/Neuron_Rebuttal_IMS_Batch50_MaskFix_20261009/_PV_selected_figures/panelA_kw/' \
  "$HOME/Desktop/PV_KW_Final/"
```

## Integrity checks

The script checks each HDF5's matched cell count against `pv_nuclear_positive_cells` in the input CSV; also checks the nuclear-associated PV union area is no greater than the total PV union area or organoid tissue mask area. It will refuse to generate misleading statistics if any check fails. If your original PV segmentation includes outside-mask pixels, the area comparison can fail; then the appropriate remedy is intersection with the Hoechst organoid mask, not disabling validation.
