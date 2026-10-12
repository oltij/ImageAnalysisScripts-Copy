# Full package: see README_ALL.md for one command to regenerate publication figures and all variance figures.

# PV figures — V1/V2, including BC-only and separate legends

This package **replaces only the plotting/export script**, preserving your already-extracted `per_acquisition_metrics.csv` and all original scratch analysis results. It builds on the previous manifest-validated code rather than an independent extraction pipeline.

## Output styles

- **Version 1:** all dots are circles, with no outlines; red color indicates experimental batch (**17** light red, **24** medium red, **RA02** dark red).
- **Version 2:** same color assignments, means, axes, and data; cell line is additionally encoded by **H1 circle**, **H1 CL32 square**, **H1 CL23 triangle**. The shape legend is neutral gray.
- **Area and PV cell-density:** boxplots with the conventional **Tukey 1.5×IQR whiskers**, median and quartiles, and all individual organoid dots (including outliers), no dot outlines.
- **B and C:** closely grouped specificity and coverage dots; each group has a **thick black horizontal bar at the arithmetic mean of the organoid-level percentages**.
- **NEW:** `BC_only_specificity_coverage_publication_v1/v2`: compact B/C side-by-side figure. V1 includes the experiment-color legend; V2 includes **both experiment-color and clone-shape legends**.
- **NEW:** standalone A area and supplementary PV density figures each include an experiment-color legend in **both versions**. Version 2 additionally includes the shape legend.
- No statistics or normality figures are produced with the recommended command.

Each of **six figures × two versions** is exported as an editable-text vector SVG, vector PDF, and **600-DPI PNG** = **36 files**. Matplotlib saves SVG text as actual `<text>` elements (`svg.fonttype=none`); Illustrator should have the selected typeface installed (DejaVu Sans by default). The SVG is scalable without intrinsic DPI constraints; embedded raster content, if any, is rendered at 300 DPI.

## 1. Mac → Great Lakes

Download the ZIP and run in your **Mac Terminal**:

```bash
scp ~/Downloads/pv_figures_bc_legends_full.zip \
  oltij@greatlakes-xfer.arc-ts.umich.edu:~/
```

## 2. Great Lakes installation

Run in the Great Lakes shell:

```bash
conda activate ims-mosaic
mkdir -p ~/pv_figure_analysis
if [ -f ~/pv_figure_analysis/pv_figures.py ]; then
  cp ~/pv_figure_analysis/pv_figures.py \
    ~/pv_figure_analysis/pv_figures.py.backup_$(date +%Y%m%d_%H%M%S)
fi
cd ~
unzip -oq pv_figures_bc_legends_full.zip
cd ~/pv_figure_analysis
```

## 3. Redraw from your existing CSV (fast; no CellProfiler run)

```bash
python pv_figures.py \
  --batch-dir ~/neuron_ims_batch50_maskfix \
  --strict --redraw-only --no-normality
```

The script obtains the **actual scratch root** from `~/neuron_ims_batch50_maskfix/prepared.json`, checks the 35 selected IMS indices, verifies block/experiment/clone/enhancer metadata, and reads `per_acquisition_metrics.csv` without altering it.

**Only if the CSV is not already present**, reconstruct it from completed scratch analyses (still no segmentation rerun):

```bash
python pv_figures.py \
  --batch-dir ~/neuron_ims_batch50_maskfix \
  --strict --no-normality
```

## 4. Verify outputs

```bash
python verify_figure_exports.py \
  --figures /gpfs/accounts/parent_root/parent0/oltij/Neuron_Rebuttal_IMS_Batch50_MaskFix_20261009/_PV_selected_figures
```

Expected: **12 figure sets / 36 files**, all SVGs containing editable text.

Figure stems are:

```
A_PV_area_publication_v1/v2
A_supp_PV_density_publication_v1/v2
B_BiVe3_specificity_coverage_publication_v1/v2
C_BiVe4_specificity_coverage_publication_v1/v2
ABC_PV_area_specificity_coverage_publication_v1/v2
BC_only_specificity_coverage_publication_v1/v2
```

Each has `.svg`, `.pdf`, `.png` extensions.

## 5. Great Lakes → Mac

Run in **Mac Terminal**, not the SSH session:

```bash
mkdir -p "$HOME/Desktop/PV_Figures_BC_Legends"
rsync -avh --progress \
  --include='*/' \
  --include='*_publication_v1.svg' --include='*_publication_v1.pdf' --include='*_publication_v1.png' \
  --include='*_publication_v2.svg' --include='*_publication_v2.pdf' --include='*_publication_v2.png' \
  --exclude='*' \
  'oltij@greatlakes-xfer.arc-ts.umich.edu:/gpfs/accounts/parent_root/parent0/oltij/Neuron_Rebuttal_IMS_Batch50_MaskFix_20261009/_PV_selected_figures/' \
  "$HOME/Desktop/PV_Figures_BC_Legends/"
open "$HOME/Desktop/PV_Figures_BC_Legends"
```

**Caution:** Each observation is one selected IMS acquisition. Confirm that acquisitions correspond to independent organoids before statistical comparisons. The script does not report p-values and does not alter any quantitative measurements.
