# Experiment 17 ventral: readable per-field stable-matched QC PDFs

`build_per_organoid_reports.py` is a **drop-in replacement** for the previous per-organoid report builder in `FigureAnalysis/mNeonGreenReporterAnalysis/`. It uses the adjacent existing `plot_igf_exp17_ventral_with_gfp_plots.py` and `igf_exp17_ventral_block_mapping.csv` unchanged.

## Core design

One PDF **per microscope field**, grouped by sample-section (`151-3/`, etc.). It is a read-only visualization/export layer. It does not rerun or modify alignment, CellProfiler, ROI filtering, stable matching, or prior analyses. New PDFs go to `mentor_field_plots/per_organoid_reports_quartiles/` and use `_organoid_report_quartiles.pdf` filenames. To regenerate any existing short PDF, pass `--overwrite`.

Marker-pair pages use the original unconstrained Stage-08 tree by default.
Use `--marker-pair-variant nucleus-constrained` to build a separate report
from Stage 09. That mode labels every marker pair `NUCLEUS-CONSTRAINED`, uses a
separate default output directory, and adds `_nucleus_constrained` to every PDF
and displayed-exemplar CSV name. It never mixes Stage-08 and Stage-09
marker-pair HDF5 files. Stage-07 marker-nuclear pages remain available in both
report variants and retain their original meaning.

### What changed in this replacement

This version is designed to be much easier to read:

- **Final stable matches only** for colocalization pages.
- **Quartiles** (Q1–Q4) instead of deciles.
- **Up to 10 displayed exemplars per quartile by default** (`--max-examples 10`).
- **All measurements are printed directly beside the exact ROI or ROI pair they describe.**
- **Only displayed exemplars** are written to the companion CSV files.
- If a quartile has more exemplars than fit comfortably, the script automatically uses **multiple pages**.

## Per-field pages

1. **Cover**: field, biological labels, anti-GFP status, final-match availability and cautions.
2. **Aligned full-field fluorescence**: Hoechst, GFP, LHX6, PV (when present), with Hoechst-defined tissue outline.
3. **Standalone final-stable-matched ROI quartiles**: for each selected channel × ranking metric, pages are organized by **Q1–Q4**. Each displayed row is one exact ROI and includes:
   - full-field location marker,
   - fluorescence crop,
   - exact ROI ID,
   - quartile membership,
   - ranking metric value,
   - area,
   - mean fluorescence,
   - integrated fluorescence,
   - circularity,
   - solidity,
   - eccentricity.
4. **Colocalization pages**: for each pair HDF5 with final stable matches, pages are organized by **Q1–Q4** using the upstream `A_Percent_Inside_B` ranking. Each displayed row is one exact final stable-matched pair and includes:
   - Object B fluorescence with blue outline,
   - Object A fluorescence with red outline,
   - A-red / B-green composite,
   - exact A/B ROI IDs,
   - directional overlap percentages,
   - full-ROI mean and integrated intensity for both ROIs,
   - overlap-only integrated intensities for both channels.
5. **Eight key IGF comparison charts** (PV/LHX6 and GFP) showing **the exact current field outlined in black**; use `--all-summary-plots` for all 20 charts.
6. **Methods, provenance and legends**.

Unmatched, best-partial-overlap, below-threshold, before/after and exhaustive pair-intensity pages are no longer assembled into the new PDF. `--append-existing-qc-pdfs` is deliberately rejected because it defeats the final-only scope.

## What “quartiles” means here

Quartiles are **equal-count rank quartiles**, sorted ascending by the selected ROI feature (for standalone pages) or by the upstream **`A_Percent_Inside_B`** (for matched pair pages). Within each quartile, the script chooses up to the requested number of real exemplars spaced across the rank order so the displayed set spans that quartile. The **actual observed measurement range and number of candidates within that quartile** are printed on the page.

Thus Q1–Q4 are **not** absolute overlap ranges 0–25%, 25–50%, etc. Upstream match thresholds remain unchanged.

For standalone quartiles, the allowed pool consists of **final filtered ROI IDs** that also occur in at least one existing final stable-match pair. This is only a display filter; it does not redefine any biological class.

## Integrated fluorescence intensity

All displayed intensities are computed using the **actual original aligned TIFF pixels**:

- for single-channel pages: the aligned channel image,
- for matched-pair pages: the original embedded HDF5 fluorescence images used by the upstream colocalization run.

Full-ROI integrated intensity is:

`sum(image[y, x] for all ROI pixels)`

Mean fluorescence is the corresponding raw average. Overlap-only intensity is the raw sum on the intersection of the two displayed partner masks. Display brightness scaling is never used for these measurements.

These values are **uncorrected** for camera/background, gain, or anti-GFP amplification and should not be interpreted as molecule counts or absolute normalized expression. The supplied ROI metrics CSV determines area/shape annotations. The numbers reflect only the **displayed** exemplars.

## Sample selection

Same as the existing `plot_igf_exp17_ventral_with_gfp_plots.py`: Experiment 17 ventral blocks 149 (CL32 No IGF), 150 (CL32 +IGF), 151 (H1 +IGF), 152 (CL23 +IGF), 153 (CL23 No IGF), 157 (H1 No IGF); `-2` anti-GFP stained and `-3` no anti-GFP; verified H1 selection and exclusion of failed 12/13.

## Install / run

The files should live together:

```text
FigureAnalysis/mNeonGreenReporterAnalysis/
  build_per_organoid_reports.py
  plot_igf_exp17_ventral_with_gfp_plots.py
  igf_exp17_ventral_block_mapping.csv
  README_per_organoid_reports.md
  tests/test_per_organoid_reports_quartiles.py
```

Python 3.10+ packages: `numpy`, `pandas`, `scipy`, `tifffile`, `matplotlib`, `h5py`. Running the synthetic test additionally requires `pytest` and PyMuPDF (imported as `fitz`). Run from the existing `ims-mosaic` environment on Great Lakes.

### Pilot: one real field

```bash
conda activate ims-mosaic
REPO="$HOME/ImageAnalysisScripts-Copy"
W1="/scratch/parent_root/parent0/$USER/ims_folder_batches/igf_94_production_v2"
W2="/scratch/parent_root/parent0/$USER/igf_organoids/production_94_remaining_cleanup_20261016"
python "$REPO/FigureAnalysis/mNeonGreenReporterAnalysis/build_per_organoid_reports.py" \
  --h1-workspace "$W1" --remaining-workspace "$W2" \
  --only-sample 151-3 --max-fields 1 --max-examples 10
```

### Every eligible field

```bash
python "$REPO/FigureAnalysis/mNeonGreenReporterAnalysis/build_per_organoid_reports.py" \
  --h1-workspace "$W1" --remaining-workspace "$W2" --max-examples 10
```

### Separate nucleus-constrained marker-pair reports

```bash
python "$REPO/FigureAnalysis/mNeonGreenReporterAnalysis/build_per_organoid_reports.py" \
  --h1-workspace "$W1" --remaining-workspace "$W2" \
  --marker-pair-variant nucleus-constrained --max-examples 10
```

This requires
`analysis/09_marker_colocalization_nucleus_constrained/*/analysis.h5` for each
field being reported and fails explicitly if no Stage-09 marker-pair result is
present. The default remains `--marker-pair-variant unconstrained`.

Useful options:

- `--max-examples 10` (default): up to 10 displayed exemplars per quartile.
- `--marker-pair-variant unconstrained|nucleus-constrained`: select one
  marker-pair result tree; default is unconstrained.
- `--roi-metrics integrated_intensity,area`: adds standalone quartile pages for additional ranking metrics.
- `--all-summary-plots`: include all 20 summary charts.
- `--only-field "Field_1"`
- `--out /scratch/.../custom_folder`
- `--overwrite`
- `--skip-summary-plots`
- `--skip-single-rois`
- `--skip-colocalization`
- `--max-fields 1`

**Run substantial batches on a Great Lakes compute node via Slurm rather than on a login node.** Test one field first.

## CSV outputs and identity matching

For every PDF, two CSV files are saved in the same sample-section folder:

- `<field>_displayed_single_roi_quartiles.csv`: one row for each exact standalone ROI shown, including quartile label, display order, ROI ID and measurements.
- `<field>_displayed_stable_matched_pair_quartiles.csv`: one row for each exact A/B pair shown, including quartile label, display order, directional overlaps and fluorescence totals.

Only **displayed** objects occur in either CSV. The report folder also contains `report_manifest.csv` for the generated field reports.

## Tests

```bash
python -m py_compile FigureAnalysis/mNeonGreenReporterAnalysis/build_per_organoid_reports.py
python FigureAnalysis/mNeonGreenReporterAnalysis/build_per_organoid_reports.py --help
python -m pytest -q FigureAnalysis/mNeonGreenReporterAnalysis/tests/test_per_organoid_reports_quartiles.py
```

The tests use **synthetic** microscopy arrays and HDF5; they are not biological data.
