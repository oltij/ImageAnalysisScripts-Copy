# Experiment 17 ventral: concise per-field stable-matched QC PDFs

`build_per_organoid_reports.py` is a **drop-in replacement** for the previous per-organoid report builder in `FigureAnalysis/mNeonGreenReporterAnalysis/`. It uses the adjacent existing `plot_igf_exp17_ventral_with_gfp_plots.py` and `igf_exp17_ventral_block_mapping.csv` unchanged.

## Core design

One PDF **per microscope field**, grouped by sample-section (`151-3/`, etc.). It is a read-only visualization/export layer. It does not rerun or modify alignment, CellProfiler, ROI filtering, stable matching, or prior analyses. The original 400-page PDFs are left alone by default: new PDFs go to `mentor_field_plots/per_organoid_reports_quartiles/` and use `_organoid_report_quartiles.pdf` filenames. To regenerate any existing short PDF, pass `--overwrite`.

### Per-field pages

1. Cover: field, biological labels, anti-GFP status, final-match availability and cautions.
2. Aligned full-field fluorescence: Hoechst, GFP, LHX6, PV (when present), with Hoechst-defined tissue outline.
3. Standalone *final-stable-matched* ROI quartiles: one **landscape contact sheet per selected channel × ranking metric**. The default ranking is `Integrated intensity`, with one real ROI representing each nonempty rank quartile (Q1–Q4). Each row includes the field-wide position locator, the selected ROI's fluorescence crop and colored contour, the exact ROI ID, quartile range and number of eligible ROIs, and that SAME ROI's area, mean fluorescence, integrated fluorescence, circularity, solidity and eccentricity.
4. Colocalization: one **four-row contact sheet per existing pair HDF5 with final stable matches**. Only `stable_matching` pairs are eligible, and there is exactly one representative final pair in each nonempty directional overlap rank quartile. Each row includes **Object B fluorescence with blue outline**, **Object A fluorescence with red outline**, **two-channel A-red/B-green composite**, exact **A/B ObjectIDs**, directional overlap percentages, matched area in pixels, full-ROI mean/integrated intensity for BOTH ROIs, and overlap-only integrated intensities for BOTH channels. Every value relates to the two IDs in the SAME row.
5. Eight key IGF comparison charts (PV/LHX6 and GFP) showing **the exact current field outlined in black**; use `--all-summary-plots` for all 20 charts.
6. Methods, provenance and legends.

Unmatched, best-partial-overlap, below-threshold, before/after and exhaustive pair-intensity pages are no longer assembled into the new PDF. `--append-existing-qc-pdfs` is deliberately rejected because it defeats the final-only scope. Older upstream QC files remain untouched on scratch.

### What “quartiles” means here

Quartiles are **equal-count rank quartiles**, sorted ascending by the selected ROI feature (for standalone pages) or the upstream **`A_Percent_Inside_B`** (for matched pair pages). One representative per quartile is selected deterministically at the center of each sorted bin; tied values are ordered by ObjectID. The **actual observed measurement range and number of candidates within that quartile** are printed for every shown ROI/pair. A quartile can be empty if there are fewer than four eligible records. Thus Q1–Q4 are **not** absolute overlap ranges 0–25%, 25–50%, etc.; upstream numeric match thresholds are still unchanged.

For standalone quartiles, the allowed pool consists of **final filtered ROI IDs** that also occur in at least one existing final stable-match pair. The channel ROI IDs and recorded metrics come from `06_final_shape_filtered/<channel>/<channel>_ROI_pixels_FILTERED.csv` and `<channel>_kept_ROI_metrics.csv` respectively. If an ROI ID isn't in both files *and* a final stable match, it isn't eligible. This is a display filter, not new biological classification. Single-channel examples and pairwise examples are separate rank selections: **do not assume an ROI in one sheet is the same as a ROI in another sheet without checking the IDs**.

### Integrated fluorescence intensity

All displayed intensities are computed using the **actual original aligned TIFF pixels** (from the original colocalization HDF5 images for pair displays, and the aligned image for single-channel displays). Full-ROI integrated intensity is `sum(image[y, x] for all ROI pixels)`; mean is the corresponding raw average. Overlap-only intensity is the raw sum on the pixel intersection of the two displayed partner masks. Figures use display-only brightness stretching, *never* for these measurements. These values are uncorrected for camera/background, gain or anti-GFP amplification and are **not** comparable as molecule counts or absolute gain-controlled quantitative expression. The supplied ROI metrics CSV determines area/shape annotations. The numbers reflect only the **displayed** exemplars, not the distribution means of all objects.

### Sample selection

Same as the existing `plot_igf_exp17_ventral_with_gfp_plots.py`: Experiment 17 ventral blocks 149 (CL32 No IGF), 150 (CL32 +IGF), 151 (H1 +IGF), 152 (CL23 +IGF), 153 (CL23 No IGF), 157 (H1 No IGF); `-2` anti-GFP stained and `-3` no anti-GFP; verified H1 selection and exclusion of failed 12/13. Untreated block labels inferred from inventory require wet-lab confirmation. Do not treat multiple fields from the same sample block as biologically independent replicates.

## Install / run

The files should live together:

```text
FigureAnalysis/mNeonGreenReporterAnalysis/
  build_per_organoid_reports.py       # replace old version
  plot_igf_exp17_ventral_with_gfp_plots.py # keep existing
  igf_exp17_ventral_block_mapping.csv     # keep existing
  README_per_organoid_reports.md     # replace old version
  tests/test_per_organoid_reports.py # replace old version
```

Python 3.10+ packages: `numpy`, `pandas`, `scipy`, `tifffile`, `matplotlib`, `h5py`. Running the synthetic test additionally requires `pytest` and PyMuPDF (imported as `fitz`). Run from the existing `ims-mosaic` environment on Great Lakes. Original analysis files must be accessible (not uploaded into GitHub).

### Pilot: one real field

```bash
conda activate ims-mosaic
REPO="$HOME/ImageAnalysisScripts-Copy"
W1="/scratch/parent_root/parent0/$USER/ims_folder_batches/igf_94_production_v2"
W2="/scratch/parent_root/parent0/$USER/igf_organoids/production_94_remaining_cleanup_20261016"
python "$REPO/FigureAnalysis/mNeonGreenReporterAnalysis/build_per_organoid_reports.py" \
  --h1-workspace "$W1" --remaining-workspace "$W2" \
  --only-sample 151-3 --max-fields 1
```

A new folder is created at `$W1/mentor_field_plots/per_organoid_reports_quartiles/151-3/` with field-named PDF and two CSVs.

### Every eligible field

```bash
python "$REPO/FigureAnalysis/mNeonGreenReporterAnalysis/build_per_organoid_reports.py" \
  --h1-workspace "$W1" --remaining-workspace "$W2"
```

Additional options: `--roi-metrics integrated_intensity,area` (adds one whole **quartile contact sheet** for each additional ranking, but still only one exemplar per quartile); `--all-summary-plots` (20 charts instead of 8); `--only-field "Field_1"`; `--out /scratch/.../custom_folder`; `--overwrite`; `--skip-summary-plots`; `--skip-single-rois`; `--skip-colocalization`; `--max-fields 1`. The old `--compact` is accepted but no longer needed; the new quartile design is always concise. `--max-examples` is accepted for old commands, but exactly one displayed exemplar per quartile is enforced to keep the PDF readable.

**Run substantial batches on a Great Lakes compute node via Slurm rather than running them unattended on a login node.** Test one field first because this implementation has not been validated on real Great Lakes source files.

## CSV outputs and identity matching

For every PDF, two CSV files are saved **under the same sample-section folder**:

- `<field>_displayed_single_roi_quartiles.csv`: one row for each exact standalone ROI shown in a PDF row, including the ranking quartile, field index, channel, ROI ID and measurements.
- `<field>_displayed_stable_matched_pair_quartiles.csv`: one row for each exact A/B pair shown, with matched ROI IDs, both directional overlaps and both raw fluorescence totals (full ROIs and overlap pixels).

Only **displayed** objects occur in either CSV. The report folder also contains `report_manifest.csv` for all field reports generated in that invocation. Existing BigStitcher/raw microscopy files and prior QC are unchanged.

## Tests

From this folder or the repo root:

```bash
python -m py_compile FigureAnalysis/mNeonGreenReporterAnalysis/build_per_organoid_reports.py
python FigureAnalysis/mNeonGreenReporterAnalysis/build_per_organoid_reports.py --help
python -m pytest -q FigureAnalysis/mNeonGreenReporterAnalysis/tests/test_per_organoid_reports.py
```

The tests use **synthetic** microscopy arrays and HDF5; the example PDF used during development is **not biological data**. In deployment, always verify one real PDF and compare displayed ROI IDs, matched pair IDs, directional overlaps and integrated intensities to the upstream HDF5.
