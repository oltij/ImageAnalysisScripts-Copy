# Experiment 17 ventral organoid: one PDF per microscopy field

`build_per_organoid_reports.py` is a **read-only report generator** for the existing
`FigureAnalysis/mNeonGreenReporterAnalysis` workflow. It uses the exact selection,
experimental metadata, and 20 quantitative measures implemented by sibling
`plot_igf_exp17_ventral_with_gfp_plots.py`. It also reads final ROI pixel/metric
CSVs and colocalization HDF5 files created by the committed general pipeline.

**It does not rerun or modify microscopy data, registration, CellProfiler,
segmentation, filtering, nuclear assignments, or stable matching.**

## Output layout

```
<out>/
  report_manifest.csv
  151-3/
    Exp17_151-3_Set3_DAPI_GFP_LHX6_PV_Field_1_organoid_report.pdf
    Exp17_151-3_Set3_DAPI_GFP_LHX6_PV_Field_1_colocalized_roi_integrated_intensity.csv
  151-2/
    ...
  152-2/
    ...
```

Each PDF is **one microscopy field**, not all fields from a block pooled together.
Fields are grouped under a sample-section folder (e.g. `151-3`) because one
sample-section may contain multiple separately acquired fields. The `report_manifest.csv`
records the exact input field and output paths.

## Input selection (identical to the existing reporter-analysis script)

| Block | Cell line | IGF condition |
|---|---|---|
| 149 | CL32 | No IGF |
| 150 | CL32 | +IGF |
| 151 | H1 | +IGF |
| 152 | CL23 | +IGF |
| 153 | CL23 | No IGF |
| 157 | H1 | No IGF |

All fields must be **Experiment 17, ventral**. `-2` sample sections have anti-GFP
staining; `-3` do not. Controls inferred from the inventory should be verified
against laboratory records. Field indices 12/13 are excluded; H1 uses only the
pre-existing verified H1 CSV. The source selector also validates that both
production manifests refer to the same tiles and that downstream colocalization
has completed. The active analysis directories are selected by the original
script (earlier fields in original workspace; remaining fields in the later
cleanup workspace).

## Data sources

Each job's `output_dir` is resolved from the production `manifest.json` and then
its `analysis/` directory is examined. The path contract is:

```
<analysis>/
  run_manifest.json
  00_max_projections/<CHANNEL>.tif                 # fallback source fluorescence
  02_registration/<CHANNEL>/aligned_<CHANNEL>_max.tif
  03_aligned_segmentation/Hoechst/OrganoidMask.tiff
  06_final_shape_filtered/<CHANNEL>/
      <CHANNEL>_ROI_pixels_FILTERED.csv
      <CHANNEL>_kept_ROI_metrics.csv
  07_nuclear_colocalization/<PAIR>/analysis.h5
  08_marker_colocalization/<PAIR>/analysis.h5
```

`<CHANNEL>` includes GFP, LHX6, PV, and the Hoechst reference when present.
The nuclear reference may be named `DAPI` or `Hoescht` in the mask directory.
Actual path selection for a channel's fluorescence prefers an aligned TIFF,
then a colocalization HDF5-embedded original registered image, and finally the
`00_max_projections` TIFF (labeled as a possible alignment fallback).

Every pair-specific `analysis.h5` is the **canonical saved analysis**, containing:
`metadata/json`, `images/object_a_original`, `images/object_b_original`,
`tables/object_a_pixels`, `tables/object_b_pixels`, `tables/stable_matching`, and
`tables/overlap_results`. All ROI identities and final matching statuses come
from these existing data; no matching is recomputed.

Existing QC overlays, when available, are located at
`<PAIR>/visualizations/<OBJECT>_ROI_before_filtering.tif` and corresponding
`_after_filtering.tif` paths. The script inserts the before/after images in the
report. With `--append-existing-qc-pdfs`, it can also append **verbatim** all
existing PDFs under the current field's `05_intensity_report/`,
`06_shape_reports/`, `07_nuclear_colocalization/`, and `08_marker_colocalization/`
folders. This appendix option can drastically increase page count and file size.

## PDF contents

1. **Cover/provenance:** field ID, cell line, condition, stain status, interpretation.
2. **Whole-fluorescence panels:** Hoechst, GFP, LHX6, PV in grayscale, with the
   Hoechst-defined tissue border in **black** when the mask exists. Each channel
   is independently display-stretched (1st–99th percentile) for visibility.
3. **Standalone ROI bins:** for GFP, LHX6 and PV, sort final retained ROIs by
   each requested metric, and divide the sorted IDs into **ten approximately
   equal-count rank groups** (0–10%, 10–20%, …, 90–100%). Each populated group
   has one whole-organoid fluorescence image overlaid with **yellow** outlines
   for all other ROIs and **blue** outlines for sampled ROIs, plus up to
   `--max-examples` individual ROI crops with grayscale fluorescence underneath.
   Source values and actual object IDs are labeled.
4. **Colocalization QC for every saved nuclear and marker pair:** original
   before/after overlays (when present); pair metadata, original threshold,
   and **absolute 10-percentage-point overlap groups** for final matched pairs,
   A-only objects with best overlapping B context, and B-only objects with best
   overlapping A context. Empty groups are omitted. The selected examples are
   displayed as **Object B grayscale with a blue contour / Object A grayscale
   with a red contour / red(A)-green(B) fluorescence composite with contours**.
5. **Colocalized integrated intensity:** per-pair scatter plot and descriptive
   statistics for every finally stable-matched pair. A CSV separately exports
   the **complete, one-row-per-matched-pair** source values; no example sampling
   affects intensity measurements.
6. **All twenty comparison figures:** 14 PV/LHX6 and 6 GFP figures, with the
   same underlying field values, conditions, cell-line panels, shared y-axis
   scale, and color-by-sample scheme as the original plotting script. GFP
   figures separately display anti-GFP `-2` and no-anti-GFP `-3` rows.
   **Only the exact current field's marker is made larger and outlined in black.**
7. **Methods/legends/provenance:** measurements, colors, source file paths,
   rank-decile versus absolute-overlap-bin definitions, and statistical caveats.
8. **Optional verbatim QC PDF appendix:** when requested, includes the existing
   individual ROI intensity/shape and colocalization QC reports as source pages.

## Colocalized integrated fluorescence intensity

For each pair in `stable_matching`, with its Object-A and Object-B ROI masks:

```
A_roi_integrated_intensity_native_sum = Σ original_A_image[y,x] over ROI A pixels
B_roi_integrated_intensity_native_sum = Σ original_B_image[y,x] over ROI B pixels
A_overlap_integrated_intensity_native_sum = Σ original_A_image[y,x] over ROI A ∩ ROI B
B_overlap_integrated_intensity_native_sum = Σ original_B_image[y,x] over ROI A ∩ ROI B
```

The export additionally includes object IDs, original image pixel counts,
measured overlap pixels, both directional overlap percentages, sample/field ID,
channel names, source HDF5, and ROI mean intensities. Numbers use the exact
image-native fluorescence intensities written to HDF5; **no background
subtraction, gain calibration, or photobleaching correction is applied here**.
Summed values are not molecular counts. They should not be compared as
absolute expression between antibody-amplified and non-amplified sections
without additional technical validation. A given cell/ROI can be present in
more than one pairwise HDF5, so do not add pair totals across analyses without
accounting for repeated objects.

## The two different “10%” definitions

- **Single-channel ROI reports:** 10% **rank bins** according to the selected
  metric, with roughly equal numbers of ROIs per bin. A bin labeled `70–80%`
  is the 70th–80th percentile of the metric, **not** 70–80% of tissue area.
- **Colocalization reports:** 10 **absolute directional overlap bins** in
  percentage points. A bin labeled `70–80%` means 70–80% of the Object A ROI
  pixels fall inside its selected Object B ROI (or B inside A for B-only QC).
  This is not the proportion of organoid area. Colocalization classifications
  are unchanged from the completed stable-matching step.

## Dependencies

Python 3.10+ and the existing `ims-mosaic` environment with `numpy`, `pandas`,
`matplotlib`, `scipy`, `tifffile`, and `h5py`; `pypdf` is only needed for
`--append-existing-qc-pdfs`. Running the synthetic test additionally requires
`pytest` and PyMuPDF (imported as `fitz`). The report requires the existing
sibling script:

`FigureAnalysis/mNeonGreenReporterAnalysis/plot_igf_exp17_ventral_with_gfp_plots.py`

and the companion `igf_exp17_ventral_block_mapping.csv`. There is no need to
edit the original plotting script or rerun any upstream pipeline.

## Run for one sample-section first (Great Lakes)

```bash
conda activate ims-mosaic
REPO="$HOME/ImageAnalysisScripts-Copy"  # set to your current checkout location
SCRIPT="$REPO/FigureAnalysis/mNeonGreenReporterAnalysis/build_per_organoid_reports.py"
W1="/scratch/parent_root/parent0/$USER/ims_folder_batches/igf_94_production_v2"
W2="/scratch/parent_root/parent0/$USER/igf_organoids/production_94_remaining_cleanup_20261016"

python "$SCRIPT" \
  --h1-workspace "$W1" \
  --remaining-workspace "$W2" \
  --only-sample 151-3 \
  --max-fields 1 \
  --max-examples 3 \
  --out "$W1/mentor_field_plots/per_organoid_reports"
```

This is a useful first-run check because even one field may generate many PDF
pages. Once its PDF looks correct, omit `--only-sample` and `--max-fields` to
process all eligible fields. If desired, use
`--roi-metrics area,integrated_intensity` for a shorter but still quantitative
ROI report; the default includes area, mean intensity, integrated intensity,
circularity, solidity, and eccentricity.

**Append the original completed ROI/colocalization PDFs, if desired:**

```bash
python "$SCRIPT" --h1-workspace "$W1" --remaining-workspace "$W2" \
  --only-sample 151-3 --max-fields 1 --overwrite \
  --append-existing-qc-pdfs
```

Only run this optional appendix mode if scratch has adequate space: it may copy
very large multipage QC PDFs into every field's report.

**Transfer one sample-section directory to a Mac** (Mac Terminal):

```bash
mkdir -p "$HOME/Desktop/IGF_151-3_report"
rsync -avhP \
  oltij@greatlakes.arc-ts.umich.edu:/scratch/parent_root/parent0/oltij/ims_folder_batches/igf_94_production_v2/mentor_field_plots/per_organoid_reports/151-3/ \
  "$HOME/Desktop/IGF_151-3_report/"
```

## CLI reference

```
--h1-workspace            Earlier original 94-field workspace (required)
--remaining-workspace     Rerun/cleanup workspace for later field outputs (required)
--sample-metadata         CSV; defaults to sibling mapping file
--h1-selection            Optional override of verified H1 fields CSV
--out                     Output root (default original workspace/mentor_field_plots/per_organoid_reports)
--only-sample             Restrict to one block-section, e.g. 151-3
--only-field              Restrict by acquisition-name substring
--max-fields              Limit fields for pilot runs
--max-examples            Up to 10 examples per populated 10% group (default 3)
--roi-metrics             Comma-separated ROI ranking metric keys
--skip-single-rois        Skip single-channel ROI decile section
--skip-colocalization     Skip colocalization images AND matched-intensity calculations
--skip-summary-plots      Skip 20 comparison figures
--append-existing-qc-pdfs Append original QC PDFs verbatim; pypdf required
--overwrite               Overwrite existing finished per-field report
```

## Testing and limitations

- `python -m py_compile build_per_organoid_reports.py`
- `python build_per_organoid_reports.py --help`
- `pytest -q tests/test_per_organoid_reports.py` (synthetic, no lab data)

The accompanying synthetic test verifies the HDF5 matching-table contract,
ROI overlays, intensity arithmetic, 10% bins, highlighted plot rendering,
and PDF extractability. **Synthetic images are not experimental results.**

The script has not been executed on a real Great Lakes field in this package;
actual field TIFF/HDF5 availability and PDF page count must be validated with
`--max-fields 1` before bulk execution. Missing optional standalone ROI assets
are explicitly labeled as unavailable. Invalid/missing required ROI identities
and malformed HDF5 matching tables cause a clear failure rather than guessing.

A report may contain hundreds of pages depending on the chosen ROI metrics,
number of colocalization pairs, and density of populated bins. This is an
intentional complete QC mode; reduce metrics, limit examples, or skip sections
for smaller review copies. Every report retains the same source-label
provenance as the original comparative analysis. Field dots in plots are nested
inside biological samples; no new independent-replicate p-values are reported.
