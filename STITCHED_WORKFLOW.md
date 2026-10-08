# Stitched-TIFF workflow: minimally invasive automation

This adds `run_stitched_pipeline.py` **without changing any analysis algorithm**. It starts from separately exported, **already Imaris-stitched** per-channel TIFFs. The original IMS/BaSiC/BigStitcher pipeline remains untouched and can be used in other projects. Existing analysis scripts are not imported or rewritten; a temporary copy of each script receives only its file-level `USER SETTINGS` assignments, and is then run as the same original program.

## Run

From the **copied** repository root, in the `ims-mosaic` environment:

```bash
conda activate ims-mosaic
cp stitched_pipeline.example.json stitched_pipeline.json
# Edit paths, thresholds, pixel sizes, and filters in stitched_pipeline.json.
python run_stitched_pipeline.py --config stitched_pipeline.json --dry-run
python run_stitched_pipeline.py --config stitched_pipeline.json
# After interruption, with UNCHANGED inputs/configuration:
python run_stitched_pipeline.py --config stitched_pipeline.json --resume
```

The `cellprofiler-native` environment must exist. Its Python and CellProfiler executable are invoked via `conda run` from the controller, so manual switching is not required. The CellProfiler `.cppipe` is used unchanged. Its expected outputs are `CellMask.tiff`, `MyExpt_FilterObjects2.csv`, and `MyExpt_FilterObjects.csv`. Check that the .cppipe actually produces these for your images: channel identification, LoadImages behavior, object/mask naming, and output module paths are **not** made generic by the controller.

## Inputs and scientific controls

* Supply **one channel per TIFF**, preferably two-dimensional grayscale max projections. For 3D ZYX TIFFs, the runner computes a **maximum intensity Z projection**; multichannel hyperstacks are deliberately rejected. Imaris may export ambiguous axes: inspect your TIFF metadata before use.
* The same organoid must appear in each channel, with matching XY image dimensions before CASTalign (its original script explicitly requires this); the controller does **not** resize/crop images or reinterpret stage coordinates.
* Choose the real nuclear channel as `nuclear_channel`. Other channel names use letters/digits/underscores, and must match how your data and CellProfiler pipeline are intended to be handled.
* Set `pixel_size_um.x/y` to the actual acquisition/export calibration; the example values are historical instrument metadata, **not universal defaults**.
* Choose overlap thresholds deliberately. `marker_nuclear` applies **marker A inside nucleus B**; `marker_marker` applies **first named marker A inside second named marker B**. Both are directional and use the original stable-matching code. The example 0.25 values reproduce one historical code setting, not validated universal choices.
* A percentile of zero disables the corresponding filter. The example intentionally disables all gates until you choose thresholds based on QC. Never assume these default values complete biological quality control.
* Set `alignment.enabled` to `false` only if your channels already share the **same pixel coordinate frame**. The pre/post-alignment CellProfiler steps are still run separately for consistent downstream inputs.
* The existing CASTalign script hard-codes `aligned_PV_max.tif` for its moving output even when the channel is LHX6/GFP. The runner creates a correctly named symlink to that file, without editing the registration algorithm. Its `aligned_Hoechst_max.tif` reference output is not used.
* All original scripts, including percentile filtering and colocalization, retain existing scientific logic and file naming.

## Processing / handoffs

```text
Already Imaris-stitched per-channel TIFF (YX or recognized ZYX)
  -> [if ZYX] independent max projection (same XY canvas)
  -> initial CellProfiler / channel (centroids used for alignment only)
  -> CASTalign marker -> common nuclear reference / channel (optional)
  -> aligned CellProfiler resegmentation / channel
  -> extract_rois.py (CellMask.tiff + object CSV + organoid CSV)
  -> intensity-report QC / all channels
  -> FinalFilter/filter.py PASS 1: intensity metrics enabled, shape OFF
  -> shape/area/circularity/eccentricity/solidity reports on PASS 1 survivors
  -> FinalFilter/filter.py PASS 2: shape metrics enabled, intensity OFF
  -> colocalizationdapiscript1.py, marker A vs nucleus B / marker
  -> extract `tables/filtered_object_a_pixels` from each nuclear HDF5
  -> colocalizationdapiscript1.py, nuclear-associated marker A vs B / pair
  -> colocalizationdapiscript2.py (visualization, all-cell summary)
```

The source `CellID` is **local to each colocalization run**. Canonical cross-marker identities use matched **nuclear ObjectID** from the marker-vs-nuclear HDF5s. The runner supplies both HDF5s to the final visualizer so it can construct the all-nuclear-cell table; it never equates unrelated run-local CellIDs. It only compares pairs specified by `pairs` (default: all unordered marker pairs). Note **pair order changes the directional overlap question**.

## Generated directory layout

`00_max_projections/`, `01_initial_segmentation/`, `02_registration/`, `03_aligned_segmentation/`, `04_roi_extraction/`, `05_intensity_report/`, `05_intensity_filtered/`, `06_shape_reports/`, `06_final_shape_filtered/`, `07_nuclear_colocalization/`, `08_marker_colocalization/`, `run_manifest.json`.

Each stage has dedicated outputs and original filenames. The two successive `filter.py` invocations are separated into two directories; the second receives only first-pass filtered ROI pixels. Each marker-vs-nuclear HDF5 gets a neighboring `<marker>_nuclear_positive_cells.csv`. Pairwise analysis stores `analysis.h5` and a `visualizations/` directory.

## Failure and provenance

The runner aborts when a stage returns a nonzero exit code or omits its expected output, and keeps previously produced files. `--resume` skips only completed stages with existing nonempty expected outputs, and refuses to resume if the config or original TIFF size/mtime have changed. It is **not** a content-hash audit and should not be used to silently reuse results after editing source scripts or the `.cppipe`; use a new `output_dir` for changed analysis code or segmentation. It does not validate microscopy biology or registration accuracy.

### Known constraints to review

1. The supplied CellProfiler `.cppipe` may contain image-specific/load-image assumptions. Verify it produces the three expected files for each channel and inspect masks before large batch runs.
2. `extract_rois.py` reconstructs object IDs from mask intensity levels, assumes 2D label TIFFs with zero-valued background, and validates centroids against the object CSV; its existing validity checks are kept intact.
3. Filter, metric-report, and colocalization scripts retain their existing PDF, TIFF, and HDF5 memory/compute costs. Multiple channels on Great Lakes may need high-memory compute nodes.
4. The original visualizer is exceptionally large and may fail on sparse/empty populations. The runner preserves existing error behavior instead of changing scientific code. Disable reports temporarily only to isolate the issue, not to claim the QC is complete.
5. The original CellProfiler, marker/nuclear matching, and alignment algorithms **have not been end-to-end tested on your Imaris images**. Review the reports, coordinate consistency, and biological plausibility before using results.
6. The controller deliberately does not automate scientific filter cutoffs from report images. You must set them in JSON; `--resume` will refuse modified configuration, so after changing cutoffs start in a fresh output directory.

## Existing scripts and documentation

The source scripts and their original README files remain authoritative for algorithms, output schemas, plots and QC details: `CellProfiler/cellprofilerdriver.py`, `Alignment/CASTalign_two_channel_registration.py`, `ExtractingROIs/extract_rois.py`, the intensity and shape reports under `ROIFiltering/`, `ROIFiltering/FinalFilter/filter.py`, `Colocalization/colocalizationdapiscript1.py`, and `Colocalization/colocalizationdapiscript2.py`. The original root `README.md` describes the full eight-step process including optional upstream independent stitching; this document describes the shortened **stitched-image entry point** without altering the original route.
