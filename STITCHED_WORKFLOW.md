# Automated workflow from Imaris-stitched TIFFs

`run_stitched_pipeline.py` connects the existing analysis programs without changing their scientific calculations. It starts with one already stitched TIFF per channel, skips IMS export/BaSiC/BigStitcher, and runs segmentation through final colocalization in the repository's established order.

The controller substitutes only each script's top-level `USER SETTINGS` in a temporary copy. The checked-in source scripts, matching logic, metric calculations, and percentile rules remain the authoritative implementations.

## Before the first run

Create the two environments documented under `CondaEnvironment/`:

```bash
conda env create -f CondaEnvironment/ims-mosaic.yml
conda env create -f CondaEnvironment/cellprofiler-native.yml
```

The repository now contains two deliberately different CellProfiler files:

- `CellProfiler/MGEOPV.cpproj` is the binary project for interactive GUI work.
- `CellProfiler/MGEOPVFinal.cppipe` is the text pipeline for headless runs.

The `.cppipe` was exported byte-for-byte from the pipeline embedded in the former project-formatted file. Its 24 modules and settings were not edited. The controller rejects HDF5 project files supplied as `cellprofiler.pipeline` so this format error fails before analysis starts.

## Inputs

Supply one single-channel TIFF for every entry in `channels`. Supported layouts are:

- two-dimensional `YX`; or
- a three-dimensional page stack reported by `tifffile` as `ZYX`, `QYX`, or `IYX`.

The runner max-projects accepted 3-D stacks across the first axis. RGB, multichannel, time-series, and ambiguous hyperstacks are rejected. Every channel must have the same final XY pixel dimensions. The images must depict the same organoid and, when alignment is disabled, must already use the same pixel coordinate frame.

Copy and edit the example:

```bash
cp stitched_pipeline.example.json stitched_pipeline.json
```

Important configuration fields:

| Field | Meaning |
|---|---|
| `nuclear_channel` | Channel used as the registration and persistent cell-identity reference |
| `channels` | Channel name to stitched TIFF path; relative paths are resolved from the JSON file |
| `output_dir` | New/empty analysis directory; relative paths are resolved from the JSON file |
| `pixel_size_um.x/y` | Actual exported-image calibration, not a universal default |
| `cellprofiler.pipeline` | Text `.cppipe`; relative paths are resolved from the repository root |
| `alignment.enabled` | `true` to run initial segmentation and CASTalign; `false` only for a shared coordinate frame |
| `alignment.extra_args` | Optional CASTalign tuning flags; runner-managed paths/calibration cannot be overridden |
| `pairs` | Ordered marker pairs; A-inside-B overlap is directional |
| `thresholds` | Explicit overlap fractions in `[0, 1]` |
| `filters.intensity/shape` | Existing percentile rules; percentile `0` disables that metric |
| `reports` | Enable or disable intensity, shape, and colocalization QC outputs |

Channel names must start with a letter and contain only ASCII letters, digits, and underscores.

## Validate, run, and resume

From the repository root:

```bash
conda activate ims-mosaic
python run_stitched_pipeline.py --config stitched_pipeline.json --dry-run
python run_stitched_pipeline.py --config stitched_pipeline.json
```

The dry run parses the complete configuration, verifies the `.cppipe`, reads every TIFF's metadata, checks axes and shared XY dimensions, and prints the exact channel/pair plan. It does not create the output directory or test installed environments.

Execution additionally checks the active environment's Python dependencies and verifies that `cellprofiler` imports through the configured Conda environment. CellProfiler is invoked with `conda run`, so no manual environment switching is required mid-run.

After an interruption, resume only with unchanged inputs and settings:

```bash
python run_stitched_pipeline.py --config stitched_pipeline.json --resume
```

`--resume` skips a stage only when the manifest marks it complete and all expected outputs are still nonempty. It refuses reuse when the configuration, input size/mtime, CellProfiler pipeline, controller, or scientific source scripts changed. Use a fresh `output_dir` after changing any of those inputs.

## Ordered handoffs

```text
stitched per-channel TIFFs
  -> max projection when needed
  -> initial CellProfiler segmentation (only when CASTalign is enabled)
  -> each marker registered to the common nuclear reference
  -> CellProfiler segmentation in the final coordinate frame
  -> ROI reconstruction and pixel-coordinate export
  -> intensity QC, then intensity-only filtering
  -> shape QC on intensity survivors, then shape-only filtering
  -> each marker matched to nuclear ROIs
  -> nuclear-associated marker populations compared pairwise
  -> colocalization QC and all-nuclear-cell summaries
```

The final coordinate-frame segmentation—not the initial registration segmentation—is used downstream. Shape filtering consumes the intensity-filtered pixels. Marker/nuclear HDF5 files preserve the canonical nuclear `ObjectID` mapping needed by the final all-cell summary; pairwise `CellID` values remain local to an individual colocalization run.

When `alignment.enabled` is `false`, the redundant initial segmentation pass is skipped. The one final-coordinate-frame CellProfiler pass is still performed.

## Outputs and diagnostics

The output directory contains:

```text
00_max_projections/             # only generated for 3-D inputs
01_initial_segmentation/        # only when alignment is enabled
02_registration/                # only when alignment is enabled
03_aligned_segmentation/
04_roi_extraction/
05_intensity_report/            # when enabled
05_intensity_filtered/
06_shape_reports/               # when enabled
06_final_shape_filtered/
07_nuclear_colocalization/
08_marker_colocalization/
run.log
run_manifest.json
```

`run.log` contains every command and its combined terminal output. `run_manifest.json` records stage completion and provenance fingerprints.

CASTalign currently writes every moving image as `aligned_PV_max.tif`. For a differently named marker, the controller adds a correctly named relative symbolic link rather than copying a multi-gigabyte TIFF or changing the registration algorithm.

The CellProfiler pipeline contract is intentionally strict. It must produce nonempty:

```text
CellMask.tiff
MyExpt_FilterObjects.csv
MyExpt_FilterObjects2.csv
```

The controller stops immediately if a subprocess fails or one of its required handoff files is missing.

## Scientific decisions that remain manual

- Confirm the real pixel calibration from the Imaris export.
- Inspect segmentation masks and registration QC on one acquisition before a batch run.
- Select intensity/shape percentiles and directional overlap thresholds for the experiment. The zero-percentile example disables filtering; it is not an optimized biological setting.
- Confirm that the supplied CellProfiler pipeline is appropriate for every channel's signal and object scale.
- Review final counts and all-nuclear classifications for biological plausibility.

The code and dry-run checks have been validated without altering the analysis algorithms, but an end-to-end biological run still requires the user's real stitched TIFFs and the two Conda environments.

See `SCRIPT_HANDOFF_REFERENCE.md` for the exact role and contract of every existing script.
