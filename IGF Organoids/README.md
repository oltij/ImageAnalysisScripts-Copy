# IGF organoid stitched-image workflow

This folder is a focused entry point for the IGF organoid workflow validated
on Great Lakes on October 13, 2026. It does not duplicate the analysis engine.
`run_igf_organoids.py` launches the repository-level
`run_stitched_pipeline.py` with the same Python interpreter, so stage ordering,
segmentation, filtering, colocalization, provenance, and resume checks remain
shared with the established workflow.

The input is one independently stitched single-channel TIFF, or stitched
max-projection TIFF, for each of DAPI, GFP, LHX6, and PV. This launcher does not
perform IMS conversion, field stitching, BaSiC, or BigStitcher processing.

## Reproducibility constraints

Every IGF configuration must enable alignment and contain exactly these
CASTalign arguments:

```text
--skip-phase-init --match-gate-um 25 --iterations 5 --max-points 1500
```

The launcher checks that profile before delegating. It does not enable affine
registration. Phase correlation remains available to other users of the root
pipeline, but it is disabled here because unconstrained phase correlation
previously produced spurious shifts of approximately 1,700 pixels for LHX6 and
PV.

The existing five-iteration fitting behavior is intentionally retained for
experimental reproducibility. In particular, the saved image transform uses
the transform fitted in the fifth iteration rather than an accumulated
composition of all five fits. Iterative residuals therefore describe the
iterative fitting state and should not be claimed to exactly describe the
single transform applied to the saved image.

Channels may have slightly different native XY dimensions while alignment is
enabled. Their pixels and CellProfiler centroid coordinates are not resized,
cropped, downsampled, or padded. With the IGF `--skip-phase-init` profile, phase
correlation is not called at all. CASTalign applies its existing final transform
once to the original moving image and writes it on the fixed DAPI canvas.

## Environments

Create the two repository environments once:

```bash
conda env create -f CondaEnvironment/ims-mosaic.yml
conda env create -f CondaEnvironment/cellprofiler-native.yml
```

Run the launcher from `ims-mosaic`. The shared controller invokes
`cellprofiler-native` automatically with `conda run` during CellProfiler stages.

```bash
conda activate ims-mosaic
```

## Prepare a configuration

Copy the template to a per-organoid configuration:

```bash
cp "IGF Organoids/configs/igf_stitched.example.json" \
   "/path/to/configs/Organoid_01.igf.json"
```

Edit all four TIFF paths and choose a new, empty `output_dir`. Enter the actual
X/Y pixel calibration for these stitched images.

The exact October 13 scientific thresholds and filter rules are not stored in
this repository. The example deliberately contains invalid `null`/text
placeholders so it cannot silently run with invented settings. Copy
`pixel_size_um`, `thresholds`, and the complete `filters.intensity` and
`filters.shape` objects from the successful Great Lakes configuration:

```text
/scratch/parent_root/parent0/oltij/ims_folder_batches/olti_0623_test1/configs/
Exp17_157-2_Set2_DAPI_GFP_LHX6_PV_Field_1_3c6dd6a0ba.no_phase_gate25_iter5.stitched.json
```

Keep DAPI as `nuclear_channel`, retain GFP/LHX6/PV, all three marker pairs, the
checked-in `CellProfiler/MGEOPVFinal.cppipe`, and the exact alignment profile.
Do not reuse an output directory from an earlier code or configuration version.

## Dry run, execution, and resume

The launcher resolves the repository engine relative to its own location, so
these commands work from any current directory:

```bash
conda activate ims-mosaic

python "/path/to/ImageAnalysisScripts-Copy/IGF Organoids/run_igf_organoids.py" \
  --config "/path/to/configs/Organoid_01.igf.json" \
  --dry-run

python "/path/to/ImageAnalysisScripts-Copy/IGF Organoids/run_igf_organoids.py" \
  --config "/path/to/configs/Organoid_01.igf.json"
```

After an interruption, resume only with unchanged code, configuration, and
inputs:

```bash
python "/path/to/ImageAnalysisScripts-Copy/IGF Organoids/run_igf_organoids.py" \
  --config "/path/to/configs/Organoid_01.igf.json" \
  --resume
```

The manifest refuses a resume if the inputs, configuration, CellProfiler
pipeline, controller, or scientific scripts changed.

## Great Lakes Slurm example

Use a CPU node; this workflow does not request a GPU. A conservative first-run
submission is four CPUs, 64 GB RAM, and 12 hours. Adjust only after checking the
pilot's elapsed time and peak memory.

```bash
#!/bin/bash
#SBATCH --job-name=igf_organoid
#SBATCH --partition=standard
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/path/to/slurm/%x-%j.out
#SBATCH --error=/path/to/slurm/%x-%j.err

set -euo pipefail

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate ims-mosaic

REPO=/home/oltij/ImageAnalysisScripts-Copy
CONFIG=/path/to/configs/Organoid_01.igf.json

python "$REPO/IGF Organoids/run_igf_organoids.py" --config "$CONFIG"
```

Create the Slurm log directory before submission and submit with `sbatch`.

## Ordered stages

The shared controller performs:

1. Maximum projection only for inputs that are still Z stacks.
2. Initial CellProfiler segmentation for registration anchors.
3. Rigid CASTalign registration of GFP, LHX6, and PV to DAPI.
4. Final DAPI segmentation to establish the shared organoid mask.
5. Independent marker segmentation inside the shared DAPI mask.
6. Shared-mask equality, containment, and organoid-mask QC.
7. ROI reconstruction and pixel-coordinate export.
8. Intensity QC and intensity-only filtering.
9. Shape QC and shape-only filtering.
10. Marker/DAPI association and nuclear-positive exports.
11. All three requested marker/marker comparisons and their visualizations.

## Registration and downstream outputs

Each marker directory under `02_registration/` retains the existing CASTalign
outputs, including:

```text
matched_centroids_initial.csv
matched_centroids_rigid.csv
transform_rigid.txt
transform_final.txt
aligned_PV_max.tif
aligned_Hoechst_max.tif
aligned_two_channel.tif
registration_overlay.png
registration_qc.png
anchor_qc.png
anchor_residual_histogram.png
registration_summary.txt
```

The legacy `PV`/`Hoechst` filenames are preserved by the shared registration
script; the controller adds a marker-specific aligned-image symlink. All final
marker images are checked against the DAPI XY canvas before segmentation.

Other principal output directories are:

```text
01_initial_segmentation/
03_aligned_segmentation/
04_organoid_mask_qc/
04_roi_extraction/
05_intensity_report/
05_intensity_filtered/
06_shape_reports/
06_final_shape_filtered/
07_nuclear_colocalization/
08_marker_colocalization/
```

Review registration overlays, anchor CSVs/residual plots, shared-mask QC, and
final biological counts for the first organoid before launching a batch.

## Multiple organoids

Use one configuration and one new output directory per organoid. A Slurm array
may invoke the launcher once per configuration, but tasks must never share an
`output_dir`. Keep channel paths, calibration, thresholds, and filter settings
explicit in each configuration so the manifest records the complete analysis.
