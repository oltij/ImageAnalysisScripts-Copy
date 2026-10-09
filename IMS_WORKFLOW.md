# Complete workflow directly from Imaris IMS files

`run_ims_pipeline.py` is the one-command entry point for the full pipeline. It reads `.ims` metadata, exports channel TIFFs, creates maximum projections, and then invokes `run_stitched_pipeline.py` for segmentation, registration, ROI filtering, nuclear association, and marker-pair analysis.

The scientific algorithms in the existing stitching and downstream scripts are reused unchanged. This controller supplies paths, fixed channel identities, and IMS calibration.

## Fixed channel assignment

Every input IMS file must contain exactly four channels in this ascending-wavelength order:

| IMS channel index | Biological identity | Pipeline name |
|---:|---|---|
| 0 | Hoechst | `Hoechst` |
| 1 | mNeonGreen | `mNeonGreen` |
| 2 | BiVe3 virus | `BiVe3` |
| 3 | PV | `PV` |

Assignment is based on the zero-based channel index, not the possibly inconsistent channel name stored by Imaris. The original Imaris channel name and any available wavelength attributes are preserved in `ims_metadata.json` for auditability.

Hoechst is automatically used as the nuclear reference. The shorter name `BiVe3` is used because downstream channel identifiers allow only letters, numbers, and underscores.

## Supported IMS layouts

### Mosaic fields

For independently acquired tiles, use filenames ending in a two-digit field number:

```text
Sample_F00.ims
Sample_F01.ims
Sample_F02.ims
Sample_F03.ims
```

With the default 2 × 2 grid, the controller reuses the existing processing order:

```text
IMS export -> trailing-zero-padding trim -> BaSiC
-> BigStitcher registration/fusion -> max projection
-> downstream analysis
```

Pass either the directory or any one field. When one `_FNN.ims` file is supplied, its sibling fields from the same sample are discovered automatically.

Mosaic mode requires the Fiji/BigStitcher SIF plus Singularity or Apptainer. The known working SIF can be built from `BigStitcher/fiji_latest_bigstitcher.def` as documented elsewhere in this repository.

### Already-stitched IMS

For one IMS file already stitched in Imaris, pass that file and set `input_mode` to `auto` or `stitched`. No Fiji SIF, BaSiC, or BigStitcher is used. The controller exports all four ZYX channel stacks and creates one YX maximum projection per channel before starting downstream analysis.

A directory of standalone, non-`_FNN` IMS files is interpreted as multiple already-stitched acquisitions and processed sequentially.

Do not mix `_FNN` mosaic fields and standalone IMS files in the same input directory when using `auto` mode.

## Configuration

From the repository root:

```bash
cp ims_pipeline.example.json ims_pipeline.json
```

Edit at least:

- `ims_input`: one IMS file or a directory containing acquisitions;
- `input_mode`: normally `auto`;
- `output_dir`: a new or empty directory;
- `fiji_sif`: required for mosaic fields, otherwise it may be `null`;
- filter percentiles and overlap thresholds appropriate for the experiment.

The XY and Z voxel sizes are read from `DataSetInfo/Image` extents and dimensions. `ims_unit` is used only when the IMS spatial-unit attribute is empty. The extracted X/Y calibration is passed into all downstream physical-area and shape calculations automatically.

The example's percentile value `0` disables that metric. It is a safe pass-through starting point, not a scientifically optimized threshold.

## Run

Activate the main environment:

```bash
conda activate ims-mosaic
```

Inspect all IMS files and the planned channel assignments without writing outputs:

```bash
python run_ims_pipeline.py --config ims_pipeline.json --dry-run
```

Then execute:

```bash
python run_ims_pipeline.py --config ims_pipeline.json
```

The controller invokes `cellprofiler-native` automatically during the CellProfiler stages.

After an interruption, with unchanged inputs and configuration:

```bash
python run_ims_pipeline.py --config ims_pipeline.json --resume
```

Resume is refused if the configuration, IMS file size/mtime, or workflow code changed. Use a new `output_dir` after changing scientific settings or source data.

## Preflight checks

Before output creation, the controller verifies:

- every acquisition has exactly four channels;
- all fields in a mosaic have compatible channel shapes and voxel calibration;
- mosaic field count matches `stitching.grid` and field numbers are contiguous from `F00`;
- every IMS file has exactly one timepoint at its native resolution;
- the fixed channel order is unchanged;
- the generated downstream configuration is valid.

When wavelength metadata are available for every channel, the metadata report records whether they are ascending. The declared index order remains the identity source because it is the acquisition convention specified for this project.

## Outputs

```text
<output_dir>/
├── ims_run.log
├── ims_run_manifest.json
├── generated_configs/
│   └── <sample>.stitched_pipeline.json
├── preprocessing/
│   └── <sample>/
│       ├── 01_exported/
│       ├── 02_basic_corrected/       # mosaic mode
│       ├── 03_stitched/              # mosaic mode
│       ├── 04_max_projections/
│       └── ims_metadata.json
└── analysis/
    └── <sample>/
        ├── 01_initial_segmentation/
        ├── 02_registration/
        ├── 03_aligned_segmentation/
        ├── 04_roi_extraction/
        ├── 05_intensity_filtered/
        ├── 06_final_shape_filtered/
        ├── 07_nuclear_colocalization/
        ├── 08_marker_colocalization/
        ├── run.log
        └── run_manifest.json
```

For already-stitched IMS input, the exported full Z stacks are retained under `01_exported/`. For mosaic input, all audit intermediates from the existing BaSiC/BigStitcher workflow are retained.

## First-acquisition review

Before a batch analysis, inspect:

1. `ims_metadata.json` for channel mapping and calibration;
2. the four TIFFs in `04_max_projections/`;
3. BigStitcher logs and fused images when using mosaic mode;
4. `03_aligned_segmentation/<channel>/CellMask.tiff`;
5. intensity and shape QC reports;
6. registration overlays and final cell-population counts.

The controller automates data movement and execution order; it does not choose biologically appropriate segmentation settings, filter cutoffs, or overlap thresholds.
