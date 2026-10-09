# IMS Mosaic Pipeline: Export → BaSiC → BigStitcher → Maximum Projection

## Purpose

This script implements an end-to-end microscopy mosaic pipeline for Imaris `.ims` acquisition fields.

Its documented workflow is:

```text
IMS acquisition fields
    |
    v
native-resolution TIFF stack export
    |
    v
default trailing all-zero XY padding removal
    |
    v
BaSiC illumination correction
    |
    v
BigStitcher spatial registration
    |
    v
BigStitcher overlap-based intensity adjustment
    |
    v
32-bit fused TIFF stack
    |
    v
2-D maximum-intensity projection
```

The script is designed around acquisition files named:

```text
<sample>_F00.ims
<sample>_F01.ims
...
<sample>_F99.ims
```

with exactly two digits after `F`.

A typical four-field run is:

```bash
python ims_mosaic_pipeline_illumination_bigstitcher.py /data/experiment \
    --fiji-sif /path/to/fiji_latest_bigstitcher.sif
```

---

## Main dependencies

Python-side dependencies:

```text
Python 3.10+
h5py
hdf5plugin
numpy
tifffile
```

Additional dependency required for the illumination-correction stage:

```text
basicpy
```

Install example:

```bash
pip install h5py hdf5plugin numpy tifffile basicpy
```

External dependencies:

```text
Fiji
BigStitcher
Singularity or Apptainer
```

The `--fiji-sif` argument must point to a Singularity/Apptainer image containing Fiji with BigStitcher installed.

The script explicitly masks:

```text
/opt/Fiji.app/plugins/SPIM_Registration-5.0.26.jar
```

at runtime using an empty JAR because the known working configuration uses BigStitcher without that legacy plugin.

---

# Input IMS files

## Accepted source argument

The positional argument:

```text
path
```

may be either:

```text
one .ims file
```

or:

```text
a directory containing .ims files
```

The script only searches the supplied directory itself.

It does **not** recursively search subdirectories.

---

## Required filename format

The discovery regex is:

```text
^(?P<sample>.+)_F(?P<field>\d{2})\.ims$
```

Therefore valid examples include:

```text
Exp17_151-2_Set2_DAPI_GFP_LHX6_PV_F00.ims
Exp17_151-2_Set2_DAPI_GFP_LHX6_PV_F01.ims
Exp17_151-2_Set2_DAPI_GFP_LHX6_PV_F02.ims
Exp17_151-2_Set2_DAPI_GFP_LHX6_PV_F03.ims
```

Invalid examples include:

```text
sample_F0.ims
sample_F1.ims
sample_F000.ims
sample_processed.ims
```

Files that do not match are skipped with a warning.

---

## Acquisition grouping

All files sharing the same text before:

```text
_FNN.ims
```

are grouped into one acquisition.

Example:

```text
Exp17_151-2_Set2_DAPI_GFP_LHX6_PV_F00.ims
Exp17_151-2_Set2_DAPI_GFP_LHX6_PV_F01.ims
Exp17_151-2_Set2_DAPI_GFP_LHX6_PV_F02.ims
Exp17_151-2_Set2_DAPI_GFP_LHX6_PV_F03.ims
```

are grouped under sample name:

```text
Exp17_151-2_Set2_DAPI_GFP_LHX6_PV
```

---

# IMS structural assumptions

The script expects an Imaris HDF5 structure containing:

```text
DataSet
DataSetInfo/Image
DataSetInfo/Channel N
```

Within `DataSet`, numeric groups are sorted naturally by their numeric suffix.

The script selects:

```text
first resolution level
single timepoint
all channels
```

and requires exactly one timepoint.

If an acquisition contains more than one timepoint, it raises an error.

---

# Channel names

Channel names are read from:

```text
DataSetInfo/Channel <index>
```

attribute:

```text
Name
```

The script converts bytes/NumPy string forms to text and replaces filesystem-invalid characters:

```text
\ / : * ? " < > |
```

with underscores.

If the channel name cannot be read, fallback names are:

```text
Channel0
Channel1
Channel2
...
```

Examples of expected microscope channel names might include:

```text
Confocal - Blue
Confocal - Green
Confocal - Red
Confocal - Far Red
```

---

# Physical voxel size

The script normally reads voxel calibration directly from IMS metadata.

It uses:

```text
DataSetInfo/Image
```

attributes:

```text
ExtMin0 / ExtMax0
ExtMin1 / ExtMax1
ExtMin2 / ExtMax2
X
Y
Z
Unit
```

For each axis:

```text
voxel size =
    (ExtMax - ExtMin)
    / number of voxels
```

The result is converted to micrometers.

Supported units:

```text
um
µm
micrometer
micrometers
nm
mm
m
```

If the IMS `Unit` field is blank, the default fallback is:

```text
um
```

unless changed with:

```text
--ims-unit
```

Voxel sizes across tiles are required to match within:

```text
rtol = 1e-4
atol = 1e-9
```

to tolerate small rounded acquisition differences.

---

# Manual voxel-size override

Optional argument:

```bash
--voxel-size X_UM Y_UM Z_UM
```

Example:

```bash
--voxel-size 0.312843137 0.312836345 0.494117647
```

When provided, this overrides the automatically read IMS calibration for downstream BigStitcher dataset creation.

---

# Z-plane handling

The script determines the usable Z extent from the first channel/timepoint.

It scans backward from the final Z plane until it finds a plane containing at least one non-zero pixel.

Trailing completely empty Z planes are excluded.

If the reference stack contains no non-empty Z planes, the script raises:

```text
ValueError
```

---

# XY zero-padding trimming

By default, the script removes only **trailing all-zero rows and columns** from exported IMS stacks.

This is designed for IMS files whose logical image may be, for example:

```text
2040 × 1992
```

while the underlying HDF5 dataset may be stored in a padded:

```text
2048 × 2048
```

array.

The trimming routine examines the trailing 128-pixel slab and finds the last row and column containing any non-zero data.

It removes only:

```text
bottom all-zero rows
right all-zero columns
```

It does **not** crop dark interior image regions.

Disable this behavior with:

```bash
--no-trim-zero-padding
```

This is documented in the script as not recommended for the IMS data the pipeline was designed around.

---

# Output root

For a directory input:

```text
/data/experiment
```

the pipeline output root is:

```text
/data/experiment/experiment_IMS_to_TIFF/
```

The rule is:

```python
base / f"{base.name}_IMS_to_TIFF"
```

where `base` is the supplied directory, or the parent directory when the
positional input is a single `.ims` file.

For each acquisition/sample, the script creates:

```text
<output_root>/<sample_name>/
```

---

# Per-acquisition directory structure

For an acquisition named:

```text
Exp17_151-2_Set2_DAPI_GFP_LHX6_PV
```

the script creates:

```text
experiment_IMS_to_TIFF/
└── Exp17_151-2_Set2_DAPI_GFP_LHX6_PV/
    ├── 01_exported/
    ├── 02_basic_corrected/
    ├── 03_stitched/
    └── 04_max_projections/
```

These four directories correspond to the major processing stages.

---

# Stage 1 — IMS → TIFF export

Output directory:

```text
01_exported/
```

For each field and channel, the script writes:

```text
FNN_<channel>.tif
```

Example:

```text
F00_Confocal - Blue.tif
F01_Confocal - Blue.tif
F02_Confocal - Blue.tif
F03_Confocal - Blue.tif
```

### TIFF format

Exported stacks use:

```text
ImageJ-compatible TIFF
axes = ZYX
```

The script writes:

```python
tifffile.imwrite(
    destination,
    stack,
    imagej=True,
    metadata={"axes": "ZYX"}
)
```

These exported intermediates are deliberately **uncompressed** so `tifffile.memmap()` can be used efficiently during BaSiC fitting.

### Exported stack dimensions

The expected shape is:

```text
(Z, Y, X)
```

where:

- trailing empty Z planes may be removed;
- trailing all-zero XY padding is removed by default.

---

# Stage-1 geometry audit files

For each field, the script writes:

```text
FNN_geometry.txt
```

Example:

```text
F00_geometry.txt
```

Contents:

```text
Source IMS array: Y=<source_y>, X=<source_x>
Exported logical image: Y=<valid_y>, X=<valid_x>
Trimmed trailing rows: <count>
Trimmed trailing columns: <count>
```

Example:

```text
Source IMS array: Y=2048, X=2048
Exported logical image: Y=1992, X=2040
Trimmed trailing rows: 56
Trimmed trailing columns: 8
```

This makes the automatic zero-padding trim auditable.

---

# Stage 2 — BaSiC illumination correction

Output directory:

```text
02_basic_corrected/
```

The output filenames are the same as the exported input names:

```text
FNN_<channel>.tif
```

For each channel, one BaSiC model is estimated using sampled Z planes from **all tiles for that channel**.

Default:

```text
12 sampled Z planes per tile
```

configured with:

```bash
--basic-sample-planes 12
```

The sample planes are evenly spaced through each stack using `numpy.linspace()`.

### BaSiC model

The script creates:

```python
BaSiC(get_darkfield=False)
```

so darkfield fitting is disabled.

### Corrected output dtype

BaSiC transformation produces floating-point values internally.

The script then:

1. rounds values;
2. clips to the legal `uint16` range;
3. converts to:

```text
uint16
```

Corrected TIFFs are written as:

```text
ImageJ-compatible ZYX TIFFs
```

---

# Stage 3 — BigStitcher registration and fusion

Output directory:

```text
03_stitched/
```

Each channel receives its own filesystem-safe subdirectory.

For example, a channel:

```text
Confocal - Blue
```

may use a directory such as:

```text
Confocal_-_Blue/
```

Spaces and other non-safe characters are replaced for the job directory.

The human-readable original channel name is retained for final filenames.

---

# BigStitcher input aliases

Inside each channel job directory, the BaSiC-corrected TIFFs are presented to BigStitcher as:

```text
F00.tif
F01.tif
F02.tif
F03.tif
```

These are symlinks when supported, otherwise physical copies.

This prevents broad filename patterns from accidentally loading duplicate/old versions.

---

# BigStitcher acquisition geometry

Default grid:

```text
2 × 2
```

configured with:

```bash
--grid 2 2
```

The macro uses:

```text
Snake: Up & Right
```

The script documentation states the known four-tile geometry as:

```text
F00 = bottom-left
F01 = top-left
F02 = top-right
F03 = bottom-right
```

Default overlap:

```text
10%
```

configured with:

```bash
--overlap 10
```

---

# BigStitcher processing stages

The generated ImageJ macro performs:

1. **Define Multi-View Dataset**
2. **Pairwise shifts with Phase Correlation**
3. **Filter pairwise links**
4. **Two-round global optimization**
5. **Overlap-based intensity adjustment**
6. **32-bit float image fusion**

---

# Pairwise phase correlation

Default spatial downsampling used in the macro:

```text
4 in X
4 in Y
4 in Z
```

Pairwise links are filtered by correlation.

Default minimum:

```text
r >= 0.55
```

configured with:

```bash
--min-correlation 0.55
```

---

# Global optimization

The macro uses:

```text
Two-Round using Metadata to align unconnected Tiles
```

with:

```text
relative = 2.500
absolute = 3.500
```

and fixes:

```text
group 0-0
```

---

# BigStitcher overlap-based intensity adjustment

This occurs **after spatial registration** and immediately before fusion.

Defaults:

```text
downsampling = 32
max_inliers = 10000
affine_intensity = enabled
offset_only regularization = 0.5
unmodified-intensity regularization = 0.5
```

Command-line controls:

```bash
--intensity-downsampling 32
--intensity-max-inliers 10000
--intensity-offset-only 0.5
--intensity-unmodified 0.5
```

The adjustment is derived from registered tile-overlap pixels.

---

# Fusion settings

The macro performs fusion using:

```text
Avg, Blending
Linear Interpolation
blending_range = 1000
32-bit float
adjust_image_intensities
non-rigid registration disabled
preserve_original
```

Output mode:

```text
Save as compressed TIFF stacks
```

---

# BigStitcher audit files

Each channel job directory contains:

```text
run_bigstitcher.ijm
bigstitcher.log
dataset.xml
dataset/
```

The exact BigStitcher-managed dataset files beneath `dataset/` depend on the Fiji/BigStitcher build.

The script also creates a pipeline-level masking file:

```text
.empty_SPIM_Registration-5.0.26.jar
```

under `03_stitched/`.

---

# BigStitcher completion checks

The Python wrapper expects the BigStitcher log to contain markers indicating that:

```text
intensity adjustment was reached
intensity adjustment completed
fusion was reached
```

Specifically, it checks for:

```text
BIGSTITCHER_INTENSITY_ADJUSTMENT_START
BIGSTITCHER_INTENSITY_ADJUSTMENT_DONE
BIGSTITCHER_FUSION_START
```

If required markers are missing, the script raises a runtime error and points to:

```text
bigstitcher.log
```

for debugging.

## Current marker sequence

In the current script, the generated `bigstitcher_macro()` and the Python-side
completion checks are consistent.

Immediately before the intensity-adjustment command, the macro prints:

```text
BIGSTITCHER_INTENSITY_ADJUSTMENT_START
```

It then runs:

```text
Adjust Intensities
```

After that call returns, the macro prints:

```text
BIGSTITCHER_INTENSITY_ADJUSTMENT_DONE
BIGSTITCHER_FUSION_START
```

before starting `Image Fusion`.

The Python wrapper checks for those same three markers after Fiji exits. This
provides three useful failure boundaries:

```text
START missing
    -> the macro never reached intensity adjustment

START present but DONE missing
    -> intensity adjustment was entered but did not complete

DONE present but FUSION_START missing
    -> intensity adjustment completed but fusion was not reached
```

Thus, the earlier START-marker mismatch no longer applies to the current
implementation.


---

# Final stitched 3-D TIFF

The script searches first for BigStitcher files matching:

```text
Fused_fused_tp_*_ch_*.tif
```

and then more broadly:

```text
*Fused*.tif
```

The first matching fused file is copied to:

```text
<channel job directory>/<channel>_stitched.tif
```

Example:

```text
03_stitched/
└── Confocal_-_Blue/
    └── Confocal - Blue_stitched.tif
```

This is the final per-channel **3-D fused stack** used for maximum projection.

---

# Stage 4 — maximum-intensity projection

Output directory:

```text
04_max_projections/
```

For each fused channel stack, the script runs:

```python
image = tifffile.imread(source)
projection = np.max(image, axis=0)
```

Therefore:

```text
input:  (Z, Y, X)
output: (Y, X)
```

No normalization or dtype conversion is applied by the projection function itself.

Therefore, when the BigStitcher fused input is the expected 32-bit floating-point
stack, the maximum projection preserves that floating-point intensity
representation.

The resulting file is written with:

```python
tifffile.imwrite(
    destination,
    projection,
    imagej=True,
    metadata={"axes": "YX"},
    compression="zlib"
)
```

---

# Final maximum-projection filename

The naming rule is:

```text
<sample_name>_<channel>_stitched_max_projection.tif
```

Example:

```text
Exp17_151-2_Set2_DAPI_GFP_LHX6_PV_Confocal - Blue_stitched_max_projection.tif
```

These are written to:

```text
04_max_projections/
```

---

# Complete representative output tree

```text
<experiment>_IMS_to_TIFF/
└── <sample>/
    ├── 01_exported/
    │   ├── F00_geometry.txt
    │   ├── F01_geometry.txt
    │   ├── F02_geometry.txt
    │   ├── F03_geometry.txt
    │   ├── F00_Confocal - Blue.tif
    │   ├── F01_Confocal - Blue.tif
    │   ├── F02_Confocal - Blue.tif
    │   ├── F03_Confocal - Blue.tif
    │   └── ... other channels ...
    │
    ├── 02_basic_corrected/
    │   ├── F00_Confocal - Blue.tif
    │   ├── F01_Confocal - Blue.tif
    │   ├── F02_Confocal - Blue.tif
    │   ├── F03_Confocal - Blue.tif
    │   └── ... other channels ...
    │
    ├── 03_stitched/
    │   ├── .empty_SPIM_Registration-5.0.26.jar
    │   ├── Confocal_-_Blue/
    │   │   ├── F00.tif
    │   │   ├── F01.tif
    │   │   ├── F02.tif
    │   │   ├── F03.tif
    │   │   ├── run_bigstitcher.ijm
    │   │   ├── bigstitcher.log
    │   │   ├── dataset.xml
    │   │   ├── dataset/
    │   │   ├── Fused_fused_tp_0_ch_0.tif
    │   │   └── Confocal - Blue_stitched.tif
    │   └── ... one directory per channel ...
    │
    └── 04_max_projections/
        ├── <sample>_Confocal - Blue_stitched_max_projection.tif
        ├── <sample>_Confocal - Green_stitched_max_projection.tif
        ├── <sample>_Confocal - Red_stitched_max_projection.tif
        └── <sample>_Confocal - Far Red_stitched_max_projection.tif
```

Some BigStitcher dataset/intermediate filenames may differ depending on the installed BigStitcher build.

---

# Command-line arguments

## Positional input

```bash
path
```

An `.ims` file or directory containing IMS fields.

When one valid `_FNN.ims` file is supplied, all sibling fields with the same sample prefix are discovered automatically. This makes `Sample_F00.ims` a valid shorthand for its complete mosaic acquisition.

---

## `--fiji-sif`

Required.

```bash
--fiji-sif /path/to/fiji_latest_bigstitcher.sif
```

Must point to an existing SIF containing Fiji and BigStitcher.

---

## `--output-root`

Optional destination root. When omitted, the historical `<input folder>_IMS_to_TIFF` location is used. The complete IMS controller supplies this option so preprocessing and downstream analyses remain under one requested output directory.

---

## `--channel-names`

Optional channel aliases assigned by zero-based IMS channel order. Names must be unique ASCII identifiers beginning with a letter.

The complete four-channel workflow uses:

```bash
--channel-names Hoechst mNeonGreen BiVe3 PV
```

This affects output filenames only; it does not reorder channel data or change pixel values. Original Imaris names and wavelength metadata are recorded separately by `run_ims_pipeline.py`.

---

## `--grid`

Default:

```text
2 2
```

Example:

```bash
--grid 2 2
```

The number of discovered fields must equal:

```text
grid_x × grid_y
```

or the script raises an error.

---

## `--overlap`

Default:

```text
10.0
```

Example:

```bash
--overlap 10
```

Regular-grid overlap percentage used for both X and Y.

---

## `--min-correlation`

Default:

```text
0.55
```

Pairwise phase-correlation links below this `r` are removed before optimization.

---

## `--basic-sample-planes`

Default:

```text
12
```

Number of evenly sampled Z planes from each tile used to fit the BaSiC model.

---

## `--voxel-size`

Optional manual calibration override:

```bash
--voxel-size X_UM Y_UM Z_UM
```

Normally omitted so calibration is read from IMS metadata.

---

## `--ims-unit`

Choices:

```text
um
nm
mm
m
```

Default:

```text
um
```

Used only if IMS spatial `Unit` metadata is blank.

---

## `--intensity-downsampling`

Default:

```text
32
```

BigStitcher overlap-intensity-adjustment downsampling.

---

## `--intensity-max-inliers`

Default:

```text
10000
```

Maximum overlap pixels per image pair used by BigStitcher's intensity adjustment.

---

## `--intensity-offset-only`

Default:

```text
0.5
```

BigStitcher offset-only regularization weight.

---

## `--intensity-unmodified`

Default:

```text
0.5
```

BigStitcher unmodified-intensity regularization weight.

---

## `--overwrite`

Accepted for compatibility.

The current pipeline cleans previous pipeline-generated acquisition results by default, so this flag is usually unnecessary for a normal rerun.

It becomes most relevant when `--keep-previous` is also used: without
`--overwrite`, existing exported, corrected, stitched, or projection outputs
may be reused/skipped according to the stage-specific existence checks.

---

## `--keep-previous`

By default, before reprocessing an acquisition, the script removes the previous pipeline-generated result directory for that acquisition.

Use:

```bash
--keep-previous
```

to preserve the old result directory.

Source IMS files are outside that directory and are not removed.

---

## `--keep-input-stacks`

This is retained as an intent/documentation flag.

The script already retains exported input stacks by default, and the current
implementation does not branch on this flag. In other words, supplying
`--keep-input-stacks` does not presently change pipeline behavior.

---

## `--no-trim-zero-padding`

Disables removal of trailing all-zero XY padding.

Example:

```bash
--no-trim-zero-padding
```

The source comments state that keeping this padding is not recommended for the intended IMS files because it can create artificial hard tile boundaries during stitching.

---

# Example command

For a four-field 2×2 acquisition:

```bash
python ims_mosaic_pipeline_illumination_bigstitcher.py \
    "/nfs/turbo/umms-parent/path/to/acquisition" \
    --fiji-sif "/path/to/fiji_latest_bigstitcher.sif" \
    --grid 2 2 \
    --overlap 10 \
    --min-correlation 0.55 \
    --basic-sample-planes 12
```

If IMS unit metadata are correct, no `--voxel-size` or `--ims-unit` override is needed.

---

# Failure checks

The script explicitly stops when:

```text
input path does not exist
Fiji SIF does not exist
neither Singularity nor Apptainer is available
no _FNN.ims files are discovered
field count does not match the requested grid
IMS timepoint count is not exactly one
channel names/order differ between tiles
voxel sizes differ materially between fields
reference channel contains no usable Z planes
trailing XY slab contains no data
BaSiCPy is unavailable
BigStitcher exits non-zero
required BigStitcher stage markers are absent
BigStitcher produces no fused TIFF
```

---

# Important scientific and implementation assumptions

1. The input represents a **single-timepoint** tiled acquisition.
2. Field filenames encode tile identity as exactly two digits: `F00`–`F99`.
3. All tiles in one acquisition have the same channel names and order.
4. All tiles have effectively identical physical voxel calibration.
5. The first channel can be used to determine common Z and XY storage geometry.
6. Trailing zero padding is storage padding rather than biological image content.
7. The configured grid and `Snake: Up & Right` acquisition ordering correctly describe the microscope tile layout.
8. The default BigStitcher macro is specifically tuned around the successfully tested 2×2 workflow described in the source comments.
9. BaSiC is applied before registration.
10. BigStitcher overlap-based intensity adjustment is applied after spatial registration and before fusion.
11. Fusion is 32-bit float.
12. Maximum projection assumes the fused TIFF uses Z as axis 0.
13. The maximum-projection TIFF is written as `YX` with zlib compression.
14. The previous generated acquisition directory is deleted on rerun unless `--keep-previous` is supplied.

---

# Recommended files to retain for reproducibility

At minimum, retain:

```text
source .ims files
FNN_geometry.txt audit files
BaSiC-corrected TIFFs
run_bigstitcher.ijm
bigstitcher.log
dataset.xml
final <channel>_stitched.tif
final *_stitched_max_projection.tif
exact command used to run the pipeline
Fiji/BigStitcher SIF definition or version information
```

These together document both the image-processing settings and the final image products.

---

# Summary

The script converts a tiled Imaris acquisition into analysis-ready 2-D max projections while retaining auditable intermediate products.

```text
*_FNN.ims
   |
   v
logical ZYX TIFF export
   |
   v
default trailing zero-padding trim
   |    (disable only with --no-trim-zero-padding)
   v
BaSiC correction
   |
   v
BigStitcher phase correlation
   |
   v
global optimization
   |
   v
overlap-based intensity adjustment
   |
   v
32-bit fused ZYX TIFF
   |
   v
maximum across Z
   |
   v
YX max-projection TIFF
```
