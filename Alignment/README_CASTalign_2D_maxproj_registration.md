# CASTalign 2-D Max-Projection Registration

This README documents:

```text
CASTalign_2D_maxproj_registration.py
```

uploaded as:

```text
Pasted text(20260915-172047).txt
```

This document is intentionally limited to this registration script.

# 2. Common coordinate conventions

These scripts use standard image-array and CellProfiler-style pixel coordinates, but the ordering is worth stating explicitly.

## Image arrays

NumPy/TIFF arrays are indexed as:

```text
image[Y, X]
```

For a 2-D image:

```text
shape = (height, width)
      = (Y, X)
```

Example:

```text
shape = (3840, 3840)
```

## ROI pixel-coordinate CSVs

Pixel-coordinate tables use:

```text
ROI,X,Y
```

Example:

```csv
ROI,X,Y
1,1834,1902
1,1835,1902
1,1836,1902
2,2050,1711
```

`X` is the horizontal image coordinate and `Y` is the vertical image coordinate.

## Current physical pixel size

The registration and colocalization scripts currently use:

```text
X = 0.312843137 µm/pixel
Y = 0.312836345 µm/pixel
```

Change these values if the analyzed image has a different physical scale.

---

# 4. CASTalign 2-D max-projection registration

Suggested/current script name from its usage documentation:

```text
CASTalign_2D_maxproj_registration.py
```

Uploaded source:

```text
Pasted text(20260915-172047).txt
```

---

## 4.1 Purpose

This is a **2-D max-projection registration** script.

It is designed to align a moving fluorescence channel to a fixed/reference fluorescence channel from the same organoid.

Example:

```text
fixed  = Hoechst max projection
moving = PV max projection
```

The workflow is:

1. load two 2-D max-projection TIFFs;
2. load object centroid coordinates for each channel;
3. optionally calculate a coarse phase-correlation XY shift;
4. globally match a subset of centroids using partial Hungarian assignment;
5. fit a CASTalign rigid transform;
6. iteratively rematch and reject residual outliers;
7. optionally test a planar affine transform;
8. resample the moving max projection into the fixed coordinate frame;
9. save aligned TIFFs, transforms, anchor tables, and QC.

### Scientific interpretation

The centroid pairs are **registration anchors**.

They are **not** claims that the marker-positive object and Hoechst object are the same biological cell.

The original 3-D stacks are not modified. Registration is calculated only on the 2-D max projections.

---

## 4.2 Dependencies

The source lists:

```text
castalign
tifffile
pandas
numpy
scipy
scikit-image
matplotlib
```

Example:

```bash
pip install castalign tifffile pandas numpy scipy scikit-image matplotlib
```

---

## 4.3 Command-line interface

This script is configured through command-line arguments.

### Required arguments

```text
--fixed-image
--moving-image
--fixed-csv
--moving-csv
--output
```

### Example

```bash
python CASTalign_2D_maxproj_registration.py \
    --fixed-image "/path/Hoechst_max.tif" \
    --moving-image "/path/PV_max.tif" \
    --fixed-csv "/path/Hoechst_objects.csv" \
    --moving-csv "/path/PV_objects.csv" \
    --output "/path/CASTalign_Hoechst_PV"
```

---

## 4.4 Input 1 — fixed/reference max-projection TIFF

### Argument

```text
--fixed-image
```

### Expected origin

A stitched/max-projected fluorescence image that defines the desired final coordinate frame.

Typical example:

```text
Hoechst / DAPI maximum-intensity projection
```

It may come from the IMS → TIFF → BigStitcher/max-projection workflow.

### Required format

```text
TIFF
2-D after numpy.squeeze()
```

Example:

```text
shape = (3840, 3840)
```

Any numeric TIFF dtype accepted by `tifffile.imread()` is acceptable, but the image must reduce to exactly two dimensions.

---

## 4.5 Input 2 — moving max-projection TIFF

### Argument

```text
--moving-image
```

### Expected origin

The fluorescence channel to transform into the fixed image's coordinate system.

Examples:

```text
PV max projection
LHX6 max projection
GFP max projection
```

### Required format

Same as the fixed image:

```text
TIFF
2-D after squeeze
```

### Image canvas

Fixed and moving images may have different native XY dimensions. Their pixels
and centroid coordinates remain in their original coordinate systems. The
final moving image is resampled once onto the fixed image's XY canvas.

If phase correlation is enabled for unequal images, the script temporarily
zero-pads their bottom and right edges to a common array size only for that
coarse initializer. `--skip-phase-init` uses a zero initial translation and
does not create padded arrays or call phase correlation.

---

## 4.6 Input 3 and 4 — fixed and moving centroid CSVs

### Arguments

```text
--fixed-csv
--moving-csv
```

### Expected origin

CellProfiler object-measurement CSVs for the exact max projections supplied as fixed and moving images.

Each valid centroid must fall within the corresponding image bounds.

### Default expected centroid column names

For X, the script tries:

```text
Location_Center_X
Location_Center_X_(px)
Location_Center_X_(pixels)
```

For Y:

```text
Location_Center_Y
Location_Center_Y_(px)
Location_Center_Y_(pixels)
```

It can also discover a unique fallback column containing `location_center` and ending in `_x` or `_y`.

### Important compatibility note with the ROI-export script

The ROI reconstruction script documented above expects a standard CellProfiler object CSV containing:

```text
AreaShape_Center_X
AreaShape_Center_Y
```

Those names are **not** among this registration script's default auto-detected names.

If your registration CSV uses the `AreaShape_...` names, pass the column overrides explicitly:

```bash
--fixed-x-col AreaShape_Center_X \
--fixed-y-col AreaShape_Center_Y \
--moving-x-col AreaShape_Center_X \
--moving-y-col AreaShape_Center_Y
```

### Minimum data requirement

After non-finite centroid rows are removed, each CSV must contain at least:

```text
3 valid centroids
```

The default registration fit itself requires at least:

```text
--min-anchors 8
```

matched anchors.

### Minimal example CSV

```csv
ImageNumber,ObjectNumber,AreaShape_Center_X,AreaShape_Center_Y
1,1,1900.47,2009.28
1,2,2198.16,1743.62
1,3,1437.21,1220.55
```

Run with explicit column overrides for this example.

---

## 4.7 Output directory

Specified by:

```text
--output
```

The directory is created automatically with:

```text
parents=True
exist_ok=True
```

Example:

```text
/home/oltij/Desktop/PVAligned/
```

---

## 4.8 Registration parameters

### Physical pixel size

Defaults:

```text
--pixel-size-x 0.312843137
--pixel-size-y 0.312836345
```

Units are µm/pixel.

### Phase-correlation initializer

Enabled by default.

For unequal native canvases, temporary bottom/right zero-padding is used only
inside this initializer. It does not alter either source image or its centroid
coordinates.

Disable with:

```text
--skip-phase-init
```

Subpixel upsampling:

```text
--phase-upsample 10
```

### Centroid matching gate

Default:

```text
--match-gate-um 25.0
```

The script converts this to pixels using the geometric mean of the X and Y pixel sizes.

### Iterations

```text
--iterations 5
```

### Minimum anchors

```text
--min-anchors 8
```

### Robust outlier cutoff

```text
--outlier-mad 3.5
```

### Maximum centroids used for dense matching

```text
--max-points 1500
```

Set:

```text
--max-points 0
```

to use all points.

### Random seed

```text
--seed 42
```

### Optional affine comparison

Enable:

```text
--try-affine
```

Default minimum median-residual improvement required to select affine:

```text
--affine-min-improvement 0.15
```

meaning 15%.

Force affine:

```text
--force-affine
```

Use the force option cautiously.

---

## 4.9 Output 1 — initial matched-centroid CSV

### Filename

```text
matched_centroids_initial.csv
```

### Exact columns

```text
moving_y_px
moving_x_px
fixed_y_px
fixed_x_px
initial_distance_px
initial_distance_um_approx
```

Example:

```csv
moving_y_px,moving_x_px,fixed_y_px,fixed_x_px,initial_distance_px,initial_distance_um_approx
1304.2,2115.8,1306.1,2117.0,2.247,0.703
```

This describes the initial globally assigned anchors after the phase-correlation translation.

---

## 4.10 Output 2 — final rigid-anchor CSV

### Filename

```text
matched_centroids_rigid.csv
```

### Exact columns

```text
moving_y_px
moving_x_px
fixed_y_px
fixed_x_px
transformed_moving_y_px
transformed_moving_x_px
residual_px
residual_um_approx
```

This is the principal anchor-quality table for the rigid fit.

---

## 4.11 Output 3 — rigid transform

### Filename

```text
transform_rigid.txt
```

### Format

CASTalign's serialized transform text format, written by:

```python
rigid_t.save(...)
```

This is not a generic CSV or manually edited transformation matrix.

---

## 4.12 Output 4 — optional affine transform

### Filename

```text
transform_affine_after_rigid.txt
```

### Conditional behavior

This file is saved only when:

1. `--try-affine` is enabled; and
2. affine is actually selected because either:
   - its median residual improves enough, or
   - `--force-affine` is used.

The affine component is fitted **after the rigid transform**.

---

## 4.13 Output 5 — final composed transform

### Filename

```text
transform_final.txt
```

### Meaning

This is the complete transform used for image resampling:

```text
phase-correlation translation
    +
selected CASTalign transform
```

The selected transform is normally rigid, or rigid + affine when affine is chosen.

This is the most appropriate transform file to retain if you need to reproduce the final mapping.

---

## 4.14 Output 6 — aligned moving-channel TIFF

### Current hard-coded filename

```text
aligned_PV_max.tif
```

### Format

```text
TIFF
2-D
dtype = float32
shape = same Y,X dimensions as fixed image
BigTIFF enabled
zlib compression
```

### Important naming warning

The filename is **hard-coded as PV**, even if the moving channel is LHX6 or another marker.

For example, running:

```text
moving = LHX6
```

still writes:

```text
aligned_PV_max.tif
```

unless the code is edited or the file is renamed afterward.

For a reusable repository version, this filename should ideally be made channel-configurable.

---

## 4.15 Output 7 — aligned fixed/reference TIFF

### Current hard-coded filename

```text
aligned_Hoechst_max.tif
```

### Format

```text
TIFF
2-D
same original dtype as the fixed TIFF
BigTIFF enabled
zlib compression
```

This is simply the fixed image written into the output directory so that both channels are available together in the same coordinate frame.

### Naming warning

Like the moving output, the name is hard-coded to `Hoechst`.

---

## 4.16 Output 8 — two-channel numeric TIFF

### Filename

```text
aligned_two_channel.tif
```

### Format

```text
TIFF
shape = (2, Y, X)
dtype = float32
BigTIFF
zlib compression
```

### Channel order

```text
channel 0 = fixed image
channel 1 = registered moving image
```

This is numeric image data, not an RGB visualization.

---

## 4.17 Output 9 — registration overlay PNG

### Filename

```text
registration_overlay.png
```

### Meaning

Display-normalized RGB overlay:

```text
green = fixed
red   = registered moving
yellow = intensity shared in both display channels
```

### Format

PNG, 300 dpi.

---

## 4.18 Output 10 — registration QC PNG

### Filename

```text
registration_qc.png
```

Three panels:

```text
fixed
registered moving
overlay
```

PNG, 250 dpi.

---

## 4.19 Output 11 — anchor QC PNG

### Filename

```text
anchor_qc.png
```

Shows:

```text
moving image + chosen moving anchors
fixed image + corresponding fixed anchors
registered moving over fixed
```

PNG, 250 dpi.

---

## 4.20 Output 12 — anchor residual histogram PNG

### Filename

```text
anchor_residual_histogram.png
```

Histogram of final rigid-anchor residuals in approximate µm.

PNG, 250 dpi.

---

## 4.21 Output 13 — registration summary text

### Filename

```text
registration_summary.txt
```

### Format

Plain UTF-8 text.

### Contents include

```text
fixed and moving input paths
CellProfiler CSV paths
image shape
pixel sizes
fixed/moving object counts
initial anchor count
final rigid anchor count
phase-correlation shift in pixels
phase-correlation shift in µm
correlation error
phase difference
rigid residual summary
optional affine comparison
selected model
scientific notes
```

---

## 4.22 Example registration output tree

```text
PVAligned/
├── aligned_PV_max.tif
├── aligned_Hoechst_max.tif
├── aligned_two_channel.tif
├── matched_centroids_initial.csv
├── matched_centroids_rigid.csv
├── transform_rigid.txt
├── transform_final.txt
├── transform_affine_after_rigid.txt       # only if affine selected
├── registration_overlay.png
├── registration_qc.png
├── anchor_qc.png
├── anchor_residual_histogram.png
└── registration_summary.txt
```

---

# Important implementation notes

## 26.1 Registration centroid column-name mismatch

The registration script auto-detects:

```text
Location_Center_X / Location_Center_Y
```

while the CellProfiler CSV used by the ROI verification script uses:

```text
AreaShape_Center_X / AreaShape_Center_Y
```

Use the registration command-line column overrides when necessary.

---

## 26.2 Registration output channel names are hard-coded

The files:

```text
aligned_PV_max.tif
aligned_Hoechst_max.tif
```

do not automatically adapt to different marker names.

This can become dangerous when the same registration script is used for LHX6, PV, GFP, etc.

---
