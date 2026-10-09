# General Image-Analysis Pipeline

---

## Overview

This README describes the end-to-end workflow for a tiled microscopy acquisition consisting of four unstitched fields/quadrants and one to four fluorescence channels.

The document is organized around **eight pipeline steps**. Those eight steps are the authoritative processing order. The four larger phases describe what the pipeline is doing conceptually; they do **not** add extra processing steps.

### Run the complete workflow from IMS files

`run_ims_pipeline.py` is the complete entry point when the source data are `.ims` files. It supports both `_F00.ims`…`_F03.ims` mosaic fields and one IMS file already stitched in Imaris. It assigns channel indices 0–3 to Hoechst, mNeonGreen, BiVe3/BiVe4 virus, and PV; resolves the virus label separately for each organoid; records Imaris names/wavelength metadata; retains full-resolution projections; creates explicitly configured area-mean analysis projections with automatically adjusted calibration; and then invokes every downstream stage automatically.

IMS preprocessing writes TIFF intermediates through validated temporary files
and atomic renames. Resume checks TIFF structure, dimensions, data type, and
pixel-data bounds instead of trusting a nonempty filename, so an interrupted
write is regenerated safely.

```bash
cp ims_pipeline.example.json ims_pipeline.json
# Edit IMS/output/SIF paths and scientific settings.
python run_ims_pipeline.py --config ims_pipeline.json --dry-run
python run_ims_pipeline.py --config ims_pipeline.json
```

Read [`IMS_WORKFLOW.md`](IMS_WORKFLOW.md) for input naming, the fixed channel mapping, mosaic versus already-stitched behavior, runtime requirements, outputs, and resume instructions.

### Start from Imaris-stitched TIFFs

When stitching is already complete in Imaris, `run_stitched_pipeline.py` provides a linear entry point at Step 2 and leaves the optional IMS/BaSiC/BigStitcher route intact. It validates one single-channel TIFF per channel, runs the existing scripts in order, carries every output path into the next stage, invokes the separate CellProfiler environment automatically, logs execution, and supports provenance-checked resume.

```bash
cp stitched_pipeline.example.json stitched_pipeline.json
# Edit input paths, calibration, thresholds, and filters.
python run_stitched_pipeline.py --config stitched_pipeline.json --dry-run
python run_stitched_pipeline.py --config stitched_pipeline.json
```

Read [`STITCHED_WORKFLOW.md`](STITCHED_WORKFLOW.md) before using experimental data. [`SCRIPT_HANDOFF_REFERENCE.md`](SCRIPT_HANDOFF_REFERENCE.md) lists the input/output contract for every original script. Historical pixel sizes and zero-percentile filter examples are not automatically valid for a new acquisition.

```text
ENV
 |
 v
FEED INTO → IMS LOADING
 |
 ├── STEP 1. .ims processing → stitched/max-projected images
 ├── STEP 2. Initial segmentation of ROI per marker
 └── STEP 3. Align images using ROI/centroid information
 |
 v
ROI SCRUBBING
 |
 ├── STEP 4. Resegment the aligned images
 ├── STEP 5. Apply intensity visualization + intensity filtering
 └── STEP 6. Apply shape-based filtering
 |
 v
ROI → CELL
 |
 └── STEP 7. Nuclear colocalization to convert marker ROIs
              into nuclear-anchored cell identities
 |
 v
PUTATIVE CELL TYPES
 |
 └── STEP 8. Marker colocalization to define pairwise
              marker-positive / putative cell-type relationships
```

### The eight steps, in order

1. **`.ims processing → stitched image`** — export the IMS fields, trim storage padding, apply BaSiC, stitch/fuse each channel with BigStitcher, and generate max projections.
2. **Initial segmentation of ROI per marker** — run CellProfiler on each unaligned max projection to obtain the masks/centroids needed for registration.
3. **Align images with ROI information** — align every non-nuclear marker image independently to the same nuclear reference.
4. **Resegment aligned images** — derive the final-coordinate tissue boundary from Hoechst, reuse it for every marker's otherwise independent CellProfiler detection, verify mask equality/containment, generate reference-based QC, and extract ROI pixel coordinates.
5. **Apply intensity visualization + filter** — inspect intensity distributions and perform the first ROI-scrubbing pass using intensity only.
6. **Apply shape-based filter** — inspect geometry only after the intensity gate, then apply the final shape/area gate to the intensity-passing ROI population.
7. **Apply nuclear colocalization for ROI → Cell identity** — match each marker ROI to a nuclear ROI. The persistent cross-run identity anchor is the **nuclear ObjectID**; pairwise `CellID` values remain run-local.
8. **Marker colocalization for putative cell types** — compare the nuclear-associated marker populations to one another, then generate pairwise QC and, when desired, all-nuclear-cell A+/B+, A+/B−, A−/B+, A−/B− summaries.

### Handoff rule

Each step should consume the output of the immediately preceding biological-processing step. In particular:

```text
unaligned segmentation
    is for registration

aligned segmentation
    is for downstream filtering/colocalization

intensity-filtered ROI table
    is the input to shape filtering

final intensity+shape-filtered ROI table
    is the input to nuclear colocalization

nuclear-associated marker populations
    are the input to marker-marker colocalization
```

This ordering is important because changing coordinate frame, resegmentation, or filtering changes which ROI representation is valid downstream.

---

## Terminology used throughout

---

### Field / quadrant

Each acquisition consists of four tiled microscope fields:

```text
F00
F01
F02
F03
```

For the current BigStitcher workflow, these correspond to the four positions in a 2 × 2 acquisition.

### Channel

A fluorescence channel is one biological/imaging signal, for example:

```text
Hoescht / DAPI
LHX6
PV
GFP
```

The four `.ims` field files can each contain multiple channels.

Therefore the pipeline does **not** require one separate four-file folder per channel when the channels are stored together in the IMS acquisition.

### ROI

An ROI is a segmentation object produced independently within one fluorescence channel.

Before nuclear matching, IDs such as:

```text
ROI 17
ROI 203
ROI 901
```

are channel-specific segmentation-object identifiers.

They should not yet be interpreted as global biological cell identities.

### CellID

The current colocalization script assigns `CellID` values to matched object pairs.

However, an important limitation is:

```text
CellID is local to one colocalization run.
```

A `CellID=5` generated in an LHX6-versus-Hoescht run is not automatically the same biological identifier as `CellID=5` generated in a PV-versus-Hoescht run.

The nuclear ROI/object ID is therefore the more appropriate persistent cross-run anchor if a global cell identity is needed.

This is discussed in detail in the colocalization section.

---

---

# ENV — Software and execution environment

---

The environment is setup infrastructure, not one of the eight biological/image-processing steps. Establish it first, then enter Step 1.

Two Conda environments are used:

```text
ims-mosaic
    main image-analysis / registration / filtering / colocalization environment

cellprofiler-native
    CellProfiler-only environment
```

The Fiji + BigStitcher SIF is also prepared here.

---

### Conda environments and switching

The pipeline uses **two separate Conda environments**:

```text
ims-mosaic
    = main microscopy/image-analysis environment

cellprofiler-native
    = CellProfiler environment
```

The Fiji + BigStitcher container is also built separately as a Singularity `.sif`.

#### 2.1 Create the `ims-mosaic` environment

Environment YAML:

```text
ims-mosaic.yml
```

Suggested location:

```text
/home/oltij/Desktop/ImageAnalysisScripts/CondaEnvironment/ims-mosaic.yml
```

Create:

```bash
conda env create \
    -f /home/oltij/Desktop/ImageAnalysisScripts/CondaEnvironment/ims-mosaic.yml
```

Activate:

```bash
conda activate ims-mosaic
```

Use `ims-mosaic` for all non-CellProfiler Python stages, including IMS export/stitching orchestration, BaSiC, TIFF processing, CASTalign, ROI reconstruction/export, filtering, metric reports, and colocalization.

#### 2.2 Create the `cellprofiler-native` environment

Environment YAML:

```text
cellprofiler-native.yml
```

Suggested location:

```text
/home/oltij/Desktop/ImageAnalysisScripts/CondaEnvironment/cellprofiler-native.yml
```

Create:

```bash
conda env create \
    -f /home/oltij/Desktop/ImageAnalysisScripts/CondaEnvironment/cellprofiler-native.yml
```

Activate:

```bash
conda activate cellprofiler-native
```

The provided YAML specifies:

```text
Python 3.9.25
CellProfiler 4.2.8
cellprofiler-core 4.2.8
```

Verify:

```bash
which cellprofiler
cellprofiler --version
```

#### 2.3 Switching between environments

Start the main workflow in:

```bash
conda activate ims-mosaic
```

Before any CellProfiler stage:

```bash
conda deactivate
conda activate cellprofiler-native
```

After CellProfiler is finished:

```bash
conda deactivate
conda activate ims-mosaic
```

Conceptually:

```text
ims-mosaic
    |
    | enter CellProfiler stage
    v
cellprofiler-native
    |
    | run CellProfiler
    v
ims-mosaic
```

This switch occurs for both:

```text
CellProfiler pass 1 = before alignment
CellProfiler pass 2 = after alignment
```

---

---

### Build the Fiji + BigStitcher SIF

Use the exact known-working definition:

```text
fiji_latest_bigstitcher.def
```

For example:

```text
/home/oltij/Desktop/ImageAnalysisScripts/Bigstitcher/fiji_latest_bigstitcher.def
```

Load Singularity:

```bash
module load singularity/4.4.1
```

Build:

```bash
cd /home/oltij/Desktop/ImageAnalysisScripts/Bigstitcher

singularity build --fakeroot \
    fiji_latest_bigstitcher.sif \
    fiji_latest_bigstitcher.def
```

The expected built file is:

```text
/home/oltij/Desktop/ImageAnalysisScripts/Bigstitcher/fiji_latest_bigstitcher.sif
```

The IMS pipeline does not require the SIF to have one universal filesystem location because its path is passed explicitly with:

```text
--fiji-sif
```

Therefore the important requirement is simply to know the exact SIF path.

Before starting a long run, confirm that the container launches Fiji successfully.

For example:

```bash
singularity exec \
    /home/oltij/Desktop/ImageAnalysisScripts/Bigstitcher/fiji_latest_bigstitcher.sif \
    /opt/Fiji.app/fiji-linux-x64 \
    --headless \
    --help
```

The current stitching Python script invokes Fiji inside the container using:

```text
/opt/Fiji.app/fiji
```

so it is also useful to verify that this path exists in the successfully built SIF:

```bash
singularity exec \
    /home/oltij/Desktop/ImageAnalysisScripts/Bigstitcher/fiji_latest_bigstitcher.sif \
    test -x /opt/Fiji.app/fiji
```

If your known-working SIF already passes this test, do not change the `.def`.

---

---

# Feed into → IMS loading

---

This phase converts the raw acquisition into a common set of channel images and then establishes the information needed to place every marker channel into the same nuclear-reference coordinate system.

It contains **Steps 1–3 of 8**.

---

## Step 1 of 8 — `.ims` processing → stitched image

---

**Input**

```text
four raw F00–F03 .ims acquisition fields
```

**Output handed to Step 2**

```text
one stitched/fused 2-D max projection per fluorescence channel
```

The stitching pipeline is run **once per acquisition**, not once per channel. It discovers and processes all channels automatically.

---

### Organize the raw four-field acquisition

The current IMS mosaic script discovers files named:

```text
<sample>_F00.ims
<sample>_F01.ims
<sample>_F02.ims
<sample>_F03.ims
```

Example:

```text
Exp17_157-2_Set2_DAPI_GFP_LHX6_PV_F00.ims
Exp17_157-2_Set2_DAPI_GFP_LHX6_PV_F01.ims
Exp17_157-2_Set2_DAPI_GFP_LHX6_PV_F02.ims
Exp17_157-2_Set2_DAPI_GFP_LHX6_PV_F03.ims
```

Place the four IMS fields for one acquisition in the same source directory.

The script expects exactly two digits after `F`:

```text
F00
F01
F02
F03
```

not:

```text
F0
F1
F2
F3
```

---

---

### Run the IMS → TIFF → BaSiC → BigStitcher pipeline

The script is:

```text
ims_mosaic_pipeline_illumination_bigstitcher.py
```

Run the pipeline **once for the acquisition folder**, not once per fluorescence channel.

Example:

```bash
conda activate ims-mosaic

python /path/to/ims_mosaic_pipeline_illumination_bigstitcher.py \
    "/path/to/folder/containing/the/four/IMS/files" \
    --fiji-sif \
    "/home/oltij/Desktop/ImageAnalysisScripts/Bigstitcher/fiji_latest_bigstitcher.sif"
```

The important defaults/current behaviors are:

```text
grid = 2 × 2
overlap = 10%
minimum phase-correlation r = 0.55
BaSiC sampled Z planes per tile = 12

trailing all-zero right/bottom IMS padding = removed
voxel size = read from IMS metadata
BigStitcher intensity-adjustment downsampling = 32
BigStitcher intensity-adjustment max inliers = 10000
offset-only regularization = 0.5
unmodified-intensity regularization = 0.5

previous pipeline-generated sample directory = cleaned before rerun
```

The zero-padding trim can be disabled with:

```bash
--no-trim-zero-padding
```

but this is not recommended for the padded IMS datasets for which the pipeline was designed.

The previous generated result can be preserved with:

```bash
--keep-previous
```

A more explicit command is:

```bash
python /path/to/ims_mosaic_pipeline_illumination_bigstitcher.py \
    "/path/to/folder/containing/the/four/IMS/files" \
    --fiji-sif \
    "/home/oltij/Desktop/ImageAnalysisScripts/Bigstitcher/fiji_latest_bigstitcher.sif" \
    --grid 2 2 \
    --overlap 10 \
    --min-correlation 0.55 \
    --basic-sample-planes 12 \
    --intensity-downsampling 32 \
    --intensity-max-inliers 10000 \
    --intensity-offset-only 0.5 \
    --intensity-unmodified 0.5
```

The script also accepts a manual calibration override:

```bash
--voxel-size X_UM Y_UM Z_UM
```

but this is normally unnecessary because voxel size is read directly from the IMS metadata.

---

---

### Important correction: the IMS script processes all channels automatically

The current IMS mosaic script reads the channel list from the IMS files and then loops over:

```python
for channel in channels:
```

Therefore, if the four IMS files contain:

```text
Hoescht
GFP
LHX6
PV
```

you do **not** run the script four times.

You run it once on the acquisition folder.

The pipeline automatically performs, independently for each available channel:

```text
IMS export
→
trailing zero-padding trim
→
BaSiC illumination correction
→
BigStitcher spatial registration
→
BigStitcher overlap-based intensity adjustment
→
32-bit fusion
→
maximum-intensity projection
```

If only two or three channels exist, it processes those channels.

The important condition is that every field in the acquisition must have the same channel names/order.

---

---

### What the IMS pipeline produces

For each sample, the script creates:

```text
<source-folder-name>_IMS_to_TIFF/
└── <sample-name>/
    ├── 01_exported/
    ├── 02_basic_corrected/
    ├── 03_stitched/
    └── 04_max_projections/
```

#### `01_exported`

Contains the individual channel Z-stacks exported from each IMS field.

Example:

```text
F00_Confocal - Blue.tif
F01_Confocal - Blue.tif
F02_Confocal - Blue.tif
F03_Confocal - Blue.tif
```

The pipeline trims trailing all-zero XY storage padding by default.

This directory also contains one geometry audit file per field:

```text
F00_geometry.txt
F01_geometry.txt
F02_geometry.txt
F03_geometry.txt
```

Each records the source IMS array size, exported logical image size, and the
number of trailing rows/columns removed. These files should be retained with
the run.

#### `02_basic_corrected`

Contains BaSiC-corrected Z-stacks.

#### `03_stitched`

Contains the BigStitcher job directories and final fused 3-D stack for every channel.

Conceptually:

```text
<channel>_stitched.tif
```

Each channel gets its own independently stitched 3-D image stack.

The current BigStitcher sequence is:

```text
phase correlation
→
pairwise-link filtering
→
global optimization
→
overlap-based intensity adjustment
→
32-bit float blended fusion
```

Each channel job also retains audit/debugging files such as:

```text
dataset.xml
run_bigstitcher.ijm
bigstitcher.log
```

The wrapper explicitly checks that intensity adjustment was reached and
completed before fusion is accepted.

#### `04_max_projections`

Contains the final 2-D maximum-intensity projection for each channel:

```text
<sample>_<channel>_stitched_max_projection.tif
```

These are the images that should normally enter the next 2-D analysis stages.

The maximum projection is computed directly from the final fused stack
without normalization or dtype conversion, so the expected 32-bit
floating-point fused intensity representation is preserved.

---

---

### Important correction: a separate max-projection step is normally unnecessary

The standalone max-projection script performs:

```python
projection = np.max(image, axis=0)
```

on a fused Z-stack.

However, the current IMS mosaic script already performs exactly this operation at the end of every channel run and writes the result into:

```text
04_max_projections/
```

Therefore, for the current end-to-end IMS script:

```text
DO NOT need:
stitched stack
→ manually run standalone max-projection script
```

because the current pipeline already does:

```text
stitched stack
→ max projection automatically
```

Use the standalone max-projection script only when:

```text
you have a fused stack produced outside the IMS pipeline
```

or:

```text
you intentionally want to regenerate a projection manually
```

This avoids producing duplicate projection files.

---

---

## Step 2 of 8 — Initial segmentation of ROI per marker

---

**Input**

```text
unaligned channel max projections from Step 1
```

**Purpose**

Create the initial channel-specific segmentation masks and object centroids required for image registration.

**Output handed to Step 3**

```text
unaligned fluorescence TIFF
+ CellProfiler object/centroid CSV
```

These pre-alignment ROIs are **registration support objects**. They are not the ROI pixel tables that should be used for final filtering or colocalization after the image coordinate frame changes.

---

### First CellProfiler pass: segment the unaligned channel projections

At this stage, use the max projections from:

```text
04_max_projections/
```

First switch into the CellProfiler environment:

```bash
conda deactivate
conda activate cellprofiler-native
```

There are **two supported ways to run CellProfiler**.

#### Option A — GUI

Launch:

```bash
cellprofiler
```

and open:

```text
MGEOPV.cpproj
```

Use this option when you want to inspect segmentation interactively, adjust thresholds/settings, step through modules, or visually confirm ROI identification.

#### Option B — headless

Use:

```text
MGEOPVFinal.cppipe
```

with:

```text
cellprofilerdriver.py
```

Example:

```bash
python /home/oltij/Desktop/cellprofilerdriver.py \
    --pipeline "/home/oltij/Desktop/MGEOPVFinal.cppipe" \
    --input "/path/to/channel_max_projection.tif" \
    --output "/path/to/CellProfiler_output"
```

The two options are:

```text
GUI:
    MGEOPV.cpproj

Headless:
    MGEOPVFinal.cppipe
    +
    cellprofilerdriver.py
```

Run CellProfiler independently for each channel image:

```text
nuclear projection
marker 1 projection
marker 2 projection
marker 3 projection
```

The required CellProfiler outputs for the ROI-extraction/alignment workflow include the relevant object-measurement CSVs and segmentation image.

Typical files used by the scripts are:

```text
CellMask.tiff
MyExpt_FilterObjects2.csv
MyExpt_FilterObjects.csv
```

`MyExpt_FilterObjects2.csv` contains individual object measurements such as:

```text
ObjectNumber
AreaShape_Area
AreaShape_BoundingBoxMinimum_X
AreaShape_BoundingBoxMaximum_X
AreaShape_BoundingBoxMinimum_Y
AreaShape_BoundingBoxMaximum_Y
AreaShape_Center_X
AreaShape_Center_Y
```

`MyExpt_FilterObjects.csv` can contain the organoid-level object/area information used by some downstream scripts.

After all pre-alignment CellProfiler runs are complete, switch back to:

```bash
conda deactivate
conda activate ims-mosaic
```

Continue alignment and all other non-CellProfiler stages from `ims-mosaic`.

---

---

### Are pre-alignment ROI pixel CSVs required?

There is an important distinction.

The CASTalign registration script uses:

```text
the fluorescence TIFF
+
CellProfiler centroid CSV
```

It does **not** require the extracted:

```text
<Channel>_ROI_pixels.csv
```

for registration.

Therefore:

```text
CellProfiler before alignment = required
ROI-pixel extraction before alignment = optional for registration itself
```

You may still run ROI extraction before alignment for:

```text
QC
archival purposes
verification of CellProfiler masks
pre-alignment analyses
```

but the alignment code does not consume those pixel CSVs.

The ROI extraction step becomes essential **after alignment**, because those aligned ROI coordinates are used for filtering and colocalization.

---

---

### Optional pre-alignment ROI extraction

If desired, run the ROI reconstruction/extraction script for each channel.

For each channel it uses:

```text
CellMask.tiff
channel CellProfiler CSV
organoid CellProfiler CSV
```

and produces, among other outputs:

```text
ROI_analysis/<Channel>_ROI_pixels.csv
```

Format:

```csv
ROI,X,Y
1,100,200
1,101,200
1,102,200
2,450,775
```

Each row is one pixel belonging to one segmentation object.

---

---

## Step 3 of 8 — Align images with ROI information

---

**Input**

```text
fixed image = nuclear max projection
moving image = one non-nuclear marker max projection

fixed centroid CSV = nuclear CellProfiler objects
moving centroid CSV = marker CellProfiler objects
```

**Output handed to Step 4**

```text
aligned nuclear-reference image
aligned marker images, each independently transformed into that same reference frame
```

Every non-nuclear channel is aligned **directly to the nuclear reference**. Do not build a transform chain such as marker A → marker B → nucleus.

---

### Choose the nuclear reference channel

One channel should be designated as the common spatial reference.

Examples:

```text
Hoescht
DAPI
```

From this point forward:

```text
fixed image = nuclear channel
moving image = one non-nuclear marker channel
```

Every non-nuclear channel is aligned independently to the same nuclear image.

For example:

```text
LHX6 → Hoescht
PV   → Hoescht
GFP  → Hoescht
```

Do not align LHX6 to PV and then PV to Hoescht.

All markers should independently reference the same nuclear coordinate system.

---

---

### Align every non-nuclear channel to the nuclear reference

Use:

```text
CASTalign_2D_maxproj_registration.py
```

For each non-nuclear marker, provide:

```text
fixed TIFF = nuclear max projection
moving TIFF = marker max projection

fixed CSV = nuclear CellProfiler object CSV
moving CSV = marker CellProfiler object CSV
```

The CellProfiler files in this workflow use centroid columns:

```text
AreaShape_Center_X
AreaShape_Center_Y
```

whereas the alignment script's automatic column detection is oriented toward `Location_Center_*`.

Therefore explicitly provide the CellProfiler centroid column names.

Generic example:

```bash
python CASTalign_2D_maxproj_registration.py \
    --fixed-image "/path/to/Hoescht_max_projection.tif" \
    --moving-image "/path/to/LHX6_max_projection.tif" \
    --fixed-csv "/path/to/Hoescht/MyExpt_FilterObjects2.csv" \
    --moving-csv "/path/to/LHX6/MyExpt_FilterObjects2.csv" \
    --fixed-x-col AreaShape_Center_X \
    --fixed-y-col AreaShape_Center_Y \
    --moving-x-col AreaShape_Center_X \
    --moving-y-col AreaShape_Center_Y \
    --output "/path/to/LHX6Aligned"
```

Repeat for every non-nuclear channel.

---

---

### Alignment output naming

The uploaded alignment script currently uses some hard-coded output names associated with its original Hoescht/PV use case.

Therefore, after each run, verify exactly which file is the registered moving image.

Keep the final registered marker files under clear names such as:

```text
aligned_LHX6_max.tif
aligned_PV_max.tif
aligned_GFP_max.tif
```

and preserve the nuclear reference as:

```text
aligned_Hoechst_max.tif
```

or equivalent.

The scientific requirement is more important than the literal filename:

```text
all final channel images must be in the same nuclear-reference coordinate frame
```

---

---

### The nuclear image itself does not need a transformation

The nuclear image is the fixed reference.

Therefore:

```text
nuclear aligned image = nuclear reference image
```

It does not need to be warped against itself.

However, keeping a copy in an aligned-analysis directory can make downstream path organization easier.

---

---

# ROI scrubbing

---

This phase rebuilds the segmentation in the final aligned coordinate system and then removes poor-quality ROIs in a strict two-stage order.

It contains **Steps 4–6 of 8**.

```text
aligned images
    ↓
resegment
    ↓
intensity gate
    ↓
shape/area gate
    ↓
final filtered ROI population
```

---

## Step 4 of 8 — Resegment aligned images

---

**Input**

```text
aligned images from Step 3
```

**Output handed to Step 5**

```text
unaltered aligned ROI pixel CSV for each channel
```

Resegmentation after alignment is mandatory for downstream pixel-overlap analysis because the moving images have entered a new coordinate frame.

---

### Second CellProfiler pass: segment the aligned images

After alignment, run CellProfiler again on the final aligned images.

This second CellProfiler pass is essential.

Do not use the pre-alignment ROI pixel coordinates for downstream overlap analysis because the moving fluorescence image has changed coordinate frame.

Switch back into the CellProfiler environment:

```bash
conda deactivate
conda activate cellprofiler-native
```

Again choose one of two execution modes.

#### Option A — GUI

Launch:

```bash
cellprofiler
```

and open:

```text
MGEOPV.cpproj
```

Use the aligned image for the corresponding channel.

#### Option B — headless

Use:

```text
MGEOPVFinal.cppipe
```

with:

```text
cellprofilerdriver.py
```

Example:

```bash
python /home/oltij/Desktop/cellprofilerdriver.py \
    --pipeline "/home/oltij/Desktop/MGEOPVFinal.cppipe" \
    --input "/path/to/aligned_channel.tif" \
    --output "/path/to/aligned_channel_CellProfiler_output"
```

Run CellProfiler for:

```text
aligned nuclear reference
aligned marker 1
aligned marker 2
aligned marker 3
```

The nuclear image may be unchanged geometrically, but rerunning it in the aligned-analysis workflow keeps the segmentation products and directory structure consistent.

After all aligned CellProfiler runs are complete:

```bash
conda deactivate
conda activate ims-mosaic
```

Continue ROI extraction, filtering, reports, and colocalization from `ims-mosaic`.

---

---

### Extract ROIs from all aligned channels

Run the ROI reconstruction/extraction script on each aligned channel's CellProfiler outputs.

For each channel, the important downstream result is:

```text
<Channel>/ROI_analysis/<Channel>_ROI_pixels.csv
```

Example:

```text
HoeschtAligned/ROI_analysis/Hoescht_ROI_pixels.csv
LHX6Aligned/ROI_analysis/LHX6_ROI_pixels.csv
PVAligned/ROI_analysis/PV_ROI_pixels.csv
```

These are the **unaltered aligned segmentation-object pixel tables**.

They should be preserved as the starting point for filtering.

---

---

## Step 5 of 8 — Apply intensity visualization + filter

---

**Input**

```text
unaltered aligned ROI pixel CSVs from Step 4
+ aligned fluorescence TIFFs
```

**Order inside this step**

```text
visualize/report intensity
    ↓
choose an intensity metric/cutoff
    ↓
run intensity-only filtering
```

**Output handed to Step 6**

```text
intensity-filtered ROI pixel CSV
```

Shape filters remain disabled during this step.

---

### Filtering philosophy

The intended filtering strategy is deliberately two-stage:

```text
Stage 1:
remove low-quality/dim objects using intensity

Stage 2:
inspect shape among the surviving intensity-positive objects,
then apply shape/area filtering
```

The reason for this order is to avoid letting poorly segmented low-intensity objects dominate the apparent distribution of:

```text
area
circularity
eccentricity
solidity
```

In other words:

```text
intensity quality gate first
shape quality gate second
```

---

---

### Intensity visualization

Use the ROI intensity-report script on the **aligned, unfiltered ROI pixel CSVs**.

Inputs per channel:

```text
aligned <Channel>_ROI_pixels.csv
aligned fluorescence TIFF
```

The script reports:

```text
Mean intensity
Median intensity
Maximum intensity
Integrated intensity
Intensity SD
```

and produces percentile-based visualization/QC outputs.

Inspect these manually to decide which intensity metric(s) and percentile cutoff(s) are appropriate.

Do not choose shape thresholds yet.

---

---

### First filtering pass: intensity only

Use the final metric filtering script.

For this first pass:

```text
INPUT:
unaltered aligned ROI pixel CSVs

ENABLE:
chosen intensity filter(s)

DISABLE:
all shape filters
```

In this filtering script:

```text
percentile = 0
```

means:

```text
filter OFF
```

Therefore the shape entries should remain at:

```python
{"percentile": 0, ...}
```

during the intensity-only pass.

Example conceptual configuration:

```python
FILTERS = {
    "Mean intensity":        {"percentile": 10, "keep": "above"},

    "Median intensity":      {"percentile": 0, "keep": "above"},
    "Maximum intensity":     {"percentile": 0, "keep": "above"},
    "Integrated intensity":  {"percentile": 0, "keep": "above"},
    "Intensity SD":          {"percentile": 0, "keep": "above"},

    "Area_um2":              {"percentile": 0, "keep": "above"},
    "Perimeter_um":          {"percentile": 0, "keep": "above"},
    "Solidity":              {"percentile": 0, "keep": "above"},
    "Eccentricity":          {"percentile": 0, "keep": "below"},
    "Circularity":           {"percentile": 0, "keep": "above"},
    "ConvexHullArea_um2":    {"percentile": 0, "keep": "above"},
    "MajorAxisLength_um":    {"percentile": 0, "keep": "above"},
    "MinorAxisLength_um":    {"percentile": 0, "keep": "above"},
    "AspectRatio":           {"percentile": 0, "keep": "below"},
    "EquivalentDiameter_um": {"percentile": 0, "keep": "above"},
    "Extent":                {"percentile": 0, "keep": "above"},
}
```

The important output is:

```text
<Channel>_ROI_pixels_FILTERED.csv
```

At this point, interpret it as:

```text
intensity-filtered ROI population
```

rather than the final population.

Keep these files.

---

---

## Step 6 of 8 — Apply shape-based filter

---

**Input**

```text
intensity-filtered ROI pixel CSV from Step 5
```

**Order inside this step**

```text
inspect area/shape distributions of intensity-passing objects
    ↓
choose shape/area thresholds
    ↓
run the second filtering pass on the intensity-filtered CSV
```

**Output handed to Step 7**

```text
final intensity+shape-filtered aligned ROI pixel CSVs
```

The intensity gate must remain logically upstream: shape filtering must **not** restart from the original unfiltered ROI table.

---

### Inspect shape/area only after the intensity gate

Run the separate report scripts on the **intensity-filtered ROI pixel CSVs**:

```text
Area report
Circularity report
Eccentricity report
Solidity report
```

These scripts calculate the full geometry for surviving objects and then visualize/rank the selected metric.

This means the shape distributions now describe:

```text
objects that already passed the intensity quality gate
```

rather than being contaminated by the dim objects you intentionally removed.

Manually inspect the reports and choose any desired shape/area thresholds.

---

---

### Critical rule for the second filtering pass

Your proposed logic was:

```text
turn the intensity filter off
then enable the chosen shape filters
```

That is correct **only if the input to the second filtering pass is the intensity-filtered CSV produced by pass 1**.

The second pass must therefore be:

```text
INPUT:
pass-1 <Channel>_ROI_pixels_FILTERED.csv

INTENSITY FILTERS:
off / percentile = 0

SHAPE FILTERS:
chosen values enabled
```

This produces:

```text
intensity-passed
AND
shape-passed
```

objects.

---

---

### What not to do during the second filtering pass

Do **not** do this:

```text
original unfiltered aligned ROI CSV
+
intensity filter OFF
+
shape filter ON
```

because the low-intensity objects removed in pass 1 would be eligible again.

That would break the intended intensity-first gate.

---

---

### Recommended two-pass filtering directory structure

A clean structure is:

```text
Filtering/
├── 01_intensity_reports/
├── 02_intensity_filter/
│   ├── Hoescht/
│   │   └── Hoescht_ROI_pixels_FILTERED.csv
│   ├── LHX6/
│   │   └── LHX6_ROI_pixels_FILTERED.csv
│   └── PV/
│       └── PV_ROI_pixels_FILTERED.csv
│
├── 03_shape_reports_after_intensity/
│   ├── area/
│   ├── circularity/
│   ├── eccentricity/
│   └── solidity/
│
└── 04_final_shape_filter/
    ├── Hoescht/
    │   └── Hoescht_ROI_pixels_FILTERED.csv
    ├── LHX6/
    │   └── LHX6_ROI_pixels_FILTERED.csv
    └── PV/
        └── PV_ROI_pixels_FILTERED.csv
```

The files under:

```text
04_final_shape_filter/
```

are the final filtered segmentation objects used for colocalization.

---

---

### Why the second filtering pass works sequentially

Suppose the original aligned channel has:

```text
3,000 ROIs
```

Intensity filtering retains:

```text
2,300 ROIs
```

Then the shape-filtering script receives only those 2,300 ROI pixel sets.

The shape percentiles are therefore calculated from:

```text
the 2,300 intensity-passing ROIs
```

and perhaps retain:

```text
2,050 final ROIs
```

The final population is:

```text
original ROIs
∩
intensity quality
∩
shape quality
```

which matches the intended biological/QC rationale.

---

---

### Final filtered images must all share one coordinate frame

Before beginning colocalization, verify that the following all correspond to the same aligned image canvas:

```text
final nuclear ROI pixels
final marker-1 ROI pixels
final marker-2 ROI pixels
final marker-3 ROI pixels

aligned nuclear TIFF
aligned marker-1 TIFF
aligned marker-2 TIFF
aligned marker-3 TIFF
```

This is essential because colocalization uses exact shared:

```text
X,Y
```

pixel coordinates.

---

---

# ROI → Cell

---

This phase turns channel-specific marker ROIs into **nuclear-associated cell identities**.

It contains **Step 7 of 8**.

The important identity distinction is:

```text
marker ROI ID
    = channel-specific segmentation object

pairwise CellID
    = convenience ID local to one colocalization run

nuclear ObjectID
    = persistent canonical anchor used to connect the same cell across marker runs
```

---

## Step 7 of 8 — Nuclear colocalization: ROI → Cell identity

---

For each marker, run a separate marker-versus-nucleus colocalization.

**Input**

```text
Object A = final filtered marker ROI pixels
Object B = final filtered nuclear ROI pixels

Object A TIFF = aligned marker fluorescence
Object B TIFF = aligned nuclear fluorescence
```

**Core output handed to Step 8**

```text
marker↔nucleus analysis HDF5
+ nuclear-positive marker pixel table
+ mapping from marker ObjectID to canonical nuclear ObjectID
```

This is the point where a marker ROI is promoted from a channel-specific object to a nuclear-associated cell. Preserve the marker↔nucleus HDF5 files because the later all-cell summaries need those canonical nuclear mappings.

---

### First colocalization stage: marker versus nuclear stain

For each biological marker, perform a separate colocalization against the nuclear reference.

Examples:

```text
LHX6 versus Hoescht
PV versus Hoescht
GFP versus Hoescht
```

For this stage, a useful convention is:

```text
Object A = marker
Object B = nuclear stain
```

because the current eligibility threshold is directional.

It tests:

```text
A_Percent_Inside_B
```

against the configured threshold.

Thus with:

```text
A = marker
B = nucleus
```

the threshold asks how much of the marker ROI lies inside its candidate nuclear ROI.

Choose the threshold according to the biology/segmentation logic of the marker rather than assuming the same threshold is appropriate for every marker.

---

---

### Inputs to `01_run_colocalization.py`

For a marker-versus-nuclear run:

```text
OBJECT_A_NAME = marker
OBJECT_B_NAME = nuclear stain

OBJECT_A_CSV = final filtered marker ROI-pixel CSV
OBJECT_B_CSV = final filtered nuclear ROI-pixel CSV

OBJECT_A_TIFF = aligned marker TIFF
OBJECT_B_TIFF = aligned nuclear TIFF
```

Optional:

```text
ORGANOID_CSV
```

can be supplied when density/area calculations are desired.

The script requires the two TIFFs to have the same image shape.

---

---

### Colocalization matching behavior

The current analysis:

```text
calculates exact shared pixels
calculates directional overlap percentages
applies the overlap threshold
builds preference lists
performs one-to-one stable matching
assigns CellIDs to matched pairs
```

The threshold is based on:

```text
A_Percent_Inside_B
```

not the reverse percentage.

The stable matching is one-to-one.

Therefore one marker object cannot ultimately be assigned to several nuclear objects, and one nuclear object cannot ultimately be assigned to several marker objects within that run.

---

---

### What `01_run_colocalization.py` writes

The script writes one self-contained:

```text
.h5
```

analysis file.

It contains:

```text
original fluorescence images
original object pixel tables
overlap table
eligible-pair table
preference matrix
stable matching table
CellID mapping
filtered matched Object-A pixels
filtered matched Object-B pixels
positive-object table
unmatched-object information
summary statistics
filtering-stage counts
```

Important internal tables include:

```text
tables/cell_id_mapping
tables/filtered_object_a_pixels
tables/filtered_object_b_pixels
```

---

---

### Nuclear-positive marker objects

When:

```text
Object A = marker
Object B = nuclear
```

the HDF5 table:

```text
tables/filtered_object_a_pixels
```

contains only marker objects that survived the one-to-one marker-versus-nuclear matching.

These are the marker segmentation objects that can be treated as:

```text
nuclear-associated / nuclear-positive marker cells
```

under the chosen overlap criterion.

The table is formatted using:

```text
CellID
X
Y
```

and is therefore directly compatible with the next colocalization analysis after exporting it to CSV.

---

---

### Important CellID limitation

The current script creates `CellID` using:

```text
1, 2, 3, ..., N
```

after sorting the matched pairs in that specific analysis run.

Therefore:

```text
LHX6-vs-Hoescht CellID 10
```

and:

```text
PV-vs-Hoescht CellID 10
```

are **not guaranteed to refer to the same nucleus**.

The values are only shared between Object A and Object B within one `.h5` result.

---

---

### Recommended nomenclature for global cell identity

For a stable identity across marker-vs-nuclear runs, preserve:

```text
B_ObjectID
```

from:

```text
tables/cell_id_mapping
```

when Object B is the nuclear channel.

Conceptually:

```text
Nuclear ObjectID = canonical nucleus/cell key
```

while:

```text
CellID = pairwise-run-local convenience identifier
```

If a future script is added to construct a global multi-marker table, it should join marker relationships through this nuclear `B_ObjectID`.

---

---

### Nuclear-colocalization visualization/QC

---

After each analysis HDF5 is created, run:

```text
02_generate_colocalization_visualizations.py
```

Set:

```text
RESULTS_FILE
```

to the corresponding HDF5.

The visualization script reads the self-contained HDF5 written by
`01_run_colocalization.py`; it does **not** rerun stable matching.

The ordinary visualization outputs include:

```text
matched-pair examples
below-threshold examples
filtering-stage counts
filtering Sankey
before/after ROI overlays
full-resolution ROI-overlay TIFFs
matched-object label TIFFs
RGB matched-object TIFF
5%-bin overlap-threshold QC
unmatched-object R/G/Y QC
unmatched-object geometry/intensity metrics
```

The visualization stage also now creates population-level summary outputs in:

```text
visualizations/unmatched_only_qc/
```

including:

```text
<ObjectA>_<ObjectB>_final_class_counts.csv
<ObjectA>_<ObjectB>_final_summary_metrics.csv
<ObjectA>_<ObjectB>_final_unmatched_class_5pct_bin_summary.csv
```

`final_summary_metrics.csv` contains, among other fields:

```text
total Object-A ROIs
total Object-B ROIs
double-positive count
Object-A-only count
Object-B-only count

Object A assigned to Object B %
Object A not assigned to Object B %
Object B assigned to Object A %
Object B not assigned to Object A %

A+/B+ % of the final classified marker population
A+/B- % of the final classified marker population
A-/B+ % of the final classified marker population
```

The denominator for those three marker-population percentages is:

```text
double-positive
+ Object-A-only
+ Object-B-only
```

so this three-class summary does **not** include double-negative cells.

#### Optional organoid-area summary

The visualization script can also calculate marker-positive ROI area as a
percentage of total organoid area.

Configure:

```python
ENABLE_ORGANOID_AREA_SUMMARY = True
ORGANOID_MEASUREMENTS_CSV = Path(
    "/path/to/MyExpt_FilterObjects.csv"
)
ORGANOID_AREA_COLUMN = "AreaShape_Area"
```

The CellProfiler organoid CSV supplies the denominator.

The resulting summary includes fields such as:

```text
Organoid_Area_Pixels
Organoid_Area_um2

Object_A_Positive_Area_Pixels
Object_A_Positive_Area_um2
Object_A_Positive_Area_Percent_of_Organoid

Object_B_Positive_Area_Pixels
Object_B_Positive_Area_um2
Object_B_Positive_Area_Percent_of_Organoid
```

If multiple organoid rows are present, their `AreaShape_Area` values are
summed.

---

### Export the nuclear-positive marker pixels for the next stage

The current `01_run_colocalization.py` stores the matched marker pixel table **inside the HDF5**.

It does not automatically write that table as a standalone CSV for the next analysis stage.

Therefore an export step is required before marker-vs-marker colocalization.

A minimal helper is:

```python
import io
from pathlib import Path

import h5py
import pandas as pd


results_file = Path(
    "/path/to/marker_vs_nuclear_analysis.h5"
)

output_csv = Path(
    "/path/to/Marker_nuclear_positive_cells.csv"
)


with h5py.File(results_file, "r") as h5:
    payload = h5[
        "tables/filtered_object_a_pixels"
    ][()]

    marker_cells = pd.read_csv(
        io.BytesIO(payload.tobytes())
    )


marker_cells.to_csv(
    output_csv,
    index=False
)

print(f"Saved: {output_csv}")
print(marker_cells.head())
```

The exported CSV will have the form:

```csv
CellID,X,Y
1,1200,1800
1,1201,1800
1,1201,1801
2,2340,1670
```

Because `01_run_colocalization.py` accepts either:

```text
ROI
```

or:

```text
CellID
```

as its object identifier, this exported table can be used directly in the next colocalization analysis.

---

---

### Export one nuclear-positive CSV per marker

After nuclear matching, the workflow should produce something conceptually like:

```text
NuclearPositiveCells/
├── LHX6_nuclear_positive_cells.csv
├── PV_nuclear_positive_cells.csv
└── GFP_nuclear_positive_cells.csv
```

Each contains only the pixels belonging to marker objects that passed nuclear matching.

Also retain the corresponding original marker-versus-nuclear `.h5` files.
The exported CSVs are used as marker-versus-marker inputs, while the HDF5
files preserve the canonical marker-to-nuclear ObjectID mapping needed by the
later all-nuclear-cell summary.

---

---

# Putative cell types

---

This phase compares the **nuclear-associated marker populations** produced in Step 7.

It contains **Step 8 of 8**.

The output is not a new universal cell identity namespace. Instead, marker-marker overlap/matching defines pairwise marker relationships that can be interpreted as putative multi-marker cell classes under the chosen segmentation and overlap criteria.

---

## Step 8 of 8 — Marker colocalization

---

**Input**

```text
nuclear-positive marker-A pixel table
nuclear-positive marker-B pixel table
aligned marker-A TIFF
aligned marker-B TIFF
```

**Main outputs**

```text
marker-A ↔ marker-B analysis HDF5
matched / A-only / B-only relationships
pairwise QC visualizations
summary metrics
optional organoid-area coverage
optional all-nuclear-cell A+/B+, A+/B−, A−/B+, A−/B− classification
```

For the all-nuclear-cell summary, the visualization stage reuses the two original marker↔nucleus HDF5 files from Step 7. It does **not** infer global identity from equality of pairwise `CellID` values.

---

### Second colocalization stage: marker versus marker

Now compare the nuclear-positive marker populations to each other.

Examples:

```text
LHX6-positive cells versus PV-positive cells
LHX6-positive cells versus GFP-positive cells
PV-positive cells versus GFP-positive cells
```

Inputs are:

```text
Object A CSV:
nuclear-positive marker-A CellID/X/Y table

Object B CSV:
nuclear-positive marker-B CellID/X/Y table

Object A TIFF:
aligned marker-A fluorescence image

Object B TIFF:
aligned marker-B fluorescence image
```

Because all images were aligned to the same nuclear reference, their coordinates should now be directly comparable.

After the marker-versus-marker HDF5 is created, run
`02_generate_colocalization_visualizations.py` on that result and point its
two nuclear-source settings back to the original marker-versus-nuclear HDF5
files. This allows the pairwise marker result to be summarized both as:

```text
matched / A-only / B-only marker objects
```

and as the biologically broader:

```text
A+/B+
A+/B-
A-/B+
A-/B-
```

fractions of all nuclear cells.

---

---

### Pairwise marker relationships

For a four-channel experiment containing:

```text
nuclear
A
B
C
```

the nuclear-gating analyses are:

```text
A ↔ nuclear
B ↔ nuclear
C ↔ nuclear
```

Potential downstream marker-marker analyses are:

```text
A ↔ B
A ↔ C
B ↔ C
```

Only run the pairwise comparisons relevant to the biological question.

---

---

### Pairwise marker CellIDs are again local

The second marker-marker colocalization run creates another new run-local:

```text
CellID
```

namespace.

That is acceptable for pairwise overlap analysis because the code determines matches from:

```text
pixel coordinates
```

rather than assuming that input IDs already match across channels.

However, pairwise analysis alone does not automatically create one universal:

```text
cell × all markers
```

table.

For a future global multi-marker identity table, join all results through the corresponding nuclear object IDs from the first nuclear-gating stage.

---

---

### Recommended conceptual identity hierarchy

Use the following terminology consistently:

```text
ROI
    = one segmentation object in one channel

Nuclear-positive marker ROI
    = marker ROI successfully matched to one nuclear ROI

Pairwise CellID
    = ID assigned to a matched pair inside one colocalization run

Canonical nuclear ID
    = nuclear ObjectID used to connect identities across runs
```

This avoids accidentally treating two unrelated run-local CellID numbers as the same cell.

---

---

### All-nuclear-cell classification after marker colocalization

---

This is the downstream portion of the colocalization visualization/QC stage that is specifically relevant to Step 8.

For a marker-versus-marker result, the visualization script can recover the
original marker-to-nucleus assignments from the two **pre-existing**
marker-versus-nuclear HDF5 analyses.

Configure:

```python
ENABLE_NUCLEAR_ALL_CELL_SUMMARY = True

OBJECT_A_NUCLEAR_RESULTS_FILE = Path(
    "/path/to/ObjectA_vs_Nuclear_analysis.h5"
)

OBJECT_B_NUCLEAR_RESULTS_FILE = Path(
    "/path/to/ObjectB_vs_Nuclear_analysis.h5"
)
```

This is important because pairwise `CellID` values are run-local. The
visualization script therefore maps the current marker objects back to the
canonical nuclear `ObjectID` stored in the original marker↔nucleus matching
results.

It does **not** spatially rematch the markers to nuclei.

Once both marker populations are linked back to the canonical nuclear
namespace, every nuclear ROI is classified into exactly one of:

```text
A+/B+
A+/B-
A-/B+
A-/B-
```

The corresponding outputs are:

```text
<ObjectA>_<ObjectB>_all_nuclear_cells_classification.csv
<ObjectA>_<ObjectB>_all_nuclear_cells_class_counts.csv
```

and the four class counts/percentages are also incorporated into:

```text
<ObjectA>_<ObjectB>_final_summary_metrics.csv
```

The denominator for these percentages is **all canonical nuclear ROIs**,
including nuclei negative for both markers.

This is the appropriate summary when the biological question is:

```text
What percentage of all cells are A+/B+, A+/B-, A-/B+, or A-/B-?
```

rather than merely:

```text
How are the marker-positive objects partitioned?
```

#### Nuclear-denominator compatibility check

The two source marker↔nucleus HDF5 files should refer to the same nuclear
segmentation/ObjectID namespace.

With:

```python
ALLOW_NESTED_NUCLEAR_OBJECT_SETS = True
```

the script permits:

```text
identical nuclear ObjectID sets
```

or:

```text
one nuclear set being a strict subset of the other
```

in which case the union/larger set is used.

Partially different, non-nested nuclear ObjectID sets cause the script to
stop because they suggest incompatible nuclear segmentations or ID
namespaces.

Therefore, preserve the original marker↔nucleus HDF5 files. They are not only
QC records; they provide the canonical nuclear links required for the later
four-class all-cell summary.

Before accepting a colocalization result, inspect both the visual QC and the
relevant summary CSVs.


---

---

# Pipeline-wide reference

---

## The complete eight-step flow

---

```text
ENV
│
├── create/activate ims-mosaic
├── create/activate cellprofiler-native when needed
└── build/locate Fiji + BigStitcher SIF
     │
     v
FEED INTO → IMS LOADING
│
├── STEP 1 — .ims processing → stitched/max-projected images
│      raw F00–F03 IMS
│      → logical TIFF export
│      → zero-padding trim
│      → BaSiC
│      → BigStitcher registration
│      → overlap-based intensity adjustment
│      → 32-bit fusion
│      → max projection
│
├── STEP 2 — initial ROI segmentation per marker
│      unaligned max projection
│      → CellProfiler mask + centroid measurements
│
└── STEP 3 — image alignment using ROI/centroid information
       every marker independently → nuclear reference
       │
       v
ROI SCRUBBING
│
├── STEP 4 — resegment aligned images
│      aligned TIFF
│      → aligned ROI pixel table
│
├── STEP 5 — intensity visualization + filtering
│      aligned ROI table
│      → intensity QC
│      → intensity-filtered ROI table
│
└── STEP 6 — shape-based filtering
       intensity-filtered ROI table
       → shape/area QC
       → final filtered ROI table
       │
       v
ROI → CELL
│
└── STEP 7 — nuclear colocalization
       final marker ROI + final nuclear ROI
       → marker↔nucleus HDF5
       → nuclear-positive marker pixels
       → canonical nuclear ObjectID mapping
       │
       v
PUTATIVE CELL TYPES
│
└── STEP 8 — marker colocalization
       nuclear-positive marker A ↔ nuclear-positive marker B
       → pairwise marker relationships
       → QC / summaries
       → optional all-nuclear-cell four-class table
```

---

## Which steps happen once versus per channel

---

### Once per software setup

```text
Create ims-mosaic Conda environment
Create cellprofiler-native Conda environment
Build Fiji + BigStitcher SIF
```

### Once per four-field acquisition

```text
Run IMS mosaic pipeline
```

The IMS pipeline automatically processes all channels in the acquisition.

### Once per channel before alignment

```text
CellProfiler pass 1
```

### Once per non-nuclear channel

```text
CASTalign to nuclear reference
```

### Once per aligned channel

```text
CellProfiler pass 2
ROI extraction
intensity filtering
shape filtering
```

### Once per marker

```text
marker-versus-nuclear colocalization
colocalization QC
export nuclear-positive marker pixels
```

### Once per desired marker pair

```text
marker-versus-marker colocalization
colocalization QC
population summary
optional organoid-area summary
all-nuclear-cell four-class summary
```

The all-nuclear-cell summary additionally reuses the two corresponding
marker-versus-nuclear HDF5 files from the first colocalization stage.

---

---

## Suggested project directory

---

A practical analysis directory could be:

```text
Experiment/
├── 00_raw_ims/
│   ├── Sample_F00.ims
│   ├── Sample_F01.ims
│   ├── Sample_F02.ims
│   └── Sample_F03.ims
│
├── 01_mosaic_pipeline/
│   └── ... IMS_to_TIFF outputs ...
│
├── 02_cellprofiler_prealignment/
│   ├── Nuclear/
│   ├── MarkerA/
│   ├── MarkerB/
│   └── MarkerC/
│
├── 03_aligned/
│   ├── Nuclear/
│   ├── MarkerA/
│   ├── MarkerB/
│   └── MarkerC/
│
├── 04_cellprofiler_aligned/
│   ├── Nuclear/
│   ├── MarkerA/
│   ├── MarkerB/
│   └── MarkerC/
│
├── 05_aligned_roi_pixels/
│   ├── Nuclear/
│   ├── MarkerA/
│   ├── MarkerB/
│   └── MarkerC/
│
├── 06_filtering/
│   ├── 01_intensity_reports/
│   ├── 02_intensity_filter/
│   ├── 03_shape_reports_after_intensity/
│   └── 04_final_shape_filter/
│
├── 07_nuclear_colocalization/
│   ├── MarkerA_vs_Nuclear/
│   ├── MarkerB_vs_Nuclear/
│   └── MarkerC_vs_Nuclear/
│
├── 08_nuclear_positive_marker_cells/
│   ├── MarkerA_nuclear_positive_cells.csv
│   ├── MarkerB_nuclear_positive_cells.csv
│   └── MarkerC_nuclear_positive_cells.csv
│
└── 09_marker_marker_colocalization/
    ├── MarkerA_vs_MarkerB/
    │   └── visualizations/
    │       └── unmatched_only_qc/
    │           ├── *_final_summary_metrics.csv
    │           ├── *_all_nuclear_cells_classification.csv
    │           └── *_all_nuclear_cells_class_counts.csv
    ├── MarkerA_vs_MarkerC/
    └── MarkerB_vs_MarkerC/
```

The exact directory names are optional, but keeping each stage separate helps prevent accidentally feeding pre-alignment, unfiltered, or first-pass-filtered objects into a later stage.

---

---

## Files that should never be overwritten accidentally

---

Preserve separately:

```text
raw IMS files

unaligned max projections

CellProfiler pass-1 measurements

aligned TIFFs

CellProfiler pass-2 measurements

unfiltered aligned ROI pixel CSVs

intensity-filtered ROI pixel CSVs

final intensity+shape-filtered ROI pixel CSVs

marker-versus-nuclear HDF5 files

marker-versus-marker HDF5 files

final summary-metrics CSVs

all-nuclear-cell classification/count CSVs

all visualization/QC directories
```

These represent different biological/processing stages and should not share ambiguous filenames in one directory.

---

---

## QC checkpoints

---

Do not treat the pipeline as one completely blind automated chain.

Recommended manual QC checkpoints are:

```text
after stitching
    inspect seams and gross spatial failures

after cross-channel alignment
    inspect nuclear/marker overlays and residuals

after CellProfiler segmentation
    inspect ROI masks

after intensity report
    select/justify intensity filtering

after intensity-only filtering
    verify low-intensity segmentation artifacts were removed

after shape reports
    choose shape thresholds

after final filtering
    inspect retained/removed populations

after marker-versus-nuclear colocalization
    inspect matched and unmatched objects

after marker-versus-marker colocalization
    inspect representative matched pairs

after all-nuclear-cell classification
    verify the two source marker↔nucleus HDF5 files use a compatible
    canonical nuclear ObjectID set and inspect the four-class counts
```

---

---

## Key implementation corrections

---

### Correction 1 — stitching is not run separately per channel

Run the IMS mosaic script once on the four-field acquisition.

It automatically loops over all available channels.

### Correction 2 — the IMS script already creates max projections

The standalone max-projection script is optional when using the current IMS mosaic pipeline.

### Correction 3 — ROI extraction before alignment is optional

The alignment script consumes CellProfiler centroid CSVs, not ROI-pixel CSVs.

Aligned ROI extraction is the critical extraction stage for filtering/colocalization.

### Correction 4 — intensity must remain logically upstream of shape filtering

If the intensity filter is disabled for the second pass, the second pass must read the **intensity-filtered CSV**, not the original ROI CSV.

### Correction 5 — current CellIDs are not globally stable

The current colocalization script assigns sequential CellIDs independently for each pairwise run.

Use the nuclear ObjectID mapping when a persistent cross-marker cell identity is required.

### Correction 6 — nuclear-positive marker pixels are stored in HDF5

The first colocalization script does not automatically emit those pixels as a standalone downstream CSV.

Export:

```text
tables/filtered_object_a_pixels
```

before feeding nuclear-positive marker cells into another colocalization run.

---

### Correction 7 — CellProfiler uses a separate environment and has two execution modes

The main pipeline runs in:

```text
ims-mosaic
```

CellProfiler runs in:

```text
cellprofiler-native
```

At each CellProfiler stage:

```text
ims-mosaic
→
cellprofiler-native
→
run CellProfiler
→
ims-mosaic
```

CellProfiler may be run either:

```text
GUI: MGEOPV.cpproj
```

or:

```text
Headless: MGEOPVFinal.cppipe + cellprofilerdriver.py
```

Both are valid execution options for the CellProfiler stages of the same larger workflow.

---


### Correction 8 — the current colocalization visualization script is also a summary stage

`02_generate_colocalization_visualizations.py` is no longer only a plotting
script.

In addition to visual QC, the current version can produce:

```text
final marker-population class counts
final summary metrics
marker-positive area as % of organoid area
per-nucleus A/B classification
all-nuclear-cell four-class counts and percentages
```

For a true all-cell denominator, it must reuse the original marker↔nucleus
HDF5 mappings rather than infer nuclear identity from the later pairwise
`CellID` values.

---

### Correction 9 — IMS reruns clean the previous generated sample result by default

The current IMS mosaic script removes the previous pipeline-generated
directory for the same acquisition before rerunning unless:

```bash
--keep-previous
```

is supplied.

This cleanup affects only pipeline-generated output under the
`*_IMS_to_TIFF/<sample>/` result directory; it does not delete the original
source `.ims` files.

This behavior is intentional because stale BigStitcher aliases, XML/HDF5
state, corrected TIFFs, or projections can otherwise contaminate a rerun.

---

## Final interpretation

---

The pipeline has four useful conceptual identity levels:

```text
1. Segmentation ROI
   independent object in one channel

2. Nuclear-associated marker cell
   marker ROI that successfully matches a nucleus

3. Pairwise multi-marker relationship
   nuclear-associated marker cells that overlap/match another
   nuclear-associated marker population

4. Canonical all-cell classification
   each nuclear ObjectID classified as A+/B+, A+/B-, A-/B+, or A-/B-
   by reusing the original marker↔nucleus assignments
```

Keeping these levels separate prevents channel-specific ROI labels or
run-local `CellID` values from being mistaken for universal biological cell
identities.

The nuclear reference provides the shared spatial coordinate system, while
the canonical nuclear `ObjectID` provides the persistent identity anchor
across separate marker analyses. The final all-cell summary should therefore
be built from those preserved nuclear assignments, not from numerical
equality of pairwise `CellID` values.

---

# Eight-step order — final reminder

The pipeline should be read and executed in this order:

```text
IMS LOADING
1. .ims processing → stitched image
2. Initial segmentation of ROI per marker
3. Align images with ROI information

ROI SCRUBBING
4. Resegment aligned images
5. Apply intensity visualization + filter
6. Apply shape-based filter

ROI → CELL
7. Apply nuclear colocalization for ROI → Cell identity

PUTATIVE CELL TYPES
8. Marker colocalization
```

The phase labels organize the biology and image-processing logic; the numbered steps define the actual execution order.
