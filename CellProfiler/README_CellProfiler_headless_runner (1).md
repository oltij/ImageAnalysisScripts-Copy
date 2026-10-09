# Headless CellProfiler Single-Image Runner

## Purpose

This script runs a CellProfiler `.cppipe` pipeline headlessly on **one fluorescence image**, with an optional aligned organoid-mask image, and collects all outputs into one specified directory.

It is designed for:

```text
CellProfiler pipeline
+
one microscopy image
+
optional aligned organoid mask
+
one output directory
```

The script also protects against a common CellProfiler behavior: some `SaveImages` or `ExportToSpreadsheet` modules may write files beside the input image rather than only into the command-line output directory. To avoid losing those files, the script uses a temporary input directory, symlinks the image into it, runs CellProfiler, then copies any generated files from that temporary directory into the final output directory before cleanup.

## Suggested script name

```text
cellprofilerdriver.py
```

or:

```text
run_cellprofiler_headless.py
```

The examples below use:

```text
cellprofilerdriver.py
```

## Dependencies

The Python imports are all from the standard library:

```text
argparse
subprocess
tempfile
shutil
pathlib
```

The external requirement is that the `cellprofiler` executable is available in the active environment.

For example:

```bash
conda activate cellprofiler-native
```

Then verify:

```bash
which cellprofiler
```

and optionally:

```bash
cellprofiler --version
```


## CellProfiler can be run in two ways

There are two supported ways to run the CellProfiler portion of this workflow.

### Option 1 — CellProfiler GUI

Use the CellProfiler graphical interface with the project:

```text
MGEOPV.cpproj
```

This option is useful when you want to:

```text
inspect images interactively
review segmentation
adjust thresholds or module settings
step through the pipeline
visually confirm that the ROIs are being identified correctly
```

In this mode, open CellProfiler normally and load:

```text
MGEOPV.cpproj
```

Then select or load the appropriate input image and run the project through the GUI.

### Option 2 — headless CellProfiler

Use the headless runner documented in this README together with:

```text
MGEOPVFinal.cppipe
```

The checked-in `MGEOPVFinal.cppipe` is a genuine text pipeline. It was exported unchanged from the pipeline embedded in the previously project-formatted file so CellProfiler's `-p` command-line option can load it reliably. Do not substitute the binary `MGEOPV.cpproj` in a headless command.

Example:

```bash
python /home/oltij/Desktop/cellprofilerdriver.py \
    --pipeline "/home/oltij/Desktop/MGEOPVFinal.cppipe" \
    --input "/path/to/input_image.tif" \
    --output "/path/to/CellProfiler_output"
```

The headless option is useful when:

```text
the CellProfiler settings have already been finalized
you want reproducible command-line execution
you are running on Great Lakes
you do not need to interactively inspect each module during the run
```

The two options are therefore:

```text
GUI:
    MGEOPV.cpproj

Headless:
    MGEOPVFinal.cppipe
    +
    cellprofilerdriver.py
```

## Important: switch Conda environments before running CellProfiler

The broader image-analysis pipeline commonly runs in:

```text
ims-mosaic
```

but CellProfiler is run from:

```text
cellprofiler-native
```

Therefore, if you are currently in:

```text
(ims-mosaic)
```

switch environments before running CellProfiler.

A safe sequence is:

```bash
conda deactivate
conda activate cellprofiler-native
```

Your prompt should then begin with something similar to:

```text
(cellprofiler-native)
```

Run CellProfiler either through the GUI or with the headless runner.

After the CellProfiler stage is finished, leave the CellProfiler environment:

```bash
conda deactivate
```

and return to the main image-analysis environment:

```bash
conda activate ims-mosaic
```

Your prompt should again begin with:

```text
(ims-mosaic)
```

Conceptually:

```text
ims-mosaic
    |
    | switch before CellProfiler
    v
cellprofiler-native
    |
    | run MGEOPV.cpproj in GUI
    | OR
    | run MGEOPVFinal.cppipe headlessly
    |
    | switch back afterward
    v
ims-mosaic
```

This environment switch should be performed each time the workflow enters or leaves a CellProfiler stage.

## Command-line arguments

The script requires:

```text
--pipeline
--input
--output
```

It also accepts optional `--organoid-mask`. This is used by the automated
stitched-image controller for final marker segmentation and requires the
controller-generated two-input pipeline.

### `--pipeline`

Path to a CellProfiler pipeline:

```text
.cppipe
```

Example:

```text
/home/oltij/Desktop/MGEOPV.cppipe
```

The script resolves the path and checks that it exists. If not, it raises `FileNotFoundError`.

### `--input`

Path to one input image.

Typical example:

```text
/home/oltij/Desktop/PVImage/Exp17_157-2_Set2_DAPI_GFP_LHX6_PV_FusionStitcher_F0_Confocal - Far Red_MaxProjection.tif
```

The script itself does not inspect image dimensions or dtype. The input simply needs to be readable by the CellProfiler pipeline.

### `--output`

Destination directory for all CellProfiler outputs.

Example:

```text
/home/oltij/Desktop/CellProfiler_output
```

The script creates the directory automatically if needed.

### `--organoid-mask`

Optional path to the final-coordinate Hoechst `OrganoidMask.tiff`. When
supplied, the driver links it beside the fluorescence input under the stable
name expected by the generated CellProfiler pipeline. The initial
pre-registration passes do not use this argument.

## Basic usage

```bash
python cellprofilerdriver.py \
    --pipeline "/path/to/pipeline.cppipe" \
    --input "/path/to/input_image.tif" \
    --output "/path/to/output_directory"
```

## Example

```bash
conda activate cellprofiler-native

python /home/oltij/Desktop/cellprofilerdriver.py \
    --pipeline "/home/oltij/Desktop/MGEOPVFinal.cppipe" \
    --input "/home/oltij/Desktop/PVImage/Exp17_157-2_Set2_DAPI_GFP_LHX6_PV_FusionStitcher_F0_Confocal - Far Red_MaxProjection.tif" \
    --output "/home/oltij/Desktop/CellProfiler_output"
```

Quotes are important when paths contain spaces.

## Workflow

```text
pipeline.cppipe
+
one input image
+
requested output directory
        |
        v
validate input paths
        |
        v
create temporary input directory
        |
        v
symlink input image into temporary directory
        |
        v
run CellProfiler headlessly
        |
        v
collect files written to requested output directory
and files written beside temporary input
        |
        v
copy recovered files into final output directory
        |
        v
delete temporary directory
```

## Temporary input directory

The script creates:

```python
with tempfile.TemporaryDirectory(
    prefix="cellprofiler_input_"
) as temp_dir_string:
```

A temporary directory may look like:

```text
/tmp/cellprofiler_input_abcd1234/
```

The exact location depends on the operating system and environment.

## Why a temporary directory is used

The CellProfiler `-i` option takes an input directory, not a single image path.

This script is intended to process exactly one image. Therefore it creates a temporary directory containing only that image.

This avoids accidentally processing other images that might be present in the original image's directory.

## The input TIFF is symlinked, not copied

The script creates:

```python
temporary_input.symlink_to(input_file)
```

Conceptually:

```text
/tmp/cellprofiler_input_xxxx/input.tif
        |
        | symbolic link
        v
/home/oltij/.../actual_input.tif
```

This avoids copying a potentially very large microscopy image.

The original input file is not modified.

## CellProfiler command executed

The command is:

```python
command = [
    "cellprofiler",
    "-c",
    "-r",
    "-p", str(pipeline_path),
    "-i", str(temp_dir),
    "-o", str(output_dir),
    "-L", "INFO",
]
```

Equivalent shell command:

```bash
cellprofiler \
    -c \
    -r \
    -p "/path/to/pipeline.cppipe" \
    -i "/tmp/cellprofiler_input_xxxx" \
    -o "/path/to/output_directory" \
    -L INFO
```

### `-c`

Run CellProfiler without the graphical user interface.

### `-r`

Run the pipeline.

### `-p`

Specify the `.cppipe` pipeline.

### `-i`

Specify the temporary single-image input directory.

### `-o`

Specify the final output directory.

### `-L INFO`

Use INFO-level CellProfiler logging.

## Running CellProfiler

The script executes:

```python
result = subprocess.run(command)
```

If CellProfiler returns:

```text
0
```

the script continues.

If CellProfiler returns a nonzero status, the script raises:

```text
RuntimeError
```

for example:

```text
CellProfiler exited with code 1
```

## Output recovery

Some pipelines may be configured so that `SaveImages` or `ExportToSpreadsheet` writes beside the input image.

Because the input directory in this script is temporary, those files would disappear unless copied elsewhere.

After CellProfiler finishes, the script scans:

```python
temp_dir.iterdir()
```

and skips only the original input-image symlink.

Anything else is treated as CellProfiler output.

## Recovering ordinary files

Files are copied with:

```python
shutil.copy2(
    item,
    destination
)
```

The destination is:

```text
<output directory>/<same filename>
```

## Recovering directories

Directories are copied recursively.

If a directory with the same name already exists in the output directory, the existing directory is removed first:

```python
shutil.rmtree(destination)
```

and then replaced with:

```python
shutil.copytree(item, destination)
```

## Typical CellProfiler outputs

Exact outputs depend on the `.cppipe`.

For the microscopy workflow in this repository, common outputs may include:

```text
CellMask.tiff
MyExpt_FilterObjects.csv
MyExpt_FilterObjects2.csv
```

Other pipelines may generate:

```text
CSV files
TIFF files
PNG files
subdirectories
other CellProfiler exports
```

The driver does not require fixed output names.

## Final report

After CellProfiler completes and any temporary outputs are recovered, the script prints:

```text
CELLPROFILER COMPLETED SUCCESSFULLY
```

followed by the final output directory and every item currently present in it.

If no files exist, it prints:

```text
[NO FILES FOUND]
```

If files were recovered from the temporary input directory, it also reports the number of recovered items.

## Temporary-directory cleanup

The temporary input directory is deleted automatically when the `TemporaryDirectory` context exits.

The symlink is therefore removed automatically.

The original microscopy image remains unchanged.

## Use in the larger microscopy pipeline

This driver can be used for both CellProfiler stages.

### Before alignment

Run CellProfiler on the stitched/max-projected channels to obtain object centroids and masks for registration.

Conceptually:

```text
max-projected channel TIFF
        |
        v
CellProfiler
        |
        +--> object CSV
        +--> mask TIFF
```

The centroid CSV can then be used by the registration script.

### After alignment

The automated controller first runs CellProfiler on final-frame Hoechst to
create the reference organoid mask. It then runs every aligned marker with
that mask as the optional second input. Marker fluorescence still drives cell
detection; the Hoechst mask supplies the common tissue boundary.

Conceptually:

```text
aligned channel TIFF
        |
        v
CellProfiler
        |
        +--> aligned object CSVs
        +--> aligned mask TIFF
        |
        v
ROI extraction
```

These post-alignment ROIs are used for downstream filtering and colocalization.

## One fluorescence image per invocation

The script is intentionally designed to process one fluorescence image at a
time. Final marker runs may additionally receive one shared mask.

For multiple channels, run it once per channel.

Example:

```bash
python cellprofilerdriver.py \
    --pipeline "/path/to/Hoescht.cppipe" \
    --input "/path/to/aligned_Hoechst_max.tif" \
    --output "/path/to/HoeschtAligned/CellProfiler_output"
```

```bash
python cellprofilerdriver.py \
    --pipeline "/path/to/LHX6.cppipe" \
    --input "/path/to/aligned_LHX6_max.tif" \
    --output "/path/to/LHX6Aligned/CellProfiler_output"
```

```bash
python cellprofilerdriver.py \
    --pipeline "/path/to/PV.cppipe" \
    --input "/path/to/aligned_PV_max.tif" \
    --output "/path/to/PVAligned/CellProfiler_output"
```

Use the correct `.cppipe` for each channel if different channel-specific pipelines are used.

## Important assumptions

1. `cellprofiler` is available in the active shell environment.
2. The `.cppipe` file is already configured correctly.
3. The script processes one image per run.
4. The input image is symlinked rather than copied.
5. The filesystem must support symbolic links.
6. CellProfiler must be able to follow the symlink.
7. The output directory is not cleared before the run.
8. Existing ordinary files with the same name may be overwritten during recovery.
9. Existing directories with the same name are removed and replaced during recovery.
10. Exact output filenames depend on the `.cppipe`.
11. The script does not validate image dimensions, dtype, or channel identity.
12. A nonzero CellProfiler exit code causes the script to stop.

## Common error: `cellprofiler` not found

Activate the CellProfiler environment:

```bash
conda activate cellprofiler-native
```

Then verify:

```bash
which cellprofiler
```

## Common error: pipeline not found

Check:

```bash
ls -lh "/path/to/pipeline.cppipe"
```

Then correct the `--pipeline` path.

## Common error: input image not found

Check:

```bash
ls -lh "/path/to/input_image.tif"
```

Make sure paths containing spaces are quoted.

## Common error: CellProfiler exits nonzero

If you see:

```text
RuntimeError:
CellProfiler exited with code 1
```

inspect the CellProfiler terminal output above the Python exception. The underlying problem is usually inside the CellProfiler pipeline or one of its modules.

## Which CellProfiler file should I use?

Use:

```text
MGEOPV.cpproj
```

for the graphical CellProfiler workflow.

Use:

```text
MGEOPVFinal.cppipe
```

for the headless command-line workflow described by this script.

The `.cpproj` and `.cppipe` files serve different usage modes and should not be treated as interchangeable command-line inputs.

The stitched-image controller also checks the file signature and required output settings before starting; a project/HDF5 file renamed with a `.cppipe` extension is rejected with an actionable error.

## Minimal command

```bash
conda activate cellprofiler-native

python /home/oltij/Desktop/cellprofilerdriver.py \
    --pipeline "/home/oltij/Desktop/MGEOPVFinal.cppipe" \
    --input "/path/to/input_image.tif" \
    --output "/path/to/CellProfiler_output"
```

## Data-flow summary

```text
.cppipe pipeline
       +
one TIFF/image
       |
       v
temporary single-image directory
       |
       v
CellProfiler headless
       |
       ├── outputs written to -o directory
       │
       └── outputs written beside temporary input
                  |
                  v
               recovered
                  |
                  v
       final requested output directory
```
