# General microscopy pipeline

`GeneralPipeline/` is the experiment-independent interface to the repository's
existing microscopy analysis. It supports two starting points and sends both to
the same committed downstream controller, `run_stitched_pipeline.py`.

```text
raw_ims_tiles                         stitched_tiffs
      |                                     |
discover complete 2x2 fields         validate exact sample/channel mapping
validate IMS channels/calibration    validate YX or ZYX single-channel TIFFs
IMS export -> BaSiC -> BigStitcher    maximum projection when ZYX
      |                                     |
      +---------- run_stitched_pipeline.py -+
                         |
 initial segmentation -> CASTalign -> final registered segmentation
 -> shared nuclear-reference tissue mask -> ROI/filtering/matching/QC
```

The general wrappers orchestrate existing scripts. They do not implement a
second copy of registration, segmentation, filtering, or colocalization.

## Installation and runtime

Use the repository's existing `ims-mosaic` and `cellprofiler-native` Conda
environments. Raw tiled IMS processing additionally needs the existing
Fiji/BigStitcher Singularity image. On Great Lakes the supplied Slurm wrapper
loads `singularity/4.4.1`, binds `/gpfs` and `/scratch`, activates
`ims-mosaic`, and invokes CellProfiler through the configured
`cellprofiler-native` environment.

Run commands from a clean, committed checkout. `prepare` records the Git
revision and scientific-file hashes and refuses a dirty checkout by default.
`--allow-dirty-repository` is only for local development.

## Route A: original tiled IMS acquisitions

Copy and edit the example:

```bash
cd /path/to/ImageAnalysisScripts-Copy
cp GeneralPipeline/configs/raw_ims.example.json raw_ims.json
```

At minimum, replace `preprocessing.fiji_sif`. Confirm all scientific settings;
the values in an example are not automatically valid for a new experiment.

Supported source layouts include a field directly below the input root or any
nested experiment/date/session layers up to `raw_ims_tiles.max_depth`:

```text
INPUT_ROOT/
  Experiment_A/                    # optional; no Olti_* name is required
    2026-06-18/                    # optional
      Sample_Set3_DAPI_GFP_Field_1/
        Sample_Set3_DAPI_GFP_F02.ims
        Sample_Set3_DAPI_GFP_F04.ims
        Sample_Set3_DAPI_GFP_F06.ims
        Sample_Set3_DAPI_GFP_F08.ims
```

One additional tile-directory layer is also supported, for example
`tile_F02/..._F02.ims`. Original tile numbers need not start at zero. The four
tiles are sorted numerically and staged as `F00` through `F03` symlinks in the
workspace; the manifest retains the original numbers and paths.

The validated grid is 2×2. A different grid is rejected instead of being
advertised without corresponding stitching validation. Incomplete fields may
be skipped and recorded. Duplicate indices, extra tiles, unnumbered IMS files,
malformed metadata, or ambiguous channels are fatal.

Biological channel identity comes from either:

- a recognized suffix in every IMS filename when `channel_names` is `auto`; or
- an explicit `channel_names` list in physical IMS index order.

`Set3` is stored as experiment metadata, not a channel. Marker digits and
underscores in names such as `ASCL1`, `BiVe3`, and `NKX2_1` are preserved. An
optional standalone acquisition number after the final marker is supported,
but its directory and filename values must agree. Stored labels such as
`Confocal - Green` are provenance only and never become a biological mapping.

Every tile must have the same physical channel indices, stored names,
wavelength metadata, ZYX dimensions, and voxel calibration. A physical phase
channel or a physical/biological count mismatch is rejected. Phase exclusion
is not implemented because silently deleting an index could shift every
downstream biological channel.

Run the lifecycle:

```bash
REPO=/path/to/ImageAnalysisScripts-Copy
INPUT=/path/to/experiments
WORKSPACE=/scratch/parent_root/parent0/$USER/my_experiment_001

python "$REPO/GeneralPipeline/run_pipeline.py" scan \
  --mode raw_ims_tiles --input "$INPUT" --config "$REPO/raw_ims.json"

python "$REPO/GeneralPipeline/run_pipeline.py" prepare \
  --mode raw_ims_tiles --input "$INPUT" --config "$REPO/raw_ims.json" \
  --workspace "$WORKSPACE"

python "$REPO/GeneralPipeline/run_pipeline.py" validate \
  --manifest "$WORKSPACE/manifest.json"

python "$REPO/GeneralPipeline/run_pipeline.py" preflight \
  --manifest "$WORKSPACE/manifest.json" --index 0

python "$REPO/GeneralPipeline/run_pipeline.py" submit \
  --manifest "$WORKSPACE/manifest.json" --max-parallel 2 --print-only

# Remove --print-only only after reviewing the command.
python "$REPO/GeneralPipeline/run_pipeline.py" submit \
  --manifest "$WORKSPACE/manifest.json" --max-parallel 2
```

To run locally instead of Slurm:

```bash
python "$REPO/GeneralPipeline/run_pipeline.py" run \
  --manifest "$WORKSPACE/manifest.json" --index 0
```

After an interrupted task, use the explicit resume command only with the same
manifest, inputs, repository revision, and parameters:

```bash
python "$REPO/GeneralPipeline/run_pipeline.py" resume \
  --manifest "$WORKSPACE/manifest.json" --index 0

python "$REPO/GeneralPipeline/run_pipeline.py" inspect \
  --manifest "$WORKSPACE/manifest.json"
```

## Route B: TIFFs already stitched in Imaris

This route accepts one channel-specific, single-channel TIFF per biological
channel. Each input may be a 2-D `YX` image or a 3-D `ZYX` page stack. The
existing downstream runner makes a maximum projection for a supported stack.
RGB, time series, multi-position files, and other hyperstacks are rejected.

Copy and edit the example:

```bash
cp GeneralPipeline/configs/stitched_tiff.example.json stitched_tiffs.json
```

Replace every TIFF path and both pixel sizes. The pipeline deliberately
requires explicit positive `pixel_size_um.x` and `.y`; it does not substitute a
historical or arbitrary calibration. TIFF paths in `mapping.samples` are
relative to `--input` unless absolute.

The example contains a single explicit sample. For a separate explicit sample
manifest, replace the mapping with:

```json
{"sample_manifest": "samples.json"}
```

`samples.json` is resolved relative to the configuration file and has this
shape:

```json
{
  "schema": 1,
  "samples": [
    {
      "name": "sample_001",
      "channels": {
        "DAPI": "sample_001_DAPI.tif",
        "GFP": "sample_001_GFP.tif",
        "LHX6": "sample_001_LHX6.tif"
      }
    }
  ]
}
```

Channel paths inside the sample manifest are still resolved relative to
`--input`. Every sample must contain exactly the configured channels.

For deterministic filename discovery, use one anchored rule with named
`sample` and `channel` groups, for example:

```json
{
  "filename_regex": "(?P<sample>.+)__(?P<channel>DAPI|GFP|LHX6)\\.tif"
}
```

The rule is applied as a full filename match. Missing and duplicate channels
stop preparation. There is no fuzzy matching or arbitrary-sort pairing.

Run the same lifecycle with a different mode:

```bash
REPO=/path/to/ImageAnalysisScripts-Copy
INPUT=/path/to/imaris_stitched_tiffs
WORKSPACE=/scratch/parent_root/parent0/$USER/my_stitched_experiment_001

python "$REPO/GeneralPipeline/run_pipeline.py" scan \
  --mode stitched_tiffs --input "$INPUT" --config "$REPO/stitched_tiffs.json"

python "$REPO/GeneralPipeline/run_pipeline.py" prepare \
  --mode stitched_tiffs --input "$INPUT" --config "$REPO/stitched_tiffs.json" \
  --workspace "$WORKSPACE"

python "$REPO/GeneralPipeline/run_pipeline.py" validate \
  --manifest "$WORKSPACE/manifest.json"

python "$REPO/GeneralPipeline/run_pipeline.py" preflight \
  --manifest "$WORKSPACE/manifest.json" --index 0

python "$REPO/GeneralPipeline/run_pipeline.py" submit \
  --manifest "$WORKSPACE/manifest.json" --max-parallel 2 --print-only
```

This mode bypasses IMS export, BaSiC, and BigStitcher entirely. If native
channel canvases differ, alignment must remain enabled. CASTalign registers
each marker onto the nuclear reference canvas before final segmentation.

## Commands and resume behavior

- `scan` is read-only and reports valid acquisitions and skipped incomplete raw
  fields.
- `prepare` requires a new/empty workspace, stages read-only symlinks, writes
  per-acquisition configs, and atomically writes `manifest.json`.
- `validate` checks source stats, symlinks, repository revision, configuration,
  generic runner, Slurm runner, and every scientific-file hash.
- `preflight` reopens one raw field's IMS metadata or validates one stitched
  configuration without running scientific stages.
- `run` executes one acquisition and refuses tracked downstream outputs unless
  the underlying runner is explicitly resumed.
- `resume` executes one acquisition with provenance-checked resume. The raw
  preprocessing path uses the existing keep/validation behavior. IMS
  export/stitching and general-pipeline area-mean TIFF intermediates retain the
  repository's atomic validated-write behavior.
- `inspect` reports whether each downstream run manifest exists.
- `submit` creates one Slurm array task per acquisition. `--print-only` never
  calls `sbatch`. A recorded submission blocks accidental duplicates.

All job ordering is deterministic. A resume is refused after source file
statistics, configuration, Git revision, runner, or scientific code changes.
Use a fresh workspace for changed inputs or parameters.

## Workspace layout

```text
WORKSPACE/
  manifest.json
  configs/                  # immutable per-acquisition inputs
  input_staging/            # symlinks only; original sources remain untouched
  controller_compat/        # raw IMS only; preserves staged F00-F03 aliases
  logs/
  outputs/
    ACQUISITION/
      preprocessing/        # raw IMS only
      general_pipeline_provenance.json
      analysis/             # shared downstream outputs and run_manifest.json
```

For sources below the protected Great Lakes Turbo acquisition root, the
workspace must resolve below `/scratch`. The source tree is never used for
intermediates, logs, manifests, or outputs. Existing nonempty workspaces are
never overwritten.

## Shared downstream behavior

Both modes ultimately invoke the same top-level file:

```text
run_stitched_pipeline.py
```

It retains the established order: initial CellProfiler segmentation,
CASTalign, final segmentation in registered coordinates, one final nuclear
tissue mask shared by marker segmentations, ROI extraction, intensity and
shape filtering, nuclear-anchored matching, marker-pair colocalization, and QC.
The generic wrapper does not alter those algorithms or numerical settings.

## Already-stitched IMS limitation

`run_ims_pipeline.py` remains available for an IMS file already stitched in
Imaris. That controller currently assumes exactly four reporters in the fixed
Hoechst, mNeonGreen, BiVe3/BiVe4, PV layout. It is therefore preserved and
documented in `IMS_WORKFLOW.md`, but it is not presented as a generalized
`GeneralPipeline` input mode. Use it only when that fixed mapping is correct.

## Known limitations

- Raw tiled input is validated only for a 2×2/four-tile acquisition.
- Physical phase-channel exclusion is unsupported; such files are rejected.
- TIFF calibration must be explicit even when TIFF tags appear to contain a
  resolution, because no tag convention has yet been validated for this
  workflow.
- TIFF input supports single-channel `YX` and simple `ZYX`/page stacks only.
- A filename discovery rule must provide exact named sample/channel groups.
- The general wrapper is additive. The active IGF production wrappers retain
  their original implementation and paths; common-code migration is deferred
  until it can be compatibility-tested against completed production manifests.

For the lower-level handoff contract, see `../SCRIPT_HANDOFF_REFERENCE.md`. For
the direct downstream controller, see `../STITCHED_WORKFLOW.md`.
