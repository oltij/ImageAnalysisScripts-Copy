# IGF organoid IMS folder batch

This directory packages the working Great Lakes folder-batch controller around
the repository's existing scientific pipeline. It starts with the original
four-tile `.ims` acquisitions and runs the established sequence:

1. IMS export
2. BaSiC correction
3. BigStitcher mosaic assembly and maximum projections
4. initial CellProfiler segmentation for registration anchors
5. CASTalign registration
6. segmentation again in the aligned coordinate system
7. shared nuclear-tissue-mask QC, ROI filtering, and colocalization

The wrapper does not replace or optimize those algorithms. It calls the
existing repository controllers and preserves their settings.

## Files

- `run_folder_batch.py`: recursive discovery, validation, symlink staging,
  manifests, per-field execution, inspection, and resume support.
- `submit_folder_batch.py`: one-command preparation and Slurm array submission.
- `run_folder_batch.sbatch`: the working Great Lakes CPU job profile.
- `configs/ims_mge_generic.example.json`: portable batch profile derived from
  the supplied working configuration.

All commands may be launched from any current directory. The scripts determine
the repository root from their own absolute paths. `submit_folder_batch.py`
passes the absolute `run_folder_batch.py` path into the Slurm job because `$0`
inside a job refers to Slurm's spooled copy, not reliably to the checkout.

## Source and scratch safety

For production, keep the two roots separate:

```bash
IMS_INPUT=/nfs/turbo/umms-parent/andor_microscope_data_dump/path/to/session
IGF_WORKSPACE=/scratch/parent_root/parent0/$USER/igf_organoids/run_001
```

The source tree is scanned, stat-checked, and read, but never written. The
wrapper creates only symbolic links to original IMS files. The links,
`controller_compat`, generated configs, manifests, preprocessing, logs, and all
analysis results are created below `IGF_WORKSPACE`. Preparation rejects a
workspace inside the source tree and rejects reuse of a nonempty workspace.

For inputs below the protected Turbo root, preparation rejects any workspace
that is not below `/scratch`. Each scientific run should use a new scratch
workspace. Original `.ims` files under
`/nfs/turbo/umms-parent/andor_microscope_data_dump/` remain untouched.

## Great Lakes prerequisites

Use the existing environments:

```bash
conda env create -f CondaEnvironment/ims-mosaic.yml
conda env create -f CondaEnvironment/cellprofiler-native.yml
conda activate ims-mosaic
```

The Slurm runner preserves the working settings:

- `singularity/4.4.1`
- `SINGULARITY_BINDPATH=/gpfs:/gpfs,/scratch:/scratch`
- `ims-mosaic` for orchestration, IMS export, and stitching
- `cellprofiler-native` through the downstream CellProfiler configuration
- four CPUs, 64 GB memory, standard CPU partition, and a one-day time limit

No GPU is requested.

## Acquisition layout and naming

The input may be a complete session, a group directory, or one `Field_N`
directory. Discovery is recursive (five levels by default) and supports an
optional date-folder layer. Each complete field must contain four IMS tiles,
either directly or in tile directories. Derived `FusionStitcher` IMS files are
ignored.

Both naming styles are supported:

```text
Sample_DAPI_GFP_LHX6_PV_Field_1/
  Sample_DAPI_GFP_LHX6_PV_F00.ims

Sample_DAPI_GFP_LHX6_PV_7_Field_1/
  Sample_DAPI_GFP_LHX6_PV_7_F00.ims
```

The standalone integer after the final marker is an optional acquisition
identifier. If present, the directory and every tile filename must use the same
identifier. Numbers intrinsic to marker names, including `BiVe3`, `BiVe4`,
`ASCL1`, and `NKX2_1`, are retained as marker names and are not mistaken for
acquisition identifiers.

Source tile IDs need not begin at zero. Four discovered tiles are sorted by
their original numeric IDs and staged as `F00` through `F03`. The manifest
records both identities. Incomplete fields are skipped and recorded by default;
other malformed or ambiguous inputs stop preparation.

## Configure

Copy the example and edit it only when the acquisition requires different
channel or scientific settings:

```bash
cp "IGF Organoids/configs/ims_mge_generic.example.json" \
  /scratch/parent_root/parent0/$USER/igf_organoids/ims_mge_generic.json
```

The supplied Great Lakes Fiji image remains configured at:

```text
/nfs/turbo/umms-parent/olti_ims2tiffscript/fiji_latest_bigstitcher.sif
```

Channel names default to automatic parsing from the filename in stored
wavelength order. The nuclear channel is identified from DAPI or Hoechst.
Uncommon markers can be appended to `filename_markers`, or `channel_names` can
be set explicitly in wavelength order.

The IGF alignment profile is fixed to the successful settings:

```text
--skip-phase-init --match-gate-um 25 --iterations 5 --max-points 1500
```

The example otherwise preserves the supplied IMS export, BaSiC, BigStitcher,
CellProfiler, ROI-filtering, and colocalization settings.

## Scan without writing

From the repository checkout:

```bash
python "IGF Organoids/run_folder_batch.py" scan \
  --input "$IMS_INPUT" \
  --spec "IGF Organoids/configs/ims_mge_generic.example.json"
```

Review every field, original tile order, staged order, and detected channel
mapping. `scan` does not create files.

Useful selectors include:

```text
--group 'Exp17_157-2_Set2*'
--field 1 --field 3
--max-depth 6
--require-contiguous
```

## Prepare and inspect the Slurm command

Use `--print-only` first. It creates the scratch workspace and manifests, but
does not call `sbatch`:

```bash
python "IGF Organoids/submit_folder_batch.py" \
  --input "$IMS_INPUT" \
  --spec "IGF Organoids/configs/ims_mge_generic.example.json" \
  --workspace "$IGF_WORKSPACE" \
  --max-parallel 2 \
  --print-only
```

Preparation creates:

```text
manifest.json
skipped_fields.csv
controller_compat/
input_staging/
configs/
logs/
outputs/
```

`controller_compat` is intentional and must remain. The upstream IMS entry
points call `Path.resolve()`, which would dereference a staged `F00` link back
to a differently numbered source tile before validating tile order. The
adapter changes only those path-canonicalization calls to `Path.absolute()`.
It symlinks the unmodified downstream runner and TIFF-safety module; it does not
change any scientific algorithm.

## Submit

Because preparation requires a new empty workspace, submit the already prepared
manifest after reviewing the print-only command:

```bash
python "IGF Organoids/submit_folder_batch.py" \
  --manifest "$IGF_WORKSPACE/manifest.json" \
  --max-parallel 2
```

Each array task processes exactly one biological `Field_N` acquisition. The
parallel limit controls how many full field jobs run concurrently. The helper
records the Slurm job ID in `slurm_submission.json` and refuses accidental
duplicate submissions.

Monitor with the printed job ID:

```bash
squeue -j JOB_ID
sacct -j JOB_ID --format=JobID,State,ExitCode,Elapsed,MaxRSS
tail -F "$IGF_WORKSPACE/logs/ims_JOB_ID_0.out"
tail -F "$IGF_WORKSPACE/logs/ims_JOB_ID_0.err"
```

## Resume and status

The Slurm runner always requests `--resume`. Existing stage manifests are
checked before outputs are reused, and source IMS size and modification time
must still match the preparation manifest.

```bash
python "IGF Organoids/run_folder_batch.py" inspect \
  --manifest "$IGF_WORKSPACE/manifest.json"
```

After confirming an earlier job has ended, an intentional resubmission is:

```bash
python "IGF Organoids/submit_folder_batch.py" \
  --manifest "$IGF_WORKSPACE/manifest.json" \
  --max-parallel 2 \
  --allow-duplicate
```

## CASTalign changes from the repository baseline

The supplied working Great Lakes versions of
`Alignment/CASTalign_two_channel_registration.py` and
`run_stitched_pipeline.py` are used directly.

The relevant changes are narrowly scoped:

1. The old early error requiring identical fixed and moving XY dimensions was
   removed from CASTalign.
2. The stitched runner accepts unequal native dimensions only when alignment is
   enabled. Alignment-disabled workflows still require one common XY canvas.
3. If phase correlation is enabled elsewhere, unequal images are temporarily
   zero-padded on the bottom and right to a common size only for that coarse
   initializer. Source images and centroid coordinate systems are unchanged.
4. IGF sets `--skip-phase-init`, so no phase-correlation call or padding occurs
   for this workflow; the initial translation is zero.
5. The final moving image is transformed once onto the fixed nuclear image's
   output dimensions.
6. The existing five-iteration transformation behavior is unchanged.

Phase correlation was therefore disabled for IGF, not deleted from the general
CASTalign command. This keeps other existing workflows backward-compatible.

## Intentional differences from the supplied Great Lakes wrapper

Repository integration required only these non-scientific changes:

- the three batch files and example configuration live under `IGF Organoids/`;
- their default repository root is the parent of this directory;
- the hardcoded `/home/oltij/ImageAnalysisScripts-Copy` Slurm path is replaced
  by absolute runner paths recorded in the manifest and passed as job arguments;
- the optional directory/filename acquisition identifiers are explicitly
  cross-validated;
- the successful phase-disabled CASTalign arguments are present in the example;
- documentation and focused regression tests were added.

The supplied CASTalign and stitched-pipeline Python files themselves are
otherwise reproduced exactly. No IMS export, BaSiC, BigStitcher, CellProfiler,
ROI-filtering, or colocalization algorithm was rewritten.
