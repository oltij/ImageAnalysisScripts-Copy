# IGF organoid IMS folder batch

This directory contains the production Great Lakes orchestration around the
repository's existing microscopy pipeline. It starts from each field's four
original `.ims` tiles and runs the established workflow:

1. IMS export
2. BaSiC correction
3. BigStitcher 2×2 mosaic assembly and maximum projections
4. initial CellProfiler segmentation for registration anchors
5. CASTalign registration
6. segmentation again in the aligned coordinate system
7. shared DAPI-tissue-mask QC, ROI filtering, matching, and colocalization

The wrappers do not replace, optimize, or change those scientific algorithms.

## Committed files

- `run_folder_batch.py`: session-scoped discovery, metadata validation, scratch
  staging, per-field manifests, execution, inspection, and resume support.
- `submit_folder_batch.py`: generic preparation and single-array submission.
- `submit_production_batches.py`: the production 94-field submission as two
  dependent arrays of 47.
- `run_folder_batch.sbatch`: the validated Great Lakes per-field job profile.
- `configs/ims_mge_generic.example.json`: the unchanged production scientific
  configuration.
- `configs/production_94_fields.json`: the audited sessions, expected counts,
  channel mapping, and configuration checksum for the 94-field run.

The scripts resolve the repository root from their own absolute locations.
Absolute manifest and runner paths are passed to Slurm because a job executes a
spooled copy of the `.sbatch` file.

## Exact production procedure

No command in the inspection, preparation, validation, preflight, or preview
steps submits a Slurm job. Only step 7 calls `sbatch`.

### 1. Pull a clean GitHub checkout

For a new checkout:

```bash
cd /home/$USER
git clone https://github.com/oltij/ImageAnalysisScripts-Copy.git
cd ImageAnalysisScripts-Copy
git switch main
git pull --ff-only origin main
git status --short
```

`git status --short` must print nothing. Preparation with the production run
plan records `git rev-parse HEAD` and refuses a dirty checkout. The escape hatch
`--allow-dirty-repository` exists for development only and should not be used
for this production run.

Activate the existing orchestration environment:

```bash
conda activate ims-mosaic
```

Set shell variables from the repository root:

```bash
REPO_ROOT="$PWD"
IMS_INPUT_ROOT=/nfs/turbo/umms-parent/andor_microscope_data_dump
IGF_CONFIG="$REPO_ROOT/IGF Organoids/configs/ims_mge_generic.example.json"
IGF_RUN_PLAN="$REPO_ROOT/IGF Organoids/configs/production_94_fields.json"
IGF_WORKSPACE=/scratch/parent_root/parent0/$USER/igf_organoids/production_94_001
```

Use a new workspace name if that path already exists.

### 2. Inspect exactly the 12 production sessions, read-only

```bash
python "$REPO_ROOT/IGF Organoids/run_folder_batch.py" scan \
  --input "$IMS_INPUT_ROOT" \
  --spec "$IGF_CONFIG" \
  --run-plan "$IGF_RUN_PLAN"
```

The committed plan expands to these sessions only:

| Session | Complete fields |
|---|---:|
| `Olti_06_23_2026` | 12 |
| `Olti_06_25_2026` | 8 |
| `Olti_06_29_2026` | 7 |
| `Olti_06_30_2026` | 6 |
| `Olti_06_30_2026-secondslide` | 6 |
| `Olti_07_01_2026` | 6 |
| `Olti_07_02_2026` | 15 |
| `Olti_07_03_2026` | 6 |
| `Olti_07_06_2026` | 6 |
| `Olti_07_08_2026` | 7 |
| `Olti_07_17_2026` | 7 |
| `Olti_07_24_2026` | 8 |
| **Total** | **94** |

`Olti_07_16_2026` and every unlisted session are excluded. The scan must end
with `RUN PLAN VALIDATED (read-only)` and confirm:

- 94 complete fields;
- 19 incomplete fields skipped and reported;
- 376 original IMS tiles;
- four physical channels per tile, indices 0–3;
- biological channels `DAPI`, `GFP`, `LHX6`, `PV` in that order;
- no phase channel and no cross-tile metadata inconsistency.

`--session NAME` is also repeatable for generic scans and preparations. Session
selectors must be unique direct child directory names; missing names, symlinks,
duplicates, absolute paths, and traversal components are rejected.

### 3. Confirm a fresh scratch workspace

```bash
if [[ -e "$IGF_WORKSPACE" ]]; then
  echo "Choose a new IGF_WORKSPACE; this path already exists" >&2
  return 1 2>/dev/null || exit 1
fi
mkdir -p "$(dirname "$IGF_WORKSPACE")"
df -h "$(dirname "$IGF_WORKSPACE")"
```

On Great Lakes, `/scratch/parent_root/parent0/...` may canonically resolve to
`/gpfs/accounts/parent_root/parent0/...`. Validation resolves both the workspace
and `/scratch` before comparison, so this valid layout is accepted. Inputs
under the protected Turbo root are still rejected if the canonical workspace
is outside scratch.

### 4. Prepare and validate all 94 fields

```bash
python "$REPO_ROOT/IGF Organoids/run_folder_batch.py" prepare \
  --input "$IMS_INPUT_ROOT" \
  --spec "$IGF_CONFIG" \
  --run-plan "$IGF_RUN_PLAN" \
  --workspace "$IGF_WORKSPACE"

python "$REPO_ROOT/IGF Organoids/run_folder_batch.py" validate \
  --manifest "$IGF_WORKSPACE/manifest.json"
```

Preparation stops before submission unless every audited count and channel
mapping matches. The manifest records:

- the Git commit and clean/dirty status;
- the scientific configuration path and SHA-256;
- the run-plan path, SHA-256, and complete validation report;
- the batch and Slurm runner hashes;
- selected sessions, source paths and source file statistics;
- physical channel indices/names, biological names, and any `SetN` identifier;
- deterministic job indices 0–93 and original-to-staged tile mappings.

The original Turbo tree is only read and stat-checked. The wrapper writes all
symlinks, compatibility controllers, generated configurations, manifests,
logs, preprocessing, and analysis outputs below `IGF_WORKSPACE`.

The generated scratch `controller_compat` files are derived from the committed
checkout. They change only path canonicalization from `Path.resolve()` to
`Path.absolute()` in the two IMS entry points so a staged `F00` link is not
renamed back to its original tile number. The downstream runner and TIFF-safety
module remain symlinks to the committed repository. No scientific algorithm is
changed.

### 5. Run representative metadata preflights

Run both array boundaries and one representative field in each batch:

```bash
for IGF_INDEX in 0 46 47 93; do
  python "$REPO_ROOT/IGF Organoids/run_folder_batch.py" preflight \
    --manifest "$IGF_WORKSPACE/manifest.json" \
    --index "$IGF_INDEX"
done
```

Preflight validates the original source links, source file statistics, physical
IMS metadata, calibration, channel ordering, and intended commands without
performing stitching or analysis.

### 6. Preview the two Slurm arrays

```bash
python "$REPO_ROOT/IGF Organoids/submit_production_batches.py" \
  --manifest "$IGF_WORKSPACE/manifest.json" \
  --print-only
```

The preview must show:

```text
Batch 1: --array=0-46%47
Batch 2: --array=47-93%47 --dependency=afterok:<BATCH1_JOB_ID>
```

Preview mode does not call `sbatch` and does not write a submission record.

### 7. Submit both arrays

```bash
python "$REPO_ROOT/IGF Organoids/submit_production_batches.py" \
  --manifest "$IGF_WORKSPACE/manifest.json"
```

Batch 1 permits at most 47 simultaneous field tasks. Batch 2 also permits 47,
but Slurm holds it under `afterok` until every Batch 1 array task succeeds. The
helper atomically records both IDs and exact commands in:

```text
$IGF_WORKSPACE/slurm_two_batch_submission.json
```

An existing record blocks accidental duplicate submission. If Batch 1 was
accepted but the Batch 2 `sbatch` call failed, repair the submission issue and
run exactly:

```bash
python "$REPO_ROOT/IGF Organoids/submit_production_batches.py" \
  --manifest "$IGF_WORKSPACE/manifest.json" \
  --resume-second-batch
```

This reuses the recorded Batch 1 job ID instead of duplicating Batch 1.

### 8. Monitor and diagnose

Read the IDs:

```bash
python -c 'import json,sys; d=json.load(open(sys.argv[1])); print(d["batch1_job_id"], d["batch2_job_id"])' \
  "$IGF_WORKSPACE/slurm_two_batch_submission.json"
```

Then substitute them below:

```bash
squeue -j BATCH1_JOB_ID,BATCH2_JOB_ID
sacct -j BATCH1_JOB_ID,BATCH2_JOB_ID \
  --format=JobID,JobName,State,ExitCode,Elapsed,MaxRSS,NodeList,Reason

tail -F "$IGF_WORKSPACE/logs/ims_BATCH1_JOB_ID_0.out"
tail -F "$IGF_WORKSPACE/logs/ims_BATCH1_JOB_ID_0.err"

python "$REPO_ROOT/IGF Organoids/run_folder_batch.py" inspect \
  --manifest "$IGF_WORKSPACE/manifest.json"
```

Each log filename is `ims_<array-job-id>_<array-index>.out` or `.err`. Use the
failed array index reported by `sacct` to inspect the corresponding log. Batch 2
will not run if any Batch 1 task fails. Per-field execution always uses
`--resume`; unchanged completed stages can therefore be reused after a failure.
Do not use a duplicate-submission override until the recorded jobs have ended
and their state is understood.

## Source and scratch safety

The protected source root is:

```text
/nfs/turbo/umms-parent/andor_microscope_data_dump
```

Original `.ims` files are never renamed, copied over, modified, or deleted.
Staged inputs are scratch symbolic links. Preparation rejects a workspace
inside the source tree, rejects a nonempty workspace, and—for protected Turbo
inputs—rejects any workspace outside canonical scratch.

The per-field Slurm profile remains:

- account `parent0`;
- partition `standard`;
- one node and one task;
- four CPUs;
- 64 GB memory;
- one-day time limit;
- no GPU;
- `singularity/4.4.1`;
- `/gpfs:/gpfs,/scratch:/scratch` container binds;
- `ims-mosaic` for orchestration/export/stitching;
- `cellprofiler-native` through the downstream configuration.

## Naming and biological-channel safety

Discovery supports an optional date directory and field layouts such as:

```text
Sample_Set2_DAPI_GFP_LHX6_PV_Field_1/
  Sample_Set2_DAPI_GFP_LHX6_PV_F00.ims

Sample_Set2_DAPI_GFP_LHX6_PV_7_Field_1/
  Sample_Set2_DAPI_GFP_LHX6_PV_7_F00.ims
```

`Set1`, `Set2`, `Set3`, and similar tokens are experiment metadata, never
biological channels. The set number is stored separately in each manifest job.
An optional standalone acquisition integer after the final marker must agree
between the field directory and every tile. Digits intrinsic to `BiVe3`,
`BiVe4`, `ASCL1`, and multi-token markers such as `NKX2_1` remain part of their
marker names.

Physical IMS channel count, indices, and stored names are checked across all
four tiles. A phase channel or a mismatch between physical channel count and
filename biological markers fails before processing. In particular, a
five-channel phase-plus-DAPI/GFP/LHX6/PV acquisition is rejected; the wrapper
never fabricates `Set3` as a fifth marker. Partial phase-channel export is not
implemented here.

Source tile IDs need not start at zero. Four complete tiles are sorted by their
numeric IDs and staged as `F00`–`F03`; both identities remain in the manifest.
Incomplete fields are skipped and listed in `skipped_fields.csv`. Extra,
duplicate, ambiguous, or inconsistent inputs remain fatal.

## Preserved scientific configuration

`configs/ims_mge_generic.example.json` remains the production configuration.
It preserves native IMS resolution and voxel calibration,
`analysis_downsampling.factor=1`, 2×2 stitching, existing BaSiC and BigStitcher
parameters, CellProfiler segmentation and shared DAPI mask, ROI filters,
matching, colocalization, and report generation.

IGF registration remains exactly:

```text
--skip-phase-init --match-gate-um 25 --iterations 5 --max-points 1500
```

CASTalign's existing five-iteration behavior is unchanged. The previous early
equal-dimension rejection remains removed: when alignment is enabled, markers
with different native stitched XY canvases are registered onto the fixed DAPI
canvas. Phase correlation is disabled for IGF, not deleted from the general
CASTalign tool.
