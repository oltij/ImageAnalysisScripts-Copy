# Stage 09: nucleus-constrained marker matching

Stage 09 is an additional marker-marker analysis. It does not replace or
modify either existing analysis:

- Stage 07 independently matches each marker ROI to the nuclear channel.
- Stage 08 performs the original unconstrained marker-marker overlap and
  stable matching.
- Stage 09 repeats marker-marker matching independently after requiring both
  marker objects to have the same verified Stage-07 canonical nuclear
  `ObjectID`.

Both Stage 08 and Stage 09 run by default for every configured marker pair.
The all-nuclei A+/B+, A+/B-, A-/B+, and A-/B- classification remains a
separate Stage-07-derived outcome; it is not redefined as a Stage-09 count.

## Scientific definition

A candidate A-B pair is eligible in Stage 09 only when all of the following
are true:

1. It passes the existing directional threshold, A pixels inside B.
2. The current A and B ROIs each resolve unambiguously to their respective
   Stage-07 marker object by an exact pixel-coordinate fingerprint.
3. The two Stage-07 assignments point to the same original nuclear Object B
   `ObjectID`.

This restriction is applied before preference lists are constructed. The
existing A-proposing Gale-Shapley implementation, overlap ranking, tie
behavior, thresholds, images, and calibration are then reused without change
on the reduced candidate set. Stage 09 is recomputed from the nuclear-positive
marker pixel tables; it is not produced by deleting pairs from Stage 08.
Consequently a constrained final pair can be absent from the unconstrained
final matching.

`CellID` is local to a single matching run and is never treated as a canonical
identity. Stage 09 records and distinguishes:

- the current marker-marker `ObjectID`;
- the original Stage-07 marker `ObjectID`;
- the Stage-07 run-local source `CellID`;
- the canonical Stage-07 nuclear Object B `ObjectID`; and
- the new Stage-09 output `CellID`.

Current-to-Stage-07 marker identity is verified from exact `X,Y` pixel sets.
ROI and CellID values are cross-checks only. Missing, coordinate-ambiguous, or
identifier-conflicting mappings are excluded with an explicit audit reason.
Before matching, the two Stage-07 sources must also have the same nuclear
channel name, image canvas, pixel calibration, and exact nuclear
`ObjectID/X/Y` segmentation table. An incompatibility stops Stage 09 rather
than guessing.

## Configuration and execution

The stitched controller accepts this optional top-level setting:

```json
"nucleus_constrained_marker_matching": true
```

It defaults to `true`. Set it to `false` only when intentionally requesting
the legacy Stage-07/08 workflow without Stage 09. Pair directions and
thresholds still come from the existing `pairs` and
`thresholds.marker_marker` settings.

On Great Lakes, use the normal committed controller and the existing Conda and
Slurm wrappers. No special environment is required:

```bash
conda activate ims-mosaic
REPO="$HOME/ImageAnalysisScripts-Copy"
python "$REPO/run_stitched_pipeline.py" \
  --config /scratch/parent_root/parent0/$USER/my_run/stitched_pipeline.json \
  --dry-run

# Run inside the existing Slurm workflow, or resume the same tracked output:
python "$REPO/run_stitched_pipeline.py" \
  --config /scratch/parent_root/parent0/$USER/my_run/stitched_pipeline.json \
  --resume
```

The general and IGF wrappers pass this setting through their normal analysis
configuration. Continue to use their documented Slurm submission commands;
do not run a microscopy batch on a login node.

Resume labels are independent. A verified completed Stage 08 is skipped while
a missing Stage 09, its visualization stage, or its comparison stage runs.
As elsewhere in the repository, resume requires an unchanged manifest
fingerprint; use a new output directory after changing scientific code or
configuration.

## Output tree

The original output remains at:

```text
analysis/08_marker_colocalization/<A>_vs_<B>/
```

The additional output is separate:

```text
analysis/09_marker_colocalization_nucleus_constrained/<A>_vs_<B>/
  analysis.h5
  nuclear_constraint_pair_audit.csv
  visualizations/
  comparison_to_unconstrained/
    matching_variant_comparison.csv
    matching_variant_pair_details.csv
    discordant_unconstrained_final_pairs.csv
    matching_variant_comparison.pdf
    matching_variant_comparison.png
```

`analysis.h5` retains the standard images, metadata, pixel tables, overlap
results, eligible pairs, preference matrix, stable matching, CellID mapping,
filtered pixels, summaries, and filtering counts. It additionally stores the
two verified nuclear-assignment tables and `nuclear_constraint_audit`.
Metadata records both Stage-07 source HDF5 paths and the identity-validation
policy.

The audit has one row per spatially overlapping candidate and records original
A/B object IDs, Stage-07 source CellIDs, canonical nuclear ObjectIDs, mapping
status, coordinate verification, both directional overlap percentages,
original threshold status, same-nucleus/constraint status, final eligibility,
final constrained-match status, and an exclusion reason. Zero eligible or
zero final matches are valid: empty tables and an informative matched-pair
figure are emitted instead of crashing.

The comparison outputs report input ROI counts, overlapping candidates,
threshold-passing candidates, different-nucleus and unresolved exclusions,
constraint-passing candidates, final counts in both variants, common and
variant-only final pairs, and the fraction of unconstrained final pairs whose
resolved nuclear identities disagree. Pair-detail and discordance CSVs retain
the exact provenance identifiers.

Visualization filenames retain their established names inside the Stage-09
parent directory. Titles explicitly say `NUCLEUS-CONSTRAINED`. The nuclear
classification tables copied into visualization output continue to describe
the same Stage-07 all-nuclei classification and are not Stage-09 match counts.

## Per-organoid reports

The existing reporter remains unconstrained by default. Request the separate
Stage-09 report explicitly:

```bash
python FigureAnalysis/mNeonGreenReporterAnalysis/build_per_organoid_reports.py \
  --h1-workspace "$W1" \
  --remaining-workspace "$W2" \
  --marker-pair-variant nucleus-constrained \
  --max-examples 10
```

It reads marker-pair HDF5 files only from the Stage-09 tree, labels the pair as
`NUCLEUS-CONSTRAINED`, and adds `_nucleus_constrained` to report and CSV names.
The unconstrained files and output directory are not overwritten. The
established maximum of 10 displayed exemplars per quartile and ROI-level
traceability are unchanged.

## Focused validation

The synthetic tests exercise same-nucleus eligibility, different-nucleus
exclusion despite strong overlap, below-threshold exclusion, missing and
ambiguous provenance, incompatible nuclear segmentations, source CellID versus
original ObjectID, zero matches, independent rematching, output creation,
Stage-08 immutability, reporter discovery, and resume labels:

```bash
python -m py_compile \
  Colocalization/colocalizationdapiscript1.py \
  Colocalization/colocalizationdapiscript2.py \
  Colocalization/compare_matching_variants.py \
  run_stitched_pipeline.py

python -m pytest -q \
  tests/test_nucleus_constrained_matching.py \
  tests/test_stitched_pipeline.py \
  FigureAnalysis/mNeonGreenReporterAnalysis/tests/test_per_organoid_reports_quartiles.py
```

These tests use synthetic arrays only and do not submit Slurm jobs or create
experimental results.
