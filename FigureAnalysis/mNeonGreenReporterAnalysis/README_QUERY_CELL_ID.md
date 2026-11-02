# Query a Stage-07 marker CellID, canonical nucleus, or Stage-08/09 PairID

This is a **read-only per-field investigation tool**. It reads an existing field's `analysis/` folder and outputs labeled fluorescence images, a PDF, per-marker native intensity/size measurements, and **all Stage-08/09 spatial candidate pairs involving the queried nucleus**. Both final matching outcomes are explicitly identified.

## Files to integrate

- `FigureAnalysis/mNeonGreenReporterAnalysis/query_cell_id.py` — standalone executable script
- `FigureAnalysis/mNeonGreenReporterAnalysis/tests/test_query_cell_id.py` — synthetic regression tests

Integrate at these exact paths in `oltij/ImageAnalysisScripts-Copy`. The script uses the real Stage-07 source HDF5 `tables/cell_id_mapping`, Stage-08/09 `tables/stable_matching` and `tables/overlap_results`, stored pixel tables and native fluorescence images. It also checks the Stage-09 nuclear constraint audit.

## IDs: keep the scopes explicit

**Stage-07 CellID is not globally unique.** For example, LHX6 CellID 12 and PV CellID 12 are independently assigned numbers and may belong to different nuclei. Always identify a Stage-07 CellID with field + marker + number. The canonical shared identity is the **Hoechst nuclear ObjectID**. Neither the Stage-08 nor Stage-09 final pairwise CellID is a Stage-07 marker CellID.

Recommended labels:

- `Field <field>, LHX6 Stage-07 CellID 12 (ROI 247), Hoechst nucleus 501`
- `Field <field>, PV Stage-07 CellID 29 (ROI 83), Hoechst nucleus 501`
- `Stage 08: LHX6 C12 ↔ PV C29; same nucleus: YES (501 ↔ 501)`

The script uses long `Field_Analysis_Path::Marker::C<CellID>` keys in exported CSV to prevent accidental cross-field or cross-marker joins. These keys identify objects within a **specific processing run**; after segmentation or matching is recomputed, CellID numbering can change. Use original ROI and nuclear IDs plus pixel provenance to compare runs.

## Great Lakes usage

```bash
conda activate ims-mosaic
REPO="$HOME/ImageAnalysisScripts-Copy-reports"
git -C "$REPO" pull --ff-only

# Replace with an actual completed field directory containing 07_... / 08_... / 09_...
ANALYSIS="/scratch/parent_root/parent0/$USER/igf_organoids/<workspace>/outputs/<field>/analysis"

# Query one marker-associated CellID
python "$REPO/FigureAnalysis/mNeonGreenReporterAnalysis/query_cell_id.py" \
  --analysis "$ANALYSIS" --marker LHX6 --cell-id 12

# Or query by a canonical Hoechst nuclear ID directly
python "$REPO/FigureAnalysis/mNeonGreenReporterAnalysis/query_cell_id.py" \
  --analysis "$ANALYSIS" --nuclear-id 501

# Or query a deterministic final colocalization pair
python "$REPO/FigureAnalysis/mNeonGreenReporterAnalysis/query_cell_id.py" \
  --analysis "$ANALYSIS" --pair-id S08_LHX6_PV_P000017
```

If `--cell-id` is supplied without `--marker`, the script accepts it only when that numeric CellID occurs in **exactly one** Stage-07 marker source; otherwise it rejects the ambiguous input and lists the alternatives. `--overwrite` regenerates the output for the selected nucleus. The default output is `<analysis>/cell_id_queries/nucleus_<id>/`. Set `--max-pair-visuals 8` (default) to limit the number of *visual* examples **per stage**; the overlap-candidates CSV remains complete.

Exactly one of `--cell-id`, `--nuclear-id`, or `--pair-id` is required. A
PairID query detects its variant, recovers both marker identities and native
measurements, and checks the exact `(A_ObjectID, B_ObjectID)` pair in the other
variant. PairID numbers are field- and variant-scoped; number equality is never
used as cross-variant identity.

### Outputs

- `marker_measurements.csv`: each available marker's Stage-07 CellID, original ROI ID, nucleus ID, native mean and integrated marker/Hoeschst intensity, ROI pixel area, Stage-07 marker/nucleus shared pixel count, and any final-ROI metrics (Area_um2, circularity, solidity, etc.). Unassigned markers remain visible as negative rows.
- `stage08_stage09_overlap_candidates.csv`: every candidate spatially overlapping the selected nucleus's Stage-07-associated marker ROIs in available Stages 08 and 09, with BOTH marker-specific Stage-07 CellIDs, original and current ROI IDs, both canonical nuclear IDs, YES/NO/unknown shared-nucleus status, directional overlap, areas, raw native fluorescence, whether the pair is a **final** stable match, the pair-run local CellID, and Stage-09 rejection reason if applicable.
- `stage08_stage09_variant_summary.csv`: counts of candidate and final matches involving that nucleus for each processed marker pair and stage. Absent stages are explicitly marked.
- `marker_fluorescence.png`: labeled marker/Hoeschst/multi-channel crops with segmentation contours.
- `cell_id_report.pdf`: marker/nucleus crop overview followed by representative Stage-08 and Stage-09 candidate pair crops (final matches prioritized) with measured overlap and both nuclear identities.
- `cell_identity_summary.json`: nucleus ID, selection provenance, available analyses and per-variant matching counts.

PairID queries instead write under `<analysis>/pair_id_queries/<PairID>/`:

- `pair_query_summary.json`;
- `pair_query_measurements.csv` with full-ROI and overlap fluorescence;
- `pair_query_stage08_stage09_cross_reference.csv`;
- `pair_query_report.pdf` with both fluorescence crops, composite, IDs,
  same-nucleus status, and cross-variant presence.

### Testing

```bash
python -m py_compile FigureAnalysis/mNeonGreenReporterAnalysis/query_cell_id.py
python -m pytest -q FigureAnalysis/mNeonGreenReporterAnalysis/tests/test_query_cell_id.py
```

The synthetic tests cover ambiguous marker CellID rejection, two marker runs using the same numeric CellID for different nuclei, complete per-marker ID tracing, discordant Stage-08 versus same-nucleus Stage-09 final pairing, independent matching statuses, missing Stage-09 results, Stage-08 and Stage-09 PairID queries, common-pair cross-referencing, mixed-selector rejection, and report generation. **The new script has not been tested on actual Great Lakes microscopy fields.**

## Important distinction for the report-generator integration

`query_cell_id.py` is a **separate standalone query tool**. The requested changes to the existing 10-exemplar quartile report and all-field counts are described in `CODEX_INSTRUCTIONS.md` and should be implemented by Codex; do not replace the existing report with this single-nucleus tool.
