# Source-script interface and handoff reference

This catalog documents the role of **every existing script** in the repository. `run_ims_pipeline.py` covers the complete IMS entry point and hands projections/calibration to `run_stitched_pipeline.py`. The automation does not alter source algorithms, percentile math, matching strategy, or segmentation settings. It routes outputs and generates an ephemeral source copy containing only top-level `USER SETTINGS` substitutions for scripts that lack a CLI.

| Script | Position | Inputs | Output and handoff | Runner behavior |
|---|---|---|---|---|
| `Stitching/IMSTIFF_trim_zero_padding.py` | IMS mosaic preprocessing | Original `_FNN.ims` fields | Exported, BaSiC-corrected, fused and max-projected channel TIFFs | Run by `run_ims_pipeline.py` in mosaic mode with fixed channel names supplied by index; not run by the stitched-TIFF entry point |
| `Stitching/maxprojectscript.py` | Optional upstream example | Fixed-path TIFF Z-stack | Fixed-path max projection | **Not run**: the controller validates axes and performs the same axis-0 max projection directly. The example's hard-coded paths and dangling `PY` token cannot affect this workflow. |
| `BigStitcher/fiji_latest_bigstitcher.def` | Optional upstream infrastructure | Singularity definition | BigStitcher environment | Required by `run_ims_pipeline.py` for mosaic fields; not used for an already-stitched IMS/TIFF |
| `CellProfiler/cellprofilerdriver.py` | Initial and final-frame segmentation | `--pipeline`, `--input`, `--output` | `OrganoidMask.tiff`, `CellMask.tiff`, `MyExpt_FilterObjects2.csv`, `MyExpt_FilterObjects.csv` | Runs via `conda run -n cellprofiler-native`; the initial pass is omitted when alignment is disabled |
| `Alignment/CASTalign_two_channel_registration.py` | Registration | fixed/moving images + CellProfiler object CSVs, centroid CLI flags | `aligned_PV_max.tif`, `aligned_Hoechst_max.tif`, transform and QC | `--output` set per marker, link generic filename to true moving marker name; no modification to registration code |
| `QC/organoid_mask_qc.py` | Final-frame segmentation QC | Aligned fluorescence, `OrganoidMask.tiff`, and `CellMask.tiff` for every channel | Actual-mask panels, containment overlay, and `organoid_mask_statistics.csv` | Runs after aligned segmentation and before ROI extraction; reports only and does not filter objects |
| `ExtractingROIs/extract_rois.py` | Aligned ROI reconstruction | CP label TIFF, marker object CSV, organoid object CSV | `<Marker>_ROI_pixels.csv`, verification/summary and reconstructed labels | Temporary settings replacement; source unchanged |
| `ROIFiltering/Intensity/intensityfilter.py` | Intensity QC | aligned ROI-pixel CSV(s) + aligned fluorescence image(s) | `All_Channel_ROI_Intensity_Report.pdf`, CSV reports | All channels in one run; no filtering done here |
| `ROIFiltering/FinalFilter/filter.py` | Two distinct passes | ROI-pixel CSV(s), aligned image(s), `FILTERS` | `<Marker>_ROI_pixels_FILTERED.csv`, cutoffs, audit tables | Pass 1 disables all shape metrics; pass 2 consumes pass-1 CSVs and disables all intensity metrics |
| `ROIFiltering/Area/areafilter.py` | Shape QC | intensity-positive pixel CSV(s) + aligned images | shape PDFs and CSVs | Runs on intensity-gated ROIs only |
| `ROIFiltering/Circularity/circularityfilter.py` | Shape QC | same as above | shape PDFs and CSVs | Runs independently on same intensity-gated ROIs |
| `ROIFiltering/Eccentricity/eccentricityfilter.py` | Shape QC | same as above | shape PDFs and CSVs | Runs independently on same intensity-gated ROIs |
| `ROIFiltering/Solidity/solidityfilter.py` | Shape QC | same as above | shape PDFs and CSVs | Runs independently on same intensity-gated ROIs |
| `Colocalization/colocalizationdapiscript1.py` | Nuclear association and pair analysis | A/B pixel CSV + A/B aligned TIFF + optional organoid CSV + overlap threshold | `analysis.h5` containing `tables/filtered_object_a_pixels`, mapping, overlap and summary tables | Run once per marker vs nucleus, then once per requested marker pair |
| `Colocalization/colocalizationdapiscript2.py` | QC, plots, all-nuclear summary | `analysis.h5`, optional two marker-vs-nucleus HDF5 files and organoid-area CSV | matched/unmatched PDFs, pair label TIFFs, population and all-cell tables | Runs after each HDF5; nuclear-summary flag off for marker-vs-nucleus, on for marker-vs-marker |

## Source documentation retained, but see the companion orchestration documentation

- `README.md`: original eight-step workflow, including upstream IMS/BigStitcher route.
- `IMS_WORKFLOW.md`: one-command IMS channel assignment, preprocessing, metadata, and downstream handoff.
- `Stitching/README_IMS_Mosaic_Pipeline_updated(1).md`, `Stitching/README_TIFF_max_projection.md`: optional upstream TIFF processing.
- `BigStitcher/README_BigStitcher_DEF_exact.md`: optional container construction.
- `CellProfiler/README_CellProfiler_headless_runner (1).md`: source-specific CP semantics and outputs.
- `Alignment/README_CASTalign_2D_maxproj_registration.md`: CASTalign CLI and QC.
- `ExtractingROIs/README_ROI_reconstruction_export.md`: ROI reconstruction checks.
- `ROIFiltering/Intensity/README_ROI_intensity_report.md`: intensity distribution and metrics.
- `ROIFiltering/Area/README_ROI_area_report.md`, `ROIFiltering/Circularity/README_ROI_circularity_report.md`, `ROIFiltering/Eccentricity/README_ROI_eccentricity_report.md`, `ROIFiltering/Solidity/README_ROI_solidity_report.md`: shape diagnostics.
- `ROIFiltering/FinalFilter/README_ROI_final_metric_filtering.md`: two-pass filter settings.
- `Colocalization/README_01_run_colocalization.md`, `Colocalization/README_02_generate_colocalization_visualizations_updated.md`: stable matching, HDF5 contracts, and visualizations.

## Intentional boundaries

The controller validates that a `.cppipe` is text and contains the output names needed by the handoffs, but it does not rewrite pipeline modules, reinterpret microscope channel naming, auto-select thresholds, or infer pixel calibration. The source scripts were built for particular microscopy data rather than arbitrary image formats; manual inspection of first-sample masks, alignment, and counts remains mandatory.
