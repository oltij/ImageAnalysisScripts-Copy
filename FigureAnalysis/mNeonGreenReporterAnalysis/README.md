# Experiment 17 ventral IGF/mNeonGreen reporter analysis

This directory contains the read-only plotting analysis for the Experiment 17
ventral organoids. It compares untreated (`No IGF`) and IGF1-treated (`+IGF`)
fields for H1, CL23, and CL32 and produces:

- 14 PV/LHX6 figures;
- 6 mNeonGreen/GFP figures split by anti-GFP staining status;
- one combined multipage PDF;
- field-, object-, sample-, metadata-, and selection-audit CSV files.

The script reads completed outputs from the existing IGF microscopy pipeline.
It does not run segmentation, registration, filtering, or colocalization and
does not modify either Great Lakes workspace.

## Files in this directory

- `plot_igf_exp17_ventral_with_gfp_plots.py`: the supplied plotting and data
  validation script. Scientific calculations, filtering rules, experiment
  assignments, plotting layouts, and figure legends are preserved.
- `igf_exp17_ventral_block_mapping.csv`: the supplied authoritative mapping of
  the six Experiment 17 ventral blocks.
- `.gitignore`: excludes generated figures/tables and large microscopy data
  from this source directory.

The mapping CSV is resolved relative to the script by default. Scratch
workspaces, the verified H1 selection, and the output directory remain
configurable command-line paths.

## Experimental design

Only Experiment 17 ventral organoids are eligible. Experiment 16, dorsal
organoids, unknown blocks, and known failed manifest indices 12 and 13 are
excluded. Both Great Lakes manifests must describe the same 94 fields and must
agree on manifest indices, source fields, original tile IDs, and source paths.

The block-to-condition assignment is:

| Block | Cell line | Condition | Basis recorded in mapping | Inventory note |
|---:|---|---|---|---|
| 149 | CL32 | No IGF | Inferred untreated because IGF1 is absent from the inventory | PV32 Ruku |
| 150 | CL32 | +IGF | Explicit IGF1 | PV32 Ruku IGF1; stored after drying one day late |
| 151 | H1 | +IGF | Explicit IGF1 | H1 IGF1 Ruku |
| 152 | CL23 | +IGF | Explicit IGF1 | PV23 Ruku IGF1 |
| 153 | CL23 | No IGF | Inferred untreated because IGF1 is absent from the inventory | PV23 Ruku |
| 157 | H1 | No IGF | Inferred untreated because IGF1 is absent from the inventory | H1 Ruku |

The inferred `No IGF` assignments should be confirmed against wet-lab records.
The block 150 storage note is retained as a potential QC covariate; it does not
change any calculation.

H1 uses the previously verified field indices listed in
`PV_density_per_field_verified_H1.csv`. CL23 and CL32 use all complete eligible
Experiment 17 ventral fields. Workspace selection is preserved from the
supplied script:

- indices 0–46 are read from the original/H1 workspace;
- index 47 may be read from either workspace, preferring the original one;
- indices 48–93 are read from the remaining-fields workspace.

A field is complete only when its `analysis/run_manifest.json` contains a
completed key matching `08_pair_*_viz`.

### Anti-GFP staining

The acquisition section suffix determines anti-GFP status:

| Sample suffix | Interpretation |
|---|---|
| `-2` | Anti-GFP stained |
| `-3` | Unstained/no anti-GFP |

Any included sample without one of these suffixes is rejected. The PV/LHX6
figures combine the eligible fields without making anti-GFP a panel dimension.
The six mNeonGreen/GFP figures use separate rows for stained and unstained
fields.

## Plot conventions and legends

Every dot represents one microscopy field. Dot color identifies the
sample/block-section, not an independent experiment. Within each IGF condition,
the black horizontal line is the arithmetic mean of the displayed field-level
values. The x-axis order is `No IGF`, then `+IGF`.

PV/LHX6 figures have one panel each for H1, CL23, and CL32. GFP figures have
those cell-line columns and two anti-GFP rows. Each figure uses a common y-axis
range across its visible panels, calculated from all finite values for that
metric.

Nuclear classes come from the canonical Hoechst-associated LHX6/PV
classification:

- `LHX6+/PV+` — double-positive nucleus;
- `LHX6+/PV-` — LHX6-only nucleus;
- `LHX6-/PV+` — PV-only nucleus;
- `LHX6-/PV-` — double-negative nucleus.

“PV ROI associated with a nucleus” does not mean that PV signal is physically
restricted to the nucleus. It means that the final filtered PV object was
linked to that nucleus by the upstream analysis. GFP measurements use all
final kept GFP ROIs and are not restricted to nuclei.

The script performs descriptive plotting only. It reports no p-values.
Multiple fields from one sample/block-section should not automatically be
treated as independent biological replicates.

## The 14 PV/LHX6 figures

Each metric is exported as `<stem>_side_by_side.pdf` and
`<stem>_side_by_side.png`.

| # | Output stem | Biological meaning | Calculation and normalization |
|---:|---|---|---|
| 1 | `pv_nuclear_area_percent` | Nuclear-associated PV-positive ROI burden within the organoid | Sum of areas of PV ROIs associated with PV-positive nuclei ÷ Hoechst organoid area × 100 |
| 2 | `pv_nuclear_mean_area` | Average PV-associated segmented area per PV-positive nucleus | Arithmetic mean PV ROI area across `LHX6+/PV+` and `LHX6-/PV+` nuclei; µm²/PV-positive nucleus |
| 3 | `pv_nuclear_median_area` | Typical PV-associated segmented area per PV-positive nucleus, less sensitive to large ROIs | Median PV ROI area across `LHX6+/PV+` and `LHX6-/PV+` nuclei; µm²/PV-positive nucleus |
| 4 | `pv_area_double_positive` | PV ROI size specifically in LHX6/PV double-positive cells | Mean PV ROI area among `LHX6+/PV+` nuclei; µm²/double-positive nucleus |
| 5 | `pv_area_pv_only` | PV ROI size in PV-positive, LHX6-negative cells | Mean PV ROI area among `LHX6-/PV+` nuclei; µm²/PV-only nucleus |
| 6 | `double_positive_density` | Abundance of LHX6/PV double-positive cells | `LHX6+/PV+` nuclei ÷ Hoechst organoid area in mm² |
| 7 | `lhx6_only_density` | Abundance of LHX6-only cells | `LHX6+/PV-` nuclei ÷ Hoechst organoid area in mm² |
| 8 | `pv_only_density` | Abundance of PV-only cells | `LHX6-/PV+` nuclei ÷ Hoechst organoid area in mm² |
| 9 | `lhx6_total_density` | Abundance of all LHX6-positive cells regardless of PV | (`LHX6+/PV+` + `LHX6+/PV-`) nuclei ÷ Hoechst organoid area in mm² |
| 10 | `pv_within_lhx6` | Fraction of the LHX6-positive population that is also PV-positive | `LHX6+/PV+` ÷ all LHX6-positive nuclei × 100 |
| 11 | `pv_total_density` | Abundance of all PV-positive cells regardless of LHX6 | (`LHX6+/PV+` + `LHX6-/PV+`) nuclei ÷ Hoechst organoid area in mm² |
| 12 | `double_positive_fraction` | Double-positive nuclei as a share of the complete nuclear population | `LHX6+/PV+` ÷ all four Hoechst nuclear classes × 100 |
| 13 | `lhx6_only_fraction` | LHX6-only nuclei as a share of the complete nuclear population | `LHX6+/PV-` ÷ all four Hoechst nuclear classes × 100 |
| 14 | `pv_only_fraction` | PV-only nuclei as a share of the complete nuclear population | `LHX6-/PV+` ÷ all four Hoechst nuclear classes × 100 |

## The 6 mNeonGreen/GFP figures

Each metric is exported as `<stem>_anti_gfp_split.pdf` and
`<stem>_anti_gfp_split.png`. The panels compare IGF conditions separately for
anti-GFP-stained (`-2`) and unstained (`-3`) fields.

| # | Output stem | Biological meaning | Calculation and normalization |
|---:|---|---|---|
| 1 | `gfp_density` | Abundance of final retained mNeonGreen/GFP objects | Number of final kept GFP ROIs ÷ Hoechst organoid area in mm² |
| 2 | `gfp_area_percent` | Total retained reporter-positive area within the organoid | Sum of final kept GFP ROI areas ÷ Hoechst organoid area × 100 |
| 3 | `gfp_area_mean` | Average segmented reporter-object size | Arithmetic mean `Area_um2` across final kept GFP ROIs; µm²/ROI |
| 4 | `gfp_area_median` | Typical reporter-object size, less sensitive to unusually large objects | Median `Area_um2` across final kept GFP ROIs; µm²/ROI |
| 5 | `gfp_intensity_mean` | Average reporter intensity across retained objects | Arithmetic mean of the per-ROI mean-intensity values |
| 6 | `gfp_intensity_median` | Typical reporter intensity across retained objects | Median of the per-ROI mean-intensity values |

For GFP intensity, the script preserves this input-column priority:
`Mean intensity`, `MeanIntensity`, `Mean_Intensity`, `Intensity_Mean`,
`MeanIntensity_GFP`, then `IntegratedIntensity`. If none exists, both GFP
intensity metrics are `NaN`; the other GFP calculations are still retained.

## Required inputs and directory structure

The analysis expects two independently prepared scratch workspaces whose
manifests both contain the same 94 ordered production jobs:

```text
H1_WORKSPACE/
├── manifest.json
├── mentor_field_plots/
│   └── PV_density_per_field_verified_H1.csv
└── outputs/
    └── <job>/analysis/...

REMAINING_WORKSPACE/
├── manifest.json
└── outputs/
    └── <job>/analysis/...
```

The script follows each manifest job's recorded `output_dir`, then expects:

```text
<job-output>/analysis/
├── run_manifest.json
├── 06_final_shape_filtered/
│   ├── GFP/GFP_kept_ROI_metrics.csv
│   └── PV/PV_kept_ROI_metrics.csv
└── 08_marker_colocalization/
    └── LHX6_vs_PV/visualizations/unmatched_only_qc/
        ├── LHX6_PV_all_nuclear_cells_class_counts.csv
        ├── LHX6_PV_all_nuclear_cells_classification.csv
        └── LHX6_PV_final_summary_metrics.csv
```

The verified H1 selection CSV must contain an `index` column. If it also
contains `sample`, `name`, `organoid_area_um2`, or `pv_positive_cells`, those
values are used for consistency checks against the manifest and current
outputs.

The tracked block-mapping CSV must contain:

```text
block_id,cell_line,condition,experiment_id,organoid_type
```

`condition_basis` and `inventory_note` are retained as optional audit columns.
The script requires exactly blocks 149, 150, 151, 152, 153, and 157 with the
assignments shown above.

## Dependencies

- Python 3.10 or newer;
- NumPy;
- pandas;
- Matplotlib with PDF support.

The existing Great Lakes `ims-mosaic` Conda environment is expected to provide
these packages. To create a small standalone environment instead:

```bash
conda create -n igf-reporter-plots -c conda-forge \
  python=3.11 numpy pandas matplotlib
conda activate igf-reporter-plots
```

## Command-line usage

From any directory:

```bash
python /path/to/ImageAnalysisScripts-Copy/FigureAnalysis/mNeonGreenReporterAnalysis/plot_igf_exp17_ventral_with_gfp_plots.py \
  --h1-workspace /path/to/original_workspace \
  --remaining-workspace /path/to/remaining_workspace \
  --out /path/to/new_output_directory
```

The repository mapping CSV is loaded automatically from beside the script.
The defaults for the other optional arguments are:

- `--h1-selection`:
  `<h1-workspace>/mentor_field_plots/PV_density_per_field_verified_H1.csv`;
- `--out`:
  `<h1-workspace>/mentor_field_plots/igf_exp17_ventral_only_side_by_side`.

To override either one, or to audit an alternate mapping explicitly:

```bash
python plot_igf_exp17_ventral_with_gfp_plots.py \
  --h1-workspace "$H1_WORKSPACE" \
  --remaining-workspace "$REMAINING_WORKSPACE" \
  --h1-selection /scratch/path/verified_H1.csv \
  --sample-metadata /path/to/audited_mapping.csv \
  --out "$OUTPUT_DIR"
```

Use `python plot_igf_exp17_ventral_with_gfp_plots.py --help` for the complete
CLI reference. No output directory is written by `--help`.

## Expected outputs

The output directory contains:

- `all_cell_lines_side_by_side.pdf`: 22 pages—design/assignment page, legend
  page, 14 PV/LHX6 pages, and 6 GFP pages;
- 14 PV/LHX6 PDFs and 14 matching 230-DPI PNGs using the stems in the table
  above;
- 6 GFP PDFs and 6 matching 230-DPI PNGs using the stems in the table above;
- `combined_field_cell_classes_pv_nuclear_area.csv`: one row per included
  field with nuclear counts, PV measurements, GFP measurements, mapping, and
  workspace provenance;
- `combined_pv_nuclear_positive_single_roi_areas.csv`: one row per
  nuclear-associated PV ROI, including nucleus class and source object ID;
- `sample_level_descriptive_summary.csv`: field count and selected arithmetic
  means grouped by cell line, IGF condition, experiment, and sample;
- `sample_metadata_used.csv`: exact block metadata used for the run;
- `field_selection_audit.csv`: included and excluded manifest indices, selected
  workspace, and reason.

That is 21 PDF files, 20 PNG files, and 5 CSV files. Generated outputs are not
tracked by Git.

## Running on Great Lakes

After this change is merged, use a clean checkout and keep outputs on scratch:

```bash
cd "$HOME"
if [ ! -d ImageAnalysisScripts-Copy/.git ]; then
  git clone https://github.com/oltij/ImageAnalysisScripts-Copy.git
fi
cd ImageAnalysisScripts-Copy
git switch main
git pull --ff-only origin main

conda activate ims-mosaic

REPO_ROOT="$PWD"
H1_WORKSPACE="/scratch/parent_root/parent0/$USER/<original-workspace>"
REMAINING_WORKSPACE="/scratch/parent_root/parent0/$USER/<remaining-workspace>"
OUTPUT_DIR="/scratch/parent_root/parent0/$USER/igf_exp17_ventral_reporter_figures"

python "$REPO_ROOT/FigureAnalysis/mNeonGreenReporterAnalysis/plot_igf_exp17_ventral_with_gfp_plots.py" \
  --h1-workspace "$H1_WORKSPACE" \
  --remaining-workspace "$REMAINING_WORKSPACE" \
  --out "$OUTPUT_DIR"
```

This plotting job does not need a GPU. If Great Lakes policy or load requires a
compute node, request a short standard interactive allocation before running
the same command:

```bash
salloc --account=parent_root --partition=standard \
  --nodes=1 --ntasks=1 --cpus-per-task=2 --mem=16G --time=01:00:00
srun --pty bash
conda activate ims-mosaic
```

The source workspaces are read-only inputs. Choose a new `OUTPUT_DIR` if you
need to preserve a previous plotting run.

## Retrieving results

Run this from a Mac terminal, not from the Great Lakes shell. Replace the
remote output path with the exact value used above:

```bash
LOCAL_OUTPUT="$HOME/Downloads/igf_exp17_ventral_reporter_figures"
REMOTE_OUTPUT="/scratch/parent_root/parent0/oltij/igf_exp17_ventral_reporter_figures"
mkdir -p "$LOCAL_OUTPUT"

rsync -avh --progress \
  --include='*.pdf' \
  --include='*.png' \
  --include='*.csv' \
  --exclude='*' \
  "oltij@greatlakes-xfer.arc-ts.umich.edu:${REMOTE_OUTPUT}/" \
  "$LOCAL_OUTPUT/"
```

Confirm that the five CSVs, 21 PDFs, and 20 PNGs arrived. Scratch is temporary
storage, so retrieve or archive the results promptly.
