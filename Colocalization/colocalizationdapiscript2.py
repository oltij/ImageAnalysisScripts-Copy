#!/usr/bin/env python3
"""
02_generate_colocalization_visualizations.py

Visualization/QC-only stage.

This script reads ONE self-contained HDF5 analysis file created by
01_run_colocalization.py.

Change the USER SETTINGS in this file freely without rerunning the
colocalization analysis.

Current visual outputs preserved from the original script:
    - matched-pair PDF
    - below-threshold PDF
    - filtering-stage CSV
    - filtering Sankey PDF
    - full-resolution before/after ROI overlay TIFFs
    - ROI color mapping CSVs
    - before/after ROI filtering PDF
    - kept Object-A label TIFF
    - kept Object-B label TIFF
    - two-channel label TIFF
    - RGB visualization TIFF
"""

from pathlib import Path
import io
import json

import h5py
import numpy as np
import pandas as pd
import tifffile
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Rectangle, Polygon


# ============================================================
# USER SETTINGS
# ============================================================

# THIS IS THE ONLY ANALYSIS INPUT NEEDED.
RESULTS_FILE = Path(
    "/home/oltij/Desktop/Aligned_LHX6PV_colocalization_results_25_new/"
    "Aligned_LHX6PV_colocalization_analysis.h5"
)

# Output directory for all visualization/QC files.
OUTPUT_DIR = RESULTS_FILE.parent / "visualizations"

# ------------------------------------------------------------
# Optional marker-positive AREA as a percentage of organoid area
# ------------------------------------------------------------
# This adds area-based summary metrics such as:
#     total LHX6-positive ROI area / total organoid area * 100
#     total PV-positive ROI area   / total organoid area * 100
#
# ORGANOID_MEASUREMENTS_CSV should be the CellProfiler object-level CSV
# containing the organoid ROI measurement(s). CellProfiler AreaShape_Area
# is interpreted as an area in image pixels; if multiple organoid rows are
# present, their areas are summed.
ENABLE_ORGANOID_AREA_SUMMARY = True
ORGANOID_MEASUREMENTS_CSV = Path(
    "/home/oltij/Desktop/LHX6Image/MyExpt_FilterObjects.csv"
)
ORGANOID_AREA_COLUMN = "AreaShape_Area"

# ------------------------------------------------------------
# Optional TRUE all-cell summary using PRE-EXISTING nuclear assignments
# ------------------------------------------------------------
# Keep this enabled to calculate A+/B+, A+/B-, A-/B+, and A-/B-
# as percentages of ALL nuclear ROIs, WITHOUT re-matching either marker
# to Hoechst/DAPI in this visualization script.
ENABLE_NUCLEAR_ALL_CELL_SUMMARY = True

# IMPORTANT:
# Point these to the TWO original marker-vs-nuclear HDF5 files produced by
# 01_run_colocalization.py before the current marker-vs-marker run.
#
# Example conceptual runs:
#     OBJECT_A_NUCLEAR_RESULTS_FILE = LHX6 <-> Hoechst/DAPI analysis HDF5
#     OBJECT_B_NUCLEAR_RESULTS_FILE = PV   <-> Hoechst/DAPI analysis HDF5
#
# Leave either as None only if you want the all-cell summary skipped.
OBJECT_A_NUCLEAR_RESULTS_FILE = Path(
    "/home/oltij/Desktop/Aligned_DAPILHX6_colocalization_results_25_new/"
    "Aligned_DAPILHX6_colocalization_analysis.h5"
)
OBJECT_B_NUCLEAR_RESULTS_FILE = Path(
    "/home/oltij/Desktop/Aligned_DAPIPV_colocalization_results_25_new/"
    "Aligned_DAPIPV_colocalization_analysis.h5"
)

# Nuclear denominator handling:
# - identical nuclear ObjectID sets -> use that shared set
# - one set is a strict subset of the other -> use the union/larger set
# - partially different, non-nested sets -> stop, because that suggests
#   incompatible nuclear segmentations / ID namespaces
ALLOW_NESTED_NUCLEAR_OBJECT_SETS = True

# ------------------------------------------------------------
# Matched-pair / below-threshold visualization
# ------------------------------------------------------------

N_VISUALIZATION_PAIRS = 15
RANDOM_SEED = 42
CROP_PADDING = 25

DISPLAY_LOW_PERCENTILE = 1
DISPLAY_HIGH_PERCENTILE = 99

OTHER_ROI_LINEWIDTH = 0.65
SELECTED_ROI_LINEWIDTH = 1.0

# ------------------------------------------------------------
# Full-resolution TIFF ROI overlay settings
# ------------------------------------------------------------

ROI_TIFF_OUTLINE_WIDTH = 3

# 1.0 = pure ROI color
# 0.5 = half ROI color / half fluorescence
ROI_OUTLINE_ALPHA = 1.0

ROI_COLOR_SATURATION = 0.85
ROI_COLOR_VALUE = 0.95

# ------------------------------------------------------------
# Visualization colors
# ------------------------------------------------------------

OBJECT_A_COLOR = "red"
OBJECT_B_COLOR = "blue"

OBJECT_A_OTHER_COLOR = "orange"
OBJECT_B_OTHER_COLOR = "cyan"

FILTERING_OVERLAY_MAX_FIGURE_SIZE = (18, 14)


# ============================================================
# HELPERS: HDF5
# ============================================================

def dataframe_from_h5(
    h5: h5py.File,
    name: str
) -> pd.DataFrame:
    payload = h5[f"tables/{name}"][()]
    return pd.read_csv(io.BytesIO(payload.tobytes()))


def load_results(results_file: Path):
    with h5py.File(results_file, "r") as h5:
        metadata = json.loads(
            h5["metadata/json"][()].decode("utf-8")
        )

        a_original = h5["images/object_a_original"][()]
        b_original = h5["images/object_b_original"][()]

        object_a = dataframe_from_h5(
            h5, "object_a_pixels"
        )
        object_b = dataframe_from_h5(
            h5, "object_b_pixels"
        )

        overlap_df = dataframe_from_h5(
            h5, "overlap_results"
        )
        eligible_df = dataframe_from_h5(
            h5, "eligible_pairs"
        )
        matching_df = dataframe_from_h5(
            h5, "stable_matching"
        )
        cell_id_mapping = dataframe_from_h5(
            h5, "cell_id_mapping"
        )
        unmatched_df = dataframe_from_h5(
            h5, "unmatched_object_a"
        )
        filtering_stage_df = dataframe_from_h5(
            h5, "filtering_stage_counts"
        )

        filtered_a = dataframe_from_h5(
            h5, "filtered_object_a_pixels"
        )
        filtered_b = dataframe_from_h5(
            h5, "filtered_object_b_pixels"
        )

        return {
            "metadata": metadata,
            "a_original": a_original,
            "b_original": b_original,
            "object_a": object_a,
            "object_b": object_b,
            "overlap_df": overlap_df,
            "eligible_df": eligible_df,
            "matching_df": matching_df,
            "cell_id_mapping": cell_id_mapping,
            "unmatched_df": unmatched_df,
            "filtering_stage_df": filtering_stage_df,
            "filtered_a": filtered_a,
            "filtered_b": filtered_b
        }


# ============================================================
# LOAD ANALYSIS RESULTS
# ============================================================

if not RESULTS_FILE.exists():
    raise FileNotFoundError(
        f"Analysis results file not found:\n{RESULTS_FILE}"
    )

data = load_results(RESULTS_FILE)

metadata = data["metadata"]
a_original = data["a_original"]
b_original = data["b_original"]
object_a = data["object_a"]
object_b = data["object_b"]
overlap_df = data["overlap_df"]
eligible_df = data["eligible_df"]
matching_df = data["matching_df"]
cell_id_mapping = data["cell_id_mapping"]
unmatched_df = data["unmatched_df"]
filtering_stage_df = data["filtering_stage_df"]

OBJECT_A_NAME = metadata["objects"]["object_a_name"]
OBJECT_B_NAME = metadata["objects"]["object_b_name"]
THRESHOLD = metadata["matching"]["threshold"]
ANALYSIS_VARIANT = metadata.get("analysis_variant", "unconstrained_marker_marker")
IS_NUCLEUS_CONSTRAINED = ANALYSIS_VARIANT == "nucleus_constrained_marker_marker"

PIXEL_SIZE_X_UM = metadata["pixel_size_um"]["x"]
PIXEL_SIZE_Y_UM = metadata["pixel_size_um"]["y"]

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True
)

print("=" * 70)
print("COLocalization visualization stage")
print("=" * 70)
print(f"\nAnalysis file:\n{RESULTS_FILE}")
if IS_NUCLEUS_CONSTRAINED:
    print("\nMatching variant: NUCLEUS-CONSTRAINED")
print(f"\nVisualization output directory:\n{OUTPUT_DIR}")


# ============================================================
# OBJECT COUNTS / IDS
# ============================================================

all_a_ids = set(
    object_a["ObjectID"].astype(int)
)

all_b_ids = set(
    object_b["ObjectID"].astype(int)
)

matched_a_ids = set(
    matching_df["A_ObjectID"].astype(int)
) if len(matching_df) else set()

matched_b_ids = set(
    matching_df["B_ObjectID"].astype(int)
) if len(matching_df) else set()

below_threshold_a_ids = set(
    unmatched_df["A_ObjectID"].astype(int)
) if len(unmatched_df) else set()

total_a_objects = len(all_a_ids)
total_b_objects = len(all_b_ids)


# ============================================================
# OUTPUT PATHS
# ============================================================

output_prefix = f"{OBJECT_A_NAME}_{OBJECT_B_NAME}"

visualization_pdf = OUTPUT_DIR / (
    f"{output_prefix}_matched_pair_visualization.pdf"
)

below_threshold_pdf = OUTPUT_DIR / (
    f"{output_prefix}_below_threshold_visualization.pdf"
)

filtering_sankey_pdf = OUTPUT_DIR / (
    f"{output_prefix}_filtering_sankey.pdf"
)

roi_filtering_overlay_pdf = OUTPUT_DIR / (
    f"{output_prefix}_ROI_filtering_before_after.pdf"
)

filtering_stage_csv = OUTPUT_DIR / (
    f"{output_prefix}_filtering_stage_counts.csv"
)

object_a_tiff = OUTPUT_DIR / (
    f"kept_{OBJECT_A_NAME}.tif"
)

object_b_tiff = OUTPUT_DIR / (
    f"kept_{OBJECT_B_NAME}.tif"
)

combined_tiff = OUTPUT_DIR / (
    f"kept_{output_prefix}_two_channel.tif"
)

rgb_visualization_tiff = OUTPUT_DIR / (
    f"kept_{output_prefix}_RGB.tif"
)

roi_filtering_overlay_a_before_tiff = OUTPUT_DIR / (
    f"{OBJECT_A_NAME}_ROI_before_filtering.tif"
)

roi_filtering_overlay_a_after_tiff = OUTPUT_DIR / (
    f"{OBJECT_A_NAME}_ROI_after_filtering.tif"
)

roi_filtering_overlay_b_before_tiff = OUTPUT_DIR / (
    f"{OBJECT_B_NAME}_ROI_before_filtering.tif"
)

roi_filtering_overlay_b_after_tiff = OUTPUT_DIR / (
    f"{OBJECT_B_NAME}_ROI_after_filtering.tif"
)

a_color_mapping_csv = OUTPUT_DIR / (
    f"{OBJECT_A_NAME}_ROI_color_mapping.csv"
)

b_color_mapping_csv = OUTPUT_DIR / (
    f"{OBJECT_B_NAME}_ROI_color_mapping.csv"
)


# ============================================================
# IMAGE SIZE / DISPLAY CONTRAST
# ============================================================

image_height, image_width = a_original.shape

a_vmin = np.percentile(
    a_original,
    DISPLAY_LOW_PERCENTILE
)
a_vmax = np.percentile(
    a_original,
    DISPLAY_HIGH_PERCENTILE
)

b_vmin = np.percentile(
    b_original,
    DISPLAY_LOW_PERCENTILE
)
b_vmax = np.percentile(
    b_original,
    DISPLAY_HIGH_PERCENTILE
)

print("\nDisplay contrast:")
print(
    f"{OBJECT_A_NAME}: "
    f"{a_vmin:.3f} to {a_vmax:.3f}"
)
print(
    f"{OBJECT_B_NAME}: "
    f"{b_vmin:.3f} to {b_vmax:.3f}"
)


# ============================================================
# PIXEL GROUPS
# ============================================================

a_pixel_groups = {
    int(object_id): group[["X", "Y"]].to_numpy()
    for object_id, group in object_a.groupby("ObjectID")
}

b_pixel_groups = {
    int(object_id): group[["X", "Y"]].to_numpy()
    for object_id, group in object_b.groupby("ObjectID")
}


# ============================================================
# ROI OUTLINE HELPER FOR CROPPED PDF VIEWS
# ============================================================

def make_roi_outline(
    roi_pixels,
    crop_xmin,
    crop_ymin,
    crop_xmax,
    crop_ymax
):
    crop_width = crop_xmax - crop_xmin
    crop_height = crop_ymax - crop_ymin

    mask = np.zeros(
        (crop_height, crop_width),
        dtype=bool
    )

    for x, y in roi_pixels:
        local_x = int(x) - crop_xmin
        local_y = int(y) - crop_ymin

        if (
            0 <= local_x < crop_width
            and 0 <= local_y < crop_height
        ):
            mask[local_y, local_x] = True

    boundary = np.zeros_like(mask)

    boundary[1:, :] |= (
        mask[1:, :] & ~mask[:-1, :]
    )
    boundary[:-1, :] |= (
        mask[:-1, :] & ~mask[1:, :]
    )
    boundary[:, 1:] |= (
        mask[:, 1:] & ~mask[:, :-1]
    )
    boundary[:, :-1] |= (
        mask[:, :-1] & ~mask[:, 1:]
    )

    return boundary


# ============================================================
# FIND OBJECTS PRESENT IN CROP
# ============================================================

def get_rois_in_crop(
    pixel_groups,
    xmin,
    ymin,
    xmax,
    ymax
):
    roi_ids = []

    for roi_id, pixels in pixel_groups.items():
        xs = pixels[:, 0]
        ys = pixels[:, 1]

        overlaps_x = (
            (xs >= xmin) &
            (xs < xmax)
        )
        overlaps_y = (
            (ys >= ymin) &
            (ys < ymax)
        )

        if np.any(overlaps_x & overlaps_y):
            roi_ids.append(roi_id)

    return roi_ids


# ============================================================
# MATCHED-PAIR PDF
# ============================================================

if len(matching_df) > 0:

    rng = np.random.default_rng(RANDOM_SEED)

    n_to_show = min(
        N_VISUALIZATION_PAIRS,
        len(matching_df)
    )

    selected_indices = rng.choice(
        len(matching_df),
        size=n_to_show,
        replace=False
    )

    visualization_pairs = (
        matching_df
        .iloc[selected_indices]
        .reset_index(drop=True)
    )

    print(
        f"\nCreating matched visualization "
        f"for {n_to_show} random pairs."
    )

    with PdfPages(visualization_pdf) as pdf:

        fig = plt.figure(figsize=(30, 14))

        gs = fig.add_gridspec(
            3,
            n_to_show,
            hspace=0.35,
            wspace=0.08
        )

        for column_number, (_, match) in enumerate(
            visualization_pairs.iterrows()
        ):
            a_id = int(match["A_ObjectID"])
            b_id = int(match["B_ObjectID"])

            overlap_percent = float(
                match["A_Percent_Inside_B"]
            )

            cell_id = int(
                cell_id_mapping.loc[
                    (
                        cell_id_mapping["A_ObjectID"]
                        == a_id
                    )
                    &
                    (
                        cell_id_mapping["B_ObjectID"]
                        == b_id
                    ),
                    "CellID"
                ].iloc[0]
            )

            a_pixels = a_pixel_groups[a_id]
            b_pixels = b_pixel_groups[b_id]

            all_x = np.concatenate(
                [a_pixels[:, 0], b_pixels[:, 0]]
            )
            all_y = np.concatenate(
                [a_pixels[:, 1], b_pixels[:, 1]]
            )

            xmin = max(
                0,
                int(np.min(all_x)) - CROP_PADDING
            )
            xmax = min(
                image_width,
                int(np.max(all_x))
                + CROP_PADDING + 1
            )
            ymin = max(
                0,
                int(np.min(all_y)) - CROP_PADDING
            )
            ymax = min(
                image_height,
                int(np.max(all_y))
                + CROP_PADDING + 1
            )

            a_crop = a_original[
                ymin:ymax,
                xmin:xmax
            ]
            b_crop = b_original[
                ymin:ymax,
                xmin:xmax
            ]

            a_rois_in_crop = get_rois_in_crop(
                a_pixel_groups,
                xmin,
                ymin,
                xmax,
                ymax
            )

            b_rois_in_crop = get_rois_in_crop(
                b_pixel_groups,
                xmin,
                ymin,
                xmax,
                ymax
            )

            a_other_outlines = []

            for other_a_id in a_rois_in_crop:
                if other_a_id == a_id:
                    continue

                a_other_outlines.append(
                    make_roi_outline(
                        a_pixel_groups[other_a_id],
                        xmin,
                        ymin,
                        xmax,
                        ymax
                    )
                )

            selected_a_outline = make_roi_outline(
                a_pixels,
                xmin,
                ymin,
                xmax,
                ymax
            )

            b_other_outlines = []

            for other_b_id in b_rois_in_crop:
                if other_b_id == b_id:
                    continue

                b_other_outlines.append(
                    make_roi_outline(
                        b_pixel_groups[other_b_id],
                        xmin,
                        ymin,
                        xmax,
                        ymax
                    )
                )

            selected_b_outline = make_roi_outline(
                b_pixels,
                xmin,
                ymin,
                xmax,
                ymax
            )

            # ------------------------------------------------
            # ROW 1 — B
            # ------------------------------------------------

            ax_b = fig.add_subplot(
                gs[0, column_number]
            )

            ax_b.imshow(
                b_crop,
                cmap="gray",
                vmin=b_vmin,
                vmax=b_vmax
            )

            for outline in b_other_outlines:
                ax_b.contour(
                    outline,
                    levels=[0.5],
                    linewidths=OTHER_ROI_LINEWIDTH,
                    colors=OBJECT_B_OTHER_COLOR
                )

            ax_b.contour(
                selected_b_outline,
                levels=[0.5],
                linewidths=SELECTED_ROI_LINEWIDTH,
                colors=OBJECT_B_COLOR
            )

            ax_b.set_title(
                f"{OBJECT_B_NAME} {b_id}",
                fontsize=10,
                fontweight="bold"
            )
            ax_b.axis("off")

            # ------------------------------------------------
            # ROW 2 — A
            # ------------------------------------------------

            ax_a = fig.add_subplot(
                gs[1, column_number]
            )

            ax_a.imshow(
                a_crop,
                cmap="gray",
                vmin=a_vmin,
                vmax=a_vmax
            )

            for outline in a_other_outlines:
                ax_a.contour(
                    outline,
                    levels=[0.5],
                    linewidths=OTHER_ROI_LINEWIDTH,
                    colors=OBJECT_A_OTHER_COLOR
                )

            ax_a.contour(
                selected_a_outline,
                levels=[0.5],
                linewidths=SELECTED_ROI_LINEWIDTH,
                colors=OBJECT_A_COLOR
            )

            ax_a.set_title(
                f"{OBJECT_A_NAME} {a_id}",
                fontsize=10,
                fontweight="bold"
            )
            ax_a.axis("off")

            # ------------------------------------------------
            # ROW 3 — COMBINED
            # ------------------------------------------------

            ax_combined = fig.add_subplot(
                gs[2, column_number]
            )

            a_range = a_vmax - a_vmin
            b_range = b_vmax - b_vmin

            if a_range == 0:
                a_range = 1
            if b_range == 0:
                b_range = 1

            a_norm = (
                a_crop.astype(float) - a_vmin
            ) / a_range

            b_norm = (
                b_crop.astype(float) - b_vmin
            ) / b_range

            a_norm = np.clip(a_norm, 0, 1)
            b_norm = np.clip(b_norm, 0, 1)

            combined = np.zeros(
                (*a_crop.shape, 3),
                dtype=float
            )

            combined[..., 0] = a_norm
            combined[..., 1] = b_norm

            ax_combined.imshow(combined)

            for outline in a_other_outlines:
                ax_combined.contour(
                    outline,
                    levels=[0.5],
                    linewidths=OTHER_ROI_LINEWIDTH,
                    colors=OBJECT_A_OTHER_COLOR
                )

            for outline in b_other_outlines:
                ax_combined.contour(
                    outline,
                    levels=[0.5],
                    linewidths=OTHER_ROI_LINEWIDTH,
                    colors=OBJECT_B_OTHER_COLOR
                )

            ax_combined.contour(
                selected_a_outline,
                levels=[0.5],
                linewidths=SELECTED_ROI_LINEWIDTH,
                colors=OBJECT_A_COLOR
            )

            ax_combined.contour(
                selected_b_outline,
                levels=[0.5],
                linewidths=SELECTED_ROI_LINEWIDTH,
                colors=OBJECT_B_COLOR
            )

            center_x = np.mean(all_x) - xmin
            center_y = np.mean(all_y) - ymin

            ax_combined.text(
                center_x,
                center_y,
                f"Cell {cell_id}",
                color="white",
                fontsize=8,
                fontweight="bold",
                ha="center",
                va="center",
                bbox=dict(
                    facecolor="black",
                    alpha=0.55,
                    edgecolor="none",
                    pad=1.5
                )
            )

            ax_combined.set_title(
                f"{overlap_percent:.1f}% A inside B",
                fontsize=9,
                fontweight="bold"
            )

            ax_combined.axis("off")

        fig.suptitle(
            (
                "NUCLEUS-CONSTRAINED: "
                if IS_NUCLEUS_CONSTRAINED else ""
            )
            + f"Random {OBJECT_A_NAME}-{OBJECT_B_NAME} Matched Cells",
            fontsize=18,
            fontweight="bold",
            y=0.995
        )

        fig.text(
            0.01,
            0.83,
            OBJECT_B_NAME,
            fontsize=13,
            fontweight="bold",
            rotation=90,
            va="center"
        )

        fig.text(
            0.01,
            0.50,
            OBJECT_A_NAME,
            fontsize=13,
            fontweight="bold",
            rotation=90,
            va="center"
        )

        fig.text(
            0.01,
            0.17,
            "Combined",
            fontsize=13,
            fontweight="bold",
            rotation=90,
            va="center"
        )

        fig.text(
            0.5,
            0.015,
            f"Selected {OBJECT_A_NAME} = {OBJECT_A_COLOR} | "
            f"Other {OBJECT_A_NAME} = {OBJECT_A_OTHER_COLOR} | "
            f"Selected {OBJECT_B_NAME} = {OBJECT_B_COLOR} | "
            f"Other {OBJECT_B_NAME} = {OBJECT_B_OTHER_COLOR}",
            ha="center",
            fontsize=10
        )

        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

else:
    print("\nNo matched pairs available.")
    if IS_NUCLEUS_CONSTRAINED:
        with PdfPages(visualization_pdf) as pdf:
            fig = plt.figure(figsize=(11, 6.5))
            fig.text(
                0.08,
                0.78,
                f"{OBJECT_A_NAME}–{OBJECT_B_NAME}: NUCLEUS-CONSTRAINED matching",
                fontsize=18,
                fontweight="bold",
            )
            fig.text(
                0.08,
                0.56,
                "No final stable-matched pairs were produced.",
                fontsize=14,
            )
            fig.text(
                0.08,
                0.36,
                "See analysis.h5 and nuclear_constraint_pair_audit.csv for explicit\n"
                "overlap-threshold, canonical-nuclear-identity, and exclusion provenance.",
                fontsize=11,
            )
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)


# ============================================================
# BELOW-THRESHOLD VISUALIZATION
# ============================================================

print("\nCreating below-threshold visualization...")

n_failed_to_show = min(
    N_VISUALIZATION_PAIRS,
    len(below_threshold_a_ids)
)

if n_failed_to_show > 0:

    rng_failed = np.random.default_rng(
        RANDOM_SEED + 1
    )

    selected_failed_a_ids = rng_failed.choice(
        list(below_threshold_a_ids),
        size=n_failed_to_show,
        replace=False
    )

    selected_failed_a_ids = [
        int(x) for x in selected_failed_a_ids
    ]

    # Re-create the best partial match lookup directly from analysis data.
    best_partial = {}
    if len(overlap_df) > 0:
        best_rows = (
            overlap_df
            .sort_values(
                ["A_ObjectID", "A_Percent_Inside_B", "B_ObjectID"],
                ascending=[True, False, True]
            )
            .drop_duplicates(
                subset=["A_ObjectID"],
                keep="first"
            )
        )
        for _, row in best_rows.iterrows():
            best_partial[int(row["A_ObjectID"])] = row

    print(
        f"Creating visualization for "
        f"{n_failed_to_show} below-threshold "
        f"{OBJECT_A_NAME} objects."
    )

    with PdfPages(below_threshold_pdf) as pdf:

        fig = plt.figure(figsize=(30, 14))

        gs = fig.add_gridspec(
            3,
            n_failed_to_show,
            hspace=0.35,
            wspace=0.08
        )

        for column_number, a_id in enumerate(
            selected_failed_a_ids
        ):
            best_row = best_partial.get(a_id)

            if best_row is not None:
                b_id = int(best_row["B_ObjectID"])
                overlap_percent = float(
                    best_row["A_Percent_Inside_B"]
                )
                has_b_overlap = True
            else:
                b_id = None
                overlap_percent = 0.0
                has_b_overlap = False

            a_pixels = a_pixel_groups[a_id]

            if has_b_overlap:
                b_pixels = b_pixel_groups[b_id]

                all_x = np.concatenate(
                    [a_pixels[:, 0], b_pixels[:, 0]]
                )
                all_y = np.concatenate(
                    [a_pixels[:, 1], b_pixels[:, 1]]
                )
            else:
                b_pixels = np.empty((0, 2), dtype=int)
                all_x = a_pixels[:, 0]
                all_y = a_pixels[:, 1]

            xmin = max(
                0,
                int(np.min(all_x)) - CROP_PADDING
            )
            xmax = min(
                image_width,
                int(np.max(all_x)) + CROP_PADDING + 1
            )
            ymin = max(
                0,
                int(np.min(all_y)) - CROP_PADDING
            )
            ymax = min(
                image_height,
                int(np.max(all_y)) + CROP_PADDING + 1
            )

            a_crop = a_original[
                ymin:ymax,
                xmin:xmax
            ]
            b_crop = b_original[
                ymin:ymax,
                xmin:xmax
            ]

            a_rois_in_crop = get_rois_in_crop(
                a_pixel_groups,
                xmin,
                ymin,
                xmax,
                ymax
            )

            b_rois_in_crop = get_rois_in_crop(
                b_pixel_groups,
                xmin,
                ymin,
                xmax,
                ymax
            )

            a_other_outlines = []

            for other_a_id in a_rois_in_crop:
                if other_a_id == a_id:
                    continue

                a_other_outlines.append(
                    make_roi_outline(
                        a_pixel_groups[other_a_id],
                        xmin,
                        ymin,
                        xmax,
                        ymax
                    )
                )

            selected_a_outline = make_roi_outline(
                a_pixels,
                xmin,
                ymin,
                xmax,
                ymax
            )

            b_other_outlines = []

            for other_b_id in b_rois_in_crop:
                if has_b_overlap and other_b_id == b_id:
                    continue

                b_other_outlines.append(
                    make_roi_outline(
                        b_pixel_groups[other_b_id],
                        xmin,
                        ymin,
                        xmax,
                        ymax
                    )
                )

            if has_b_overlap:
                selected_b_outline = make_roi_outline(
                    b_pixels,
                    xmin,
                    ymin,
                    xmax,
                    ymax
                )
            else:
                selected_b_outline = np.zeros(
                    a_crop.shape,
                    dtype=bool
                )

            # Row 1 — B
            ax_b = fig.add_subplot(
                gs[0, column_number]
            )

            ax_b.imshow(
                b_crop,
                cmap="gray",
                vmin=b_vmin,
                vmax=b_vmax
            )

            for outline in b_other_outlines:
                ax_b.contour(
                    outline,
                    levels=[0.5],
                    linewidths=OTHER_ROI_LINEWIDTH,
                    colors=OBJECT_B_OTHER_COLOR
                )

            if has_b_overlap:
                ax_b.contour(
                    selected_b_outline,
                    levels=[0.5],
                    linewidths=SELECTED_ROI_LINEWIDTH,
                    colors=OBJECT_B_COLOR
                )

                ax_b.set_title(
                    f"{OBJECT_B_NAME} {b_id}\n"
                    f"best partial overlap",
                    fontsize=9,
                    fontweight="bold"
                )
            else:
                ax_b.set_title(
                    f"No {OBJECT_B_NAME} overlap",
                    fontsize=9,
                    fontweight="bold"
                )

            ax_b.axis("off")

            # Row 2 — A
            ax_a = fig.add_subplot(
                gs[1, column_number]
            )

            ax_a.imshow(
                a_crop,
                cmap="gray",
                vmin=a_vmin,
                vmax=a_vmax
            )

            for outline in a_other_outlines:
                ax_a.contour(
                    outline,
                    levels=[0.5],
                    linewidths=OTHER_ROI_LINEWIDTH,
                    colors=OBJECT_A_OTHER_COLOR
                )

            ax_a.contour(
                selected_a_outline,
                levels=[0.5],
                linewidths=SELECTED_ROI_LINEWIDTH,
                colors=OBJECT_A_COLOR
            )

            ax_a.set_title(
                f"{OBJECT_A_NAME} {a_id}",
                fontsize=10,
                fontweight="bold"
            )

            ax_a.axis("off")

            # Row 3 — Combined
            ax_combined = fig.add_subplot(
                gs[2, column_number]
            )

            a_range = a_vmax - a_vmin
            b_range = b_vmax - b_vmin

            if a_range == 0:
                a_range = 1
            if b_range == 0:
                b_range = 1

            a_norm = np.clip(
                (a_crop.astype(float) - a_vmin) / a_range,
                0,
                1
            )

            b_norm = np.clip(
                (b_crop.astype(float) - b_vmin) / b_range,
                0,
                1
            )

            combined = np.zeros(
                (*a_crop.shape, 3),
                dtype=float
            )

            combined[..., 0] = a_norm
            combined[..., 1] = b_norm

            ax_combined.imshow(combined)

            for outline in a_other_outlines:
                ax_combined.contour(
                    outline,
                    levels=[0.5],
                    linewidths=OTHER_ROI_LINEWIDTH,
                    colors=OBJECT_A_OTHER_COLOR
                )

            for outline in b_other_outlines:
                ax_combined.contour(
                    outline,
                    levels=[0.5],
                    linewidths=OTHER_ROI_LINEWIDTH,
                    colors=OBJECT_B_OTHER_COLOR
                )

            ax_combined.contour(
                selected_a_outline,
                levels=[0.5],
                linewidths=SELECTED_ROI_LINEWIDTH,
                colors=OBJECT_A_COLOR
            )

            if has_b_overlap:
                ax_combined.contour(
                    selected_b_outline,
                    levels=[0.5],
                    linewidths=SELECTED_ROI_LINEWIDTH,
                    colors=OBJECT_B_COLOR
                )

            center_x = np.mean(all_x) - xmin
            center_y = np.mean(all_y) - ymin

            if has_b_overlap:
                label_text = (
                    f"{OBJECT_A_NAME} {a_id}\n"
                    f"{OBJECT_B_NAME} {b_id}\n"
                    f"{overlap_percent:.1f}% A in B"
                )

                title_text = (
                    f"Below {THRESHOLD * 100:.0f}%: "
                    f"{overlap_percent:.1f}%"
                )
            else:
                label_text = (
                    f"{OBJECT_A_NAME} {a_id}\n"
                    f"No {OBJECT_B_NAME} overlap"
                )

                title_text = "No overlap"

            ax_combined.text(
                center_x,
                center_y,
                label_text,
                color="white",
                fontsize=7,
                fontweight="bold",
                ha="center",
                va="center",
                bbox=dict(
                    facecolor="black",
                    alpha=0.55,
                    edgecolor="none",
                    pad=1.5
                )
            )

            ax_combined.set_title(
                title_text,
                fontsize=9,
                fontweight="bold"
            )

            ax_combined.axis("off")

        fig.suptitle(
            f"Random Below-Threshold {OBJECT_A_NAME} Objects",
            fontsize=18,
            fontweight="bold",
            y=0.995
        )

        fig.text(
            0.01,
            0.83,
            OBJECT_B_NAME,
            fontsize=13,
            fontweight="bold",
            rotation=90,
            va="center"
        )

        fig.text(
            0.01,
            0.50,
            OBJECT_A_NAME,
            fontsize=13,
            fontweight="bold",
            rotation=90,
            va="center"
        )

        fig.text(
            0.01,
            0.17,
            "Combined",
            fontsize=13,
            fontweight="bold",
            rotation=90,
            va="center"
        )

        fig.text(
            0.5,
            0.015,
            f"Selected {OBJECT_A_NAME} = {OBJECT_A_COLOR} | "
            f"Other {OBJECT_A_NAME} = {OBJECT_A_OTHER_COLOR} | "
            f"Best {OBJECT_B_NAME} = {OBJECT_B_COLOR} | "
            f"Other {OBJECT_B_NAME} = {OBJECT_B_OTHER_COLOR}",
            ha="center",
            fontsize=10
        )

        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

else:
    print(
        f"No below-threshold {OBJECT_A_NAME} objects "
        f"available for visualization."
    )



# ============================================================
# 5%-INCREMENT HOECHST-AREA-OVERLAP THRESHOLD QC
# ============================================================
#
# This module is specifically for evaluating the colocalization
# threshold based on the percentage of Object A (Hoechst) area
# that lies inside Object B (PV).
#
# Each bin is 5 percentage points wide:
#     0–5%, 5–10%, ..., 95–100%
#
# ABOVE-THRESHOLD:
#     Uses the FINAL STABLE-MATCHED pairs.
#
# BELOW-THRESHOLD:
#     Uses each failed Hoechst ROI's BEST PARTIAL PV overlap.
#
# For each 5% bin, up to MAX_OVERLAP_QC_PAIRS examples are shown.
# Each example is a 3-row panel:
#     row 1 = PV
#     row 2 = Hoechst
#     row 3 = combined
#
# This section is intentionally visualization-only. Changing any
# setting below does NOT rerun the colocalization analysis.
# ============================================================

# ------------------------------------------------------------
# USER SETTINGS FOR THRESHOLD QC
# ------------------------------------------------------------

OVERLAP_QC_BIN_WIDTH = 5
MAX_OVERLAP_QC_PAIRS = 5
OVERLAP_QC_RANDOM_SEED = 42

OVERLAP_QC_FIGSIZE = (24, 11)

# Keep the same crop/contrast behavior as the ordinary pair PDF.
OVERLAP_QC_CROP_PADDING = CROP_PADDING

# ------------------------------------------------------------
# Validate bin width
# ------------------------------------------------------------

N_OVERLAP_QC_BINS = int(
    round(100 / OVERLAP_QC_BIN_WIDTH)
)

if (
    N_OVERLAP_QC_BINS <= 0
    or 100 % OVERLAP_QC_BIN_WIDTH != 0
):
    raise ValueError(
        "OVERLAP_QC_BIN_WIDTH must divide evenly into 100."
    )


# ------------------------------------------------------------
# Output files
# ------------------------------------------------------------

overlap_qc_above_pdf = OUTPUT_DIR / (
    f"{output_prefix}_Hoechst_area_overlap_5pct_above_threshold.pdf"
)

overlap_qc_below_pdf = OUTPUT_DIR / (
    f"{output_prefix}_Hoechst_area_overlap_5pct_below_threshold.pdf"
)

overlap_qc_bin_summary_csv = OUTPUT_DIR / (
    f"{output_prefix}_Hoechst_area_overlap_5pct_bin_summary.csv"
)


# ------------------------------------------------------------
# Helper: create one matched-pair figure page
# ------------------------------------------------------------

def make_overlap_qc_page(
    pair_rows,
    pdf,
    page_title,
    seed_offset=0
):
    """
    Create one PDF page containing up to MAX_OVERLAP_QC_PAIRS
    pair examples.

    Each pair occupies one column with:
        PV / Hoechst / Combined

    pair_rows must contain:
        A_ObjectID
        B_ObjectID
        A_Percent_Inside_B
    """

    if len(pair_rows) == 0:
        return

    n_to_show = min(
        MAX_OVERLAP_QC_PAIRS,
        len(pair_rows)
    )

    # Deterministically sample from the bin.
    rng = np.random.default_rng(
        OVERLAP_QC_RANDOM_SEED + seed_offset
    )

    if len(pair_rows) > n_to_show:
        selected_indices = rng.choice(
            len(pair_rows),
            size=n_to_show,
            replace=False
        )

        page_pairs = (
            pair_rows
            .iloc[selected_indices]
            .reset_index(drop=True)
        )
    else:
        page_pairs = (
            pair_rows
            .reset_index(drop=True)
        )

    fig = plt.figure(
        figsize=OVERLAP_QC_FIGSIZE
    )

    gs = fig.add_gridspec(
        3,
        n_to_show,
        hspace=0.38,
        wspace=0.10
    )

    for column_number, (_, pair) in enumerate(
        page_pairs.iterrows()
    ):
        a_id = int(pair["A_ObjectID"])
        b_id = int(pair["B_ObjectID"])

        overlap_percent = float(
            pair["A_Percent_Inside_B"]
        )

        a_pixels = a_pixel_groups.get(a_id)
        b_pixels = b_pixel_groups.get(b_id)

        if a_pixels is None or b_pixels is None:
            continue

        all_x = np.concatenate(
            [
                a_pixels[:, 0],
                b_pixels[:, 0]
            ]
        )

        all_y = np.concatenate(
            [
                a_pixels[:, 1],
                b_pixels[:, 1]
            ]
        )

        xmin = max(
            0,
            int(np.min(all_x))
            - OVERLAP_QC_CROP_PADDING
        )

        xmax = min(
            image_width,
            int(np.max(all_x))
            + OVERLAP_QC_CROP_PADDING
            + 1
        )

        ymin = max(
            0,
            int(np.min(all_y))
            - OVERLAP_QC_CROP_PADDING
        )

        ymax = min(
            image_height,
            int(np.max(all_y))
            + OVERLAP_QC_CROP_PADDING
            + 1
        )

        a_crop = a_original[
            ymin:ymax,
            xmin:xmax
        ]

        b_crop = b_original[
            ymin:ymax,
            xmin:xmax
        ]

        a_rois_in_crop = get_rois_in_crop(
            a_pixel_groups,
            xmin,
            ymin,
            xmax,
            ymax
        )

        b_rois_in_crop = get_rois_in_crop(
            b_pixel_groups,
            xmin,
            ymin,
            xmax,
            ymax
        )

        # ----------------------------------------------------
        # Other ROI outlines
        # ----------------------------------------------------

        a_other_outlines = []

        for other_a_id in a_rois_in_crop:
            if other_a_id == a_id:
                continue

            a_other_outlines.append(
                make_roi_outline(
                    a_pixel_groups[other_a_id],
                    xmin,
                    ymin,
                    xmax,
                    ymax
                )
            )

        b_other_outlines = []

        for other_b_id in b_rois_in_crop:
            if other_b_id == b_id:
                continue

            b_other_outlines.append(
                make_roi_outline(
                    b_pixel_groups[other_b_id],
                    xmin,
                    ymin,
                    xmax,
                    ymax
                )
            )

        selected_a_outline = make_roi_outline(
            a_pixels,
            xmin,
            ymin,
            xmax,
            ymax
        )

        selected_b_outline = make_roi_outline(
            b_pixels,
            xmin,
            ymin,
            xmax,
            ymax
        )

        # ----------------------------------------------------
        # ROW 1 — PV
        # ----------------------------------------------------

        ax_b = fig.add_subplot(
            gs[0, column_number]
        )

        ax_b.imshow(
            b_crop,
            cmap="gray",
            vmin=b_vmin,
            vmax=b_vmax
        )

        for outline in b_other_outlines:
            ax_b.contour(
                outline,
                levels=[0.5],
                linewidths=OTHER_ROI_LINEWIDTH,
                colors=OBJECT_B_OTHER_COLOR
            )

        ax_b.contour(
            selected_b_outline,
            levels=[0.5],
            linewidths=SELECTED_ROI_LINEWIDTH,
            colors=OBJECT_B_COLOR
        )

        ax_b.set_title(
            f"{OBJECT_B_NAME} {b_id}",
            fontsize=9,
            fontweight="bold"
        )

        ax_b.axis("off")

        # ----------------------------------------------------
        # ROW 2 — HOECHST
        # ----------------------------------------------------

        ax_a = fig.add_subplot(
            gs[1, column_number]
        )

        ax_a.imshow(
            a_crop,
            cmap="gray",
            vmin=a_vmin,
            vmax=a_vmax
        )

        for outline in a_other_outlines:
            ax_a.contour(
                outline,
                levels=[0.5],
                linewidths=OTHER_ROI_LINEWIDTH,
                colors=OBJECT_A_OTHER_COLOR
            )

        ax_a.contour(
            selected_a_outline,
            levels=[0.5],
            linewidths=SELECTED_ROI_LINEWIDTH,
            colors=OBJECT_A_COLOR
        )

        ax_a.set_title(
            f"{OBJECT_A_NAME} {a_id}",
            fontsize=9,
            fontweight="bold"
        )

        ax_a.axis("off")

        # ----------------------------------------------------
        # ROW 3 — COMBINED
        # ----------------------------------------------------

        ax_combined = fig.add_subplot(
            gs[2, column_number]
        )

        a_range = a_vmax - a_vmin
        b_range = b_vmax - b_vmin

        if a_range == 0:
            a_range = 1

        if b_range == 0:
            b_range = 1

        a_norm = (
            a_crop.astype(float)
            - a_vmin
        ) / a_range

        b_norm = (
            b_crop.astype(float)
            - b_vmin
        ) / b_range

        a_norm = np.clip(
            a_norm,
            0,
            1
        )

        b_norm = np.clip(
            b_norm,
            0,
            1
        )

        combined = np.zeros(
            (
                a_crop.shape[0],
                a_crop.shape[1],
                3
            ),
            dtype=float
        )

        # Preserve the existing combined display:
        # Hoechst = red, PV = green.
        combined[..., 0] = a_norm
        combined[..., 1] = b_norm

        ax_combined.imshow(
            combined
        )

        for outline in a_other_outlines:
            ax_combined.contour(
                outline,
                levels=[0.5],
                linewidths=OTHER_ROI_LINEWIDTH,
                colors=OBJECT_A_OTHER_COLOR
            )

        for outline in b_other_outlines:
            ax_combined.contour(
                outline,
                levels=[0.5],
                linewidths=OTHER_ROI_LINEWIDTH,
                colors=OBJECT_B_OTHER_COLOR
            )

        ax_combined.contour(
            selected_a_outline,
            levels=[0.5],
            linewidths=SELECTED_ROI_LINEWIDTH,
            colors=OBJECT_A_COLOR
        )

        ax_combined.contour(
            selected_b_outline,
            levels=[0.5],
            linewidths=SELECTED_ROI_LINEWIDTH,
            colors=OBJECT_B_COLOR
        )

        center_x = np.mean(all_x) - xmin
        center_y = np.mean(all_y) - ymin

        ax_combined.text(
            center_x,
            center_y,
            f"{overlap_percent:.1f}%",
            color="white",
            fontsize=9,
            fontweight="bold",
            ha="center",
            va="center",
            bbox=dict(
                facecolor="black",
                alpha=0.55,
                edgecolor="none",
                pad=1.5
            )
        )

        ax_combined.set_title(
            f"{overlap_percent:.1f}% "
            f"{OBJECT_A_NAME} area in {OBJECT_B_NAME}",
            fontsize=8.5,
            fontweight="bold"
        )

        ax_combined.axis("off")

    fig.suptitle(
        page_title,
        fontsize=17,
        fontweight="bold",
        y=0.995
    )

    fig.text(
        0.01,
        0.83,
        OBJECT_B_NAME,
        fontsize=12,
        fontweight="bold",
        rotation=90,
        va="center"
    )

    fig.text(
        0.01,
        0.50,
        OBJECT_A_NAME,
        fontsize=12,
        fontweight="bold",
        rotation=90,
        va="center"
    )

    fig.text(
        0.01,
        0.17,
        "Combined",
        fontsize=12,
        fontweight="bold",
        rotation=90,
        va="center"
    )

    fig.text(
        0.5,
        0.018,
        (
            f"Selected {OBJECT_A_NAME} = {OBJECT_A_COLOR} | "
            f"Selected {OBJECT_B_NAME} = {OBJECT_B_COLOR} | "
            f"Other {OBJECT_A_NAME} = {OBJECT_A_OTHER_COLOR} | "
            f"Other {OBJECT_B_NAME} = {OBJECT_B_OTHER_COLOR}"
        ),
        ha="center",
        fontsize=9.5
    )

    pdf.savefig(
        fig,
        bbox_inches="tight"
    )

    plt.close(fig)


# ------------------------------------------------------------
# Build ABOVE-THRESHOLD bin table
# ------------------------------------------------------------

above_qc_rows = []

if len(matching_df) > 0:

    above_qc = matching_df.copy()

    above_qc["A_Percent_Inside_B"] = pd.to_numeric(
        above_qc["A_Percent_Inside_B"],
        errors="coerce"
    )

    above_qc = above_qc[
        above_qc["A_Percent_Inside_B"].notna()
    ].copy()

    above_qc["OverlapBinIndex"] = np.floor(
        above_qc["A_Percent_Inside_B"]
        / OVERLAP_QC_BIN_WIDTH
    ).astype(int)

    # 100% belongs to the final 95–100% bin.
    above_qc.loc[
        above_qc["OverlapBinIndex"] >= N_OVERLAP_QC_BINS,
        "OverlapBinIndex"
    ] = N_OVERLAP_QC_BINS - 1

    above_qc["OverlapBinStart"] = (
        above_qc["OverlapBinIndex"]
        * OVERLAP_QC_BIN_WIDTH
    )

    above_qc["OverlapBinEnd"] = (
        (
            above_qc["OverlapBinIndex"] + 1
        )
        * OVERLAP_QC_BIN_WIDTH
    )

    above_qc["OverlapBinLabel"] = (
        above_qc["OverlapBinStart"].astype(str)
        + "–"
        + above_qc["OverlapBinEnd"].astype(str)
        + "%"
    )

    for bin_index in range(N_OVERLAP_QC_BINS):

        bin_start = (
            bin_index * OVERLAP_QC_BIN_WIDTH
        )

        bin_end = (
            bin_start + OVERLAP_QC_BIN_WIDTH
        )

        in_bin = above_qc[
            above_qc["OverlapBinIndex"] == bin_index
        ]

        above_qc_rows.append(
            {
                "View": "Above threshold",
                "BinIndex": bin_index,
                "BinStartPercent": bin_start,
                "BinEndPercent": bin_end,
                "AvailablePairs": len(in_bin),
                "SelectedPairs": min(
                    MAX_OVERLAP_QC_PAIRS,
                    len(in_bin)
                )
            }
        )

else:
    above_qc = pd.DataFrame(
        columns=[
            "A_ObjectID",
            "B_ObjectID",
            "A_Percent_Inside_B"
        ]
    )


# ------------------------------------------------------------
# Build BELOW-THRESHOLD bin table
# ------------------------------------------------------------

below_qc_rows = []

best_partial_rows = []

if len(overlap_df) > 0:

    best_partial = (
        overlap_df
        .copy()
        .sort_values(
            [
                "A_ObjectID",
                "A_Percent_Inside_B",
                "B_ObjectID"
            ],
            ascending=[True, False, True]
        )
        .drop_duplicates(
            subset=["A_ObjectID"],
            keep="first"
        )
    )

    # Only the A objects that failed the eligibility criterion.
    best_partial = best_partial[
        best_partial["A_ObjectID"].isin(
            below_threshold_a_ids
        )
    ].copy()

    best_partial["A_Percent_Inside_B"] = pd.to_numeric(
        best_partial["A_Percent_Inside_B"],
        errors="coerce"
    )

    best_partial = best_partial[
        best_partial["A_Percent_Inside_B"].notna()
    ].copy()

    best_partial["OverlapBinIndex"] = np.floor(
        best_partial["A_Percent_Inside_B"]
        / OVERLAP_QC_BIN_WIDTH
    ).astype(int)

    best_partial.loc[
        best_partial["OverlapBinIndex"] >= N_OVERLAP_QC_BINS,
        "OverlapBinIndex"
    ] = N_OVERLAP_QC_BINS - 1

    for bin_index in range(N_OVERLAP_QC_BINS):

        bin_start = (
            bin_index * OVERLAP_QC_BIN_WIDTH
        )

        bin_end = (
            bin_start + OVERLAP_QC_BIN_WIDTH
        )

        # Only show below-threshold bins here.
        # This keeps this PDF explicitly useful for evaluating the
        # current threshold.
        if bin_start >= THRESHOLD * 100:
            continue

        in_bin = best_partial[
            best_partial["OverlapBinIndex"] == bin_index
        ]

        below_qc_rows.append(
            {
                "View": "Below threshold",
                "BinIndex": bin_index,
                "BinStartPercent": bin_start,
                "BinEndPercent": bin_end,
                "AvailablePairs": len(in_bin),
                "SelectedPairs": min(
                    MAX_OVERLAP_QC_PAIRS,
                    len(in_bin)
                )
            }
        )


# ------------------------------------------------------------
# Save bin summary
# ------------------------------------------------------------

overlap_qc_summary_df = pd.DataFrame(
    above_qc_rows + below_qc_rows
)

overlap_qc_summary_df.to_csv(
    overlap_qc_bin_summary_csv,
    index=False
)


# ------------------------------------------------------------
# ABOVE-THRESHOLD PDF
# ------------------------------------------------------------

print()
print("=" * 70)
print(
    f"CREATING {OVERLAP_QC_BIN_WIDTH}% HOECHST-AREA "
    "OVERLAP THRESHOLD QC"
)
print("=" * 70)

print(
    f"\nAbove-threshold PDF:\n"
    f"{overlap_qc_above_pdf}"
)

with PdfPages(overlap_qc_above_pdf) as pdf:

    if len(above_qc) == 0:

        print(
            "No final stable-matched pairs available "
            "for above-threshold overlap QC."
        )

    else:

        for bin_index in range(
            N_OVERLAP_QC_BINS
        ):

            bin_start = (
                bin_index * OVERLAP_QC_BIN_WIDTH
            )

            bin_end = (
                bin_start + OVERLAP_QC_BIN_WIDTH
            )

            # Respect the actual analysis threshold.
            # Example: with a 25% threshold, the above-threshold
            # PDF starts at the 25–30% bin.
            if bin_end <= THRESHOLD * 100:
                continue

            bin_df = above_qc[
                above_qc["OverlapBinIndex"]
                == bin_index
            ].copy()

            if len(bin_df) == 0:
                continue

            page_title = (
                f"FINAL MATCHED PAIRS — "
                f"{bin_start:g}–{bin_end:g}% "
                f"{OBJECT_A_NAME} area inside {OBJECT_B_NAME}"
            )

            make_overlap_qc_page(
                bin_df,
                pdf,
                page_title,
                seed_offset=10_000 + bin_index
            )

            print(
                f"Above threshold "
                f"{bin_start:>3g}–{bin_end:<3g}%: "
                f"{len(bin_df):,} available; "
                f"{min(MAX_OVERLAP_QC_PAIRS, len(bin_df)):,} shown"
            )


# ------------------------------------------------------------
# BELOW-THRESHOLD PDF
# ------------------------------------------------------------

print(
    f"\nBelow-threshold PDF:\n"
    f"{overlap_qc_below_pdf}"
)

with PdfPages(overlap_qc_below_pdf) as pdf:

    if len(best_partial) == 0:

        print(
            "No below-threshold Hoechst objects with "
            "overlapping PV partners are available."
        )

    else:

        for bin_index in range(
            N_OVERLAP_QC_BINS
        ):

            bin_start = (
                bin_index * OVERLAP_QC_BIN_WIDTH
            )

            bin_end = (
                bin_start + OVERLAP_QC_BIN_WIDTH
            )

            # Explicitly only bins below the current threshold.
            if bin_start >= THRESHOLD * 100:
                continue

            bin_df = best_partial[
                best_partial["OverlapBinIndex"]
                == bin_index
            ].copy()

            if len(bin_df) == 0:
                continue

            page_title = (
                f"BEST PARTIAL OVERLAPS — "
                f"{bin_start:g}–{bin_end:g}% "
                f"{OBJECT_A_NAME} area inside {OBJECT_B_NAME}"
            )

            make_overlap_qc_page(
                bin_df,
                pdf,
                page_title,
                seed_offset=20_000 + bin_index
            )

            print(
                f"Below threshold "
                f"{bin_start:>3g}–{bin_end:<3g}%: "
                f"{len(bin_df):,} available; "
                f"{min(MAX_OVERLAP_QC_PAIRS, len(bin_df)):,} shown"
            )


print()
print(
    "5% Hoechst-area overlap QC summary saved to:"
)
print(overlap_qc_bin_summary_csv)


# ============================================================
# FILTERING STAGE COUNTS
# ============================================================

filtering_stage_df.to_csv(
    filtering_stage_csv,
    index=False
)


# ============================================================
# SANKEY HELPER
# ============================================================

def draw_filtering_sankey(
    ax,
    marker_name,
    counts
):
    """
    counts = [all, overlap, eligible, matched]
    """

    all_n, overlap_n, eligible_n, matched_n = [
        int(x) for x in counts
    ]

    stages = [
        ("All\nobjects", all_n),
        ("Any\noverlap", overlap_n),
        ("Eligible\npair", eligible_n),
        ("Final\nmatched", matched_n)
    ]

    x_positions = [0.03, 0.35, 0.67, 0.95]

    max_count = max(all_n, 1)
    top_y = 0.78
    max_height = 0.50
    scale = max_height / max_count

    node_width = 0.065
    node_edge = 1.2

    node_bottoms = []
    node_tops = []

    for x, (label, count) in zip(
        x_positions,
        stages
    ):
        height = max(
            count * scale,
            0.015 if count > 0 else 0.0
        )

        bottom = top_y - height

        rect = Rectangle(
            (x - node_width / 2, bottom),
            node_width,
            height,
            facecolor="white",
            edgecolor="black",
            linewidth=node_edge,
            zorder=5
        )

        ax.add_patch(rect)

        ax.text(
            x,
            bottom - 0.035,
            label,
            ha="center",
            va="top",
            fontsize=10,
            fontweight="bold"
        )

        ax.text(
            x,
            bottom + height + 0.018,
            f"{count:,}",
            ha="center",
            va="bottom",
            fontsize=11,
            fontweight="bold"
        )

        node_bottoms.append(bottom)
        node_tops.append(bottom + height)

    stage_counts = [
        all_n,
        overlap_n,
        eligible_n,
        matched_n
    ]

    for i in range(3):
        current_n = stage_counts[i]
        next_n = stage_counts[i + 1]

        current_height = max(
            current_n * scale,
            0.015 if current_n > 0 else 0.0
        )

        next_height = max(
            next_n * scale,
            0.015 if next_n > 0 else 0.0
        )

        current_center = (
            node_bottoms[i] + node_tops[i]
        ) / 2

        next_center = (
            node_bottoms[i + 1] + node_tops[i + 1]
        ) / 2

        flow_height_current = min(
            current_height,
            next_height
        )

        flow_height_next = flow_height_current

        x0 = x_positions[i] + node_width / 2
        x1 = x_positions[i + 1] - node_width / 2

        y0_top = (
            current_center
            + flow_height_current / 2
        )

        y0_bottom = (
            current_center
            - flow_height_current / 2
        )

        y1_top = (
            next_center
            + flow_height_next / 2
        )

        y1_bottom = (
            next_center
            - flow_height_next / 2
        )

        polygon = Polygon(
            [
                (x0, y0_top),
                (x1, y1_top),
                (x1, y1_bottom),
                (x0, y0_bottom)
            ],
            closed=True,
            facecolor="0.72",
            edgecolor="0.35",
            linewidth=0.7,
            alpha=0.75,
            zorder=2
        )

        ax.add_patch(polygon)

        removed_n = current_n - next_n

        if removed_n > 0:
            removed_height = removed_n * scale

            reject_y0 = (
                node_bottoms[i] - removed_height
            )

            branch_offset = (
                0.06 + 0.075 * i
            )

            x0_reject = (
                x_positions[i] + node_width / 2
            )

            x1_reject = (
                x_positions[i + 1] - node_width / 2
            )

            polygon_reject = Polygon(
                [
                    (
                        x0_reject,
                        node_bottoms[i]
                    ),
                    (
                        x0_reject,
                        reject_y0
                    ),
                    (
                        x1_reject,
                        reject_y0 - branch_offset
                    ),
                    (
                        x1_reject,
                        reject_y0
                    )
                ],
                closed=True,
                facecolor="0.88",
                edgecolor="0.55",
                linewidth=0.6,
                alpha=0.9,
                zorder=1
            )

            ax.add_patch(polygon_reject)

            ax.text(
                (
                    x0_reject + x1_reject
                ) / 2,
                (
                    reject_y0
                    - branch_offset
                    - 0.03
                ),
                f"filtered: {removed_n:,}",
                ha="center",
                va="top",
                fontsize=8.5
            )

    ax.set_xlim(-0.03, 1.02)
    ax.set_ylim(0.02, 0.95)

    ax.set_title(
        marker_name,
        fontsize=15,
        fontweight="bold"
    )

    ax.axis("off")


# ============================================================
# FILTERING SANKEY PDF
# ============================================================

print("\nCreating filtering Sankey diagrams...")

a_sankey = filtering_stage_df[
    filtering_stage_df["Marker"] == OBJECT_A_NAME
]["Count"].tolist()

b_sankey = filtering_stage_df[
    filtering_stage_df["Marker"] == OBJECT_B_NAME
]["Count"].tolist()

with PdfPages(filtering_sankey_pdf) as pdf:

    fig = plt.figure(figsize=(18, 10))

    ax_a_sankey = fig.add_axes(
        [0.02, 0.10, 0.46, 0.76]
    )

    ax_b_sankey = fig.add_axes(
        [0.52, 0.10, 0.46, 0.76]
    )

    draw_filtering_sankey(
        ax_a_sankey,
        OBJECT_A_NAME,
        a_sankey
    )

    draw_filtering_sankey(
        ax_b_sankey,
        OBJECT_B_NAME,
        b_sankey
    )

    fig.suptitle(
        f"{OBJECT_A_NAME}/{OBJECT_B_NAME} Cell Filtering Flow",
        fontsize=19,
        fontweight="bold",
        y=0.96
    )

    fig.text(
        0.5,
        0.035,
        (
            f"Flow stages: all segmented objects → "
            f"any overlap → eligible pair → "
            f"final one-to-one stable match. "
            f"Eligibility is based on {OBJECT_A_NAME} having >= "
            f"{THRESHOLD * 100:.0f}% of its area inside "
            f"{OBJECT_B_NAME}."
        ),
        ha="center",
        va="center",
        fontsize=10.5
    )

    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


# ============================================================
# FULL-RESOLUTION ROI OVERLAY HELPERS
# ============================================================

def make_roi_boundary(
    roi_pixels,
    outline_width=ROI_TIFF_OUTLINE_WIDTH
):
    roi_pixels = np.asarray(
        roi_pixels,
        dtype=int
    )

    if len(roi_pixels) == 0:
        return (
            np.empty(0, dtype=int),
            np.empty(0, dtype=int)
        )

    xs = roi_pixels[:, 0]
    ys = roi_pixels[:, 1]

    xmin = int(xs.min())
    xmax = int(xs.max())
    ymin = int(ys.min())
    ymax = int(ys.max())

    outline_width = max(
        1,
        int(outline_width)
    )

    if outline_width % 2 == 0:
        outline_width += 1

    radius = (outline_width - 1) // 2

    local_width = xmax - xmin + 1
    local_height = ymax - ymin + 1

    padded_width = (
        local_width + 2 * radius
    )
    padded_height = (
        local_height + 2 * radius
    )

    mask = np.zeros(
        (padded_height, padded_width),
        dtype=bool
    )

    local_x = (
        xs - xmin + radius
    )

    local_y = (
        ys - ymin + radius
    )

    mask[local_y, local_x] = True

    boundary = np.zeros_like(mask)

    boundary[1:, :] |= (
        mask[1:, :] & ~mask[:-1, :]
    )
    boundary[:-1, :] |= (
        mask[:-1, :] & ~mask[1:, :]
    )
    boundary[:, 1:] |= (
        mask[:, 1:] & ~mask[:, :-1]
    )
    boundary[:, :-1] |= (
        mask[:, :-1] & ~mask[:, 1:]
    )

    if radius > 0:
        original_boundary = boundary.copy()
        thick_boundary = boundary.copy()

        for dy in range(-radius, radius + 1):
            for dx in range(-radius, radius + 1):
                if dx == 0 and dy == 0:
                    continue

                shifted = np.zeros_like(
                    original_boundary
                )

                y_src_start = max(0, -dy)
                y_src_end = min(
                    padded_height,
                    padded_height - dy
                )

                x_src_start = max(0, -dx)
                x_src_end = min(
                    padded_width,
                    padded_width - dx
                )

                y_dst_start = max(0, dy)
                y_dst_end = min(
                    padded_height,
                    padded_height + dy
                )

                x_dst_start = max(0, dx)
                x_dst_end = min(
                    padded_width,
                    padded_width + dx
                )

                if (
                    y_src_end > y_src_start
                    and x_src_end > x_src_start
                ):
                    shifted[
                        y_dst_start:y_dst_end,
                        x_dst_start:x_dst_end
                    ] = original_boundary[
                        y_src_start:y_src_end,
                        x_src_start:x_src_end
                    ]

                    thick_boundary |= shifted

        boundary = thick_boundary

    by, bx = np.nonzero(boundary)

    return (
        bx + xmin - radius,
        by + ymin - radius
    )


def hsv_roi_color(
    ordinal,
    n_rois
):
    n_rois = max(
        int(n_rois),
        1
    )

    hue = (
        float(ordinal)
        / float(n_rois)
    ) % 1.0

    raw = np.asarray(
        plt.cm.hsv(hue)
    )[:3]

    hsv_rgb = raw.copy()

    hsv_rgb = (
        hsv_rgb * ROI_COLOR_VALUE
    )

    mean_value = np.mean(hsv_rgb)

    hsv_rgb = (
        mean_value
        +
        ROI_COLOR_SATURATION
        * (hsv_rgb - mean_value)
    )

    return np.clip(
        hsv_rgb,
        0,
        1
    )


def build_roi_outline_overlay(
    original_image,
    pixel_groups,
    roi_ids,
    vmin,
    vmax,
    color_lookup=None
):
    image = np.asarray(
        original_image,
        dtype=float
    )

    if vmax <= vmin:
        vmax = vmin + 1.0

    normalized = (
        image - vmin
    ) / (
        vmax - vmin
    )

    normalized = np.clip(
        normalized,
        0,
        1
    )

    gray = (
        normalized * 255
    ).astype(np.uint8)

    rgb_image = np.stack(
        [gray, gray, gray],
        axis=-1
    )

    sorted_ids = sorted(
        int(x) for x in roi_ids
    )

    if color_lookup is None:
        color_lookup = {
            roi_id: hsv_roi_color(
                index,
                len(sorted_ids)
            )
            for index, roi_id in enumerate(sorted_ids)
        }

    for roi_id in sorted_ids:
        if roi_id not in pixel_groups:
            continue

        bx, by = make_roi_boundary(
            pixel_groups[roi_id],
            outline_width=ROI_TIFF_OUTLINE_WIDTH
        )

        if len(bx) == 0:
            continue

        valid = (
            (bx >= 0)
            & (bx < rgb_image.shape[1])
            & (by >= 0)
            & (by < rgb_image.shape[0])
        )

        bx = bx[valid]
        by = by[valid]

        if len(bx) == 0:
            continue

        color = (
            np.asarray(
                color_lookup[roi_id]
            ) * 255.0
        )

        base = (
            rgb_image[
                by,
                bx,
                :
            ].astype(float)
        )

        blended = (
            (1.0 - ROI_OUTLINE_ALPHA) * base
            +
            ROI_OUTLINE_ALPHA * color
        )

        rgb_image[
            by,
            bx,
            :
        ] = np.clip(
            blended,
            0,
            255
        ).astype(np.uint8)

    return rgb_image, color_lookup


def save_color_mapping(
    color_lookup,
    csv_path
):
    rows = []

    for roi_id, color in sorted(
        color_lookup.items()
    ):
        rows.append(
            {
                "ObjectID": int(roi_id),
                "Red": int(round(color[0] * 255)),
                "Green": int(round(color[1] * 255)),
                "Blue": int(round(color[2] * 255))
            }
        )

    pd.DataFrame(
        rows,
        columns=[
            "ObjectID",
            "Red",
            "Green",
            "Blue"
        ]
    ).to_csv(
        csv_path,
        index=False
    )


def add_roi_color_figure(
    ax,
    image,
    title
):
    ax.imshow(
        image,
        interpolation="nearest"
    )

    ax.set_title(
        title,
        fontsize=13,
        fontweight="bold"
    )

    ax.axis("off")


# ============================================================
# FULL-RESOLUTION BEFORE/AFTER TIFF OVERLAYS
# ============================================================

print()
print("=" * 70)
print("CREATING FULL-RESOLUTION ROI OVERLAY TIFFS")
print("=" * 70)

print(
    f"\nTIFF outline width: "
    f"{ROI_TIFF_OUTLINE_WIDTH} pixels"
)

print(
    f"ROI outline alpha: "
    f"{ROI_OUTLINE_ALPHA:.2f}"
)

print(
    f"\nBuilding {OBJECT_A_NAME} before-filtering overlay..."
)

a_before_overlay, a_color_lookup = (
    build_roi_outline_overlay(
        a_original,
        a_pixel_groups,
        all_a_ids,
        a_vmin,
        a_vmax
    )
)

print(
    f"Building {OBJECT_A_NAME} after-filtering overlay..."
)

a_after_overlay, _ = (
    build_roi_outline_overlay(
        a_original,
        a_pixel_groups,
        matched_a_ids,
        a_vmin,
        a_vmax,
        color_lookup=a_color_lookup
    )
)

print(
    f"Building {OBJECT_B_NAME} before-filtering overlay..."
)

b_before_overlay, b_color_lookup = (
    build_roi_outline_overlay(
        b_original,
        b_pixel_groups,
        all_b_ids,
        b_vmin,
        b_vmax
    )
)

print(
    f"Building {OBJECT_B_NAME} after-filtering overlay..."
)

b_after_overlay, _ = (
    build_roi_outline_overlay(
        b_original,
        b_pixel_groups,
        matched_b_ids,
        b_vmin,
        b_vmax,
        color_lookup=b_color_lookup
    )
)


# ============================================================
# WRITE FULL-RESOLUTION RGB TIFFS
# ============================================================

print("\nWriting full-resolution RGB TIFFs...")

tifffile.imwrite(
    roi_filtering_overlay_a_before_tiff,
    a_before_overlay,
    photometric="rgb",
    compression="zlib"
)

tifffile.imwrite(
    roi_filtering_overlay_a_after_tiff,
    a_after_overlay,
    photometric="rgb",
    compression="zlib"
)

tifffile.imwrite(
    roi_filtering_overlay_b_before_tiff,
    b_before_overlay,
    photometric="rgb",
    compression="zlib"
)

tifffile.imwrite(
    roi_filtering_overlay_b_after_tiff,
    b_after_overlay,
    photometric="rgb",
    compression="zlib"
)

save_color_mapping(
    a_color_lookup,
    a_color_mapping_csv
)

save_color_mapping(
    b_color_lookup,
    b_color_mapping_csv
)


# ============================================================
# BEFORE/AFTER PDF OVERVIEW
# ============================================================

print("\nCreating before/after ROI outline PDF...")

with PdfPages(roi_filtering_overlay_pdf) as pdf:

    fig, axes = plt.subplots(
        2,
        2,
        figsize=FILTERING_OVERLAY_MAX_FIGURE_SIZE
    )

    add_roi_color_figure(
        axes[0, 0],
        a_before_overlay,
        (
            f"{OBJECT_A_NAME} — BEFORE FILTERING\n"
            f"All {total_a_objects:,} ROIs"
        )
    )

    add_roi_color_figure(
        axes[0, 1],
        a_after_overlay,
        (
            f"{OBJECT_A_NAME} — AFTER FILTERING\n"
            f"Final matched: {len(matched_a_ids):,}"
        )
    )

    add_roi_color_figure(
        axes[1, 0],
        b_before_overlay,
        (
            f"{OBJECT_B_NAME} — BEFORE FILTERING\n"
            f"All {total_b_objects:,} ROIs"
        )
    )

    add_roi_color_figure(
        axes[1, 1],
        b_after_overlay,
        (
            f"{OBJECT_B_NAME} — AFTER FILTERING\n"
            f"Final matched: {len(matched_b_ids):,}"
        )
    )

    fig.suptitle(
        (
            f"{OBJECT_A_NAME}/{OBJECT_B_NAME} "
            "ROI Filtering — Full-Image QC"
        ),
        fontsize=19,
        fontweight="bold",
        y=0.995
    )

    fig.text(
        0.5,
        0.012,
        (
            "Each ROI is outlined in a distinct color. "
            "The fluorescence image underneath is unchanged; "
            "only the displayed ROI set changes after filtering. "
            "For detailed zooming, use the full-resolution TIFFs."
        ),
        ha="center",
        fontsize=10.5
    )

    plt.tight_layout(
        rect=[0, 0.03, 1, 0.97]
    )

    pdf.savefig(
        fig,
        bbox_inches="tight"
    )

    plt.close(fig)


# ============================================================
# KEPT OBJECT LABEL TIFFS / RGB VISUALIZATION TIFF
# ============================================================

height, width = a_original.shape

a_image = np.zeros(
    (height, width),
    dtype=np.uint32
)

matched_a_pixels = object_a[
    object_a["ObjectID"].isin(matched_a_ids)
]

for object_id, group in matched_a_pixels.groupby(
    "ObjectID"
):
    xs = group["X"].to_numpy()
    ys = group["Y"].to_numpy()

    a_image[ys, xs] = int(object_id)

tifffile.imwrite(
    object_a_tiff,
    a_image
)

b_image = np.zeros(
    (height, width),
    dtype=np.uint32
)

matched_b_pixels = object_b[
    object_b["ObjectID"].isin(matched_b_ids)
]

for object_id, group in matched_b_pixels.groupby(
    "ObjectID"
):
    xs = group["X"].to_numpy()
    ys = group["Y"].to_numpy()

    b_image[ys, xs] = int(object_id)

tifffile.imwrite(
    object_b_tiff,
    b_image
)

two_channel = np.stack(
    [a_image, b_image],
    axis=0
)

tifffile.imwrite(
    combined_tiff,
    two_channel
)

rgb = np.zeros(
    (height, width, 3),
    dtype=np.uint8
)

a_mask = a_image > 0
b_mask = b_image > 0

# Preserve original script:
# Object A = red, Object B = green.
rgb[a_mask, 0] = 255
rgb[b_mask, 1] = 255

# Object B outlines = blue
b_outline = np.zeros(
    (height, width),
    dtype=bool
)

for object_id in matched_b_ids:
    object_mask = (
        b_image == int(object_id)
    )

    boundary = np.zeros(
        (height, width),
        dtype=bool
    )

    boundary[1:, :] |= (
        object_mask[1:, :]
        &
        ~object_mask[:-1, :]
    )

    boundary[:-1, :] |= (
        object_mask[:-1, :]
        &
        ~object_mask[1:, :]
    )

    boundary[:, 1:] |= (
        object_mask[:, 1:]
        &
        ~object_mask[:, :-1]
    )

    boundary[:, :-1] |= (
        object_mask[:, :-1]
        &
        ~object_mask[:, 1:]
    )

    b_outline |= boundary

# Object A outlines = red
a_outline = np.zeros(
    (height, width),
    dtype=bool
)

for object_id in matched_a_ids:
    object_mask = (
        a_image == int(object_id)
    )

    boundary = np.zeros(
        (height, width),
        dtype=bool
    )

    boundary[1:, :] |= (
        object_mask[1:, :]
        &
        ~object_mask[:-1, :]
    )

    boundary[:-1, :] |= (
        object_mask[:-1, :]
        &
        ~object_mask[1:, :]
    )

    boundary[:, 1:] |= (
        object_mask[:, 1:]
        &
        ~object_mask[:, :-1]
    )

    boundary[:, :-1] |= (
        object_mask[:, :-1]
        &
        ~object_mask[:, 1:]
    )

    a_outline |= boundary

rgb[b_outline, 0] = 0
rgb[b_outline, 1] = 0
rgb[b_outline, 2] = 255

rgb[a_outline, 0] = 255
rgb[a_outline, 1] = 0
rgb[a_outline, 2] = 0

tifffile.imwrite(
    rgb_visualization_tiff,
    rgb
)



# ============================================================
# FINAL CLASS-SPECIFIC R/G/Y + METRIC / INTENSITY QC
# ============================================================
#
# This additional visualization-only module divides the final objects into:
#
#   1. POSITIVE:
#      Object-A / Object-B pairs retained by the final one-to-one matching.
#
#   2. A-ONLY:
#      Object-A ROIs that are NOT present in the final matching.
#      Their best overlapping Object-B ROI (when one exists) is shown only
#      as context; it is not reclassified as a positive pair.
#
#   3. B-ONLY:
#      Object-B ROIs that are NOT present in the final matching.
#      Their best overlapping Object-A ROI (when one exists) is shown only
#      as context.
#
# IMPORTANT:
# "A-only" and "B-only" mean "not retained in the FINAL stable match".
# They do not necessarily mean literally zero spatial overlap. An object can
# be unmatched because it has no overlap, only sub-threshold overlap, or an
# eligible partner that was lost during the one-to-one matching competition.
#
# R/G/Y mask TIFF convention used below:
#      red    = selected/context Object A mask
#      green  = selected/context Object B mask
#      yellow = spatial overlap of those two masks
#
# The class-specific 5%-bin PDFs add two rows to the image-style QC:
#      row 1 = Object B grayscale
#      row 2 = Object A grayscale
#      row 3 = red/green/yellow fluorescence composite
#      row 4 = normalized A/B intensity distributions in the FOCAL ROI
#      row 5 = geometry + raw-intensity metrics
#
# This section does NOT change the colocalization analysis.
# ============================================================

print()
print("=" * 70)
print("CREATING UNMATCHED-OBJECT R/G/Y + METRIC/INTENSITY QC")
print("=" * 70)

# ------------------------------------------------------------
# USER SETTINGS FOR THE NEW CLASS-SPECIFIC QC
# ------------------------------------------------------------

CLASS_QC_BIN_WIDTH = OVERLAP_QC_BIN_WIDTH
CLASS_QC_MAX_PER_BIN = MAX_OVERLAP_QC_PAIRS
CLASS_QC_RANDOM_SEED = OVERLAP_QC_RANDOM_SEED + 100_000
CLASS_QC_CROP_PADDING = OVERLAP_QC_CROP_PADDING
CLASS_QC_FIGSIZE = (24, 17)
CLASS_QC_HIST_BINS = 30

# ------------------------------------------------------------
# Final class ID sets
# ------------------------------------------------------------

a_only_ids = set(all_a_ids) - set(matched_a_ids)
b_only_ids = set(all_b_ids) - set(matched_b_ids)

eligible_a_ids = (
    set(eligible_df["A_ObjectID"].astype(int))
    if len(eligible_df) and "A_ObjectID" in eligible_df.columns
    else set()
)

eligible_b_ids = (
    set(eligible_df["B_ObjectID"].astype(int))
    if len(eligible_df) and "B_ObjectID" in eligible_df.columns
    else set()
)

# ------------------------------------------------------------
# New output paths
# ------------------------------------------------------------
# Keep the added unmatched-object QC in its own clean subfolder so
# stale outputs from older versions are not mixed with these results.
UNMATCHED_QC_DIR = OUTPUT_DIR / "unmatched_only_qc"
UNMATCHED_QC_DIR.mkdir(parents=True, exist_ok=True)

a_only_strict_rgy_tiff = UNMATCHED_QC_DIR / (
    f"{output_prefix}_{OBJECT_A_NAME}_ONLY_mask_RGY.tif"
)

b_only_strict_rgy_tiff = UNMATCHED_QC_DIR / (
    f"{output_prefix}_{OBJECT_B_NAME}_ONLY_mask_RGY.tif"
)

a_only_rgy_tiff = UNMATCHED_QC_DIR / (
    f"{output_prefix}_{OBJECT_A_NAME}_ONLY_with_best_{OBJECT_B_NAME}_context_RGY.tif"
)

b_only_rgy_tiff = UNMATCHED_QC_DIR / (
    f"{output_prefix}_{OBJECT_B_NAME}_ONLY_with_best_{OBJECT_A_NAME}_context_RGY.tif"
)

a_only_metrics_csv = UNMATCHED_QC_DIR / (
    f"{output_prefix}_{OBJECT_A_NAME}_ONLY_all_metrics.csv"
)

b_only_metrics_csv = UNMATCHED_QC_DIR / (
    f"{output_prefix}_{OBJECT_B_NAME}_ONLY_all_metrics.csv"
)

unmatched_class_metrics_csv = UNMATCHED_QC_DIR / (
    f"{output_prefix}_unmatched_classes_metrics.csv"
)

class_counts_csv = UNMATCHED_QC_DIR / (
    f"{output_prefix}_final_class_counts.csv"
)

# One-row population-level summary containing marker-relative assignment
# percentages and the final A+/B+, A+/B-, A-/B+ class composition.
summary_metrics_csv = UNMATCHED_QC_DIR / (
    f"{output_prefix}_final_summary_metrics.csv"
)

# Per-nucleus classification and four-class all-cell summary.
nuclear_cell_classification_csv = UNMATCHED_QC_DIR / (
    f"{output_prefix}_all_nuclear_cells_classification.csv"
)

all_cell_class_counts_csv = UNMATCHED_QC_DIR / (
    f"{output_prefix}_all_nuclear_cells_class_counts.csv"
)

class_bin_summary_csv = UNMATCHED_QC_DIR / (
    f"{output_prefix}_final_unmatched_class_5pct_bin_summary.csv"
)

a_only_metric_intensity_pdf = UNMATCHED_QC_DIR / (
    f"{output_prefix}_{OBJECT_A_NAME}_ONLY_metrics_intensity_5pct.pdf"
)

b_only_metric_intensity_pdf = UNMATCHED_QC_DIR / (
    f"{output_prefix}_{OBJECT_B_NAME}_ONLY_metrics_intensity_5pct.pdf"
)

# ------------------------------------------------------------
# Pixel-index caches
# ------------------------------------------------------------

def _flat_pixel_group(pixel_group):
    """Return sorted unique flat indices for an X/Y pixel group."""
    pixels = np.asarray(pixel_group, dtype=np.int64)

    if len(pixels) == 0:
        return np.empty(0, dtype=np.int64)

    flat = (
        pixels[:, 1] * int(image_width)
        + pixels[:, 0]
    )

    return np.unique(flat)


a_flat_groups = {
    int(object_id): _flat_pixel_group(pixels)
    for object_id, pixels in a_pixel_groups.items()
}

b_flat_groups = {
    int(object_id): _flat_pixel_group(pixels)
    for object_id, pixels in b_pixel_groups.items()
}

_intersection_count_cache = {}


def _intersection_count(a_id, b_id):
    """Number of pixels shared by one A ROI and one B ROI."""
    if a_id is None or b_id is None:
        return 0

    key = (int(a_id), int(b_id))

    if key in _intersection_count_cache:
        return _intersection_count_cache[key]

    a_flat = a_flat_groups.get(int(a_id))
    b_flat = b_flat_groups.get(int(b_id))

    if a_flat is None or b_flat is None:
        count = 0
    else:
        count = int(
            np.intersect1d(
                a_flat,
                b_flat,
                assume_unique=True
            ).size
        )

    _intersection_count_cache[key] = count
    return count


def _a_percent_inside_b(a_id, b_id):
    a_flat = a_flat_groups.get(int(a_id)) if a_id is not None else None

    if a_flat is None or len(a_flat) == 0 or b_id is None:
        return 0.0

    return (
        100.0
        * _intersection_count(a_id, b_id)
        / len(a_flat)
    )


def _b_percent_inside_a(a_id, b_id):
    b_flat = b_flat_groups.get(int(b_id)) if b_id is not None else None

    if b_flat is None or len(b_flat) == 0 or a_id is None:
        return 0.0

    return (
        100.0
        * _intersection_count(a_id, b_id)
        / len(b_flat)
    )


# ------------------------------------------------------------
# Best contextual partner lookups for unmatched objects
# ------------------------------------------------------------

best_b_for_a = {}
best_a_for_b = {}

if len(overlap_df) > 0:

    overlap_pairs_for_context = (
        overlap_df[["A_ObjectID", "B_ObjectID"]]
        .dropna()
        .drop_duplicates()
        .copy()
    )

    overlap_pairs_for_context["A_ObjectID"] = (
        overlap_pairs_for_context["A_ObjectID"].astype(int)
    )

    overlap_pairs_for_context["B_ObjectID"] = (
        overlap_pairs_for_context["B_ObjectID"].astype(int)
    )

    # Best B for each A = largest percentage of A inside B.
    for a_id, rows in overlap_pairs_for_context.groupby("A_ObjectID"):
        candidates = []

        for b_id in rows["B_ObjectID"].astype(int):
            candidates.append(
                (
                    _a_percent_inside_b(int(a_id), int(b_id)),
                    int(b_id)
                )
            )

        if candidates:
            # Highest overlap first; lowest ID breaks ties reproducibly.
            candidates.sort(key=lambda x: (-x[0], x[1]))
            best_b_for_a[int(a_id)] = candidates[0][1]

    # Best A for each B = largest percentage of B inside A.
    for b_id, rows in overlap_pairs_for_context.groupby("B_ObjectID"):
        candidates = []

        for a_id in rows["A_ObjectID"].astype(int):
            candidates.append(
                (
                    _b_percent_inside_a(int(a_id), int(b_id)),
                    int(a_id)
                )
            )

        if candidates:
            candidates.sort(key=lambda x: (-x[0], x[1]))
            best_a_for_b[int(b_id)] = candidates[0][1]


# ------------------------------------------------------------
# Geometry + intensity metrics
# ------------------------------------------------------------

PIXEL_AREA_UM2 = float(PIXEL_SIZE_X_UM) * float(PIXEL_SIZE_Y_UM)


def _intensity_stats(image, pixels, prefix):
    """Raw intensity summary for an image sampled at X/Y pixels."""
    result = {
        f"{prefix}_Mean": np.nan,
        f"{prefix}_Median": np.nan,
        f"{prefix}_Std": np.nan,
        f"{prefix}_Min": np.nan,
        f"{prefix}_P05": np.nan,
        f"{prefix}_P95": np.nan,
        f"{prefix}_Max": np.nan,
    }

    if pixels is None:
        return result

    pixels = np.asarray(pixels, dtype=int)

    if len(pixels) == 0:
        return result

    xs = pixels[:, 0]
    ys = pixels[:, 1]

    valid = (
        (xs >= 0)
        & (xs < image.shape[1])
        & (ys >= 0)
        & (ys < image.shape[0])
    )

    xs = xs[valid]
    ys = ys[valid]

    if len(xs) == 0:
        return result

    values = np.asarray(image[ys, xs], dtype=float)
    values = values[np.isfinite(values)]

    if len(values) == 0:
        return result

    result.update(
        {
            f"{prefix}_Mean": float(np.mean(values)),
            f"{prefix}_Median": float(np.median(values)),
            f"{prefix}_Std": float(np.std(values)),
            f"{prefix}_Min": float(np.min(values)),
            f"{prefix}_P05": float(np.percentile(values, 5)),
            f"{prefix}_P95": float(np.percentile(values, 95)),
            f"{prefix}_Max": float(np.max(values)),
        }
    )

    return result


def _intersection_pixels(a_id, b_id):
    """Return X/Y coordinates of the actual A/B mask intersection."""
    if a_id is None or b_id is None:
        return np.empty((0, 2), dtype=int)

    a_flat = a_flat_groups.get(int(a_id))
    b_flat = b_flat_groups.get(int(b_id))

    if a_flat is None or b_flat is None:
        return np.empty((0, 2), dtype=int)

    overlap_flat = np.intersect1d(
        a_flat,
        b_flat,
        assume_unique=True
    )

    if len(overlap_flat) == 0:
        return np.empty((0, 2), dtype=int)

    ys = overlap_flat // int(image_width)
    xs = overlap_flat % int(image_width)

    return np.column_stack([xs, ys]).astype(int)


def _pair_all_metrics(a_id=None, b_id=None):
    """Compute geometry and raw-intensity metrics for one A/B pair/context."""

    a_id_int = int(a_id) if a_id is not None and not pd.isna(a_id) else None
    b_id_int = int(b_id) if b_id is not None and not pd.isna(b_id) else None

    a_pixels = (
        a_pixel_groups.get(a_id_int)
        if a_id_int is not None
        else None
    )

    b_pixels = (
        b_pixel_groups.get(b_id_int)
        if b_id_int is not None
        else None
    )

    a_n = int(len(a_pixels)) if a_pixels is not None else 0
    b_n = int(len(b_pixels)) if b_pixels is not None else 0

    intersection_n = (
        _intersection_count(a_id_int, b_id_int)
        if a_id_int is not None and b_id_int is not None
        else 0
    )

    union_n = a_n + b_n - intersection_n

    a_in_b = (
        100.0 * intersection_n / a_n
        if a_n > 0
        else np.nan
    )

    b_in_a = (
        100.0 * intersection_n / b_n
        if b_n > 0
        else np.nan
    )

    iou_percent = (
        100.0 * intersection_n / union_n
        if union_n > 0
        else np.nan
    )

    dice_percent = (
        100.0 * (2.0 * intersection_n) / (a_n + b_n)
        if (a_n + b_n) > 0
        else np.nan
    )

    centroid_distance_pixels = np.nan
    centroid_distance_um = np.nan
    centroid_dx_pixels = np.nan
    centroid_dy_pixels = np.nan
    centroid_dx_um = np.nan
    centroid_dy_um = np.nan

    if a_n > 0 and b_n > 0:
        a_center_x = float(np.mean(a_pixels[:, 0]))
        a_center_y = float(np.mean(a_pixels[:, 1]))
        b_center_x = float(np.mean(b_pixels[:, 0]))
        b_center_y = float(np.mean(b_pixels[:, 1]))

        centroid_dx_pixels = b_center_x - a_center_x
        centroid_dy_pixels = b_center_y - a_center_y

        centroid_distance_pixels = float(
            np.hypot(
                centroid_dx_pixels,
                centroid_dy_pixels
            )
        )

        centroid_dx_um = (
            centroid_dx_pixels * float(PIXEL_SIZE_X_UM)
        )

        centroid_dy_um = (
            centroid_dy_pixels * float(PIXEL_SIZE_Y_UM)
        )

        centroid_distance_um = float(
            np.hypot(
                centroid_dx_um,
                centroid_dy_um
            )
        )

    intersection_pixels = _intersection_pixels(
        a_id_int,
        b_id_int
    )

    metrics = {
        "A_ObjectID": a_id_int if a_id_int is not None else np.nan,
        "B_ObjectID": b_id_int if b_id_int is not None else np.nan,
        "A_Pixel_Count": a_n if a_n > 0 else np.nan,
        "B_Pixel_Count": b_n if b_n > 0 else np.nan,
        "Intersection_Pixels": intersection_n,
        "Union_Pixels": union_n if union_n > 0 else np.nan,
        "A_Area_um2": a_n * PIXEL_AREA_UM2 if a_n > 0 else np.nan,
        "B_Area_um2": b_n * PIXEL_AREA_UM2 if b_n > 0 else np.nan,
        "Intersection_Area_um2": intersection_n * PIXEL_AREA_UM2,
        "Union_Area_um2": union_n * PIXEL_AREA_UM2 if union_n > 0 else np.nan,
        "A_Percent_Inside_B": a_in_b,
        "B_Percent_Inside_A": b_in_a,
        "IoU_Percent": iou_percent,
        "Dice_Percent": dice_percent,
        "Centroid_DX_Pixels": centroid_dx_pixels,
        "Centroid_DY_Pixels": centroid_dy_pixels,
        "Centroid_Distance_Pixels": centroid_distance_pixels,
        "Centroid_DX_um": centroid_dx_um,
        "Centroid_DY_um": centroid_dy_um,
        "Centroid_Distance_um": centroid_distance_um,
    }

    # A-channel and B-channel intensity inside the A ROI.
    metrics.update(
        _intensity_stats(
            a_original,
            a_pixels,
            f"{OBJECT_A_NAME}_Channel_In_{OBJECT_A_NAME}_ROI"
        )
    )

    metrics.update(
        _intensity_stats(
            b_original,
            a_pixels,
            f"{OBJECT_B_NAME}_Channel_In_{OBJECT_A_NAME}_ROI"
        )
    )

    # A-channel and B-channel intensity inside the B ROI.
    metrics.update(
        _intensity_stats(
            a_original,
            b_pixels,
            f"{OBJECT_A_NAME}_Channel_In_{OBJECT_B_NAME}_ROI"
        )
    )

    metrics.update(
        _intensity_stats(
            b_original,
            b_pixels,
            f"{OBJECT_B_NAME}_Channel_In_{OBJECT_B_NAME}_ROI"
        )
    )

    # Channel intensity inside the true mask intersection.
    metrics.update(
        _intensity_stats(
            a_original,
            intersection_pixels,
            f"{OBJECT_A_NAME}_Channel_In_Intersection"
        )
    )

    metrics.update(
        _intensity_stats(
            b_original,
            intersection_pixels,
            f"{OBJECT_B_NAME}_Channel_In_Intersection"
        )
    )

    return metrics


# ------------------------------------------------------------
# Why an unmatched object is unmatched
# ------------------------------------------------------------

def _a_only_reason(a_id):
    a_id = int(a_id)

    if a_id not in best_b_for_a:
        return "No spatial overlap with any Object B ROI"

    if a_id in eligible_a_ids:
        return "Eligible pair existed, but A was not retained in final one-to-one match"

    return "Best Object B overlap remained below the A-inside-B threshold"


def _b_only_reason(b_id):
    b_id = int(b_id)

    if b_id not in best_a_for_b:
        return "No spatial overlap with any Object A ROI"

    if b_id in eligible_b_ids:
        return "Eligible pair existed, but B was not retained in final one-to-one match"

    return "Only non-eligible/sub-threshold Object A partners overlapped B"


# ------------------------------------------------------------
# Build the two unmatched-class metric tables
# ------------------------------------------------------------

a_only_rows = []

for a_id in sorted(a_only_ids):
    b_id = best_b_for_a.get(int(a_id))

    row = {
        "Class": f"{OBJECT_A_NAME} only / not final matched",
        "Reason": _a_only_reason(a_id),
        "Focal_Object": OBJECT_A_NAME,
        "Context_Partner_Type": (
            f"Best overlapping {OBJECT_B_NAME} context"
            if b_id is not None
            else f"No overlapping {OBJECT_B_NAME} ROI"
        ),
    }

    row.update(_pair_all_metrics(a_id, b_id))
    row["Focal_Overlap_Percent"] = (
        row["A_Percent_Inside_B"]
        if b_id is not None and np.isfinite(row["A_Percent_Inside_B"])
        else 0.0
    )

    a_only_rows.append(row)


b_only_rows = []

for b_id in sorted(b_only_ids):
    a_id = best_a_for_b.get(int(b_id))

    row = {
        "Class": f"{OBJECT_B_NAME} only / not final matched",
        "Reason": _b_only_reason(b_id),
        "Focal_Object": OBJECT_B_NAME,
        "Context_Partner_Type": (
            f"Best overlapping {OBJECT_A_NAME} context"
            if a_id is not None
            else f"No overlapping {OBJECT_A_NAME} ROI"
        ),
    }

    row.update(_pair_all_metrics(a_id, b_id))
    row["Focal_Overlap_Percent"] = (
        row["B_Percent_Inside_A"]
        if a_id is not None and np.isfinite(row["B_Percent_Inside_A"])
        else 0.0
    )

    b_only_rows.append(row)


a_only_metrics_df = pd.DataFrame(a_only_rows)
b_only_metrics_df = pd.DataFrame(b_only_rows)

unmatched_class_metrics_df = pd.concat(
    [
        a_only_metrics_df,
        b_only_metrics_df,
    ],
    ignore_index=True,
    sort=False
)

a_only_metrics_df.to_csv(
    a_only_metrics_csv,
    index=False
)

b_only_metrics_df.to_csv(
    b_only_metrics_csv,
    index=False
)

unmatched_class_metrics_df.to_csv(
    unmatched_class_metrics_csv,
    index=False
)

# ------------------------------------------------------------
# Population-level class / assignment summary
# ------------------------------------------------------------
# Because the final matching is one-to-one, every retained match is one
# inferred double-positive cell (Object A+ / Object B+). The remaining
# Object-A and Object-B ROIs form the two single-positive classes.

n_a_only = len(a_only_ids)
n_b_only = len(b_only_ids)
n_double_positive = len(matching_df)
n_final_classified_cells = (
    n_double_positive + n_a_only + n_b_only
)

def _safe_percent(numerator, denominator):
    if denominator == 0:
        return np.nan
    return 100.0 * float(numerator) / float(denominator)


def _unique_positive_area_pixels(pixel_groups):
    """Return the union area, in pixels, occupied by all ROIs in one marker.

    Using a union mask prevents accidental double-counting if two ROI pixel
    lists ever overlap. Coordinates outside the image are ignored.
    """
    positive_mask = np.zeros(
        (image_height, image_width),
        dtype=bool,
    )

    for pixels in pixel_groups.values():
        if pixels is None or len(pixels) == 0:
            continue

        pixels = np.asarray(pixels, dtype=int)
        xs = pixels[:, 0]
        ys = pixels[:, 1]

        valid = (
            (xs >= 0)
            & (xs < image_width)
            & (ys >= 0)
            & (ys < image_height)
        )

        positive_mask[ys[valid], xs[valid]] = True

    return int(np.count_nonzero(positive_mask))


def _load_organoid_area_pixels(csv_path, area_column):
    """Load total CellProfiler organoid AreaShape_Area in image pixels."""
    csv_path = Path(csv_path)

    if not csv_path.exists():
        raise FileNotFoundError(
            "Organoid measurement CSV not found:\n"
            f"{csv_path}"
        )

    organoid_df = pd.read_csv(csv_path)

    if area_column not in organoid_df.columns:
        area_like_columns = [
            c for c in organoid_df.columns
            if "area" in str(c).lower()
        ]
        raise ValueError(
            f"Organoid area column '{area_column}' was not found in:\n"
            f"{csv_path}\n"
            f"Area-like columns found: {area_like_columns}"
        )

    areas = pd.to_numeric(
        organoid_df[area_column],
        errors="coerce",
    )
    areas = areas[np.isfinite(areas) & (areas > 0)]

    if len(areas) == 0:
        raise ValueError(
            f"No positive finite values were found in '{area_column}' of:\n"
            f"{csv_path}"
        )

    return float(areas.sum()), int(len(areas))


# ------------------------------------------------------------
# Marker-positive area / total organoid area summary
# ------------------------------------------------------------
# The numerator is the UNION of all current marker ROI pixels for that
# marker. The denominator is the total CellProfiler organoid AreaShape_Area.
# These calculations are independent of A/B matching and therefore include
# both matched and unmatched positive marker ROIs.
object_a_positive_area_pixels = np.nan
object_b_positive_area_pixels = np.nan
object_a_positive_area_um2 = np.nan
object_b_positive_area_um2 = np.nan
organoid_area_pixels = np.nan
organoid_area_um2 = np.nan
object_a_positive_area_percent_of_organoid = np.nan
object_b_positive_area_percent_of_organoid = np.nan
organoid_roi_count = np.nan

if ENABLE_ORGANOID_AREA_SUMMARY:
    organoid_area_pixels, organoid_roi_count = _load_organoid_area_pixels(
        ORGANOID_MEASUREMENTS_CSV,
        ORGANOID_AREA_COLUMN,
    )

    object_a_positive_area_pixels = _unique_positive_area_pixels(
        a_pixel_groups
    )
    object_b_positive_area_pixels = _unique_positive_area_pixels(
        b_pixel_groups
    )

    organoid_area_um2 = organoid_area_pixels * PIXEL_AREA_UM2
    object_a_positive_area_um2 = (
        object_a_positive_area_pixels * PIXEL_AREA_UM2
    )
    object_b_positive_area_um2 = (
        object_b_positive_area_pixels * PIXEL_AREA_UM2
    )

    object_a_positive_area_percent_of_organoid = _safe_percent(
        object_a_positive_area_pixels,
        organoid_area_pixels,
    )
    object_b_positive_area_percent_of_organoid = _safe_percent(
        object_b_positive_area_pixels,
        organoid_area_pixels,
    )

    print()
    print("=" * 70)
    print("MARKER-POSITIVE AREA / ORGANOID AREA")
    print("=" * 70)
    print(
        f"Organoid area: {organoid_area_pixels:,.0f} px | "
        f"{organoid_area_um2:,.2f} um^2 "
        f"({int(organoid_roi_count)} organoid ROI row(s))"
    )
    print(
        f"{OBJECT_A_NAME}+ area: "
        f"{object_a_positive_area_pixels:,.0f} px | "
        f"{object_a_positive_area_um2:,.2f} um^2 | "
        f"{object_a_positive_area_percent_of_organoid:.3f}% of organoid"
    )
    print(
        f"{OBJECT_B_NAME}+ area: "
        f"{object_b_positive_area_pixels:,.0f} px | "
        f"{object_b_positive_area_um2:,.2f} um^2 | "
        f"{object_b_positive_area_percent_of_organoid:.3f}% of organoid"
    )

# Final three-class composition. This file now contains the double-positive
# class as well as both single-positive classes.
class_counts_df = pd.DataFrame(
    [
        {
            "Class": f"{OBJECT_A_NAME}+/{OBJECT_B_NAME}+",
            "Count": n_double_positive,
            "Percent_of_Final_Classified_Cells": _safe_percent(
                n_double_positive, n_final_classified_cells
            ),
        },
        {
            "Class": f"{OBJECT_A_NAME}+/{OBJECT_B_NAME}-",
            "Count": n_a_only,
            "Percent_of_Final_Classified_Cells": _safe_percent(
                n_a_only, n_final_classified_cells
            ),
        },
        {
            "Class": f"{OBJECT_A_NAME}-/{OBJECT_B_NAME}+",
            "Count": n_b_only,
            "Percent_of_Final_Classified_Cells": _safe_percent(
                n_b_only, n_final_classified_cells
            ),
        },
    ]
)

class_counts_df.to_csv(
    class_counts_csv,
    index=False
)

# One-row headline summary. The first four percentages use each marker's own
# ROI count as the denominator. The three marker-positive composition
# percentages use A+/B+ + A+/B- + A-/B+ as the denominator.
summary_metrics_record = {
    "Object_A_Name": OBJECT_A_NAME,
    "Object_B_Name": OBJECT_B_NAME,
    "Object_A_Total_ROIs": total_a_objects,
    "Object_B_Total_ROIs": total_b_objects,
    "Double_Positive_Count": n_double_positive,
    "Object_A_Only_Count": n_a_only,
    "Object_B_Only_Count": n_b_only,
    "Final_Classified_Cell_Count": n_final_classified_cells,
    "Object_A_Assigned_to_Object_B_Percent": _safe_percent(
        n_double_positive, total_a_objects
    ),
    "Object_A_Not_Assigned_to_Object_B_Percent": _safe_percent(
        n_a_only, total_a_objects
    ),
    "Object_B_Assigned_to_Object_A_Percent": _safe_percent(
        n_double_positive, total_b_objects
    ),
    "Object_B_Not_Assigned_to_Object_A_Percent": _safe_percent(
        n_b_only, total_b_objects
    ),
    "Double_Positive_Label": f"{OBJECT_A_NAME}+/{OBJECT_B_NAME}+",
    "Double_Positive_Percent_of_Final_Classified_Cells": _safe_percent(
        n_double_positive, n_final_classified_cells
    ),
    "Object_A_Only_Label": f"{OBJECT_A_NAME}+/{OBJECT_B_NAME}-",
    "Object_A_Only_Percent_of_Final_Classified_Cells": _safe_percent(
        n_a_only, n_final_classified_cells
    ),
    "Object_B_Only_Label": f"{OBJECT_A_NAME}-/{OBJECT_B_NAME}+",
    "Object_B_Only_Percent_of_Final_Classified_Cells": _safe_percent(
        n_b_only, n_final_classified_cells
    ),
    # Area-based marker coverage of the organoid.
    "Organoid_Area_Source_CSV": (
        str(ORGANOID_MEASUREMENTS_CSV)
        if ENABLE_ORGANOID_AREA_SUMMARY
        else None
    ),
    "Organoid_Area_Column": (
        ORGANOID_AREA_COLUMN
        if ENABLE_ORGANOID_AREA_SUMMARY
        else None
    ),
    "Organoid_ROI_Count": organoid_roi_count,
    "Organoid_Area_Pixels": organoid_area_pixels,
    "Organoid_Area_um2": organoid_area_um2,
    "Object_A_Positive_Area_Pixels": object_a_positive_area_pixels,
    "Object_A_Positive_Area_um2": object_a_positive_area_um2,
    "Object_A_Positive_Area_Percent_of_Organoid": (
        object_a_positive_area_percent_of_organoid
    ),
    "Object_B_Positive_Area_Pixels": object_b_positive_area_pixels,
    "Object_B_Positive_Area_um2": object_b_positive_area_um2,
    "Object_B_Positive_Area_Percent_of_Organoid": (
        object_b_positive_area_percent_of_organoid
    ),
}


# ============================================================
# TRUE ALL-CELL SUMMARY USING PRE-EXISTING NUCLEAR ASSIGNMENTS
# ============================================================
# This section DOES NOT perform any new marker-to-nucleus spatial matching.
# Instead it reuses the first-stage marker-vs-nuclear results already stored
# by 01_run_colocalization.py.
#
# Why this matters:
#   - Pairwise CellID is run-local and is assigned sequentially.
#   - The canonical cross-run identity is the nuclear ObjectID from the
#     original marker-vs-nuclear stable match.
#
# For each current marker object, this section attempts to recover its
# canonical nuclear ObjectID from the corresponding source HDF5. It supports
# either of the two normal pipeline situations:
#   1. current ObjectID is still the original marker ROI/ObjectID; or
#   2. the current input table carries the first-stage run-local CellID.
#
# Once Object A and Object B are both mapped to canonical nuclear ObjectIDs,
# every nucleus is classified as exactly one of:
#     A+/B+, A+/B-, A-/B+, A-/B-
# and those four percentages use ALL nuclei as the denominator.


def _normalize_name(value):
    return str(value).strip().lower().replace(" ", "")


def _load_source_nuclear_assignment_h5(results_file, expected_marker_name):
    """Load one marker-vs-nuclear HDF5 and expose canonical nucleus mapping."""
    results_file = Path(results_file)
    if not results_file.exists():
        raise FileNotFoundError(
            f"Source marker-vs-nuclear analysis file not found:\n{results_file}"
        )

    with h5py.File(results_file, "r") as h5:
        source_metadata = json.loads(
            h5["metadata/json"][()].decode("utf-8")
        )
        source_mapping = dataframe_from_h5(h5, "cell_id_mapping")
        source_object_a = dataframe_from_h5(h5, "object_a_pixels")
        source_object_b = dataframe_from_h5(h5, "object_b_pixels")

    source_a_name = source_metadata["objects"]["object_a_name"]
    source_b_name = source_metadata["objects"]["object_b_name"]

    expected_norm = _normalize_name(expected_marker_name)
    a_norm = _normalize_name(source_a_name)
    b_norm = _normalize_name(source_b_name)

    if a_norm == expected_norm and b_norm != expected_norm:
        marker_side = "A"
        marker_id_col = "A_ObjectID"
        nuclear_id_col = "B_ObjectID"
        nuclear_table = source_object_b
        nuclear_name = source_b_name
    elif b_norm == expected_norm and a_norm != expected_norm:
        marker_side = "B"
        marker_id_col = "B_ObjectID"
        nuclear_id_col = "A_ObjectID"
        nuclear_table = source_object_a
        nuclear_name = source_a_name
    else:
        raise ValueError(
            "Could not uniquely identify the marker side in source HDF5.\n"
            f"Expected marker: {expected_marker_name}\n"
            f"Source objects: A={source_a_name}, B={source_b_name}\n"
            f"File: {results_file}"
        )

    required_cols = {"CellID", marker_id_col, nuclear_id_col}
    missing = required_cols - set(source_mapping.columns)
    if missing:
        raise ValueError(
            f"Source cell_id_mapping is missing {sorted(missing)} in {results_file}"
        )

    mapping = source_mapping[["CellID", marker_id_col, nuclear_id_col]].copy()
    mapping = mapping.rename(
        columns={
            "CellID": "Source_CellID",
            marker_id_col: "Source_Marker_ObjectID",
            nuclear_id_col: "Nuclear_ObjectID",
        }
    )

    for col in ["Source_CellID", "Source_Marker_ObjectID", "Nuclear_ObjectID"]:
        mapping[col] = pd.to_numeric(mapping[col], errors="raise").astype(int)

    nuclear_ids = set(
        pd.to_numeric(nuclear_table["ObjectID"], errors="raise").astype(int)
    )

    if not nuclear_ids:
        raise ValueError(
            f"No nuclear ObjectIDs were found in source HDF5: {results_file}"
        )

    return {
        "results_file": results_file,
        "metadata": source_metadata,
        "marker_side": marker_side,
        "marker_name": expected_marker_name,
        "nuclear_name": nuclear_name,
        "mapping": mapping,
        "nuclear_ids": nuclear_ids,
    }


def _unique_current_object_to_original_column(current_table, column_name):
    """Return current ObjectID -> selected preserved input-ID column."""
    if column_name not in current_table.columns:
        return None

    pairs = current_table[["ObjectID", column_name]].dropna().drop_duplicates()
    pairs["ObjectID"] = pd.to_numeric(
        pairs["ObjectID"], errors="raise"
    ).astype(int)
    pairs[column_name] = pd.to_numeric(
        pairs[column_name], errors="raise"
    ).astype(int)

    counts = pairs.groupby("ObjectID")[column_name].nunique()
    bad = counts[counts != 1]
    if len(bad):
        raise ValueError(
            f"Current table has ObjectIDs linked to multiple {column_name} values: "
            f"{bad.index.astype(int).tolist()[:10]}"
        )

    return dict(zip(pairs["ObjectID"], pairs[column_name]))


def _map_current_objects_to_canonical_nuclei(
    current_table,
    current_object_ids,
    source_info,
    marker_name,
):
    """
    Map current marker objects back to source canonical nuclear ObjectIDs.

    Preference order:
      1. Direct original marker ObjectID match.
      2. Preserved source CellID column in current object table.
      3. Current ObjectID itself as source CellID (when prior exported data
         contained only CellID/X/Y and the loader promoted CellID to ObjectID).
    """
    source_mapping = source_info["mapping"].copy()

    by_marker_id = dict(
        zip(
            source_mapping["Source_Marker_ObjectID"].astype(int),
            source_mapping["Nuclear_ObjectID"].astype(int),
        )
    )
    by_source_cellid = dict(
        zip(
            source_mapping["Source_CellID"].astype(int),
            source_mapping["Nuclear_ObjectID"].astype(int),
        )
    )

    current_ids = sorted(int(x) for x in current_object_ids)
    current_id_set = set(current_ids)

    # Route 1: current ObjectID is the original marker ROI/ObjectID.
    if current_id_set.issubset(set(by_marker_id)):
        rows = [
            {
                "Current_ObjectID": object_id,
                "Source_Link_Method": "current_ObjectID_to_source_marker_ObjectID",
                "Source_Marker_ObjectID": object_id,
                "Source_CellID": int(
                    source_mapping.loc[
                        source_mapping["Source_Marker_ObjectID"] == object_id,
                        "Source_CellID",
                    ].iloc[0]
                ),
                "Nuclear_ObjectID": int(by_marker_id[object_id]),
            }
            for object_id in current_ids
        ]
        return pd.DataFrame(rows)

    # Route 2: source CellID was preserved as a column in the current table.
    current_to_source_cellid = _unique_current_object_to_original_column(
        current_table, "CellID"
    )
    if current_to_source_cellid is not None:
        source_cellids_needed = {
            int(current_to_source_cellid[x]) for x in current_ids
        }
        if source_cellids_needed.issubset(set(by_source_cellid)):
            source_cellid_to_marker = dict(
                zip(
                    source_mapping["Source_CellID"].astype(int),
                    source_mapping["Source_Marker_ObjectID"].astype(int),
                )
            )
            rows = []
            for object_id in current_ids:
                source_cellid = int(current_to_source_cellid[object_id])
                rows.append(
                    {
                        "Current_ObjectID": object_id,
                        "Source_Link_Method": "preserved_CellID_to_source_CellID",
                        "Source_Marker_ObjectID": int(
                            source_cellid_to_marker[source_cellid]
                        ),
                        "Source_CellID": source_cellid,
                        "Nuclear_ObjectID": int(
                            by_source_cellid[source_cellid]
                        ),
                    }
                )
            return pd.DataFrame(rows)

    # Route 3: current ObjectID itself may be the prior run's CellID.
    if current_id_set.issubset(set(by_source_cellid)):
        source_cellid_to_marker = dict(
            zip(
                source_mapping["Source_CellID"].astype(int),
                source_mapping["Source_Marker_ObjectID"].astype(int),
            )
        )
        rows = [
            {
                "Current_ObjectID": object_id,
                "Source_Link_Method": "current_ObjectID_to_source_CellID",
                "Source_Marker_ObjectID": int(
                    source_cellid_to_marker[object_id]
                ),
                "Source_CellID": object_id,
                "Nuclear_ObjectID": int(by_source_cellid[object_id]),
            }
            for object_id in current_ids
        ]
        return pd.DataFrame(rows)

    missing_marker = sorted(current_id_set - set(by_marker_id))[:10]
    missing_cellid = sorted(current_id_set - set(by_source_cellid))[:10]
    raise ValueError(
        f"Could not link current {marker_name} objects to the stored nuclear "
        "assignment. The current objects matched neither the source marker "
        "ObjectIDs nor source CellIDs.\n"
        f"First unmatched current IDs vs marker IDs: {missing_marker}\n"
        f"First unmatched current IDs vs source CellIDs: {missing_cellid}\n"
        f"Source HDF5: {source_info['results_file']}"
    )


if ENABLE_NUCLEAR_ALL_CELL_SUMMARY:
    if (
        OBJECT_A_NUCLEAR_RESULTS_FILE is None
        or OBJECT_B_NUCLEAR_RESULTS_FILE is None
    ):
        print(
            "\nNuclear all-cell summary requested, but one or both source "
            "marker-vs-nuclear HDF5 paths are None. Skipping all-cell summary.\n"
            "Set OBJECT_A_NUCLEAR_RESULTS_FILE and "
            "OBJECT_B_NUCLEAR_RESULTS_FILE to the original first-stage HDF5s."
        )
        nuclear_summary_completed = False

    else:
        print()
        print("=" * 70)
        print("CREATING ALL-CELL SUMMARY FROM STORED NUCLEAR ASSIGNMENTS")
        print("=" * 70)

        source_a = _load_source_nuclear_assignment_h5(
            OBJECT_A_NUCLEAR_RESULTS_FILE,
            OBJECT_A_NAME,
        )
        source_b = _load_source_nuclear_assignment_h5(
            OBJECT_B_NUCLEAR_RESULTS_FILE,
            OBJECT_B_NAME,
        )

        print(
            f"{OBJECT_A_NAME} nuclear source: "
            f"{source_a['results_file']}"
        )
        print(
            f"{OBJECT_B_NAME} nuclear source: "
            f"{source_b['results_file']}"
        )

        nuclear_ids_a = set(source_a["nuclear_ids"])
        nuclear_ids_b = set(source_b["nuclear_ids"])

        # ------------------------------------------------------------
        # Establish one canonical ALL-NUCLEUS denominator.
        # ------------------------------------------------------------
        # The two first-stage runs may have used differently filtered
        # versions of the SAME nuclear segmentation.  If one nuclear-ID
        # set is a strict subset of the other, their ID namespace is still
        # compatible and the union (equivalently the larger set) is the
        # correct all-cell denominator.
        #
        # If the two sets are non-nested, stop rather than silently mixing
        # potentially different nuclear segmentations.
        if nuclear_ids_a == nuclear_ids_b:
            canonical_nuclear_ids = nuclear_ids_a
            nuclear_denominator_method = "identical_source_nuclear_sets"
        elif nuclear_ids_a.issubset(nuclear_ids_b) or nuclear_ids_b.issubset(nuclear_ids_a):
            if not ALLOW_NESTED_NUCLEAR_OBJECT_SETS:
                raise ValueError(
                    "The marker-vs-nuclear source HDF5 files contain nested "
                    "but non-identical nuclear ObjectID populations, and "
                    "ALLOW_NESTED_NUCLEAR_OBJECT_SETS is False."
                )
            canonical_nuclear_ids = nuclear_ids_a | nuclear_ids_b
            nuclear_denominator_method = "nested_source_nuclear_sets_union"
            smaller_name = (
                OBJECT_A_NAME if len(nuclear_ids_a) < len(nuclear_ids_b)
                else OBJECT_B_NAME
            )
            larger_name = (
                OBJECT_B_NAME if len(nuclear_ids_a) < len(nuclear_ids_b)
                else OBJECT_A_NAME
            )
            print(
                "NOTE: source nuclear populations are nested rather than "
                "identical; using the union/larger canonical nuclear set "
                f"as the all-cell denominator ({len(canonical_nuclear_ids):,} nuclei)."
            )
            print(
                f"  {smaller_name} source nuclear IDs are a subset of "
                f"{larger_name} source nuclear IDs."
            )
        else:
            only_a = sorted(nuclear_ids_a - nuclear_ids_b)
            only_b = sorted(nuclear_ids_b - nuclear_ids_a)
            raise ValueError(
                "The two marker-vs-nuclear source HDF5 files contain "
                "non-nested nuclear ObjectID populations. This may indicate "
                "different nuclear segmentations or incompatible ID namespaces.\n"
                f"{OBJECT_A_NAME} source nuclear count: {len(nuclear_ids_a):,}\n"
                f"{OBJECT_B_NAME} source nuclear count: {len(nuclear_ids_b):,}\n"
                f"Only in {OBJECT_A_NAME} source (first 10): {only_a[:10]}\n"
                f"Only in {OBJECT_B_NAME} source (first 10): {only_b[:10]}"
            )

        # ------------------------------------------------------------
        # IMPORTANT: use the ORIGINAL stored marker->nucleus assignments
        # directly.  Do NOT require every object in the current LHX6/PV
        # analysis to map back to a nucleus.
        # ------------------------------------------------------------
        # The current marker-vs-marker run can contain ROIs that were not
        # retained in the earlier marker-vs-nuclear gate.  Those objects
        # correctly have no stored canonical nuclear assignment and should
        # not cause the all-cell summary to fail.
        #
        # Positivity for the all-cell table therefore means:
        #     A+ = nucleus appears in the final A<->nuclear source mapping
        #     B+ = nucleus appears in the final B<->nuclear source mapping
        # ------------------------------------------------------------
        a_source_mapping = source_a["mapping"].copy()
        b_source_mapping = source_b["mapping"].copy()

        if a_source_mapping["Nuclear_ObjectID"].duplicated().any():
            raise ValueError(
                f"The stored {OBJECT_A_NAME}-vs-nuclear mapping contains "
                "multiple marker assignments to the same canonical nucleus."
            )
        if b_source_mapping["Nuclear_ObjectID"].duplicated().any():
            raise ValueError(
                f"The stored {OBJECT_B_NAME}-vs-nuclear mapping contains "
                "multiple marker assignments to the same canonical nucleus."
            )

        a_positive_nuclear_ids = set(
            a_source_mapping["Nuclear_ObjectID"].astype(int)
        )
        b_positive_nuclear_ids = set(
            b_source_mapping["Nuclear_ObjectID"].astype(int)
        )

        if not a_positive_nuclear_ids.issubset(canonical_nuclear_ids):
            raise ValueError(
                f"Stored {OBJECT_A_NAME} positive nuclear IDs extend outside "
                "the canonical nuclear population."
            )
        if not b_positive_nuclear_ids.issubset(canonical_nuclear_ids):
            raise ValueError(
                f"Stored {OBJECT_B_NAME} positive nuclear IDs extend outside "
                "the canonical nuclear population."
            )

        a_by_nucleus = {
            int(row["Nuclear_ObjectID"]): row
            for _, row in a_source_mapping.iterrows()
        }
        b_by_nucleus = {
            int(row["Nuclear_ObjectID"]): row
            for _, row in b_source_mapping.iterrows()
        }

        nuclear_rows = []
        for nuclear_id in sorted(int(x) for x in canonical_nuclear_ids):
            a_row = a_by_nucleus.get(nuclear_id)
            b_row = b_by_nucleus.get(nuclear_id)

            a_positive = a_row is not None
            b_positive = b_row is not None
            class_label = (
                f"{OBJECT_A_NAME}{'+' if a_positive else '-'}"
                f"/{OBJECT_B_NAME}{'+' if b_positive else '-'}"
            )

            nuclear_rows.append(
                {
                    "Nuclear_ObjectID": nuclear_id,
                    "Nuclear_Marker": source_a["nuclear_name"],
                    "Nuclear_Denominator_Method": nuclear_denominator_method,
                    "Object_A_Name": OBJECT_A_NAME,
                    "Object_A_Positive": bool(a_positive),
                    "Object_A_Current_ObjectID": np.nan,
                    "Object_A_Source_Marker_ObjectID": (
                        int(a_row["Source_Marker_ObjectID"])
                        if a_row is not None else np.nan
                    ),
                    "Object_A_Source_CellID": (
                        int(a_row["Source_CellID"])
                        if a_row is not None else np.nan
                    ),
                    "Object_A_Source_Link_Method": (
                        "direct_first_stage_marker_to_nucleus_assignment"
                        if a_row is not None else ""
                    ),
                    "Object_B_Name": OBJECT_B_NAME,
                    "Object_B_Positive": bool(b_positive),
                    "Object_B_Current_ObjectID": np.nan,
                    "Object_B_Source_Marker_ObjectID": (
                        int(b_row["Source_Marker_ObjectID"])
                        if b_row is not None else np.nan
                    ),
                    "Object_B_Source_CellID": (
                        int(b_row["Source_CellID"])
                        if b_row is not None else np.nan
                    ),
                    "Object_B_Source_Link_Method": (
                        "direct_first_stage_marker_to_nucleus_assignment"
                        if b_row is not None else ""
                    ),
                    "Class": class_label,
                }
            )

        nuclear_cell_df = pd.DataFrame(nuclear_rows)
        nuclear_cell_df.to_csv(
            nuclear_cell_classification_csv,
            index=False,
        )

        all_cell_labels = [
            f"{OBJECT_A_NAME}+/{OBJECT_B_NAME}+",
            f"{OBJECT_A_NAME}+/{OBJECT_B_NAME}-",
            f"{OBJECT_A_NAME}-/{OBJECT_B_NAME}+",
            f"{OBJECT_A_NAME}-/{OBJECT_B_NAME}-",
        ]

        all_cell_counts = {
            label: int((nuclear_cell_df["Class"] == label).sum())
            for label in all_cell_labels
        }
        n_all_nuclear_cells = int(len(nuclear_cell_df))

        all_cell_class_counts_df = pd.DataFrame(
            [
                {
                    "Class": label,
                    "Count": all_cell_counts[label],
                    "Percent_of_All_Nuclear_Cells": _safe_percent(
                        all_cell_counts[label], n_all_nuclear_cells
                    ),
                }
                for label in all_cell_labels
            ]
        )
        all_cell_class_counts_df.to_csv(
            all_cell_class_counts_csv,
            index=False,
        )

        summary_metrics_record.update(
            {
                "All_Cell_Assignment_Source": (
                    "direct_preexisting_marker_vs_nuclear_HDF5_assignments"
                ),
                "All_Cell_Object_A_Nuclear_Source_H5": str(
                    source_a["results_file"]
                ),
                "All_Cell_Object_B_Nuclear_Source_H5": str(
                    source_b["results_file"]
                ),
                "All_Cell_Nuclear_Marker": source_a["nuclear_name"],
                "All_Cell_Nuclear_Denominator_Method": nuclear_denominator_method,
                "All_Cell_Total_Nuclear_Cells": n_all_nuclear_cells,
                "All_Cell_Double_Positive_Label": (
                    f"{OBJECT_A_NAME}+/{OBJECT_B_NAME}+"
                ),
                "All_Cell_Double_Positive_Count": all_cell_counts[
                    f"{OBJECT_A_NAME}+/{OBJECT_B_NAME}+"
                ],
                "All_Cell_Double_Positive_Percent": _safe_percent(
                    all_cell_counts[f"{OBJECT_A_NAME}+/{OBJECT_B_NAME}+"],
                    n_all_nuclear_cells,
                ),
                "All_Cell_Object_A_Only_Label": (
                    f"{OBJECT_A_NAME}+/{OBJECT_B_NAME}-"
                ),
                "All_Cell_Object_A_Only_Count": all_cell_counts[
                    f"{OBJECT_A_NAME}+/{OBJECT_B_NAME}-"
                ],
                "All_Cell_Object_A_Only_Percent": _safe_percent(
                    all_cell_counts[f"{OBJECT_A_NAME}+/{OBJECT_B_NAME}-"],
                    n_all_nuclear_cells,
                ),
                "All_Cell_Object_B_Only_Label": (
                    f"{OBJECT_A_NAME}-/{OBJECT_B_NAME}+"
                ),
                "All_Cell_Object_B_Only_Count": all_cell_counts[
                    f"{OBJECT_A_NAME}-/{OBJECT_B_NAME}+"
                ],
                "All_Cell_Object_B_Only_Percent": _safe_percent(
                    all_cell_counts[f"{OBJECT_A_NAME}-/{OBJECT_B_NAME}+"],
                    n_all_nuclear_cells,
                ),
                "All_Cell_Double_Negative_Label": (
                    f"{OBJECT_A_NAME}-/{OBJECT_B_NAME}-"
                ),
                "All_Cell_Double_Negative_Count": all_cell_counts[
                    f"{OBJECT_A_NAME}-/{OBJECT_B_NAME}-"
                ],
                "All_Cell_Double_Negative_Percent": _safe_percent(
                    all_cell_counts[f"{OBJECT_A_NAME}-/{OBJECT_B_NAME}-"],
                    n_all_nuclear_cells,
                ),
            }
        )

        nuclear_summary_completed = True

        print(
            f"\nCanonical nuclear cells: {n_all_nuclear_cells:,}"
        )
        for label in all_cell_labels:
            count = all_cell_counts[label]
            pct = _safe_percent(count, n_all_nuclear_cells)
            print(f"  {label}: {count:,} ({pct:.2f}%)")

        print(
            f"Per-nucleus classification CSV: "
            f"{nuclear_cell_classification_csv}"
        )
        print(
            f"All-cell class-count CSV: {all_cell_class_counts_csv}"
        )
else:
    nuclear_summary_completed = False


summary_metrics_df = pd.DataFrame([summary_metrics_record])
summary_metrics_df.to_csv(
    summary_metrics_csv,
    index=False
)

print(
    f"\nFinal marker-positive classes: "
    f"{OBJECT_A_NAME}+/{OBJECT_B_NAME}+={n_double_positive:,}, "
    f"{OBJECT_A_NAME}+/{OBJECT_B_NAME}-={n_a_only:,}, "
    f"{OBJECT_A_NAME}-/{OBJECT_B_NAME}+={n_b_only:,}"
)
print(f"Summary metrics CSV: {summary_metrics_csv}")


# ------------------------------------------------------------
# Full-resolution red / green / yellow MASK TIFFs
# ------------------------------------------------------------

def _mask_from_ids(pixel_groups, roi_ids):
    mask = np.zeros(
        (image_height, image_width),
        dtype=bool
    )

    for roi_id in sorted(int(x) for x in roi_ids):
        pixels = pixel_groups.get(roi_id)

        if pixels is None or len(pixels) == 0:
            continue

        pixels = np.asarray(pixels, dtype=int)
        xs = pixels[:, 0]
        ys = pixels[:, 1]

        valid = (
            (xs >= 0)
            & (xs < image_width)
            & (ys >= 0)
            & (ys < image_height)
        )

        mask[ys[valid], xs[valid]] = True

    return mask


def _make_rgy_mask_tiff(a_ids, b_ids, output_path):
    a_mask_local = _mask_from_ids(
        a_pixel_groups,
        a_ids
    )

    b_mask_local = _mask_from_ids(
        b_pixel_groups,
        b_ids
    )

    rgb_local = np.zeros(
        (image_height, image_width, 3),
        dtype=np.uint8
    )

    # Red Object A + green Object B = yellow true mask overlap.
    rgb_local[a_mask_local, 0] = 255
    rgb_local[b_mask_local, 1] = 255

    tifffile.imwrite(
        output_path,
        rgb_local,
        photometric="rgb",
        compression="zlib"
    )


# Context partner IDs for the unmatched-object RGB views.
# These come from the best overlapping partner stored in each unmatched
# object's metrics row. Objects with no overlapping partner simply do not
# contribute an opposite-channel context ROI.
a_only_context_b_ids = {
    int(x)
    for x in a_only_metrics_df.get(
        "B_ObjectID",
        pd.Series(dtype=float)
    ).dropna()
}

b_only_context_a_ids = {
    int(x)
    for x in b_only_metrics_df.get(
        "A_ObjectID",
        pd.Series(dtype=float)
    ).dropna()
}


# Strict negative-class masks: only the final-unmatched focal objects.
_make_rgy_mask_tiff(
    a_only_ids,
    set(),
    a_only_strict_rgy_tiff
)

_make_rgy_mask_tiff(
    set(),
    b_only_ids,
    b_only_strict_rgy_tiff
)

# Context versions additionally show each unmatched object's best spatially
# overlapping counterpart, so partial yellow overlap can be inspected.
_make_rgy_mask_tiff(
    a_only_ids,
    a_only_context_b_ids,
    a_only_rgy_tiff
)

_make_rgy_mask_tiff(
    b_only_context_a_ids,
    b_only_ids,
    b_only_rgy_tiff
)


# ------------------------------------------------------------
# Helpers for the new 5%-increment class PDFs
# ------------------------------------------------------------

def _display_normalize(image, vmin, vmax):
    if vmax <= vmin:
        vmax = vmin + 1.0

    result = (
        np.asarray(image, dtype=float) - float(vmin)
    ) / (
        float(vmax) - float(vmin)
    )

    return np.clip(result, 0, 1)


def _crop_bounds_for_pair(a_id, b_id, focal_object):
    point_groups = []

    if a_id is not None and int(a_id) in a_pixel_groups:
        point_groups.append(a_pixel_groups[int(a_id)])

    if b_id is not None and int(b_id) in b_pixel_groups:
        point_groups.append(b_pixel_groups[int(b_id)])

    if not point_groups:
        return 0, image_width, 0, image_height

    all_points = np.vstack(point_groups)
    all_x = all_points[:, 0]
    all_y = all_points[:, 1]

    xmin = max(
        0,
        int(np.min(all_x)) - CLASS_QC_CROP_PADDING
    )

    xmax = min(
        image_width,
        int(np.max(all_x)) + CLASS_QC_CROP_PADDING + 1
    )

    ymin = max(
        0,
        int(np.min(all_y)) - CLASS_QC_CROP_PADDING
    )

    ymax = min(
        image_height,
        int(np.max(all_y)) + CLASS_QC_CROP_PADDING + 1
    )

    return xmin, xmax, ymin, ymax


def _selected_outline(pixel_groups, object_id, xmin, ymin, xmax, ymax, shape):
    if object_id is None or int(object_id) not in pixel_groups:
        return np.zeros(shape, dtype=bool)

    return make_roi_outline(
        pixel_groups[int(object_id)],
        xmin,
        ymin,
        xmax,
        ymax
    )


def _focal_intensity_values(focal_object, a_id, b_id):
    if focal_object == OBJECT_A_NAME:
        pixels = (
            a_pixel_groups.get(int(a_id))
            if a_id is not None
            else None
        )
    else:
        pixels = (
            b_pixel_groups.get(int(b_id))
            if b_id is not None
            else None
        )

    if pixels is None or len(pixels) == 0:
        return np.empty(0), np.empty(0)

    pixels = np.asarray(pixels, dtype=int)
    xs = pixels[:, 0]
    ys = pixels[:, 1]

    valid = (
        (xs >= 0)
        & (xs < image_width)
        & (ys >= 0)
        & (ys < image_height)
    )

    xs = xs[valid]
    ys = ys[valid]

    if len(xs) == 0:
        return np.empty(0), np.empty(0)

    a_values = _display_normalize(
        a_original[ys, xs],
        a_vmin,
        a_vmax
    ).ravel()

    b_values = _display_normalize(
        b_original[ys, xs],
        b_vmin,
        b_vmax
    ).ravel()

    return a_values, b_values


def _fmt_metric(value, decimals=2):
    if value is None or pd.isna(value):
        return "NA"

    return f"{float(value):.{decimals}f}"


def _metrics_text(row):
    a_id = (
        int(row["A_ObjectID"])
        if "A_ObjectID" in row and pd.notna(row["A_ObjectID"])
        else None
    )

    b_id = (
        int(row["B_ObjectID"])
        if "B_ObjectID" in row and pd.notna(row["B_ObjectID"])
        else None
    )

    a_a_mean = row.get(
        f"{OBJECT_A_NAME}_Channel_In_{OBJECT_A_NAME}_ROI_Mean",
        np.nan
    )
    b_a_mean = row.get(
        f"{OBJECT_B_NAME}_Channel_In_{OBJECT_A_NAME}_ROI_Mean",
        np.nan
    )
    a_b_mean = row.get(
        f"{OBJECT_A_NAME}_Channel_In_{OBJECT_B_NAME}_ROI_Mean",
        np.nan
    )
    b_b_mean = row.get(
        f"{OBJECT_B_NAME}_Channel_In_{OBJECT_B_NAME}_ROI_Mean",
        np.nan
    )
    a_i_mean = row.get(
        f"{OBJECT_A_NAME}_Channel_In_Intersection_Mean",
        np.nan
    )
    b_i_mean = row.get(
        f"{OBJECT_B_NAME}_Channel_In_Intersection_Mean",
        np.nan
    )

    a_a_med = row.get(
        f"{OBJECT_A_NAME}_Channel_In_{OBJECT_A_NAME}_ROI_Median",
        np.nan
    )
    b_a_med = row.get(
        f"{OBJECT_B_NAME}_Channel_In_{OBJECT_A_NAME}_ROI_Median",
        np.nan
    )
    a_b_med = row.get(
        f"{OBJECT_A_NAME}_Channel_In_{OBJECT_B_NAME}_ROI_Median",
        np.nan
    )
    b_b_med = row.get(
        f"{OBJECT_B_NAME}_Channel_In_{OBJECT_B_NAME}_ROI_Median",
        np.nan
    )

    a_a_p95 = row.get(
        f"{OBJECT_A_NAME}_Channel_In_{OBJECT_A_NAME}_ROI_P95",
        np.nan
    )
    b_a_p95 = row.get(
        f"{OBJECT_B_NAME}_Channel_In_{OBJECT_A_NAME}_ROI_P95",
        np.nan
    )
    a_b_p95 = row.get(
        f"{OBJECT_A_NAME}_Channel_In_{OBJECT_B_NAME}_ROI_P95",
        np.nan
    )
    b_b_p95 = row.get(
        f"{OBJECT_B_NAME}_Channel_In_{OBJECT_B_NAME}_ROI_P95",
        np.nan
    )

    lines = [
        f"A={a_id if a_id is not None else 'NA'} | B={b_id if b_id is not None else 'NA'}",
        f"A area: {_fmt_metric(row.get('A_Pixel_Count'), 0)} px | {_fmt_metric(row.get('A_Area_um2'))} um²",
        f"B area: {_fmt_metric(row.get('B_Pixel_Count'), 0)} px | {_fmt_metric(row.get('B_Area_um2'))} um²",
        f"Overlap: {_fmt_metric(row.get('Intersection_Pixels'), 0)} px | {_fmt_metric(row.get('Intersection_Area_um2'))} um²",
        f"A in B: {_fmt_metric(row.get('A_Percent_Inside_B'))}%",
        f"B in A: {_fmt_metric(row.get('B_Percent_Inside_A'))}%",
        f"IoU: {_fmt_metric(row.get('IoU_Percent'))}% | Dice: {_fmt_metric(row.get('Dice_Percent'))}%",
        f"Centroid dist: {_fmt_metric(row.get('Centroid_Distance_Pixels'))} px | {_fmt_metric(row.get('Centroid_Distance_um'))} um",
        "",
        f"A ch in A ROI mean/med/P95: {_fmt_metric(a_a_mean)}/{_fmt_metric(a_a_med)}/{_fmt_metric(a_a_p95)}",
        f"B ch in A ROI mean/med/P95: {_fmt_metric(b_a_mean)}/{_fmt_metric(b_a_med)}/{_fmt_metric(b_a_p95)}",
        f"A ch in B ROI mean/med/P95: {_fmt_metric(a_b_mean)}/{_fmt_metric(a_b_med)}/{_fmt_metric(a_b_p95)}",
        f"B ch in B ROI mean/med/P95: {_fmt_metric(b_b_mean)}/{_fmt_metric(b_b_med)}/{_fmt_metric(b_b_p95)}",
        f"A/B ch mean in intersection: {_fmt_metric(a_i_mean)} / {_fmt_metric(b_i_mean)}",
        "",
        f"Reason: {row.get('Reason', '')}",
    ]

    return "\n".join(lines)


def _make_class_qc_page(
    page_rows,
    pdf,
    page_title,
    focal_object,
    seed_offset=0
):
    if len(page_rows) == 0:
        return

    n_to_show = min(
        CLASS_QC_MAX_PER_BIN,
        len(page_rows)
    )

    rng = np.random.default_rng(
        CLASS_QC_RANDOM_SEED + int(seed_offset)
    )

    if len(page_rows) > n_to_show:
        selected_indices = rng.choice(
            len(page_rows),
            size=n_to_show,
            replace=False
        )

        shown = (
            page_rows
            .iloc[selected_indices]
            .reset_index(drop=True)
        )
    else:
        shown = page_rows.reset_index(drop=True)

    fig = plt.figure(figsize=CLASS_QC_FIGSIZE)

    gs = fig.add_gridspec(
        5,
        n_to_show,
        height_ratios=[1.0, 1.0, 1.0, 0.85, 1.45],
        hspace=0.40,
        wspace=0.16
    )

    for column_number, (_, row) in enumerate(shown.iterrows()):
        a_id = (
            int(row["A_ObjectID"])
            if pd.notna(row.get("A_ObjectID", np.nan))
            else None
        )

        b_id = (
            int(row["B_ObjectID"])
            if pd.notna(row.get("B_ObjectID", np.nan))
            else None
        )

        xmin, xmax, ymin, ymax = _crop_bounds_for_pair(
            a_id,
            b_id,
            focal_object
        )

        a_crop = a_original[ymin:ymax, xmin:xmax]
        b_crop = b_original[ymin:ymax, xmin:xmax]

        a_outline_local = _selected_outline(
            a_pixel_groups,
            a_id,
            xmin,
            ymin,
            xmax,
            ymax,
            a_crop.shape
        )

        b_outline_local = _selected_outline(
            b_pixel_groups,
            b_id,
            xmin,
            ymin,
            xmax,
            ymax,
            b_crop.shape
        )

        # -----------------------------
        # Row 1: B grayscale
        # -----------------------------
        ax_b = fig.add_subplot(gs[0, column_number])
        ax_b.imshow(
            b_crop,
            cmap="gray",
            vmin=b_vmin,
            vmax=b_vmax
        )

        if np.any(b_outline_local):
            ax_b.contour(
                b_outline_local,
                levels=[0.5],
                linewidths=SELECTED_ROI_LINEWIDTH,
                colors="green"
            )

        ax_b.set_title(
            f"{OBJECT_B_NAME} {b_id if b_id is not None else 'no ROI'}",
            fontsize=9,
            fontweight="bold"
        )
        ax_b.axis("off")

        # -----------------------------
        # Row 2: A grayscale
        # -----------------------------
        ax_a = fig.add_subplot(gs[1, column_number])
        ax_a.imshow(
            a_crop,
            cmap="gray",
            vmin=a_vmin,
            vmax=a_vmax
        )

        if np.any(a_outline_local):
            ax_a.contour(
                a_outline_local,
                levels=[0.5],
                linewidths=SELECTED_ROI_LINEWIDTH,
                colors="red"
            )

        ax_a.set_title(
            f"{OBJECT_A_NAME} {a_id if a_id is not None else 'no ROI'}",
            fontsize=9,
            fontweight="bold"
        )
        ax_a.axis("off")

        # -----------------------------
        # Row 3: fluorescence R/G/Y
        # -----------------------------
        ax_rgb = fig.add_subplot(gs[2, column_number])

        a_norm = _display_normalize(
            a_crop,
            a_vmin,
            a_vmax
        )

        b_norm = _display_normalize(
            b_crop,
            b_vmin,
            b_vmax
        )

        rgb_crop = np.zeros(
            (*a_crop.shape, 3),
            dtype=float
        )

        rgb_crop[..., 0] = a_norm
        rgb_crop[..., 1] = b_norm

        ax_rgb.imshow(rgb_crop)

        if np.any(a_outline_local):
            ax_rgb.contour(
                a_outline_local,
                levels=[0.5],
                linewidths=SELECTED_ROI_LINEWIDTH,
                colors="red"
            )

        if np.any(b_outline_local):
            ax_rgb.contour(
                b_outline_local,
                levels=[0.5],
                linewidths=SELECTED_ROI_LINEWIDTH,
                colors="green"
            )

        focal_overlap = float(
            row.get("Focal_Overlap_Percent", 0.0)
        )

        ax_rgb.set_title(
            f"R={OBJECT_A_NAME}, G={OBJECT_B_NAME}, Y=overlap\n"
            f"focal overlap={focal_overlap:.1f}%",
            fontsize=8.5,
            fontweight="bold"
        )
        ax_rgb.axis("off")

        # -----------------------------
        # Row 4: focal ROI intensity distributions
        # -----------------------------
        ax_hist = fig.add_subplot(gs[3, column_number])

        a_values, b_values = _focal_intensity_values(
            focal_object,
            a_id,
            b_id
        )

        bins = np.linspace(0, 1, CLASS_QC_HIST_BINS + 1)

        if len(a_values) > 0:
            ax_hist.hist(
                a_values,
                bins=bins,
                histtype="step",
                density=True,
                linewidth=1.4,
                color="red",
                label=OBJECT_A_NAME
            )

        if len(b_values) > 0:
            ax_hist.hist(
                b_values,
                bins=bins,
                histtype="step",
                density=True,
                linewidth=1.4,
                color="green",
                label=OBJECT_B_NAME
            )

        ax_hist.set_xlim(0, 1)
        ax_hist.set_xlabel(
            "Display-normalized intensity",
            fontsize=7
        )
        ax_hist.set_ylabel("Density", fontsize=7)
        ax_hist.tick_params(labelsize=6)
        ax_hist.set_title(
            f"Intensity inside focal {focal_object} ROI",
            fontsize=8,
            fontweight="bold"
        )

        if len(a_values) > 0 or len(b_values) > 0:
            ax_hist.legend(fontsize=6, frameon=False)

        # -----------------------------
        # Row 5: metrics
        # -----------------------------
        ax_metrics = fig.add_subplot(gs[4, column_number])
        ax_metrics.axis("off")

        ax_metrics.text(
            0.0,
            1.0,
            _metrics_text(row),
            transform=ax_metrics.transAxes,
            ha="left",
            va="top",
            fontsize=6.7,
            family="monospace"
        )

    fig.suptitle(
        page_title,
        fontsize=17,
        fontweight="bold",
        y=0.995
    )

    fig.text(
        0.5,
        0.006,
        (
            f"Red = {OBJECT_A_NAME}; green = {OBJECT_B_NAME}; "
            "yellow = channel/mask overlap. "
            "Intensity histogram values are normalized using the same "
            "display-percentile limits as the image panels; metric text/CSV "
            "retains raw image intensities."
        ),
        ha="center",
        fontsize=9
    )

    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


def _write_empty_pdf(pdf_path, title, message):
    with PdfPages(pdf_path) as pdf:
        fig = plt.figure(figsize=(11, 8.5))
        ax = fig.add_subplot(111)
        ax.axis("off")
        ax.text(
            0.5,
            0.58,
            title,
            ha="center",
            va="center",
            fontsize=18,
            fontweight="bold"
        )
        ax.text(
            0.5,
            0.45,
            message,
            ha="center",
            va="center",
            fontsize=12
        )
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)


def _write_class_5pct_pdf(
    class_df,
    pdf_path,
    class_label,
    focal_object,
    metric_label,
    seed_base
):
    if len(class_df) == 0:
        _write_empty_pdf(
            pdf_path,
            class_label,
            "No objects/pairs were available in this final class."
        )
        return []

    working = class_df.copy()

    working["Focal_Overlap_Percent"] = pd.to_numeric(
        working["Focal_Overlap_Percent"],
        errors="coerce"
    ).fillna(0.0)

    working["OverlapBinIndex"] = np.floor(
        working["Focal_Overlap_Percent"]
        / CLASS_QC_BIN_WIDTH
    ).astype(int)

    n_bins = int(round(100 / CLASS_QC_BIN_WIDTH))

    working.loc[
        working["OverlapBinIndex"] < 0,
        "OverlapBinIndex"
    ] = 0

    working.loc[
        working["OverlapBinIndex"] >= n_bins,
        "OverlapBinIndex"
    ] = n_bins - 1

    bin_rows = []

    with PdfPages(pdf_path) as pdf:
        for bin_index in range(n_bins):
            bin_start = bin_index * CLASS_QC_BIN_WIDTH
            bin_end = bin_start + CLASS_QC_BIN_WIDTH

            in_bin = working[
                working["OverlapBinIndex"] == bin_index
            ].copy()

            bin_rows.append(
                {
                    "Class": class_label,
                    "FocalObject": focal_object,
                    "OverlapMetric": metric_label,
                    "BinIndex": bin_index,
                    "BinStartPercent": bin_start,
                    "BinEndPercent": bin_end,
                    "AvailableObjectsOrPairs": len(in_bin),
                    "SelectedForPDF": min(
                        CLASS_QC_MAX_PER_BIN,
                        len(in_bin)
                    ),
                }
            )

            if len(in_bin) == 0:
                continue

            page_title = (
                f"{class_label} — {bin_start:g}–{bin_end:g}% {metric_label}"
            )

            _make_class_qc_page(
                in_bin,
                pdf,
                page_title,
                focal_object,
                seed_offset=seed_base + bin_index
            )

    return bin_rows


# ------------------------------------------------------------
# Write the two 5%-increment metric/intensity PDFs
# ------------------------------------------------------------

class_bin_rows = []

class_bin_rows.extend(
    _write_class_5pct_pdf(
        a_only_metrics_df,
        a_only_metric_intensity_pdf,
        f"{OBJECT_A_NAME}-ONLY / NOT FINAL MATCHED",
        OBJECT_A_NAME,
        f"{OBJECT_A_NAME} area inside best overlapping {OBJECT_B_NAME}",
        20_000
    )
)

class_bin_rows.extend(
    _write_class_5pct_pdf(
        b_only_metrics_df,
        b_only_metric_intensity_pdf,
        f"{OBJECT_B_NAME}-ONLY / NOT FINAL MATCHED",
        OBJECT_B_NAME,
        f"{OBJECT_B_NAME} area inside best overlapping {OBJECT_A_NAME}",
        30_000
    )
)

pd.DataFrame(class_bin_rows).to_csv(
    class_bin_summary_csv,
    index=False
)

print("\nNew unmatched-object files created:")
for _path in [
    a_only_strict_rgy_tiff,
    b_only_strict_rgy_tiff,
    a_only_rgy_tiff,
    b_only_rgy_tiff,
    a_only_metrics_csv,
    b_only_metrics_csv,
    unmatched_class_metrics_csv,
    class_counts_csv,
    summary_metrics_csv,
    nuclear_cell_classification_csv,
    all_cell_class_counts_csv,
    class_bin_summary_csv,
    a_only_metric_intensity_pdf,
    b_only_metric_intensity_pdf,
]:
    print(_path)


# ============================================================
# FINAL OUTPUT SUMMARY
# ============================================================

print()
print("=" * 70)
print("VISUALIZATION FILES CREATED")
print("=" * 70)

for path in [
    visualization_pdf,
    below_threshold_pdf,
    filtering_stage_csv,
    filtering_sankey_pdf,
    roi_filtering_overlay_pdf,
    roi_filtering_overlay_a_before_tiff,
    roi_filtering_overlay_a_after_tiff,
    roi_filtering_overlay_b_before_tiff,
    roi_filtering_overlay_b_after_tiff,
    a_color_mapping_csv,
    b_color_mapping_csv,
    object_a_tiff,
    object_b_tiff,
    combined_tiff,
    rgb_visualization_tiff,
    a_only_strict_rgy_tiff,
    b_only_strict_rgy_tiff,
    a_only_rgy_tiff,
    b_only_rgy_tiff,
    a_only_metrics_csv,
    b_only_metrics_csv,
    unmatched_class_metrics_csv,
    class_counts_csv,
    summary_metrics_csv,
    nuclear_cell_classification_csv,
    all_cell_class_counts_csv,
    class_bin_summary_csv,
    a_only_metric_intensity_pdf,
    b_only_metric_intensity_pdf
]:
    print(path)

print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
