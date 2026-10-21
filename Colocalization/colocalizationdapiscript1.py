#!/usr/bin/env python3
"""
01_run_colocalization.py

Analysis-only stage.

This script:
    1. Loads the two ROI pixel CSVs.
    2. Calculates shared pixels and directional overlap percentages.
    3. Applies the overlap threshold.
    4. Builds preference lists.
    5. Runs the same Object-A-proposing stable matching used in the
       original script.
    6. Assigns shared CellIDs.
    7. Calculates summary statistics and optional organoid densities.
    8. Saves EVERYTHING needed by the visualization/QC stage into ONE
       self-contained HDF5 file.

The visualization script does not need the original CSVs or TIFFs after
this script has successfully created the HDF5 results file.

Dependencies:
    pandas
    numpy
    tifffile
    h5py
"""

from pathlib import Path
from datetime import datetime
import io
import json

import h5py
import numpy as np
import pandas as pd
import tifffile


# ============================================================
# USER SETTINGS
# ============================================================

OBJECT_A_NAME = "LHX6"
OBJECT_B_NAME = "PV"

OBJECT_A_CSV = Path(
    "/home/oltij/Desktop/ROI_intensitymetric_filtering/LHX6/LHX6_ROI_pixels_FILTERED.csv"
)

OBJECT_B_CSV = Path(
    "/home/oltij/Desktop/ROI_intensitymetric_filtering/PV/PV_ROI_pixels_FILTERED.csv"
)

OBJECT_A_TIFF = Path(
    "/home/oltij/Desktop/LHX6Aligned/aligned_LHX6_max.tif"
)

OBJECT_B_TIFF = Path(
    "/home/oltij/Desktop/PVAligned/aligned_PV_max.tif"
)

# Set to None to skip organoid-area/density calculations.
ORGANOID_CSV = Path(
    "/home/oltij/Desktop/HoeschtAligned/MyExpt_FilterObjects.csv"
)
# ORGANOID_CSV = None

# ONE file is produced by this script.
RESULTS_FILE = Path(
    "/home/oltij/Desktop/Aligned_LHX6PV_colocalization_results_25_new/"
    "Aligned_LHX6PV_colocalization_analysis.h5"
)

THRESHOLD = 0.25

PIXEL_SIZE_X_UM = 0.312843137
PIXEL_SIZE_Y_UM = 0.312836345

# Optional downstream marker-marker constraint. The default False path is the
# historical unconstrained analysis and intentionally preserves its behavior.
NUCLEAR_IDENTITY_CONSTRAINT = False
OBJECT_A_NUCLEAR_RESULTS_FILE = None
OBJECT_B_NUCLEAR_RESULTS_FILE = None
NUCLEAR_CONSTRAINT_AUDIT_CSV = None


# ============================================================
# HELPERS
# ============================================================

def find_id_column(df: pd.DataFrame, object_name: str) -> str:
    if "ROI" in df.columns:
        return "ROI"
    if "CellID" in df.columns:
        return "CellID"
    raise ValueError(
        f"{object_name} CSV must contain either 'ROI' or 'CellID'. "
        f"Found columns: {list(df.columns)}"
    )


def load_object_csv(csv_path: Path, object_name: str):
    print(f"\nLoading {object_name} pixel file...")
    df = pd.read_csv(csv_path)
    print(f"{object_name} pixel rows: {len(df):,}")

    id_column = find_id_column(df, object_name)

    required = {id_column, "X", "Y"}
    if not required.issubset(df.columns):
        raise ValueError(
            f"{object_name} CSV must contain {required}. "
            f"Found columns: {list(df.columns)}"
        )

    df = df.copy()
    df["ObjectID"] = (
        pd.to_numeric(df[id_column], errors="raise").astype(int)
    )
    df["X"] = pd.to_numeric(df["X"], errors="raise").astype(int)
    df["Y"] = pd.to_numeric(df["Y"], errors="raise").astype(int)

    df = df.drop_duplicates(
        subset=["ObjectID", "X", "Y"]
    ).copy()

    return df, id_column


def dataframe_to_bytes(df: pd.DataFrame) -> bytes:
    """Serialize a DataFrame losslessly enough for this pipeline."""
    buffer = io.StringIO()
    df.to_csv(buffer, index=False)
    return buffer.getvalue().encode("utf-8")


def store_dataframe(h5: h5py.File, name: str, df: pd.DataFrame):
    payload = np.frombuffer(
        dataframe_to_bytes(df),
        dtype=np.uint8
    )
    h5.create_dataset(
        f"tables/{name}",
        data=payload,
        compression="gzip",
        compression_opts=6
    )


def store_image(h5: h5py.File, name: str, image: np.ndarray):
    h5.create_dataset(
        f"images/{name}",
        data=image,
        compression="gzip",
        compression_opts=4,
        shuffle=True
    )


def dataframe_from_bytes(h5: h5py.File, name: str) -> pd.DataFrame:
    payload = h5[f"tables/{name}"][()]
    return pd.read_csv(io.BytesIO(payload.tobytes()))


def require_columns(df: pd.DataFrame, columns: set[str], label: str):
    missing = columns - set(df.columns)
    if missing:
        raise ValueError(f"{label} is missing required columns: {sorted(missing)}")


def normalized_xy_groups(df: pd.DataFrame, id_column: str, label: str):
    require_columns(df, {id_column, "X", "Y"}, label)
    work = df[[id_column, "X", "Y"]].copy()
    for column in (id_column, "X", "Y"):
        work[column] = pd.to_numeric(work[column], errors="raise").astype(int)
    work = work.drop_duplicates().sort_values([id_column, "Y", "X"])
    return {
        int(object_id): tuple(map(tuple, group[["X", "Y"]].to_numpy(dtype=int)))
        for object_id, group in work.groupby(id_column, sort=True)
    }


def unique_optional_int(group: pd.DataFrame, column: str):
    if column not in group.columns:
        return None, None
    values = pd.to_numeric(group[column], errors="coerce").dropna().astype(int).unique()
    if len(values) == 0:
        return None, None
    if len(values) > 1:
        return None, f"multiple_{column}_values"
    return int(values[0]), None


def load_nuclear_assignment_context(
    results_file: Path,
    marker_name: str,
    current_pixels: pd.DataFrame,
    current_id_column: str,
):
    """Resolve current marker ROIs to canonical nuclei using coordinates.

    Stage-local CellIDs and original marker ROI identifiers are supporting
    provenance only. Exact Stage-07 marker pixel coordinates are the
    authoritative identity relationship.
    """
    if results_file is None or not Path(results_file).is_file():
        raise FileNotFoundError(
            f"Missing Stage-07 nuclear assignment HDF5 for {marker_name}: {results_file}"
        )
    results_file = Path(results_file)
    with h5py.File(results_file, "r") as h5:
        metadata = json.loads(h5["metadata/json"][()].decode("utf-8"))
        source_marker = dataframe_from_bytes(h5, "object_a_pixels")
        source_nuclear = dataframe_from_bytes(h5, "object_b_pixels")
        mapping = dataframe_from_bytes(h5, "cell_id_mapping")
        filtered_marker = dataframe_from_bytes(h5, "filtered_object_a_pixels")
        marker_image_shape = tuple(h5["images/object_a_original"].shape)
        nuclear_image_shape = tuple(h5["images/object_b_original"].shape)

    source_a_name = metadata["objects"]["object_a_name"]
    nuclear_name = metadata["objects"]["object_b_name"]
    if source_a_name != marker_name:
        raise ValueError(
            f"{results_file}: expected Stage-07 Object A {marker_name!r}, found {source_a_name!r}"
        )
    require_columns(
        mapping,
        {"CellID", "A_ObjectID", "B_ObjectID"},
        f"{results_file} cell_id_mapping",
    )
    mapping = mapping[["CellID", "A_ObjectID", "B_ObjectID"]].copy()
    for column in mapping.columns:
        mapping[column] = pd.to_numeric(mapping[column], errors="raise").astype(int)
    for column in mapping.columns:
        if mapping[column].duplicated().any():
            raise ValueError(
                f"{results_file}: Stage-07 {column} is not one-to-one; canonical assignment is ambiguous"
            )

    source_groups = normalized_xy_groups(
        source_marker, "ObjectID", f"{results_file} marker pixels"
    )
    nuclear_groups = normalized_xy_groups(
        source_nuclear, "ObjectID", f"{results_file} nuclear pixels"
    )
    filtered_groups = normalized_xy_groups(
        filtered_marker, "CellID", f"{results_file} filtered marker pixels"
    )
    mapping_by_cell = mapping.set_index("CellID")
    mapping_by_marker = mapping.set_index("A_ObjectID")

    for row in mapping.itertuples(index=False):
        marker_id = int(row.A_ObjectID)
        cell_id = int(row.CellID)
        nuclear_id = int(row.B_ObjectID)
        if marker_id not in source_groups or cell_id not in filtered_groups:
            raise ValueError(
                f"{results_file}: mapping row CellID={cell_id}, marker={marker_id} lacks source pixels"
            )
        if source_groups[marker_id] != filtered_groups[cell_id]:
            raise ValueError(
                f"{results_file}: Stage-07 filtered CellID {cell_id} does not reproduce marker ObjectID {marker_id} coordinates"
            )
        if nuclear_id not in nuclear_groups:
            raise ValueError(
                f"{results_file}: canonical nuclear ObjectID {nuclear_id} is absent from nuclear pixels"
            )

    matched_marker_ids = set(mapping["A_ObjectID"].astype(int))
    coordinate_to_marker_ids: dict[tuple, list[int]] = {}
    for marker_id in sorted(matched_marker_ids):
        coordinate_to_marker_ids.setdefault(source_groups[marker_id], []).append(marker_id)

    current_groups = normalized_xy_groups(
        current_pixels, "ObjectID", f"current {marker_name} pixels"
    )
    assignment_rows = []
    for current_id, coordinates in current_groups.items():
        current_group = current_pixels[current_pixels["ObjectID"].astype(int) == current_id]
        input_roi, roi_error = unique_optional_int(current_group, "ROI")
        input_cell, cell_error = unique_optional_int(current_group, "CellID")
        coordinate_candidates = coordinate_to_marker_ids.get(coordinates, [])

        identifier_candidates = set()
        if input_roi is not None and input_roi in mapping_by_marker.index:
            identifier_candidates.add(int(input_roi))
        if input_cell is not None and input_cell in mapping_by_cell.index:
            identifier_candidates.add(int(mapping_by_cell.loc[input_cell, "A_ObjectID"]))
        if current_id_column == "ROI" and current_id in mapping_by_marker.index:
            identifier_candidates.add(int(current_id))
        if current_id_column == "CellID" and current_id in mapping_by_cell.index:
            identifier_candidates.add(int(mapping_by_cell.loc[current_id, "A_ObjectID"]))

        source_marker_id = None
        source_cell_id = None
        canonical_nuclear_id = None
        status = "resolved"
        if roi_error or cell_error:
            status = roi_error or cell_error
        elif len(coordinate_candidates) == 0:
            status = "missing_coordinate_identity"
        elif len(coordinate_candidates) > 1:
            status = "ambiguous_coordinate_identity"
        else:
            source_marker_id = int(coordinate_candidates[0])
            if identifier_candidates and identifier_candidates != {source_marker_id}:
                status = "identifier_coordinate_mismatch"
                source_marker_id = None
            else:
                source_row = mapping_by_marker.loc[source_marker_id]
                source_cell_id = int(source_row["CellID"])
                canonical_nuclear_id = int(source_row["B_ObjectID"])

        assignment_rows.append(
            {
                "Current_ObjectID": int(current_id),
                "Input_ID_Column": current_id_column,
                "Input_ROI": input_roi,
                "Input_CellID": input_cell,
                "Stage07_Marker_ObjectID": source_marker_id,
                "Stage07_Source_CellID": source_cell_id,
                "Canonical_Nuclear_ObjectID": canonical_nuclear_id,
                "Coordinate_Verified": status == "resolved",
                "Mapping_Status": status,
            }
        )

    assignments = pd.DataFrame(
        assignment_rows,
        columns=[
            "Current_ObjectID", "Input_ID_Column", "Input_ROI", "Input_CellID",
            "Stage07_Marker_ObjectID", "Stage07_Source_CellID",
            "Canonical_Nuclear_ObjectID", "Coordinate_Verified", "Mapping_Status",
        ],
    )
    return {
        "results_file": results_file,
        "metadata": metadata,
        "marker_name": marker_name,
        "nuclear_name": nuclear_name,
        "marker_image_shape": marker_image_shape,
        "nuclear_image_shape": nuclear_image_shape,
        "pixel_size_um": metadata.get("pixel_size_um", {}),
        "nuclear_pixels": source_nuclear[["ObjectID", "X", "Y"]].copy(),
        "assignments": assignments,
    }


def validate_nuclear_contexts(context_a: dict, context_b: dict):
    if context_a["nuclear_name"] != context_b["nuclear_name"]:
        raise ValueError(
            "Stage-07 sources use different nuclear channels: "
            f"{context_a['nuclear_name']} vs {context_b['nuclear_name']}"
        )
    if context_a["nuclear_image_shape"] != context_b["nuclear_image_shape"]:
        raise ValueError("Stage-07 sources use different nuclear image dimensions")
    if context_a["marker_image_shape"] != context_a["nuclear_image_shape"] or context_b["marker_image_shape"] != context_b["nuclear_image_shape"]:
        raise ValueError("Stage-07 marker and nuclear images do not share one coordinate canvas")

    for axis in ("x", "y"):
        a_value = float(context_a["pixel_size_um"].get(axis, np.nan))
        b_value = float(context_b["pixel_size_um"].get(axis, np.nan))
        if not np.isfinite(a_value) or not np.isfinite(b_value) or not np.isclose(a_value, b_value):
            raise ValueError(f"Stage-07 sources use incompatible {axis}-pixel calibration")

    def canonical_nuclear_table(context):
        table = context["nuclear_pixels"].copy()
        for column in ("ObjectID", "X", "Y"):
            table[column] = pd.to_numeric(table[column], errors="raise").astype(int)
        return table.drop_duplicates().sort_values(["ObjectID", "Y", "X"]).reset_index(drop=True)

    if not canonical_nuclear_table(context_a).equals(canonical_nuclear_table(context_b)):
        raise ValueError(
            "Stage-07 sources do not contain the same canonical nuclear segmentation/ObjectID coordinate table"
        )


def build_nuclear_constraint_audit(
    overlap_df: pd.DataFrame,
    assignments_a: pd.DataFrame,
    assignments_b: pd.DataFrame,
) -> pd.DataFrame:
    columns = [
        "A_ObjectID", "B_ObjectID", "A_Current_ObjectID", "B_Current_ObjectID",
        "A_Stage07_Marker_ObjectID", "B_Stage07_Marker_ObjectID",
        "A_Stage07_Source_CellID", "B_Stage07_Source_CellID",
        "A_Canonical_Nuclear_ObjectID", "B_Canonical_Nuclear_ObjectID",
        "A_Mapping_Status", "B_Mapping_Status", "A_Coordinate_Verified", "B_Coordinate_Verified",
        "Same_Nucleus", "A_Percent_Inside_B", "B_Percent_Inside_A",
        "Original_Overlap_Threshold_Pass", "Nuclear_Constraint_Pass",
        "Final_Constrained_Eligibility", "Final_Constrained_Stable_Match",
        "Exclusion_Reason",
    ]
    if overlap_df.empty:
        return pd.DataFrame(columns=columns)

    a_lookup = assignments_a.set_index("Current_ObjectID").to_dict("index")
    b_lookup = assignments_b.set_index("Current_ObjectID").to_dict("index")
    rows = []
    for overlap in overlap_df.itertuples(index=False):
        a_id = int(overlap.A_ObjectID)
        b_id = int(overlap.B_ObjectID)
        a_assignment = a_lookup.get(a_id)
        b_assignment = b_lookup.get(b_id)
        a_status = "missing_current_object" if a_assignment is None else a_assignment["Mapping_Status"]
        b_status = "missing_current_object" if b_assignment is None else b_assignment["Mapping_Status"]
        a_resolved = a_status == "resolved"
        b_resolved = b_status == "resolved"
        resolved = a_resolved and b_resolved
        a_nucleus = (
            int(a_assignment["Canonical_Nuclear_ObjectID"])
            if a_resolved else None
        )
        b_nucleus = (
            int(b_assignment["Canonical_Nuclear_ObjectID"])
            if b_resolved else None
        )
        same_nucleus = bool(a_nucleus == b_nucleus) if resolved else pd.NA
        threshold_pass = bool(overlap.Eligible)
        nuclear_pass = bool(resolved and same_nucleus)
        final_eligible = bool(threshold_pass and nuclear_pass)
        if not threshold_pass:
            reason = "below_directional_overlap_threshold"
        elif not resolved:
            reason = f"unresolved_nuclear_assignment:A={a_status};B={b_status}"
        elif not same_nucleus:
            reason = "different_canonical_nuclear_objectids"
        else:
            reason = ""
        rows.append(
            {
                "A_ObjectID": a_id,
                "B_ObjectID": b_id,
                "A_Current_ObjectID": a_id,
                "B_Current_ObjectID": b_id,
                "A_Stage07_Marker_ObjectID": None if a_assignment is None else a_assignment["Stage07_Marker_ObjectID"],
                "B_Stage07_Marker_ObjectID": None if b_assignment is None else b_assignment["Stage07_Marker_ObjectID"],
                "A_Stage07_Source_CellID": None if a_assignment is None else a_assignment["Stage07_Source_CellID"],
                "B_Stage07_Source_CellID": None if b_assignment is None else b_assignment["Stage07_Source_CellID"],
                "A_Canonical_Nuclear_ObjectID": a_nucleus,
                "B_Canonical_Nuclear_ObjectID": b_nucleus,
                "A_Mapping_Status": a_status,
                "B_Mapping_Status": b_status,
                "A_Coordinate_Verified": False if a_assignment is None else bool(a_assignment["Coordinate_Verified"]),
                "B_Coordinate_Verified": False if b_assignment is None else bool(b_assignment["Coordinate_Verified"]),
                "Same_Nucleus": same_nucleus,
                "A_Percent_Inside_B": float(overlap.A_Percent_Inside_B),
                "B_Percent_Inside_A": float(overlap.B_Percent_Inside_A),
                "Original_Overlap_Threshold_Pass": threshold_pass,
                "Nuclear_Constraint_Pass": nuclear_pass,
                "Final_Constrained_Eligibility": final_eligible,
                "Final_Constrained_Stable_Match": False,
                "Exclusion_Reason": reason,
            }
        )
    return pd.DataFrame(rows, columns=columns)


# ============================================================
# LOAD OBJECTS
# ============================================================

object_a, object_a_original_id_column = load_object_csv(
    OBJECT_A_CSV,
    OBJECT_A_NAME
)

object_b, object_b_original_id_column = load_object_csv(
    OBJECT_B_CSV,
    OBJECT_B_NAME
)

all_a_ids = set(object_a["ObjectID"].astype(int))
all_b_ids = set(object_b["ObjectID"].astype(int))

total_a_objects = len(all_a_ids)
total_b_objects = len(all_b_ids)

print(f"\nTotal {OBJECT_A_NAME} objects: {total_a_objects:,}")
print(f"Total {OBJECT_B_NAME} objects: {total_b_objects:,}")

constraint_context_a = None
constraint_context_b = None
constraint_audit_df = pd.DataFrame()

if NUCLEAR_IDENTITY_CONSTRAINT:
    print("\nResolving canonical nuclear identities from Stage-07 analyses...")
    constraint_context_a = load_nuclear_assignment_context(
        OBJECT_A_NUCLEAR_RESULTS_FILE,
        OBJECT_A_NAME,
        object_a,
        object_a_original_id_column,
    )
    constraint_context_b = load_nuclear_assignment_context(
        OBJECT_B_NUCLEAR_RESULTS_FILE,
        OBJECT_B_NAME,
        object_b,
        object_b_original_id_column,
    )
    validate_nuclear_contexts(constraint_context_a, constraint_context_b)
    for label, context in ((OBJECT_A_NAME, constraint_context_a), (OBJECT_B_NAME, constraint_context_b)):
        counts = context["assignments"]["Mapping_Status"].value_counts().to_dict()
        print(f"{label} Stage-07 provenance status: {counts}")


# ============================================================
# OBJECT AREAS
# ============================================================

a_areas = (
    object_a.groupby("ObjectID")
    .size()
    .rename("A_Area")
)

b_areas = (
    object_b.groupby("ObjectID")
    .size()
    .rename("B_Area")
)


# ============================================================
# CALCULATE SHARED PIXELS
# ============================================================

print()
print("=" * 70)
print(f"CALCULATING {OBJECT_A_NAME}/{OBJECT_B_NAME} OVERLAPS")
print("=" * 70)

print("\nFinding shared pixels...")

shared_pixels = pd.merge(
    object_a[["ObjectID", "X", "Y"]],
    object_b[["ObjectID", "X", "Y"]],
    on=["X", "Y"],
    how="inner",
    suffixes=("_A", "_B")
)

print(f"Shared pixels: {len(shared_pixels):,}")

if len(shared_pixels) > 0:
    overlap_counts = (
        shared_pixels
        .groupby(["ObjectID_A", "ObjectID_B"])
        .size()
        .reset_index(name="Shared_Area")
    )

    overlap_counts = (
        overlap_counts
        .merge(
            a_areas,
            left_on="ObjectID_A",
            right_index=True,
            how="left"
        )
        .merge(
            b_areas,
            left_on="ObjectID_B",
            right_index=True,
            how="left"
        )
    )

    overlap_counts["A_Percent_Inside_B"] = (
        overlap_counts["Shared_Area"]
        / overlap_counts["A_Area"]
        * 100
    )

    overlap_counts["B_Percent_Inside_A"] = (
        overlap_counts["Shared_Area"]
        / overlap_counts["B_Area"]
        * 100
    )

    overlap_counts["Eligible"] = (
        overlap_counts["A_Percent_Inside_B"]
        >= THRESHOLD * 100
    )

    overlap_df = (
        overlap_counts
        .rename(
            columns={
                "ObjectID_A": "A_ObjectID",
                "ObjectID_B": "B_ObjectID"
            }
        )[
            [
                "A_ObjectID",
                "B_ObjectID",
                "A_Area",
                "B_Area",
                "Shared_Area",
                "A_Percent_Inside_B",
                "B_Percent_Inside_A",
                "Eligible"
            ]
        ]
    )
else:
    overlap_df = pd.DataFrame(
        columns=[
            "A_ObjectID",
            "B_ObjectID",
            "A_Area",
            "B_Area",
            "Shared_Area",
            "A_Percent_Inside_B",
            "B_Percent_Inside_A",
            "Eligible"
        ]
    )

if len(overlap_df) == 0 and not NUCLEAR_IDENTITY_CONSTRAINT:
    raise RuntimeError(
        f"No overlapping pixels were found between "
        f"{OBJECT_A_NAME} and {OBJECT_B_NAME}."
    )

if NUCLEAR_IDENTITY_CONSTRAINT:
    constraint_audit_df = build_nuclear_constraint_audit(
        overlap_df,
        constraint_context_a["assignments"],
        constraint_context_b["assignments"],
    )
    if len(overlap_df):
        constraint_columns = [
            "A_ObjectID", "B_ObjectID",
            "A_Stage07_Marker_ObjectID", "B_Stage07_Marker_ObjectID",
            "A_Stage07_Source_CellID", "B_Stage07_Source_CellID",
            "A_Canonical_Nuclear_ObjectID", "B_Canonical_Nuclear_ObjectID",
            "A_Mapping_Status", "B_Mapping_Status",
            "A_Coordinate_Verified", "B_Coordinate_Verified", "Same_Nucleus",
            "Original_Overlap_Threshold_Pass", "Nuclear_Constraint_Pass",
            "Final_Constrained_Eligibility", "Exclusion_Reason",
        ]
        overlap_df = overlap_df.drop(columns=["Eligible"]).merge(
            constraint_audit_df[constraint_columns],
            on=["A_ObjectID", "B_ObjectID"],
            how="left",
            validate="one_to_one",
        )
        overlap_df["Eligible"] = overlap_df["Final_Constrained_Eligibility"].astype(bool)

eligible_df = overlap_df[
    overlap_df["Eligible"]
].copy()

print(f"Total overlapping pairs: {len(overlap_df):,}")
print(
    f"Eligible pairs (>= {THRESHOLD * 100:.0f}% of "
    f"{OBJECT_A_NAME} inside {OBJECT_B_NAME}): "
    f"{len(eligible_df):,}"
)

if len(overlap_df):
    print(
        f"\nMaximum {OBJECT_A_NAME} percentage inside "
        f"{OBJECT_B_NAME}: "
        f"{overlap_df['A_Percent_Inside_B'].max():.2f}%"
    )

    print(
        f"Median {OBJECT_A_NAME} percentage inside "
        f"{OBJECT_B_NAME}: "
        f"{overlap_df['A_Percent_Inside_B'].median():.2f}%"
    )
else:
    print("\nNo spatially overlapping marker pairs were found; constrained outputs will be empty.")


# ============================================================
# PREFERENCE LISTS
# ============================================================

a_preferences = {}
b_preferences = {}

if len(eligible_df) > 0:
    print("\nBuilding preference lists...")

    for a_id, group in eligible_df.groupby("A_ObjectID"):
        ranked = group.sort_values(
            by=["A_Percent_Inside_B", "B_ObjectID"],
            ascending=[False, True]
        )
        a_preferences[int(a_id)] = (
            ranked["B_ObjectID"].astype(int).tolist()
        )

    for b_id, group in eligible_df.groupby("B_ObjectID"):
        # PRESERVED FROM THE ORIGINAL SCRIPT:
        # B-side preferences are ranked by A_Percent_Inside_B.
        ranked = group.sort_values(
            by=["A_Percent_Inside_B", "A_ObjectID"],
            ascending=[False, True]
        )
        b_preferences[int(b_id)] = (
            ranked["A_ObjectID"].astype(int).tolist()
        )


# ============================================================
# PREFERENCE MATRIX
# ============================================================

preference_rows = []

if len(eligible_df) > 0:
    for _, row in eligible_df.iterrows():
        a_id = int(row["A_ObjectID"])
        b_id = int(row["B_ObjectID"])

        a_rank = a_preferences[a_id].index(b_id) + 1
        b_rank = b_preferences[b_id].index(a_id) + 1

        preference_rows.append(
            {
                "A_ObjectID": a_id,
                "B_ObjectID": b_id,
                "A_Percent_Inside_B": row["A_Percent_Inside_B"],
                "B_Percent_Inside_A": row["B_Percent_Inside_A"],
                "A_Preference_Rank": a_rank,
                "B_Preference_Rank": b_rank
            }
        )

preference_df = pd.DataFrame(
    preference_rows,
    columns=[
        "A_ObjectID",
        "B_ObjectID",
        "A_Percent_Inside_B",
        "B_Percent_Inside_A",
        "A_Preference_Rank",
        "B_Preference_Rank"
    ]
)


# ============================================================
# STABLE MATCHING
# ============================================================

b_current_match = {}
a_current_match = {}

if len(eligible_df) > 0:
    print()
    print("=" * 70)
    print(
        f"STARTING {OBJECT_A_NAME}-PROPOSING "
        f"STABLE MATCHING"
    )
    print("=" * 70)

    next_preference_index = {
        a_id: 0
        for a_id in a_preferences
    }

    proposal_queue = list(a_preferences.keys())

    overlap_lookup = {
        (
            int(row["A_ObjectID"]),
            int(row["B_ObjectID"])
        ): float(row["A_Percent_Inside_B"])
        for _, row in eligible_df.iterrows()
    }

    while proposal_queue:
        a_id = proposal_queue.pop(0)

        preferences = a_preferences[a_id]
        current_index = next_preference_index[a_id]

        if current_index >= len(preferences):
            continue

        b_id = preferences[current_index]
        next_preference_index[a_id] += 1

        if b_id not in b_current_match:
            b_current_match[b_id] = a_id
            a_current_match[a_id] = b_id
            continue

        current_a = b_current_match[b_id]

        new_score = overlap_lookup[(a_id, b_id)]
        current_score = overlap_lookup[(current_a, b_id)]

        if new_score > current_score:
            del a_current_match[current_a]

            b_current_match[b_id] = a_id
            a_current_match[a_id] = b_id

            if (
                current_a in a_preferences
                and
                next_preference_index[current_a]
                < len(a_preferences[current_a])
            ):
                proposal_queue.append(current_a)
        else:
            if (
                next_preference_index[a_id]
                < len(a_preferences[a_id])
            ):
                proposal_queue.append(a_id)


# ============================================================
# FINAL MATCHING TABLE
# ============================================================

matching_rows = []

if len(a_current_match) > 0:
    for a_id, b_id in sorted(a_current_match.items()):
        key = (a_id, b_id)
        score = overlap_lookup[key]

        pair_row = (
            eligible_df[
                (eligible_df["A_ObjectID"] == a_id)
                &
                (eligible_df["B_ObjectID"] == b_id)
            ].iloc[0]
        )

        a_rank = a_preferences[a_id].index(b_id) + 1
        b_rank = b_preferences[b_id].index(a_id) + 1

        matching_rows.append(
            {
                "A_ObjectID": a_id,
                "B_ObjectID": b_id,
                "A_Area": int(pair_row["A_Area"]),
                "B_Area": int(pair_row["B_Area"]),
                "Shared_Area": int(pair_row["Shared_Area"]),
                "A_Percent_Inside_B": score,
                "B_Percent_Inside_A": float(
                    pair_row["B_Percent_Inside_A"]
                ),
                "A_Preference_Rank": a_rank,
                "B_Preference_Rank": b_rank,
                "Pass_Threshold": True
            }
        )

matching_df = pd.DataFrame(
    matching_rows,
    columns=[
        "A_ObjectID",
        "B_ObjectID",
        "A_Area",
        "B_Area",
        "Shared_Area",
        "A_Percent_Inside_B",
        "B_Percent_Inside_A",
        "A_Preference_Rank",
        "B_Preference_Rank",
        "Pass_Threshold"
    ]
)

matched_a_ids = set(
    matching_df["A_ObjectID"].astype(int)
) if len(matching_df) else set()

matched_b_ids = set(
    matching_df["B_ObjectID"].astype(int)
) if len(matching_df) else set()

if NUCLEAR_IDENTITY_CONSTRAINT:
    matched_keys = set(
        zip(
            matching_df["A_ObjectID"].astype(int),
            matching_df["B_ObjectID"].astype(int),
        )
    ) if len(matching_df) else set()
    if len(constraint_audit_df):
        constraint_audit_df["Final_Constrained_Stable_Match"] = [
            (int(row.A_ObjectID), int(row.B_ObjectID)) in matched_keys
            for row in constraint_audit_df.itertuples(index=False)
        ]
        eligible_not_matched = (
            constraint_audit_df["Final_Constrained_Eligibility"].astype(bool)
            & ~constraint_audit_df["Final_Constrained_Stable_Match"].astype(bool)
        )
        constraint_audit_df.loc[
            eligible_not_matched, "Exclusion_Reason"
        ] = "not_selected_by_constrained_stable_matching"
    if NUCLEAR_CONSTRAINT_AUDIT_CSV is not None:
        audit_path = Path(NUCLEAR_CONSTRAINT_AUDIT_CSV)
        audit_path.parent.mkdir(parents=True, exist_ok=True)
        audit_temp = audit_path.with_suffix(audit_path.suffix + ".tmp")
        constraint_audit_df.to_csv(audit_temp, index=False)
        audit_temp.replace(audit_path)


# ============================================================
# CELL ID ASSIGNMENT
# ============================================================

cell_id_mapping = pd.DataFrame(
    columns=["CellID", "A_ObjectID", "B_ObjectID"]
)

if len(matching_df) > 0:
    print()
    print("=" * 70)
    print("ASSIGNING SHARED CELL IDS")
    print("=" * 70)

    cell_id_mapping = (
        matching_df[["A_ObjectID", "B_ObjectID"]]
        .sort_values(["A_ObjectID", "B_ObjectID"])
        .reset_index(drop=True)
    )

    cell_id_mapping["CellID"] = np.arange(
        1,
        len(cell_id_mapping) + 1
    )

    cell_id_mapping = cell_id_mapping[
        ["CellID", "A_ObjectID", "B_ObjectID"]
    ]


# ============================================================
# FILTERED PIXEL TABLES
# ============================================================

a_id_to_cell = dict(
    zip(
        cell_id_mapping["A_ObjectID"].astype(int),
        cell_id_mapping["CellID"].astype(int)
    )
)

b_id_to_cell = dict(
    zip(
        cell_id_mapping["B_ObjectID"].astype(int),
        cell_id_mapping["CellID"].astype(int)
    )
)

filtered_a = object_a[
    object_a["ObjectID"].isin(a_id_to_cell.keys())
].copy()

if len(filtered_a) > 0:
    filtered_a["CellID"] = (
        filtered_a["ObjectID"].astype(int).map(a_id_to_cell)
    )
    filtered_a = filtered_a.drop(columns=["ObjectID"])
    original_cols = [
        col for col in filtered_a.columns if col != "CellID"
    ]
    filtered_a = filtered_a[
        ["CellID"] + original_cols
    ].sort_values(["CellID", "Y", "X"])
else:
    filtered_a = pd.DataFrame(
        columns=["CellID", "X", "Y"]
    )

filtered_b = object_b[
    object_b["ObjectID"].isin(b_id_to_cell.keys())
].copy()

if len(filtered_b) > 0:
    filtered_b["CellID"] = (
        filtered_b["ObjectID"].astype(int).map(b_id_to_cell)
    )
    filtered_b = filtered_b.drop(columns=["ObjectID"])
    original_cols = [
        col for col in filtered_b.columns if col != "CellID"
    ]
    filtered_b = filtered_b[
        ["CellID"] + original_cols
    ].sort_values(["CellID", "Y", "X"])
else:
    filtered_b = pd.DataFrame(
        columns=["CellID", "X", "Y"]
    )


# ============================================================
# POSITIVE OBJECT A TABLE
# ============================================================

positive_df = matching_df.copy()


# ============================================================
# BELOW-THRESHOLD A OBJECTS + BEST PARTIAL MATCH
# ============================================================

all_eligible_a_ids = set(
    eligible_df["A_ObjectID"].astype(int)
)

below_threshold_a_ids = (
    all_a_ids - all_eligible_a_ids
)

best_b_for_a = {}

if len(overlap_df) > 0:
    best_overlap_rows = (
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

    for _, row in best_overlap_rows.iterrows():
        best_b_for_a[int(row["A_ObjectID"])] = row

unmatched_rows = []

for a_id in sorted(below_threshold_a_ids):
    best_row = best_b_for_a.get(a_id)

    if best_row is not None:
        best_b = int(best_row["B_ObjectID"])
        best_shared_area = int(best_row["Shared_Area"])
        best_a_percent = float(
            best_row["A_Percent_Inside_B"]
        )
        best_b_percent = float(
            best_row["B_Percent_Inside_A"]
        )
    else:
        best_b = np.nan
        best_shared_area = 0
        best_a_percent = 0.0
        best_b_percent = 0.0

    unmatched_rows.append(
        {
            "A_ObjectID": a_id,
            "Best_B_ObjectID": best_b,
            "Best_Shared_Area": best_shared_area,
            "Best_A_Percent_Inside_B": best_a_percent,
            "Best_B_Percent_Inside_A": best_b_percent
        }
    )

unmatched_df = pd.DataFrame(unmatched_rows)


# ============================================================
# SUMMARY / ORGANOID AREA / DENSITIES
# ============================================================

n_matched_a = len(matched_a_ids)
n_matched_b = len(matched_b_ids)

a_proportion = (
    n_matched_a / total_a_objects
    if total_a_objects > 0 else np.nan
)

b_proportion = (
    n_matched_b / total_b_objects
    if total_b_objects > 0 else np.nan
)

organoid_area_pixels = np.nan
organoid_area_mm2 = np.nan

if ORGANOID_CSV is not None:
    print("\nLoading organoid mask CSV...")
    organoid_df = pd.read_csv(ORGANOID_CSV)

    if "AreaShape_Area" not in organoid_df.columns:
        raise ValueError(
            "Organoid CSV does not contain 'AreaShape_Area'."
        )

    organoid_area_pixels = (
        pd.to_numeric(
            organoid_df["AreaShape_Area"],
            errors="coerce"
        )
        .dropna()
        .sum()
    )

    pixel_area_um2 = (
        PIXEL_SIZE_X_UM * PIXEL_SIZE_Y_UM
    )

    pixel_area_mm2 = (
        pixel_area_um2 / 1_000_000
    )

    organoid_area_mm2 = (
        organoid_area_pixels * pixel_area_mm2
    )

if (
    pd.notna(organoid_area_mm2)
    and organoid_area_mm2 > 0
):
    a_density = total_a_objects / organoid_area_mm2
    b_density = total_b_objects / organoid_area_mm2
    matched_a_density = n_matched_a / organoid_area_mm2
    matched_b_density = n_matched_b / organoid_area_mm2
else:
    a_density = np.nan
    b_density = np.nan
    matched_a_density = np.nan
    matched_b_density = np.nan

summary_rows = [
    (f"Total {OBJECT_A_NAME} objects", total_a_objects, "objects"),
    (f"Matched {OBJECT_A_NAME} objects", n_matched_a, "objects"),
    (
        f"Proportion of {OBJECT_A_NAME} objects matched",
        a_proportion,
        "fraction"
    ),
    (f"Total {OBJECT_B_NAME} objects", total_b_objects, "objects"),
    (f"Matched {OBJECT_B_NAME} objects", n_matched_b, "objects"),
    (
        f"Proportion of {OBJECT_B_NAME} objects matched",
        b_proportion,
        "fraction"
    ),
    ("Matched CellIDs", len(cell_id_mapping), "cells"),
    (
        (
            f"{OBJECT_A_NAME} objects without a constrained-eligible partner"
            if NUCLEAR_IDENTITY_CONSTRAINT
            else f"{OBJECT_A_NAME} objects below threshold"
        ),
        len(below_threshold_a_ids),
        "objects"
    ),
    ("Organoid area", organoid_area_pixels, "pixels"),
    ("Organoid area", organoid_area_mm2, "mm^2"),
    (f"{OBJECT_A_NAME} density", a_density, "objects/mm^2"),
    (f"{OBJECT_B_NAME} density", b_density, "objects/mm^2"),
    (
        f"Matched {OBJECT_A_NAME} density",
        matched_a_density,
        "objects/mm^2"
    ),
    (
        f"Matched {OBJECT_B_NAME} density",
        matched_b_density,
        "objects/mm^2"
    )
]

summary_df = pd.DataFrame(
    summary_rows,
    columns=["Metric", "Value", "Units"]
)


# ============================================================
# FILTERING STAGE COUNTS
# ============================================================

a_overlap_ids = set(overlap_df["A_ObjectID"].astype(int))
b_overlap_ids = set(overlap_df["B_ObjectID"].astype(int))

a_eligible_ids = set(eligible_df["A_ObjectID"].astype(int))
b_eligible_ids = set(eligible_df["B_ObjectID"].astype(int))

filtering_stage_rows = []

for marker_name, total_n, overlap_n, eligible_n, matched_n in [
    (
        OBJECT_A_NAME,
        total_a_objects,
        len(a_overlap_ids),
        len(a_eligible_ids),
        len(matched_a_ids)
    ),
    (
        OBJECT_B_NAME,
        total_b_objects,
        len(b_overlap_ids),
        len(b_eligible_ids),
        len(matched_b_ids)
    )
]:
    filtering_stage_rows.extend(
        [
            {
                "Marker": marker_name,
                "Stage": "All segmented objects",
                "Count": int(total_n)
            },
            {
                "Marker": marker_name,
                "Stage": "At least one overlapping partner",
                "Count": int(overlap_n)
            },
            {
                "Marker": marker_name,
                "Stage": (
                    (
                        f"Participating in eligible pair "
                        f"(A >= {THRESHOLD * 100:.0f}% inside B and same canonical nucleus)"
                        if NUCLEAR_IDENTITY_CONSTRAINT
                        else f"Participating in eligible pair "
                        f"(A >= {THRESHOLD * 100:.0f}% inside B)"
                    )
                ),
                "Count": int(eligible_n)
            },
            {
                "Marker": marker_name,
                "Stage": "Final one-to-one matched objects",
                "Count": int(matched_n)
            }
        ]
    )

filtering_stage_df = pd.DataFrame(
    filtering_stage_rows
)


# ============================================================
# LOAD ORIGINAL FLUORESCENCE
# ============================================================

print("\nLoading original fluorescence images...")

a_original = np.squeeze(tifffile.imread(OBJECT_A_TIFF))
b_original = np.squeeze(tifffile.imread(OBJECT_B_TIFF))

if a_original.ndim != 2:
    raise ValueError(
        f"{OBJECT_A_NAME} original TIFF is not 2D after squeezing."
    )

if b_original.ndim != 2:
    raise ValueError(
        f"{OBJECT_B_NAME} original TIFF is not 2D after squeezing."
    )

if a_original.shape != b_original.shape:
    raise ValueError(
        "The two fluorescence images have different dimensions:\n"
        f"{OBJECT_A_NAME}: {a_original.shape}\n"
        f"{OBJECT_B_NAME}: {b_original.shape}"
    )

if NUCLEAR_IDENTITY_CONSTRAINT:
    if tuple(a_original.shape) != tuple(constraint_context_a["marker_image_shape"]):
        raise ValueError(
            f"Current {OBJECT_A_NAME} fluorescence canvas {a_original.shape} does not match "
            f"its Stage-07 source {constraint_context_a['marker_image_shape']}"
        )
    if tuple(b_original.shape) != tuple(constraint_context_b["marker_image_shape"]):
        raise ValueError(
            f"Current {OBJECT_B_NAME} fluorescence canvas {b_original.shape} does not match "
            f"its Stage-07 source {constraint_context_b['marker_image_shape']}"
        )

height, width = a_original.shape


# ============================================================
# WRITE SINGLE SELF-CONTAINED HDF5 ANALYSIS FILE
# ============================================================

RESULTS_FILE.parent.mkdir(
    parents=True,
    exist_ok=True
)

metadata = {
    "run": {
        "date_time": datetime.now().astimezone().isoformat(),
        "script": Path(__file__).name
    },
    "objects": {
        "object_a_name": OBJECT_A_NAME,
        "object_b_name": OBJECT_B_NAME
    },
    "id_columns": {
        "object_a": object_a_original_id_column,
        "object_b": object_b_original_id_column
    },
    "input_files": {
        "object_a_csv": str(OBJECT_A_CSV),
        "object_b_csv": str(OBJECT_B_CSV),
        "object_a_tiff": str(OBJECT_A_TIFF),
        "object_b_tiff": str(OBJECT_B_TIFF),
        "organoid_csv": (
            None if ORGANOID_CSV is None
            else str(ORGANOID_CSV)
        )
    },
    "matching": {
        "threshold": THRESHOLD,
        "threshold_percent": THRESHOLD * 100,
        "matching_direction": f"{OBJECT_A_NAME}-proposing",
        "algorithm": "Stable matching / deferred acceptance"
    },
    "pixel_size_um": {
        "x": PIXEL_SIZE_X_UM,
        "y": PIXEL_SIZE_Y_UM,
        "pixel_area_um2": (
            PIXEL_SIZE_X_UM * PIXEL_SIZE_Y_UM
        )
    },
    "image": {
        "height": int(height),
        "width": int(width),
        "object_a_dtype": str(a_original.dtype),
        "object_b_dtype": str(b_original.dtype)
    }
}

if NUCLEAR_IDENTITY_CONSTRAINT:
    metadata["analysis_variant"] = "nucleus_constrained_marker_marker"
    metadata["nuclear_identity_constraint"] = {
        "enabled": True,
        "applied_before_preference_lists": True,
        "canonical_identity": "Stage-07 nuclear Object B ObjectID",
        "unresolved_assignment_policy": "exclude with explicit audit reason",
        "object_a_stage07_results_file": str(OBJECT_A_NUCLEAR_RESULTS_FILE),
        "object_b_stage07_results_file": str(OBJECT_B_NUCLEAR_RESULTS_FILE),
        "audit_csv": (
            None if NUCLEAR_CONSTRAINT_AUDIT_CSV is None
            else str(NUCLEAR_CONSTRAINT_AUDIT_CSV)
        ),
        "same_nuclear_segmentation_verified_by": (
            "exact ObjectID/X/Y nuclear pixel table, image dimensions, and pixel calibration"
        ),
        "marker_identity_verified_by": (
            "exact current-to-Stage-07 marker ROI pixel coordinates; ROI and CellID are provenance checks only"
        ),
    }
    metadata["identifier_definitions"] = {
        "current_marker_objectid": "ObjectID used by this marker-marker run",
        "stage07_marker_objectid": "original marker ROI ObjectID in the marker-vs-nuclear run",
        "stage07_source_cellid": "run-local CellID assigned by that marker-vs-nuclear run",
        "canonical_nuclear_objectid": "Stage-07 nuclear Object B ObjectID; cross-marker identity anchor",
        "output_cellid": "run-local CellID assigned only to this constrained marker-marker result",
    }

print("\nWriting single analysis results file:")
print(RESULTS_FILE)

with h5py.File(RESULTS_FILE, "w") as h5:
    h5.attrs["format"] = "colocalization_analysis_v1"
    h5.attrs["created_by"] = Path(__file__).name

    h5.create_dataset(
        "metadata/json",
        data=np.bytes_(
            json.dumps(metadata, indent=4)
        )
    )

    store_image(h5, "object_a_original", a_original)
    store_image(h5, "object_b_original", b_original)

    # Raw, deduplicated object pixel tables.
    store_dataframe(h5, "object_a_pixels", object_a)
    store_dataframe(h5, "object_b_pixels", object_b)

    # Core analysis tables.
    store_dataframe(h5, "overlap_results", overlap_df)
    store_dataframe(h5, "eligible_pairs", eligible_df)
    store_dataframe(h5, "preference_matrix", preference_df)
    store_dataframe(h5, "stable_matching", matching_df)
    store_dataframe(h5, "cell_id_mapping", cell_id_mapping)
    store_dataframe(h5, "filtered_object_a_pixels", filtered_a)
    store_dataframe(h5, "filtered_object_b_pixels", filtered_b)
    store_dataframe(h5, "positive_objects", positive_df)
    store_dataframe(h5, "unmatched_object_a", unmatched_df)
    store_dataframe(h5, "summary_statistics", summary_df)
    store_dataframe(h5, "filtering_stage_counts", filtering_stage_df)
    if NUCLEAR_IDENTITY_CONSTRAINT:
        store_dataframe(
            h5,
            "nuclear_assignments_object_a",
            constraint_context_a["assignments"],
        )
        store_dataframe(
            h5,
            "nuclear_assignments_object_b",
            constraint_context_b["assignments"],
        )
        store_dataframe(h5, "nuclear_constraint_audit", constraint_audit_df)


# ============================================================
# PRINT SUMMARY
# ============================================================

print()
print("=" * 70)
print("SUMMARY STATISTICS")
print("=" * 70)

print(f"\nTotal {OBJECT_A_NAME} objects: {total_a_objects:,}")
print(f"Matched {OBJECT_A_NAME} objects: {n_matched_a:,}")
print(
    f"Proportion of {OBJECT_A_NAME}: "
    f"{a_proportion:.4f} ({a_proportion * 100:.2f}%)"
)

print(f"\nTotal {OBJECT_B_NAME} objects: {total_b_objects:,}")
print(f"Matched {OBJECT_B_NAME} objects: {n_matched_b:,}")
print(
    f"Proportion of {OBJECT_B_NAME}: "
    f"{b_proportion:.4f} ({b_proportion * 100:.2f}%)"
)

print(f"\nMatched CellIDs: {len(cell_id_mapping):,}")
print(
    (
        f"{OBJECT_A_NAME} without a constrained-eligible partner: "
        if NUCLEAR_IDENTITY_CONSTRAINT
        else f"{OBJECT_A_NAME} below threshold: "
    )
    + f"{len(below_threshold_a_ids):,}"
)

if pd.notna(organoid_area_mm2):
    print(f"\nOrganoid area: {organoid_area_pixels:,.0f} pixels")
    print(f"Organoid area: {organoid_area_mm2:.4f} mm^2")
    print(
        f"{OBJECT_A_NAME} density: "
        f"{a_density:.2f} objects/mm^2"
    )
    print(
        f"{OBJECT_B_NAME} density: "
        f"{b_density:.2f} objects/mm^2"
    )
    print(
        f"Matched {OBJECT_A_NAME} density: "
        f"{matched_a_density:.2f} objects/mm^2"
    )
    print(
        f"Matched {OBJECT_B_NAME} density: "
        f"{matched_b_density:.2f} objects/mm^2"
    )

print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
print("\nThe visualization script now needs ONLY:")
print(RESULTS_FILE)
