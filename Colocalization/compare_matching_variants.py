#!/usr/bin/env python3
"""Compare unconstrained and nucleus-constrained marker matching results."""

from __future__ import annotations

import argparse
import io
import json
from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


def read_table(h5: h5py.File, name: str) -> pd.DataFrame:
    key = f"tables/{name}"
    if key not in h5:
        raise ValueError(f"{h5.filename}: missing {key}")
    return pd.read_csv(io.BytesIO(h5[key][()].tobytes()))


def load_analysis(path: Path) -> dict:
    with h5py.File(path, "r") as h5:
        result = {
            "metadata": json.loads(h5["metadata/json"][()].decode("utf-8")),
            "object_a": read_table(h5, "object_a_pixels"),
            "object_b": read_table(h5, "object_b_pixels"),
            "overlap": read_table(h5, "overlap_results"),
            "eligible": read_table(h5, "eligible_pairs"),
            "matching": read_table(h5, "stable_matching"),
        }
        if "tables/nuclear_constraint_audit" in h5:
            result["audit"] = read_table(h5, "nuclear_constraint_audit")
    return result


def canonical_pixels(frame: pd.DataFrame, label: str) -> pd.DataFrame:
    required = {"ObjectID", "X", "Y"}
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"{label} missing {sorted(missing)}")
    result = frame[["ObjectID", "X", "Y"]].copy()
    for column in result.columns:
        result[column] = pd.to_numeric(result[column], errors="raise").astype(int)
    return result.drop_duplicates().sort_values(["ObjectID", "Y", "X"]).reset_index(drop=True)


def bool_series(frame: pd.DataFrame, column: str) -> pd.Series:
    if column not in frame.columns:
        raise ValueError(f"Missing boolean column {column}")
    values = frame[column]
    if values.dtype == bool:
        return values
    normalized = values.astype(str).str.strip().str.lower()
    invalid = ~normalized.isin({"true", "false", "1", "0"})
    if invalid.any():
        raise ValueError(f"Column {column} contains non-boolean values")
    return normalized.isin({"true", "1"})


def pair_keys(frame: pd.DataFrame) -> set[tuple[int, int]]:
    if frame.empty:
        return set()
    return set(
        zip(
            pd.to_numeric(frame["A_ObjectID"], errors="raise").astype(int),
            pd.to_numeric(frame["B_ObjectID"], errors="raise").astype(int),
        )
    )


def build_comparison(unconstrained: dict, constrained: dict):
    unconstrained_meta = unconstrained["metadata"]
    constrained_meta = constrained["metadata"]
    if constrained_meta.get("analysis_variant") != "nucleus_constrained_marker_marker":
        raise ValueError("Constrained HDF5 lacks nucleus-constrained analysis metadata")
    if unconstrained_meta["objects"] != constrained_meta["objects"]:
        raise ValueError("Matching variants use different marker names/directions")
    for key in ("object_a", "object_b"):
        if not canonical_pixels(unconstrained[key], f"unconstrained {key}").equals(
            canonical_pixels(constrained[key], f"constrained {key}")
        ):
            raise ValueError(f"Matching variants do not use identical {key} ROI pixels")
    if unconstrained_meta["matching"]["threshold"] != constrained_meta["matching"]["threshold"]:
        raise ValueError("Matching variants use different directional overlap thresholds")

    audit = constrained.get("audit")
    if audit is None:
        raise ValueError("Constrained HDF5 lacks nuclear_constraint_audit")
    required = {
        "A_ObjectID", "B_ObjectID", "A_Stage07_Source_CellID", "B_Stage07_Source_CellID",
        "A_Canonical_Nuclear_ObjectID", "B_Canonical_Nuclear_ObjectID",
        "A_Mapping_Status", "B_Mapping_Status", "Original_Overlap_Threshold_Pass",
        "Final_Constrained_Eligibility", "Final_Constrained_Stable_Match", "Exclusion_Reason",
    }
    missing = required - set(audit.columns)
    if missing:
        raise ValueError(f"Constraint audit missing {sorted(missing)}")
    audit_pairs = pair_keys(audit)
    unconstrained_overlap_pairs = pair_keys(unconstrained["overlap"])
    constrained_overlap_pairs = pair_keys(constrained["overlap"])
    if audit_pairs != unconstrained_overlap_pairs or audit_pairs != constrained_overlap_pairs:
        raise ValueError(
            "The two variants and the constraint audit do not contain the same "
            "spatial-overlap candidate pairs"
        )

    spatial_pass = bool_series(audit, "Original_Overlap_Threshold_Pass")
    constrained_eligible = bool_series(audit, "Final_Constrained_Eligibility")
    resolved = (
        audit["A_Mapping_Status"].eq("resolved")
        & audit["B_Mapping_Status"].eq("resolved")
    )
    same_nucleus = (
        pd.to_numeric(audit["A_Canonical_Nuclear_ObjectID"], errors="coerce")
        == pd.to_numeric(audit["B_Canonical_Nuclear_ObjectID"], errors="coerce")
    )
    different_nuclei = spatial_pass & resolved & ~same_nucleus
    unresolved = spatial_pass & ~resolved

    unconstrained_final = pair_keys(unconstrained["matching"])
    constrained_final = pair_keys(constrained["matching"])
    common = unconstrained_final & constrained_final
    only_unconstrained = unconstrained_final - constrained_final
    only_constrained = constrained_final - unconstrained_final
    audit_lookup = audit.set_index(["A_ObjectID", "B_ObjectID"], drop=False)

    pair_rows = []
    for a_id, b_id in sorted(unconstrained_final | constrained_final):
        if (a_id, b_id) not in audit_lookup.index:
            raise ValueError(f"Final pair {(a_id, b_id)} is absent from constraint audit")
        row = audit_lookup.loc[(a_id, b_id)]
        if isinstance(row, pd.DataFrame):
            raise ValueError(f"Duplicate constraint audit row for pair {(a_id, b_id)}")
        record = row.to_dict()
        in_unconstrained = (a_id, b_id) in unconstrained_final
        in_constrained = (a_id, b_id) in constrained_final
        record["Unconstrained_Final_Stable_Match"] = in_unconstrained
        record["Constrained_Final_Stable_Match"] = in_constrained
        record["Pair_Category"] = (
            "common" if in_unconstrained and in_constrained
            else "unconstrained_only" if in_unconstrained
            else "constrained_only"
        )
        record["Nuclear_Identity_Discordant"] = bool(
            row["A_Mapping_Status"] == "resolved"
            and row["B_Mapping_Status"] == "resolved"
            and float(row["A_Canonical_Nuclear_ObjectID"])
            != float(row["B_Canonical_Nuclear_ObjectID"])
        )
        pair_rows.append(record)
    pair_detail_columns = list(audit.columns) + [
        "Unconstrained_Final_Stable_Match",
        "Constrained_Final_Stable_Match",
        "Pair_Category",
        "Nuclear_Identity_Discordant",
    ]
    pair_details = pd.DataFrame(pair_rows, columns=pair_detail_columns)

    discordant_unconstrained = (
        pair_details[
            pair_details["Unconstrained_Final_Stable_Match"]
            & pair_details["Nuclear_Identity_Discordant"]
        ].copy()
        if len(pair_details)
        else pair_details.copy()
    )
    discordant_count = len(discordant_unconstrained)
    discordant_percent = (
        100.0 * discordant_count / len(unconstrained_final)
        if unconstrained_final else 0.0
    )

    metrics = [
        ("Marker A input ROIs", unconstrained["object_a"]["ObjectID"].nunique(), "ROIs"),
        ("Marker B input ROIs", unconstrained["object_b"]["ObjectID"].nunique(), "ROIs"),
        ("Candidate spatial-overlap pairs", len(unconstrained["overlap"]), "pairs"),
        ("Pairs passing original overlap threshold", int(spatial_pass.sum()), "pairs"),
        ("Threshold-passing pairs excluded for different nuclei", int(different_nuclei.sum()), "pairs"),
        ("Threshold-passing pairs with unresolved nuclear identity", int(unresolved.sum()), "pairs"),
        ("Pairs passing overlap and nuclear constraints", int(constrained_eligible.sum()), "pairs"),
        ("Final unconstrained stable matches", len(unconstrained_final), "pairs"),
        ("Final nucleus-constrained stable matches", len(constrained_final), "pairs"),
        ("Common final pairs", len(common), "pairs"),
        ("Final pairs unique to unconstrained", len(only_unconstrained), "pairs"),
        ("Final pairs unique to nucleus-constrained", len(only_constrained), "pairs"),
        ("Unconstrained final pairs with discordant nuclear identity", discordant_count, "pairs"),
        ("Percent unconstrained final pairs with discordant nuclear identity", discordant_percent, "percent"),
    ]
    summary = pd.DataFrame(metrics, columns=["Metric", "Value", "Units"])
    return summary, pair_details, discordant_unconstrained


def save_figure(summary: pd.DataFrame, object_a: str, object_b: str, pdf_path: Path, png_path: Path):
    count_metrics = summary[summary["Units"] == "pairs"].copy()
    fig, ax = plt.subplots(figsize=(13, 8.5))
    ax.barh(count_metrics["Metric"], count_metrics["Value"], color="#4472C4")
    ax.invert_yaxis()
    ax.set_xlabel("Number of marker pairs")
    ax.set_title(
        f"{object_a}–{object_b}: unconstrained vs NUCLEUS-CONSTRAINED stable matching",
        weight="bold",
    )
    for index, value in enumerate(count_metrics["Value"]):
        ax.text(float(value), index, f"  {float(value):g}", va="center", fontsize=9)
    percent = float(
        summary.loc[
            summary["Metric"].eq("Percent unconstrained final pairs with discordant nuclear identity"),
            "Value",
        ].iloc[0]
    )
    fig.text(
        0.01,
        0.01,
        f"Discordant nuclear identity among unconstrained final pairs: {percent:.2f}%. "
        "Constrained matching is recomputed independently after candidate restriction; it is not a post-filtered subset.",
        fontsize=9,
    )
    fig.tight_layout(rect=(0, 0.055, 1, 1))
    fig.savefig(pdf_path, bbox_inches="tight")
    fig.savefig(png_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--unconstrained-h5", type=Path, required=True)
    parser.add_argument("--constrained-h5", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args(argv)

    unconstrained = load_analysis(args.unconstrained_h5)
    constrained = load_analysis(args.constrained_h5)
    summary, details, discordant = build_comparison(unconstrained, constrained)
    output = args.output_dir
    output.mkdir(parents=True, exist_ok=True)
    summary.to_csv(output / "matching_variant_comparison.csv", index=False)
    details.to_csv(output / "matching_variant_pair_details.csv", index=False)
    discordant.to_csv(output / "discordant_unconstrained_final_pairs.csv", index=False)
    names = unconstrained["metadata"]["objects"]
    save_figure(
        summary,
        names["object_a_name"],
        names["object_b_name"],
        output / "matching_variant_comparison.pdf",
        output / "matching_variant_comparison.png",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
