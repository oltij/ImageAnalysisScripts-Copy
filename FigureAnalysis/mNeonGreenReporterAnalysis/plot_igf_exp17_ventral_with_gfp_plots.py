#!/usr/bin/env python3
"""Experiment 17 ventral-only H1 / CL23 / CL32 IGF-noIGF comparisons.

Preserves the previously verified H1 fields; selects all complete Experiment 17
ventral CL23/CL32 fields from both production workspaces using the explicit
block-ID inventory. Automatically includes block 153 (untreated CL23).
Never includes experiment 16, dorsal, or excluded failed indices 12/13.
Produces 14 PV/LHX6 and 6 mNeonGreen/GFP comparison figures plus a multipage PDF.
Read-only: no upstream reruns.
"""
from __future__ import annotations

import argparse
import csv
import json
import textwrap
from pathlib import Path
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_SAMPLE_METADATA = SCRIPT_DIR / "igf_exp17_ventral_block_mapping.csv"

CLASSES = ("LHX6+/PV+", "LHX6+/PV-", "LHX6-/PV+", "LHX6-/PV-")
CONDITIONS = ("No IGF", "+IGF")
CELL_LINE_ORDER = ("H1", "CL23", "CL32")
EXPECTED_BLOCKS = {149, 150, 151, 152, 153, 157}
EXCLUDED_INDICES = {12, 13}

ANTI_GFP_ORDER = ("Anti-GFP stained (-2)", "No anti-GFP (-3)")

GFP_PLOT_METRICS = (
    ("gfp_density_per_mm2", "mNeonGreen / GFP ROI density", "GFP ROIs / mm² Hoechst area", "gfp_density"),
    ("gfp_area_percent_of_organoid", "mNeonGreen / GFP ROI area", "% of Hoechst organoid area", "gfp_area_percent"),
    ("gfp_area_mean_um2", "Mean mNeonGreen / GFP ROI area", "µm² / GFP ROI", "gfp_area_mean"),
    ("gfp_area_median_um2", "Median mNeonGreen / GFP ROI area", "µm² / GFP ROI", "gfp_area_median"),
    ("gfp_mean_intensity_mean", "Mean mNeonGreen / GFP ROI intensity", "Mean intensity / GFP ROI", "gfp_intensity_mean"),
    ("gfp_mean_intensity_median", "Median mNeonGreen / GFP ROI intensity", "Mean intensity / GFP ROI", "gfp_intensity_median"),
)

PLOT_METRICS = (
    ("pv_area_percent_of_organoid_nuclear_only", "PV+ ROI area, nuclear-associated only", "% of Hoechst organoid area", "pv_nuclear_area_percent"),
    ("pv_area_mean_um2_per_positive_nucleus", "Mean PV+ ROI area per PV+ nucleus", "µm² / PV+ nucleus", "pv_nuclear_mean_area"),
    ("pv_area_median_um2_per_positive_nucleus", "Median PV+ ROI area per PV+ nucleus", "µm² / PV+ nucleus", "pv_nuclear_median_area"),
    ("pv_area_mean_um2_double_positive", "Mean PV ROI area in LHX6+/PV+ nuclei", "µm² / double-positive nucleus", "pv_area_double_positive"),
    ("pv_area_mean_um2_pv_only", "Mean PV ROI area in LHX6−/PV+ nuclei", "µm² / PV-only nucleus", "pv_area_pv_only"),
    ("lhx6_pv_double_positive_density_per_mm2", "LHX6+/PV+ cells", "Cells / mm² Hoechst area", "double_positive_density"),
    ("lhx6_only_density_per_mm2", "LHX6+/PV− cells", "Cells / mm² Hoechst area", "lhx6_only_density"),
    ("pv_only_density_per_mm2", "LHX6−/PV+ cells", "Cells / mm² Hoechst area", "pv_only_density"),
    ("lhx6_positive_density_per_mm2", "All LHX6+ cells", "Cells / mm² Hoechst area", "lhx6_total_density"),
    ("pv_within_lhx6_percent", "PV+ among LHX6+ cells", "% of LHX6+ nuclei", "pv_within_lhx6"),
    ("pv_positive_density_per_mm2", "All PV+ cells", "Cells / mm² Hoechst area", "pv_total_density"),
    ("double_positive_percent_of_all_nuclei", "LHX6+/PV+ nuclei", "% of all Hoechst nuclei", "double_positive_fraction"),
    ("lhx6_only_percent_of_all_nuclei", "LHX6+/PV− nuclei", "% of all Hoechst nuclei", "lhx6_only_fraction"),
    ("pv_only_percent_of_all_nuclei", "LHX6−/PV+ nuclei", "% of all Hoechst nuclei", "pv_only_fraction"),
)


METRIC_LEGENDS = {
    "pv_area_percent_of_organoid_nuclear_only": "Total area of PV-positive ROIs associated with nuclei, expressed as a percentage of the Hoechst-defined organoid area.",
    "pv_area_mean_um2_per_positive_nucleus": "Average PV ROI area for PV-positive nuclei; larger values mean larger PV-associated segmented regions per PV-positive cell.",
    "pv_area_median_um2_per_positive_nucleus": "Median PV ROI area for PV-positive nuclei; this is less sensitive to unusually large PV ROIs than the mean.",
    "pv_area_mean_um2_double_positive": "Average PV ROI area specifically among LHX6+/PV+ nuclei.",
    "pv_area_mean_um2_pv_only": "Average PV ROI area among LHX6-/PV+ nuclei (PV-only cells).",
    "lhx6_pv_double_positive_density_per_mm2": "Density of nuclei positive for both LHX6 and PV, normalized to Hoechst-defined organoid area.",
    "lhx6_only_density_per_mm2": "Density of LHX6+/PV- nuclei, normalized to Hoechst-defined organoid area.",
    "pv_only_density_per_mm2": "Density of LHX6-/PV+ nuclei, normalized to Hoechst-defined organoid area.",
    "lhx6_positive_density_per_mm2": "Density of all LHX6+ nuclei, regardless of PV status.",
    "pv_within_lhx6_percent": "Percentage of LHX6+ nuclei that are also PV+, i.e. LHX6+/PV+ divided by all LHX6+ nuclei.",
    "pv_positive_density_per_mm2": "Density of all PV+ nuclei, regardless of LHX6 status.",
    "double_positive_percent_of_all_nuclei": "Percentage of all Hoechst nuclei classified as LHX6+/PV+.",
    "lhx6_only_percent_of_all_nuclei": "Percentage of all Hoechst nuclei classified as LHX6+/PV-.",
    "pv_only_percent_of_all_nuclei": "Percentage of all Hoechst nuclei classified as LHX6-/PV+.",
    "gfp_density_per_mm2": "Density of final kept mNeonGreen / GFP ROIs per mm² of Hoechst-defined organoid area.",
    "gfp_area_percent_of_organoid": "Total area of final kept mNeonGreen / GFP ROIs expressed as a percentage of the Hoechst-defined organoid area.",
    "gfp_area_mean_um2": "Average area of final kept mNeonGreen / GFP ROIs.",
    "gfp_area_median_um2": "Median area of final kept mNeonGreen / GFP ROIs.",
    "gfp_mean_intensity_mean": "Average ROI mean intensity for final kept mNeonGreen / GFP ROIs.",
    "gfp_mean_intensity_median": "Median ROI mean intensity for final kept mNeonGreen / GFP ROIs."
}

GENERAL_LEGEND_TEXT = (
    "Legend: each dot = one microscopy field; dot color = one sample/block-section; black horizontal line = condition mean. "
    "Panels show H1, CL23, and CL32 side by side. Cell classes are defined by nuclear colocalization with Hoechst: "
    "LHX6+/PV+, LHX6+/PV-, LHX6-/PV+, LHX6-/PV-. PV area means the area of a segmented PV-positive ROI associated with a nucleus; "
    "it is not PV signal physically restricted to the nucleus. Density plots are normalized to Hoechst-defined organoid area; "
    "fraction plots are normalized to all Hoechst nuclei. The mNeonGreen/GFP comparison pages use final kept GFP ROIs after filtering; they are not restricted to nuclei. "
    "Anti-GFP status is inferred from the sample suffix: -2 = anti-GFP stained, -3 = no anti-GFP. All included samples are Experiment 17 ventral organoids."
)


def require_file(path: Path) -> Path:
    if not path.is_file():
        raise FileNotFoundError(f"Required file not found: {path}")
    return path


def scalar_float(value, label: str) -> float:
    try:
        v = float(value)
    except (ValueError, TypeError) as e:
        raise ValueError(f"Missing or nonnumeric {label}: {value!r}") from e
    if not np.isfinite(v):
        raise ValueError(f"Invalid {label}: {value!r}")
    return v


def identify_sample(value: str) -> str:
    found = re.findall(r"(?<!\d)(\d{3}-\d)(?!\d)", str(value))
    unique = sorted(set(found))
    if len(unique) != 1:
        raise ValueError(f"Cannot identify unique sample in {value!r}: {unique}")
    return unique[0]


def load_sample_info(path: Path) -> dict[int, dict[str, str]]:
    """Inventory is keyed by block ID; suffixes (-2/-3) are acquisition sections."""
    frame = pd.read_csv(require_file(path), dtype=str).fillna("")
    needed = {"block_id", "cell_line", "condition", "experiment_id", "organoid_type"}
    if not needed.issubset(frame.columns):
        raise ValueError(f"Block metadata CSV must contain columns {sorted(needed)}")
    info: dict[int, dict[str, str]] = {}
    for row in frame.to_dict("records"):
        block = int(row["block_id"])
        if block in info:
            raise ValueError(f"Duplicate block ID in metadata: {block}")
        cell_line = str(row["cell_line"]).strip()
        condition = str(row["condition"]).strip()
        experiment = str(row["experiment_id"]).strip()
        tissue = str(row["organoid_type"]).strip()
        if (condition not in CONDITIONS or cell_line not in CELL_LINE_ORDER
                or experiment != "17" or tissue.lower() != "ventral"):
            raise ValueError(f"Invalid exp17 ventral cell line / condition for block {block}: {row}")
        info[block] = {
            "cell_line": cell_line, "condition": condition,
            "experiment_id": experiment, "organoid_type": tissue,
            "condition_basis": str(row.get("condition_basis", "")).strip(),
            "inventory_note": str(row.get("inventory_note", "")).strip(),
        }
    if set(info) != EXPECTED_BLOCKS:
        raise ValueError(f"Expected exactly Experiment 17 ventral blocks {sorted(EXPECTED_BLOCKS)}, got {sorted(info)}")
    expected = {149: ("CL32", "No IGF"), 150: ("CL32", "+IGF"),
                151: ("H1", "+IGF"), 152: ("CL23", "+IGF"),
                153: ("CL23", "No IGF"), 157: ("H1", "No IGF")}
    for block, (line, condition) in expected.items():
        if (info[block]["cell_line"], info[block]["condition"]) != (line, condition):
            raise ValueError(f"Block {block} contradicts inventory mapping: expected {line}/{condition}")
    return info


def source_identifiers(job: dict) -> tuple[int, int, str]:
    source = str(job.get("source_field", ""))
    # Example: Exp17_153-2_SetN_DAPI_GFP_LHX6_PV_Field_1
    match = re.search(r"(?:^|/)Exp(?P<exp>\d+)_(?P<block>\d{3})-(?P<section>\d)(?!\d)", source)
    if not match:
        raise ValueError(f"Cannot identify experiment/block in source_field: {source!r}")
    return int(match["exp"]), int(match["block"]), f"{match['block']}-{match['section']}"


def load_h1_selection(path: Path):
    with require_file(path).open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    if not rows or "index" not in rows[0]:
        raise ValueError(f"No rows or 'index' missing in {path}")
    indices = [int(row["index"]) for row in rows]
    if len(indices) != len(set(indices)):
        raise ValueError("Duplicate verified indices in H1 selection CSV")
    if any(i in (12, 13) or i < 0 or i > 47 for i in indices):
        raise ValueError("H1 selection must stay within verified original H1 fields (0-47 excluding 12/13)")
    return rows


def read_counts(path: Path) -> dict[str, int]:
    frame = pd.read_csv(require_file(path))
    required = {"Class", "Count"}
    if not required.issubset(frame.columns):
        raise ValueError(f"Missing {required - set(frame.columns)} in {path}")
    if frame["Class"].duplicated().any():
        raise ValueError(f"Duplicate class labels in {path}")
    counts = {}
    for name in CLASSES:
        subset = frame.loc[frame["Class"] == name, "Count"]
        if len(subset) != 1:
            raise ValueError(f"Exactly one count required for {name}: {path}")
        v = scalar_float(subset.iloc[0], f"{name} count")
        if v < 0 or not v.is_integer():
            raise ValueError(f"Invalid cell count {name}={v}")
        counts[name] = int(v)
    if len(frame) != 4:
        raise ValueError(f"Expected exactly four nuclear classes in {path}; got {len(frame)}")
    return counts


def anti_gfp_status(sample: str) -> str:
    if str(sample).endswith("-2"):
        return "Anti-GFP stained (-2)"
    if str(sample).endswith("-3"):
        return "No anti-GFP (-3)"
    raise ValueError(f"Cannot infer anti-GFP status from sample {sample!r}; expected -2 or -3 suffix")


def choose_intensity_column(frame: pd.DataFrame) -> str | None:
    candidates = [
        "Mean intensity", "MeanIntensity", "Mean_Intensity", "Intensity_Mean",
        "MeanIntensity_GFP", "IntegratedIntensity"
    ]
    for c in candidates:
        if c in frame.columns:
            return c
    return None


def extract_gfp_metrics(analysis: Path, organoid_um2: float) -> dict[str, float]:
    path = analysis / "06_final_shape_filtered" / "GFP" / "GFP_kept_ROI_metrics.csv"
    frame = pd.read_csv(require_file(path))
    if "Area_um2" not in frame.columns:
        raise ValueError(f"GFP metrics missing Area_um2: {path}")
    areas = pd.to_numeric(frame["Area_um2"], errors="coerce")
    if areas.isna().any():
        raise ValueError(f"Non-numeric GFP Area_um2 values in {path}")
    n = int(len(frame))
    area_sum = float(areas.sum())
    organoid_mm2 = organoid_um2 / 1e6
    metrics = {
        "gfp_roi_count": n,
        "gfp_density_per_mm2": (n / organoid_mm2 if organoid_mm2 else np.nan),
        "gfp_area_sum_um2": area_sum,
        "gfp_area_percent_of_organoid": (100.0 * area_sum / organoid_um2 if organoid_um2 else np.nan),
        "gfp_area_mean_um2": float(areas.mean()) if n else np.nan,
        "gfp_area_median_um2": float(areas.median()) if n else np.nan,
    }
    intensity_col = choose_intensity_column(frame)
    if intensity_col is not None:
        intens = pd.to_numeric(frame[intensity_col], errors="coerce")
        metrics["gfp_mean_intensity_mean"] = float(intens.mean()) if len(intens) else np.nan
        metrics["gfp_mean_intensity_median"] = float(intens.median()) if len(intens) else np.nan
        metrics["gfp_intensity_column"] = intensity_col
    else:
        metrics["gfp_mean_intensity_mean"] = np.nan
        metrics["gfp_mean_intensity_median"] = np.nan
        metrics["gfp_intensity_column"] = ""
    return metrics


def extract_pv_area(analysis: Path, nuclear_csv: Path, expected_pv: int):
    nuclei = pd.read_csv(require_file(nuclear_csv))
    wanted = {"Nuclear_ObjectID", "Class", "Object_B_Name", "Object_B_Source_Marker_ObjectID"}
    if not wanted.issubset(nuclei.columns):
        raise ValueError(f"Missing all-nuclei identity columns {wanted - set(nuclei.columns)} in {nuclear_csv}")
    if nuclei["Nuclear_ObjectID"].duplicated().any():
        raise ValueError(f"Duplicate nucleus IDs: {nuclear_csv}")
    pv_rows = nuclei[nuclei["Class"].isin(("LHX6+/PV+", "LHX6-/PV+"))].copy()
    if len(pv_rows) != expected_pv:
        raise ValueError(f"PV nuclear count mismatch: {len(pv_rows)} != {expected_pv}")
    if (pv_rows["Object_B_Name"] != "PV").any():
        raise ValueError("Expected Object B to be PV in the marker-pair all-cell classification")
    source_ids = pd.to_numeric(pv_rows["Object_B_Source_Marker_ObjectID"], errors="coerce")
    if source_ids.isna().any() or (source_ids % 1 != 0).any():
        raise ValueError("PV+ nuclei lack valid original PV marker object IDs")
    if source_ids.duplicated().any():
        raise ValueError("Multiple PV+ nuclei refer to the same PV marker ROI")
    pv_rows["PV_Source_ObjectID"] = source_ids.astype(int)

    metrics_path = analysis / "06_final_shape_filtered" / "PV" / "PV_kept_ROI_metrics.csv"
    metrics = pd.read_csv(require_file(metrics_path))
    if not {"ObjectID", "Area_um2"}.issubset(metrics.columns):
        raise ValueError(f"PV final metrics missing ObjectID or Area_um2: {metrics_path}")
    if metrics["ObjectID"].duplicated().any():
        raise ValueError(f"Duplicate PV ObjectID in {metrics_path}")

    merged = pv_rows.merge(
        metrics[["ObjectID", "Area_um2"]],
        how="left",
        left_on="PV_Source_ObjectID",
        right_on="ObjectID",
        validate="one_to_one",
    )
    if merged["Area_um2"].isna().any() or (merged["Area_um2"] <= 0).any():
        bad = merged.loc[merged["Area_um2"].isna() | (merged["Area_um2"] <= 0), "PV_Source_ObjectID"].tolist()
        raise ValueError(f"PV nuclear-positive source ROI IDs missing positive final area metrics: {bad[:10]}")
    return merged[["Nuclear_ObjectID", "Class", "PV_Source_ObjectID", "Area_um2"]].rename(columns={"Area_um2": "PV_ROI_Area_um2"})


def get_field(job: dict, index: int, sample_info: dict[int, dict[str, str]], original_h1_row: dict | None = None):
    source = str(job.get("source_field", ""))
    exp, block, sample = source_identifiers(job)
    anti_gfp = anti_gfp_status(sample)
    if exp != 17 or block not in sample_info:
        raise ValueError(f"Refusing non-Exp17-ventral or unknown block: {source}")
    cell_line = sample_info[block]["cell_line"]
    condition = sample_info[block]["condition"]
    experiment_id = "17"
    condition_basis = sample_info[block].get("condition_basis", "")
    inventory_note = sample_info[block].get("inventory_note", "")

    if original_h1_row is not None:
        if "sample" in original_h1_row and str(original_h1_row["sample"]).strip():
            old_sample = identify_sample(str(original_h1_row["sample"]))
            if sample != old_sample:
                raise ValueError(f"Index {index}: plotted sample {old_sample} != manifest sample {sample}")
        if "name" in original_h1_row and str(original_h1_row["name"]).strip() and str(original_h1_row["name"]).strip() != job.get("name"):
            raise ValueError(f"Index {index}: manifest name differs from verified H1 selection")

    analysis = Path(job["output_dir"]).resolve() / "analysis"
    prefix = analysis / "08_marker_colocalization" / "LHX6_vs_PV" / "visualizations" / "unmatched_only_qc"
    counts = read_counts(prefix / "LHX6_PV_all_nuclear_cells_class_counts.csv")
    n_all = sum(counts.values())
    dp, lx, pv, dn = (counts[label] for label in CLASSES)
    pv_total = dp + pv
    lhx6_total = dp + lx

    summary_file = prefix / "LHX6_PV_final_summary_metrics.csv"
    summary = pd.read_csv(require_file(summary_file))
    if len(summary) != 1:
        raise ValueError(f"Expected one row in {summary_file}")
    summary = summary.iloc[0]
    if "All_Cell_Total_Nuclear_Cells" in summary and int(scalar_float(summary["All_Cell_Total_Nuclear_Cells"], "total nuclei")) != n_all:
        raise ValueError(f"Index {index}: class counts != canonical nuclei total")
    organoid_um2 = scalar_float(summary["Organoid_Area_um2"], "organoid area")
    if organoid_um2 <= 0:
        raise ValueError("Organoid area must be positive")
    organoid_mm2 = organoid_um2 / 1e6

    if original_h1_row is not None:
        if "organoid_area_um2" in original_h1_row and str(original_h1_row["organoid_area_um2"]).strip():
            old_area = scalar_float(original_h1_row["organoid_area_um2"], "previous plotted organoid area")
            if not np.isclose(organoid_um2, old_area, rtol=1e-5, atol=0.1):
                raise ValueError(f"Index {index}: area mismatch with previous PV plot: {organoid_um2} versus {old_area}")
        if "pv_positive_cells" in original_h1_row and str(original_h1_row["pv_positive_cells"]).strip():
            old_pv = int(float(original_h1_row["pv_positive_cells"]))
            if old_pv != pv_total:
                raise ValueError(f"Index {index}: PV count differs from previous plot: {pv_total} versus {old_pv}")

    area_df = extract_pv_area(analysis, prefix / "LHX6_PV_all_nuclear_cells_classification.csv", pv_total)
    area_df.insert(0, "index", index)
    area_df.insert(1, "sample", sample)
    area_df.insert(2, "block_id", block)
    area_df.insert(3, "cell_line", cell_line)
    area_df.insert(4, "condition", condition)
    area_df.insert(5, "experiment_id", experiment_id)
    area_df.insert(6, "anti_gfp_status", anti_gfp)
    area_sum = float(area_df["PV_ROI_Area_um2"].sum())
    gfp_metrics = extract_gfp_metrics(analysis, organoid_um2)

    result = {
        "index": index,
        "sample": sample,
        "block_id": block,
        "organoid_type": "Ventral",
        "cell_line": cell_line,
        "condition": condition,
        "experiment_id": experiment_id,
        "anti_gfp_status": anti_gfp,
        "condition_basis": condition_basis,
        "inventory_note": inventory_note,
        "source_field": source,
        "output_dir": str(analysis),
        "organoid_area_um2": organoid_um2,
        "organoid_area_mm2": organoid_mm2,
        "all_hoechst_nuclei": n_all,
        "lhx6_pv_double_positive_cells": dp,
        "lhx6_only_cells": lx,
        "pv_only_cells": pv,
        "double_negative_cells": dn,
        "lhx6_positive_cells": lhx6_total,
        "pv_positive_cells": pv_total,
        "pv_nuclear_associated_roi_count": len(area_df),
        "pv_area_sum_um2_nuclear_only": area_sum,
        "pv_area_percent_of_organoid_nuclear_only": 100.0 * area_sum / organoid_um2,
        "pv_area_mean_um2_per_positive_nucleus": (float(area_df["PV_ROI_Area_um2"].mean()) if pv_total else np.nan),
        "pv_area_median_um2_per_positive_nucleus": (float(area_df["PV_ROI_Area_um2"].median()) if pv_total else np.nan),
        "pv_area_mean_um2_double_positive": float(area_df.loc[area_df["Class"] == "LHX6+/PV+", "PV_ROI_Area_um2"].mean()),
        "pv_area_mean_um2_pv_only": float(area_df.loc[area_df["Class"] == "LHX6-/PV+", "PV_ROI_Area_um2"].mean()),
        "pv_within_lhx6_percent": (100.0 * dp / lhx6_total if lhx6_total else np.nan),
    }
    result.update(gfp_metrics)
    for key, n in (("lhx6_pv_double_positive", dp), ("lhx6_only", lx), ("pv_only", pv), ("lhx6_positive", lhx6_total), ("pv_positive", pv_total)):
        result[f"{key}_density_per_mm2"] = n / organoid_mm2
    for key, n in (("double_positive", dp), ("lhx6_only", lx), ("pv_only", pv), ("lhx6_positive", lhx6_total), ("pv_positive", pv_total)):
        result[f"{key}_percent_of_all_nuclei"] = 100.0 * n / n_all if n_all else np.nan
    return result, area_df


def verified_complete(job: dict) -> bool:
    """Check durable downstream completion marker (manifest alone is insufficient)."""
    analysis = Path(job["output_dir"]).resolve() / "analysis"
    path = analysis / "run_manifest.json"
    if not path.is_file():
        return False
    completed = json.loads(path.read_text()).get("completed", {})
    return any(k.startswith("08_pair_") and k.endswith("_viz") for k in completed)


def collect_fields(original_workspace: Path, remaining_workspace: Path,
                   selection_csv: Path, sample_info: dict[int, dict[str, str]]):
    """All finished CL23/CL32 Exp17 ventral fields; exact verified H1 selection."""
    original = json.loads(require_file(original_workspace / "manifest.json").read_text())
    remaining = json.loads(require_file(remaining_workspace / "manifest.json").read_text())
    old_jobs, new_jobs = original["jobs"], remaining["jobs"]
    if len(old_jobs) != 94 or len(new_jobs) != 94:
        raise ValueError("Both production manifests must describe the same 94 fields")
    h1_rows = load_h1_selection(selection_csv)
    h1_by_index = {int(row["index"]): row for row in h1_rows}

    all_fields, all_pv, audit = [], [], []
    seen_h1 = set()
    for i in range(94):
        old, new = old_jobs[i], new_jobs[i]
        if old.get("index") != i or new.get("index") != i:
            raise ValueError(f"Manifest index mismatch at {i}")
        if str(old.get("source_field")) != str(new.get("source_field")):
            raise ValueError(f"Source field differs between original/new manifests at {i}")
        old_tiles = [(t["source_index"], t["source"]) for t in old.get("tile_map", [])]
        new_tiles = [(t["source_index"], t["source"]) for t in new.get("tile_map", [])]
        if old_tiles != new_tiles:
            raise ValueError(f"Original tiles differ between workspaces at index {i}")

        exp, block, sample = source_identifiers(old)
        reason = ""
        if i in EXCLUDED_INDICES:
            reason = "known-failed excluded index"
        elif exp != 17:
            reason = "different experiment"
        elif block not in sample_info:
            reason = "not one of the six Experiment 17 ventral blocks"
        elif sample_info[block]["cell_line"] == "H1" and i not in h1_by_index:
            reason = "not in original verified H1 plot selection"

        entry = {"index": i, "source_field": old["source_field"], "sample": sample,
                 "block_id": block, "experiment_id": exp,
                 "cell_line": sample_info[block]["cell_line"] if block in sample_info else "",
                 "condition": sample_info[block]["condition"] if block in sample_info else "",
                 "included": False, "source_workspace": "", "reason": reason}
        if reason:
            audit.append(entry)
            continue

        # Earlier fields use original workspace. Index 47 succeeded on retry,
        # and can be found in either workspace. The later array 48..93 is new.
        if i <= 46:
            candidates = [("original", old)]
        elif i == 47:
            candidates = [("original", old), ("remaining", new)]
        else:
            candidates = [("remaining", new)]
        selected = next(((tag, job) for tag, job in candidates if verified_complete(job)), None)
        if selected is None:
            entry["reason"] = "no complete downstream marker-colocalization stage"
            audit.append(entry)
            print(f"[SKIP {i:02d}] {sample}: no complete analysis in expected workspace")
            continue
        tag, job = selected
        old_h1 = h1_by_index.get(i) if sample_info[block]["cell_line"] == "H1" else None
        result, pv = get_field(job, i, sample_info, original_h1_row=old_h1)
        result["source_workspace"] = tag
        result["verified_h1_selection"] = bool(old_h1)
        all_fields.append(result)
        all_pv.append(pv)
        if old_h1:
            seen_h1.add(i)
        entry.update({"included": True, "source_workspace": tag, "reason": "complete analysis"})
        audit.append(entry)
        print(f"[OK {i:02d}] {result['cell_line']} {sample}: {result['condition']} ({tag})")

    missing_h1 = sorted(set(h1_by_index) - seen_h1)
    if missing_h1:
        raise ValueError(f"Verified H1 selection indices were not processed: {missing_h1}; inspect selection audit and original workspace")
    if not all_fields:
        raise ValueError("No finished Exp17 ventral fields were found")
    return all_fields, all_pv, audit


def draw(ax, df: pd.DataFrame, metric: str, ylabel: str, palette: dict[str, tuple], *, title: str):
    for ci, condition in enumerate(CONDITIONS):
        group = df[df["condition"] == condition]
        if group.empty:
            continue
        for sample in sorted(group["sample"].unique()):
            subset = group[group["sample"] == sample].sort_values("index")
            vals = subset[metric].to_numpy(dtype=float)
            offsets = np.array([((int(i) * 37 % 101) / 100 - .5) * .40 for i in subset["index"]])
            finite = np.isfinite(vals)
            ax.scatter(
                np.full(finite.sum(), ci) + offsets[finite],
                vals[finite],
                s=26,
                color=palette[sample],
                alpha=.82,
                edgecolors="white",
                linewidths=.3,
                zorder=3,
            )
        good = group[metric].dropna().to_numpy(dtype=float)
        if good.size:
            ax.plot([ci - .19, ci + .19], [good.mean(), good.mean()], color="black", lw=1.7, zorder=4)
    ax.set_xticks(range(len(CONDITIONS)), CONDITIONS)
    ax.set_xlim(-.55, len(CONDITIONS) - .45)
    ax.set_title(title, fontsize=10, weight="semibold")
    ax.set_ylabel(ylabel, fontsize=9)
    ax.grid(axis="y", alpha=.2, linewidth=.6)
    ax.spines[["top", "right"]].set_visible(False)


def write_plots(df: pd.DataFrame, dest: Path):
    samples = sorted(df["sample"].unique())
    cmap = plt.get_cmap("tab20")
    palette = {sample: cmap(i % 20) for i, sample in enumerate(samples)}
    cell_lines = [c for c in CELL_LINE_ORDER if c in set(df["cell_line"]) ]
    extras = [c for c in sorted(set(df["cell_line"])) if c not in cell_lines]
    cell_lines.extend(extras)

    def add_footer(fig, metric_key: str | None = None):
        wrapped_general = textwrap.fill(GENERAL_LEGEND_TEXT, width=175)
        fig.text(.5, .105, wrapped_general, ha="center", va="center", fontsize=8.2, color="black")
        if metric_key is not None:
            specific = "What this plot shows: " + METRIC_LEGENDS[metric_key]
            wrapped_specific = textwrap.fill(specific, width=165)
            fig.text(.5, .045, wrapped_specific, ha="center", va="center", fontsize=8.4, color="black")

    with PdfPages(dest / "all_cell_lines_side_by_side.pdf") as pdf:
        # Document sample assignments and key limitations before the actual plots.
        fig = plt.figure(figsize=(15, 8.5))
        fig.patch.set_facecolor("white")
        fig.text(.055, .94, "Experiment 17 ventral only · IGF vs No IGF", fontsize=20, weight="semibold")
        fig.text(.055, .885, "Subsequent pages include PV/LHX6 measurements and mNeonGreen/GFP comparisons.", fontsize=12)
        fig.text(.055, .85, "Each colored dot = one processed field; colors identify sample IDs, not independent experiments.", fontsize=11)
        groups = df.groupby(["cell_line", "sample", "condition", "experiment_id"], observed=True).size().reset_index(name="Fields")
        rows = [[str(v) for v in record] for record in groups[["cell_line", "sample", "condition", "experiment_id", "Fields"]].itertuples(index=False, name=None)]
        ax = fig.add_axes([.05, .30, .9, .47]); ax.axis("off")
        table = ax.table(cellText=rows, colLabels=["Cell line", "Sample", "IGF condition", "Experiment", "Fields"], loc="center", cellLoc="left", colLoc="left")
        table.auto_set_font_size(False); table.set_fontsize(10); table.scale(1, 1.23)
        fig.text(.055, .235, "Selected blocks: H1 = 151 (+IGF), 157 (No IGF); CL23 = 152 (+IGF), 153 (No IGF).", fontsize=11)
        fig.text(.055, .195, "CL32 = 150 (+IGF), 149 (No IGF). All groups are ventral and Experiment 17.", fontsize=11, weight="semibold")
        fig.text(.055, .155, "CL32 block 150 was stored after drying one day late (per inventory); review as a potential QC covariate.", fontsize=11)
        fig.text(.055, .115, "No inferential p-values: fields within the same sample are not independent biological replicates.", fontsize=11)
        add_footer(fig)
        pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)

        # Dedicated legend page
        fig = plt.figure(figsize=(15, 8.5))
        fig.patch.set_facecolor("white")
        fig.text(.05, .95, "Figure legend and plot definitions", fontsize=18, weight="semibold")
        legend_lines = [
            "1. PV+ ROI area, nuclear-associated only: total area of PV-positive ROIs associated with nuclei, as % of Hoechst organoid area.",
            "2. Mean PV+ ROI area per PV+ nucleus: average PV ROI area for PV-positive nuclei.",
            "3. Median PV+ ROI area per PV+ nucleus: median PV ROI area for PV-positive nuclei.",
            "4. Mean PV ROI area in LHX6+/PV+ nuclei: average PV ROI area among double-positive cells.",
            "5. Mean PV ROI area in LHX6-/PV+ nuclei: average PV ROI area among PV-only cells.",
            "6. LHX6+/PV+ cells: density of double-positive nuclei per mm² Hoechst area.",
            "7. LHX6+/PV- cells: density of LHX6-only nuclei per mm² Hoechst area.",
            "8. LHX6-/PV+ cells: density of PV-only nuclei per mm² Hoechst area.",
            "9. All LHX6+ cells: density of all LHX6-positive nuclei per mm² Hoechst area.",
            "10. PV+ among LHX6+ cells: percentage of LHX6+ nuclei that are also PV+.",
            "11. All PV+ cells: density of all PV-positive nuclei per mm² Hoechst area.",
            "12. LHX6+/PV+ nuclei: percentage of all Hoechst nuclei that are double-positive.",
            "13. LHX6+/PV- nuclei: percentage of all Hoechst nuclei that are LHX6-only.",
            "14. LHX6-/PV+ nuclei: percentage of all Hoechst nuclei that are PV-only."
        ]
        y=.90
        for line in legend_lines:
            fig.text(.055, y, textwrap.fill(line, width=145), fontsize=10.3, va="top")
            y -= .055
        add_footer(fig)
        pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)

        for metric, heading, unit, stem in PLOT_METRICS:
            fig, axes = plt.subplots(1, len(cell_lines), figsize=(5.3 * len(cell_lines), 4.8), squeeze=False, sharey=False)
            axes = axes[0]
            all_y = []
            for ax, cell_line in zip(axes, cell_lines):
                sub = df[df["cell_line"] == cell_line].copy()
                if sub.empty:
                    ax.set_axis_off()
                    continue
                draw(ax, sub, metric, unit, palette, title=f"{cell_line} (n={len(sub)} fields)")
                vals = pd.to_numeric(sub[metric], errors="coerce")
                vals = vals[np.isfinite(vals)]
                if len(vals):
                    all_y.extend(vals.tolist())
            if all_y:
                ymin, ymax = min(all_y), max(all_y)
                if np.isfinite(ymin) and np.isfinite(ymax):
                    if ymin == ymax:
                        pad = abs(ymin) * 0.08 if ymin else 1.0
                    else:
                        pad = (ymax - ymin) * 0.08
                    for ax in axes:
                        if ax.axison:
                            ax.set_ylim(ymin - pad, ymax + pad)
            handles = [plt.Line2D([0], [0], marker="o", linestyle="none", markerfacecolor=palette[s], markeredgecolor="none", label=s) for s in samples]
            fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 0.155), ncol=min(6, len(handles)), frameon=False, title="Sample (block-section)")
            fig.suptitle(heading + " · each dot = one field", fontsize=13, weight="semibold")
            add_footer(fig, metric)
            fig.subplots_adjust(top=.84, bottom=.31, wspace=.30)
            fig.savefig(dest / f"{stem}_side_by_side.png", dpi=230, bbox_inches="tight")
            fig.savefig(dest / f"{stem}_side_by_side.pdf", bbox_inches="tight")
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)

        # mNeonGreen / GFP pages split by anti-GFP stain status (-2 vs -3)
        for metric, heading, unit, stem in GFP_PLOT_METRICS:
            fig, axes = plt.subplots(len(ANTI_GFP_ORDER), len(cell_lines), figsize=(5.3 * len(cell_lines), 3.85 * len(ANTI_GFP_ORDER)), squeeze=False, sharey=False)
            all_y = []
            for ri, stain_status in enumerate(ANTI_GFP_ORDER):
                for ci, cell_line in enumerate(cell_lines):
                    ax = axes[ri][ci]
                    sub = df[(df["cell_line"] == cell_line) & (df["anti_gfp_status"] == stain_status)].copy()
                    if sub.empty:
                        ax.set_axis_off()
                        continue
                    draw(ax, sub, metric, unit, palette, title=f"{cell_line} · {stain_status} (n={len(sub)} fields)")
                    vals = pd.to_numeric(sub[metric], errors="coerce")
                    vals = vals[np.isfinite(vals)]
                    if len(vals):
                        all_y.extend(vals.tolist())
                axes[ri][0].text(-0.55, 0.5, stain_status, transform=axes[ri][0].transAxes, rotation=90, va="center", ha="center", fontsize=10, weight="semibold")
            if all_y:
                ymin, ymax = min(all_y), max(all_y)
                if np.isfinite(ymin) and np.isfinite(ymax):
                    if ymin == ymax:
                        pad = abs(ymin) * 0.08 if ymin else 1.0
                    else:
                        pad = (ymax - ymin) * 0.08
                    for row in axes:
                        for ax in row:
                            if ax.axison:
                                ax.set_ylim(ymin - pad, ymax + pad)
            handles = [plt.Line2D([0], [0], marker="o", linestyle="none", markerfacecolor=palette[s], markeredgecolor="none", label=s) for s in samples]
            fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 0.145), ncol=min(6, len(handles)), frameon=False, title="Sample (block-section)")
            fig.suptitle(heading + " · IGF vs No IGF, split by anti-GFP staining", fontsize=13, weight="semibold")
            add_footer(fig, metric)
            fig.subplots_adjust(top=.90, bottom=.27, left=.08, wspace=.30, hspace=.42)
            fig.savefig(dest / f"{stem}_anti_gfp_split.png", dpi=230, bbox_inches="tight")
            fig.savefig(dest / f"{stem}_anti_gfp_split.pdf", bbox_inches="tight")
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--h1-workspace", type=Path, required=True, help="Original IGF workspace containing the verified H1 selection CSV and earlier fields")
    ap.add_argument("--remaining-workspace", type=Path, required=True, help="New completed IGF workspace containing the remaining finished fields")
    ap.add_argument("--sample-metadata", type=Path, help="CSV mapping Experiment 17 ventral block IDs -> cell_line and condition; defaults to the tracked CSV beside this script")
    ap.add_argument("--h1-selection", type=Path, help="Defaults to h1-workspace/mentor_field_plots/PV_density_per_field_verified_H1.csv")
    ap.add_argument("--out", type=Path, help="Output folder for combined results")
    args = ap.parse_args()

    h1_workspace = args.h1_workspace.expanduser().resolve()
    remaining_workspace = args.remaining_workspace.expanduser().resolve()
    sample_metadata = (args.sample_metadata or DEFAULT_SAMPLE_METADATA).expanduser().resolve()
    sample_info = load_sample_info(sample_metadata)
    h1_selection = (args.h1_selection or h1_workspace / "mentor_field_plots" / "PV_density_per_field_verified_H1.csv").expanduser().resolve()
    output = (args.out or h1_workspace / "mentor_field_plots" / "igf_exp17_ventral_only_side_by_side").expanduser().resolve()
    all_fields, all_pv, audit = collect_fields(h1_workspace, remaining_workspace, h1_selection, sample_info)

    frame = pd.DataFrame(all_fields).sort_values(["cell_line", "condition", "sample", "index"]).reset_index(drop=True)
    if frame["index"].duplicated().any():
        raise ValueError(f"Duplicate field indices in results: {frame.loc[frame['index'].duplicated(), 'index'].tolist()}")
    if set(frame["experiment_id"]) != {"17"} or set(frame["organoid_type"]) != {"Ventral"}:
        raise ValueError("Non-Experiment 17 or non-ventral field found in results")
    for line in CELL_LINE_ORDER:
        for condition in CONDITIONS:
            if not ((frame["cell_line"] == line) & (frame["condition"] == condition)).any():
                raise ValueError(f"Missing measured {line} / {condition} fields; cannot draw a valid IGF comparison")
    print("\nIncluded fields by cell line, condition, block and sample:")
    print(frame.groupby(["cell_line", "condition", "block_id", "sample"], observed=True).size().to_string())
    print("\nExcluded and incomplete indices:")
    print(pd.DataFrame(audit).loc[lambda d: ~d["included"], ["index", "sample", "reason"]].to_string(index=False))
    pv_frame = pd.concat(all_pv, ignore_index=True)
    output.mkdir(parents=True, exist_ok=True)
    frame.to_csv(output / "combined_field_cell_classes_pv_nuclear_area.csv", index=False)
    pv_frame.to_csv(output / "combined_pv_nuclear_positive_single_roi_areas.csv", index=False)

    summary = frame.groupby(["cell_line", "condition", "experiment_id", "sample"], observed=True).agg(
        n_fields=("index", "count"),
        mean_pv_within_lhx6_percent=("pv_within_lhx6_percent", "mean"),
        mean_lhx6_density_per_mm2=("lhx6_positive_density_per_mm2", "mean"),
        mean_double_positive_density_per_mm2=("lhx6_pv_double_positive_density_per_mm2", "mean"),
        mean_pv_area_percent_of_organoid=("pv_area_percent_of_organoid_nuclear_only", "mean"),
        mean_pv_area_per_cell_um2=("pv_area_mean_um2_per_positive_nucleus", "mean"),
    ).reset_index()
    summary.to_csv(output / "sample_level_descriptive_summary.csv", index=False)
    sample_metadata_out = pd.DataFrame([{"block_id": k, **v} for k, v in sorted(sample_info.items())])
    sample_metadata_out.to_csv(output / "sample_metadata_used.csv", index=False)
    pd.DataFrame(audit).to_csv(output / "field_selection_audit.csv", index=False)

    write_plots(frame, output)
    print(f"\nCompleted: {len(frame)} fields total; {len(pv_frame)} nuclear-positive PV objects")
    print(f"Outputs: {output}")
    print("IMPORTANT: condition labels without explicit IGF1 in the box inventory are inferred controls and require wet-lab confirmation.")
    print("Only Experiment 17 ventral block IDs 149,150,151,152,153,157 included. H1 uses previously verified fields.")


if __name__ == "__main__":
    try:
        main()
    except (ValueError, KeyError, FileNotFoundError) as exc:
        sys.exit(f"ERROR: {exc}")
