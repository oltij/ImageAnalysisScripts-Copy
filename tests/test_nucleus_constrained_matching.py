"""Synthetic regression tests for independent Stage-09 marker matching."""

from __future__ import annotations

import hashlib
import importlib.util
import io
import json
from pathlib import Path
import subprocess
import sys

import h5py
import numpy as np
import pandas as pd
import tifffile


ROOT = Path(__file__).resolve().parents[1]
SCRIPT1 = ROOT / "Colocalization/colocalizationdapiscript1.py"
COMPARE = ROOT / "Colocalization/compare_matching_variants.py"
SCRIPT2 = ROOT / "Colocalization/colocalizationdapiscript2.py"
RUNNER_SPEC = importlib.util.spec_from_file_location(
    "stage09_test_runner", ROOT / "run_stitched_pipeline.py"
)
runner = importlib.util.module_from_spec(RUNNER_SPEC)
assert RUNNER_SPEC.loader is not None
RUNNER_SPEC.loader.exec_module(runner)


def _store_table(h5: h5py.File, name: str, frame: pd.DataFrame) -> None:
    payload = np.frombuffer(frame.to_csv(index=False).encode("utf-8"), dtype=np.uint8)
    h5.create_dataset(f"tables/{name}", data=payload)


def _read_table(path: Path, name: str) -> pd.DataFrame:
    with h5py.File(path, "r") as h5:
        return pd.read_csv(io.BytesIO(h5[f"tables/{name}"][()].tobytes()))


def _pixels(object_id: int, coordinates: list[tuple[int, int]]) -> list[dict]:
    return [
        {"ROI": object_id, "ObjectID": object_id, "X": x, "Y": y}
        for x, y in coordinates
    ]


def _source_inputs():
    a_coordinates = {
        101: [(x, 1) for x in range(1, 6)],
        102: [(x, 1) for x in range(6, 11)],
        103: [(x, 2) for x in range(1, 6)],
        # Deliberately identical Stage-07 objects exercise fail-safe ambiguity.
        105: [(x, 4) for x in range(1, 6)],
        106: [(x, 4) for x in range(1, 6)],
    }
    b_coordinates = {
        201: [(x, 1) for x in range(1, 4)] + [(x, 3) for x in range(1, 6)],
        202: [(x, 1) for x in range(1, 9)],
        203: [(x, 2) for x in range(1, 3)],
        204: [(x, 4) for x in range(1, 6)],
    }
    a_mapping = [(11, 101, 9001), (12, 102, 9002), (13, 103, 9003),
                 (15, 105, 9004), (16, 106, 9005)]
    b_mapping = [(21, 201, 9001), (22, 202, 9002), (23, 203, 9003),
                 (24, 204, 9004)]
    return a_coordinates, b_coordinates, a_mapping, b_mapping


def _write_stage07(
    path: Path,
    marker: str,
    coordinates: dict[int, list[tuple[int, int]]],
    mapping_rows: list[tuple[int, int, int]],
    *,
    alter_nuclear_canvas: bool = False,
) -> None:
    marker_pixels = pd.DataFrame(
        [record for object_id, xy in coordinates.items() for record in _pixels(object_id, xy)]
    )
    mapping = pd.DataFrame(
        mapping_rows, columns=["CellID", "A_ObjectID", "B_ObjectID"]
    )
    cell_lookup = mapping.set_index("A_ObjectID")["CellID"].to_dict()
    filtered_rows = []
    for row in marker_pixels.itertuples(index=False):
        filtered_rows.append(
            {"CellID": cell_lookup[int(row.ObjectID)], "ROI": int(row.ROI),
             "X": int(row.X), "Y": int(row.Y)}
        )
    nuclear_rows = [
        {"ObjectID": nucleus, "X": index + 1, "Y": 6}
        for index, nucleus in enumerate((9001, 9002, 9003, 9004, 9005))
    ]
    if alter_nuclear_canvas:
        nuclear_rows[0]["X"] = 9
    metadata = {
        "objects": {"object_a_name": marker, "object_b_name": "Hoechst"},
        "pixel_size_um": {"x": 0.5, "y": 0.5},
        "matching": {"threshold": 0.5},
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(path, "w") as h5:
        h5.create_dataset("metadata/json", data=np.bytes_(json.dumps(metadata)))
        h5.create_dataset("images/object_a_original", data=np.zeros((8, 12), np.uint16))
        h5.create_dataset("images/object_b_original", data=np.zeros((8, 12), np.uint16))
        _store_table(h5, "object_a_pixels", marker_pixels)
        _store_table(h5, "object_b_pixels", pd.DataFrame(nuclear_rows))
        _store_table(h5, "cell_id_mapping", mapping)
        _store_table(h5, "filtered_object_a_pixels", pd.DataFrame(filtered_rows))


def _current_pixels(
    coordinates: dict[int, list[tuple[int, int]]],
    source_cells: dict[int, int],
) -> pd.DataFrame:
    rows = []
    for object_id, xy in coordinates.items():
        for x, y in xy:
            rows.append(
                {"ROI": object_id, "CellID": source_cells[object_id], "X": x, "Y": y}
            )
    return pd.DataFrame(rows)


def _run_analysis(
    work: Path,
    label: str,
    a_csv: Path,
    b_csv: Path,
    a_tiff: Path,
    b_tiff: Path,
    organoid_csv: Path,
    *,
    constrained: bool,
    a_stage07: Path | None = None,
    b_stage07: Path | None = None,
) -> tuple[Path, Path | None, subprocess.CompletedProcess[str]]:
    result = work / label / "analysis.h5"
    audit = work / label / "nuclear_constraint_pair_audit.csv" if constrained else None
    replacements = {
        "OBJECT_A_NAME": repr("LHX6"),
        "OBJECT_B_NAME": repr("PV"),
        "OBJECT_A_CSV": runner.path_value(a_csv),
        "OBJECT_B_CSV": runner.path_value(b_csv),
        "OBJECT_A_TIFF": runner.path_value(a_tiff),
        "OBJECT_B_TIFF": runner.path_value(b_tiff),
        "ORGANOID_CSV": runner.path_value(organoid_csv),
        "RESULTS_FILE": runner.path_value(result),
        "THRESHOLD": "0.5",
        "PIXEL_SIZE_X_UM": "0.5",
        "PIXEL_SIZE_Y_UM": "0.5",
    }
    if constrained:
        replacements.update(
            {
                "NUCLEAR_IDENTITY_CONSTRAINT": "True",
                "OBJECT_A_NUCLEAR_RESULTS_FILE": runner.path_value(a_stage07),
                "OBJECT_B_NUCLEAR_RESULTS_FILE": runner.path_value(b_stage07),
                "NUCLEAR_CONSTRAINT_AUDIT_CSV": runner.path_value(audit),
            }
        )
    patched = work / f"{label}.py"
    patched.write_bytes(runner.patched_source(SCRIPT1, replacements))
    completed = subprocess.run(
        [sys.executable, str(patched)], text=True, capture_output=True, check=False
    )
    return result, audit, completed


def _run_visualization(
    work: Path, h5_path: Path, paths: dict
) -> tuple[Path, subprocess.CompletedProcess[str]]:
    output = work / "visualizations"
    replacements = {
        "RESULTS_FILE": runner.path_value(h5_path),
        "OUTPUT_DIR": runner.path_value(output),
        "ENABLE_ORGANOID_AREA_SUMMARY": "True",
        "ORGANOID_MEASUREMENTS_CSV": runner.path_value(paths["organoid"]),
        "ENABLE_NUCLEAR_ALL_CELL_SUMMARY": "True",
        "OBJECT_A_NUCLEAR_RESULTS_FILE": runner.path_value(paths["a_stage07"]),
        "OBJECT_B_NUCLEAR_RESULTS_FILE": runner.path_value(paths["b_stage07"]),
    }
    patched = work / "stage09_visualization.py"
    patched.write_bytes(runner.patched_source(SCRIPT2, replacements))
    completed = subprocess.run(
        [sys.executable, str(patched)], text=True, capture_output=True, check=False
    )
    return output, completed


def _prepare_fixture(work: Path, *, incompatible_b_nuclei: bool = False) -> dict:
    a_coordinates, b_coordinates, a_mapping, b_mapping = _source_inputs()
    a_stage07 = work / "stage07_lhx6.h5"
    b_stage07 = work / "stage07_pv.h5"
    _write_stage07(a_stage07, "LHX6", a_coordinates, a_mapping)
    _write_stage07(
        b_stage07,
        "PV",
        b_coordinates,
        b_mapping,
        alter_nuclear_canvas=incompatible_b_nuclei,
    )

    # ROI 104 has no Stage-07 record; ROI 105 has two coordinate-identical
    # Stage-07 sources. Both must remain unresolved, never matched by ID alone.
    current_a_coordinates = dict(a_coordinates)
    current_a_coordinates.pop(106)
    current_a_coordinates[104] = [(x, 3) for x in range(1, 6)]
    a_csv = work / "lhx6_nuclear_positive.csv"
    b_csv = work / "pv_nuclear_positive.csv"
    _current_pixels(
        current_a_coordinates,
        {101: 11, 102: 12, 103: 13, 104: 999, 105: 15},
    ).to_csv(a_csv, index=False)
    _current_pixels(
        b_coordinates, {201: 21, 202: 22, 203: 23, 204: 24}
    ).to_csv(b_csv, index=False)
    a_tiff = work / "lhx6.tif"
    b_tiff = work / "pv.tif"
    tifffile.imwrite(a_tiff, np.arange(96, dtype=np.uint16).reshape(8, 12))
    tifffile.imwrite(b_tiff, np.arange(96, dtype=np.uint16).reshape(8, 12) * 2)
    organoid = work / "organoid.csv"
    pd.DataFrame({"AreaShape_Area": [80]}).to_csv(organoid, index=False)
    return {
        "a_csv": a_csv, "b_csv": b_csv, "a_tiff": a_tiff, "b_tiff": b_tiff,
        "organoid": organoid, "a_stage07": a_stage07, "b_stage07": b_stage07,
    }


def _pair_set(frame: pd.DataFrame) -> set[tuple[int, int]]:
    return set(zip(frame["A_ObjectID"].astype(int), frame["B_ObjectID"].astype(int)))


def test_stage09_restricts_before_matching_and_records_provenance(tmp_path: Path):
    paths = _prepare_fixture(tmp_path)
    stage08, _, run08 = _run_analysis(tmp_path, "stage08", **paths_without_sources(paths), constrained=False)
    assert run08.returncode == 0, run08.stdout + run08.stderr
    before_hash = hashlib.sha256(stage08.read_bytes()).hexdigest()

    stage09, audit_csv, run09 = _run_analysis(
        tmp_path,
        "stage09",
        **paths_without_sources(paths),
        constrained=True,
        a_stage07=paths["a_stage07"],
        b_stage07=paths["b_stage07"],
    )
    assert run09.returncode == 0, run09.stdout + run09.stderr
    assert audit_csv is not None and audit_csv.is_file()
    assert hashlib.sha256(stage08.read_bytes()).hexdigest() == before_hash

    unconstrained = _read_table(stage08, "stable_matching")
    constrained = _read_table(stage09, "stable_matching")
    assert _pair_set(unconstrained) == {(101, 202), (104, 201), (105, 204)}
    assert _pair_set(constrained) == {(101, 201), (102, 202)}
    # This proves rematching was recomputed; neither result is a filtered
    # subset of the other.
    assert (101, 201) not in _pair_set(unconstrained)

    with h5py.File(stage08, "r") as h5:
        stage08_meta = json.loads(h5["metadata/json"][()].decode("utf-8"))
        assert "analysis_variant" not in stage08_meta
        assert "tables/nuclear_constraint_audit" not in h5
    with h5py.File(stage09, "r") as h5:
        stage09_meta = json.loads(h5["metadata/json"][()].decode("utf-8"))
        assert stage09_meta["analysis_variant"] == "nucleus_constrained_marker_marker"
        assert stage09_meta["nuclear_identity_constraint"]["applied_before_preference_lists"]
        assert Path(stage09_meta["nuclear_identity_constraint"]["object_a_stage07_results_file"]) == paths["a_stage07"]
        assert "tables/nuclear_assignments_object_a" in h5

    audit = pd.read_csv(audit_csv).set_index(["A_ObjectID", "B_ObjectID"])
    same = audit.loc[(101, 201)]
    assert bool(same["Original_Overlap_Threshold_Pass"])
    assert bool(same["Nuclear_Constraint_Pass"])
    assert bool(same["Final_Constrained_Stable_Match"])
    assert int(same["A_Stage07_Marker_ObjectID"]) == 101
    assert int(same["A_Stage07_Source_CellID"]) == 11
    assert int(same["A_Canonical_Nuclear_ObjectID"]) == 9001
    different = audit.loc[(101, 202)]
    assert different["Exclusion_Reason"] == "different_canonical_nuclear_objectids"
    below = audit.loc[(103, 203)]
    assert below["Exclusion_Reason"] == "below_directional_overlap_threshold"
    missing = audit.loc[(104, 201)]
    assert "missing_coordinate_identity" in missing["Exclusion_Reason"]
    assert int(missing["B_Canonical_Nuclear_ObjectID"]) == 9001
    ambiguous = audit.loc[(105, 204)]
    assert "ambiguous_coordinate_identity" in ambiguous["Exclusion_Reason"]

    comparison_dir = tmp_path / "comparison"
    comparison = subprocess.run(
        [sys.executable, str(COMPARE), "--unconstrained-h5", str(stage08),
         "--constrained-h5", str(stage09), "--output-dir", str(comparison_dir)],
        text=True, capture_output=True, check=False,
    )
    assert comparison.returncode == 0, comparison.stdout + comparison.stderr
    for filename in (
        "matching_variant_comparison.csv", "matching_variant_pair_details.csv",
        "discordant_unconstrained_final_pairs.csv", "matching_variant_comparison.pdf",
        "matching_variant_comparison.png",
    ):
        assert (comparison_dir / filename).is_file()
        assert (comparison_dir / filename).stat().st_size > 0
    details = pd.read_csv(comparison_dir / "matching_variant_pair_details.csv")
    discordant = details[details["Nuclear_Identity_Discordant"]]
    assert _pair_set(discordant) == {(101, 202)}


def paths_without_sources(paths: dict) -> dict:
    return {
        "a_csv": paths["a_csv"],
        "b_csv": paths["b_csv"],
        "a_tiff": paths["a_tiff"],
        "b_tiff": paths["b_tiff"],
        "organoid_csv": paths["organoid"],
    }


def test_zero_eligible_constrained_result_is_valid(tmp_path: Path):
    paths = _prepare_fixture(tmp_path)
    a = pd.read_csv(paths["a_csv"])
    b = pd.read_csv(paths["b_csv"])
    a[a["ROI"] == 103].to_csv(paths["a_csv"], index=False)
    b[b["ROI"] == 203].to_csv(paths["b_csv"], index=False)
    stage09, audit_csv, completed = _run_analysis(
        tmp_path, "stage09_zero", **paths_without_sources(paths), constrained=True,
        a_stage07=paths["a_stage07"], b_stage07=paths["b_stage07"],
    )
    assert completed.returncode == 0, completed.stdout + completed.stderr
    assert _read_table(stage09, "stable_matching").empty
    audit = pd.read_csv(audit_csv)
    assert len(audit) == 1
    assert audit.iloc[0]["Exclusion_Reason"] == "below_directional_overlap_threshold"
    visualization, visualized = _run_visualization(tmp_path, stage09, paths)
    assert visualized.returncode == 0, visualized.stdout + visualized.stderr
    assert (visualization / "LHX6_PV_matched_pair_visualization.pdf").stat().st_size > 0
    assert (visualization / "LHX6_PV_filtering_stage_counts.csv").stat().st_size > 0
    all_nuclei = (
        visualization
        / "unmatched_only_qc"
        / "LHX6_PV_all_nuclear_cells_classification.csv"
    )
    assert all_nuclei.stat().st_size > 0


def test_incompatible_stage07_nuclear_segmentation_fails_safely(tmp_path: Path):
    paths = _prepare_fixture(tmp_path, incompatible_b_nuclei=True)
    stage09, audit_csv, completed = _run_analysis(
        tmp_path, "stage09_incompatible", **paths_without_sources(paths), constrained=True,
        a_stage07=paths["a_stage07"], b_stage07=paths["b_stage07"],
    )
    assert completed.returncode != 0
    assert "same canonical nuclear segmentation" in (completed.stdout + completed.stderr)
    assert not stage09.exists()
    assert audit_csv is not None and not audit_csv.exists()
