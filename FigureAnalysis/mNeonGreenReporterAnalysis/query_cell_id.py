#!/usr/bin/env python3
"""Query one Stage-07 marker-specific CellID across all available markers and Stages 08/09.

IDs are scoped to a FIELD and a STAGE-07 MARKER-NUCLEAR RUN. They are not
cross-marker identities: canonical shared identity is the Hoechst nuclear ObjectID.
Read-only relative to source analysis. Writes a new CSV/JSON/PNG/PDF report.

Requires numpy, pandas, h5py, matplotlib. No imports from other repo scripts.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
import io
import json
from pathlib import Path
import re
import sys

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd


@dataclass
class MarkerSource:
    name: str
    nucleus_name: str
    file: Path
    mapping: pd.DataFrame
    pixels: dict[int, np.ndarray] | None = None
    nuclear_pixels: dict[int, np.ndarray] | None = None


def table(h5: h5py.File, name: str) -> pd.DataFrame:
    key = f"tables/{name}"
    if key not in h5:
        raise ValueError(f"{h5.filename}: missing {key}")
    return pd.read_csv(io.BytesIO(h5[key][()].tobytes()))


def pair_token(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", str(value)).strip("_")


def optional_int(value):
    return None if pd.isna(value) else int(value)


def xy_groups(df: pd.DataFrame, *, label: str) -> dict[int, np.ndarray]:
    if not {"ObjectID", "X", "Y"}.issubset(df.columns):
        raise ValueError(f"{label}: missing ObjectID/X/Y")
    out = {}
    for oid, group in df.groupby("ObjectID", sort=False):
        points = group[["X", "Y"]].to_numpy(dtype=np.int64)
        if len(np.unique(points, axis=0)) != len(points):
            raise ValueError(f"{label}: ROI {oid} has duplicated XY pixels")
        out[int(oid)] = points
    return out


def same_points(a: np.ndarray, b: np.ndarray) -> bool:
    if len(a) != len(b):
        return False
    return np.array_equal(a[np.lexsort((a[:, 1], a[:, 0]))],
                          b[np.lexsort((b[:, 1], b[:, 0]))])


def collect_sources(analysis: Path) -> dict[str, MarkerSource]:
    found = {}
    for path in sorted((analysis / "07_nuclear_colocalization").glob("*/analysis.h5")):
        with h5py.File(path, "r") as h5:
            meta = json.loads(h5["metadata/json"][()].decode())
            a = meta["objects"]["object_a_name"]
            b = meta["objects"]["object_b_name"]
            mapping = table(h5, "cell_id_mapping")
            needed = {"CellID", "A_ObjectID", "B_ObjectID"}
            if not needed.issubset(mapping.columns):
                raise ValueError(f"{path}: incomplete Stage-07 mapping")
            mapping = mapping[list(("CellID", "A_ObjectID", "B_ObjectID"))].copy()
            for c in needed:
                mapping[c] = pd.to_numeric(mapping[c], errors="raise").astype("int64")
                if mapping[c].duplicated().any():
                    raise ValueError(f"{path}: duplicate {c} in one-to-one Stage-07 mapping")
        if a in found:
            raise ValueError(f"Ambiguous Stage-07 sources for {a}: {found[a].file} and {path}")
        found[a] = MarkerSource(a, b, path, mapping)
    if not found:
        raise FileNotFoundError(f"No Stage-07 HDF5 files in {analysis}/07_nuclear_colocalization")
    nuclear_names = {src.nucleus_name for src in found.values()}
    if len(nuclear_names) != 1:
        raise ValueError(f"Stage-07 sources have incompatible nuclear names: {sorted(nuclear_names)}")
    # Reusing a numeric nuclear ObjectID across independently segmented images
    # would produce a false same-cell conclusion. Require identical source
    # nuclear segmentation and canvas, just as Stage 09 does.
    reference = None
    for src in found.values():
        with h5py.File(src.file, "r") as h5:
            size = tuple(h5["images/object_b_original"].shape)
            pixels = table(h5, "object_b_pixels")[["ObjectID", "X", "Y"]].copy()
            for col in ("ObjectID", "X", "Y"):
                pixels[col] = pd.to_numeric(pixels[col], errors="raise").astype("int64")
            pixels = pixels.drop_duplicates().sort_values(["ObjectID", "Y", "X"]).reset_index(drop=True)
            metadata = json.loads(h5["metadata/json"][()].decode())
            cal = metadata.get("pixel_size_um")
        if reference is None:
            reference = (size, pixels, cal)
        elif size != reference[0] or not pixels.equals(reference[1]) or cal != reference[2]:
            raise ValueError("Stage-07 sources do not share the same canonical nuclear segmentation, image size and calibration")
    return found


def select_nucleus(sources: dict[str, MarkerSource], cell_id: int | None,
                   marker: str | None, nuclear_id: int | None) -> tuple[int, str]:
    if nuclear_id is not None:
        found = any(int(nuclear_id) in set(src.mapping.B_ObjectID.astype(int)) for src in sources.values())
        if not found:
            # Negative results are ambiguous: may be a real nucleus with no matched markers.
            raise ValueError(f"Nuclear ObjectID {nuclear_id} has no Stage-07 positive marker assignment; use a matched CellID instead")
        return nuclear_id, f"{next(iter(sources.values())).nucleus_name} nucleus {nuclear_id}"
    if cell_id is None:
        raise ValueError("Supply --cell-id or --nuclear-id")
    if marker is not None and marker not in sources:
        raise ValueError(f"Marker {marker!r} has no Stage-07 source. Available: {', '.join(sources)}")
    matched = []
    for name, src in sources.items():
        if marker and name != marker:
            continue
        row = src.mapping[src.mapping.CellID == cell_id]
        if len(row):
            matched.append((name, int(row.iloc[0].B_ObjectID)))
    if not matched:
        raise ValueError(f"Stage-07 CellID {cell_id} was not found" + (f" for {marker}" if marker else ""))
    if len(matched) != 1:
        options = ", ".join(f"{name} CellID {cell_id} -> nucleus {n}" for name, n in matched)
        raise ValueError("CellID is ambiguous across marker runs; specify --marker. Candidates: " + options)
    name, nucleus = matched[0]
    return nucleus, f"{name} Stage-07 CellID {cell_id}"


def get_pixels(source: MarkerSource, *, nuclear=False) -> dict[int, np.ndarray]:
    attribute = "nuclear_pixels" if nuclear else "pixels"
    value = getattr(source, attribute)
    if value is None:
        with h5py.File(source.file, "r") as h5:
            value = xy_groups(table(h5, "object_b_pixels" if nuclear else "object_a_pixels"), label=str(source.file))
        setattr(source, attribute, value)
    return value


def crop_bounds(xy: np.ndarray, shape: tuple[int, int], padding: int = 24) -> tuple[int, int, int, int]:
    h, w = shape
    return (max(0, int(xy[:, 0].min()) - padding), min(w, int(xy[:, 0].max()) + padding + 1),
            max(0, int(xy[:, 1].min()) - padding), min(h, int(xy[:, 1].max()) + padding + 1))


def values_for_roi(img: h5py.Dataset, xy: np.ndarray) -> np.ndarray:
    if not len(xy):
        raise ValueError("ROI cannot contain zero pixels")
    h, w = img.shape
    if np.any(xy[:, 0] < 0) or np.any(xy[:, 0] >= w) or np.any(xy[:, 1] < 0) or np.any(xy[:, 1] >= h):
        raise ValueError(f"ROI outside image {img.shape}")
    x0, x1, y0, y1 = crop_bounds(xy, img.shape, 0)
    crop = img[y0:y1, x0:x1]
    return crop[xy[:, 1] - y0, xy[:, 0] - x0].astype(np.float64)


def norm(arr: np.ndarray) -> np.ndarray:
    a = np.asarray(arr, dtype=np.float64)
    values = a[np.isfinite(a)]
    if len(values) == 0:
        return np.zeros(a.shape, dtype=np.float32)
    lo, hi = np.percentile(values, [1, 99])
    return np.clip((a - lo) / max(float(hi - lo), 1e-10), 0, 1).astype(np.float32)


def outline(xy: np.ndarray, shape: tuple[int,int], x0: int, y0: int) -> np.ndarray:
    mask = np.zeros(shape, dtype=bool)
    mask[xy[:, 1] - y0, xy[:, 0] - x0] = True
    inner = mask.copy()
    inner[1:-1, 1:-1] &= mask[:-2, 1:-1] & mask[2:, 1:-1] & mask[1:-1, :-2] & mask[1:-1, 2:]
    return mask & ~inner


def source_row_for_nucleus(src: MarkerSource, nuclear_id: int) -> pd.Series | None:
    selected = src.mapping[src.mapping.B_ObjectID == nuclear_id]
    if len(selected) > 1:
        raise ValueError(f"{src.file}: more than one marker ROI assigned to nucleus {nuclear_id}")
    return selected.iloc[0] if len(selected) else None


def measure_markers(analysis: Path, sources: dict[str, MarkerSource], nucleus: int):
    records, crops = [], []
    for name, src in sources.items():
        r = source_row_for_nucleus(src, nucleus)
        rec = {"Field_Analysis_Path": str(analysis), "Marker": name,
               "Canonical_Hoechst_Nuclear_ObjectID": nucleus,
               "Stage07_Source_H5": str(src.file), "Stage07_Assigned": r is not None,
               "Stage07_CellID": None, "Stage07_Cell_Key": None,
               "Original_Marker_ObjectID": None}
        if r is None:
            records.append(rec)
            continue
        roi_id, cell_id = int(r.A_ObjectID), int(r.CellID)
        rec.update(Stage07_CellID=cell_id, Stage07_Cell_Key=f"{analysis}::{name}::C{cell_id}",
                   Original_Marker_ObjectID=roi_id)
        coords = get_pixels(src).get(roi_id)
        nuc_coords = get_pixels(src, nuclear=True).get(nucleus)
        if coords is None or nuc_coords is None:
            raise ValueError(f"{src.file}: missing matched source marker or nuclear pixels")
        with h5py.File(src.file, "r") as h5:
            a = h5["images/object_a_original"]
            b = h5["images/object_b_original"]
            if a.shape != b.shape or len(a.shape) != 2:
                raise ValueError(f"Incompatible source images in {src.file}")
            vals_a = values_for_roi(a, coords)
            vals_b = values_for_roi(b, nuc_coords)
            combined = np.vstack((coords, nuc_coords))
            x0, x1, y0, y1 = crop_bounds(combined, a.shape)
            ca = np.asarray(a[y0:y1, x0:x1]); cb = np.asarray(b[y0:y1, x0:x1])
        flat_marker = set((int(x), int(y)) for x, y in coords)
        flat_nuclear = set((int(x), int(y)) for x, y in nuc_coords)
        rec.update(Marker_Area_Pixels=len(coords), Marker_Mean_Native=float(vals_a.mean()),
                   Marker_Integrated_Native=float(vals_a.sum()),
                   Hoechst_Area_Pixels=len(nuc_coords), Hoechst_Mean_Native=float(vals_b.mean()),
                   Hoechst_Integrated_Native=float(vals_b.sum()),
                   Marker_Nuclear_Shared_Pixels=len(flat_marker & flat_nuclear))
        metrics_file = analysis / "06_final_shape_filtered" / name / f"{name}_kept_ROI_metrics.csv"
        if metrics_file.is_file():
            metrics = pd.read_csv(metrics_file)
            if "ObjectID" not in metrics:
                raise ValueError(f"{metrics_file} lacks ObjectID")
            matched = metrics[metrics["ObjectID"] == roi_id]
            if len(matched) != 1:
                raise ValueError(f"{metrics_file}: expected exactly one measurement row for ROI {roi_id}, found {len(matched)}")
            rec.update({f"FinalROI_{k}": v for k, v in matched.iloc[0].items()})
        records.append(rec)
        crops.append(dict(marker=name, cell_id=cell_id, roi_id=roi_id, ca=ca, cb=cb,
                          xy=coords, nxy=nuc_coords, x0=x0, y0=y0))
    return pd.DataFrame(records), crops


def resolve_pair_roi(raw: pd.DataFrame, source: MarkerSource) -> tuple[int | None, int | None, int | None, str]:
    """Get Stage-07 (original ROI, source CellID, nuclear ID) using XY as identity.

    Caller passes exactly one current ObjectID's full pixel table. Never infer
    same-cell status from numeric ID alone. Unresolved identities are explicit.
    """
    if raw.empty:
        return None, None, None, "missing_current_roi"
    points = raw[["X", "Y"]].to_numpy(dtype=np.int64)
    hints = set()
    if "ROI" in raw:
        ids = pd.to_numeric(raw.ROI, errors="raise").unique()
        if len(ids) != 1:
            return None, None, None, "multiple_original_ROI_ids"
        hints.add(("ROI", int(ids[0])))
    if "CellID" in raw:
        ids = pd.to_numeric(raw.CellID, errors="raise").unique()
        if len(ids) != 1:
            return None, None, None, "multiple_source_CellIDs"
        hints.add(("CellID", int(ids[0])))
    proposals = []
    for col, value in hints:
        q = source.mapping[source.mapping.A_ObjectID.eq(value) if col == "ROI" else source.mapping.CellID.eq(value)]
        if len(q) == 1:
            proposals.append(int(q.iloc[0].A_ObjectID))
    # If numeric ID hints absent, coordinate matching is possible but expensive;
    # query only candidates mapped to the currently selected nucleus is not enough
    # for determining a partner's independent nucleus, so fail safely.
    if not proposals:
        return None, None, None, "no_authoritative_source_id_hint"
    if len(set(proposals)) != 1:
        return None, None, None, "conflicting_source_id_hints"
    oid = proposals[0]
    ref = get_pixels(source).get(oid)
    if ref is None or not same_points(points, ref):
        return None, None, None, "source_pixels_mismatch"
    row = source.mapping[source.mapping.A_ObjectID == oid].iloc[0]
    return oid, int(row.CellID), int(row.B_ObjectID), "verified"


def pair_mapping_from_h5(path: Path, stage: int, sources: dict[str, MarkerSource]) -> pd.DataFrame:
    """Read pair provenance, reconstructing deterministic IDs for legacy HDF5."""
    with h5py.File(path, "r") as h5:
        meta = json.loads(h5["metadata/json"][()].decode())
        a, b = meta["objects"]["object_a_name"], meta["objects"]["object_b_name"]
        matching = table(h5, "stable_matching").sort_values(
            ["A_ObjectID", "B_ObjectID"], kind="mergesort"
        ).reset_index(drop=True)
        if "tables/pair_id_mapping" in h5:
            mapping = table(h5, "pair_id_mapping")
            required = {"PairID", "A_ObjectID", "B_ObjectID", "Same_Nucleus"}
            if not required.issubset(mapping.columns):
                raise ValueError(f"{path}: pair_id_mapping missing {sorted(required - set(mapping.columns))}")
            if set(zip(mapping.A_ObjectID.astype(int), mapping.B_ObjectID.astype(int))) != set(
                zip(matching.A_ObjectID.astype(int), matching.B_ObjectID.astype(int))
            ):
                raise ValueError(f"{path}: pair_id_mapping does not match final stable pairs")
            expected_ids = [
                f"S{stage:02d}_{pair_token(a)}_{pair_token(b)}_P{i:06d}"
                for i in range(1, len(matching) + 1)
            ]
            actual = mapping.sort_values(["A_ObjectID", "B_ObjectID"], kind="mergesort")
            if actual.PairID.astype(str).tolist() != expected_ids:
                raise ValueError(f"{path}: PairIDs are not deterministic for sorted final pairs")
            if stage == 9 and len(mapping) and not mapping.Same_Nucleus.astype(str).eq("YES").all():
                raise ValueError(f"{path}: Stage-09 final PairID is not SAME NUCLEUS: YES")
            if "A_Marker" not in mapping:
                mapping["A_Marker"] = a
            if "B_Marker" not in mapping:
                mapping["B_Marker"] = b
            return mapping
        if a not in sources or b not in sources:
            raise ValueError(f"{path}: cannot reconstruct provenance without Stage-07 sources")
        a_pixels = table(h5, "object_a_pixels")
        b_pixels = table(h5, "object_b_pixels")
        overlaps = table(h5, "overlap_results").set_index(["A_ObjectID", "B_ObjectID"])
    rows = []
    for sequence, match in enumerate(matching.itertuples(index=False), 1):
        aid, bid = int(match.A_ObjectID), int(match.B_ObjectID)
        ar, ac, an, ast = resolve_pair_roi(a_pixels[a_pixels.ObjectID == aid], sources[a])
        br, bc, bn, bst = resolve_pair_roi(b_pixels[b_pixels.ObjectID == bid], sources[b])
        same = "YES" if ast == bst == "verified" and an == bn else (
            "NO" if ast == bst == "verified" else "UNRESOLVED"
        )
        overlap = overlaps.loc[(aid, bid)]
        rows.append({
            "PairID": f"S{stage:02d}_{pair_token(a)}_{pair_token(b)}_P{sequence:06d}",
            "Analysis_Variant": "stage08_unconstrained" if stage == 8 else "stage09_nucleus_constrained",
            "A_Marker": a, "B_Marker": b,
            "A_ObjectID": aid, "B_ObjectID": bid,
            "A_Original_ROI_ObjectID": ar, "B_Original_ROI_ObjectID": br,
            "A_Stage07_Source_CellID": ac, "B_Stage07_Source_CellID": bc,
            "A_Canonical_Nuclear_ObjectID": an, "B_Canonical_Nuclear_ObjectID": bn,
            "A_Mapping_Status": ast, "B_Mapping_Status": bst,
            "Same_Nucleus": same,
            "Shared_Area": int(overlap.Shared_Area),
            "A_Percent_Inside_B": float(overlap.A_Percent_Inside_B),
            "B_Percent_Inside_A": float(overlap.B_Percent_Inside_A),
            "Final_Stable_Match": True,
        })
    result = pd.DataFrame(rows)
    if stage == 9 and len(result) and not result.Same_Nucleus.eq("YES").all():
        raise ValueError(f"{path}: Stage-09 final pair has unresolved/different nuclear identity")
    return result


def collect_pair_mappings(analysis: Path, sources: dict[str, MarkerSource]) -> pd.DataFrame:
    records = []
    roots = {
        8: "08_marker_colocalization",
        9: "09_marker_colocalization_nucleus_constrained",
    }
    for stage, root in roots.items():
        for path in sorted((analysis / root).glob("*/analysis.h5")):
            mapping = pair_mapping_from_h5(path, stage, sources).copy()
            mapping["Stage"] = stage
            mapping["Source_H5"] = str(path)
            records.append(mapping)
    if not records:
        return pd.DataFrame(columns=["PairID", "Stage", "Source_H5"])
    result = pd.concat(records, ignore_index=True)
    if result.PairID.duplicated().any():
        duplicates = result.loc[result.PairID.duplicated(False), "PairID"].tolist()
        raise ValueError(f"PairID is ambiguous within this field: {duplicates}")
    return result


def pair_measurement(path: Path, row: pd.Series) -> tuple[dict, dict]:
    aid, bid = int(row.A_ObjectID), int(row.B_ObjectID)
    with h5py.File(path, "r") as h5:
        meta = json.loads(h5["metadata/json"][()].decode())
        a, b = meta["objects"]["object_a_name"], meta["objects"]["object_b_name"]
        a_df, b_df = table(h5, "object_a_pixels"), table(h5, "object_b_pixels")
        aa = a_df[a_df.ObjectID == aid][["X", "Y"]].to_numpy(dtype=np.int64)
        bb = b_df[b_df.ObjectID == bid][["X", "Y"]].to_numpy(dtype=np.int64)
        if not len(aa) or not len(bb):
            raise ValueError(f"{path}: PairID {row.PairID} lacks ROI pixels")
        im_a, im_b = h5["images/object_a_original"], h5["images/object_b_original"]
        av, bv = values_for_roi(im_a, aa), values_for_roi(im_b, bb)
        a_flat = aa[:, 1] * im_a.shape[1] + aa[:, 0]
        b_flat = bb[:, 1] * im_a.shape[1] + bb[:, 0]
        common = np.intersect1d(a_flat, b_flat)
        cy, cx = np.divmod(common, im_a.shape[1])
        common_xy = np.column_stack((cx, cy)).astype(np.int64)
        a_overlap = values_for_roi(im_a, common_xy) if len(common_xy) else np.array([], dtype=float)
        b_overlap = values_for_roi(im_b, common_xy) if len(common_xy) else np.array([], dtype=float)
    record = row.to_dict()
    record.update({
        "A_Marker": a, "B_Marker": b,
        "A_ROI_Area_Pixels": len(aa), "B_ROI_Area_Pixels": len(bb),
        "Overlap_Area_Pixels": len(common),
        "A_Mean_Fluorescence_Native": float(av.mean()),
        "B_Mean_Fluorescence_Native": float(bv.mean()),
        "A_Integrated_Fluorescence_Native": float(av.sum()),
        "B_Integrated_Fluorescence_Native": float(bv.sum()),
        "A_Overlap_Integrated_Fluorescence_Native": float(a_overlap.sum()),
        "B_Overlap_Integrated_Fluorescence_Native": float(b_overlap.sum()),
    })
    view = {
        **record,
        "file": path,
        "A_XY": aa,
        "B_XY": bb,
        "Pair": f"{a}_vs_{b}",
        "A_Stage07_CellID": row.get("A_Stage07_Source_CellID"),
        "B_Stage07_CellID": row.get("B_Stage07_Source_CellID"),
        "A_Original_Marker_ObjectID": row.get("A_Original_ROI_ObjectID", aid),
        "B_Original_Marker_ObjectID": row.get("B_Original_ROI_ObjectID", bid),
        "Same_Hoechst_Nucleus": {"YES": True, "NO": False}.get(str(row.Same_Nucleus)),
        "Shared_Area_Pixels": len(common),
        "Stage_Local_Pairwise_CellID": row.get("Stage_Local_Pairwise_CellID"),
    }
    return record, view


def analyze_variant(analysis: Path, stage: int, sources: dict[str, MarkerSource], nucleus: int):
    root = ("08_marker_colocalization" if stage == 8 else "09_marker_colocalization_nucleus_constrained")
    files = sorted((analysis / root).glob("*/analysis.h5"))
    details, summaries, views = [], [], []
    for path in files:
        with h5py.File(path, "r") as h5:
            meta = json.loads(h5["metadata/json"][()].decode())
            a, b = meta["objects"]["object_a_name"], meta["objects"]["object_b_name"]
            if a not in sources or b not in sources:
                raise ValueError(f"{path}: pair uses marker without Stage-07 identity: {a}, {b}")
            if stage == 9 and meta.get("analysis_variant") != "nucleus_constrained_marker_marker":
                raise ValueError(f"{path}: marked Stage 09 but missing constrained metadata")
            a_df = table(h5, "object_a_pixels"); b_df = table(h5, "object_b_pixels")
            overlaps = table(h5, "overlap_results")
            matching = table(h5, "stable_matching")
            maps = table(h5, "cell_id_mapping")
            if "tables/pair_id_mapping" in h5:
                pair_ids = {
                    (int(r.A_ObjectID), int(r.B_ObjectID)): str(r.PairID)
                    for r in table(h5, "pair_id_mapping").itertuples(index=False)
                }
            else:
                ordered_matches = matching.sort_values(
                    ["A_ObjectID", "B_ObjectID"], kind="mergesort"
                ).reset_index(drop=True)
                pair_ids = {
                    (int(r.A_ObjectID), int(r.B_ObjectID)): f"S{stage:02d}_{pair_token(a)}_{pair_token(b)}_P{i:06d}"
                    for i, r in enumerate(ordered_matches.itertuples(index=False), 1)
                }
            if not {"A_ObjectID", "B_ObjectID"}.issubset(overlaps) or not {"A_ObjectID", "B_ObjectID"}.issubset(matching):
                raise ValueError(f"{path}: pair HDF5 missing ID columns")
            audit = table(h5, "nuclear_constraint_audit") if stage == 9 else None
            if stage == 9 and set(zip(audit.A_ObjectID.astype(int), audit.B_ObjectID.astype(int))) != set(zip(overlaps.A_ObjectID.astype(int), overlaps.B_ObjectID.astype(int))):
                raise ValueError(f"{path}: Stage-09 audit does not cover precisely the overlap candidates")
            # Find only the requested nucleus' ROIs first. Resolve other ROIs
            # only when they occur in a spatial candidate involving that nucleus.
            idmaps: dict[str, dict[int, tuple]] = {"A": {}, "B": {}}
            ids: dict[str, set[int]] = {"A": set(), "B": set()}
            pair_pixel_tables = {"A": a_df, "B": b_df}
            for side, name in (("A", a), ("B", b)):
                source_row = source_row_for_nucleus(sources[name], nucleus)
                if source_row is None:
                    continue
                original, source_cell = int(source_row.A_ObjectID), int(source_row.CellID)
                px = pair_pixel_tables[side]
                possible = pd.Series(False, index=px.index)
                if "ROI" in px:
                    possible |= pd.to_numeric(px.ROI, errors="coerce").eq(original)
                if "CellID" in px:
                    possible |= pd.to_numeric(px.CellID, errors="coerce").eq(source_cell)
                # Handle CSVs originally identified only by ROI or CellID.
                possible |= px.ObjectID.eq(original) | px.ObjectID.eq(source_cell)
                candidate_ids = set(px.loc[possible, "ObjectID"].astype(int))
                for oid, group in px[px.ObjectID.isin(candidate_ids)].groupby("ObjectID", sort=False):
                    resolved = resolve_pair_roi(group, sources[name])
                    idmaps[side][int(oid)] = resolved
                    if resolved[3] == "verified" and resolved[2] == nucleus:
                        ids[side].add(int(oid))
            candidates = overlaps[
                overlaps.A_ObjectID.isin(ids["A"]) | overlaps.B_ObjectID.isin(ids["B"])
            ].copy()
            # Resolve the opposite marker only for candidates actually examined.
            for side, name, col in (("A", a, "A_ObjectID"), ("B", b, "B_ObjectID")):
                needed = set(candidates[col].astype(int)) - set(idmaps[side])
                if needed:
                    for oid, group in pair_pixel_tables[side][pair_pixel_tables[side].ObjectID.isin(needed)].groupby("ObjectID", sort=False):
                        idmaps[side][int(oid)] = resolve_pair_roi(group, sources[name])
            final_set = set(zip(matching.A_ObjectID.astype(int), matching.B_ObjectID.astype(int)))
            output_ids = {(int(row.A_ObjectID), int(row.B_ObjectID)): int(row.CellID) for row in maps.itertuples(index=False)}
            audit_by_pair = {(int(row.A_ObjectID), int(row.B_ObjectID)): row for row in audit.itertuples(index=False)} if audit is not None else {}
            candidate_a = set(candidates.A_ObjectID.astype(int))
            candidate_b = set(candidates.B_ObjectID.astype(int))
            candidate_xy = {
                "A": {int(oid): grp[["X", "Y"]].to_numpy(dtype=np.int64)
                      for oid, grp in a_df[a_df.ObjectID.isin(candidate_a)].groupby("ObjectID", sort=False)},
                "B": {int(oid): grp[["X", "Y"]].to_numpy(dtype=np.int64)
                      for oid, grp in b_df[b_df.ObjectID.isin(candidate_b)].groupby("ObjectID", sort=False)},
            }
            # Cache native ROI fluorescence statistics once per object.
            pair_img_a, pair_img_b = h5["images/object_a_original"], h5["images/object_b_original"]
            native_stats = {}
            for side, img in (("A", pair_img_a), ("B", pair_img_b)):
                for oid, pixels in candidate_xy[side].items():
                    values = values_for_roi(img, pixels)
                    native_stats[(side, oid)] = (int(len(pixels)), float(values.mean()), float(values.sum()))
            n_final_for_nucleus = 0
            n_candidates_for_nucleus = 0
            for row in candidates.itertuples(index=False):
                aid, bid = int(row.A_ObjectID), int(row.B_ObjectID)
                if aid not in ids["A"] and bid not in ids["B"]:
                    continue
                n_candidates_for_nucleus += 1
                ar, ac, an, ast = idmaps["A"].get(aid, (None,None,None,"missing"))
                br, bc, bn, bst = idmaps["B"].get(bid, (None,None,None,"missing"))
                same = (an == bn) if ast == "verified" and bst == "verified" else None
                matched = (aid, bid) in final_set
                if stage == 9 and matched and same is not True:
                    raise ValueError(f"{path}: Stage-09 final pair {(aid, bid)} does not have verified same-nucleus identity")
                n_final_for_nucleus += int(matched)
                aud = audit_by_pair.get((aid, bid))
                record = {
                    "Field_Analysis_Path": str(analysis),
                    "Stage": stage, "Pair": f"{a}_vs_{b}", "A_Marker": a, "B_Marker": b,
                    "A_Current_ObjectID": aid, "B_Current_ObjectID": bid,
                    "A_Original_Marker_ObjectID": ar, "B_Original_Marker_ObjectID": br,
                    "A_Stage07_CellID": ac, "B_Stage07_CellID": bc,
                    "A_Stage07_Cell_Key": f"{analysis}::{a}::C{ac}" if ac is not None else None,
                    "B_Stage07_Cell_Key": f"{analysis}::{b}::C{bc}" if bc is not None else None,
                    "A_Canonical_Nuclear_ObjectID": an, "B_Canonical_Nuclear_ObjectID": bn,
                    "A_Identity_Status": ast, "B_Identity_Status": bst,
                    "Same_Hoechst_Nucleus": same,
                    "A_Percent_Inside_B": float(row.A_Percent_Inside_B),
                    "B_Percent_Inside_A": float(row.B_Percent_Inside_A),
                    "Shared_Area_Pixels": int(row.Shared_Area),
                    "A_ROI_Area_Pixels": native_stats[("A", aid)][0],
                    "B_ROI_Area_Pixels": native_stats[("B", bid)][0],
                    "A_Mean_Fluorescence_Native": native_stats[("A", aid)][1],
                    "B_Mean_Fluorescence_Native": native_stats[("B", bid)][1],
                    "A_Integrated_Fluorescence_Native": native_stats[("A", aid)][2],
                    "B_Integrated_Fluorescence_Native": native_stats[("B", bid)][2],
                    "Original_Overlap_Threshold_Pass": bool(row.A_Percent_Inside_B >= float(meta["matching"]["threshold_percent"]) - 1e-8),
                    "Stage09_Constrained_Eligibility": bool(aud.Final_Constrained_Eligibility) if aud is not None else None,
                    "Stage09_Exclusion_Reason": str(aud.Exclusion_Reason) if aud is not None and pd.notna(aud.Exclusion_Reason) else "",
                    "Final_Stable_Match": matched,
                    "PairID": pair_ids.get((aid, bid)),
                    "Stage_Local_Pairwise_CellID": output_ids.get((aid, bid)),
                    "Source_H5": str(path),
                }
                details.append(record)
                views.append({**record, "file": path, "rank": (int(matched), float(row.A_Percent_Inside_B)),
                              "A_XY": candidate_xy["A"][aid], "B_XY": candidate_xy["B"][bid]})
            summaries.append({"Stage": stage, "Pair": f"{a}_vs_{b}",
                              "Query_Nucleus": nucleus, "N_Overlapping_Candidates_Involving_Nucleus": n_candidates_for_nucleus,
                              "N_Final_Matches_Involving_Nucleus": n_final_for_nucleus,
                              "Status": "complete"})
    if not files:
        summaries.append({"Stage": stage, "Pair": "(no result HDF5)", "Query_Nucleus": nucleus,
                          "N_Overlapping_Candidates_Involving_Nucleus": 0,
                          "N_Final_Matches_Involving_Nucleus": 0,
                          "Status": "stage_not_available"})
    return details, summaries, views


def make_marker_figure(crops: list[dict], nucleus: int):
    n = max(1, len(crops))
    fig, ax = plt.subplots(n, 3, figsize=(13.6, 3.45 * n + 1.0), squeeze=False)
    for row in ax:
        for a in row:
            a.axis("off")
    for i, entry in enumerate(crops):
        a, b = entry["ca"], entry["cb"]
        m, nm = outline(entry["xy"], a.shape, entry["x0"], entry["y0"]), outline(entry["nxy"], b.shape, entry["x0"], entry["y0"])
        ca, cb = norm(a), norm(b)
        ia, ib = np.repeat(ca[:, :, None], 3, axis=2), np.repeat(cb[:, :, None], 3, axis=2)
        ia[m] = (0.12, .77, 1.); ib[nm] = (1., .3, .2)
        overlay = np.stack((ca, cb, np.zeros_like(ca)), axis=-1)
        overlay[m] = (.12, .77, 1.); overlay[nm] = (1., .3, .2)
        for col, img, title in ((0,ia,f'{entry["marker"]} · CellID {entry["cell_id"]} · ROI {entry["roi_id"]}'),
                                (1,ib,f'Hoechst · nucleus {nucleus}'), (2,overlay,'Shared XY crop')):
            ax[i, col].imshow(img, interpolation="nearest")
            ax[i, col].set_title(title, fontsize=10)
    if not crops:
        ax[0, 1].text(.5,.5,"No Stage-07 marker assigned to this nucleus", ha="center")
    fig.suptitle(f'Field-local canonical Hoechst nucleus {nucleus} — Stage-07 marker associations', fontsize=15, weight='bold')
    fig.text(.5, .008, 'Images contrast-scaled for display only; numerical measurements use raw pixel values.', ha='center', fontsize=9)
    fig.tight_layout(rect=(0,.035,1,.95))
    return fig


def make_pair_figure(item: dict, nucleus: int):
    with h5py.File(item["file"], "r") as h5:
        aa, bb = item["A_XY"], item["B_XY"]
        if not len(aa) or not len(bb):
            raise ValueError(f'{item["file"]}: final/candidate pair lacks ROI pixel coordinates')
        im_a, im_b = h5["images/object_a_original"], h5["images/object_b_original"]
        if im_a.shape != im_b.shape:
            raise ValueError(f'{item["file"]}: source image coordinate mismatch')
        x0,x1,y0,y1 = crop_bounds(np.vstack([aa,bb]), im_a.shape)
        ca, cb = np.asarray(im_a[y0:y1,x0:x1]), np.asarray(im_b[y0:y1,x0:x1])
    va, vb = norm(ca), norm(cb)
    ma, mb = outline(aa, ca.shape, x0,y0), outline(bb, ca.shape, x0,y0)
    img_a, img_b = np.repeat(va[:,:,None],3,axis=2), np.repeat(vb[:,:,None],3,axis=2)
    img_a[ma] = (1.,.2,.2); img_b[mb] = (.2,.5,1.)
    composite = np.stack((va,vb,np.zeros_like(va)),axis=-1)
    composite[ma] = (1.,.2,.2); composite[mb] = (.2,.5,1.)
    fig, axs = plt.subplots(1,4,figsize=(16,4.5),gridspec_kw={"width_ratios":[1,1,1,1.55]})
    for a,im,t in zip(axs[:3], (img_a,img_b,composite),(f'{item["A_Marker"]} Stage-07 CellID {item["A_Stage07_CellID"]}',
                       f'{item["B_Marker"]} Stage-07 CellID {item["B_Stage07_CellID"]}','A red / B green (B outline blue)')):
        a.imshow(im,interpolation='nearest');a.axis('off');a.set_title(t,fontsize=10)
    same = item["Same_Hoechst_Nucleus"]
    relation = 'YES' if same is True else ('NO' if same is False else 'UNRESOLVED')
    txt = (f'Stage {item["Stage"]} · {item["Pair"]}\n'
           f'PairID: {item.get("PairID") or "(not a final pair)"}\n'
           f'SAME NUCLEUS: {relation}\n'
           f'A Stage-07 CellID: {item["A_Stage07_CellID"]}   B: {item["B_Stage07_CellID"]}\n'
           f'A nuclear ID: {item["A_Canonical_Nuclear_ObjectID"]}\n'
           f'B nuclear ID: {item["B_Canonical_Nuclear_ObjectID"]}\n'
           f'A ROI: {item["A_Original_Marker_ObjectID"]}   B ROI: {item["B_Original_Marker_ObjectID"]}\n'
           f'Final stable match: {item["Final_Stable_Match"]}\n'
           f'A inside B: {item["A_Percent_Inside_B"]:.2f}%\n'
           f'B inside A: {item["B_Percent_Inside_A"]:.2f}%\n'
           f'Overlap pixels: {item["Shared_Area_Pixels"]}\n'
           f'A/B area (pixels): {item["A_ROI_Area_Pixels"]}/{item["B_ROI_Area_Pixels"]}\n'
           f'A/B native mean: {item["A_Mean_Fluorescence_Native"]:.4g}/{item["B_Mean_Fluorescence_Native"]:.4g}\n'
           f'A/B native sum: {item["A_Integrated_Fluorescence_Native"]:.4g}/{item["B_Integrated_Fluorescence_Native"]:.4g}\n'
           f'Output pair CellID: {item["Stage_Local_Pairwise_CellID"]}')
    axs[3].axis('off');axs[3].text(0,.95,txt,va='top',fontsize=10,linespacing=1.4)
    fig.suptitle(f'Field nucleus {nucleus} · Stage {item["Stage"]} (matched pairs first)', fontsize=13)
    fig.tight_layout(rect=(0,.02,1,.9))
    return fig


def run_pair_query(
    analysis: Path,
    sources: dict[str, MarkerSource],
    pair_id: str,
    out: Path | None,
    *,
    overwrite: bool,
) -> int:
    mappings = collect_pair_mappings(analysis, sources)
    selected = mappings[mappings.PairID.eq(pair_id)]
    if len(selected) != 1:
        available = ", ".join(mappings.PairID.astype(str).head(12))
        raise ValueError(
            f"PairID {pair_id!r} was not found exactly once in this field. "
            f"Available examples: {available or '(none)'}"
        )
    query = selected.iloc[0]
    key_mask = (
        mappings.A_Marker.eq(query.A_Marker)
        & mappings.B_Marker.eq(query.B_Marker)
        & mappings.A_ObjectID.astype(int).eq(int(query.A_ObjectID))
        & mappings.B_ObjectID.astype(int).eq(int(query.B_ObjectID))
    )
    linked = mappings[key_mask].sort_values("Stage")
    measurement_rows, views = [], []
    by_stage = {int(row.Stage): row for _, row in linked.iterrows()}
    if 8 in by_stage and 9 in by_stage:
        for column in (
            "A_Stage07_Source_CellID", "B_Stage07_Source_CellID",
            "A_Canonical_Nuclear_ObjectID", "B_Canonical_Nuclear_ObjectID",
            "Same_Nucleus",
        ):
            left, right = by_stage[8].get(column, pd.NA), by_stage[9].get(column, pd.NA)
            if pd.isna(left) and pd.isna(right):
                continue
            equal = (
                str(left) == str(right)
                if column == "Same_Nucleus"
                else pd.notna(left) and pd.notna(right) and float(left) == float(right)
            )
            if not equal:
                raise ValueError(
                    f"Cross-variant provenance mismatch for {pair_id}, column {column}: "
                    f"{left!r} vs {right!r}"
                )
    for stage, row in sorted(by_stage.items()):
        measured, view = pair_measurement(Path(row.Source_H5), row)
        measurement_rows.append(measured)
        views.append(view)
    crosswalk_rows = []
    for stage in (8, 9):
        row = by_stage.get(stage)
        provenance = query if row is None else row
        crosswalk_rows.append({
            "Query_PairID": pair_id,
            "Stage": stage,
            "PairID": None if row is None else row.PairID,
            "Present": row is not None,
            "A_Marker": query.A_Marker,
            "B_Marker": query.B_Marker,
            "A_ObjectID": int(query.A_ObjectID),
            "B_ObjectID": int(query.B_ObjectID),
            "A_Original_ROI_ObjectID": provenance.get("A_Original_ROI_ObjectID", pd.NA),
            "B_Original_ROI_ObjectID": provenance.get("B_Original_ROI_ObjectID", pd.NA),
            "A_Stage07_Source_CellID": provenance.get("A_Stage07_Source_CellID", pd.NA),
            "B_Stage07_Source_CellID": provenance.get("B_Stage07_Source_CellID", pd.NA),
            "A_Canonical_Nuclear_ObjectID": provenance.get("A_Canonical_Nuclear_ObjectID", pd.NA),
            "B_Canonical_Nuclear_ObjectID": provenance.get("B_Canonical_Nuclear_ObjectID", pd.NA),
            "Same_Nucleus": "ABSENT" if row is None else row.Same_Nucleus,
            "Tracking_Category": (
                "common_to_both" if len(by_stage) == 2
                else "stage08_only" if 8 in by_stage
                else "stage09_only"
            ),
        })
    folder = (out or analysis / "pair_id_queries").expanduser().resolve() / pair_id
    names = [
        "pair_query_summary.json",
        "pair_query_measurements.csv",
        "pair_query_stage08_stage09_cross_reference.csv",
        "pair_query_report.pdf",
    ]
    if not overwrite and any((folder / name).exists() for name in names):
        raise FileExistsError(f"Pair query already exists: {folder}. Use --overwrite to regenerate")
    folder.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(measurement_rows).to_csv(folder / names[1], index=False)
    pd.DataFrame(crosswalk_rows).to_csv(folder / names[2], index=False)
    with PdfPages(folder / names[3]) as pdf:
        for view in views:
            title_nucleus = view.get("A_Canonical_Nuclear_ObjectID")
            fig = make_pair_figure(view, title_nucleus)
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)
        fig = plt.figure(figsize=(11, 7), facecolor="white")
        fig.text(.06, .92, f"Pair cross-variant tracking · {pair_id}", fontsize=18, weight="bold")
        lines = [
            f"Marker pair: {query.A_Marker} vs {query.B_Marker}",
            f"Underlying final-pair key: A ObjectID {int(query.A_ObjectID)} / B ObjectID {int(query.B_ObjectID)}",
        ]
        for row in crosswalk_rows:
            lines.append(
                f"Stage {row['Stage']}: "
                + (f"PRESENT as {row['PairID']} · SAME NUCLEUS: {row['Same_Nucleus']}" if row["Present"] else "ABSENT")
            )
        fig.text(.07, .78, "\n\n".join(lines), fontsize=12, va="top")
        fig.text(.07, .10, "Pair identity is field- and variant-scoped. Cross-variant tracking uses the exact A/B ObjectID pair, not pair-number equality.", fontsize=9)
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)
    summary = {
        "analysis": str(analysis),
        "query_type": "pair_id",
        "queried_pair_id": pair_id,
        "queried_stage": int(query.Stage),
        "A_marker": str(query.A_Marker),
        "B_marker": str(query.B_Marker),
        "A_ObjectID": int(query.A_ObjectID),
        "B_ObjectID": int(query.B_ObjectID),
        "A_Stage07_Source_CellID": optional_int(query.A_Stage07_Source_CellID),
        "B_Stage07_Source_CellID": optional_int(query.B_Stage07_Source_CellID),
        "A_Canonical_Nuclear_ObjectID": optional_int(query.A_Canonical_Nuclear_ObjectID),
        "B_Canonical_Nuclear_ObjectID": optional_int(query.B_Canonical_Nuclear_ObjectID),
        "same_nucleus": str(query.Same_Nucleus),
        "present_in_stage08": 8 in by_stage,
        "present_in_stage09": 9 in by_stage,
        "stage08_pair_id": None if 8 not in by_stage else str(by_stage[8].PairID),
        "stage09_pair_id": None if 9 not in by_stage else str(by_stage[9].PairID),
        "report_files": names,
    }
    (folder / names[0]).write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    print("Saved:", folder)
    return 0


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis", required=True, type=Path, help="Existing field analysis directory")
    selector = parser.add_mutually_exclusive_group(required=True)
    selector.add_argument("--cell-id", type=int, help="STAGE-07 marker-specific CellID, not pairwise CellID")
    selector.add_argument("--nuclear-id", type=int, help="Canonical Hoechst nuclear ObjectID")
    selector.add_argument("--pair-id", help="Stage-08 or Stage-09 deterministic colocalization PairID")
    parser.add_argument("--marker", help="Marker name when using --cell-id, e.g. LHX6, PV, BiVe3")
    parser.add_argument("--out", type=Path, help="Explicit output directory; otherwise <analysis>/cell_id_queries")
    parser.add_argument("--max-pair-visuals", type=int, default=8, help="Maximum pair example pages per stage (default 8)")
    parser.add_argument("--overwrite", action="store_true", help="Replace report files in the chosen query folder")
    args = parser.parse_args(argv)
    if args.max_pair_visuals < 0:
        parser.error("--max-pair-visuals cannot be negative")
    analysis = args.analysis.expanduser().resolve()
    if not analysis.is_dir():
        parser.error(f"Missing analysis folder: {analysis}")
    sources = collect_sources(analysis)
    if args.pair_id is not None:
        if args.marker:
            parser.error("--marker is only valid with --cell-id")
        return run_pair_query(
            analysis,
            sources,
            args.pair_id,
            args.out,
            overwrite=args.overwrite,
        )
    nucleus, selected_as = select_nucleus(sources, args.cell_id, args.marker, args.nuclear_id)
    markers, crops = measure_markers(analysis, sources, nucleus)
    all_rows, variant_rows, image_rows = [], [], []
    for stage in (8,9):
        rows, sums, imgs = analyze_variant(analysis, stage, sources, nucleus)
        all_rows.extend(rows);variant_rows.extend(sums);image_rows.extend(imgs)
    folder = (args.out or analysis / "cell_id_queries").expanduser().resolve() / f"nucleus_{nucleus}"
    names = ["marker_measurements.csv", "stage08_stage09_overlap_candidates.csv", "stage08_stage09_variant_summary.csv",
             "cell_id_report.pdf", "marker_fluorescence.png", "cell_identity_summary.json"]
    if not args.overwrite and any((folder / filename).exists() for filename in names):
        raise FileExistsError(f"Query already exists: {folder}. Use --overwrite to regenerate")
    folder.mkdir(parents=True, exist_ok=True)
    markers.to_csv(folder / names[0],index=False)
    df = pd.DataFrame(all_rows)
    if df.empty:
        df = pd.DataFrame(columns=["Stage", "Pair", "A_Stage07_CellID", "B_Stage07_CellID", "A_Canonical_Nuclear_ObjectID", "B_Canonical_Nuclear_ObjectID", "Same_Hoechst_Nucleus", "Final_Stable_Match", "Source_H5"])
    df.to_csv(folder / names[1],index=False)
    pd.DataFrame(variant_rows).to_csv(folder / names[2],index=False)
    fig = make_marker_figure(crops,nucleus)
    fig.savefig(folder / names[4],dpi=160,bbox_inches='tight')
    with PdfPages(folder / names[3]) as pdf:
        pdf.savefig(fig,bbox_inches='tight')
        plt.close(fig)
        for stage in (8,9):
            selected = sorted((i for i in image_rows if i["Stage"] == stage),
                              key=lambda v: v["rank"],reverse=True)[:args.max_pair_visuals]
            for item in selected:
                fig = make_pair_figure(item,nucleus)
                pdf.savefig(fig,bbox_inches='tight');plt.close(fig)
        if not image_rows:
            fig = plt.figure(figsize=(9,4))
            fig.text(.08,.7,"No spatial-overlap candidates involving the query nucleus were found\nin available Stage-08 or Stage-09 HDF5 results",fontsize=13)
            pdf.savefig(fig);plt.close(fig)
    summary = {"analysis":str(analysis), "queried_as":selected_as,"canonical_Hoechst_nuclear_ObjectID":nucleus,
               "available_Stage07_markers":list(sources),
               "stage08_available":any(x["Stage"] == 8 and x["Status"] == "complete" for x in variant_rows),
               "stage09_available":any(x["Stage"] == 9 and x["Status"] == "complete" for x in variant_rows),
               "number_assigned_markers":int(markers.Stage07_Assigned.sum()),
               "stage08_final_pairs_involving_query_nucleus":sum(bool(r["Final_Stable_Match"]) for r in all_rows if r["Stage"] == 8),
               "stage09_final_pairs_involving_query_nucleus":sum(bool(r["Final_Stable_Match"]) for r in all_rows if r["Stage"] == 9),
               "report_files":names}
    (folder / names[5]).write_text(json.dumps(summary,indent=2)+"\n")
    print(json.dumps(summary,indent=2)); print("Saved:",folder)
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (ValueError, OSError, KeyError) as error:
        sys.exit(f"ERROR: {error}")
