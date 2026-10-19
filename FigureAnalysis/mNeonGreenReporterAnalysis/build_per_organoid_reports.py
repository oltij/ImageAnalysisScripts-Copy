#!/usr/bin/env python3
"""Read-only Experiment 17 ventral per-field microscopy QC PDF reports.

Uses the existing mNeonGreenReporterAnalysis selector and metrics, stored
final ROI pixels/metrics, and the committed colocalization analysis.h5 format.
Does not rerun segmentation, filtering, registration, or stable matching.
"""
from __future__ import annotations

import argparse
import importlib.util
import io
import json
import logging
import math
import re
import sys
import textwrap
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from scipy.ndimage import binary_erosion
import tifffile

HERE = Path(__file__).resolve().parent
LOGGER = logging.getLogger("organoid-reports")
CHANNELS = ("Hoechst", "GFP", "LHX6", "PV")
ROI_METRICS = {
    "area": ("Area_um2", "ROI area (µm²)"),
    "mean_intensity": ("Mean intensity", "ROI mean fluorescence intensity (a.u.)"),
    "integrated_intensity": ("Integrated intensity", "ROI integrated fluorescence intensity (a.u.·pixel)"),
    "circularity": ("Circularity", "ROI circularity"),
    "solidity": ("Solidity", "ROI solidity"),
    "eccentricity": ("Eccentricity", "ROI eccentricity"),
}


def read_table(h5: h5py.File, key: str) -> pd.DataFrame:
    name = f"tables/{key}"
    if name not in h5:
        raise ValueError(f"HDF5 missing {name}: {h5.filename}")
    return pd.read_csv(io.BytesIO(h5[name][()].tobytes()))


def to_pairs(df: pd.DataFrame, label: str) -> dict[int, np.ndarray]:
    id_col = next((x for x in ("ObjectID", "ROI", "CellID") if x in df.columns), None)
    if not id_col or not {"X", "Y"}.issubset(df.columns):
        raise ValueError(f"{label}: expected ObjectID/ROI/CellID, X, Y, found {list(df.columns)}")
    if df.empty:
        return {}
    pixels: dict[int, np.ndarray] = {}
    for oid, g in df.groupby(id_col, sort=True):
        pixels[int(oid)] = g[["X", "Y"]].to_numpy(dtype=np.int32)
    return pixels


def check_pixels(groups: dict[int, np.ndarray], shape: tuple[int, int], label: str) -> None:
    h, w = shape
    for object_id, xy in groups.items():
        if xy.size and ((xy[:, 0] < 0).any() or (xy[:, 0] >= w).any() or (xy[:, 1] < 0).any() or (xy[:, 1] >= h).any()):
            raise ValueError(f"{label}: ROI {object_id} contains out-of-bounds coordinates")


def load_image(path: Path) -> np.ndarray:
    arr = np.squeeze(tifffile.imread(path))
    if arr.ndim != 2:
        raise ValueError(f"Expected 2-D fluorescence image, got {arr.shape}: {path}")
    return arr


def fluorescence_image(analysis: Path, channel: str) -> tuple[np.ndarray, str] | None:
    # Prefer geometrically aligned files; the nuclear reference is the fixed canvas.
    registered = analysis / "02_registration" / channel / f"aligned_{channel}_max.tif"
    if registered.is_file():
        return load_image(registered), str(registered)
    # The H5 from each colocalization run contains the exact registered TIFF
    # used for the matching calculation, so prefer it over a raw projection.
    for p in sorted(list((analysis / "07_nuclear_colocalization").glob("*/analysis.h5")) +
                    list((analysis / "08_marker_colocalization").glob("*/analysis.h5"))):
        with h5py.File(p, "r") as h5:
            meta = json.loads(h5["metadata/json"][()].decode())
            objects = meta["objects"]
            for suffix in ("a", "b"):
                if objects[f"object_{suffix}_name"].lower() == channel.lower():
                    return h5[f"images/object_{suffix}_original"][()], f"{p} (embedded object {suffix.upper()})"
    # Fallback for runs without alignment, or with no relevant pair H5.
    projection = analysis / "00_max_projections" / f"{channel}.tif"
    if projection.is_file():
        return load_image(projection), f"{projection} (projection fallback; verify geometric alignment)"
    return None


def norm_gray(image: np.ndarray, *, lo: float = 1, hi: float = 99) -> np.ndarray:
    finite = image[np.isfinite(image)]
    if finite.size == 0:
        return np.zeros(image.shape, dtype=np.float32)
    vmin, vmax = np.percentile(finite, [lo, hi])
    if vmax <= vmin:
        vmax = vmin + 1
    return np.clip((image.astype(np.float32) - vmin) / (vmax - vmin), 0, 1)


def edge(mask: np.ndarray) -> np.ndarray:
    return mask & ~binary_erosion(mask, structure=np.ones((3, 3), dtype=bool), border_value=0)


def make_label(groups: dict[int, np.ndarray], shape: tuple[int, int]) -> np.ndarray:
    label = np.zeros(shape, dtype=np.int32)
    check_pixels(groups, shape, "ROI label image")
    for oid, xy in groups.items():
        if xy.size:
            label[xy[:, 1], xy[:, 0]] = oid
    return label


def overlay(image: np.ndarray, labels: np.ndarray | None = None,
            selected: set[int] | None = None, organoid: np.ndarray | None = None,
            max_side: int = 1250) -> np.ndarray:
    """Rasterized overlays keep multipage PDFs relatively small."""
    step = max(1, int(math.ceil(max(image.shape) / max_side)))
    gray = norm_gray(image)[::step, ::step]
    result = np.repeat(gray[..., None], 3, axis=2)
    if labels is not None:
        label = labels[::step, ::step]
        boundaries = (label > 0) & ((label != np.roll(label, 1, axis=0)) |
                                    (label != np.roll(label, -1, axis=0)) |
                                    (label != np.roll(label, 1, axis=1)) |
                                    (label != np.roll(label, -1, axis=1)))
        selection = np.isin(label, list(selected)) if selected is not None else np.zeros_like(boundaries)
        result[boundaries & ~selection] = (0.99, 0.78, 0.06)  # all others: yellow
        result[boundaries & selection] = (0.0, 0.7, 1.0)  # selected: blue
    if organoid is not None:
        outline = edge(organoid.astype(bool))[::step, ::step]
        result[outline] = (0., 0., 0.)
    return result


def note_page(pdf: PdfPages, title: str, paragraphs: list[str], *, subtitle="") -> None:
    fig = plt.figure(figsize=(13, 8.2), facecolor="white")
    fig.text(.06, .93, title, fontsize=20, weight="bold", va="top")
    if subtitle:
        fig.text(.06, .865, subtitle, fontsize=11, color="dimgray")
    y = .78
    for paragraph in paragraphs:
        for line in textwrap.wrap(str(paragraph), width=126) or [""]:
            fig.text(.065, y, line, fontsize=10.3, va="top")
            y -= .034
        y -= .035
        if y < .08:
            break
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


def graph_page(pdf: PdfPages, fig, footer: str, title: str | None = None) -> None:
    if title:
        fig.suptitle(title, fontsize=14, weight="bold")
    fig.text(.5, .025, "\n".join(textwrap.wrap(footer, width=166)), ha="center", va="bottom", fontsize=7.7)
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


def get_reference_mask(analysis: Path, shape) -> tuple[np.ndarray | None, str | None]:
    base = analysis / "03_aligned_segmentation"
    for channel in ("Hoechst", "DAPI", "Hoescht"):
        path = base / channel / "OrganoidMask.tiff"
        if path.is_file():
            mask = load_image(path) > 0
            if mask.shape != shape:
                raise ValueError(f"Organoid mask {mask.shape} does not match aligned image {shape}: {path}")
            return mask, str(path)
    return None, None


def whole_page(pdf: PdfPages, analysis: Path, provenance: list[str]) -> None:
    images = [(channel, fluorescence_image(analysis, channel)) for channel in CHANNELS]
    images = [(label, found) for label, found in images if found is not None]
    if not images:
        note_page(pdf, "Whole-organoid fluorescence unavailable", ["No aligned source image or embedded colocalization image could be located."])
        return
    shape = images[0][1][0].shape
    mask, mask_path = get_reference_mask(analysis, shape)
    fig, axes = plt.subplots(2, 2, figsize=(15, 11.2))
    for ax in axes.flat:
        ax.axis("off")
    for ax, (channel, (image, src)) in zip(axes.flat, images):
        if image.shape != shape:
            raise ValueError(f"Misaligned channel {channel} {image.shape} != {shape}")
        ax.imshow(overlay(image, organoid=mask), interpolation="nearest")
        ax.set_title(channel + " · aligned whole field", fontsize=13)
    footer = "Actual aligned grayscale fluorescence; black organoid perimeter = Hoechst-defined tissue mask (when present). Per-channel brightness is independently stretched to the 1st–99th percentiles for display; no contrast adjustment is applied to quantitative intensities."
    graph_page(pdf, fig, footer, "Whole-organoid fluorescence and tissue outline")
    provenance.extend([f"Fluorescence {ch}: {src}" for ch, (_, src) in images])
    if mask_path:
        provenance.append("Tissue mask: " + mask_path)


def resolve_roi_files(analysis: Path, channel: str) -> tuple[Path, Path] | None:
    directory = analysis / "06_final_shape_filtered" / channel
    pixels = directory / f"{channel}_ROI_pixels_FILTERED.csv"
    metrics = directory / f"{channel}_kept_ROI_metrics.csv"
    if pixels.is_file() and metrics.is_file():
        return pixels, metrics
    return None


def single_roi_pages(pdf: PdfPages, analysis: Path, channel: str, metric_keys: list[str],
                     max_examples: int, seed: int, provenance: list[str]) -> None:
    pair = resolve_roi_files(analysis, channel)
    found = fluorescence_image(analysis, channel)
    if pair is None or found is None:
        note_page(pdf, f"{channel}: ROI report unavailable", ["Final filtered ROI pixels/metrics or aligned fluorescence are missing. No ROIs or numeric values were invented."])
        return
    image, source = found
    pixels = pd.read_csv(pair[0], usecols=lambda x: x in {"ROI", "ObjectID", "CellID", "X", "Y"})
    groups = to_pairs(pixels, channel)
    metrics = pd.read_csv(pair[1])
    if "ObjectID" not in metrics.columns:
        raise ValueError(f"Metrics missing ObjectID: {pair[1]}")
    if metrics["ObjectID"].duplicated().any():
        raise ValueError(f"Duplicate ObjectID in {pair[1]}")
    check_pixels(groups, image.shape, channel)
    label = make_label(groups, image.shape)
    all_metrics = metrics.set_index("ObjectID")
    present = set(all_metrics.index.astype(int)) & set(groups)
    if not present:
        note_page(pdf, f"{channel}: no matched final ROI IDs", ["Pixel and metric tables contain no common ROI IDs."])
        return
    if present != set(groups) or present != set(all_metrics.index.astype(int)):
        raise ValueError(f"{channel}: final ROI pixel IDs and ROI metric IDs disagree. Refusing mislabeled decile panels.")
    provenance.extend([f"{channel} ROI pixels: {pair[0]}", f"{channel} ROI metrics: {pair[1]}", f"{channel} fluorescence: {source}"])
    rng = np.random.default_rng(seed)
    # Each decile is equal-rank (not a percentage of tissue area).
    for metric_key in metric_keys:
        col, description = ROI_METRICS[metric_key]
        if col not in metrics.columns:
            note_page(pdf, f"{channel}: {description} unavailable", [f"The final ROI metrics CSV does not contain column {col!r}. This metric is skipped, rather than being estimated."])
            continue
        values = pd.to_numeric(all_metrics[col], errors="coerce")
        finite = values[np.isfinite(values)].sort_values(kind="mergesort")
        if finite.empty:
            note_page(pdf, f"{channel}: {description} unavailable", ["There are no finite values for this measurement."])
            continue
        bins = np.array_split(finite.index.to_numpy(dtype=int), 10)
        for decile, ids in enumerate(bins):
            if not len(ids):
                continue
            ids = ids.astype(int)
            chosen = rng.choice(ids, size=min(max_examples, len(ids)), replace=False)
            fig = plt.figure(figsize=(14.4, 9.5))
            gs = fig.add_gridspec(2, max(3, max_examples), height_ratios=[2.4, 1.3], hspace=.34)
            full_ax = fig.add_subplot(gs[0, :])
            full_ax.imshow(overlay(image, label, set(map(int, chosen))), interpolation="nearest")
            full_ax.set_title(f"Whole organoid · {len(ids)} ROIs in this 10% rank bin · {len(chosen)} displayed in blue")
            full_ax.axis("off")
            for k, oid in enumerate(chosen):
                ax = fig.add_subplot(gs[1, k])
                xy = groups[int(oid)]
                x0 = max(0, int(xy[:, 0].min()) - 18); x1 = min(image.shape[1], int(xy[:, 0].max()) + 19)
                y0 = max(0, int(xy[:, 1].min()) - 18); y1 = min(image.shape[0], int(xy[:, 1].max()) + 19)
                patch = image[y0:y1, x0:x1]
                roi_mask = label[y0:y1, x0:x1] == oid
                gray = np.repeat(norm_gray(patch)[..., None], 3, axis=2)
                gray[edge(roi_mask)] = (0., .7, 1.)
                ax.imshow(gray, interpolation="nearest")
                ax.set_title(f"ROI {oid}\n{float(values.loc[oid]):.4g}", fontsize=9)
                ax.axis("off")
            for k in range(len(chosen), max_examples):
                fig.add_subplot(gs[1, k]).axis("off")
            fig.suptitle(f"{channel}: {description} · rank {decile*10}–{(decile+1)*10}%", fontsize=15, weight="bold")
            graph_page(pdf, fig, "Blue = sampled ROI boundary; yellow = all other ROI boundaries; grayscale = original aligned fluorescence. Rank bins contain approximately equal numbers of ROIs, sorted by metric. Values shown on crop titles use the unmodified measurements. The whole-image view shows only the sampled ROIs in blue.")


@dataclass
class PairData:
    name: str
    a: str
    b: str
    img_a: np.ndarray
    img_b: np.ndarray
    a_pixels: dict[int, np.ndarray]
    b_pixels: dict[int, np.ndarray]
    matching: pd.DataFrame
    overlaps: pd.DataFrame
    threshold: float
    source: Path


def load_pair(path: Path) -> PairData:
    with h5py.File(path, "r") as h5:
        meta = json.loads(h5["metadata/json"][()].decode())
        a = meta["objects"]["object_a_name"]
        b = meta["objects"]["object_b_name"]
        ia = h5["images/object_a_original"][()]
        ib = h5["images/object_b_original"][()]
        aa = to_pairs(read_table(h5, "object_a_pixels"), f"{a} H5")
        bb = to_pairs(read_table(h5, "object_b_pixels"), f"{b} H5")
        matching = read_table(h5, "stable_matching")
        overlaps = read_table(h5, "overlap_results")
    if ia.shape != ib.shape or ia.ndim != 2:
        raise ValueError(f"{path}: incompatible image shapes {ia.shape}, {ib.shape}")
    check_pixels(aa, ia.shape, a)
    check_pixels(bb, ib.shape, b)
    expected_cols = {"A_ObjectID", "B_ObjectID"}
    if not expected_cols.issubset(matching.columns) or not expected_cols.issubset(overlaps.columns):
        raise ValueError(f"{path}: missing A_ObjectID/B_ObjectID in matching or overlap table")
    return PairData(f"{a} vs {b}", a, b, ia, ib, aa, bb, matching, overlaps,
                    float(meta["matching"]["threshold_percent"]), path)


def intensity_rows(pair: PairData, field_index: int, sample: str) -> pd.DataFrame:
    """Use native source fluorescence, not segmentation masks or display scaling."""
    records = []
    overlap_lookup = pair.overlaps.set_index(["A_ObjectID", "B_ObjectID"])
    for row in pair.matching.itertuples(index=False):
        aid, bid = int(row.A_ObjectID), int(row.B_ObjectID)
        if aid not in pair.a_pixels or bid not in pair.b_pixels:
            raise ValueError(f"Matched ROI IDs missing from pixels: {pair.name}: {aid}/{bid}")
        ap, bp = pair.a_pixels[aid], pair.b_pixels[bid]
        a_int = pair.img_a[ap[:, 1], ap[:, 0]].astype(np.float64)
        b_int = pair.img_b[bp[:, 1], bp[:, 0]].astype(np.float64)
        a_flat = ap[:, 1].astype(np.int64) * pair.img_a.shape[1] + ap[:, 0]
        b_flat = bp[:, 1].astype(np.int64) * pair.img_a.shape[1] + bp[:, 0]
        common = np.intersect1d(a_flat, b_flat, assume_unique=False)
        y, x = np.divmod(common, pair.img_a.shape[1])
        overlap = overlap_lookup.loc[(aid, bid)] if (aid, bid) in overlap_lookup.index else None
        records.append({
            "field_index": field_index, "sample": sample,
            "pair": pair.name, "object_a": pair.a, "object_b": pair.b,
            "A_ObjectID": aid, "B_ObjectID": bid,
            "A_roi_integrated_intensity_native_sum": float(a_int.sum()),
            "B_roi_integrated_intensity_native_sum": float(b_int.sum()),
            "A_roi_mean_intensity_native": float(a_int.mean()) if a_int.size else np.nan,
            "B_roi_mean_intensity_native": float(b_int.mean()) if b_int.size else np.nan,
            "A_roi_area_pixels": len(ap), "B_roi_area_pixels": len(bp),
            "overlap_area_pixels": len(common),
            "A_overlap_integrated_intensity_native_sum": float(pair.img_a[y, x].astype(np.float64).sum()),
            "B_overlap_integrated_intensity_native_sum": float(pair.img_b[y, x].astype(np.float64).sum()),
            "A_percent_inside_B": float(overlap["A_Percent_Inside_B"]) if overlap is not None else np.nan,
            "B_percent_inside_A": float(overlap["B_Percent_Inside_A"]) if overlap is not None else np.nan,
            "source_h5": str(pair.source),
        })
    return pd.DataFrame(records)


def pair_panel(axs, pair: PairData, aid: int, bid: int, title: str) -> None:
    ap, bp = pair.a_pixels[aid], pair.b_pixels[bid]
    xy = np.concatenate([ap, bp])
    x0 = max(0, int(xy[:, 0].min()) - 18); x1 = min(pair.img_a.shape[1], int(xy[:, 0].max()) + 19)
    y0 = max(0, int(xy[:, 1].min()) - 18); y1 = min(pair.img_a.shape[0], int(xy[:, 1].max()) + 19)
    am = np.zeros((y1-y0,x1-x0), dtype=bool)
    bm = np.zeros_like(am)
    am[ap[:, 1]-y0, ap[:, 0]-x0] = True
    bm[bp[:, 1]-y0, bp[:, 0]-x0] = True
    ga = norm_gray(pair.img_a[y0:y1,x0:x1])
    gb = norm_gray(pair.img_b[y0:y1,x0:x1])
    arr_a = np.repeat(ga[..., None], 3, axis=2)
    arr_b = np.repeat(gb[..., None], 3, axis=2)
    arr_a[edge(am)] = (1., .08, .08)
    arr_b[edge(bm)] = (.12, .42, 1.)
    composite = np.stack((ga, gb, np.zeros_like(ga)), axis=-1)
    composite[edge(am)] = (1., 0., 0.)
    composite[edge(bm)] = (0., 0., 1.)
    for ax, panel, name in zip(axs, (arr_b, arr_a, composite), (pair.b, pair.a, "A=red/B=green")):
        ax.imshow(panel, interpolation="nearest")
        ax.set_title(f"{name}: {title}", fontsize=8)
        ax.axis("off")


def colocalization_pages(pdf: PdfPages, pair: PairData, max_examples: int, seed: int,
                         field_index: int, sample: str) -> pd.DataFrame:
    intensities = intensity_rows(pair, field_index, sample)
    note_page(pdf, f"Colocalization: {pair.name}", [
        f"Upstream stable-matching threshold: {pair.threshold:g}% of Object A ({pair.a}) area inside Object B ({pair.b}). Matching identities are read from analysis.h5 unchanged.",
        "Matched: final one-to-one retained pairs. A-only: A objects not finally matched; B-only: B objects not finally matched. For unmatched examples, the best partial partner is DISPLAY CONTEXT, not a positive match.",
        "Each 10% bin is an ABSOLUTE DIRECTIONAL OVERLAP percentage (0–10%, …, 90–100%), not rank-based deciles or a percentage of organoid area. Empty bins are skipped and listed in the index CSV."
    ])
    matching_keys = {(int(x.A_ObjectID), int(x.B_ObjectID)) for x in pair.matching.itertuples(index=False)}
    matched_a = {a for a, _ in matching_keys}; matched_b = {b for _, b in matching_keys}
    overlaps = pair.overlaps.copy()
    for col in ("A_Percent_Inside_B", "B_Percent_Inside_A"):
        overlaps[col] = pd.to_numeric(overlaps[col], errors="coerce")
    categories = [("Final stable-matched A/B", pair.matching.merge(overlaps, on=["A_ObjectID", "B_ObjectID"], how="left", suffixes=("", "_overlap")), "A_Percent_Inside_B"),
                  ("A-only · best B context", overlaps[~overlaps["A_ObjectID"].isin(matched_a)].sort_values("A_Percent_Inside_B", ascending=False).drop_duplicates("A_ObjectID"), "A_Percent_Inside_B"),
                  ("B-only · best A context", overlaps[~overlaps["B_ObjectID"].isin(matched_b)].sort_values("B_Percent_Inside_A", ascending=False).drop_duplicates("B_ObjectID"), "B_Percent_Inside_A")]
    rng = np.random.default_rng(seed)
    for title, data, metric in categories:
        if metric not in data.columns:
            raise ValueError(f"{pair.source}: missing directional overlap metric {metric}")
        if data.empty:
            note_page(pdf, f"{pair.name}: {title}", ["No overlapping ROI pairs in this category; zero-overlap unmatched ROIs have no partner to show."])
            continue
        for bin_idx in range(10):
            start, stop = bin_idx*10, (bin_idx+1)*10
            values = pd.to_numeric(data[metric], errors="coerce")
            members = data[(values >= start) & ((values < stop) if bin_idx < 9 else (values <= stop))]
            if members.empty:
                continue
            chosen_idx = rng.choice(len(members), size=min(max_examples, len(members)), replace=False)
            selected = members.iloc[chosen_idx]
            fig, axes = plt.subplots(3, len(selected), figsize=(max(5.0, 4.2*len(selected)), 9.5), squeeze=False)
            for j, row in enumerate(selected.itertuples(index=False)):
                aid, bid = int(row.A_ObjectID), int(row.B_ObjectID)
                pair_panel(axes[:,j], pair, aid, bid, f"A {aid} / B {bid}")
            graph_page(pdf, fig,
                       f"{title}. {start}–{stop}% of {pair.a if metric.startswith('A_') else pair.b} ROI area contained in {pair.b if metric.startswith('A_') else pair.a}. "
                       "Top: Object B grayscale + blue selected outline. Middle: Object A grayscale + red selected outline. Bottom: red A / green B fluorescence composite, with red/blue contours. "
                       "Unmatched counterpart is best-overlap context ONLY; positive status is determined exclusively by upstream stable matching.",
                       f"{pair.name} · {title} · overlap {start}–{stop}% ({len(members)} candidates)")
    return intensities


def intensity_pages(pdf: PdfPages, pair_name: str, values: pd.DataFrame) -> None:
    if values.empty:
        note_page(pdf, f"Integrated fluorescence: {pair_name}", ["No finally colocalized ROI pairs; no intensities to summarize."])
        return
    fig, ax = plt.subplots(figsize=(12.4, 6.8))
    x = values["A_roi_integrated_intensity_native_sum"].to_numpy(float)
    y = values["B_roi_integrated_intensity_native_sum"].to_numpy(float)
    ax.scatter(x, y, alpha=.6, s=14)
    ax.set(xlabel=f"{values['object_a'].iloc[0]}: integrated native intensity sum (a.u.·pixel)",
           ylabel=f"{values['object_b'].iloc[0]}: integrated native intensity sum (a.u.·pixel)")
    ax.set_title(f"{pair_name}: all {len(values)} finally matched pairs")
    ax.grid(alpha=.2)
    stats = []
    for side in ("A", "B"):
        col = f"{side}_roi_integrated_intensity_native_sum"
        vals = values[col].astype(float)
        stats.append(f"{values['object_a' if side == 'A' else 'object_b'].iloc[0]}: n={len(vals)}, mean={vals.mean():.5g}, median={vals.median():.5g}, total={vals.sum():.6g}")
    graph_page(pdf, fig,
               " | ".join(stats) + ". Whole ROI integrated intensity is the raw sum of original aligned fluorescence pixels belonging to each finally matched ROI; overlap-only sums and one row per matched pair are exported to CSV. "
               "No background subtraction or gain normalization is applied; values cannot be treated as molecule counts or directly compared across anti-GFP protocols.")
    # Every matched pair's raw numerical measurement is also inside the PDF,
    # not solely in the companion CSV or scatter plot.
    shown_columns = [
        ("A_ObjectID", "A ROI"), ("B_ObjectID", "B ROI"),
        ("A_roi_integrated_intensity_native_sum", "A full ROI\nΣ fluorescence"),
        ("B_roi_integrated_intensity_native_sum", "B full ROI\nΣ fluorescence"),
        ("A_overlap_integrated_intensity_native_sum", "A overlap\nΣ fluorescence"),
        ("B_overlap_integrated_intensity_native_sum", "B overlap\nΣ fluorescence"),
    ]
    page_size = 32
    for start in range(0, len(values), page_size):
        chunk = values.iloc[start:start+page_size]
        formatted = []
        for row in chunk.to_dict("records"):
            formatted.append([str(int(row[k])) if k in ("A_ObjectID", "B_ObjectID") else f"{float(row[k]):.5g}" for k, _ in shown_columns])
        fig = plt.figure(figsize=(13.4, 9.2), facecolor="white")
        fig.text(.04, .94, f"{pair_name}: integrated fluorescence for every matched pair", fontsize=16, weight="bold")
        fig.text(.04, .90, f"Pairs {start+1}–{start+len(chunk)} of {len(values)}; native TIFF units, full ROI and spatial intersection", fontsize=10)
        ax = fig.add_axes([.035, .13, .93, .73]); ax.axis("off")
        tab = ax.table(cellText=formatted, colLabels=[x[1] for x in shown_columns],
                       loc="upper center", cellLoc="right", colLoc="center")
        tab.auto_set_font_size(False); tab.set_fontsize(7.4); tab.scale(1.0, 1.25)
        fig.text(.04, .055, "A/B ROI IDs are the original object IDs from the pair's analysis.h5. Σ fluorescence = raw sum of pixels in the corresponding ROI (or intersection); not background corrected.", fontsize=8)
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)


def existing_pair_overlay_page(pdf: PdfPages, pair: PairData, analysis: Path) -> None:
    vis = pair.source.parent / "visualizations"
    paths = [vis / f"{pair.a}_ROI_before_filtering.tif",
             vis / f"{pair.a}_ROI_after_filtering.tif",
             vis / f"{pair.b}_ROI_before_filtering.tif",
             vis / f"{pair.b}_ROI_after_filtering.tif"]
    available = [p.is_file() for p in paths]
    if not any(available):
        note_page(pdf, f"{pair.name}: prior ROI overlay TIFFs absent", [
            "New 10%-increment colocalization QC crops are included in this report; the optional existing full-image before/after overlays were not found."])
        return
    fig, axes = plt.subplots(2, 2, figsize=(13.8, 11))
    titles = [f"{pair.a} · all ROIs", f"{pair.a} · finally matched",
              f"{pair.b} · all ROIs", f"{pair.b} · finally matched"]
    for ax, path, title in zip(axes.flat, paths, titles):
        ax.axis("off")
        if path.is_file():
            arr = tifffile.imread(path)
            if arr.ndim != 3 or arr.shape[2] not in (3, 4):
                raise ValueError(f"Expected RGB overlay image: {path}")
            skip = max(1, math.ceil(max(arr.shape[:2])/1350))
            ax.imshow(arr[::skip, ::skip, :3], interpolation="nearest")
            ax.set_title(title)
        else:
            ax.set_title(title + " (unavailable)")
    graph_page(pdf, fig, "Existing colocalization visualization outputs, not newly computed. Each ROI has its own distinct outline color over its own channel's grayscale fluorescence. The 'after' panels retain the final one-to-one matched ROI set, which is not equivalent to intensity/shape filtering alone.", f"{pair.name}: upstream whole-field before/after overlays")


def appendix_pdfs(analysis: Path) -> list[Path]:
    """Optional verbatim appendix: retain exact upstream QC pages and captions."""
    candidates = []
    for directory in (analysis / "05_intensity_report", analysis / "06_shape_reports",
                      analysis / "07_nuclear_colocalization", analysis / "08_marker_colocalization"):
        if directory.is_dir():
            candidates.extend(directory.rglob("*.pdf"))
    return sorted(set(p.resolve() for p in candidates if p.is_file()))


def merge_existing_appendix(generated: Path, output: Path, analysis: Path, appendix: list[Path]) -> None:
    """Append existing PDF pages losslessly; never rasterize original PDFs."""
    from pypdf import PdfWriter
    writer = PdfWriter()
    writer.append(str(generated), import_outline=False)
    for p in appendix:
        try:
            writer.append(str(p), import_outline=False)
        except Exception as exc:
            LOGGER.warning("Skipping unreadable source QC PDF %s: %s", p, exc)
    with output.open("wb") as handle:
        writer.write(handle)
    writer.close()


def metric_palette(samples: list[str]):
    cmap = plt.get_cmap("tab20")
    return {sample: cmap(i % 20) for i, sample in enumerate(samples)}


def one_metric_page(pdf: PdfPages, df: pd.DataFrame, plotmod, metric: str, title: str,
                    unit: str, field_index: int, *, gfp: bool = False) -> None:
    samples = sorted(df["sample"].unique()); palette = metric_palette(samples)
    cell_lines = [c for c in plotmod.CELL_LINE_ORDER if c in set(df["cell_line"])]
    stains = list(plotmod.ANTI_GFP_ORDER) if gfp else [None]
    fig, axes = plt.subplots(len(stains), len(cell_lines), figsize=(5.2*len(cell_lines), 4.3*len(stains)+1), squeeze=False)
    y_values = pd.to_numeric(df[metric], errors="coerce")
    y_values = y_values[np.isfinite(y_values)]
    if y_values.empty:
        plt.close(fig)
        note_page(pdf, title + ": numeric data unavailable", [
            f"The source reporter-analysis table contains no finite values for {metric}. This plot cannot be drawn without inventing values."])
        return
    ylo, yhi = float(y_values.min()), float(y_values.max())
    pad = .07*(yhi-ylo) if yhi > ylo else max(1, abs(ylo)*.08)
    for ri, stain in enumerate(stains):
        for ci, cell_line in enumerate(cell_lines):
            ax = axes[ri, ci]
            sub = df[(df["cell_line"] == cell_line) & ((df["anti_gfp_status"] == stain) if gfp else True)]
            for k, condition in enumerate(plotmod.CONDITIONS):
                group = sub[sub["condition"] == condition]
                for sample in sorted(group["sample"].unique()):
                    records = group[group["sample"] == sample].sort_values("index")
                    x = np.array([k + ((int(i) * 37 % 101) / 100 - .5)*.40 for i in records["index"]])
                    yy = records[metric].to_numpy(float)
                    finite = np.isfinite(yy)
                    ax.scatter(x[finite], yy[finite], s=26, color=palette[sample], alpha=.82,
                               edgecolors="white", linewidths=.3, zorder=3)
                    target = records[(records["index"] == field_index) & np.isfinite(records[metric])]
                    if not target.empty:
                        v = target.iloc[0]
                        xx = k + ((field_index * 37 % 101) / 100 - .5)*.40
                        ax.scatter([xx], [float(v[metric])], color=[palette[sample]], s=130,
                                   edgecolors="black", linewidths=2.6, zorder=8)
                vals = pd.to_numeric(group[metric], errors="coerce").dropna()
                if not vals.empty:
                    ax.plot([k-.19, k+.19], [float(vals.mean())]*2, color="black", lw=1.7, zorder=5)
            ax.set_xticks([0,1], plotmod.CONDITIONS)
            ax.set_xlim(-.55,1.55)
            ax.set_ylim(ylo-pad, yhi+pad)
            ax.set_ylabel(unit, fontsize=9)
            ax.set_title(f"{cell_line}" + (f" · {stain}" if stain else "") + f" (n={len(sub)} fields)", fontsize=10)
            ax.grid(axis="y", alpha=.2)
            ax.spines[["top", "right"]].set_visible(False)
    handles = [Line2D([0], [0], marker="o", linestyle="none", markerfacecolor=palette[s], markeredgecolor="none", label=s) for s in samples]
    handles.append(Line2D([0], [0], marker="o", linestyle="none", markerfacecolor="white", markeredgecolor="black", markeredgewidth=2.5, markersize=10, label="THIS FIELD"))
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(.5, .11), ncol=min(7, len(handles)), fontsize=8, frameon=False)
    fig.subplots_adjust(top=.88, bottom=.28 if not gfp else .18, wspace=.28, hspace=.40)
    label = plotmod.METRIC_LEGENDS.get(metric, title)
    graph_page(pdf, fig, label + " Black ring = exact current FIELD dot (not all fields in the same sample); other dots = other fields. Black bar = condition-level mean across fields. All panels share the same y scale. Fields from one biological block are not independent biological replicates.", title)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--h1-workspace", type=Path, required=True)
    ap.add_argument("--remaining-workspace", type=Path, required=True)
    ap.add_argument("--sample-metadata", type=Path, default=HERE/"igf_exp17_ventral_block_mapping.csv")
    ap.add_argument("--h1-selection", type=Path, help="Override existing verified H1 fields CSV")
    ap.add_argument("--out", type=Path, help="Default: original H1 workspace/mentor_field_plots/per_organoid_reports")
    ap.add_argument("--only-sample", help="Block-section filter, e.g. 151-3")
    ap.add_argument("--only-field", help="Exact or substring field acquisition name")
    ap.add_argument("--max-fields", type=int, help="Limit number of fields after filtering, for testing")
    ap.add_argument("--max-examples", type=int, default=3, help="ROI/pair examples per nonempty bin (default 3)")
    ap.add_argument("--roi-metrics", default="area,mean_intensity,integrated_intensity,circularity,solidity,eccentricity",
                    help="Comma-separated ROI ranking metrics: " + ",".join(ROI_METRICS))
    ap.add_argument("--skip-single-rois", action="store_true")
    ap.add_argument("--skip-colocalization", action="store_true")
    ap.add_argument("--skip-summary-plots", action="store_true")
    ap.add_argument("--append-existing-qc-pdfs", action="store_true",
                    help="Append ALL existing intensity, shape and colocalization QC PDFs verbatim after the newly generated report; may greatly increase file size")
    ap.add_argument("--overwrite", action="store_true")
    args = ap.parse_args(argv)
    if args.max_examples < 1 or args.max_examples > 10 or (args.max_fields is not None and args.max_fields < 1):
        ap.error("--max-examples must be 1..10 and --max-fields must be positive")
    keys = [x.strip() for x in args.roi_metrics.split(",") if x.strip()]
    if not keys or any(x not in ROI_METRICS for x in keys) or len(set(keys)) != len(keys):
        ap.error("Invalid --roi-metrics. Choose unique names from: " + ", ".join(ROI_METRICS))
    plotpath = HERE / "plot_igf_exp17_ventral_with_gfp_plots.py"
    if not plotpath.is_file():
        ap.error(f"Expected sibling main analysis script: {plotpath}")
    spec = importlib.util.spec_from_file_location("_reporter_metrics", plotpath)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to import {plotpath}")
    plotmod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(plotmod)
    info = plotmod.load_sample_info(args.sample_metadata.expanduser().resolve())
    original = args.h1_workspace.expanduser().resolve()
    other = args.remaining_workspace.expanduser().resolve()
    selection = args.h1_selection or original / "mentor_field_plots/PV_density_per_field_verified_H1.csv"
    fields, _, _ = plotmod.collect_fields(original, other, selection, info)
    data = pd.DataFrame(fields)
    out = (args.out or original / "mentor_field_plots/per_organoid_reports").expanduser().resolve()
    chosen = data.copy()
    if args.only_sample:
        chosen = chosen[chosen["sample"] == args.only_sample]
    if args.only_field:
        chosen = chosen[chosen["source_field"].str.contains(args.only_field, regex=False)]
    chosen = chosen.sort_values("index")
    if args.max_fields:
        chosen = chosen.head(args.max_fields)
    if chosen.empty:
        raise ValueError("No selected verified Experiment 17 ventral fields match your filters")
    LOGGER.info("Generating %d organoid field reports", len(chosen))
    report_records = []
    for entry in chosen.to_dict("records"):
        name = re.sub(r"[^A-Za-z0-9_.-]+", "_", Path(entry["source_field"]).name).strip("_")
        target_dir = out / entry["sample"]
        dest = target_dir / f"{name}_organoid_report.pdf"
        if dest.exists() and not args.overwrite:
            LOGGER.info("Skip existing report (use --overwrite): %s", dest)
            continue
        LOGGER.info("Field %s (index %s): %s", entry["sample"], entry["index"], entry["source_field"])
        target_dir.mkdir(parents=True, exist_ok=True)
        partial = dest.with_name(dest.stem + ".incomplete.pdf")
        analysis = Path(entry["output_dir"])
        provenance = []
        intensity_dfs = []
        pair_candidates = sorted((analysis / "07_nuclear_colocalization").glob("*/analysis.h5")) + sorted((analysis / "08_marker_colocalization").glob("*/analysis.h5"))
        try:
            with PdfPages(partial) as pdf:
                note_page(pdf, f"Exp 17 · {entry['sample']} · Field {entry['index']}", [
                    f"Source: {entry['source_field']}",
                    f"Cell line: {entry['cell_line']} | Condition: {entry['condition']} | {entry['anti_gfp_status']} | Ventral organoid",
                    "One PDF represents one microscopy FIELD. The sample-section may contain multiple fields, and those fields are not independent biological replicates.",
                    "Sections: whole-field fluorescence; individual ROI rank deciles; pairwise colocalization absolute-overlap bins; integrated fluorescence for every finally matched pair; all 20 comparison figures with THIS FIELD circled in black.",
                    "Displays are contrast adjusted for viewing, but raw pixels are used for integrated intensity. Colocalization membership and thresholds are read unchanged from the previously completed analysis.",
                ])
                whole_page(pdf, analysis, provenance)
                if not args.skip_single_rois:
                    for channel in ("GFP", "LHX6", "PV"):
                        single_roi_pages(pdf, analysis, channel, keys, args.max_examples, int(entry["index"])+31, provenance)
                if not args.skip_colocalization:
                    if not pair_candidates:
                        note_page(pdf, "Colocalization outputs unavailable", ["No completed 07_nuclear_colocalization or 08_marker_colocalization analysis.h5 files found."])
                    for n, pairfile in enumerate(pair_candidates):
                        pair = load_pair(pairfile)
                        provenance.append("Colocalization HDF5: " + str(pairfile))
                        existing_pair_overlay_page(pdf, pair, analysis)
                        intensity = colocalization_pages(pdf, pair, args.max_examples, int(entry["index"])+n+456, int(entry["index"]), entry["sample"])
                        intensity_pages(pdf, pair.name, intensity)
                        if not intensity.empty:
                            intensity_dfs.append(intensity)
                        del pair
                if not args.skip_summary_plots:
                    note_page(pdf, "Comparison plots: black-ring field", [
                        "The large marker with a black outline is exactly the current acquisition field. Its location and value match the source comparative analysis. Other fields and sample colors are retained, as are original condition means.",
                        "The GFP plots separately display -2 antibody-amplified and -3 native-only fluorescence acquisitions. These must not be treated as directly comparable absolute intensity scales without technical validation.",
                    ])
                    for metric, title, unit, _ in plotmod.PLOT_METRICS:
                        one_metric_page(pdf, data, plotmod, metric, title, unit, int(entry["index"]))
                    for metric, title, unit, _ in plotmod.GFP_PLOT_METRICS:
                        one_metric_page(pdf, data, plotmod, metric, title, unit, int(entry["index"]), gfp=True)
                note_page(pdf, "Methods, legends, and provenance", [
                    "Fluorescence underlays: grayscale originals. Blue=selected single-channel ROI contour; yellow=other ROI contours. Organoid outline is the Hoechst tissue mask when available.",
                    "Colocalization: Object A red, Object B blue contours; composite fluorescence A red / B green. Matched pairs are final upstream stable matches. Alone/unmatched pair contexts are best overlapping alternatives, never treated as true matched positives.",
                    "Single ROI 10% groups are equal-rank deciles per measured feature. Colocalization 10% groups are absolute percent overlap of A inside B (or B inside A for B-only). Empty bins are omitted, not merged.",
                    "Integrated intensity = sum of native fluorescence pixel values across each full matched ROI. Overlap-only integrated intensity = sum inside the shared ROI pixels. These are image-native units, NOT background subtracted or exposure normalized.",
                    "Each report includes all metrics used by the source summary figure script (when finite), with a black ring identifying THIS FIELD, not all fields of the same organoid block-section.",
                ] + provenance)
            if args.append_existing_qc_pdfs:
                sources = appendix_pdfs(analysis)
                if sources:
                    indexed = partial.with_name(partial.stem + ".indexed.pdf")
                    combined = partial.with_name(partial.stem + ".combined.pdf")
                    with PdfPages(indexed) as pdf:
                        for i in range(0, len(sources), 13):
                            note_page(pdf, "Original upstream QC PDFs · appendix index",
                                      [str(p.relative_to(analysis)) for p in sources[i:i+13]])
                    from pypdf import PdfWriter
                    with PdfWriter() as writer:
                        writer.append(str(partial))
                        writer.append(str(indexed))
                        for source in sources:
                            writer.append(str(source))
                        with combined.open("wb") as handle:
                            writer.write(handle)
                    combined.replace(dest)
                    partial.unlink(missing_ok=True)
                    indexed.unlink(missing_ok=True)
                else:
                    partial.replace(dest)
            else:
                partial.replace(dest)
        except Exception:
            partial.unlink(missing_ok=True)
            partial.with_name(partial.stem + ".indexed.pdf").unlink(missing_ok=True)
            partial.with_name(partial.stem + ".combined.pdf").unlink(missing_ok=True)
            raise
        intensity_path = target_dir / f"{name}_colocalized_roi_integrated_intensity.csv"
        if intensity_dfs:
            pd.concat(intensity_dfs, ignore_index=True).to_csv(intensity_path, index=False)
        else:
            pd.DataFrame(columns=["field_index", "sample", "pair", "A_ObjectID", "B_ObjectID", "A_roi_integrated_intensity_native_sum", "B_roi_integrated_intensity_native_sum"]).to_csv(intensity_path, index=False)
        report_records.append({"index": entry["index"], "sample": entry["sample"], "source_field": entry["source_field"], "pdf": str(dest), "intensity_csv": str(intensity_path), "n_pair_hdf5": len(pair_candidates)})
        LOGGER.info("Saved %s", dest)
    if report_records:
        pd.DataFrame(report_records).to_csv(out / "report_manifest.csv", index=False)
    LOGGER.info("Finished %d reports in %s", len(report_records), out)
    return 0


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    try:
        sys.exit(main())
    except (ValueError, FileNotFoundError, KeyError) as exc:
        sys.exit(f"ERROR: {exc}")
