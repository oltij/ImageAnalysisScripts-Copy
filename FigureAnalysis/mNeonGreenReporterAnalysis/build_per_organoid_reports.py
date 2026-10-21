#!/usr/bin/env python3
"""Concise, ROI-traceable Experiment 17 ventral per-field microscopy QC PDFs.

Uses the existing mNeonGreenReporterAnalysis selector and metrics, stored
final ROI pixels/metrics, and the committed colocalization analysis.h5 format.
Shows only ROI objects participating in existing final stable matching.
Uses equal-count rank quartiles for ROI metrics and final matched pair overlap.
Every selected exemplar has its own images and native fluorescence measurements
on the same page, and the companion CSV contains only displayed exemplars.
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
                    list((analysis / "08_marker_colocalization").glob("*/analysis.h5")) +
                    list((analysis / "09_marker_colocalization_nucleus_constrained").glob("*/analysis.h5"))):
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


def quartile_examples(frame: pd.DataFrame, metric: str, id_columns: list[str],
                      max_examples: int = 10) -> list[dict]:
    """Return up to max_examples real records from each equal-count quartile.

    Records are ranked within eligible objects, not absolute 0-25% bins. Within a
    quartile, examples are chosen at approximately even rank intervals so the
    displayed set spans that quartile instead of repeating near-identical medians.
    """
    work = frame.copy()
    work[metric] = pd.to_numeric(work[metric], errors="coerce")
    work = work[np.isfinite(work[metric])].sort_values([metric] + id_columns, kind="mergesort").reset_index(drop=True)
    chunks = np.array_split(np.arange(len(work)), 4)
    out: list[dict] = []
    for q, inds in enumerate(chunks, 1):
        if not len(inds):
            out.append({"quartile": q, "count": 0, "range": None, "examples": []})
            continue
        span = (float(work.iloc[inds[0]][metric]), float(work.iloc[inds[-1]][metric]))
        if len(inds) <= max_examples:
            chosen_positions = list(range(len(inds)))
        else:
            chosen_positions = sorted(set(int(x) for x in np.linspace(0, len(inds)-1, num=max_examples)))
        examples = []
        for pos in chosen_positions:
            idx = int(inds[pos])
            row = work.iloc[idx].copy()
            examples.append({"row": row, "rank_in_quartile": pos + 1, "quartile_size": len(inds)})
        out.append({"quartile": q, "count": len(inds), "range": span, "examples": examples})
    return out


def img_crop(image: np.ndarray, xy: np.ndarray, padding: int = 18) -> tuple[int, int, int, int]:
    h, w = image.shape
    return (max(0, int(xy[:,0].min())-padding), min(w, int(xy[:,0].max())+padding+1),
            max(0, int(xy[:,1].min())-padding), min(h, int(xy[:,1].max())+padding+1))


def gray_with_limits(image: np.ndarray, bounds: tuple[float, float]) -> np.ndarray:
    lo, hi = bounds
    return np.clip((image.astype(np.float32) - lo) / max(hi-lo, 1e-12), 0, 1)


def image_limits(image: np.ndarray) -> tuple[float, float]:
    valid = image[np.isfinite(image)]
    if not valid.size:
        return (0.0, 1.0)
    lo, hi = np.percentile(valid, (1, 99))
    return float(lo), float(hi if hi > lo else lo + 1)


def roi_crop_rgba(image: np.ndarray, xy: np.ndarray, limits: tuple[float,float],
                  *, color=(0.05, .67, .96)) -> np.ndarray:
    x0,x1,y0,y1 = img_crop(image, xy, padding=20)
    g = gray_with_limits(image[y0:y1, x0:x1], limits)
    rgb = np.repeat(g[:,:,None],3,axis=2)
    m = np.zeros(g.shape, bool)
    m[xy[:,1]-y0, xy[:,0]-x0]=True
    rgb[edge(m)] = color
    return rgb


def roi_quartile_pages(pdf: PdfPages, analysis: Path, channel: str, metric_keys: list[str],
                       allowed_ids: set[int], field_index: int, sample: str,
                       provenance: list[str], *, examples_per_quartile: int = 10,
                       rows_per_page: int = 5) -> list[dict]:
    """Display up to N exemplar ROIs per quartile, with each row fully traceable.

    Every displayed row contains exactly one ROI's locator image, crop, and
    quantitative values side-by-side. Pages are split by metric and quartile, with
    multiple pages when needed for readability.
    """
    files = resolve_roi_files(analysis, channel)
    img = fluorescence_image(analysis, channel)
    if files is None or img is None:
        LOGGER.warning("No final ROI metric/pixel files or fluorescence for %s", channel)
        return []
    image, image_source = img
    pixels = pd.read_csv(files[0], usecols=lambda c: c in {"ROI", "ObjectID", "CellID", "X", "Y"})
    groups = to_pairs(pixels, channel)
    metrics = pd.read_csv(files[1])
    if "ObjectID" not in metrics:
        raise ValueError(f"{files[1]}: missing ObjectID")
    if metrics["ObjectID"].duplicated().any():
        raise ValueError(f"{files[1]}: duplicate ROI ObjectID")
    check_pixels(groups, image.shape, channel)
    eligible = (set(groups) & set(metrics["ObjectID"].astype(int)) & allowed_ids)
    if not eligible:
        LOGGER.warning("No final-stable-matched %s ROI IDs also present in final metric/pixel tables", channel)
        return []
    d = metrics[metrics["ObjectID"].isin(eligible)].copy()
    provenance.extend([f"{channel} final ROI pixels: {files[0]}", f"{channel} final ROI metrics: {files[1]}", f"{channel} image: {image_source}"])
    limits = image_limits(image)
    skip = max(1, math.ceil(max(image.shape) / 590))
    thumb = gray_with_limits(image[::skip, ::skip], limits)
    mask, _ = get_reference_mask(analysis, image.shape)
    if mask is not None:
        orgedge = edge(mask)[::skip, ::skip]
        thumb_rgb = np.repeat(thumb[:, :, None], 3, axis=2)
        thumb_rgb[orgedge] = (.05, .05, .05)
    else:
        thumb_rgb = np.repeat(thumb[:, :, None], 3, axis=2)
    output = []
    for key in metric_keys:
        col, display_name = ROI_METRICS[key]
        if col not in d:
            LOGGER.warning("%s: missing metric %s (omitting this quartile sheet)", channel, col)
            continue
        quartiles = quartile_examples(d, col, ["ObjectID"], max_examples=examples_per_quartile)
        for qpack in quartiles:
            q = qpack["quartile"]
            qcount = qpack["count"]
            qrange = qpack["range"]
            examples = qpack["examples"]
            if not examples:
                fig = plt.figure(figsize=(16.6, 4.6), facecolor='white')
                fig.text(.055, .85, f'{channel} · {display_name} · Q{q}', fontsize=17, weight='bold')
                fig.text(.055, .62, f'No eligible final-stable-matched ROI present in Q{q}.', fontsize=13)
                fig.text(.055, .18, 'Quartiles are equal-count rank groups among final stable-matched ROI IDs only.', fontsize=9)
                pdf.savefig(fig, bbox_inches='tight'); plt.close(fig)
                continue
            for page_start in range(0, len(examples), rows_per_page):
                page_examples = examples[page_start:page_start + rows_per_page]
                page_no = page_start // rows_per_page + 1
                page_total = math.ceil(len(examples) / rows_per_page)
                fig = plt.figure(figsize=(18.2, 2.35 * len(page_examples) + 1.9), facecolor='white')
                gs = fig.add_gridspec(len(page_examples), 3, left=.045, right=.97, top=.86, bottom=.08, wspace=.16, hspace=.38, width_ratios=[1.08, 1.18, 2.55])
                for ridx, example in enumerate(page_examples):
                    axloc = fig.add_subplot(gs[ridx, 0]); axcrop = fig.add_subplot(gs[ridx, 1]); axinfo = fig.add_subplot(gs[ridx, 2])
                    for ax in (axloc, axcrop, axinfo):
                        ax.axis('off')
                    chosen = example['row']
                    oid = int(chosen['ObjectID']); xy = groups[oid]
                    x, y = float(xy[:, 0].mean()), float(xy[:, 1].mean())
                    axloc.imshow(thumb_rgb, interpolation='nearest'); axloc.scatter([x/skip], [y/skip], s=120, marker='o', facecolors='none', edgecolors='#00bfff', linewidths=2.0)
                    axloc.set_title(f'Field locator', fontsize=9)
                    axcrop.imshow(roi_crop_rgba(image, xy, limits), interpolation='nearest')
                    axcrop.set_title(f'{channel} ROI {oid} · fluorescence', fontsize=9)
                    native = image[xy[:, 1], xy[:, 0]].astype(np.float64)
                    values = {"Area (µm²)": chosen.get('Area_um2', float('nan')),
                              "Mean intensity (native)": float(native.mean()),
                              "Integrated intensity (native sum)": float(native.sum()),
                              "Circularity": chosen.get('Circularity', float('nan')),
                              "Solidity": chosen.get('Solidity', float('nan')),
                              "Eccentricity": chosen.get('Eccentricity', float('nan'))}
                    lines = [
                        f'Q{q} · ROI #{oid} · example {page_start + ridx + 1} of {len(examples)} displayed ({example["quartile_size"]} total in quartile)',
                        f'Rank metric: {display_name} = {float(chosen[col]):.5g}    Quartile range: {qrange[0]:.5g}–{qrange[1]:.5g}',
                        f'Area: {float(values["Area (µm²)"]):.5g} µm²    Mean fluorescence: {values["Mean intensity (native)"]:.5g}    Integrated fluorescence Σ: {values["Integrated intensity (native sum)"]:.6g}',
                        f'Circularity: {float(values["Circularity"]):.4g}    Solidity: {float(values["Solidity"]):.4g}    Eccentricity: {float(values["Eccentricity"]):.4g}',
                    ]
                    axinfo.text(.01, .99, lines[0], weight='bold', fontsize=10.4, va='top')
                    axinfo.text(.01, .76, "\n".join(lines[1:]), fontsize=9.5, va='top', linespacing=1.5)
                    output.append({"field_index": field_index, "sample": sample, "channel": channel, "quartile": f'Q{q}',
                                   "rank_metric": col, "rank_value": float(chosen[col]), "quartile_n": qcount,
                                   "quartile_range_min": qrange[0], "quartile_range_max": qrange[1], "ROI_ObjectID": oid,
                                   "display_order_within_quartile": page_start + ridx + 1, "displayed_examples_in_quartile": len(examples),
                                   "roi_pixel_count": len(xy), "Area_um2": values["Area (µm²)"],
                                   "mean_native_intensity": values["Mean intensity (native)"],
                                   "integrated_native_intensity_sum": values["Integrated intensity (native sum)"],
                                   "Circularity": values["Circularity"], "Solidity": values["Solidity"], "Eccentricity": values["Eccentricity"]})
                fig.text(.045, .965, f'{channel} · {display_name} · Q{q} final-matched ROI exemplars', fontsize=16.5, weight='bold')
                fig.text(.045, .928, f'Sample {sample} | Field {field_index} | {len(eligible)} ROI IDs eligible through final stable matches | displaying {len(page_examples)} of {len(examples)} requested exemplars (page {page_no}/{page_total})', fontsize=10.2)
                fig.text(.045, .03, 'Each row is one ROI only: left = whole-field location, middle = the fluorescence crop for that same ROI ID, right = values for that same ROI ID. Cyan circle/outline marks the displayed ROI. Quartiles are equal-count rank groups, not 0–25% absolute bins.', fontsize=8.8)
                pdf.savefig(fig, bbox_inches='tight'); plt.close(fig)
    return output


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
    variant: str


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
    variant = meta.get("analysis_variant", "unconstrained")
    display_name = (
        f"{a} vs {b} [NUCLEUS-CONSTRAINED]"
        if variant == "nucleus_constrained_marker_marker"
        else f"{a} vs {b}"
    )
    return PairData(display_name, a, b, ia, ib, aa, bb, matching, overlaps,
                    float(meta["matching"]["threshold_percent"]), path, variant)


def colocalization_files(analysis: Path, marker_pair_variant: str) -> list[Path]:
    nuclear = sorted((analysis / "07_nuclear_colocalization").glob("*/analysis.h5"))
    if marker_pair_variant == "unconstrained":
        marker = sorted((analysis / "08_marker_colocalization").glob("*/analysis.h5"))
    elif marker_pair_variant == "nucleus-constrained":
        marker = sorted(
            (analysis / "09_marker_colocalization_nucleus_constrained").glob("*/analysis.h5")
        )
        if not marker:
            raise FileNotFoundError(
                f"No Stage-09 nucleus-constrained marker-pair analysis.h5 files under {analysis}"
            )
    else:
        raise ValueError(f"Unknown marker-pair variant: {marker_pair_variant}")
    return nuclear + marker


def intensity_rows(pair: PairData, field_index: int, sample: str, selected_matches: pd.DataFrame | None = None) -> pd.DataFrame:
    """Use native source fluorescence, not segmentation masks or display scaling."""
    records = []
    overlap_lookup = pair.overlaps.set_index(["A_ObjectID", "B_ObjectID"])
    matches = selected_matches if selected_matches is not None else pair.matching
    for row in matches.itertuples(index=False):
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


def matched_pair_quartile_pages(pdf: PdfPages, pair: PairData, field_index: int, sample: str, *,
                              examples_per_quartile: int = 10, rows_per_page: int = 5) -> pd.DataFrame:
    """Display up to N final stable-matched pairs per overlap-rank quartile.

    Each displayed row contains the B image, A image, composite, and exact
    measurements for those same A/B ROI IDs. Only final stable matches are used.
    """
    if pair.matching.empty:
        LOGGER.info("%s has zero final stable matches; omit quartile image page", pair.name)
        return pd.DataFrame()
    merged = pair.matching[["A_ObjectID", "B_ObjectID"]].merge(
        pair.overlaps[["A_ObjectID", "B_ObjectID", "A_Percent_Inside_B", "B_Percent_Inside_A"]],
        on=["A_ObjectID", "B_ObjectID"], how='left', validate='one_to_one')
    if merged['A_Percent_Inside_B'].isna().any() or merged['B_Percent_Inside_A'].isna().any():
        raise ValueError(f'{pair.source}: final matches missing overlap metrics')
    quartiles = quartile_examples(merged, 'A_Percent_Inside_B', ["A_ObjectID", "B_ObjectID"], max_examples=examples_per_quartile)
    selected_rows = []
    for qpack in quartiles:
        q = qpack['quartile']
        for pos, ex in enumerate(qpack['examples'], start=1):
            row = ex['row'].to_dict()
            row.update(quartile=f'Q{q}', quartile_n=qpack['count'], quartile_range_min=qpack['range'][0], quartile_range_max=qpack['range'][1],
                       display_order_within_quartile=pos, displayed_examples_in_quartile=len(qpack['examples']))
            selected_rows.append(row)
    selected = pd.DataFrame(selected_rows)
    if selected.empty:
        return pd.DataFrame()
    measured = intensity_rows(pair, field_index, sample, selected_matches=selected)
    if len(measured) != len(selected):
        raise ValueError('Selected matched pair intensity count mismatch')
    measured = measured.merge(selected, on=['A_ObjectID', 'B_ObjectID'], how='left', validate='one_to_one', suffixes=('', '_selected'))
    bounds_a = image_limits(pair.img_a); bounds_b = image_limits(pair.img_b)
    for q in range(1, 5):
        rows = measured[measured['quartile'] == f'Q{q}'].reset_index(drop=True)
        if rows.empty:
            fig = plt.figure(figsize=(16.8, 4.4), facecolor='white')
            fig.text(.05, .83, f'{pair.name} · Q{q} final stable-matched exemplars', fontsize=16.5, weight='bold')
            fig.text(.05, .60, f'No final stable-matched pairs are present in Q{q}.', fontsize=12.5)
            fig.text(.05, .18, f'Pair source: {pair.source}', fontsize=9)
            pdf.savefig(fig, bbox_inches='tight'); plt.close(fig)
            continue
        for page_start in range(0, len(rows), rows_per_page):
            chunk = rows.iloc[page_start:page_start + rows_per_page].copy()
            page_no = page_start // rows_per_page + 1
            page_total = math.ceil(len(rows) / rows_per_page)
            fig = plt.figure(figsize=(18.2, 2.45 * len(chunk) + 2.0), facecolor='white')
            gs = fig.add_gridspec(len(chunk), 4, left=.042, right=.975, top=.86, bottom=.08, wspace=.12, hspace=.34, width_ratios=[1.04, 1.04, 1.04, 2.9])
            for ridx, r in enumerate(chunk.itertuples(index=False)):
                ax_b = fig.add_subplot(gs[ridx, 0]); ax_a = fig.add_subplot(gs[ridx, 1]); ax_c = fig.add_subplot(gs[ridx, 2]); ax_t = fig.add_subplot(gs[ridx, 3])
                for ax in (ax_b, ax_a, ax_c, ax_t): ax.axis('off')
                aid = int(r.A_ObjectID); bid = int(r.B_ObjectID)
                ap = pair.a_pixels[aid]; bp = pair.b_pixels[bid]; coords = np.vstack((ap, bp))
                x0, x1, y0, y1 = img_crop(pair.img_a, coords, padding=22)
                ca = gray_with_limits(pair.img_a[y0:y1, x0:x1], bounds_a)
                cb = gray_with_limits(pair.img_b[y0:y1, x0:x1], bounds_b)
                am = np.zeros(ca.shape, bool); bm = np.zeros(cb.shape, bool)
                am[ap[:, 1]-y0, ap[:, 0]-x0] = True; bm[bp[:, 1]-y0, bp[:, 0]-x0] = True
                ra = np.repeat(ca[:, :, None], 3, axis=2); rb = np.repeat(cb[:, :, None], 3, axis=2)
                ra[edge(am)] = (.98, .15, .16); rb[edge(bm)] = (.14, .44, 1.)
                composite = np.stack([ca, cb, np.zeros_like(ca)], axis=-1)
                composite[edge(am)] = (.98, .15, .16); composite[edge(bm)] = (.13, .35, 1.)
                for ax, img, title in ((ax_b, rb, f'{pair.b} ROI #{bid}'), (ax_a, ra, f'{pair.a} ROI #{aid}'), (ax_c, composite, 'A red · B green')):
                    ax.imshow(img, interpolation='nearest'); ax.set_title(title, fontsize=8.8, pad=3)
                lines = [
                    f'Q{q} · pair {page_start + ridx + 1} of {len(rows)} displayed · FINAL STABLE MATCH',
                    f'{pair.a} #{aid} ↔ {pair.b} #{bid}    {pair.a} in {pair.b}: {float(r.A_percent_inside_B):.2f}%    {pair.b} in {pair.a}: {float(r.B_percent_inside_A):.2f}%',
                    f'Quartile size: {int(r.quartile_n)}    overlap range in Q{q}: {float(r.quartile_range_min):.2f}–{float(r.quartile_range_max):.2f}%',
                    f'{pair.a}: ROI pixels {int(r.A_roi_area_pixels)}    mean {float(r.A_roi_mean_intensity_native):.5g}    Σ {float(r.A_roi_integrated_intensity_native_sum):.6g}',
                    f'{pair.b}: ROI pixels {int(r.B_roi_area_pixels)}    mean {float(r.B_roi_mean_intensity_native):.5g}    Σ {float(r.B_roi_integrated_intensity_native_sum):.6g}',
                    f'Overlap pixels: {int(r.overlap_area_pixels)}    overlap Σ: {pair.a}={float(r.A_overlap_integrated_intensity_native_sum):.5g}    {pair.b}={float(r.B_overlap_integrated_intensity_native_sum):.5g}',
                ]
                ax_t.text(.01, .99, lines[0], weight='bold', fontsize=10, va='top')
                ax_t.text(.01, .80, "\n".join(lines[1:]), fontsize=9.2, va='top', linespacing=1.4)
            fig.text(.042, .965, f'{pair.name} · Q{q} final stable-matched exemplars', fontsize=16.2, weight='bold')
            fig.text(.042, .928, f'Sample {sample} | Field {field_index} | Matched pairs: {len(pair.matching)} | A-proposing threshold: {pair.threshold:g}% | displaying {len(chunk)} of {len(rows)} requested exemplars (page {page_no}/{page_total})', fontsize=10)
            fig.text(.042, .03, 'Each row is one exact final matched pair only: left = B fluorescence, center-left = A fluorescence, center-right = composite, right = values for these same A/B ROI IDs. Blue outline = B ROI; red outline = A ROI. Intensities are native uncorrected fluorescence.', fontsize=8.6)
            pdf.savefig(fig, bbox_inches='tight'); plt.close(fig)
    return measured


def stable_ids_by_channel(pair_files: list[Path]) -> dict[str,set[int]]:
    """Read ONLY stable-matching IDs for pre-selection, without loading image arrays.

    The final IDs are not recomputed or derived from mere overlap.
    """
    eligible: dict[str,set[int]]=defaultdict(set)
    for path in pair_files:
        with h5py.File(path,'r') as h5:
            meta=json.loads(h5['metadata/json'][()].decode())
            matching=read_table(h5,'stable_matching')
            if matching.empty:continue
            a=meta['objects']['object_a_name'];b=meta['objects']['object_b_name']
            eligible[a].update(matching['A_ObjectID'].astype(int))
            eligible[b].update(matching['B_ObjectID'].astype(int))
    return eligible


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
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--h1-workspace',type=Path,required=True)
    ap.add_argument('--remaining-workspace',type=Path,required=True)
    ap.add_argument('--sample-metadata',type=Path,default=HERE/'igf_exp17_ventral_block_mapping.csv')
    ap.add_argument('--h1-selection',type=Path)
    ap.add_argument('--out',type=Path,help='Default: original workspace/mentor_field_plots/per_organoid_reports_quartiles')
    ap.add_argument('--only-sample',help='Sample-section, e.g. 151-3')
    ap.add_argument('--only-field',help='Source field name substring')
    ap.add_argument('--max-fields',type=int,help='Only this many sorted eligible fields (for testing)')
    ap.add_argument('--roi-metrics',default='integrated_intensity',
                    help='Ranking metric(s) for standalone matched-ROI quartile sheets. Choose: '+','.join(ROI_METRICS))
    ap.add_argument('--max-examples',type=int,default=10,
                    help='Number of displayed exemplar ROIs or matched pairs per quartile (default: 10).')
    ap.add_argument('--skip-single-rois',action='store_true')
    ap.add_argument('--skip-colocalization',action='store_true')
    ap.add_argument('--skip-summary-plots',action='store_true')
    ap.add_argument('--all-summary-plots',action='store_true',help='Include all 20 comparison figures instead of 8 core figures')
    ap.add_argument(
        '--marker-pair-variant',
        choices=('unconstrained', 'nucleus-constrained'),
        default='unconstrained',
        help='Marker-pair result tree to report; variants are written to separate output files',
    )
    ap.add_argument('--compact',action='store_true',help='Backward-compatible alias; concise quartiles are now the default')
    ap.add_argument('--append-existing-qc-pdfs',action='store_true',help='Unsupported in quartile-only reports to avoid unmatched/below-threshold pages')
    ap.add_argument('--overwrite',action='store_true')
    args=ap.parse_args(argv)
    if args.append_existing_qc_pdfs:ap.error('Appending legacy QC PDFs is deliberately disabled: they include non-stable matched examples and very large appendices')
    if args.max_fields is not None and args.max_fields<1:ap.error('--max-fields must be positive')
    if args.max_examples < 1 or args.max_examples > 25:
        ap.error('--max-examples must be 1..25')
    keys=[k.strip() for k in args.roi_metrics.split(',') if k.strip()]
    if not keys or len(set(keys))!=len(keys) or any(k not in ROI_METRICS for k in keys):
        ap.error('Invalid --roi-metrics; allowed: '+','.join(ROI_METRICS))
    script=HERE/'plot_igf_exp17_ventral_with_gfp_plots.py'
    if not script.is_file():ap.error(f'Requires existing sibling script: {script}')
    spec=importlib.util.spec_from_file_location('_reporter_metrics',script)
    if spec is None or spec.loader is None:raise RuntimeError(f'Cannot load {script}')
    plotmod=importlib.util.module_from_spec(spec);spec.loader.exec_module(plotmod)
    original=args.h1_workspace.expanduser().resolve();other=args.remaining_workspace.expanduser().resolve()
    selection=(args.h1_selection or original/'mentor_field_plots/PV_density_per_field_verified_H1.csv')
    fields,_,_=plotmod.collect_fields(original,other,selection,plotmod.load_sample_info(args.sample_metadata.expanduser().resolve()))
    data=pd.DataFrame(fields)
    chosen=data.copy()
    if args.only_sample:chosen=chosen[chosen['sample']==args.only_sample]
    if args.only_field:chosen=chosen[chosen['source_field'].str.contains(args.only_field,regex=False)]
    chosen=chosen.sort_values('index')
    if args.max_fields:chosen=chosen.head(args.max_fields)
    if chosen.empty:raise ValueError('No eligible Experiment 17 ventral fields match requested filters')
    variant_suffix = '_nucleus_constrained' if args.marker_pair_variant == 'nucleus-constrained' else ''
    out=(args.out or original/f'mentor_field_plots/per_organoid_reports_quartiles{variant_suffix}').expanduser().resolve()
    LOGGER.info('Creating %d field reports (final stable matched, Q1..Q4)',len(chosen))
    report_records=[]
    for entry in chosen.to_dict('records'):
        name=re.sub(r'[^A-Za-z0-9_.-]+','_',Path(entry['source_field']).name).strip('_')
        folder=out/entry['sample'];dest=folder/f'{name}_organoid_report_quartiles{variant_suffix}.pdf'
        if dest.exists() and not args.overwrite:
            LOGGER.info('Existing report, skipping: %s (use --overwrite)',dest)
            continue
        folder.mkdir(parents=True,exist_ok=True)
        analysis=Path(entry['output_dir'])
        files=colocalization_files(analysis, args.marker_pair_variant)
        # Read each pair's stored final match identities. No new matching occurs.
        stable=stable_ids_by_channel(files)
        single_rows=[];matched_rows=[];provenance=[]
        partial=dest.with_name(dest.stem+'.incomplete.pdf')
        try:
            with PdfPages(partial) as pdf:
                note_page(pdf,f'Exp 17 · {entry["sample"]} · Field {entry["index"]}',[
                    f'Source: {entry["source_field"]}',
                    f'{entry["cell_line"]} | {entry["condition"]} | {entry["anti_gfp_status"]} | Ventral organoid',
                    f'Marker-pair report variant: {args.marker_pair_variant.upper()}. Stage-07 marker–nuclear pages retain their original definition.',
                    'This PDF represents ONE microscope field, not all fields in its sample-section.',
                    'Final stable matches ONLY: standalone ROI examples must participate in at least one final upstream stable match. Colocalization examples are final one-to-one matched pairs only. No unmatched or below-threshold displays.',
                    'A row of images, IDs and numeric measurements always refers to exactly the SAME selected ROI or A/B PAIR. Each rank quartile Q1–Q4 displays up to the requested number of real exemplars, chosen to span that quartile.',
                    'Native integrated intensity is the sum of ORIGINAL fluorescence pixel values in that exact ROI; measurements are not background or antibody normalized.',
                    'Comparison plots circle this exact field in black. Other dots remain visible. Fields from the same block are not independent biological replicates.',
                    'Final stable-matched IDs eligible for single-channel quartiles: '+', '.join(f'{ch}={len(stable.get(ch,set()))}' for ch in ('GFP','LHX6','PV'))+'. A channel with no eligible final matches is omitted, never populated with unmatched objects.',
                    f'Display policy: up to {args.max_examples} exemplar ROIs or stable-matched pairs per quartile, with multiple pages when needed for readability.',
                    f'Available colocalization analysis files: {len(files)}. Each pair with final matches contributes quartile pages, split automatically when needed.'
                ])
                whole_page(pdf,analysis,provenance)
                if not args.skip_single_rois:
                    for channel in ('GFP','LHX6','PV'):
                        single_rows.extend(roi_quartile_pages(pdf, analysis, channel, keys, stable.get(channel, set()), int(entry['index']), entry['sample'], provenance, examples_per_quartile=args.max_examples))
                if not args.skip_colocalization:
                    for pairfile in files:
                        pair=load_pair(pairfile)
                        provenance.append('Pair HDF5: '+str(pair.source))
                        if not pair.matching.empty:
                            records = matched_pair_quartile_pages(pdf, pair, int(entry['index']), entry['sample'], examples_per_quartile=args.max_examples)
                            if not records.empty:matched_rows.append(records)
                        del pair
                if not args.skip_summary_plots:
                    key_pv={'pv_area_percent_of_organoid_nuclear_only','lhx6_pv_double_positive_density_per_mm2',
                            'lhx6_positive_density_per_mm2','pv_within_lhx6_percent','pv_positive_density_per_mm2'}
                    key_gfp={'gfp_density_per_mm2','gfp_area_percent_of_organoid','gfp_mean_intensity_mean'}
                    for metric,title,unit,_ in plotmod.PLOT_METRICS:
                        if args.all_summary_plots or metric in key_pv:
                            one_metric_page(pdf,data,plotmod,metric,title,unit,int(entry['index']))
                    for metric,title,unit,_ in plotmod.GFP_PLOT_METRICS:
                        if args.all_summary_plots or metric in key_gfp:
                            one_metric_page(pdf,data,plotmod,metric,title,unit,int(entry['index']),gfp=True)
                note_page(pdf,'Measurement key, traceability and source files',[
                    'Q1–Q4 are equal-count quartiles ranked by the named measured feature, or by directional overlap A-inside-B for a matched pair. Quartiles ARE NOT 0–25%, 25–50%, etc. absolute area fractions. If fewer than 4 eligible records exist, remaining quartiles are empty.',
                    'Standalone displays: cyan circle = location of the exact selected ROI in the full aligned field; cyan cropped contour = the same ROI ID. ONLY ROIs in a final stable-matched pair are eligible. Raw pixel integrated intensity is summed over the COMPLETE selected ROI. Up to the requested number of exemplar ROIs are shown within each quartile.',
                    'Colocalization displays: each HORIZONTAL ROW has its own exact A ID and B ID; blue boundary = B ROI, red boundary = A ROI. Composite fluorescence A=red/B=green. Percentages are directional ROI pixel-area overlap (A in B vs B in A), copied from upstream HDF5. Stable matching is not rerun.',
                    'Integrated intensity (Σ) = direct sum of the original stored fluorescence pixels within the FULL matched ROI. Overlap Σ = direct sum of those pixel values at the intersection. No background correction or cross-condition gain normalization. Do not interpret fluorescence units as molecules.',
                    'The CSV exports contain ONLY the displayed quartile exemplars, including the same ROI IDs and quartile labels shown inside this PDF. Summary charts still use the complete previously verified field metrics and highlight THIS field in black. If a quartile contains fewer than the requested examples, all available members are shown.',
                ]+provenance[:12])
            partial.replace(dest)
        except Exception:
            partial.unlink(missing_ok=True)
            raise
        solo_csv=folder/f'{name}_displayed_single_roi_quartiles{variant_suffix}.csv'
        matched_csv=folder/f'{name}_displayed_stable_matched_pair_quartiles{variant_suffix}.csv'
        pd.DataFrame(single_rows, columns=['field_index', 'sample', 'channel', 'quartile', 'rank_metric', 'rank_value', 'quartile_n', 'quartile_range_min', 'quartile_range_max', 'display_order_within_quartile', 'displayed_examples_in_quartile', 'ROI_ObjectID', 'roi_pixel_count', 'Area_um2', 'mean_native_intensity', 'integrated_native_intensity_sum', 'Circularity', 'Solidity', 'Eccentricity']).to_csv(solo_csv, index=False)
        matched_frame = pd.concat(matched_rows, ignore_index=True) if matched_rows else pd.DataFrame(
            columns=['field_index','sample','pair','A_ObjectID','B_ObjectID','quartile','display_order_within_quartile','displayed_examples_in_quartile','A_roi_integrated_intensity_native_sum','B_roi_integrated_intensity_native_sum','A_overlap_integrated_intensity_native_sum','B_overlap_integrated_intensity_native_sum']
        )
        matched_frame.to_csv(matched_csv,index=False)
        report_records.append({'index':entry['index'],'sample':entry['sample'],'source_field':entry['source_field'],
            'marker_pair_variant':args.marker_pair_variant,
            'pdf':str(dest),'single_roi_displayed_csv':str(solo_csv),'stable_matched_displayed_csv':str(matched_csv),
            'n_pair_hdf5':len(files),'displayed_single_roi':len(single_rows),'displayed_final_pairs':len(matched_frame)})
        LOGGER.info('Saved %s | standalone ROI exemplars=%d | final matched pair exemplars=%d',dest,len(single_rows),len(matched_frame))
    if report_records:
        out.mkdir(parents=True,exist_ok=True)
        pd.DataFrame(report_records).to_csv(out/'report_manifest.csv',index=False)
    return 0


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    try:
        sys.exit(main())
    except (ValueError, FileNotFoundError, KeyError) as exc:
        sys.exit(f"ERROR: {exc}")
