#!/usr/bin/env python3
"""Summarize aligned CellProfiler organoid masks and containment failures."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import tifffile
from scipy import ndimage


STAT_COLUMNS = (
    "channel",
    "height_pixels",
    "width_pixels",
    "organoid_mask_area_pixels",
    "connected_components",
    "largest_component_area_pixels",
    "touches_image_border",
    "segmented_cell_pixels",
    "cell_pixels_outside_organoid",
    "cell_pixels_outside_organoid_percent",
)


def read_2d(path: Path, label: str) -> np.ndarray:
    array = np.squeeze(tifffile.imread(path))
    if array.ndim != 2:
        raise ValueError(f"{label} must be a 2-D TIFF, got {array.shape}: {path}")
    return array


def normalize_fluorescence(image: np.ndarray) -> np.ndarray:
    finite = np.asarray(image, dtype=np.float32)
    valid = finite[np.isfinite(finite)]
    if valid.size == 0:
        return np.zeros(finite.shape, dtype=np.float32)
    low, high = np.percentile(valid, (1, 99.8))
    if high <= low:
        low = float(valid.min())
        high = float(valid.max())
    if high <= low:
        return np.zeros(finite.shape, dtype=np.float32)
    return np.clip((finite - low) / (high - low), 0, 1)


def mask_statistics(
    channel: str, organoid_mask: np.ndarray, cell_mask: np.ndarray
) -> tuple[dict[str, object], np.ndarray]:
    organoid = organoid_mask > 0
    cells = cell_mask > 0
    if organoid.shape != cells.shape:
        raise ValueError(
            f"{channel}: organoid and cell masks have different shapes "
            f"({organoid.shape} and {cells.shape})"
        )
    labels, component_count = ndimage.label(organoid)
    component_areas = np.bincount(labels.ravel())[1:]
    largest = int(component_areas.max()) if component_areas.size else 0
    outside = cells & ~organoid
    cell_pixels = int(cells.sum())
    outside_pixels = int(outside.sum())
    touches_border = bool(
        organoid[0, :].any()
        or organoid[-1, :].any()
        or organoid[:, 0].any()
        or organoid[:, -1].any()
    )
    statistics: dict[str, object] = {
        "channel": channel,
        "height_pixels": organoid.shape[0],
        "width_pixels": organoid.shape[1],
        "organoid_mask_area_pixels": int(organoid.sum()),
        "connected_components": int(component_count),
        "largest_component_area_pixels": largest,
        "touches_image_border": touches_border,
        "segmented_cell_pixels": cell_pixels,
        "cell_pixels_outside_organoid": outside_pixels,
        "cell_pixels_outside_organoid_percent": (
            100.0 * outside_pixels / cell_pixels if cell_pixels else 0.0
        ),
    }
    return statistics, outside


def run_qc(
    entries: list[tuple[str, Path, Path, Path]], output_dir: Path
) -> tuple[Path, Path]:
    if not entries:
        raise ValueError("At least one channel entry is required")
    output_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = []
    panels: list[tuple[str, np.ndarray, np.ndarray, np.ndarray]] = []

    for channel, fluorescence_path, organoid_path, cell_path in entries:
        fluorescence = read_2d(fluorescence_path, f"{channel} fluorescence")
        organoid = read_2d(organoid_path, f"{channel} organoid mask") > 0
        cells = read_2d(cell_path, f"{channel} cell mask") > 0
        if fluorescence.shape != organoid.shape or organoid.shape != cells.shape:
            raise ValueError(
                f"{channel}: fluorescence, organoid mask, and cell mask must share "
                f"one shape; got {fluorescence.shape}, {organoid.shape}, {cells.shape}"
            )
        statistics, outside = mask_statistics(channel, organoid, cells)
        rows.append(statistics)
        panels.append((channel, fluorescence, organoid, outside))

    csv_path = output_dir / "organoid_mask_statistics.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=STAT_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)

    figure, axes = plt.subplots(
        len(panels),
        3,
        figsize=(15, 4.5 * len(panels)),
        squeeze=False,
        constrained_layout=True,
    )
    for row_index, (channel, fluorescence, organoid, outside) in enumerate(panels):
        normalized = normalize_fluorescence(fluorescence)
        boundary = organoid & ~ndimage.binary_erosion(organoid)
        overlay = np.repeat(normalized[..., None], 3, axis=2)
        overlay[boundary] = (0.0, 1.0, 0.0)
        overlay[outside] = (1.0, 0.0, 0.0)

        axes[row_index, 0].imshow(normalized, cmap="gray", vmin=0, vmax=1)
        axes[row_index, 0].set_title(f"{channel}: fluorescence")
        axes[row_index, 1].imshow(organoid, cmap="gray", vmin=0, vmax=1)
        axes[row_index, 1].set_title(f"{channel}: OrganoidMask.tiff")
        axes[row_index, 2].imshow(overlay, vmin=0, vmax=1)
        axes[row_index, 2].set_title(f"{channel}: green boundary, red cells outside")
        for axis in axes[row_index]:
            axis.set_axis_off()

    figure_path = output_dir / "organoid_mask_qc.png"
    figure.savefig(figure_path, dpi=160)
    plt.close(figure)

    for row in rows:
        print(
            f"{row['channel']}: mask_area={row['organoid_mask_area_pixels']}, "
            f"components={row['connected_components']}, "
            f"largest={row['largest_component_area_pixels']}, "
            f"touches_border={row['touches_image_border']}, "
            f"outside_cell_pixels={row['cell_pixels_outside_organoid']} "
            f"({float(row['cell_pixels_outside_organoid_percent']):.4f}%)"
        )
    print(f"Wrote {csv_path}")
    print(f"Wrote {figure_path}")
    return figure_path, csv_path


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--entry",
        nargs=4,
        action="append",
        metavar=("CHANNEL", "FLUORESCENCE", "ORGANOID_MASK", "CELL_MASK"),
        required=True,
        help="Repeat once per aligned fluorescence channel.",
    )
    args = parser.parse_args(argv)
    entries = [
        (channel, Path(fluorescence), Path(organoid), Path(cells))
        for channel, fluorescence, organoid, cells in args.entry
    ]
    run_qc(entries, args.output_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
