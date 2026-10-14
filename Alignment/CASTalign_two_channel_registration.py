#!/usr/bin/env python3
"""
Headless 2-D max-projection registration with CellProfiler centroids + CASTalign.

Designed for two fluorescence channels from the SAME organoid, e.g.:
    fixed  = Hoechst max projection
    moving = PV max projection

Workflow
--------
1. Load the two 2-D max-projection TIFFs.
2. Read CellProfiler object centroids (Location_Center_X/Y).
3. Optionally run image-based phase correlation only as a coarse XY initializer.
4. Globally assign a subset of centroids one-to-one using a partial Hungarian
   assignment, allowing unmatched objects.
5. Fit a CASTalign point-based Rigid transform to the accepted anchors.
6. Iterate transform -> rematch -> robust outlier rejection -> refit.
7. Optionally fit an Affine CASTalign transform as a comparison, but keep Rigid
   as the default unless affine clearly improves the anchor residuals.
8. Apply the final XY transform to the moving max projection.
9. Save the transform, anchor tables, registered image, overlay/QC images, and
   a text summary.

Important scientific note
-------------------------
The centroid pairs are registration anchors. They are NOT assertions that a
Hoechst nucleus and a PV-positive object are the same biological cell.

This script operates ONLY on 2-D max projections. It does not alter or resample
your original 3-D stacks. CASTalign's core API is 3-D volume-coordinate based, so
2-D anchor points are represented internally as (z=0, y, x), and the resulting
rigid transform is constrained by those coplanar anchors to an in-plane XY motion.

Install
-------
    pip install castalign tifffile pandas numpy scipy scikit-image matplotlib

CASTalign's current official repository recommends `pip install castalign` and
optionally CuPy for GPU acceleration. The registration in this script is small
and 2-D, so CPU execution is generally sufficient for the first test.

Example
-------
python CASTalign_2D_maxproj_registration.py \
    --fixed-image /path/Hoechst_max.tif \
    --moving-image /path/PV_max.tif \
    --fixed-csv /path/Hoechst_objects.csv \
    --moving-csv /path/PV_objects.csv \
    --output /path/CASTalign_Hoechst_PV

If you want an affine comparison on the first test:
    --try-affine

The default physical scale is the XY pixel size from your IMS data:
    X = 0.312843137 um
    Y = 0.312836345 um
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Optional, Tuple

import numpy as np
import pandas as pd
import tifffile
from scipy.optimize import linear_sum_assignment
from skimage.registration import phase_cross_correlation

import castalign as ca


def log(msg: str) -> None:
    print(msg, flush=True)


def load_tiff_2d(path: Path) -> np.ndarray:
    arr = np.squeeze(tifffile.imread(path))
    if arr.ndim != 2:
        raise ValueError(
            f"{path} has shape {arr.shape}. This script expects a 2-D max projection TIFF."
        )
    return arr


def find_xy_columns(df: pd.DataFrame, x_override: Optional[str], y_override: Optional[str]) -> Tuple[str, str]:
    if x_override:
        if x_override not in df.columns:
            raise ValueError(f"X column '{x_override}' not found in CSV.")
        x_col = x_override
    else:
        candidates = [
            "Location_Center_X",
            "Location_Center_X_(px)",
            "Location_Center_X_(pixels)",
        ]
        x_col = next((c for c in candidates if c in df.columns), None)
        if x_col is None:
            fallback = [c for c in df.columns if "location_center" in str(c).lower() and str(c).lower().endswith("_x")]
            if len(fallback) != 1:
                raise ValueError(
                    "Could not identify the CellProfiler X centroid column. "
                    f"Available columns include: {list(df.columns[:50])}"
                )
            x_col = fallback[0]

    if y_override:
        if y_override not in df.columns:
            raise ValueError(f"Y column '{y_override}' not found in CSV.")
        y_col = y_override
    else:
        candidates = [
            "Location_Center_Y",
            "Location_Center_Y_(px)",
            "Location_Center_Y_(pixels)",
        ]
        y_col = next((c for c in candidates if c in df.columns), None)
        if y_col is None:
            fallback = [c for c in df.columns if "location_center" in str(c).lower() and str(c).lower().endswith("_y")]
            if len(fallback) != 1:
                raise ValueError(
                    "Could not identify the CellProfiler Y centroid column. "
                    f"Available columns include: {list(df.columns[:50])}"
                )
            y_col = fallback[0]

    return x_col, y_col


def read_centroids(csv_path: Path, x_override: Optional[str], y_override: Optional[str]) -> Tuple[np.ndarray, pd.DataFrame, Tuple[str, str]]:
    df = pd.read_csv(csv_path)
    x_col, y_col = find_xy_columns(df, x_override, y_override)

    x = pd.to_numeric(df[x_col], errors="coerce").to_numpy(dtype=float)
    y = pd.to_numeric(df[y_col], errors="coerce").to_numpy(dtype=float)
    points = np.column_stack([y, x])  # CASTalign point convention for 2-D image: (y, x)

    good = np.all(np.isfinite(points), axis=1)
    points = points[good]
    clean_df = df.loc[good].reset_index(drop=True)

    if len(points) < 3:
        raise ValueError(f"{csv_path}: fewer than 3 valid centroids.")

    return points, clean_df, (x_col, y_col)


def partial_hungarian_match(moving: np.ndarray, fixed: np.ndarray, gate_px: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Globally assign moving -> fixed with one-to-one matching and explicit
    unmatched-object handling.

    Any pair farther than gate_px is treated as invalid. Both unmatched moving
    and unmatched fixed points are allowed.
    """
    n_m = len(moving)
    n_f = len(fixed)
    if n_m == 0 or n_f == 0:
        return np.empty(0, int), np.empty(0, int), np.empty(0, float)

    diff = moving[:, None, :] - fixed[None, :, :]
    dist2 = np.einsum("ijk,ijk->ij", diff, diff)
    gate2 = gate_px * gate_px

    # Dense assignment is intentionally used for a preliminary script because
    # it is easy to audit. --max-points protects against very large CSVs.
    dummy_cost = gate2 * 1.01
    invalid_cost = dummy_cost * 4.0

    size = n_m + n_f
    cost = np.full((size, size), dummy_cost, dtype=float)
    real = dist2.copy()
    real[real > gate2] = invalid_cost
    cost[:n_m, :n_f] = real
    cost[n_m:, n_f:] = 0.0

    rows, cols = linear_sum_assignment(cost)

    mi, fi, d = [], [], []
    for r, c in zip(rows, cols):
        if r < n_m and c < n_f:
            distance = math.sqrt(dist2[r, c])
            if distance <= gate_px:
                mi.append(r)
                fi.append(c)
                d.append(distance)

    return np.asarray(mi, int), np.asarray(fi, int), np.asarray(d, float)


def spatial_subsample(points: np.ndarray, max_points: int, seed: int) -> np.ndarray:
    """Keep a spatially distributed subset rather than a purely random one."""
    if max_points <= 0 or len(points) <= max_points:
        return np.arange(len(points))

    rng = np.random.default_rng(seed)
    n_side = max(2, int(round(max_points ** 0.5)))
    mins = np.min(points, axis=0)
    maxs = np.max(points, axis=0)
    span = np.maximum(maxs - mins, 1e-9)

    ij = np.floor((points - mins) / span * n_side).astype(int)
    ij = np.clip(ij, 0, n_side - 1)
    keys = ij[:, 0] * n_side + ij[:, 1]
    _, first = np.unique(keys, return_index=True)

    selected = list(first)
    if len(selected) > max_points:
        selected = rng.choice(selected, size=max_points, replace=False).tolist()
    else:
        remaining = np.setdiff1d(np.arange(len(points)), np.asarray(selected), assume_unique=False)
        n_extra = min(max_points - len(selected), len(remaining))
        if n_extra:
            selected.extend(rng.choice(remaining, size=n_extra, replace=False).tolist())

    return np.asarray(selected, dtype=int)


def phase_correlation_shift(fixed: np.ndarray, moving: np.ndarray, upsample_factor: int = 10) -> Tuple[np.ndarray, float, float]:
    """Return moving-to-fixed coarse shift; tolerate different 2-D image shapes.

    Images are padded on bottom/right *temporarily* only to make their shapes
    compatible with phase_cross_correlation. Their (0,0) pixel origins stay
    unchanged; no registered output or source TIFF is padded or resampled here.
    The original identical-shape code path is retained exactly.
    """
    f = fixed.astype(np.float32, copy=False)
    m = moving.astype(np.float32, copy=False)
    if f.shape != m.shape:
        common_y = max(f.shape[0], m.shape[0])
        common_x = max(f.shape[1], m.shape[1])
        temporary_fixed = np.zeros((common_y, common_x), dtype=np.float32)
        temporary_moving = np.zeros((common_y, common_x), dtype=np.float32)
        temporary_fixed[:f.shape[0], :f.shape[1]] = f
        temporary_moving[:m.shape[0], :m.shape[1]] = m
        f, m = temporary_fixed, temporary_moving

    shift, error, phasediff = phase_cross_correlation(
        f, m, upsample_factor=upsample_factor
    )
    return np.asarray(shift, dtype=float), float(error), float(phasediff)


def robust_cutoff(residuals: np.ndarray, minimum: float, multiplier: float) -> float:
    if len(residuals) == 0:
        return minimum
    med = float(np.median(residuals))
    mad = float(np.median(np.abs(residuals - med)))
    robust_sigma = 1.4826 * mad
    return max(minimum, med + multiplier * robust_sigma)


def lift_points_to_castalign(points_2d: np.ndarray) -> np.ndarray:
    """Represent 2-D (Y,X) points in CASTalign's required (Z,Y,X) space."""
    points_2d = np.asarray(points_2d, dtype=float)
    if points_2d.ndim != 2 or points_2d.shape[1] != 2:
        raise ValueError(f"Expected points with shape (N,2), got {points_2d.shape}")
    z = np.zeros((len(points_2d), 1), dtype=float)
    return np.hstack([z, points_2d])


def drop_z_from_castalign(points_3d: np.ndarray) -> np.ndarray:
    """Convert CASTalign (Z,Y,X) points back to 2-D (Y,X)."""
    arr = np.asarray(points_3d, dtype=float)
    if arr.ndim != 2 or arr.shape[1] != 3:
        raise ValueError(f"Expected points with shape (N,3), got {arr.shape}")
    return arr[:, 1:3]


def transform_2d_points_with_castalign(transform, points_2d: np.ndarray) -> np.ndarray:
    """Apply a CASTalign transform to 2-D points via the z=0 plane."""
    return drop_z_from_castalign(transform.transform(lift_points_to_castalign(points_2d)))


def fit_rigid_iteratively(
    fixed_points: np.ndarray,
    moving_points: np.ndarray,
    initial_shift: np.ndarray,
    gate_px: float,
    iterations: int,
    min_anchors: int,
    outlier_mad: float,
) -> dict:
    """
    Iterative CASTalign Rigid registration for 2-D max projections.

    CASTalign's Transform.transform() requires 3-D (Z,Y,X) coordinates.
    Because this is a true 2-D registration, all anchors are embedded in the
    z=0 plane. The fitted rigid transform is therefore constrained by the
    coplanar source/target points to the in-plane XY geometry we want.
    """
    current_moving = moving_points + initial_shift
    fixed_cast = lift_points_to_castalign(fixed_points)
    final = None

    for it in range(1, iterations + 1):
        mi, fi, _ = partial_hungarian_match(current_moving, fixed_points, gate_px)
        if len(mi) < min_anchors:
            raise RuntimeError(
                f"Iteration {it}: only {len(mi)} anchors found; need {min_anchors}. "
                f"Try increasing --match-gate-um or checking the CellProfiler CSVs."
            )

        moving_cast = lift_points_to_castalign(current_moving[mi])
        fixed_cast_matches = fixed_cast[fi]
        t = ca.Rigid(points_start=moving_cast, points_end=fixed_cast_matches)
        transformed_matches_3d = t.transform(moving_cast)
        transformed_matches = drop_z_from_castalign(transformed_matches_3d)
        residuals = np.linalg.norm(transformed_matches - fixed_points[fi], axis=1)

        cutoff = min(gate_px, robust_cutoff(residuals, gate_px * 0.25, outlier_mad))
        keep = residuals <= cutoff
        if np.count_nonzero(keep) < min_anchors:
            keep = np.ones_like(residuals, dtype=bool)

        mi_keep = mi[keep]
        fi_keep = fi[keep]
        t_refit = ca.Rigid(
            points_start=lift_points_to_castalign(current_moving[mi_keep]),
            points_end=lift_points_to_castalign(fixed_points[fi_keep]),
        )

        transformed_all_3d = t_refit.transform(lift_points_to_castalign(current_moving))
        transformed_all = drop_z_from_castalign(transformed_all_3d)
        final_matches = drop_z_from_castalign(
            t_refit.transform(lift_points_to_castalign(current_moving[mi_keep]))
        )
        final_residuals = np.linalg.norm(
            final_matches - fixed_points[fi_keep], axis=1
        )

        # Verify that CASTalign is effectively preserving the z=0 plane.
        z_values = transformed_all_3d[:, 0]
        if np.max(np.abs(z_values)) > 1e-6:
            raise RuntimeError(
                "CASTalign Rigid produced a non-zero Z component from coplanar "
                f"2-D anchors (max |Z|={np.max(np.abs(z_values)):.3g}). "
                "This preliminary 2-D implementation will not continue because "
                "an out-of-plane transform would violate the intended workflow."
            )

        final = {
            "transform": t_refit,
            "moving_indices": mi_keep,
            "fixed_indices": fi_keep,
            "transformed_all": transformed_all,
            "transformed_all_3d": transformed_all_3d,
            "residuals": final_residuals,
            "iteration": it,
            "cutoff": cutoff,
        }

        log(
            f"  Iteration {it}: anchors={len(mi_keep)}, "
            f"median residual={np.median(final_residuals):.3f} px, "
            f"95th percentile={np.percentile(final_residuals, 95):.3f} px"
        )
        current_moving = transformed_all

    if final is None:
        raise RuntimeError("Rigid fitting did not produce a result.")
    return final


def save_anchor_csv(
    path: Path,
    moving_orig: np.ndarray,
    fixed_orig: np.ndarray,
    transformed_moving: np.ndarray,
    residual_px: np.ndarray,
    pixel_x_um: float,
    pixel_y_um: float,
) -> None:
    residual_um = residual_px * math.sqrt(pixel_x_um * pixel_y_um)
    df = pd.DataFrame(
        {
            "moving_y_px": moving_orig[:, 0],
            "moving_x_px": moving_orig[:, 1],
            "fixed_y_px": fixed_orig[:, 0],
            "fixed_x_px": fixed_orig[:, 1],
            "transformed_moving_y_px": transformed_moving[:, 0],
            "transformed_moving_x_px": transformed_moving[:, 1],
            "residual_px": residual_px,
            "residual_um_approx": residual_um,
        }
    )
    df.to_csv(path, index=False)


def save_image_qc(fixed: np.ndarray, moving_reg: np.ndarray, output: Path) -> None:
    import matplotlib.pyplot as plt

    def norm(a):
        a = a.astype(np.float32)
        lo, hi = np.percentile(a[np.isfinite(a)], [1, 99.5])
        if hi <= lo:
            hi = lo + 1.0
        return np.clip((a - lo) / (hi - lo), 0, 1)

    f = norm(fixed)
    m = norm(moving_reg)

    # Save numeric two-channel overlay data.
    tifffile.imwrite(
        output / "aligned_two_channel.tif",
        np.stack([fixed.astype(np.float32), moving_reg.astype(np.float32)], axis=0),
        bigtiff=True,
        compression="zlib",
    )

    # Simple RGB visualization: fixed -> green, moving -> red.
    rgb = np.zeros((*f.shape, 3), dtype=np.float32)
    rgb[..., 1] = f
    rgb[..., 0] = m

    plt.figure(figsize=(9, 9))
    plt.imshow(rgb)
    plt.title("Hoechst (green) + registered moving channel (red)")
    plt.axis("off")
    plt.tight_layout(pad=0)
    plt.savefig(output / "registration_overlay.png", dpi=300, bbox_inches="tight", pad_inches=0)
    plt.close()

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    axes[0].imshow(f, cmap="gray")
    axes[0].set_title("Fixed / Hoechst")
    axes[1].imshow(m, cmap="gray")
    axes[1].set_title("Registered moving")
    axes[2].imshow(rgb)
    axes[2].set_title("Overlay")
    for ax in axes:
        ax.axis("off")
    fig.tight_layout()
    fig.savefig(output / "registration_qc.png", dpi=250)
    plt.close(fig)


def save_anchor_visualization(
    fixed: np.ndarray,
    moving: np.ndarray,
    moving_reg: np.ndarray,
    moving_anchor_orig: np.ndarray,
    fixed_anchor: np.ndarray,
    transformed_anchor: np.ndarray,
    output: Path,
) -> None:
    import matplotlib.pyplot as plt

    def norm(a):
        a = a.astype(np.float32)
        vals = a[np.isfinite(a)]
        lo, hi = np.percentile(vals, [1, 99.5])
        if hi <= lo:
            hi = lo + 1.0
        return np.clip((a - lo) / (hi - lo), 0, 1)

    f = norm(fixed)
    m = norm(moving)
    mr = norm(moving_reg)

    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    axes[0].imshow(m, cmap="gray")
    axes[0].scatter(moving_anchor_orig[:, 1], moving_anchor_orig[:, 0], s=8)
    axes[0].set_title("Moving + chosen anchors")

    axes[1].imshow(f, cmap="gray")
    axes[1].scatter(fixed_anchor[:, 1], fixed_anchor[:, 0], s=8)
    axes[1].set_title("Fixed + corresponding anchors")

    axes[2].imshow(mr, cmap="gray")
    axes[2].imshow(f, alpha=0.35, cmap="gray")
    axes[2].set_title("Registered moving / fixed")

    for ax in axes:
        ax.axis("off")
    fig.tight_layout()
    fig.savefig(output / "anchor_qc.png", dpi=250)
    plt.close(fig)


def save_residual_histogram(residual_um: np.ndarray, output: Path) -> None:
    import matplotlib.pyplot as plt
    plt.figure(figsize=(7, 5))
    plt.hist(residual_um, bins=30)
    plt.xlabel("Anchor residual (µm, approximate)")
    plt.ylabel("Anchor count")
    plt.title("CASTalign rigid-registration residuals")
    plt.tight_layout()
    plt.savefig(output / "anchor_residual_histogram.png", dpi=250)
    plt.close()


def main() -> int:
    ap = argparse.ArgumentParser(description="Headless 2-D max-projection CASTalign registration")
    ap.add_argument("--fixed-image", required=True, type=Path, help="Fixed/reference max-projection TIFF, e.g. Hoechst")
    ap.add_argument("--moving-image", required=True, type=Path, help="Moving max-projection TIFF, e.g. PV")
    ap.add_argument("--fixed-csv", required=True, type=Path, help="CellProfiler fixed-channel object CSV")
    ap.add_argument("--moving-csv", required=True, type=Path, help="CellProfiler moving-channel object CSV")
    ap.add_argument("--output", required=True, type=Path, help="Output directory")

    ap.add_argument("--pixel-size-x", type=float, default=0.312843137)
    ap.add_argument("--pixel-size-y", type=float, default=0.312836345)

    ap.add_argument("--fixed-x-col", default=None)
    ap.add_argument("--fixed-y-col", default=None)
    ap.add_argument("--moving-x-col", default=None)
    ap.add_argument("--moving-y-col", default=None)

    ap.add_argument("--skip-phase-init", action="store_true", help="Skip coarse phase-correlation initialization")
    ap.add_argument("--phase-upsample", type=int, default=10)

    ap.add_argument("--match-gate-um", type=float, default=25.0,
                    help="Maximum centroid matching distance in microns (default 25)")
    ap.add_argument("--iterations", type=int, default=5)
    ap.add_argument("--min-anchors", type=int, default=8)
    ap.add_argument("--outlier-mad", type=float, default=3.5)
    ap.add_argument("--max-points", type=int, default=1500,
                    help="Maximum centroids per channel used in matching. 0 = all.")
    ap.add_argument("--seed", type=int, default=42)

    ap.add_argument("--try-affine", action="store_true",
                    help="Fit an affine CASTalign transform after rigid registration for comparison")
    ap.add_argument("--affine-min-improvement", type=float, default=0.15,
                    help="Minimum median-residual improvement needed to select affine (default 15%%)")
    ap.add_argument("--force-affine", action="store_true",
                    help="Select affine if it fits successfully; use cautiously")

    args = ap.parse_args()
    out = args.output
    out.mkdir(parents=True, exist_ok=True)

    # ---------------------------------------------------------------
    # Load images
    # ---------------------------------------------------------------
    log("Loading max projections...")
    fixed = load_tiff_2d(args.fixed_image)
    moving = load_tiff_2d(args.moving_image)
    log(f"  Fixed shape : {fixed.shape}")
    log(f"  Moving shape: {moving.shape}")

    if fixed.shape != moving.shape:
        log(
            "  Different native XY shapes are permitted; coarse phase correlation "
            "uses temporary zero-padded arrays, and the final CASTalign output "
            "keeps the fixed image canvas."
        )

    # ---------------------------------------------------------------
    # CellProfiler centroids
    # ---------------------------------------------------------------
    log("Reading CellProfiler centroids...")
    fixed_points, fixed_df, fixed_cols = read_centroids(args.fixed_csv, args.fixed_x_col, args.fixed_y_col)
    moving_points, moving_df, moving_cols = read_centroids(args.moving_csv, args.moving_x_col, args.moving_y_col)
    log(f"  Fixed centroids : {len(fixed_points)}")
    log(f"  Moving centroids: {len(moving_points)}")
    log(f"  Fixed columns   : X={fixed_cols[0]}, Y={fixed_cols[1]}")
    log(f"  Moving columns  : X={moving_cols[0]}, Y={moving_cols[1]}")

    # Confirm that the CellProfiler coordinates actually refer to these max projections.
    for label, pts, shape in (("fixed", fixed_points, fixed.shape), ("moving", moving_points, moving.shape)):
        bad = (pts[:, 0] < 0) | (pts[:, 0] >= shape[0]) | (pts[:, 1] < 0) | (pts[:, 1] >= shape[1])
        if np.any(bad):
            raise ValueError(
                f"{label} CellProfiler CSV contains {int(np.sum(bad))} centroid(s) outside "
                f"the corresponding image bounds {shape}. This usually means the CSV and TIFF "
                "are not the exact same max projection/crop."
            )

    fi = spatial_subsample(fixed_points, args.max_points, args.seed)
    mi = spatial_subsample(moving_points, args.max_points, args.seed + 1)
    fixed_use = fixed_points[fi]
    moving_use = moving_points[mi]
    if len(fi) < len(fixed_points) or len(mi) < len(moving_points):
        log(f"  Spatial subset used: fixed={len(fi)}, moving={len(mi)}")

    # ---------------------------------------------------------------
    # Coarse image-based initialization
    # ---------------------------------------------------------------
    shift_px = np.zeros(2, dtype=float)
    phase_error = float("nan")
    phase_diff = float("nan")

    if not args.skip_phase_init:
        log("Running coarse 2-D phase correlation (initializer only)...")
        shift_px, phase_error, phase_diff = phase_correlation_shift(
            fixed, moving, upsample_factor=args.phase_upsample
        )
        log(f"  Shift to apply to moving (Y,X px): {shift_px}")
        log(f"  Correlation error: {phase_error:.6g}")

    match_gate_px = args.match_gate_um / math.sqrt(args.pixel_size_x * args.pixel_size_y)
    log(f"Centroid matching gate: {args.match_gate_um:.2f} µm ≈ {match_gate_px:.2f} px")

    initial_moving = moving_use + shift_px
    mi0, fi0, d0 = partial_hungarian_match(initial_moving, fixed_use, match_gate_px)
    if len(mi0) < args.min_anchors:
        raise RuntimeError(
            f"Only {len(mi0)} initial anchors found; need {args.min_anchors}. "
            "Try increasing --match-gate-um, checking the CellProfiler CSVs, or using --skip-phase-init."
        )
    log(f"Initial global centroid anchors: {len(mi0)}")

    initial_anchor_df = pd.DataFrame({
        "moving_y_px": moving_use[mi0, 0],
        "moving_x_px": moving_use[mi0, 1],
        "fixed_y_px": fixed_use[fi0, 0],
        "fixed_x_px": fixed_use[fi0, 1],
        "initial_distance_px": d0,
        "initial_distance_um_approx": d0 * math.sqrt(args.pixel_size_x * args.pixel_size_y),
    })
    initial_anchor_df.to_csv(out / "matched_centroids_initial.csv", index=False)

    # ---------------------------------------------------------------
    # Iterative CASTalign rigid
    # ---------------------------------------------------------------
    log("Fitting iterative CASTalign Rigid transform...")
    rigid = fit_rigid_iteratively(
        fixed_points=fixed_use,
        moving_points=moving_use,
        initial_shift=shift_px,
        gate_px=match_gate_px,
        iterations=args.iterations,
        min_anchors=args.min_anchors,
        outlier_mad=args.outlier_mad,
    )

    rigid_t = rigid["transform"]
    rmi = rigid["moving_indices"]
    rfi = rigid["fixed_indices"]
    rigid_residual_px = rigid["residuals"]

    # Original moving/fixed points corresponding to selected subsets.
    moving_anchor_orig = moving_use[rmi]
    fixed_anchor_orig = fixed_use[rfi]
    moving_anchor_transformed = transform_2d_points_with_castalign(rigid_t, moving_anchor_orig + shift_px)

    residual_um = rigid_residual_px * math.sqrt(args.pixel_size_x * args.pixel_size_y)
    save_anchor_csv(
        out / "matched_centroids_rigid.csv",
        moving_anchor_orig,
        fixed_anchor_orig,
        moving_anchor_transformed,
        rigid_residual_px,
        args.pixel_size_x,
        args.pixel_size_y,
    )

    # Save reconstructible CASTalign rigid transform.
    rigid_t.save(out / "transform_rigid.txt")

    # ---------------------------------------------------------------
    # Optional affine comparison
    # ---------------------------------------------------------------
    chosen_model = "Rigid"
    final_t = rigid_t
    affine_t = None
    affine_residual_px = None
    affine_improvement = None

    if args.try_affine:
        log("Fitting optional CASTalign Affine comparison...")
        # Fit affine to the rigidly transformed anchor positions.
        rigid_anchor_positions = transform_2d_points_with_castalign(
            rigid_t, moving_anchor_orig + shift_px
        )
        affine_t = ca.LaminarAffine(
            points_start=lift_points_to_castalign(rigid_anchor_positions),
            points_end=lift_points_to_castalign(fixed_anchor_orig),
        )
        affine_anchor_positions = transform_2d_points_with_castalign(
            affine_t, rigid_anchor_positions
        )
        affine_residual_px = np.linalg.norm(
            affine_anchor_positions - fixed_anchor_orig, axis=1
        )

        rigid_median = float(np.median(rigid_residual_px))
        affine_median = float(np.median(affine_residual_px))
        affine_improvement = 0.0 if rigid_median <= 0 else (rigid_median - affine_median) / rigid_median

        log(f"  Rigid median residual : {rigid_median:.3f} px")
        log(f"  Affine median residual: {affine_median:.3f} px")
        log(f"  Improvement           : {100 * affine_improvement:.2f}%")

        if args.force_affine or affine_improvement >= args.affine_min_improvement:
            chosen_model = "Affine"
            final_t = rigid_t + affine_t
            affine_t.save(out / "transform_affine_after_rigid.txt")

    # Final transform including the initial phase-correlation translation.
    # Translation is applied first to the moving image coordinates, then the
    # CASTalign fitted transform acts on the translated coordinates.
    # CASTalign's 2-D Rigid transform is therefore composed after translation.
    initializer_t = ca.Translate(points_start=lift_points_to_castalign(np.zeros((1, 2))), points_end=lift_points_to_castalign(np.array([[shift_px[0], shift_px[1]]])))
    full_final_t = initializer_t + final_t
    full_final_t.save(out / "transform_final.txt")

    # ---------------------------------------------------------------
    # Apply to the max projection
    # ---------------------------------------------------------------
    log(f"Applying final {chosen_model} transform to moving max projection...")
    registered = full_final_t.transform_image(
        moving,
        output_size=(1, fixed.shape[0], fixed.shape[1]),
        force_size=True,
        labels=False,
    )
    registered = np.asarray(registered, dtype=np.float32)
    if registered.ndim == 3 and registered.shape[0] == 1:
        registered = registered[0]

    # Save the two channels as separate TIFFs in the SAME final coordinate frame.
    # Hoechst is the fixed/reference image; PV is the transformed/moving image.
    tifffile.imwrite(
        out / "aligned_PV_max.tif",
        registered,
        bigtiff=True,
        compression="zlib",
    )
    tifffile.imwrite(
        out / "aligned_Hoechst_max.tif",
        fixed,
        bigtiff=True,
        compression="zlib",
    )

    # ---------------------------------------------------------------
    # QC
    # ---------------------------------------------------------------
    log("Writing QC...")
    save_image_qc(fixed, registered, out)
    save_anchor_visualization(
        fixed,
        moving,
        registered,
        moving_anchor_orig,
        fixed_anchor_orig,
        transform_2d_points_with_castalign(full_final_t, moving_anchor_orig),
        out,
    )
    save_residual_histogram(residual_um, out)

    # ---------------------------------------------------------------
    # Summary
    # ---------------------------------------------------------------
    with open(out / "registration_summary.txt", "w", encoding="utf-8") as f:
        f.write("CASTalign 2-D max-projection registration\n")
        f.write("=" * 68 + "\n\n")
        f.write(f"Fixed image: {args.fixed_image}\n")
        f.write(f"Moving image: {args.moving_image}\n")
        f.write(f"Fixed CellProfiler CSV: {args.fixed_csv}\n")
        f.write(f"Moving CellProfiler CSV: {args.moving_csv}\n\n")
        f.write(f"Image shape (Y,X): {fixed.shape}\n")
        f.write(f"Pixel size X: {args.pixel_size_x:.9f} um\n")
        f.write(f"Pixel size Y: {args.pixel_size_y:.9f} um\n\n")
        f.write(f"Fixed object count: {len(fixed_points)}\n")
        f.write(f"Moving object count: {len(moving_points)}\n")
        f.write(f"Initial anchors: {len(mi0)}\n")
        f.write(f"Final rigid anchors: {len(rigid_residual_px)}\n\n")
        f.write("Phase-correlation initializer\n")
        f.write(f"Shift applied to moving (Y,X px): {shift_px.tolist()}\n")
        f.write(f"Shift applied to moving (Y,X um): {(shift_px * np.array([args.pixel_size_y, args.pixel_size_x])).tolist()}\n")
        f.write(f"Correlation error: {phase_error}\n")
        f.write(f"Phase difference: {phase_diff}\n\n")
        f.write("Rigid CASTalign fit\n")
        f.write(f"Median residual: {np.median(residual_um):.4f} um (approx)\n")
        f.write(f"Mean residual: {np.mean(residual_um):.4f} um (approx)\n")
        f.write(f"95th percentile: {np.percentile(residual_um,95):.4f} um (approx)\n")
        f.write(f"Maximum residual: {np.max(residual_um):.4f} um (approx)\n\n")
        if affine_t is not None:
            f.write("Affine comparison\n")
            f.write(f"Median residual: {np.median(affine_residual_px) * math.sqrt(args.pixel_size_x * args.pixel_size_y):.4f} um (approx)\n")
            f.write(f"Median improvement over rigid: {100 * affine_improvement:.2f}%\n")
            f.write(f"Selected: {chosen_model}\n\n")
        f.write("Scientific note\n")
        f.write("Centroid pairs are registration anchors, not assertions of biological cell identity between channels.\n")
        f.write("The registration is calculated in 2-D max-projection coordinates. The original 3-D stacks are not modified.\n")
        f.write("CASTalign receives 2-D anchors embedded as z=0 in its required (Z,Y,X) volume coordinate system.\n")
        f.write("Rigid is the default model. Affine uses CASTalign LaminarAffine because the max projection is rice-paper/planar geometry.\n")
        f.write("Affine should only be used if the anchor distribution and residual improvement justify it.\n")

    log("")
    log("=" * 72)
    log("DONE")
    log(f"Final model: {chosen_model}")
    log(f"Output: {out}")
    log(f"Final anchors: {len(rigid_residual_px)}")
    log(f"Median residual: {np.median(residual_um):.3f} um (approx)")
    log("=" * 72)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
