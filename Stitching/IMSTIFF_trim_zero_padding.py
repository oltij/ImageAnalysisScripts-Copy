#!/usr/bin/env python3
"""IMS mosaic pipeline: export -> BaSiC -> BigStitcher -> max projection.

Example (four fields, Fiji/BigStitcher inside Singularity):
    python ims_mosaic_pipeline_illumination_bigstitcher.py /data/experiment \
        --fiji-sif /path/to/fiji_latest_bigstitcher.sif

Dependencies: Python 3.10+, h5py, numpy, tifffile.  For illumination
correction install BaSiCPy: ``pip install basicpy``.  Fiji must have
BigStitcher installed and up to date.  The program retains intermediates within the current run for auditability, but
cleans the prior pipeline-generated result directory for the same acquisition
before starting a new run. Use --keep-previous to preserve the old result.
"""
from __future__ import annotations

import argparse
import logging
import re
import shutil
import subprocess
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import h5py
# Registers optional HDF5 compression filters used by some IMS files.
try:
    import hdf5plugin  # noqa: F401
except ImportError:
    hdf5plugin = None
import numpy as np
import tifffile


# Acquisition fields use exactly two digits (F00 through F99).  The anchors
# intentionally reject similarly named processing outputs and F0-style names.
FIELD_RE = re.compile(r"^(?P<sample>.+)_F(?P<field>\d{2})\.ims$", re.IGNORECASE)
INVALID_FILENAME_CHARS = re.compile(r'[\\/:*?"<>|]')


@dataclass(frozen=True)
class Tile:
    path: Path
    field: int


def numeric_groups(group: h5py.Group) -> list[str]:
    """Return Imaris groups in their natural numeric order."""
    return sorted(group.keys(), key=lambda item: int(item.rsplit(" ", 1)[-1]))


def channel_name(handle: h5py.File, index: int) -> str:
    """Return a filesystem-safe Imaris channel name."""
    try:
        value = handle[f"DataSetInfo/Channel {index}"].attrs.get("Name")
        if isinstance(value, np.ndarray):
            if value.dtype.kind == "S":
                value = b"".join(value.tolist()).decode("utf-8", errors="replace")
            elif value.dtype.kind == "U":
                value = "".join(value.tolist())
            else:
                value = str(value)
        elif isinstance(value, bytes):
            value = value.decode("utf-8", errors="replace")
        value = INVALID_FILENAME_CHARS.sub("_", str(value)).strip()
        if value:
            return value
    except (KeyError, TypeError, OSError):
        pass
    return f"Channel{index}"


def usable_z(stack: h5py.Dataset) -> int:
    """Exclude trailing, entirely empty planes occasionally present in IMS data."""
    for index in range(stack.shape[0] - 1, -1, -1):
        if np.any(stack[index]):
            return index + 1
    raise ValueError("the first channel/timepoint has no non-empty z planes")


def numeric_attribute(group: h5py.Group, name: str) -> float:
    """Read an IMS numeric attribute, which may be stored as a C-string."""
    value = group.attrs[name]
    if isinstance(value, np.ndarray):
        value = value.tobytes()
    if isinstance(value, bytes):
        value = value.decode("utf-8", errors="strict").rstrip("\x00")
    return float(value)


def voxel_size_from_ims(source: Path, fallback_unit: str = "um") -> tuple[float, float, float]:
    """Read native physical voxel spacing from IMS extents, in micrometers.

    Imaris records the physical lower/upper bounds (ExtMin*, ExtMax*) and
    image dimensions (X, Y, Z) in DataSetInfo/Image.  Their difference divided
    by the corresponding number of voxels is the authoritative calibration.
    """
    with h5py.File(source, "r") as handle:
        image = handle["DataSetInfo/Image"]
        unit = image.attrs.get("Unit", b"um")
        if isinstance(unit, np.ndarray):
            unit = unit.tobytes()
        if isinstance(unit, bytes):
            unit = unit.decode("utf-8", errors="replace").rstrip("\x00")
        factors = {"um": 1.0, "µm": 1.0, "micrometer": 1.0, "micrometers": 1.0, "nm": 0.001, "mm": 1000.0, "m": 1_000_000.0}
        normalized_unit = str(unit).strip().lower()
        if not normalized_unit:
            # Some acquisition software writes valid extents but leaves this
            # optional IMS attribute empty.  The microscope coordinates used
            # here are conventionally micrometers; retain an explicit escape
            # hatch for datasets known to use another unit.
            normalized_unit = fallback_unit
            logging.warning("%s has no IMS spatial unit; assuming %s", source.name, fallback_unit)
        try:
            factor = factors[normalized_unit]
        except KeyError as error:
            raise ValueError(f"unsupported IMS spatial unit {unit!r}") from error
        result = []
        for axis, axis_name in enumerate(("X", "Y", "Z")):
            extent = numeric_attribute(image, f"ExtMax{axis}") - numeric_attribute(image, f"ExtMin{axis}")
            pixels = numeric_attribute(image, axis_name)
            if extent <= 0 or pixels <= 0:
                raise ValueError(f"invalid {axis_name} extent or size in {source.name}")
            result.append(extent / pixels * factor)
    return tuple(result)  # type: ignore[return-value]


def discover_groups(input_path: Path) -> dict[str, list[Tile]]:
    """Group ``*_FNN.ims`` (two-digit) fields belonging to one acquisition."""
    # Deliberately inspect only the supplied folder itself.  Recursive search
    # could accidentally combine source tiles with files in output or archive
    # subdirectories.
    def directory_ims_files(directory: Path) -> list[Path]:
        return sorted(
            path for path in directory.iterdir()
            if path.is_file() and path.suffix.lower() == ".ims"
        )

    if input_path.is_file():
        selected = FIELD_RE.match(input_path.name)
        if selected:
            # Supplying one field is shorthand for the complete sibling set of
            # that acquisition. This avoids a surprising one-tile grid error.
            sample = selected.group("sample")
            candidates = [
                path
                for path in directory_ims_files(input_path.parent)
                if (match := FIELD_RE.match(path.name)) and match.group("sample") == sample
            ]
        else:
            candidates = [input_path]
    else:
        candidates = directory_ims_files(input_path)
    groups: dict[str, list[Tile]] = defaultdict(list)
    for path in candidates:
        match = FIELD_RE.match(path.name)
        if not match:
            logging.warning("Skipping %s: expected a name ending in _F00.ims through _F99.ims", path.name)
            continue
        groups[match.group("sample")].append(Tile(path.resolve(), int(match.group("field"))))
    return {name: sorted(tiles, key=lambda tile: tile.field) for name, tiles in groups.items()}


def output_root(input_path: Path) -> Path:
    base = input_path.parent if input_path.is_file() else input_path
    return base / f"{base.name}_IMS_to_TIFF"


def trailing_zero_padding_xy(stack: h5py.Dataset, z_count: int) -> tuple[int, int]:
    """Return the logical (Y, X) extent by removing only trailing all-zero rows/columns.

    Some of the microscope IMS files used here store logical 2040 x 1992 data
    in a 2048 x 2048 HDF5 array, with an exactly-zero 8-pixel right pad and
    56-pixel bottom pad.  Those pixels contain no image data and should not be
    exported to TIFF, because doing so creates artificial hard tile boundaries
    for downstream registration/fusion.

    The search is deliberately restricted to the trailing edges.  It never
    removes dark/black regions from the interior of the image.  Only the first
    channel is used to determine the geometry because all channels in an IMS
    acquisition share the same spatial dimensions/padding.
    """
    z_count = min(z_count, stack.shape[0])
    y_total, x_total = stack.shape[1], stack.shape[2]

    # Inspect a modest trailing slab rather than loading the whole volume.
    slab = np.asarray(stack[:z_count, max(0, y_total - 128):y_total,
                             max(0, x_total - 128):x_total])
    row_has_data = np.any(slab, axis=(0, 2))
    col_has_data = np.any(slab, axis=(0, 1))

    if np.any(row_has_data):
        y_start = max(0, y_total - 128)
        valid_y = y_start + int(np.flatnonzero(row_has_data)[-1]) + 1
    else:
        raise ValueError("the reference channel has no non-zero pixels in the trailing XY slab")

    if np.any(col_has_data):
        x_start = max(0, x_total - 128)
        valid_x = x_start + int(np.flatnonzero(col_has_data)[-1]) + 1
    else:
        raise ValueError("the reference channel has no non-zero pixels in the trailing XY slab")

    return valid_y, valid_x


def export_tile_stacks(
    tile: Tile,
    exported: Path,
    overwrite: bool,
    trim_zero_padding: bool = True,
    channel_names: list[str] | None = None,
) -> list[str]:
    """Write logical native-resolution ZYX TIFFs, trimming only trailing zero padding."""
    with h5py.File(tile.path, "r") as handle:
        dataset = handle["DataSet"]
        resolution = numeric_groups(dataset)[0]
        timepoints = numeric_groups(dataset[resolution])
        if len(timepoints) != 1:
            raise ValueError(f"{tile.path.name} has {len(timepoints)} timepoints; this mosaic pipeline requires one")
        channels = numeric_groups(dataset[resolution][timepoints[0]])
        reference = dataset[resolution][timepoints[0]][channels[0]]["Data"]
        z_count = usable_z(reference)

        source_y, source_x = reference.shape[1], reference.shape[2]
        if trim_zero_padding:
            valid_y, valid_x = trailing_zero_padding_xy(reference, z_count)
        else:
            valid_y, valid_x = source_y, source_x
        trim_bottom = source_y - valid_y
        trim_right = source_x - valid_x

        metadata_names = [channel_name(handle, index) for index in range(len(channels))]
        if channel_names is not None and len(channel_names) != len(channels):
            raise ValueError(
                f"{tile.path.name} has {len(channels)} channels, but "
                f"{len(channel_names)} --channel-names were supplied"
            )
        names = list(channel_names) if channel_names is not None else metadata_names
        logging.info(
            "Exporting F%02d: %d channels, %d Z planes; IMS storage=%dx%d YX; "
            "logical image=%dx%d YX; trimming bottom=%d px right=%d px",
            tile.field, len(channels), z_count, source_x, source_y,
            valid_x, valid_y, trim_bottom, trim_right,
        )

        # Record the geometry used for this tile so the trim is auditable.
        geometry_path = exported / f"F{tile.field:02d}_geometry.txt"
        geometry_path.write_text(
            f"Source IMS array: Y={source_y}, X={source_x}\n"
            f"Exported logical image: Y={valid_y}, X={valid_x}\n"
            f"Trimmed trailing rows: {trim_bottom}\n"
            f"Trimmed trailing columns: {trim_right}\n",
            encoding="utf-8",
        )

        for index, channel in enumerate(channels):
            destination = exported / f"F{tile.field:02d}_{names[index]}.tif"
            if destination.exists() and not overwrite:
                continue

            # Read only the logical image extent.  This is the key change: the
            # padded 2048x2048 storage array is never written to TIFF.
            stack = dataset[resolution][timepoints[0]][channel]["Data"][:z_count, :valid_y, :valid_x]

            # These intermediates are deliberately uncompressed: tifffile can
            # memory-map them during BaSiC fitting, keeping RAM use bounded.
            tifffile.imwrite(destination, stack, imagej=True, metadata={"axes": "ZYX"})
    return names


def export_acquisition(
    tiles: list[Tile],
    exported: Path,
    overwrite: bool,
    fallback_unit: str,
    trim_zero_padding: bool = True,
    channel_names: list[str] | None = None,
) -> tuple[list[str], tuple[float, float, float]]:
    exported.mkdir(parents=True, exist_ok=True)
    expected: list[str] | None = None
    calibration: tuple[float, float, float] | None = None
    for tile in tiles:
        tile_calibration = voxel_size_from_ims(tile.path, fallback_unit)
        if calibration is None:
            calibration = tile_calibration
        # Physical extents are often rounded independently for each field.
        # Permit sub-0.01% differences while still rejecting genuinely
        # mismatched acquisition calibrations.
        elif not np.allclose(calibration, tile_calibration, rtol=1e-4, atol=1e-9):
            raise ValueError(f"voxel size in {tile.path.name} differs from the first field: {tile_calibration} versus {calibration}")
        names = export_tile_stacks(
            tile,
            exported,
            overwrite,
            trim_zero_padding=trim_zero_padding,
            channel_names=channel_names,
        )
        if expected is None:
            expected = names
        elif names != expected:
            raise ValueError(f"channel names/order in {tile.path.name} differ from the first field")
    assert expected is not None and calibration is not None
    return expected, calibration


def correct_basicpy(input_files: list[Path], corrected: Path, overwrite: bool, sample_planes: int) -> None:
    """Estimate one BaSiC flat-field per channel, then correct every tile.

    BaSiCPy is the Python implementation of BaSiC.  Sampling evenly spaced Z
    planes from every field makes the estimate feasible without loading an
    entire mosaic into memory; correction is then written one tile at a time.
    """
    try:
        from basicpy import BaSiC
    except ImportError as error:
        raise RuntimeError("BaSiCPy is required: pip install basicpy") from error
    corrected.mkdir(parents=True, exist_ok=True)
    destinations = [corrected / path.name for path in input_files]
    if all(path.exists() for path in destinations) and not overwrite:
        return
    training: list[np.ndarray] = []
    for path in input_files:
        data = tifffile.memmap(path)
        indices = np.unique(np.linspace(0, data.shape[0] - 1, min(sample_planes, data.shape[0]), dtype=int))
        training.extend(np.asarray(data[index], dtype=np.float32) for index in indices)
    model = BaSiC(get_darkfield=False)
    model.fit(np.stack(training))
    for source, destination in zip(input_files, destinations):
        if destination.exists() and not overwrite:
            continue
        source_data = tifffile.memmap(source)
        corrected_stack = model.transform(np.asarray(source_data, dtype=np.float32))
        info = np.iinfo(np.uint16)
        result = np.clip(np.rint(corrected_stack), info.min, info.max).astype(np.uint16)
        tifffile.imwrite(destination, result, imagej=True, metadata={"axes": "ZYX"})


def ij_quote(path: Path) -> str:
    """Return a filesystem path for BigStitcher macro options.

    Do not escape spaces here: in this BigStitcher/ImageJ invocation, backslashes
    become literal characters and cause file-not-found errors. The channel job
    directory itself is made space-free in ``stitch_channel``.
    """
    return str(path).replace("\\", "/")


def bigstitcher_macro(directory: Path, grid_x: int, grid_y: int, overlap: float,
                      min_correlation: float, voxel_size: tuple[float, float, float],
                      channel: str, downsample: int = 4,
                      intensity_downsampling: int = 32,
                      intensity_max_inliers: int = 10000,
                      affine_intensity: bool = True,
                      offset_only: float = 0.5,
                      unmodified: float = 0.5) -> str:
    """Build the proven BigStitcher 3.0.8 macro used in the successful A/B test.

    Apply the known 2x2 acquisition geometry: F00 bottom-left, F01 top-left, F02 top-right, F03 bottom-right, using Snake: Up & Right with the requested overlap. The Automatic Loader receives exactly one FNN.tif per field from a freshly cleaned job directory.
    The legacy SPIM Registration plugin is masked at runtime by stitch_channel. After registration, BigStitcher's overlap-based intensity adjustment is computed from registered overlaps and stored in the dataset XML, then applied during 32-bit float fusion. The macro emits explicit START/DONE markers around the adjustment call so the Python wrapper can verify that the stage was reached.
    """
    base = ij_quote(directory)
    voxel_x, voxel_y, voxel_z = voxel_size

    return f"""setBatchMode(true);

// ============================================================
// 1. DEFINE DATASET
// ============================================================
run(\"Define Multi-View Dataset\",
    \"define_dataset=[Automatic Loader (Bioformats based)] \" +
    \"project_filename=dataset.xml \" +
    \"path={base}/F??.tif \" +
    \"pattern_0=Tiles \" +
    \"modify_voxel_size? voxel_size_x={voxel_x} voxel_size_y={voxel_y} voxel_size_z={voxel_z} voxel_size_unit=um \" +
    \"move_tiles_to_grid_(per_angle)?=[Move Tile to Grid (Macro-scriptable)] \" +
    "grid_type=[Snake: Up & Right] tiles_x={grid_x} tiles_y={grid_y} tiles_z=1 overlap_x_(%)={overlap} overlap_y_(%)={overlap} overlap_z_(%)=0 " +
    \"how_to_load_images=[Re-save as multiresolution HDF5] \" +
    \"dataset_save_path={base}/dataset \" +
    \"check_stack_sizes \" +
    \"subsampling_factors=[{{ {{1,1,1}}, {{2,2,2}}, {{4,4,4}} }}] \" +
    \"hdf5_chunk_sizes=[{{ {{16,16,16}}, {{16,16,16}}, {{16,16,16}} }}] \" +
    \"timepoints_per_partition=1 setups_per_partition=0 use_deflate_compression \" +
    \"export_path={base}/dataset");

if (!File.exists(\"{base}/dataset.xml\"))
    exit(\"BigStitcher dataset creation failed: {base}/dataset.xml\");

// ============================================================
// 2. PHASE CORRELATION (4x in XYZ)
// ============================================================
run(\"Calculate pairwise shifts ...\",
    \"select={base}/dataset.xml process_angle=[All angles] process_channel=[All channels] \" +
    \"process_illumination=[All illuminations] process_tile=[All tiles] process_timepoint=[All Timepoints] \" +
    \"method=[Phase Correlation] channels=[Average Channels] \" +
    \"downsample_in_x={downsample} downsample_in_y={downsample} downsample_in_z={downsample}");

// ============================================================
// 3. FILTER PAIRWISE LINKS
//    Keep links at or above the configured correlation threshold.
// ============================================================
run("Filter pairwise shifts ...",
    "select={base}/dataset.xml filter_by_link_quality min_r={min_correlation} max_r=1");

// 4. GLOBAL OPTIMIZATION
// ============================================================
run(\"Optimize globally and apply shifts ...\",
    \"select={base}/dataset.xml process_angle=[All angles] process_channel=[All channels] \" +
    \"process_illumination=[All illuminations] process_tile=[All tiles] process_timepoint=[All Timepoints] \" +
    \"relative=2.500 absolute=3.500 \" +
    \"global_optimization_strategy=[Two-Round using Metadata to align unconnected Tiles] fix_group_0-0,\");

// ============================================================
// 5. OVERLAP-BASED INTENSITY ADJUSTMENT
//    BigStitcher computes a per-view linear intensity transform from
//    registered overlap pixels. This is intentionally done AFTER
//    BaSiC and AFTER spatial registration, immediately before fusion.
// ============================================================
print(\"BIGSTITCHER_INTENSITY_ADJUSTMENT_START\");
run(\"Adjust Intensities\",
    \"select={base}/dataset.xml \" +
    \"process_angle=[All angles] process_channel=[All channels] \" +
    \"process_illumination=[All illuminations] process_tile=[All tiles] \" +
    \"process_timepoint=[All Timepoints] bounding_box=[All Views] \" +
    \"downsampling={intensity_downsampling} max_inliers={intensity_max_inliers} \" +
    \"affine_intensity offset_only={offset_only} unmodified={unmodified}\");

// ============================================================
// 6. FUSION
//    Apply the computed intensity transforms during 32-bit float fusion.
// ============================================================
print(\"BIGSTITCHER_INTENSITY_ADJUSTMENT_DONE\");
print(\"BIGSTITCHER_FUSION_START\");
run(\"Image Fusion\",
    \"browse={base}/dataset.xml select={base}/dataset.xml \" +
    \"process_angle=[All angles] process_channel=[All channels] \" +
    \"process_illumination=[All illuminations] process_tile=[All tiles] process_timepoint=[All Timepoints] \" +
    \"bounding_box=[All Views] downsampling=1 interpolation=[Linear Interpolation] blending_range=1000 \" +
    \"fusion_type=[Avg, Blending] pixel_type=[32-bit float] adjust_image_intensities \" +
    \"interest_points_for_non_rigid=[-= Disable Non-Rigid =-] preserve_original \" +
    \"produce=[Each timepoint & channel] fused_image=[Save as (compressed) TIFF stacks] \" +
    \"define_input=[Auto-load from input data (values shown below)] output_file_directory={base}/ \" +
    \"filename_addition=Fused\");

eval(\"script\", \"System.exit(0);\");
"""

def make_empty_jar(path: Path) -> None:
    """Create a valid empty JAR used to mask the legacy SPIM plugin."""
    if path.exists():
        return
    import zipfile
    with zipfile.ZipFile(path, "w"):
        pass


def stitch_channel(sif: Path, channel: str, corrected: Path, stitched: Path,
                   grid_x: int, grid_y: int, overlap: float,
                   min_correlation: float, voxel_size: tuple[float, float, float],
                   overwrite: bool,
                   intensity_downsampling: int, intensity_max_inliers: int,
                   affine_intensity: bool, offset_only: float, unmodified: float) -> Path:
    """Run BigStitcher using the original SIF while masking SPIM Registration."""
    # BigStitcher was successfully tested with filesystem-safe channel directories
    # such as ``Confocal_-_Blue``. Do not escape spaces in ImageJ macro option
    # strings: this build interprets the backslashes literally and then cannot
    # find F*.tif. Keep the human-readable channel name for filenames/logs, but
    # use a space-free directory for all BigStitcher paths.
    channel_dir_name = re.sub(r"[^A-Za-z0-9._-]+", "_", channel).strip("_") or "Channel"
    job = stitched / channel_dir_name
    job.mkdir(parents=True, exist_ok=True)
    output = job / f"{channel}_stitched.tif"
    if output.exists() and not overwrite:
        return output

    # Clean up stale BigStitcher input aliases from previous runs.
    # Keep the final fused TIFFs and other non-input files untouched.
    for stale in job.glob("F[0-9][0-9]*.tif"):
        if stale.is_symlink():
            stale.unlink()

    # Present BigStitcher with exactly one input per acquisition field:
    # F00.tif, F01.tif, ...  The underlying files remain the BaSiC-corrected
    # FNN_<channel>.tif files. This prevents old and new aliases from being
    # accidentally included together by a broad glob.
    for source in sorted(corrected.glob(f"F*_{channel}.tif")):
        match = re.match(r"^(F\d{2})_", source.name)
        if not match:
            continue
        target = job / f"{match.group(1)}.tif"
        try:
            target.symlink_to(source)
        except FileExistsError:
            pass
        except OSError:
            shutil.copy2(source, target)

    macro_path = job / "run_bigstitcher.ijm"
    macro_path.write_text(
        bigstitcher_macro(
            job, grid_x, grid_y, overlap, min_correlation, voxel_size,
            channel=channel, intensity_downsampling=intensity_downsampling,
            intensity_max_inliers=intensity_max_inliers,
            affine_intensity=affine_intensity, offset_only=offset_only,
            unmodified=unmodified,
        ),
        encoding="utf-8",
    )

    # Mask only the known-problematic legacy plugin. The successful A/B test
    # showed that the rest of the original SIF works without it.
    empty_jar = stitched / ".empty_SPIM_Registration-5.0.26.jar"
    make_empty_jar(empty_jar)
    masked_target = "/opt/Fiji.app/plugins/SPIM_Registration-5.0.26.jar"

    logging.info(
        "BigStitcher: %s (SPIM_Registration masked; min_r=%.2f; BaSiC retained; "
        "post-registration overlap-based intensity adjustment ON; 32-bit fusion)",
        channel, min_correlation,
    )
    command = [
        "singularity" if shutil.which("singularity") else "apptainer",
        "exec", "--cleanenv",
        "--bind", f"{empty_jar}:{masked_target}",
        str(sif),
        "/opt/Fiji.app/fiji",
        "--headless",
        "-macro", str(macro_path),
    ]
    completed = subprocess.run(command, text=True, capture_output=True)
    (job / "bigstitcher.log").write_text(
        completed.stdout + "\n" + completed.stderr, encoding="utf-8"
    )
    if completed.returncode:
        raise RuntimeError(f"BigStitcher failed for {channel}; see {job / 'bigstitcher.log'}")

    log_text = completed.stdout + "\n" + completed.stderr
    if "BIGSTITCHER_INTENSITY_ADJUSTMENT_START" not in log_text:
        raise RuntimeError(f"BigStitcher intensity-adjustment stage was not reached for {channel}; see {job / 'bigstitcher.log'}")
    if "BIGSTITCHER_INTENSITY_ADJUSTMENT_DONE" not in log_text:
        raise RuntimeError(f"BigStitcher intensity-adjustment stage did not complete for {channel}; see {job / 'bigstitcher.log'}")
    if "BIGSTITCHER_FUSION_START" not in log_text:
        raise RuntimeError(f"BigStitcher fusion stage was not reached for {channel}; see {job / 'bigstitcher.log'}")

    fused_candidates = sorted(job.glob("Fused_fused_tp_*_ch_*.tif"))
    if not fused_candidates:
        fused_candidates = sorted(job.glob("*Fused*.tif"))
    if not fused_candidates:
        raise RuntimeError(
            f"BigStitcher completed but produced no fused TIFF for {channel}; "
            f"see {job / 'bigstitcher.log'}"
        )

    shutil.copy2(fused_candidates[0], output)
    return output

def max_project(source: Path, destination: Path, overwrite: bool) -> None:
    if destination.exists() and not overwrite:
        return
    # BigStitcher may produce compressed 32-bit float TIFF stacks here; read
    # the final fused volume directly so the intensity-adjusted values are preserved.
    image = tifffile.imread(source)
    projection = np.max(image, axis=0)
    tifffile.imwrite(destination, projection, imagej=True, metadata={"axes": "YX"}, compression="zlib")


def clean_previous_results(root: Path, sample_name: str) -> None:
    """Remove all outputs from a previous run for this acquisition only.

    The IMS source files live outside ``root``.  Only the pipeline-generated
    sample directory is removed, including exported TIFFs, BaSiC-corrected
    TIFFs, BigStitcher XML/Zarr state, symlinks, logs, and projections.
    """
    acquisition = root / sample_name
    if acquisition.exists() or acquisition.is_symlink():
        logging.warning("Cleaning previous pipeline results: %s", acquisition)
        if acquisition.is_symlink() or acquisition.is_file():
            acquisition.unlink()
        else:
            shutil.rmtree(acquisition)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Convert IMS mosaic fields to BaSiC-corrected, BigStitcher-registered, intensity-adjusted, fused max projections.")
    parser.add_argument("path", type=Path, help="an IMS file or a directory containing IMS files")
    parser.add_argument("--fiji-sif", type=Path, required=True, help="Singularity/Apptainer SIF containing Fiji + BigStitcher")
    parser.add_argument("--output-root", type=Path, help="optional output root; default is <input folder>_IMS_to_TIFF beside the inputs")
    parser.add_argument("--channel-names", nargs="+", metavar="NAME", help="optional channel names assigned by zero-based IMS channel order")
    parser.add_argument("--grid", nargs=2, type=int, metavar=("X", "Y"), default=(2, 2), help="regular grid dimensions (default: 2 2)")
    parser.add_argument("--overlap", type=float, default=10.0, help="regular-grid overlap percentage in X/Y (default: 10)")
    parser.add_argument("--min-correlation", type=float, default=0.55, help="Minimum pairwise Phase-Correlation r retained for optimization (default: 0.55)")
    parser.add_argument("--basic-sample-planes", type=int, default=12, help="Z planes per tile used to fit BaSiC (default: 12)")
    parser.add_argument("--voxel-size", nargs=3, type=float, metavar=("X_UM", "Y_UM", "Z_UM"), help="optional manual calibration override in micrometers; normally read automatically from IMS")
    parser.add_argument("--ims-unit", choices=("um", "nm", "mm", "m"), default="um", help="unit to assume only when an IMS file leaves its Unit metadata blank (default: um)")
    parser.add_argument("--intensity-downsampling", type=int, default=32, help="downsampling used by BigStitcher overlap-based intensity adjustment (default: 32)")
    parser.add_argument("--intensity-max-inliers", type=int, default=10000, help="maximum overlap pixels per image pair for intensity adjustment (default: 10000)")
    parser.add_argument("--intensity-offset-only", type=float, default=0.5, help="BigStitcher offset-only regularization weight lambda1 (default: 0.5)")
    parser.add_argument("--intensity-unmodified", type=float, default=0.5, help="BigStitcher unmodified-intensity regularization weight lambda2 (default: 0.5)")
    parser.add_argument("--overwrite", action="store_true", help="accepted for compatibility; outputs are cleaned by default on each run")
    parser.add_argument("--keep-previous", action="store_true", help="do not remove the previous pipeline result directory before rerunning")
    parser.add_argument("--keep-input-stacks", action="store_true", help="document intent; inputs are retained by default")
    parser.add_argument("--no-trim-zero-padding", action="store_true", help="keep padded all-zero TIFF borders (not recommended for these IMS files)")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    input_path = args.path.expanduser().resolve()
    if not input_path.exists():
        raise SystemExit(f"Input does not exist: {input_path}")
    if not args.fiji_sif.is_file():
        raise SystemExit(f"Fiji SIF does not exist: {args.fiji_sif}")
    if hdf5plugin is None:
        raise SystemExit(
            "hdf5plugin is unavailable; activate the ims-mosaic environment before reading IMS files"
        )
    if args.channel_names:
        if len(set(args.channel_names)) != len(args.channel_names):
            raise SystemExit("--channel-names must be unique")
        invalid = [name for name in args.channel_names if not re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*", name)]
        if invalid:
            raise SystemExit(
                "--channel-names must start with a letter and contain only letters, digits, "
                f"and underscores; invalid: {invalid}"
            )
    if shutil.which("singularity") is None and shutil.which("apptainer") is None:
        raise SystemExit("Neither singularity nor apptainer is available; load the Singularity module first.")
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    groups = discover_groups(input_path)
    if not groups:
        raise SystemExit("No IMS files with names ending _FNN.ims were found.")
    root = args.output_root.expanduser().resolve() if args.output_root else output_root(input_path)
    root.mkdir(parents=True, exist_ok=True)
    for name, tiles in groups.items():
        if len(tiles) != args.grid[0] * args.grid[1]:
            raise ValueError(f"{name}: found {len(tiles)} fields, but --grid {args.grid[0]} {args.grid[1]} requires {args.grid[0] * args.grid[1]}")
        logging.info("Processing %s (%d fields)", name, len(tiles))
        if not args.keep_previous:
            clean_previous_results(root, name)
        acquisition = root / name
        exported, corrected, stitched, projections = (acquisition / item for item in ("01_exported", "02_basic_corrected", "03_stitched", "04_max_projections"))
        channels, ims_voxel_size = export_acquisition(
            tiles, exported, args.overwrite, args.ims_unit,
            trim_zero_padding=not args.no_trim_zero_padding,
            channel_names=args.channel_names,
        )
        voxel_size = tuple(args.voxel_size) if args.voxel_size else ims_voxel_size
        logging.info("Using voxel size %.6g x %.6g x %.6g um (%s)", *voxel_size, "manual override" if args.voxel_size else "read from IMS")
        projections.mkdir(parents=True, exist_ok=True)
        for channel in channels:
            inputs = sorted(exported.glob(f"F*_{channel}.tif"))
            correct_basicpy(inputs, corrected, args.overwrite, args.basic_sample_planes)
            fused = stitch_channel(
                args.fiji_sif, channel, corrected, stitched,
                *args.grid, args.overlap, args.min_correlation, voxel_size,
                args.overwrite, args.intensity_downsampling,
                args.intensity_max_inliers, True,
                args.intensity_offset_only, args.intensity_unmodified,
            )
            max_project(fused, projections / f"{name}_{channel}_stitched_max_projection.tif", args.overwrite)
    logging.info("Finished. Results: %s", root)
    return 0


if __name__ == "__main__":
    sys.exit(main())
