#!/usr/bin/env python3
"""Run the complete image-analysis workflow directly from Imaris IMS files.

IMS channel indices are assigned in the user-specified fixed wavelength order:
0=Hoechst, 1=mNeonGreen, 2=BiVe3 virus, 3=PV. Mosaic fields are passed
through the existing BaSiC/BigStitcher workflow; an already-stitched IMS file
is exported and max-projected directly. The resulting channel TIFFs and IMS XY
calibration are then passed to run_stitched_pipeline.py.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import h5py
import numpy as np
import tifffile

import run_stitched_pipeline as stitched_pipeline

try:
    # Merely importing this package registers optional IMS HDF5 filters.
    import hdf5plugin  # noqa: F401
except ImportError:
    hdf5plugin = None


ROOT = Path(__file__).resolve().parent
IMS_MOSAIC_SCRIPT = ROOT / "Stitching/IMSTIFF_trim_zero_padding.py"
FIELD_RE = re.compile(r"^(?P<sample>.+)_F(?P<field>\d{2})\.ims$", re.IGNORECASE)
CHANNEL_ORDER = ("Hoechst", "mNeonGreen", "BiVe3", "PV")
CHANNEL_DISPLAY_NAMES = {
    "Hoechst": "Hoechst",
    "mNeonGreen": "mNeonGreen",
    "BiVe3": "BiVe3 virus",
    "PV": "PV",
}
UNIT_TO_UM = {
    "um": 1.0,
    "µm": 1.0,
    "micrometer": 1.0,
    "micrometers": 1.0,
    "nm": 0.001,
    "mm": 1000.0,
    "m": 1_000_000.0,
}


@dataclass(frozen=True)
class Acquisition:
    name: str
    mode: str
    files: tuple[Path, ...]


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def resolve_path(value: str, base: Path) -> Path:
    path = Path(value).expanduser()
    return (path if path.is_absolute() else base / path).resolve()


def ims_files(directory: Path) -> list[Path]:
    return sorted(
        (
            path.resolve()
            for path in directory.iterdir()
            if path.is_file() and path.suffix.lower() == ".ims"
        ),
        key=lambda path: path.name.lower(),
    )


def discover_acquisitions(input_path: Path, requested_mode: str) -> list[Acquisition]:
    """Discover tiled mosaics or standalone already-stitched IMS acquisitions."""
    if requested_mode not in {"auto", "mosaic", "stitched"}:
        raise ValueError("input_mode must be auto, mosaic, or stitched")
    if not input_path.exists():
        raise FileNotFoundError(f"IMS input does not exist: {input_path}")

    if input_path.is_file():
        if input_path.suffix.lower() != ".ims":
            raise ValueError(f"IMS input file must end in .ims: {input_path}")
        match = FIELD_RE.match(input_path.name)
        if requested_mode == "mosaic" or (requested_mode == "auto" and match):
            if not match:
                raise ValueError(
                    "Mosaic IMS filenames must end in _F00.ims through _F99.ims"
                )
            sample = match.group("sample")
            siblings = []
            for path in ims_files(input_path.parent):
                sibling = FIELD_RE.match(path.name)
                if sibling and sibling.group("sample") == sample:
                    siblings.append((int(sibling.group("field")), path))
            return [
                Acquisition(
                    sample, "mosaic", tuple(path for _, path in sorted(siblings))
                )
            ]
        return [Acquisition(input_path.stem, "stitched", (input_path.resolve(),))]

    candidates = ims_files(input_path)
    if not candidates:
        raise ValueError(f"No .ims files were found directly inside {input_path}")
    tiled = [(path, FIELD_RE.match(path.name)) for path in candidates]

    if requested_mode == "stitched":
        return [Acquisition(path.stem, "stitched", (path,)) for path in candidates]
    if requested_mode == "mosaic" and any(match is None for _, match in tiled):
        bad = [path.name for path, match in tiled if match is None]
        raise ValueError(f"Mosaic mode found filenames without an _FNN suffix: {bad}")
    if requested_mode == "auto":
        has_tiled = any(match is not None for _, match in tiled)
        has_stitched = any(match is None for _, match in tiled)
        if has_tiled and has_stitched:
            raise ValueError(
                "Auto mode will not mix _FNN mosaic fields and standalone IMS files in one folder"
            )
        if has_stitched:
            return [Acquisition(path.stem, "stitched", (path,)) for path in candidates]

    grouped: dict[str, list[tuple[int, Path]]] = {}
    for path, match in tiled:
        assert match is not None
        grouped.setdefault(match.group("sample"), []).append(
            (int(match.group("field")), path)
        )
    return [
        Acquisition(name, "mosaic", tuple(path for _, path in sorted(fields)))
        for name, fields in sorted(grouped.items())
    ]


def numeric_groups(group: h5py.Group) -> list[str]:
    try:
        return sorted(group.keys(), key=lambda item: int(item.rsplit(" ", 1)[-1]))
    except (TypeError, ValueError) as exc:
        raise ValueError(
            f"Unexpected IMS group names below {group.name}: {list(group.keys())}"
        ) from exc


def decode_attribute(value: Any) -> Any:
    """Convert common HDF5 attribute encodings into JSON-safe values."""
    if isinstance(value, np.generic):
        return decode_attribute(value.item())
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace").rstrip("\x00")
    if isinstance(value, np.ndarray):
        if value.dtype.kind == "S":
            return (
                b"".join(value.reshape(-1).tolist())
                .decode("utf-8", errors="replace")
                .rstrip("\x00")
            )
        if value.dtype.kind == "U":
            return "".join(value.reshape(-1).tolist()).rstrip("\x00")
        if value.dtype.kind in {"u", "i"} and value.dtype.itemsize == 1:
            decoded = value.tobytes().decode("utf-8", errors="replace").rstrip("\x00")
            if decoded and all(
                character.isprintable() or character.isspace() for character in decoded
            ):
                return decoded
        return [decode_attribute(item) for item in value.tolist()]
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    return str(value)


def numeric_attribute(group: h5py.Group, name: str) -> float:
    value = decode_attribute(group.attrs[name])
    if isinstance(value, list):
        if len(value) != 1:
            raise ValueError(
                f"IMS attribute {group.name}/{name} is not scalar: {value}"
            )
        value = value[0]
    return float(value)


def voxel_size_from_handle(
    handle: h5py.File, fallback_unit: str
) -> tuple[float, float, float]:
    image = handle["DataSetInfo/Image"]
    unit = str(decode_attribute(image.attrs.get("Unit", ""))).strip().lower()
    if not unit:
        unit = fallback_unit
    if unit not in UNIT_TO_UM:
        raise ValueError(f"Unsupported IMS spatial unit {unit!r} in {handle.filename}")
    factor = UNIT_TO_UM[unit]
    sizes = []
    for axis, axis_name in enumerate(("X", "Y", "Z")):
        extent = numeric_attribute(image, f"ExtMax{axis}") - numeric_attribute(
            image, f"ExtMin{axis}"
        )
        pixels = numeric_attribute(image, axis_name)
        if (
            not math.isfinite(extent)
            or not math.isfinite(pixels)
            or extent <= 0
            or pixels <= 0
        ):
            raise ValueError(f"Invalid {axis_name} extent/size in {handle.filename}")
        sizes.append(extent / pixels * factor)
    return tuple(sizes)  # type: ignore[return-value]


def first_resolution_timepoint(
    handle: h5py.File,
) -> tuple[h5py.Group, str, str, list[str]]:
    dataset = handle["DataSet"]
    resolution = numeric_groups(dataset)[0]
    timepoints = numeric_groups(dataset[resolution])
    if len(timepoints) != 1:
        raise ValueError(
            f"{handle.filename} has {len(timepoints)} timepoints; exactly one is required"
        )
    timepoint = timepoints[0]
    channels = numeric_groups(dataset[resolution][timepoint])
    return dataset[resolution][timepoint], resolution, timepoint, channels


def wavelength_from_attributes(
    attributes: dict[str, Any],
) -> tuple[str | None, float | None]:
    keys = [key for key in attributes if "wavelength" in key.lower()]
    keys.sort(
        key=lambda key: (
            "emission" not in key.lower(),
            "excitation" not in key.lower(),
            key.lower(),
        )
    )
    for key in keys:
        value = attributes[key]
        match = re.search(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)", str(value))
        if match:
            return key, float(match.group())
    return None, None


def inspect_ims(source: Path, fallback_unit: str) -> dict[str, Any]:
    try:
        handle_context = h5py.File(source, "r")
    except OSError as exc:
        raise ValueError(f"Could not open IMS file {source}: {exc}") from exc
    with handle_context as handle:
        channel_root, resolution, timepoint, channel_keys = first_resolution_timepoint(
            handle
        )
        if len(channel_keys) != len(CHANNEL_ORDER):
            raise ValueError(
                f"{source.name} has {len(channel_keys)} channels; this workflow requires exactly "
                f"{len(CHANNEL_ORDER)} in the fixed order {list(CHANNEL_ORDER)}"
            )
        voxel_size = voxel_size_from_handle(handle, fallback_unit)
        channels = []
        for index, (key, canonical) in enumerate(zip(channel_keys, CHANNEL_ORDER)):
            dataset = channel_root[key]["Data"]
            if dataset.ndim != 3:
                raise ValueError(
                    f"{source.name} channel {index} has shape {dataset.shape}; expected ZYX"
                )
            info_path = f"DataSetInfo/Channel {index}"
            attributes = {
                str(name): decode_attribute(value)
                for name, value in (
                    handle[info_path].attrs.items() if info_path in handle else []
                )
            }
            wavelength_attribute, wavelength = wavelength_from_attributes(attributes)
            channels.append(
                {
                    "index": index,
                    "canonical_name": canonical,
                    "display_name": CHANNEL_DISPLAY_NAMES[canonical],
                    "ims_name": str(attributes.get("Name", f"Channel {index}")),
                    "shape_zyx": [int(value) for value in dataset.shape],
                    "dtype": str(dataset.dtype),
                    "wavelength_attribute": wavelength_attribute,
                    "wavelength": wavelength,
                    "relevant_attributes": {
                        name: value
                        for name, value in attributes.items()
                        if name.lower()
                        in {"name", "description", "color", "coloropacity"}
                        or "wavelength" in name.lower()
                    },
                }
            )
        wavelength_values = [channel["wavelength"] for channel in channels]
        order_verified = None
        if all(value is not None for value in wavelength_values):
            numeric_values = [float(value) for value in wavelength_values]
            order_verified = numeric_values == sorted(numeric_values)
        return {
            "path": str(source),
            "resolution_group": resolution,
            "timepoint_group": timepoint,
            "voxel_size_um": {
                "x": voxel_size[0],
                "y": voxel_size[1],
                "z": voxel_size[2],
            },
            "channels": channels,
            "wavelength_order_verified": order_verified,
        }


def inspect_acquisition(acquisition: Acquisition, fallback_unit: str) -> dict[str, Any]:
    files = [inspect_ims(path, fallback_unit) for path in acquisition.files]
    first = files[0]
    first_voxel = np.array(list(first["voxel_size_um"].values()), dtype=float)
    first_shapes = [channel["shape_zyx"] for channel in first["channels"]]
    for info in files[1:]:
        voxel = np.array(list(info["voxel_size_um"].values()), dtype=float)
        if not np.allclose(first_voxel, voxel, rtol=1e-4, atol=1e-9):
            raise ValueError(
                f"{acquisition.name} contains inconsistent voxel sizes: "
                f"{first['voxel_size_um']} versus {info['voxel_size_um']}"
            )
        shapes = [channel["shape_zyx"] for channel in info["channels"]]
        if shapes != first_shapes:
            raise ValueError(
                f"{acquisition.name} contains inconsistent channel stack shapes"
            )
    return {
        "sample": acquisition.name,
        "mode": acquisition.mode,
        "channel_assignment_basis": "zero-based IMS channel index in declared ascending wavelength order",
        "fixed_channel_order": list(CHANNEL_ORDER),
        "voxel_size_um": first["voxel_size_um"],
        "files": files,
    }


def validate_analysis_settings(settings: dict[str, Any]) -> None:
    if not isinstance(settings, dict):
        raise ValueError("analysis must be a JSON object")
    reserved = {"nuclear_channel", "channels", "output_dir", "pixel_size_um"}
    if reserved.intersection(settings):
        raise ValueError(
            "analysis cannot override IMS-derived settings: "
            + ", ".join(sorted(reserved.intersection(settings)))
        )
    prototype = dict(settings)
    prototype.update(
        {
            "nuclear_channel": "Hoechst",
            "channels": {name: f"/planned/{name}.tif" for name in CHANNEL_ORDER},
            "output_dir": "/planned/analysis",
            "pixel_size_um": {"x": 1.0, "y": 1.0},
        }
    )
    with tempfile.TemporaryDirectory(prefix="ims-analysis-config-") as directory:
        path = Path(directory) / "config.json"
        path.write_text(json.dumps(prototype), encoding="utf-8")
        stitched_pipeline.read_config(path)


def read_config(config_path: Path) -> dict[str, Any]:
    raw = json.loads(config_path.read_text(encoding="utf-8"))
    if not isinstance(raw, dict):
        raise ValueError("The top level of the IMS configuration must be a JSON object")
    allowed = {
        "ims_input",
        "input_mode",
        "output_dir",
        "fiji_sif",
        "channel_order",
        "ims_unit",
        "trim_zero_padding",
        "stitching",
        "analysis",
    }
    unknown = set(raw) - allowed
    if unknown:
        raise ValueError(f"Unknown IMS configuration settings: {sorted(unknown)}")
    if raw.get("channel_order", list(CHANNEL_ORDER)) != list(CHANNEL_ORDER):
        raise ValueError(
            f"channel_order must be exactly {list(CHANNEL_ORDER)} for this workflow"
        )
    input_value = raw.get("ims_input")
    output_value = raw.get("output_dir")
    if not isinstance(input_value, str) or not input_value.strip():
        raise ValueError("ims_input must be a nonempty file or directory path")
    if not isinstance(output_value, str) or not output_value.strip():
        raise ValueError("output_dir must be a nonempty path")
    mode = raw.get("input_mode", "auto")
    if mode not in {"auto", "mosaic", "stitched"}:
        raise ValueError("input_mode must be auto, mosaic, or stitched")
    unit = raw.get("ims_unit", "um")
    if unit not in {"um", "nm", "mm", "m"}:
        raise ValueError("ims_unit must be um, nm, mm, or m")
    trim = raw.get("trim_zero_padding", True)
    if not isinstance(trim, bool):
        raise ValueError("trim_zero_padding must be true or false")
    stitching = raw.get("stitching", {})
    if not isinstance(stitching, dict):
        raise ValueError("stitching must be a JSON object")
    stitch_defaults = {
        "grid": [2, 2],
        "overlap": 10.0,
        "min_correlation": 0.55,
        "basic_sample_planes": 12,
        "intensity_downsampling": 32,
        "intensity_max_inliers": 10000,
        "intensity_offset_only": 0.5,
        "intensity_unmodified": 0.5,
    }
    unknown_stitching = set(stitching) - set(stitch_defaults)
    if unknown_stitching:
        raise ValueError(f"Unknown stitching settings: {sorted(unknown_stitching)}")
    stitch_defaults.update(stitching)
    grid = stitch_defaults["grid"]
    if (
        not isinstance(grid, list)
        or len(grid) != 2
        or any(
            isinstance(value, bool) or not isinstance(value, int) or value <= 0
            for value in grid
        )
    ):
        raise ValueError("stitching.grid must contain two positive integers")
    numeric_positive = (
        "basic_sample_planes",
        "intensity_downsampling",
        "intensity_max_inliers",
    )
    for name in numeric_positive:
        value = stitch_defaults[name]
        if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
            raise ValueError(f"stitching.{name} must be a positive integer")
    for name in (
        "overlap",
        "min_correlation",
        "intensity_offset_only",
        "intensity_unmodified",
    ):
        value = stitch_defaults[name]
        if (
            isinstance(value, bool)
            or not isinstance(value, (int, float))
            or not math.isfinite(value)
        ):
            raise ValueError(f"stitching.{name} must be a finite number")
    if not 0 <= float(stitch_defaults["overlap"]) < 100:
        raise ValueError("stitching.overlap must lie in [0, 100)")
    if not 0 <= float(stitch_defaults["min_correlation"]) <= 1:
        raise ValueError("stitching.min_correlation must lie in [0, 1]")
    analysis = raw.get("analysis")
    validate_analysis_settings(analysis)

    fiji_value = raw.get("fiji_sif")
    if fiji_value is not None and (
        not isinstance(fiji_value, str) or not fiji_value.strip()
    ):
        raise ValueError("fiji_sif must be a nonempty path or null")
    result = dict(raw)
    result["_input"] = resolve_path(input_value, config_path.parent)
    result["_output"] = resolve_path(output_value, config_path.parent)
    result["_mode"] = mode
    result["_unit"] = unit
    result["_trim"] = trim
    result["_stitching"] = stitch_defaults
    result["_analysis"] = analysis
    result["_fiji_sif"] = (
        resolve_path(fiji_value, config_path.parent) if fiji_value else None
    )
    if result["_output"].exists() and not result["_output"].is_dir():
        raise ValueError(
            f"output_dir exists but is not a directory: {result['_output']}"
        )
    return result


def source_fingerprint(
    raw_config: dict[str, Any], acquisitions: list[Acquisition]
) -> str:
    files = []
    for acquisition in acquisitions:
        for path in acquisition.files:
            stat = path.stat()
            files.append(
                {"path": str(path), "size": stat.st_size, "mtime_ns": stat.st_mtime_ns}
            )
    code_hash = hashlib.sha256()
    for path in (Path(__file__), IMS_MOSAIC_SCRIPT, ROOT / "run_stitched_pipeline.py"):
        code_hash.update(path.read_bytes())
    payload = {
        "config": {
            key: value for key, value in raw_config.items() if not key.startswith("_")
        },
        "inputs": files,
        "code_sha256": code_hash.hexdigest(),
    }
    return hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()


class IMSRun:
    def __init__(
        self, cfg: dict[str, Any], acquisitions: list[Acquisition], resume: bool
    ):
        self.output: Path = cfg["_output"]
        self.manifest = self.output / "ims_run_manifest.json"
        self.log = self.output / "ims_run.log"
        self.fingerprint = source_fingerprint(cfg, acquisitions)
        self.output.mkdir(parents=True, exist_ok=True)
        if self.manifest.exists():
            previous = json.loads(self.manifest.read_text(encoding="utf-8"))
            if not resume:
                raise RuntimeError(
                    f"Output already contains an IMS run; use --resume or a new output_dir: {self.output}"
                )
            if previous.get("run_fingerprint_sha256") != self.fingerprint:
                raise RuntimeError(
                    "--resume refused: IMS inputs, configuration, or workflow code changed"
                )
            self.started_at = previous["started_at"]
            self.completed = previous.get("completed", {})
        else:
            if any(self.output.iterdir()):
                raise RuntimeError(
                    f"Refusing to mix an IMS run with existing files: {self.output}"
                )
            self.started_at = utc_now()
            self.completed: dict[str, list[str]] = {}
            self.save()

    def save(self) -> None:
        temp = self.manifest.with_suffix(".json.tmp")
        temp.write_text(
            json.dumps(
                {
                    "manifest_version": 1,
                    "started_at": self.started_at,
                    "updated_at": utc_now(),
                    "run_fingerprint_sha256": self.fingerprint,
                    "completed": self.completed,
                },
                indent=2,
                sort_keys=True,
            )
            + "\n",
            encoding="utf-8",
        )
        temp.replace(self.manifest)

    def execute(self, label: str, command: list[str]) -> None:
        print(f"\n{'=' * 70}\nRUNNING {label}\n{'=' * 70}", flush=True)
        print(" ".join(command), flush=True)
        with self.log.open("a", encoding="utf-8") as log:
            log.write(f"\n[{utc_now()}] {label}\n{' '.join(command)}\n")
            log.flush()
            process = subprocess.Popen(
                command,
                cwd=ROOT,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                errors="replace",
                bufsize=1,
            )
            assert process.stdout is not None
            for line in process.stdout:
                print(line, end="", flush=True)
                log.write(line)
            returncode = process.wait()
        if returncode:
            raise subprocess.CalledProcessError(returncode, command)

    def record(self, label: str, paths: list[Path]) -> None:
        missing = [
            str(path) for path in paths if not path.is_file() or not path.stat().st_size
        ]
        if missing:
            raise RuntimeError(f"{label} did not produce required output(s): {missing}")
        self.completed[label] = [str(path) for path in paths]
        self.save()


def usable_z(stack: h5py.Dataset) -> int:
    if stack.ndim != 3:
        raise ValueError(
            f"Expected IMS channel data in ZYX order, got shape {stack.shape}"
        )
    for index in range(stack.shape[0] - 1, -1, -1):
        if np.any(stack[index]):
            return index + 1
    raise ValueError("The first IMS channel contains no non-empty Z planes")


def trailing_zero_padding_xy(stack: h5py.Dataset, z_count: int) -> tuple[int, int]:
    y_total, x_total = stack.shape[1], stack.shape[2]
    slab = np.asarray(
        stack[
            : min(z_count, stack.shape[0]),
            max(0, y_total - 128) :,
            max(0, x_total - 128) :,
        ]
    )
    row_has_data = np.any(slab, axis=(0, 2))
    col_has_data = np.any(slab, axis=(0, 1))
    if not np.any(row_has_data) or not np.any(col_has_data):
        raise ValueError(
            "Could not determine trailing IMS padding from the final 128-pixel XY slab; "
            "set trim_zero_padding to false if this is a real dark image border"
        )
    y_start = max(0, y_total - 128)
    x_start = max(0, x_total - 128)
    return (
        y_start + int(np.flatnonzero(row_has_data)[-1]) + 1,
        x_start + int(np.flatnonzero(col_has_data)[-1]) + 1,
    )


def projection_paths(preprocessing: Path, acquisition: Acquisition) -> dict[str, Path]:
    directory = preprocessing / acquisition.name / "04_max_projections"
    return {
        channel: directory / f"{acquisition.name}_{channel}_stitched_max_projection.tif"
        for channel in CHANNEL_ORDER
    }


def export_stitched_ims(
    acquisition: Acquisition,
    preprocessing: Path,
    trim_zero_padding: bool,
) -> dict[str, Path]:
    source = acquisition.files[0]
    acquisition_dir = preprocessing / acquisition.name
    exported = acquisition_dir / "01_exported"
    projections = acquisition_dir / "04_max_projections"
    exported.mkdir(parents=True, exist_ok=True)
    projections.mkdir(parents=True, exist_ok=True)
    outputs = projection_paths(preprocessing, acquisition)

    with h5py.File(source, "r") as handle:
        channel_root, _, _, channel_keys = first_resolution_timepoint(handle)
        if len(channel_keys) != len(CHANNEL_ORDER):
            raise ValueError(
                f"{source.name} does not contain the required four channels"
            )
        reference = channel_root[channel_keys[0]]["Data"]
        z_count = usable_z(reference)
        source_y, source_x = reference.shape[1:]
        valid_y, valid_x = (
            trailing_zero_padding_xy(reference, z_count)
            if trim_zero_padding
            else (source_y, source_x)
        )
        (exported / "geometry.txt").write_text(
            f"Source IMS array: Y={source_y}, X={source_x}\n"
            f"Exported logical image: Y={valid_y}, X={valid_x}\n"
            f"Trimmed trailing rows: {source_y - valid_y}\n"
            f"Trimmed trailing columns: {source_x - valid_x}\n",
            encoding="utf-8",
        )
        for channel_key, canonical in zip(channel_keys, CHANNEL_ORDER):
            stack_path = exported / f"{canonical}.tif"
            projection_path = outputs[canonical]
            if (
                stack_path.is_file()
                and stack_path.stat().st_size
                and projection_path.is_file()
                and projection_path.stat().st_size
            ):
                continue
            dataset = channel_root[channel_key]["Data"]
            if dataset.shape != reference.shape:
                raise ValueError(
                    f"Channel shape mismatch in {source.name}: {dataset.shape} vs {reference.shape}"
                )
            stack = dataset[:z_count, :valid_y, :valid_x]
            if not stack_path.is_file() or not stack_path.stat().st_size:
                tifffile.imwrite(
                    stack_path, stack, imagej=True, metadata={"axes": "ZYX"}
                )
            if not projection_path.is_file() or not projection_path.stat().st_size:
                tifffile.imwrite(
                    projection_path,
                    np.max(stack, axis=0),
                    imagej=True,
                    metadata={"axes": "YX"},
                    compression="zlib",
                )
    return outputs


def mosaic_command(cfg: dict[str, Any], preprocessing: Path) -> list[str]:
    stitch = cfg["_stitching"]
    command = [
        sys.executable,
        str(IMS_MOSAIC_SCRIPT),
        str(cfg["_input"]),
        "--fiji-sif",
        str(cfg["_fiji_sif"]),
        "--output-root",
        str(preprocessing),
        "--channel-names",
        *CHANNEL_ORDER,
        "--grid",
        *(str(value) for value in stitch["grid"]),
        "--overlap",
        str(stitch["overlap"]),
        "--min-correlation",
        str(stitch["min_correlation"]),
        "--basic-sample-planes",
        str(stitch["basic_sample_planes"]),
        "--ims-unit",
        cfg["_unit"],
        "--intensity-downsampling",
        str(stitch["intensity_downsampling"]),
        "--intensity-max-inliers",
        str(stitch["intensity_max_inliers"]),
        "--intensity-offset-only",
        str(stitch["intensity_offset_only"]),
        "--intensity-unmodified",
        str(stitch["intensity_unmodified"]),
        "--keep-previous",
    ]
    if not cfg["_trim"]:
        command.append("--no-trim-zero-padding")
    return command


def validate_mosaic_runtime(cfg: dict[str, Any]) -> None:
    sif = cfg["_fiji_sif"]
    if sif is None or not sif.is_file():
        raise RuntimeError(
            "Mosaic IMS input requires fiji_sif pointing to the Fiji/BigStitcher SIF"
        )
    if shutil.which("singularity") is None and shutil.which("apptainer") is None:
        raise RuntimeError("Mosaic IMS input requires Singularity or Apptainer in PATH")
    missing = []
    for module in ("basicpy",):
        try:
            __import__(module)
        except ImportError:
            missing.append(module)
    if missing:
        raise RuntimeError(
            "The active ims-mosaic environment is missing: " + ", ".join(missing)
        )


def validate_runtime(cfg: dict[str, Any], needs_mosaic: bool) -> None:
    if hdf5plugin is None:
        raise RuntimeError(
            "hdf5plugin is unavailable; activate the ims-mosaic environment before reading IMS data"
        )
    # This checks the non-CellProfiler analysis dependencies and confirms that
    # the configured cellprofiler-native environment can import CellProfiler.
    stitched_pipeline.validate_runtime(cfg["_analysis"])
    if needs_mosaic:
        validate_mosaic_runtime(cfg)


def write_metadata(
    preprocessing: Path,
    acquisition: Acquisition,
    metadata: dict[str, Any],
    outputs: dict[str, Path],
) -> Path:
    destination = preprocessing / acquisition.name / "ims_metadata.json"
    destination.parent.mkdir(parents=True, exist_ok=True)
    payload = dict(metadata)
    payload["projection_tiffs"] = {name: str(path) for name, path in outputs.items()}
    destination.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return destination


def downstream_config(
    cfg: dict[str, Any],
    acquisition: Acquisition,
    metadata: dict[str, Any],
    outputs: dict[str, Path],
) -> dict[str, Any]:
    result = dict(cfg["_analysis"])
    result.update(
        {
            "nuclear_channel": "Hoechst",
            "channels": {name: str(outputs[name]) for name in CHANNEL_ORDER},
            "output_dir": str(cfg["_output"] / "analysis" / acquisition.name),
            "pixel_size_um": {
                "x": metadata["voxel_size_um"]["x"],
                "y": metadata["voxel_size_um"]["y"],
            },
        }
    )
    return result


def print_plan(
    cfg: dict[str, Any],
    acquisitions: list[Acquisition],
    metadata: dict[str, dict[str, Any]],
) -> None:
    print("IMS preflight: OK")
    print("Fixed channel mapping by ascending-wavelength index:")
    for index, channel in enumerate(CHANNEL_ORDER):
        print(
            f"  {index}: {CHANNEL_DISPLAY_NAMES[channel]} -> downstream name {channel}"
        )
    for acquisition in acquisitions:
        info = metadata[acquisition.name]
        print(
            f"\n{acquisition.name}: mode={acquisition.mode}, files={len(acquisition.files)}"
        )
        print(f"  voxel size (um): {info['voxel_size_um']}")
        first_channels = info["files"][0]["channels"]
        for channel in first_channels:
            wavelength = (
                f", wavelength={channel['wavelength']} ({channel['wavelength_attribute']})"
                if channel["wavelength"] is not None
                else ", wavelength metadata unavailable"
            )
            print(
                f"  [{channel['index']}] IMS {channel['ims_name']!r} -> "
                f"{channel['display_name']}{wavelength}"
            )
        wavelength_check = info["files"][0]["wavelength_order_verified"]
        if wavelength_check is False:
            print(
                "  WARNING: stored wavelength attributes are not ascending; "
                "the declared channel-index convention is still being used"
            )
        elif wavelength_check is True:
            print("  stored wavelength attributes confirm ascending order")
        if acquisition.mode == "mosaic":
            print("  preprocess: IMS export -> BaSiC -> BigStitcher -> max projection")
        else:
            print(
                "  preprocess: already-stitched IMS -> channel TIFFs -> max projection"
            )
        print(
            "  downstream: stitched-TIFF runner from CellProfiler through final analyses"
        )
    print(f"\nOutput root: {cfg['_output']}")


def run(
    cfg: dict[str, Any],
    acquisitions: list[Acquisition],
    metadata: dict[str, dict[str, Any]],
    resume: bool,
) -> None:
    mosaic_acquisitions = [item for item in acquisitions if item.mode == "mosaic"]
    validate_runtime(cfg, bool(mosaic_acquisitions))
    state = IMSRun(cfg, acquisitions, resume)
    preprocessing = state.output / "preprocessing"
    if mosaic_acquisitions:
        state.execute("IMS mosaic preprocessing", mosaic_command(cfg, preprocessing))

    outputs_by_sample: dict[str, dict[str, Path]] = {}
    for acquisition in acquisitions:
        outputs = projection_paths(preprocessing, acquisition)
        if acquisition.mode == "stitched":
            outputs = export_stitched_ims(acquisition, preprocessing, cfg["_trim"])
        metadata_path = write_metadata(
            preprocessing, acquisition, metadata[acquisition.name], outputs
        )
        state.record(
            f"preprocess_{acquisition.name}", [*outputs.values(), metadata_path]
        )
        outputs_by_sample[acquisition.name] = outputs

    config_dir = state.output / "generated_configs"
    config_dir.mkdir(parents=True, exist_ok=True)
    for acquisition in acquisitions:
        generated = downstream_config(
            cfg,
            acquisition,
            metadata[acquisition.name],
            outputs_by_sample[acquisition.name],
        )
        config_path = config_dir / f"{acquisition.name}.stitched_pipeline.json"
        config_path.write_text(json.dumps(generated, indent=2) + "\n", encoding="utf-8")
        parsed = stitched_pipeline.read_config(config_path)
        stitched_pipeline.validate_inputs(parsed)
        command = [
            sys.executable,
            str(ROOT / "run_stitched_pipeline.py"),
            "--config",
            str(config_path),
        ]
        if resume:
            command.append("--resume")
        state.execute(f"downstream analysis {acquisition.name}", command)
        analysis_manifest = (
            cfg["_output"] / "analysis" / acquisition.name / "run_manifest.json"
        )
        state.record(f"analysis_{acquisition.name}", [analysis_manifest])
    print(f"\nFinished complete IMS workflow: {state.output}", flush=True)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        type=Path,
        required=True,
        help="JSON configuration; see ims_pipeline.example.json",
    )
    parser.add_argument(
        "--resume", action="store_true", help="Resume an unchanged prior run"
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Inspect IMS metadata and print the plan without writing outputs",
    )
    args = parser.parse_args(argv)
    config_path = args.config.expanduser().resolve()
    cfg = read_config(config_path)
    acquisitions = discover_acquisitions(cfg["_input"], cfg["_mode"])
    expected_tiles = int(cfg["_stitching"]["grid"][0]) * int(
        cfg["_stitching"]["grid"][1]
    )
    for acquisition in acquisitions:
        if acquisition.mode == "mosaic" and len(acquisition.files) != expected_tiles:
            raise ValueError(
                f"{acquisition.name}: found {len(acquisition.files)} fields, but grid "
                f"{cfg['_stitching']['grid']} requires {expected_tiles}"
            )
        if acquisition.mode == "mosaic":
            field_numbers = [
                int(FIELD_RE.match(path.name).group("field"))  # type: ignore[union-attr]
                for path in acquisition.files
            ]
            expected_fields = list(range(expected_tiles))
            if field_numbers != expected_fields:
                raise ValueError(
                    f"{acquisition.name}: mosaic fields must be F00 through "
                    f"F{expected_tiles - 1:02d}; found {field_numbers}"
                )
    metadata = {
        acquisition.name: inspect_acquisition(acquisition, cfg["_unit"])
        for acquisition in acquisitions
    }
    print_plan(cfg, acquisitions, metadata)
    if args.dry_run:
        if any(item.mode == "mosaic" for item in acquisitions):
            sif = cfg["_fiji_sif"]
            print(
                "Mosaic runtime will additionally require the configured SIF and "
                "Singularity/Apptainer when execution starts."
            )
            print(f"Configured SIF: {sif}")
        return 0
    run(cfg, acquisitions, metadata, args.resume)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
