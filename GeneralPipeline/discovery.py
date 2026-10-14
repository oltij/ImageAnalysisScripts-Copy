"""Input discovery for raw tiled IMS and explicitly mapped stitched TIFFs."""

from __future__ import annotations

import fnmatch
import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .channel_metadata import (
    TILE_TOKEN,
    channel_suffix,
    experimental_set_identifier,
    optional_acquisition_identifier,
    validate_tile_metadata,
)


FIELD_DIR = re.compile(r"^(?:(?P<base>.+?)[_\s-]+)?Field[_\s-]*(?P<number>\d+)$", re.I)
FUSIONSTITCHER_TOKEN = re.compile(r"(?:^|[_\s-])fusionstitcher(?:$|[_\s.-])", re.I)
IGNORE_NAMES = {
    ".git",
    "__pycache__",
    "generated_configs",
    "preprocessing",
    "analysis",
    "output",
    "outputs",
    "logs",
    "input_staging",
    "inputs",
    "results",
}
IGNORED_METADATA_FILES = {"thumbs.db", "desktop.ini", "icon\r"}


@dataclass(frozen=True)
class RawTile:
    source_index: int
    path: Path


@dataclass(frozen=True)
class RawAcquisition:
    path: Path
    name: str
    field_number: int
    tiles: tuple[RawTile, ...]
    session: str | None
    channels: tuple[str, ...]
    nuclear_channel: str
    channel_naming_basis: str
    experimental_set_number: int | None
    acquisition_identifier: int | None
    physical_metadata: dict[str, Any]


@dataclass(frozen=True)
class StitchedAcquisition:
    name: str
    channels: dict[str, Path]
    nuclear_channel: str
    pixel_size_um: dict[str, float]
    tiff_metadata: dict[str, dict[str, Any]]


def fail(message: str) -> None:
    raise ValueError(message)


def _skip_directory(name: str) -> bool:
    folded = name.casefold()
    return (
        name.startswith(".")
        or folded in IGNORE_NAMES
        or folded.endswith("_ims_to_tiff")
        or folded.endswith("_ims_to_tiff_output")
    )


def _is_metadata(path: Path) -> bool:
    return path.name.startswith(".") or path.name.casefold() in IGNORED_METADATA_FILES


def is_derived_stitched_ims(path: Path) -> bool:
    return bool(FUSIONSTITCHER_TOKEN.search(path.stem))


def _tile_index(path: Path) -> int | None:
    match = TILE_TOKEN.search(path.stem)
    return int(match.group(1)) if match else None


def _ims_in_tile_directory(tile_dir: Path) -> Path:
    matches = [
        path
        for path in tile_dir.rglob("*")
        if path.is_file()
        and path.suffix.casefold() == ".ims"
        and not is_derived_stitched_ims(path)
        and not any(
            _is_metadata(Path(part)) for part in path.relative_to(tile_dir).parts
        )
        and len(path.relative_to(tile_dir).parts) <= 3
    ]
    if len(matches) != 1:
        fail(
            f"Expected exactly one original .ims within tile directory {tile_dir}; found {len(matches)}"
        )
    return matches[0].resolve()


def discover_tiles(
    field: Path, expected_count: int, require_contiguous: bool
) -> tuple[RawTile, ...]:
    tiles: dict[int, Path] = {}
    for child in sorted(field.iterdir(), key=lambda path: path.name.casefold()):
        if _is_metadata(child) or is_derived_stitched_ims(child):
            continue
        index = _tile_index(child)
        if child.is_file() and child.suffix.casefold() == ".ims":
            if index is None:
                fail(
                    f"Unnumbered IMS candidate in {field}: {child.name}; expected a final FNN token"
                )
            source = child.resolve()
        elif child.is_dir() and not child.is_symlink() and index is not None:
            source = _ims_in_tile_directory(child)
        else:
            continue
        if index in tiles:
            fail(
                f"Duplicate tile index F{index:02d} in {field}: {tiles[index]} and {source}"
            )
        tiles[index] = source
    if len(tiles) != expected_count:
        fail(
            f"{field}: expected {expected_count} original IMS tiles; found {len(tiles)}: {sorted(tiles)}"
        )
    order = sorted(tiles)
    if require_contiguous and order != list(range(order[0], order[0] + expected_count)):
        fail(f"{field}: noncontiguous IMS tile indices {order}")
    return tuple(RawTile(index, tiles[index]) for index in order)


def _field_paths(
    root: Path, maximum_depth: int
) -> list[tuple[Path, str, int, str | None]]:
    if not root.is_dir():
        fail(f"Input must be an existing directory: {root}")
    result: list[tuple[Path, str, int, str | None]] = []
    for current, directories, _files in os.walk(root, followlinks=False):
        path = Path(current)
        relative = path.relative_to(root)
        match = FIELD_DIR.fullmatch(path.name)
        if match:
            name = match.group("base") or path.parent.name
            # Session is optional metadata: the first layer below the selected
            # root when one exists, with no name-pattern requirement.
            session = relative.parts[0] if len(relative.parts) > 1 else None
            result.append((path.resolve(), name, int(match.group("number")), session))
            directories[:] = []
            continue
        if len(relative.parts) >= maximum_depth:
            directories[:] = []
        else:
            directories[:] = sorted(
                [name for name in directories if not _skip_directory(name)],
                key=str.casefold,
            )
    return result


def _validate_identifiers(
    name: str, tiles: tuple[RawTile, ...], markers: tuple[str, ...]
) -> tuple[int | None, int | None]:
    directory_identifier = optional_acquisition_identifier(
        name, markers, filename=False
    )
    filename_identifiers = [
        optional_acquisition_identifier(tile.path.name, markers, filename=True)
        for tile in tiles
    ]
    if len(set(filename_identifiers)) != 1:
        fail(
            f"Acquisition identifier mismatch across IMS tiles in {name}: {filename_identifiers}"
        )
    if directory_identifier != filename_identifiers[0]:
        fail(
            f"Acquisition identifier mismatch for {name}: directory={directory_identifier!r}, "
            f"IMS filenames={filename_identifiers[0]!r}"
        )
    set_values = [experimental_set_identifier(name, filename=False)] + [
        experimental_set_identifier(tile.path.name, filename=True) for tile in tiles
    ]
    present = {value for value in set_values if value is not None}
    if len(present) > 1:
        fail(f"Experimental Set identifier mismatch in {name}: {set_values}")
    return directory_identifier, next(iter(present)) if present else None


def _biological_channels(
    name: str,
    tiles: tuple[RawTile, ...],
    configured: str | list[str],
    markers: tuple[str, ...],
) -> tuple[list[str], str]:
    if configured == "auto":
        detected = [channel_suffix(tile.path.name, markers) for tile in tiles]
        if any(values != detected[0] for values in detected[1:]):
            fail(
                f"Biological-channel suffix differs across tiles in {name}: {detected}"
            )
        return detected[0], "filename_suffix"
    requested = list(configured)
    try:
        detected = [channel_suffix(tile.path.name, markers) for tile in tiles]
    except ValueError:
        return requested, "explicit"
    if any(values != detected[0] for values in detected[1:]):
        fail(f"Biological-channel suffix differs across tiles in {name}: {detected}")
    if [value.casefold() for value in requested] != [
        value.casefold() for value in detected[0]
    ]:
        fail(
            f"{name}: explicit channel_names {requested} conflicts with validated filename suffix {detected[0]}"
        )
    return requested, "explicit_validated_by_filename"


def discover_raw_acquisitions(
    root: Path,
    settings: dict[str, Any],
    *,
    issues: list[dict[str, Any]] | None = None,
) -> list[RawAcquisition]:
    """Discover complete fields and record only genuinely incomplete fields."""
    grid = settings["grid"]
    if grid != [2, 2]:
        fail("Only the validated 2x2 raw-IMS grid is currently supported")
    expected_count = grid[0] * grid[1]
    session_filter = settings.get("sessions", [])
    group_pattern = settings.get("group_pattern")
    fields = _field_paths(root, settings.get("max_depth", 5))
    if not fields:
        fail(f"No Field_N acquisition directories were found under {root}")
    result: list[RawAcquisition] = []
    selected_count = 0
    for path, name, number, session in fields:
        if session_filter and session not in session_filter:
            continue
        if group_pattern and not fnmatch.fnmatchcase(
            name.casefold(), group_pattern.casefold()
        ):
            continue
        selected_count += 1
        try:
            tiles = discover_tiles(
                path, expected_count, settings.get("require_contiguous_tiles", False)
            )
        except ValueError as exc:
            match = re.search(
                rf"expected {expected_count} original IMS tiles; found (\d+)", str(exc)
            )
            if (
                match
                and int(match.group(1)) < expected_count
                and settings.get("skip_incomplete", True)
            ):
                if issues is not None:
                    issues.append(
                        {
                            "path": str(path),
                            "name": name,
                            "field_number": number,
                            "session": session,
                            "kind": "incomplete",
                            "reason": str(exc),
                        }
                    )
                continue
            raise
        acquisition_identifier, set_number = _validate_identifiers(
            name, tiles, settings["filename_markers"]
        )
        channels, naming_basis = _biological_channels(
            name,
            tiles,
            settings["channel_names"],
            settings["filename_markers"],
        )
        nuclear = settings["nuclear_channel"]
        if nuclear == "auto":
            candidates = [
                value for value in channels if value.casefold() in {"dapi", "hoechst"}
            ]
            if len(candidates) != 1:
                fail(
                    f"{path}: nuclear_channel=auto found {candidates}; configure it explicitly"
                )
            nuclear = candidates[0]
        elif nuclear not in channels:
            fail(
                f"{path}: configured nuclear channel {nuclear!r} is absent from {channels}"
            )
        physical = validate_tile_metadata(
            (tile.path for tile in tiles), channels, settings["ims_unit"]
        )
        result.append(
            RawAcquisition(
                path=path,
                name=name,
                field_number=number,
                tiles=tiles,
                session=session,
                channels=tuple(channels),
                nuclear_channel=nuclear,
                channel_naming_basis=naming_basis,
                experimental_set_number=set_number,
                acquisition_identifier=acquisition_identifier,
                physical_metadata=physical,
            )
        )
    if not selected_count:
        fail(
            "No acquisition directories remain after configured session/group selection"
        )
    if not result:
        fail("No complete raw-IMS acquisitions remain after validation")
    result.sort(
        key=lambda item: (
            str(item.path).casefold(),
            item.name.casefold(),
            item.field_number,
        )
    )
    return result


def discover_stitched_acquisitions(
    root: Path,
    settings: dict[str, Any],
    inspect_tiff: Any,
    *,
    manifest_base: Path | None = None,
) -> list[StitchedAcquisition]:
    """Resolve exact TIFF mappings from entries, a manifest, or a named regex."""
    if not root.is_dir():
        fail(f"Input must be an existing directory: {root}")
    channel_names = settings["channel_names"]
    nuclear = settings["nuclear_channel"]
    calibration = settings["pixel_size_um"]
    mapping = settings["mapping"]
    choices = [
        name
        for name in ("samples", "sample_manifest", "filename_regex")
        if name in mapping
    ]
    if len(choices) != 1:
        fail(
            "stitched_tiffs.mapping must define exactly one of samples, sample_manifest, or filename_regex"
        )
    entries: list[dict[str, Any]]
    if choices[0] == "samples":
        entries = mapping["samples"]
    elif choices[0] == "sample_manifest":
        import json

        value = Path(mapping["sample_manifest"]).expanduser()
        manifest_path = (
            value if value.is_absolute() else (manifest_base or root) / value
        )
        payload = json.loads(manifest_path.read_text(encoding="utf-8"))
        if (
            not isinstance(payload, dict)
            or payload.get("schema") != 1
            or not isinstance(payload.get("samples"), list)
        ):
            fail(
                f"Stitched sample manifest must be a schema-1 object with samples: {manifest_path}"
            )
        entries = payload["samples"]
    else:
        pattern = re.compile(mapping["filename_regex"])
        if not {"sample", "channel"}.issubset(pattern.groupindex):
            fail(
                "filename_regex must contain named groups (?P<sample>...) and (?P<channel>...)"
            )
        grouped: dict[str, dict[str, str]] = {}
        for path in sorted(root.rglob("*"), key=lambda value: str(value).casefold()):
            if not path.is_file() or path.suffix.casefold() not in {".tif", ".tiff"}:
                continue
            match = pattern.fullmatch(path.name)
            if match is None:
                continue
            sample = match.group("sample")
            channel = match.group("channel")
            if channel not in channel_names:
                fail(f"{path}: filename_regex produced unknown channel {channel!r}")
            if channel in grouped.setdefault(sample, {}):
                fail(
                    f"Duplicate stitched TIFF for sample {sample!r}, channel {channel!r}"
                )
            grouped[sample][channel] = str(path)
        entries = [
            {"name": sample, "channels": channels}
            for sample, channels in grouped.items()
        ]
    if not isinstance(entries, list) or not entries:
        fail("No stitched-TIFF samples were configured or discovered")

    result: list[StitchedAcquisition] = []
    seen: set[str] = set()
    for entry in entries:
        if not isinstance(entry, dict) or set(entry) != {"name", "channels"}:
            fail(
                f"Each stitched sample must contain exactly name and channels: {entry!r}"
            )
        name = entry["name"]
        raw_channels = entry["channels"]
        if not isinstance(name, str) or not name.strip() or "/" in name or "\\" in name:
            fail(f"Invalid stitched sample name: {name!r}")
        if name.casefold() in seen:
            fail(f"Duplicate stitched sample name: {name!r}")
        seen.add(name.casefold())
        if not isinstance(raw_channels, dict):
            fail(f"{name}: channels must be an explicit object")
        missing = set(channel_names) - set(raw_channels)
        extra = set(raw_channels) - set(channel_names)
        if missing or extra:
            fail(
                f"{name}: stitched channel mapping missing={sorted(missing)}, extra={sorted(extra)}"
            )
        channels: dict[str, Path] = {}
        metadata: dict[str, dict[str, Any]] = {}
        for channel in channel_names:
            value = raw_channels[channel]
            if not isinstance(value, str) or not value.strip():
                fail(f"{name}/{channel}: TIFF path must be a nonempty string")
            candidate = Path(value).expanduser()
            source = (
                candidate if candidate.is_absolute() else root / candidate
            ).resolve()
            if not source.is_file():
                raise FileNotFoundError(
                    f"{name}/{channel}: stitched TIFF does not exist: {source}"
                )
            channels[channel] = source
            metadata[channel] = inspect_tiff(source)
        if len(set(channels.values())) != len(channels):
            fail(f"{name}: one stitched TIFF was assigned to more than one channel")
        dimensions = {tuple(item["yx_shape"]) for item in metadata.values()}
        if len(dimensions) > 1 and not settings["analysis"]["alignment"].get(
            "enabled", True
        ):
            fail(
                f"{name}: stitched channels have different XY dimensions and alignment is disabled"
            )
        result.append(
            StitchedAcquisition(
                name=name,
                channels=channels,
                nuclear_channel=nuclear,
                pixel_size_um=dict(calibration),
                tiff_metadata=metadata,
            )
        )
    result.sort(key=lambda item: item.name.casefold())
    return result
