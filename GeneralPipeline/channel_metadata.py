"""Safe biological-channel parsing and IMS metadata validation.

Biological names come from an explicit mapping or a validated filename suffix.
Generic Imaris labels such as ``Confocal - Green`` are retained as provenance,
but are never interpreted as a biological identity.
"""

from __future__ import annotations

import math
import re
from pathlib import Path
from typing import Any, Iterable


CHANNEL_NAME = re.compile(r"^[A-Za-z][A-Za-z0-9_]*$")
TILE_TOKEN = re.compile(r"(?:^|[_\s-])F[_\s-]?(\d{1,3})$", re.I)
SET_TOKEN = re.compile(r"^Set(?P<number>\d+)$", re.I)
UNIT_TO_UM = {"um": 1.0, "µm": 1.0, "nm": 0.001, "mm": 1000.0, "m": 1_000_000.0}

# This vocabulary is only a parser dictionary. It does not define a required
# panel or map a physical channel index to a biological marker.
DEFAULT_FILENAME_MARKERS = (
    "DAPI",
    "Hoechst",
    "GFP",
    "mNeonGreen",
    "mCherry",
    "RFP",
    "YFP",
    "dTomato",
    "tdTomato",
    "BiVe3",
    "BiVe4",
    "PV",
    "PVALB",
    "LHX6",
    "ASCL1",
    "SOX2",
    "SOX1",
    "NKX2_1",
    "OLIG2",
    "OLIG1",
    "DLX1",
    "DLX2",
    "DLX5",
    "GAD1",
    "GAD2",
    "GAD65",
    "GAD67",
    "SST",
    "VIP",
    "CALB1",
    "CALB2",
    "NEUN",
    "MAP2",
    "GFAP",
    "PCNA",
    "EdU",
    "KI67",
    "MKI67",
    "DCX",
    "FOXG1",
    "TBR1",
    "TUJ1",
    "SATB2",
    "ISL1",
    "IBA1",
    "SLC17A7",
    "SLC17A6",
    "GABA",
    "CHGA",
    "CHGB",
)


def fail(message: str) -> None:
    raise ValueError(message)


def validate_channel_names(values: Iterable[str]) -> list[str]:
    names = list(values)
    if len(names) < 2 or len(names) > 12:
        fail("channel_names must contain 2-12 biological channel names")
    if any(
        not isinstance(name, str) or not CHANNEL_NAME.fullmatch(name) for name in names
    ):
        fail(
            "Channel names must use ASCII letters, digits, and underscores and start with a letter"
        )
    if len({name.casefold() for name in names}) != len(names):
        fail("Channel names must be unique ignoring case")
    return names


def marker_vocabulary(extra: Iterable[str] = ()) -> tuple[str, ...]:
    result = list(DEFAULT_FILENAME_MARKERS)
    seen = {name.casefold() for name in result}
    for marker in extra:
        if not isinstance(marker, str) or not CHANNEL_NAME.fullmatch(marker):
            fail(f"Invalid filename marker {marker!r}")
        if marker.casefold() not in seen:
            result.append(marker)
            seen.add(marker.casefold())
    return tuple(result)


def _filename_prefix(filename: str) -> tuple[list[str], bool]:
    stem = Path(filename).stem
    ending = TILE_TOKEN.search(stem)
    before = stem[: ending.start()] if ending else stem
    tokens = [part for part in before.split("_") if part]
    return tokens, ending is not None


def channel_suffix(filename: str, markers: Iterable[str]) -> list[str]:
    """Return the longest validated marker suffix before optional ID and FNN."""
    marker_map: dict[tuple[str, ...], str] = {}
    for marker in markers:
        if not isinstance(marker, str) or not CHANNEL_NAME.fullmatch(marker):
            fail(f"Invalid filename marker {marker!r}")
        key = tuple(part.casefold() for part in marker.split("_"))
        previous = marker_map.get(key)
        if previous is not None and previous != marker:
            fail(f"Duplicate filename marker differing only by case: {marker!r}")
        marker_map[key] = marker

    tokens, has_tile = _filename_prefix(filename)
    if not tokens:
        fail(f"Cannot find channel names in {filename!r}")
    numeric_is_marker = any(
        len(key) >= 2
        and len(tokens) >= len(key)
        and tuple(part.casefold() for part in tokens[-len(key) :]) == key
        for key in marker_map
    )
    if (
        has_tile
        and len(tokens) >= 2
        and tokens[-2].casefold() == "field"
        and tokens[-1].isdecimal()
    ):
        tokens = tokens[:-2]
    elif has_tile and tokens[-1].isdecimal() and not numeric_is_marker:
        tokens.pop()

    result: list[str] = []
    maximum = max((len(key) for key in marker_map), default=1)
    while tokens:
        match: tuple[int, str] | None = None
        for size in range(min(maximum, len(tokens)), 0, -1):
            key = tuple(part.casefold() for part in tokens[-size:])
            if key in marker_map:
                match = (size, marker_map[key])
                break
        if match is None:
            break
        size, marker = match
        result.append(marker)
        del tokens[-size:]
    if not result:
        fail(
            f"Cannot identify a biological-channel suffix in {filename!r}; "
            "configure channel_names explicitly or extend filename_markers"
        )
    return result[::-1]


def optional_acquisition_identifier(
    value: str, markers: Iterable[str], *, filename: bool
) -> int | None:
    """Find a standalone integer after the final marker, never marker digits."""
    stem = Path(value).stem if filename else value
    if filename:
        ending = TILE_TOKEN.search(stem)
        if ending is None:
            fail(f"Cannot find final FNN tile token in IMS filename {value!r}")
        stem = stem[: ending.start()]
    tokens = [part for part in stem.split("_") if part]
    if not tokens or not tokens[-1].isdecimal():
        return None
    vocabulary = [marker.split("_") for marker in markers]
    numeric_is_marker = any(
        len(parts) >= 2
        and len(tokens) >= len(parts)
        and [item.casefold() for item in tokens[-len(parts) :]]
        == [item.casefold() for item in parts]
        for parts in vocabulary
    )
    return None if numeric_is_marker else int(tokens[-1])


def experimental_set_identifier(value: str, *, filename: bool) -> int | None:
    """Extract SetN as experiment metadata, never as a biological channel."""
    stem = Path(value).stem if filename else value
    if filename:
        ending = TILE_TOKEN.search(stem)
        if ending is None:
            fail(f"Cannot find final FNN tile token in IMS filename {value!r}")
        stem = stem[: ending.start()]
    matches = [
        int(match.group("number"))
        for token in re.split(r"[_\s-]+", stem)
        if (match := SET_TOKEN.fullmatch(token))
    ]
    if len(set(matches)) > 1:
        fail(f"Conflicting Set identifiers in {value!r}: {matches}")
    return matches[0] if matches else None


def decode_attribute(value: Any) -> Any:
    if hasattr(value, "item") and not hasattr(value, "shape"):
        return decode_attribute(value.item())
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace").rstrip("\x00")
    dtype = getattr(value, "dtype", None)
    if dtype is not None and hasattr(value, "reshape"):
        flattened = value.reshape(-1)
        if dtype.kind == "S":
            return (
                b"".join(flattened.tolist())
                .decode("utf-8", errors="replace")
                .rstrip("\x00")
            )
        if dtype.kind == "U":
            return "".join(flattened.tolist()).rstrip("\x00")
        if dtype.kind in {"u", "i"} and dtype.itemsize == 1:
            decoded = value.tobytes().decode("utf-8", errors="replace").rstrip("\x00")
            if decoded and all(
                character.isprintable() or character.isspace() for character in decoded
            ):
                return decoded
        value = value.tolist()
    if isinstance(value, (list, tuple)):
        if len(value) == 1:
            return decode_attribute(value[0])
        return [decode_attribute(item) for item in value]
    return (
        value
        if isinstance(value, (str, int, float, bool)) or value is None
        else str(value)
    )


def _numeric_groups(group: Any) -> list[str]:
    try:
        return sorted(group.keys(), key=lambda name: int(name.rsplit(" ", 1)[-1]))
    except (TypeError, ValueError) as exc:
        fail(f"Unexpected numeric IMS groups below {group.name}: {list(group.keys())}")
        raise AssertionError from exc


def _number(group: Any, name: str) -> float:
    value = decode_attribute(group.attrs[name])
    if isinstance(value, list):
        if len(value) != 1:
            fail(f"IMS attribute {group.name}/{name} is not scalar: {value}")
        value = value[0]
    return float(value)


def _voxel_size_um(handle: Any, fallback_unit: str) -> dict[str, float]:
    image = handle["DataSetInfo/Image"]
    unit = (
        str(decode_attribute(image.attrs.get("Unit", ""))).strip().lower()
        or fallback_unit
    )
    if unit not in UNIT_TO_UM:
        fail(f"Unsupported IMS spatial unit {unit!r} in {handle.filename}")
    factor = UNIT_TO_UM[unit]
    result: dict[str, float] = {}
    for axis, label in enumerate(("X", "Y", "Z")):
        extent = _number(image, f"ExtMax{axis}") - _number(image, f"ExtMin{axis}")
        pixels = _number(image, label)
        if (
            not math.isfinite(extent)
            or not math.isfinite(pixels)
            or extent <= 0
            or pixels <= 0
        ):
            fail(f"Invalid {label} calibration in {handle.filename}")
        result[label.lower()] = extent / pixels * factor
    return result


def _wavelength(attributes: dict[str, Any]) -> float | None:
    keys = sorted(
        (key for key in attributes if "wavelength" in key.casefold()),
        key=lambda key: (
            "emission" not in key.casefold(),
            "excitation" not in key.casefold(),
            key.casefold(),
        ),
    )
    for key in keys:
        match = re.search(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)", str(attributes[key]))
        if match:
            return float(match.group())
    return None


def read_ims_metadata(source: Path, fallback_unit: str = "um") -> dict[str, Any]:
    """Read channel, calibration, and canvas metadata without loading pixels."""
    try:
        import h5py
    except ImportError as exc:  # pragma: no cover - environment-specific guidance
        raise RuntimeError("IMS metadata inspection requires h5py") from exc
    try:
        handle_context = h5py.File(source, "r")
    except OSError as exc:
        raise ValueError(f"Could not open IMS file {source}: {exc}") from exc
    try:
        with handle_context as handle:
            resolutions = _numeric_groups(handle["DataSet"])
            resolution = resolutions[0]
            timepoints = _numeric_groups(handle["DataSet"][resolution])
            if len(timepoints) != 1:
                fail(
                    f"{source}: exactly one IMS time point is required; found {len(timepoints)}"
                )
            channel_root = handle["DataSet"][resolution][timepoints[0]]
            channel_keys = _numeric_groups(channel_root)
            indices = [int(name.rsplit(" ", 1)[-1]) for name in channel_keys]
            names: list[str] = []
            wavelengths: list[float | None] = []
            shapes: list[list[int]] = []
            for index, key in zip(indices, channel_keys):
                group = handle.get(f"DataSetInfo/Channel {index}")
                attributes = (
                    {
                        str(name): decode_attribute(value)
                        for name, value in group.attrs.items()
                    }
                    if group is not None
                    else {}
                )
                names.append(str(attributes.get("Name", "")).strip())
                wavelengths.append(_wavelength(attributes))
                data = channel_root[key].get("Data")
                if data is None or data.ndim != 3:
                    fail(
                        f"{source}: physical channel {index} must contain one ZYX Data stack"
                    )
                shapes.append([int(value) for value in data.shape])
            return {
                "count": len(channel_keys),
                "indices": indices,
                "stored_names": names,
                "wavelengths": wavelengths,
                "voxel_size_um": _voxel_size_um(handle, fallback_unit),
                "zyx_shapes": shapes,
            }
    except (KeyError, OSError, TypeError, StopIteration) as exc:
        raise ValueError(
            f"Missing or malformed IMS metadata in {source}: {exc}"
        ) from exc


def _same_calibration(left: dict[str, float], right: dict[str, float]) -> bool:
    return all(
        math.isclose(left[key], right[key], rel_tol=1e-4, abs_tol=1e-9)
        for key in ("x", "y", "z")
    )


def validate_tile_metadata(
    sources: Iterable[Path], biological_channels: list[str], fallback_unit: str = "um"
) -> dict[str, Any]:
    """Require a consistent, one-to-one physical/biological mapping on all tiles."""
    paths = list(sources)
    if not paths:
        fail("An acquisition has no IMS tiles")
    metadata = [read_ims_metadata(path, fallback_unit) for path in paths]
    reference = metadata[0]
    if reference["indices"] != list(range(reference["count"])):
        fail(
            f"{paths[0]}: IMS physical channel indices are not sequential from zero: {reference['indices']}"
        )
    phase = [name for name in reference["stored_names"] if "phase" in name.casefold()]
    if phase:
        fail(
            f"{paths[0]}: unsupported physical phase channel(s) {phase}; "
            "source-index exclusion is intentionally not implemented"
        )
    if reference["count"] != len(biological_channels):
        fail(
            f"{paths[0]}: IMS contains {reference['count']} physical channels but "
            f"{len(biological_channels)} biological channels are configured ({biological_channels}); "
            "phase or other extra channels are not automatically excluded"
        )
    for path, item in zip(paths[1:], metadata[1:]):
        if item["indices"] != reference["indices"]:
            fail(f"{path}: physical channel indices differ across IMS tiles")
        if item["stored_names"] != reference["stored_names"]:
            fail(f"{path}: stored physical channel names differ across IMS tiles")
        if item["wavelengths"] != reference["wavelengths"]:
            fail(f"{path}: wavelength metadata differ across IMS tiles")
        if item["zyx_shapes"] != reference["zyx_shapes"]:
            fail(f"{path}: channel dimensions differ across IMS tiles")
        if not _same_calibration(item["voxel_size_um"], reference["voxel_size_um"]):
            fail(f"{path}: voxel calibration differs across IMS tiles")
    return reference
