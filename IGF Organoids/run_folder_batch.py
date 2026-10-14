#!/usr/bin/env python3
"""Discover, stage and run 2x2 IMS field mosaics without modifying scientific scripts.

Pure-stdlib commands: scan, prepare, inspect.
Scientific modules from the ImageAnalysisScripts-Copy repo are imported only by
run/preflight, on the compute node with the ims-mosaic environment activated.

Inputs may be an Olti_* session, a group containing Field_* directories, or a
single Field_* directory. An optional date-folder layer is supported. A field
is one biological acquisition assembled from four IMS tiles; it is NOT an IMS
file named ..._F01 (the latter is one tile).
"""
from __future__ import annotations

import argparse
import csv
import fnmatch
import hashlib
import io
import json
import os
import re
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

FIELD_DIR = re.compile(r"^(?:(?P<base>.+?)[_\s-]+)?Field[_\s-]*(?P<num>\d+)$", re.I)
TILE_TOKEN = re.compile(r"(?:^|[_\s-])F[_\s-]?(\d{1,3})$", re.I)
IGNORE_NAMES = {".git", "__pycache__", "generated_configs", "preprocessing", "analysis", "output", "outputs", "logs", "input_staging", "inputs", "results", "00_max_projections", "01_exported", "02_basic_corrected", "03_stitched", "04_max_projections", "05_analysis_projections"}
CHANNEL = re.compile(r"^[A-Za-z][A-Za-z0-9_]*$")
VALID_PROFILES = {"generic", "reporter"}
SCRIPT_PATH = Path(__file__).resolve()
SCRIPT_DIR = SCRIPT_PATH.parent
REPOSITORY_ROOT = SCRIPT_DIR.parent
SBATCH_PATH = SCRIPT_DIR / "run_folder_batch.sbatch"
TURBO_IMS_ROOT = Path("/nfs/turbo/umms-parent/andor_microscope_data_dump")
SCRATCH_ROOT = Path("/scratch")
# Built-in marker vocabulary used only to find the channel-name suffix in source
# filenames. Biological identity comes from the filename, never a fixed index
# or the Imaris channel Name attribute. Add uncommon markers through the spec.
DEFAULT_FILENAME_MARKERS = (
    'DAPI', 'Hoechst', 'GFP', 'mNeonGreen', 'mCherry', 'RFP', 'YFP',
    'dTomato', 'tdTomato', 'BiVe3', 'BiVe4', 'PV', 'PVALB',
    'LHX6', 'ASCL1', 'SOX2', 'SOX1', 'NKX2_1',
    'OLIG2', 'OLIG1', 'DLX1', 'DLX2', 'DLX5', 'GAD1', 'GAD2',
    'GAD65', 'GAD67', 'SST', 'VIP', 'CALB1', 'CALB2', 'NEUN',
    'MAP2', 'GFAP', 'PCNA', 'EdU', 'KI67', 'MKI67', 'DCX',
    'FOXG1', 'TBR1', 'TUJ1', 'SATB2', 'ISL1', 'IBA1',
    'SLC17A7', 'SLC17A6', 'GABA', 'CHGA', 'CHGB',
)
NUCLEAR_MARKERS = {'dapi', 'hoechst'}
# Files left by macOS Finder and Windows Explorer are not IMS inputs.
# Keep this separate from IMS filename validation: a genuinely malformed
# *.ims acquisition must still fail loudly rather than disappearing.
IGNORED_METADATA_FILES = {'thumbs.db', 'desktop.ini', 'icon\r'}
# FusionStitcher-created .ims exports can sit next to the original four
# acquisition tiles, sometimes with a misleading _F0.ims suffix. They are
# derived, already-stitched products and must NEVER be used as raw 2x2 tiles.
# Require an actual filename token, not an incidental substring in a marker.
FUSIONSTITCHER_TOKEN = re.compile(r"(?:^|[_\s-])fusionstitcher(?:$|[_\s.-])", re.I)


def is_derived_stitched_ims(path: Path) -> bool:
    """Recognize FusionStitcher exports without inspecting image contents."""
    return bool(FUSIONSTITCHER_TOKEN.search(path.stem))



def is_os_metadata(path: Path) -> bool:
    return path.name.startswith('.') or path.name.casefold() in IGNORED_METADATA_FILES


def channel_suffix(filename: str, markers: list[str] | tuple[str, ...]) -> list[str]:
    """Read marker suffix from filename, directly before optional ID and FNN.

    The anchor is the FINAL _FNN tile token, not an arbitrary number. Only a
    STANDALONE all-digit identifier directly before FNN may be stripped. Thus
    ASCL1, LHX6, BiVe3, NKX2_1, etc. are retained. Names can be multi-token
    (e.g. NKX2_1) if declared in the marker vocabulary. Unknown marker tokens
    cause an error instead of a silent shift of all biological channel labels.
    """
    stem = Path(filename).stem
    ending = TILE_TOKEN.search(stem)
    before = stem[:ending.start()] if ending else stem
    tokens = [part for part in before.split('_') if part]
    if not tokens:
        fail(f'Cannot find channels in IMS filename {filename!r}')
    # Index by token tuple; longest-match supports multi-token marker names.
    marker_map: dict[tuple[str, ...], str] = {}
    for marker in markers:
        if not isinstance(marker, str) or not CHANNEL.fullmatch(marker):
            fail(f'Invalid filename marker {marker!r}: use ASCII letters/digits/underscores')
        key = tuple(t.casefold() for t in marker.split('_'))
        if key in marker_map and marker_map[key] != marker:
            fail(f'Duplicate filename marker differing only in case: {marker!r}')
        marker_map[key] = marker
    # Only remove a numeric segment if it is not part of a recognized
    # multi-token marker (e.g. NKX2_1). ASCL1 is intrinsically one token.
    is_numeric_marker_suffix = any(
        len(key) >= 2 and len(tokens) >= len(key)
        and tuple(t.casefold() for t in tokens[-len(key):]) == key
        for key in marker_map
    )
    if ending and len(tokens) >= 2 and tokens[-2].casefold() == 'field' and tokens[-1].isdecimal():
        tokens = tokens[:-2]
    elif ending and tokens and tokens[-1].isdecimal() and not is_numeric_marker_suffix:
        tokens.pop()
    maxlen = max(map(len, marker_map), default=1)
    out = []
    while tokens:
        match = None
        for n in range(min(maxlen, len(tokens)), 0, -1):
            key = tuple(t.casefold() for t in tokens[-n:])
            if key in marker_map:
                match = (n, marker_map[key])
                break
        if match is None:
            break
        n, marker = match
        out.append(marker)
        del tokens[-n:]
    if not out:
        fail(f'Cannot identify channel suffix in {filename!r}. '
             'Add its marker names to filename_markers or use explicit channel_names.')
    return out[::-1]


def ims_channel_count(source: Path) -> int | None:
    """Read just the channel-group count, without loading any image pixels.

    Metadata detection requires h5py, which is present in the ims-mosaic Conda
    environment. Without it, rely on the suffix vocabulary or explicit names;
    the field-job preflight checks the count on all four tiles regardless.
    """
    try:
        import h5py
    except ImportError:
        return None
    try:
        with h5py.File(source,'r') as handle:
            resolutions = handle['DataSet']
            res = next(iter(sorted(resolutions, key=lambda name: int(name.rsplit(' ',1)[-1]))))
            timepoints = resolutions[res]
            timepoint = next(iter(sorted(timepoints,key=lambda name: int(name.rsplit(' ',1)[-1]))))
            return len(timepoints[timepoint])
    except (OSError, KeyError, ValueError, StopIteration) as exc:
        # The full scientific preflight always validates real HDF5. A dry
        # filesystem test or incomplete source can still be scanned; refuse
        # execution at the metadata preflight rather than guessing the count.
        print(f'WARNING: cannot read IMS channel count from {source}: {exc}',flush=True)
        return None


def channels_from_filename(filename: str, markers: list[str] | tuple[str,...],
                           count: int | None) -> list[str]:
    """Use longest marker suffix, or the IMS channel count for unknown labels."""
    try:
        known = channel_suffix(filename,markers)
    except ValueError:
        known = []
    if count is None or len(known) == count:
        if known:
            return known
        fail(f'Cannot determine channel labels in {filename!r}; use an explicit '
             'channel_names list or add marker names to filename_markers')
    if count < 1 or count > 12:
        fail(f'Unexpected IMS channel count {count} for {filename!r}')
    stem = Path(filename).stem
    match = TILE_TOKEN.search(stem)
    before = stem[:match.start()] if match else stem
    tokens = [part for part in before.split('_') if part]
    # Do not remove a terminal numeric token that belongs to a recognized
    # multi-token marker such as NKX2_1.
    numeric_is_marker = any(
        len(marker.split('_')) >= 2
        and [x.casefold() for x in tokens[-len(marker.split('_')):]]
            == [x.casefold() for x in marker.split('_')]
        for marker in markers
    )
    if match and len(tokens)>=2 and tokens[-2].casefold()=='field' and tokens[-1].isdecimal():
        tokens = tokens[:-2]
    elif match and tokens and tokens[-1].isdecimal() and not numeric_is_marker:
        tokens.pop()
    if len(tokens)<count:
        fail(f'{filename!r}: fewer than {count} possible channel-name tokens')
    suffix = tokens[-count:]
    if any(not CHANNEL.fullmatch(name) or name.isdecimal() or name.casefold()=='field'
           for name in suffix):
        fail(f'{filename!r}: the {count}-channel IMS metadata does not unambiguously '
             f'match filename tokens {suffix}; add names to filename_markers '
             'for markers with underscores (e.g. NKX2_1) or set channel_names explicitly')
    # If a known suffix agrees with the proposed tail, check marker spelling.
    if known and [x.casefold() for x in suffix[-len(known):]] != [x.casefold() for x in known]:
        fail(f'{filename!r}: conflicting known filename suffix {known} and '
             f'{count}-channel token extraction {suffix}; set channel_names explicitly')
    return suffix


def optional_acquisition_identifier(
    value: str, markers: list[str] | tuple[str, ...], *, filename: bool
) -> int | None:
    """Return a standalone acquisition integer after the final marker.

    A trailing number that belongs to a declared multi-token marker, such as
    NKX2_1, is not an acquisition identifier. For filenames, only the text
    immediately before the final FNN tile token is considered.
    """
    stem = Path(value).stem if filename else value
    if filename:
        ending = TILE_TOKEN.search(stem)
        if ending is None:
            fail(f'Cannot find final FNN tile token in IMS filename {value!r}')
        stem = stem[:ending.start()]
    tokens = [part for part in stem.split('_') if part]
    if not tokens or not tokens[-1].isdecimal():
        return None
    numeric_is_marker = any(
        len(marker.split('_')) >= 2
        and len(tokens) >= len(marker.split('_'))
        and [x.casefold() for x in tokens[-len(marker.split('_')):]]
            == [x.casefold() for x in marker.split('_')]
        for marker in markers
    )
    return None if numeric_is_marker else int(tokens[-1])


def validate_acquisition_identifiers(
    field: 'Field', markers: list[str] | tuple[str, ...]
) -> int | None:
    """Require the optional directory and IMS filename identifiers to agree."""
    directory_id = optional_acquisition_identifier(
        field.group, markers, filename=False
    )
    filename_ids = [
        optional_acquisition_identifier(tile.path.name, markers, filename=True)
        for tile in field.tiles
    ]
    if any(identifier != filename_ids[0] for identifier in filename_ids[1:]):
        detail = ', '.join(
            f'F{tile.index:02d}={identifier!r}'
            for tile, identifier in zip(field.tiles, filename_ids)
        )
        fail(f'Acquisition identifier mismatch across IMS tiles in {field.path}: {detail}')
    filename_id = filename_ids[0]
    if directory_id != filename_id:
        fail(
            f'Acquisition identifier mismatch for {field.path}: '
            f'directory={directory_id!r}, IMS filenames={filename_id!r}'
        )
    return directory_id


def detect_field_channels(field: 'Field', markers: list[str] | tuple[str, ...]) -> list[str]:
    validate_acquisition_identifiers(field, markers)
    count = ims_channel_count(field.tiles[0].path)
    detected = []
    for tile in field.tiles:
        channels = channels_from_filename(tile.path.name, markers, count)
        detected.append(channels)
    expected = detected[0]
    for tile, got in zip(field.tiles[1:], detected[1:]):
        if got != expected:
            fail(f'Channel-name/order mismatch across tiles in {field.path}: '
                 f'F{field.tiles[0].index:02d}={expected}, F{tile.index:02d}={got}')
    return expected


def resolve_field_channels(field: 'Field', spec: dict[str, Any]) -> tuple[list[str], str]:
    requested = spec['channel_names']
    if requested == 'auto':
        return detect_field_channels(field, spec['_filename_markers']), 'filename'
    # Explicit override is allowed for nonstandard filenames, but compare to
    # recognisable suffixes whenever possible, to catch swapped channels.
    try:
        found = detect_field_channels(field, spec['_filename_markers'])
    except ValueError:
        found = None
    if found is not None and len(found) != len(requested):
        fail(f'{field.path}: explicit channel_names has {len(requested)} entries '
             f'but filenames describe {len(found)}: {found}')
    if found is not None and [x.casefold() for x in found] != [x.casefold() for x in requested]:
        print(f'WARNING: explicit channel_names {requested} override filename '
              f'labels {found} for {field.path}',flush=True)
    return list(requested), 'explicit_override'


def resolve_nuclear_channel(names: list[str], requested: str) -> str:
    if requested == 'auto':
        possible = [name for name in names if name.casefold() in NUCLEAR_MARKERS]
        if len(possible) != 1:
            fail(f'Cannot uniquely find a nuclear channel among {names}: '
                 'set nuclear_channel explicitly in the spec')
        return possible[0]
    found = [name for name in names if name.casefold() == requested.casefold()]
    if len(found) != 1:
        fail(f'nuclear_channel {requested!r} is not in detected channels {names}')
    return found[0]


def resolve_pairs(analysis: dict[str, Any], names: list[str], nuclear: str) -> dict[str, Any]:
    from itertools import combinations
    result = json.loads(json.dumps(analysis))
    markers = [name for name in names if name != nuclear]
    pairs = result.get('pairs', 'auto')
    if pairs == 'auto':
        result['pairs'] = [list(pair) for pair in combinations(markers,2)]
    elif not isinstance(pairs,list):
        fail('analysis.pairs must be "auto" or a list of channel-name pairs')
    else:
        for pair in pairs:
            if (not isinstance(pair,list) or len(pair)!=2 or pair[0]==pair[1]
                or any(name not in markers for name in pair)):
                fail(f'analysis.pairs entry {pair!r} is not available in '
                     f'detected non-nuclear markers {markers}; use "auto" for variable panels')
    return result



@dataclass(frozen=True)
class Tile:
    index: int
    path: Path


@dataclass(frozen=True)
class Field:
    path: Path
    group: str
    field_number: int
    tiles: tuple[Tile, ...]


def fail(message: str) -> None:
    raise ValueError(message)


def abs_path(value: str | Path) -> Path:
    return Path(value).expanduser().resolve()


def validate_storage_locations(root: Path, workspace: Path) -> None:
    """Keep persistent work on scratch for the protected Turbo acquisition tree."""
    source_is_protected = root == TURBO_IMS_ROOT or TURBO_IMS_ROOT in root.parents
    workspace_is_scratch = workspace == SCRATCH_ROOT or SCRATCH_ROOT in workspace.parents
    if source_is_protected and not workspace_is_scratch:
        fail(
            f'Inputs below {TURBO_IMS_ROOT} require a workspace below '
            f'{SCRATCH_ROOT}; got {workspace}'
        )


def safe_slug(text: str, limit: int = 75) -> str:
    slug = re.sub(r"[^A-Za-z0-9_-]+", "_", text).strip("_-")
    return (slug or "field")[:limit]


def sha_text(text: str, length: int = 10) -> str:
    return hashlib.sha256(text.encode('utf-8')).hexdigest()[:length]


def tile_id(path: Path) -> int | None:
    match = TILE_TOKEN.search(path.stem)
    return int(match.group(1)) if match else None


def should_skip_dir(name: str) -> bool:
    n = name.lower()
    return n.startswith('.') or n in IGNORE_NAMES or n.endswith('_ims_to_tiff') or n.endswith('_ims_to_tiff_output')


def discover_field_paths(root: Path, max_depth: int) -> list[tuple[Path,str,int]]:
    if not root.is_dir():
        fail(f"--input must be an existing directory: {root}")
    out: list[tuple[Path,str,int]] = []
    for current, dirs, files in os.walk(root, followlinks=False):
        path = Path(current)
        depth = len(path.relative_to(root).parts)
        # root may itself be an individual Field_* folder
        match = FIELD_DIR.fullmatch(path.name)
        if match:
            group = match.group('base') or path.parent.name
            out.append((path.resolve(), group, int(match.group('num'))))
            dirs[:] = []  # do not accidentally descend into tile folders
            continue
        if depth >= max_depth:
            dirs[:] = []
        else:
            dirs[:] = sorted([d for d in dirs if not should_skip_dir(d)], key=str.lower)
    return out


def ims_within_tile_dir(tile_dir: Path) -> Path:
    if not tile_dir.is_dir():
        fail(f"Missing IMS tile directory: {tile_dir}")
    matches = [p for p in tile_dir.rglob('*') if p.is_file() and p.suffix.lower() == '.ims'
               and not is_derived_stitched_ims(p)
               and not any(part.startswith('.') or part.casefold() in IGNORED_METADATA_FILES
                           for part in p.relative_to(tile_dir).parts)
               and len(p.relative_to(tile_dir).parts) <= 3]
    if len(matches) != 1:
        fail(f"Expected exactly one .ims within tile directory {tile_dir}, found {len(matches)}: {matches[:6]}")
    return matches[0].resolve()


def discover_tiles(field: Path, expected_count: int, contiguous: bool) -> tuple[Tile, ...]:
    tiles: dict[int, Path] = {}
    for child in sorted(field.iterdir(), key=lambda p: p.name.lower()):
        if is_os_metadata(child):
            continue
        # Ignore FusionStitcher *outputs* BEFORE parsing FNN or checking for
        # duplicate tile indices. E.g. the derived ..._FusionStitcher_F0.ims
        # is not a second copy of the true ..._F00.ims source tile.
        if child.is_file() and is_derived_stitched_ims(child):
            continue
        if child.is_dir() and is_derived_stitched_ims(child):
            continue
        if child.is_file() and (child.suffix.lower() == '.ims' or
                                (not child.suffix and tile_id(child) is not None)):
            idx = tile_id(child)
            if idx is None:
                fail(f"Unnumbered IMS candidate in {field}: {child.name}. Expected ..._F02.ims, f02.ims or f02")
            ims = child.resolve()
        elif child.is_dir() and not child.is_symlink() and (idx := tile_id(child)) is not None:
            ims = ims_within_tile_dir(child)
        else:
            continue
        if idx in tiles:
            fail(f"Duplicate tile index F{idx:02d} in {field}: {tiles[idx]} and {ims}")
        tiles[idx] = ims
    if len(tiles) != expected_count:
        fail(f"{field}: expected {expected_count} IMS tiles; found {len(tiles)}: {sorted(tiles)}")
    ordered = sorted(tiles)
    if contiguous and ordered != list(range(ordered[0], ordered[0]+expected_count)):
        fail(f"{field}: noncontiguous IMS tile indexes {ordered} (required by --require-contiguous)")
    return tuple(Tile(i, tiles[i]) for i in ordered)


def discover_fields(input_path: Path, max_depth: int, tile_count: int,
                    require_contiguous: bool, group_pattern: str | None,
                    field_numbers: set[int], *,
                    issues: list[dict[str, str]] | None = None,
                    scan_only: bool = False,
                    skip_incomplete: bool = True) -> list[Field]:
    """Select fields while optionally reporting incomplete mosaics.

    Scan: report malformed/incomplete fields and keep scanning read-only.
    Prepare: automatically exclude ONLY fields with fewer than required tiles.
    --strict-incomplete restores fail-fast behavior for missing tiles.
    Other invalid raw tiles (duplicates, bad filenames, extra tiles) stay fatal.
    """
    discovered = discover_field_paths(input_path, max_depth)
    if not discovered:
        fail(f"No Field_1 or <base>_Field_1 directories found under {input_path}. Check --max-depth.")
    out: list[Field] = []
    selected_count = 0
    for path, group, number in discovered:
        if group_pattern and not fnmatch.fnmatchcase(group.casefold(), group_pattern.casefold()):
            continue
        if field_numbers and number not in field_numbers:
            continue
        selected_count += 1
        try:
            tiles = discover_tiles(path, tile_count, require_contiguous)
        except ValueError as exc:
            message = str(exc)
            count_match = re.search(
                rf": expected {tile_count} IMS tiles; found (\d+): \[", message
            )
            incomplete = (count_match is not None and
                          int(count_match.group(1)) < tile_count)
            if scan_only or (skip_incomplete and incomplete):
                if issues is not None:
                    issues.append({
                        'path': str(path), 'group': group,
                        'field_number': str(number),
                        'kind': 'incomplete' if incomplete else 'invalid',
                        'reason': message,
                    })
                continue
            raise
        out.append(Field(path, group, number, tiles))
    if selected_count == 0:
        fail(f"No fields left after selection (group={group_pattern!r}, field={sorted(field_numbers)})")
    if not out and not scan_only:
        fail(f"No complete {tile_count}-tile fields selected; all selected fields are incomplete or invalid. Review the scan report.")
    out.sort(key=lambda f: (str(f.path.parent).casefold(), f.group.casefold(), f.field_number, str(f.path)))
    return out


def job_name(field: Field) -> str:
    # Uniqueness independent of whether --input was a session, group or field.
    # Keep a human readable group/field prefix; the hash disambiguates dates.
    return f"{safe_slug(field.group, 54)}_Field_{field.field_number}_{sha_text(str(field.path))}"


def file_stats(path: Path) -> dict[str, Any]:
    s = path.stat()
    return {"path": str(path), "size": s.st_size, "mtime_ns": s.st_mtime_ns}


def atomic_json(path: Path, obj: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + f'.tmp.{os.getpid()}')
    try:
        temp.write_text(json.dumps(obj, indent=2, sort_keys=True) + '\n', encoding='utf-8')
        temp.replace(path)
    finally:
        if temp.exists():
            temp.unlink()


def load_spec(spec_file: Path, *, require_fiji: bool = True) -> dict[str, Any]:
    raw = json.loads(spec_file.read_text(encoding='utf-8'))
    if not isinstance(raw, dict):
        fail("Spec must be a JSON object")
    required = {"backend", "channel_names", "nuclear_channel", "fiji_sif", "ims_unit",
                "trim_zero_padding", "stitching", "analysis_downsampling", "analysis"}
    missing = required - set(raw)
    if missing:
        fail(f"Missing required spec keys: {sorted(missing)}")
    if raw['backend'] not in VALID_PROFILES:
        fail("backend must be 'generic' or 'reporter'")
    names = raw['channel_names']
    auto = names == 'auto'
    if not auto and (not isinstance(names,list) or not (2 <= len(names) <= 12)
                     or any(not isinstance(n,str) or not CHANNEL.fullmatch(n) for n in names)
                     or len({n.casefold() for n in names}) != len(names)):
        fail('channel_names must be "auto" or 2-12 distinct valid identifier names in IMS wavelength order')
    if not isinstance(raw['nuclear_channel'],str) or (raw['nuclear_channel'] != 'auto' and not CHANNEL.fullmatch(raw['nuclear_channel'])):
        fail('nuclear_channel must be "auto" or a valid channel identifier')
    if raw['backend'] == 'reporter' and names != ['Hoechst','mNeonGreen','Virus','PV']:
        fail('The existing reporter IMS controller only accepts four channels: ["Hoechst", "mNeonGreen", "Virus", "PV"]. Use generic backend for variable marker panels.')
    if raw['backend'] == 'reporter' and raw['nuclear_channel'] != 'Hoechst':
        fail('Reporter mode requires nuclear_channel Hoechst')
    if raw['backend'] == 'generic' and not auto and 'Virus' in names:
        fail('Generic mode cannot use placeholder "Virus"; use BiVe3/BiVe4 or reporter mode')
    extra = raw.get('filename_markers',[])
    if not isinstance(extra,list) or any(not isinstance(x,str) or not CHANNEL.fullmatch(x) for x in extra):
        fail('filename_markers must be a list of channel identifiers')
    vocabulary = list(DEFAULT_FILENAME_MARKERS)
    lookup = {name.casefold() for name in vocabulary}
    for name in extra:
        if name.casefold() not in lookup:
            vocabulary.append(name)
            lookup.add(name.casefold())
    raw['_filename_markers'] = vocabulary
    if not isinstance(raw['trim_zero_padding'], bool):
        fail('trim_zero_padding must be boolean')
    if raw['ims_unit'] not in ('um', 'nm', 'mm', 'm'):
        fail('ims_unit must be um, nm, mm or m')
    ds = raw['analysis_downsampling']
    if not isinstance(ds, dict) or set(ds) != {'factor','method'} or isinstance(ds['factor'], bool) or not isinstance(ds['factor'],int) or ds['factor'] < 1 or ds['method'] != 'area_mean':
        fail('analysis_downsampling must be {"factor": positive integer, "method": "area_mean"}')
    stitch = raw['stitching']
    if not isinstance(stitch,dict) or 'grid' not in stitch or not (isinstance(stitch['grid'],list) and len(stitch['grid']) == 2 and all(isinstance(x,int) and x>0 for x in stitch['grid'])):
        fail('stitching.grid must be [positive_x,positive_y]')
    analysis = raw['analysis']
    if not isinstance(analysis,dict) or 'cellprofiler' not in analysis or 'thresholds' not in analysis:
        fail('analysis must include downstream CellProfiler configuration and thresholds')
    pair_setting = analysis.get('pairs','auto')
    if pair_setting != 'auto' and (not isinstance(pair_setting,list) or any(
        not isinstance(pair,list) or len(pair)!=2 for pair in pair_setting)):
        fail('analysis.pairs must be "auto" or an array of two-name arrays')
    if raw['backend']=='generic' and not auto:
        nuclear = resolve_nuclear_channel(names,raw['nuclear_channel'])
        resolve_pairs(analysis,names,nuclear)
    if raw['backend']=='reporter' and raw.get('virus_channel','auto') not in ('auto','BiVe3','BiVe4'):
        fail('virus_channel must be auto, BiVe3 or BiVe4')
    sif = str(raw['fiji_sif'])
    if require_fiji and (not sif.strip() or 'path/to' in sif or '<' in sif):
        fail('Set fiji_sif to your real Fiji/BigStitcher .sif path in the spec before preparation')
    p = Path(sif).expanduser()
    if not p.is_absolute():
        p = (spec_file.parent / p).resolve()
    raw['fiji_sif'] = str(p)
    return raw


def render_report(fields: list[Field], spec: dict[str, Any] | None = None,
                  issues: list[dict[str, str]] | None = None,
                  tolerate_channel_errors: bool = False) -> None:
    print(f"Found {len(fields)} organoid fields, {sum(len(f.tiles) for f in fields)} IMS tiles")
    by_group: dict[str,int] = {}
    for i, f in enumerate(fields):
        by_group[f.group] = by_group.get(f.group,0)+1
        ids = [t.index for t in f.tiles]
        staging = ', '.join(f'F{t.index:02d}→F{j:02d}' for j,t in enumerate(f.tiles))
        print(f'  [{i:03d}] {f.group} / Field_{f.field_number}  tile IDs={ids}  {staging}')
        print(f'        {f.path}')
    print('Group totals:', ', '.join(f'{k}: {v}' for k,v in sorted(by_group.items())))
    if issues:
        print(f"\nSKIPPED / UNPROCESSABLE FIELDS ({len(issues)}):")
        for issue in issues:
            print(f"  [{issue['kind'].upper()}] {issue['path']}")
            print(f"      {issue['reason']}")
        print('Incomplete fields are skipped automatically; other invalid fields stop preparation.')
    if spec:
        print(f"Profile: {spec['backend']} | channel setting: {spec['channel_names']} | grid: {spec['stitching']['grid']}")
        for field in fields:
            try:
                names, origin = resolve_field_channels(field,spec)
                nuclear = resolve_nuclear_channel(names,spec['nuclear_channel'])
                if len(names) < 2:
                    fail(f'{field.path}: full downstream workflow needs nuclear channel and at least one marker')
            except ValueError as exc:
                if not tolerate_channel_errors:
                    raise
                print(f'    CHANNEL ISSUE {field.group}/Field_{field.field_number}: {exc}')
                continue
            print(f'    Channel mapping {field.group}/Field_{field.field_number}: '
                  f'{names}, nuclear={nuclear} ({origin})')


def choose_fields(args: argparse.Namespace,
                  issues: list[dict[str, str]] | None = None) -> list[Field]:
    path = abs_path(args.input)
    tile_count = args.grid[0]*args.grid[1]
    return discover_fields(path, args.max_depth, tile_count, args.require_contiguous,
                           args.group, set(args.field or []),
                           issues=issues, scan_only=getattr(args, 'command', None) == 'scan',
                           skip_incomplete=not getattr(args, 'strict_incomplete', False))


def create_job_config(spec: dict[str, Any], job: dict[str, Any], repo: Path) -> dict[str, Any]:
    names = job['channel_names']
    nuclear = job['nuclear_channel']
    analysis = resolve_pairs(spec['analysis'],names,nuclear) if spec['backend'] == 'generic' else spec['analysis']
    if spec['backend'] == 'reporter':
        overrides = spec.get('virus_channel_by_sample', {})
        original_keys = [job['stage_sample_name'], job['group'], Path(job['source_field']).name, job['source_field']]
        matches = {overrides[k] for k in original_keys if k in overrides}
        if len(matches) > 1:
            fail(f'Conflicting per-sample virus overrides for {job["source_field"]}: {matches}')
        by_sample = {job['stage_sample_name']: next(iter(matches))} if matches else {}
        cfg: dict[str, Any] = {
            'ims_input': job['staging_dir'], 'input_mode': 'mosaic',
            'output_dir': job['output_dir'], 'fiji_sif': spec['fiji_sif'],
            'channel_order': spec['channel_names'],
            'virus_channel': spec.get('virus_channel', 'auto'),
            'virus_channel_by_sample': by_sample,
            'ims_unit': spec['ims_unit'],
            'trim_zero_padding': spec['trim_zero_padding'],
            'stitching': spec['stitching'],
            'analysis_downsampling': spec['analysis_downsampling'],
            'analysis': analysis,
        }
        return cfg
    # Generic mode has a separate config that the wrapper will split into the
    # original mosaic process and then generic stitched-TIFF runner.
    return {
        'profile': 'generic', 'stage_sample_name': job['stage_sample_name'],
        'staging_dir': job['staging_dir'], 'output_dir': job['output_dir'],
        'fiji_sif': spec['fiji_sif'], 'channel_names': names,
        'nuclear_channel': nuclear, 'ims_unit': spec['ims_unit'],
        'trim_zero_padding': spec['trim_zero_padding'],
        'stitching': spec['stitching'],
        'analysis_downsampling': spec['analysis_downsampling'],
        'analysis': analysis,
    }


def do_scan(args: argparse.Namespace) -> int:
    issues: list[dict[str, str]] = []
    fields = choose_fields(args, issues)
    spec = load_spec(abs_path(args.spec), require_fiji=False) if args.spec else None
    render_report(fields, spec, issues, tolerate_channel_errors=True)
    if not spec:
        for field in fields:
            try:
                names = detect_field_channels(field,DEFAULT_FILENAME_MARKERS)
                print(f'    Filename channels {field.group}/Field_{field.field_number}: {names}')
            except ValueError as exc:
                print(f'    Channels unverified {field.group}/Field_{field.field_number}: {exc}')
    print('\nScan only: no files were created. Use prepare after checking tile order.')
    if issues:
        n_incomplete = sum(i['kind'] == 'incomplete' for i in issues)
        n_invalid = len(issues) - n_incomplete
        print(f'WARNING: {n_incomplete} incomplete field(s) automatically skipped from a {args.grid[0]}x{args.grid[1]} submission.')
        if n_invalid:
            print(f'WARNING: {n_invalid} other invalid field(s) require correction before submission.')
    return 0



def generate_symlink_compatibility(repo: Path, workspace: Path) -> Path:
    """Build tiny controller-only adapters that preserve staged tile aliases.

    The original entrypoints call Path.resolve() before interpreting tile names,
    which dereferences an F00 symbolic link pointing at a source F02 file. This
    breaks the upstream F00-F03 validation and changes BigStitcher's tile order.
    The adapter changes ONLY path canonicalization in those two entrypoints.
    BaSiC, BigStitcher, registration, segmentation and matching are untouched.
    """
    import shutil
    dst = workspace / 'controller_compat'
    (dst / 'Stitching').mkdir(parents=True, exist_ok=True)
    edits = (
        ('run_ims_pipeline.py',
         '            path.resolve()\n            for path in directory.iterdir()',
         '            path.absolute()\n            for path in directory.iterdir()'),
        ('Stitching/IMSTIFF_trim_zero_padding.py',
         '    input_path = args.path.expanduser().resolve()\n',
         '    input_path = args.path.expanduser().absolute()\n'),
    )
    for relative, old, new in edits:
        source = repo / relative
        body = source.read_text(encoding='utf-8')
        if body.count(old) != 1:
            fail(f'Cannot generate compatibility adapter: expected one exact source fragment in {source}; repository may have changed.')
        target = dst / relative
        target.write_text(body.replace(old, new, 1), encoding='utf-8')
    for relative in ('run_stitched_pipeline.py', 'Stitching/tiff_safety.py'):
        source = (repo / relative).resolve()
        target = dst / relative
        target.symlink_to(source)
    return dst

def do_prepare(args: argparse.Namespace) -> int:
    root = abs_path(args.input)
    workspace = abs_path(args.workspace)
    repo = abs_path(args.repository)
    spec_path = abs_path(args.spec)
    spec = load_spec(spec_path)
    if not (repo/'run_ims_pipeline.py').is_file() or not (repo/'run_stitched_pipeline.py').is_file() or not (repo/'Stitching/IMSTIFF_trim_zero_padding.py').is_file():
        fail(f'Repository is missing original IMS/mosaic/stitched pipeline scripts: {repo}')
    grid = spec['stitching']['grid']
    if args.grid != grid:
        # User can specify --grid; one source of truth is the spec.
        fail(f'--grid {args.grid} differs from stitching.grid {grid} in spec')
    issues: list[dict[str, str]] = []
    fields = choose_fields(args, issues)
    render_report(fields, spec, issues)
    if workspace.exists() and any(workspace.iterdir()):
        fail(f'Workspace exists and is nonempty: {workspace}. Use a fresh workspace (or inspect the existing manifest).')
    if workspace == root or root in workspace.parents:
        fail('Workspace must be outside the source acquisition tree (avoid recursive discovery / source writes).')
    validate_storage_locations(root, workspace)
    workspace.mkdir(parents=True,exist_ok=True)
    compat_repo = generate_symlink_compatibility(repo, workspace)
    names = set()
    jobs = []
    for index, field in enumerate(fields):
        name = job_name(field)
        if name in names:
            fail(f'Internal job naming collision: {name}')
        names.add(name)
        stage_sample_name = name
        channel_names, naming_basis = resolve_field_channels(field,spec)
        nuclear_channel = resolve_nuclear_channel(channel_names,spec['nuclear_channel'])
        if len(channel_names)<2:
            fail(f'{field.path}: full downstream analysis needs a nuclear channel and at least one marker')
        if spec['backend']=='reporter' and len(channel_names)!=4:
            fail(f'{field.path}: reporter mode requires four IMS channels; use generic backend for {len(channel_names)}')
        staging = workspace/'input_staging'/name
        staging.mkdir(parents=True)
        tilemap = []
        for ordinal, tile in enumerate(field.tiles):
            target = staging/f'{stage_sample_name}_F{ordinal:02d}.ims'
            target.symlink_to(tile.path)
            tilemap.append({'source_index': tile.index, 'source': str(tile.path),
                            'source_stats': file_stats(tile.path), 'ordinal': ordinal,
                            'staged': str(target)})
        job = {
            'index': index, 'name': name, 'group': field.group,
            'field_number': field.field_number, 'source_field': str(field.path),
            'staging_dir': str(staging), 'stage_sample_name': stage_sample_name,
            'output_dir': str(workspace/'outputs'/name),
            'config_path': str(workspace/'configs'/f'{name}.json'),
            'channel_names': channel_names,
            'nuclear_channel': nuclear_channel,
            'channel_naming_basis': naming_basis,
            'tile_map': tilemap,
        }
        atomic_json(Path(job['config_path']), create_job_config(spec, job, repo))
        jobs.append(job)
    manifest = {
        'schema': 1, 'repository': str(repo), 'compat_repository': str(compat_repo), 'input_root': str(root),
        'batch_runner': str(SCRIPT_PATH), 'slurm_runner': str(SBATCH_PATH),
        'workspace': str(workspace), 'spec_path': str(spec_path),
        'spec_sha256': hashlib.sha256(spec_path.read_bytes()).hexdigest(),
        'backend': spec['backend'], 'channel_names_setting': spec['channel_names'],
        'jobs': jobs,
        'skipped_fields': issues,
        'eligible_field_count': len(jobs),
        'skipped_field_count': len(issues),
        'incomplete_field_policy': 'skip_and_record',
    }
    atomic_json(workspace/'manifest.json',manifest)
    # Human-readable audit of every omitted incomplete field, even if the list is empty.
    skip_csv = workspace / 'skipped_fields.csv'
    skip_csv_temp = workspace / 'skipped_fields.csv.tmp'
    try:
        with skip_csv_temp.open('w', newline='', encoding='utf-8') as handle:
            writer = csv.DictWriter(handle, fieldnames=['path', 'group', 'field_number', 'kind', 'reason'])
            writer.writeheader()
            writer.writerows(issues)
        skip_csv_temp.replace(skip_csv)
    finally:
        if skip_csv_temp.exists():
            skip_csv_temp.unlink()
    (workspace/'logs').mkdir(exist_ok=True)
    print(f'\nPrepared {len(jobs)} complete-field jobs in {workspace}; automatically skipped {len(issues)} incomplete field(s).')
    print(f'Skipped-field report: {workspace / "skipped_fields.csv"}')
    print(f'Manifest: {workspace / "manifest.json"}')
    print(f'Run one metadata preflight: python {SCRIPT_PATH} preflight --manifest {workspace / "manifest.json"} --index 0')
    print(f'Slurm first job: cd {workspace} && sbatch --array=0 {SBATCH_PATH} {workspace / "manifest.json"} {SCRIPT_PATH}')
    print(f'Slurm full batch: cd {workspace} && sbatch --array=0-{len(jobs)-1}%2 {SBATCH_PATH} {workspace / "manifest.json"} {SCRIPT_PATH}')
    return 0


def load_manifest(path: Path) -> dict[str,Any]:
    manifest = json.loads(path.read_text(encoding='utf-8'))
    if manifest.get('schema') != 1 or not isinstance(manifest.get('jobs'),list):
        fail(f'Invalid manifest schema in {path}')
    return manifest


def validate_job(job: dict[str,Any], manifest: dict[str,Any]) -> None:
    for tile in job['tile_map']:
        source = Path(tile['source'])
        staged = Path(tile['staged'])
        if not staged.is_symlink() or staged.resolve() != source:
            fail(f'Staged tile symlink mismatch/missing: {staged} -> {source}')
        if file_stats(source) != tile['source_stats']:
            fail(f'Source IMS file changed since prepare: {source}. Prepare a fresh workspace.')


def run_command(command: list[str], label: str) -> None:
    print(f'\n### {label}\n' + ' '.join(repr(x) if ' ' in x else x for x in command),flush=True)
    subprocess.run(command,check=True)


def generic_run(job: dict[str,Any], config: dict[str,Any], repo: Path,
                resume: bool, dry_run: bool) -> None:
    # Heavy scientific deps imported here, never by scan or prepare.
    import h5py
    import numpy as np
    import hdf5plugin  # noqa: F401
    import tifffile
    import run_ims_pipeline as ims
    import run_stitched_pipeline as stitched
    from Stitching.tiff_safety import atomic_tiff_write, is_complete_tiff

    stage = Path(job['staging_dir'])
    tiles = job['tile_map']
    first = Path(tiles[0]['staged'])
    sample = config['stage_sample_name']
    channels = list(config['channel_names'])
    factors = config['analysis_downsampling']
    factor = factors['factor']
    output = Path(config['output_dir'])
    preprocessing = output/'preprocessing'
    stitch = config['stitching']
    grid = stitch['grid']

    # Validate every tile before starting costly BaSiC/BigStitcher work. Do not
    # pretend the source IMS channel names can be trusted as biological labels.
    spacing = None
    reference_shapes = None
    reference_indices = None
    for tile in tiles:
        path = Path(tile['staged'])
        with h5py.File(path, 'r') as h5:
            ch_root, resolution, timepoint, keys = ims.first_resolution_timepoint(h5)
            if len(keys) != len(channels):
                fail(f'{path}: {len(keys)} channels, expected {len(channels)} in configured channel-index order')
            tile_spacing = ims.voxel_size_from_handle(h5,config['ims_unit'])
            shapes = []
            names = []
            stored_wavelengths = []
            for number, key in enumerate(keys):
                data = ch_root[key]['Data']
                if data.ndim != 3:
                    fail(f'{path}: channel {number} must be a ZYX 3-D stack; got {data.shape}')
                shapes.append(tuple(data.shape))
                raw_name = h5.get(f'DataSetInfo/Channel {number}')
                names.append(str(ims.decode_attribute(raw_name.attrs.get('Name',''))) if raw_name is not None else '')
                attrs = ({str(k):ims.decode_attribute(v) for k,v in raw_name.attrs.items()}
                         if raw_name is not None else {})
                _, wavelength = ims.wavelength_from_attributes(attrs)
                stored_wavelengths.append(wavelength)
            if reference_shapes is None:
                reference_shapes = shapes
                spacing = tile_spacing
                reference_indices = keys
            else:
                if shapes != reference_shapes:
                    fail(f'{path}: channel dimensions differ across mosaic tiles: {shapes} != {reference_shapes}')
                if not np.allclose(tile_spacing, spacing, rtol=1e-4, atol=1e-9):
                    fail(f'{path}: voxel calibration differs between tiles: {tile_spacing} != {spacing}')
            print(f'  Original F{tile["source_index"]:02d}: stored channel names {names}',flush=True)
            if all(wl is not None for wl in stored_wavelengths):
                print(f'    Wavelength attributes by channel index: {stored_wavelengths}',flush=True)
                if stored_wavelengths != sorted(stored_wavelengths):
                    print('    WARNING: IMS wavelength metadata are NOT in ascending order. '
                          'Confirm true wavelength/channel identity before stitching.',flush=True)
    assert spacing is not None
    print(f'IMS voxel calibration X/Y/Z (um): {spacing}',flush=True)
    print(f'Configured biological channel mapping: {list(enumerate(channels))}',flush=True)

    cmd = [sys.executable, str(repo/'Stitching/IMSTIFF_trim_zero_padding.py'),
           str(first), '--fiji-sif', str(config['fiji_sif']), '--output-root',str(preprocessing),
           '--channel-names',*channels,
           '--grid',str(grid[0]),str(grid[1]),
           '--overlap',str(stitch.get('overlap',10.0)),
           '--min-correlation',str(stitch.get('min_correlation',0.55)),
           '--basic-sample-planes',str(stitch.get('basic_sample_planes',12)),
           '--ims-unit',config['ims_unit'],
           '--intensity-downsampling',str(stitch.get('intensity_downsampling',32)),
           '--intensity-max-inliers',str(stitch.get('intensity_max_inliers',10000)),
           '--intensity-offset-only',str(stitch.get('intensity_offset_only',0.5)),
           '--intensity-unmodified',str(stitch.get('intensity_unmodified',0.5)),
           '--keep-previous']
    if not config['trim_zero_padding']:
        cmd.append('--no-trim-zero-padding')
    print(f'Mosaic input: {first} + {len(tiles)-1} sibling tiles')
    print(f'Mosaic process: IMS -> BaSiC -> BigStitcher -> projections')
    if dry_run:
        print('Generic preflight: IMS metadata and channel order checked; no mosaics written.')
        print('Planned stitching command: '+' '.join(cmd))
        return
    # Strict protection against accidental reuse of different scientific settings.
    # jobs are prepared once; run on the same workspace or --resume after a failure.
    run_command(cmd,'Stitch and flat-field correct')
    projections = preprocessing/sample/'04_max_projections'
    analysis_projections = preprocessing/sample/'05_analysis_projections'
    outputs: dict[str,Path] = {}
    downsampling_report = {}
    for channel in channels:
        source = projections/f'{sample}_{channel}_stitched_max_projection.tif'
        if not source.is_file():
            raise FileNotFoundError(f'Missing expected BigStitcher output: {source}')
        if factor == 1:
            outputs[channel] = source
            continue
        image = tifffile.imread(source)
        averaged, details = ims.area_mean_downsample(image,factor)
        target = analysis_projections/f'{sample}_{channel}_stitched_max_projection_downsample{factor}.tif'
        target.parent.mkdir(parents=True,exist_ok=True)
        if not is_complete_tiff(target,tuple(averaged.shape),averaged.dtype):
            atomic_tiff_write(target,averaged,bigtiff=True,compression='zlib',metadata={'axes':'YX'})
        outputs[channel] = target
        downsampling_report[channel] = details
    downstream = dict(config['analysis'])
    downstream.update({
        'nuclear_channel': config['nuclear_channel'],
        'channels': {k:str(v) for k,v in outputs.items()},
        'output_dir': str(output/'analysis'),
        'pixel_size_um': {'x':float(spacing[0])*factor,'y':float(spacing[1])*factor},
    })
    downstream_config = Path(job['config_path']).with_suffix('.stitched.json')
    if downstream_config.exists() and json.loads(downstream_config.read_text()) != downstream:
        fail(f'Downstream config changed for an existing run: {downstream_config}. Use new workspace.')
    atomic_json(downstream_config,downstream)
    parsed = stitched.read_config(downstream_config)
    stitched.validate_inputs(parsed)
    atomic_json(output/'folder_batch_provenance.json',{
        'source_field':job['source_field'],'tile_map':job['tile_map'],
        'channels':channels,'channel_naming_basis':job.get('channel_naming_basis'),
        'voxel_size_um':list(spacing),
        'analysis_downsampling':factors,'downsampling_geometry':downsampling_report,
        'generated_stitched_config':str(downstream_config),
    })
    downstream_cmd=[sys.executable,str(repo/'run_stitched_pipeline.py'),'--config',str(downstream_config)]
    if resume:
        downstream_cmd.append('--resume')
    elif (output/'analysis'/'run_manifest.json').exists():
        # Automatic resume for a repeated Slurm invocation with unchanged config.
        downstream_cmd.append('--resume')
    run_command(downstream_cmd,'Existing complete downstream analysis')


def do_job(args:argparse.Namespace, dry_run:bool) -> int:
    manifest_path=abs_path(args.manifest)
    manifest=load_manifest(manifest_path)
    i=int(args.index)
    if i<0 or i>=len(manifest['jobs']):
        fail(f'Invalid index {i}; valid 0..{len(manifest["jobs"])-1}')
    repo=Path(manifest['compat_repository'])
    job=manifest['jobs'][i]
    validate_job(job,manifest)
    config=json.loads(Path(job['config_path']).read_text(encoding='utf-8'))
    print(f'\nJob [{i}/{len(manifest["jobs"])}]: {job["group"]} Field_{job["field_number"]}',flush=True)
    print('Tile map:',flush=True)
    for tile in job['tile_map']:
        print(f'  Original F{tile["source_index"]:02d} {tile["source"]} => staged F{tile["ordinal"]:02d}',flush=True)
    if manifest['backend']=='reporter':
        cmd=[sys.executable,str(repo/'run_ims_pipeline.py'),'--config',job['config_path']]
        if dry_run:
            cmd.append('--dry-run')
        elif args.resume or (Path(job['output_dir'])/'ims_run_manifest.json').exists():
            cmd.append('--resume')
        run_command(cmd,'Original full IMS controller')
    else:
        if str(repo) not in sys.path:
            sys.path.insert(0,str(repo))
        generic_run(job,config,repo,args.resume,dry_run)
    return 0


def do_inspect(args:argparse.Namespace) -> int:
    manifest=load_manifest(abs_path(args.manifest))
    print(f"Workspace: {manifest['workspace']}  Backend: {manifest['backend']}")
    for job in manifest['jobs']:
        dst=Path(job['output_dir'])
        expected = (dst/'ims_run_manifest.json') if manifest['backend']=='reporter' else (dst/'analysis'/'run_manifest.json')
        status='COMPLETE/EXISTS' if expected.is_file() else 'NOT COMPLETE'
        print(f'[{job["index"]:03d}] {job["group"]} Field_{job["field_number"]}: {status} (manifest {expected})')
    return 0


def add_selection(p:argparse.ArgumentParser) -> None:
    p.add_argument('--input',required=True,help='Olti_* directory, group directory or one Field_* directory')
    p.add_argument('--group',help='Optional group name glob, e.g. Exp17_* or *LHX6_PV*')
    p.add_argument('--field',type=int,action='append',help='Include only Field_N; repeat for multiple')
    p.add_argument('--max-depth',type=int,default=5,help='Maximum nested directory levels to examine (default 5)')
    p.add_argument('--grid',type=int,nargs=2,default=[2,2],metavar=('X','Y'))
    p.add_argument('--require-contiguous',action='store_true',help='Require original tile IDs to be contiguous; default stages any 4 ascending indexes')


def build_parser() -> argparse.ArgumentParser:
    p=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    sub=p.add_subparsers(dest='command',required=True)
    scan=sub.add_parser('scan',help='Inspect hierarchy, groups, fields and ascending IMS tiles; no writes')
    add_selection(scan)
    scan.add_argument('--spec',help='Optional profile to validate channel names during read-only scan')
    prep=sub.add_parser('prepare',help='Create staged symlinks, per-field configs, and Slurm manifest')
    add_selection(prep)
    incomplete_mode = prep.add_mutually_exclusive_group()
    incomplete_mode.add_argument('--strict-incomplete', action='store_true', help='Abort if any field has fewer than the required IMS tiles (default: skip and report)')
    incomplete_mode.add_argument('--skip-incomplete', action='store_true', help='Legacy alias; incomplete fields are already skipped and reported by default')
    prep.add_argument('--spec',required=True,help='Analysis profile JSON (generic or reporter)')
    prep.add_argument('--workspace',required=True,help='New empty work/output root, separate from source IMS tree')
    prep.add_argument('--repository',default=str(REPOSITORY_ROOT),help='ImageAnalysisScripts-Copy checkout directory')
    for kind in ('preflight','run'):
        action=sub.add_parser(kind,help=f'{kind} one manifest job by index')
        action.add_argument('--manifest',required=True)
        action.add_argument('--index',required=True,type=int)
        if kind=='run':
            action.add_argument('--resume',action='store_true')
    show=sub.add_parser('inspect',help='Inspect generated manifest & completion status')
    show.add_argument('--manifest',required=True)
    return p


def main(argv:list[str]|None=None)->int:
    parser=build_parser()
    args=parser.parse_args(argv)
    try:
        if args.command=='scan':return do_scan(args)
        if args.command=='prepare':return do_prepare(args)
        if args.command=='inspect':return do_inspect(args)
        if args.command=='preflight':
            args.resume=False
            return do_job(args,True)
        if args.command=='run':return do_job(args,False)
    except (ValueError,FileNotFoundError,KeyError,IndexError) as exc:
        parser.exit(2,f'ERROR: {exc}\n')
    return 0

if __name__=='__main__':
    raise SystemExit(main())
