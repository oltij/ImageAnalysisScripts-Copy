#!/usr/bin/env python3
"""Scan, prepare, validate, run, resume, inspect, or submit microscopy inputs."""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
import re
import subprocess
import sys
from pathlib import Path
from typing import Any

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from GeneralPipeline import batch_runner  # noqa: E402
from GeneralPipeline.channel_metadata import (  # noqa: E402
    CHANNEL_NAME,
    marker_vocabulary,
    validate_channel_names,
)
from GeneralPipeline.discovery import (  # noqa: E402
    RawAcquisition,
    StitchedAcquisition,
    discover_raw_acquisitions,
    discover_stitched_acquisitions,
)
from GeneralPipeline.manifest import (  # noqa: E402
    MANIFEST_KIND,
    SCHEMA_VERSION,
    absolute,
    atomic_json,
    copy_or_link,
    create_controller_compatibility,
    fail,
    file_stats,
    fresh_workspace,
    load_manifest,
    repository_git_state,
    require_reproducible_checkout,
    scientific_hashes,
    sha256_file,
    validate_manifest,
    validate_storage_locations,
)


SCRIPT_PATH = Path(__file__).resolve()
GENERAL_DIR = SCRIPT_PATH.parent
REPOSITORY_ROOT = GENERAL_DIR.parent
SBATCH_PATH = GENERAL_DIR / "slurm/run_pipeline.sbatch"
SUPPORTED_MODES = ("raw_ims_tiles", "stitched_tiffs")


def _downstream_module(repository: Path) -> Any:
    import importlib.util

    path = repository / "run_stitched_pipeline.py"
    spec = importlib.util.spec_from_file_location("general_pipeline_downstream", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import shared downstream runner: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _positive_pixels(value: Any, label: str) -> dict[str, float]:
    if not isinstance(value, dict) or set(value) != {"x", "y"}:
        fail(f"{label} must contain exactly x and y")
    result: dict[str, float] = {}
    for axis in ("x", "y"):
        raw = value[axis]
        if isinstance(raw, bool):
            fail(f"{label}.{axis} must be a finite positive number")
        try:
            number = float(raw)
        except (TypeError, ValueError) as exc:
            raise ValueError(
                f"{label}.{axis} must be a finite positive number"
            ) from exc
        if not math.isfinite(number) or number <= 0:
            fail(f"{label}.{axis} must be a finite positive number")
        result[axis] = number
    return result


def _validate_analysis(analysis: Any, repository: Path) -> dict[str, Any]:
    if not isinstance(analysis, dict):
        fail("analysis must be an object containing the existing downstream settings")
    required = {"cellprofiler", "alignment", "thresholds"}
    if missing := required - set(analysis):
        fail(f"analysis is missing required settings: {sorted(missing)}")
    managed = {"channels", "nuclear_channel", "pixel_size_um", "output_dir"}
    if conflict := managed & set(analysis):
        fail(f"analysis contains runner-managed settings: {sorted(conflict)}")
    cellprofiler = analysis["cellprofiler"]
    if not isinstance(cellprofiler, dict):
        fail("analysis.cellprofiler must be an object")
    pipeline = cellprofiler.get("pipeline", "CellProfiler/MGEOPVFinal.cppipe")
    if not isinstance(pipeline, str) or not pipeline.strip():
        fail("analysis.cellprofiler.pipeline must be a nonempty path")
    pipeline_path = Path(pipeline).expanduser()
    if not pipeline_path.is_absolute():
        pipeline_path = repository / pipeline_path
    downstream = _downstream_module(repository)
    downstream.validate_cellprofiler_pipeline(pipeline_path.resolve())
    unknown_cellprofiler = set(cellprofiler) - {
        "pipeline",
        "conda_env",
        "conda_executable",
    }
    if unknown_cellprofiler:
        fail(f"Unknown analysis.cellprofiler settings: {sorted(unknown_cellprofiler)}")
    alignment = analysis["alignment"]
    if not isinstance(alignment, dict) or not isinstance(
        alignment.get("enabled", True), bool
    ):
        fail("analysis.alignment must contain a boolean enabled setting")
    extra_args = alignment.get("extra_args", [])
    if not isinstance(extra_args, list) or not all(
        isinstance(value, str) for value in extra_args
    ):
        fail("analysis.alignment.extra_args must be a list of strings")
    protected_args = {
        "--fixed-image",
        "--moving-image",
        "--fixed-csv",
        "--moving-csv",
        "--output",
        "--pixel-size-x",
        "--pixel-size-y",
        "--fixed-x-col",
        "--fixed-y-col",
        "--moving-x-col",
        "--moving-y-col",
    }
    if any(
        value in protected_args
        or any(value.startswith(flag + "=") for flag in protected_args)
        for value in extra_args
    ):
        fail("analysis.alignment.extra_args cannot override runner-managed arguments")
    pairs = analysis.get("pairs", "auto")
    if pairs != "auto" and (
        not isinstance(pairs, list)
        or any(not isinstance(pair, list) or len(pair) != 2 for pair in pairs)
    ):
        fail('analysis.pairs must be "auto" or a list of two-channel lists')
    thresholds = analysis["thresholds"]
    if not isinstance(thresholds, dict):
        fail("analysis.thresholds must be an object")
    unknown_thresholds = set(thresholds) - {"marker_nuclear", "marker_marker"}
    if unknown_thresholds:
        fail(f"Unknown analysis thresholds: {sorted(unknown_thresholds)}")
    for label, value in thresholds.items():
        if (
            isinstance(value, bool)
            or not isinstance(value, (int, float))
            or not math.isfinite(float(value))
            or not 0 <= value <= 1
        ):
            fail(f"analysis.thresholds.{label} must lie in [0, 1]")
    filters = analysis.get("filters", {})
    if not isinstance(filters, dict) or set(filters) - {"intensity", "shape"}:
        fail("analysis.filters may contain only intensity and shape")
    downstream.filters_for_pass(filters, "intensity")
    downstream.filters_for_pass(filters, "shape")
    reports = analysis.get("reports", {})
    allowed_reports = {"organoid_mask", "intensity", "shape", "colocalization"}
    if not isinstance(reports, dict) or set(reports) - allowed_reports:
        fail(f"analysis.reports may contain only {sorted(allowed_reports)}")
    if any(not isinstance(value, bool) for value in reports.values()):
        fail("Every analysis.reports value must be boolean")
    constrained_matching = analysis.get("nucleus_constrained_marker_matching", True)
    if not isinstance(constrained_matching, bool):
        fail("analysis.nucleus_constrained_marker_matching must be boolean")
    return json.loads(json.dumps(analysis))


def _resolve_pairs(
    analysis: dict[str, Any], channels: list[str], nuclear: str
) -> dict[str, Any]:
    result = json.loads(json.dumps(analysis))
    markers = [name for name in channels if name != nuclear]
    pairs = result.get("pairs", "auto")
    if pairs == "auto":
        result["pairs"] = [list(pair) for pair in itertools.combinations(markers, 2)]
    else:
        for pair in pairs:
            if pair[0] == pair[1] or any(name not in markers for name in pair):
                fail(
                    f"analysis.pairs entry {pair!r} is not a pair of available non-nuclear channels {markers}"
                )
        if len({tuple(pair) for pair in pairs}) != len(pairs):
            fail("analysis.pairs contains a duplicate directional pair")
    thresholds = result["thresholds"]
    if markers and "marker_nuclear" not in thresholds:
        fail(
            "analysis.thresholds.marker_nuclear is required when marker channels exist"
        )
    if result["pairs"] and "marker_marker" not in thresholds:
        fail("analysis.thresholds.marker_marker is required when marker pairs exist")
    return result


def load_config(path: Path, mode: str, repository: Path) -> dict[str, Any]:
    raw = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(raw, dict):
        fail("Configuration must be a JSON object")
    if raw.get("input_mode") != mode:
        fail(
            f"Configuration input_mode={raw.get('input_mode')!r} does not match --mode {mode!r}"
        )
    names = raw.get("channel_names")
    if mode == "raw_ims_tiles" and names == "auto":
        normalized_names: str | list[str] = "auto"
    elif isinstance(names, list):
        normalized_names = validate_channel_names(names)
    else:
        fail("channel_names must be a list (or auto for raw_ims_tiles)")
    nuclear = raw.get("nuclear_channel")
    if not isinstance(nuclear, str) or (
        nuclear != "auto" and not CHANNEL_NAME.fullmatch(nuclear)
    ):
        fail("nuclear_channel must be auto or a valid channel name")
    if mode == "stitched_tiffs":
        assert isinstance(normalized_names, list)
        if nuclear == "auto" or nuclear not in normalized_names:
            fail(
                "stitched_tiffs requires an explicit nuclear_channel present in channel_names"
            )
    extra_markers = raw.get("filename_markers", [])
    if not isinstance(extra_markers, list):
        fail("filename_markers must be a list")
    if isinstance(normalized_names, list):
        extra_markers = [*extra_markers, *normalized_names]
    analysis = _validate_analysis(raw.get("analysis"), repository)
    result = dict(raw)
    result["channel_names"] = normalized_names
    result["nuclear_channel"] = nuclear
    result["filename_markers"] = marker_vocabulary(extra_markers)
    result["analysis"] = analysis
    pipeline_value = analysis["cellprofiler"].get(
        "pipeline", "CellProfiler/MGEOPVFinal.cppipe"
    )
    pipeline_path = Path(pipeline_value).expanduser()
    if not pipeline_path.is_absolute():
        pipeline_path = repository / pipeline_path
    result["_cellprofiler_pipeline"] = str(pipeline_path.resolve())
    result["_config_path"] = str(path)

    if mode == "raw_ims_tiles":
        discovery = raw.get("raw_ims_tiles")
        preprocessing = raw.get("preprocessing")
        if not isinstance(discovery, dict) or not isinstance(preprocessing, dict):
            fail("raw_ims_tiles mode requires raw_ims_tiles and preprocessing objects")
        grid = discovery.get("grid", [2, 2])
        if grid != [2, 2]:
            fail("Only the validated raw IMS grid [2, 2] is supported")
        if (
            not isinstance(discovery.get("max_depth", 5), int)
            or discovery.get("max_depth", 5) < 1
        ):
            fail("raw_ims_tiles.max_depth must be a positive integer")
        for flag in ("require_contiguous_tiles", "skip_incomplete"):
            if not isinstance(discovery.get(flag, flag == "skip_incomplete"), bool):
                fail(f"raw_ims_tiles.{flag} must be boolean")
        sessions = discovery.get("sessions", [])
        if not isinstance(sessions, list) or any(
            not isinstance(value, str) or not value for value in sessions
        ):
            fail("raw_ims_tiles.sessions must be a list of direct-child names")
        if len({value.casefold() for value in sessions}) != len(sessions):
            fail("raw_ims_tiles.sessions contains duplicates")
        group_pattern = discovery.get("group_pattern")
        if group_pattern is not None and not isinstance(group_pattern, str):
            fail("raw_ims_tiles.group_pattern must be a string or null")
        fiji = preprocessing.get("fiji_sif")
        if not isinstance(fiji, str) or not fiji.strip() or "path/to" in fiji:
            fail("preprocessing.fiji_sif must be the real Fiji/BigStitcher SIF path")
        fiji_path = Path(fiji).expanduser()
        if not fiji_path.is_absolute():
            fiji_path = (path.parent / fiji_path).resolve()
        unit = preprocessing.get("ims_unit")
        if unit not in {"um", "nm", "mm", "m"}:
            fail("preprocessing.ims_unit must be um, nm, mm, or m")
        if not isinstance(preprocessing.get("trim_zero_padding"), bool):
            fail("preprocessing.trim_zero_padding must be boolean")
        downsampling = preprocessing.get("analysis_downsampling")
        if (
            not isinstance(downsampling, dict)
            or set(downsampling) != {"factor", "method"}
            or isinstance(downsampling["factor"], bool)
            or not isinstance(downsampling["factor"], int)
            or downsampling["factor"] < 1
            or downsampling["method"] != "area_mean"
        ):
            fail(
                'preprocessing.analysis_downsampling must be {"factor": positive integer, "method": "area_mean"}'
            )
        stitching = preprocessing.get("stitching")
        if not isinstance(stitching, dict) or stitching.get("grid") != grid:
            fail("preprocessing.stitching.grid must exactly match raw_ims_tiles.grid")
        result["raw_ims_tiles"] = {
            "grid": grid,
            "max_depth": discovery.get("max_depth", 5),
            "require_contiguous_tiles": discovery.get(
                "require_contiguous_tiles", False
            ),
            "skip_incomplete": discovery.get("skip_incomplete", True),
            "sessions": sessions,
            "group_pattern": group_pattern,
            "channel_names": normalized_names,
            "nuclear_channel": nuclear,
            "filename_markers": result["filename_markers"],
            "ims_unit": unit,
        }
        result["preprocessing"] = {
            **preprocessing,
            "fiji_sif": str(fiji_path),
        }
    else:
        stitched = raw.get("stitched_tiffs")
        if not isinstance(stitched, dict) or not isinstance(
            stitched.get("mapping"), dict
        ):
            fail("stitched_tiffs mode requires stitched_tiffs.mapping")
        result["pixel_size_um"] = _positive_pixels(
            raw.get("pixel_size_um"), "pixel_size_um"
        )
        result["stitched_tiffs"] = {"mapping": stitched["mapping"]}
    return result


def _safe_slug(value: str, limit: int = 70) -> str:
    slug = re.sub(r"[^A-Za-z0-9_-]+", "_", value).strip("_-")
    return (slug or "acquisition")[:limit]


def _unique_name(label: str, identity: str) -> str:
    digest = hashlib.sha256(identity.encode("utf-8")).hexdigest()[:10]
    return f"{_safe_slug(label)}_{digest}"


def _discover(
    mode: str, root: Path, config: dict[str, Any], repository: Path
) -> tuple[list[RawAcquisition] | list[StitchedAcquisition], list[dict[str, Any]]]:
    issues: list[dict[str, Any]] = []
    if mode == "raw_ims_tiles":
        return discover_raw_acquisitions(
            root, config["raw_ims_tiles"], issues=issues
        ), issues
    downstream = _downstream_module(repository)
    settings = {
        "channel_names": config["channel_names"],
        "nuclear_channel": config["nuclear_channel"],
        "pixel_size_um": config["pixel_size_um"],
        "mapping": config["stitched_tiffs"]["mapping"],
        "analysis": config["analysis"],
    }
    return (
        discover_stitched_acquisitions(
            root,
            settings,
            downstream.inspect_tiff,
            manifest_base=Path(config["_config_path"]).parent,
        ),
        issues,
    )


def scan(mode: str, input_root: Path, config: dict[str, Any], repository: Path) -> int:
    acquisitions, issues = _discover(mode, input_root, config, repository)
    print(f"Found {len(acquisitions)} valid {mode} acquisition(s)")
    for index, acquisition in enumerate(acquisitions):
        if isinstance(acquisition, RawAcquisition):
            indices = [tile.source_index for tile in acquisition.tiles]
            print(
                f"[{index:03d}] {acquisition.name} / Field_{acquisition.field_number}: "
                f"tiles={indices}, channels={list(acquisition.channels)}, "
                f"nuclear={acquisition.nuclear_channel}, session={acquisition.session!r}"
            )
            physical = acquisition.physical_metadata
            print(
                f"      physical indices={physical['indices']}, "
                f"stored_names={physical['stored_names']}, "
                f"wavelengths={physical['wavelengths']}, "
                f"voxel_size_um={physical['voxel_size_um']}"
            )
        else:
            shapes = {
                name: tuple(metadata["yx_shape"])
                for name, metadata in acquisition.tiff_metadata.items()
            }
            print(
                f"[{index:03d}] {acquisition.name}: channels={list(acquisition.channels)}, "
                f"nuclear={acquisition.nuclear_channel}, XY={shapes}, "
                f"pixel_size_um={acquisition.pixel_size_um}"
            )
    for issue in issues:
        print(f"[SKIPPED {issue['kind'].upper()}] {issue['path']}: {issue['reason']}")
    print("Scan only: no files were created")
    return 0


def _raw_job(
    acquisition: RawAcquisition,
    index: int,
    workspace: Path,
    repository: Path,
    config: dict[str, Any],
) -> dict[str, Any]:
    identity = str(acquisition.path)
    name = _unique_name(
        f"{acquisition.name}_Field_{acquisition.field_number}", identity
    )
    staging = workspace / "input_staging" / name
    sources: list[dict[str, Any]] = []
    for ordinal, tile in enumerate(acquisition.tiles):
        target = staging / f"{name}_F{ordinal:02d}.ims"
        copy_or_link(tile.path, target)
        sources.append(
            {
                "source": str(tile.path),
                "source_stats": file_stats(tile.path),
                "source_tile_index": tile.source_index,
                "staged_tile_ordinal": ordinal,
                "staged": str(target),
            }
        )
    output = workspace / "outputs" / name
    generated = workspace / "configs" / f"{name}.json"
    job_config = {
        "input_mode": "raw_ims_tiles",
        "stage_sample_name": name,
        "staging_dir": str(staging),
        "output_dir": str(output),
        "fiji_sif": config["preprocessing"]["fiji_sif"],
        "channel_names": list(acquisition.channels),
        "nuclear_channel": acquisition.nuclear_channel,
        "ims_unit": config["preprocessing"]["ims_unit"],
        "trim_zero_padding": config["preprocessing"]["trim_zero_padding"],
        "stitching": config["preprocessing"]["stitching"],
        "analysis_downsampling": config["preprocessing"]["analysis_downsampling"],
        "analysis": _resolve_pairs(
            config["analysis"], list(acquisition.channels), acquisition.nuclear_channel
        ),
    }
    atomic_json(generated, job_config)
    return {
        "index": index,
        "name": name,
        "display_name": acquisition.name,
        "input_mode": "raw_ims_tiles",
        "source_acquisition": str(acquisition.path),
        "session": acquisition.session,
        "field_number": acquisition.field_number,
        "experimental_set_number": acquisition.experimental_set_number,
        "acquisition_identifier": acquisition.acquisition_identifier,
        "channel_names": list(acquisition.channels),
        "nuclear_channel": acquisition.nuclear_channel,
        "channel_naming_basis": acquisition.channel_naming_basis,
        "physical_metadata": acquisition.physical_metadata,
        "sources": sources,
        "staging_dir": str(staging),
        "output_dir": str(output),
        "config_path": str(generated),
        "config_sha256": sha256_file(generated),
        "downstream_runner": str(repository / "run_stitched_pipeline.py"),
        "preprocessing_stages": [
            "IMS export",
            "BaSiC",
            "BigStitcher",
            "maximum projection",
        ],
        "bypassed_stages": [],
    }


def _stitched_job(
    acquisition: StitchedAcquisition,
    index: int,
    workspace: Path,
    repository: Path,
    config: dict[str, Any],
) -> dict[str, Any]:
    identity = "\0".join(str(path) for path in acquisition.channels.values())
    name = _unique_name(acquisition.name, identity)
    staging = workspace / "input_staging" / name
    sources: list[dict[str, Any]] = []
    staged_channels: dict[str, str] = {}
    for channel, source in acquisition.channels.items():
        suffix = source.suffix.casefold()
        target = staging / f"{channel}{suffix}"
        copy_or_link(source, target)
        staged_channels[channel] = str(target)
        sources.append(
            {
                "source": str(source),
                "source_stats": file_stats(source),
                "channel": channel,
                "staged": str(target),
                "tiff_metadata": acquisition.tiff_metadata[channel],
            }
        )
    output = workspace / "outputs" / name
    generated = workspace / "configs" / f"{name}.json"
    job_config = _resolve_pairs(
        config["analysis"], list(acquisition.channels), acquisition.nuclear_channel
    )
    job_config.update(
        {
            "nuclear_channel": acquisition.nuclear_channel,
            "channels": staged_channels,
            "pixel_size_um": acquisition.pixel_size_um,
            "output_dir": str(output / "analysis"),
        }
    )
    atomic_json(generated, job_config)
    downstream = _downstream_module(repository)
    parsed = downstream.read_config(generated)
    downstream.validate_inputs(parsed)
    return {
        "index": index,
        "name": name,
        "display_name": acquisition.name,
        "input_mode": "stitched_tiffs",
        "source_acquisition": acquisition.name,
        "channel_names": list(acquisition.channels),
        "nuclear_channel": acquisition.nuclear_channel,
        "pixel_size_um": acquisition.pixel_size_um,
        "sources": sources,
        "staging_dir": str(staging),
        "output_dir": str(output),
        "config_path": str(generated),
        "config_sha256": sha256_file(generated),
        "downstream_runner": str(repository / "run_stitched_pipeline.py"),
        "preprocessing_stages": [],
        "bypassed_stages": ["IMS export", "BaSiC", "BigStitcher"],
    }


def prepare(
    mode: str,
    input_root: Path,
    config_path: Path,
    workspace: Path,
    repository: Path,
    *,
    allow_dirty_repository: bool,
) -> int:
    config = load_config(config_path, mode, repository)
    acquisitions, issues = _discover(mode, input_root, config, repository)
    validate_storage_locations(input_root, workspace)
    for acquisition in acquisitions:
        sources = (
            [tile.path for tile in acquisition.tiles]
            if isinstance(acquisition, RawAcquisition)
            else list(acquisition.channels.values())
        )
        for source in sources:
            validate_storage_locations(source, workspace)
    git_state = repository_git_state(repository)
    require_reproducible_checkout(git_state, allow_dirty_repository)
    fresh_workspace(workspace)
    compatibility: Path | None = None
    if mode == "raw_ims_tiles":
        compatibility = create_controller_compatibility(repository, workspace)
    jobs: list[dict[str, Any]] = []
    for index, acquisition in enumerate(acquisitions):
        if isinstance(acquisition, RawAcquisition):
            jobs.append(_raw_job(acquisition, index, workspace, repository, config))
        else:
            jobs.append(
                _stitched_job(acquisition, index, workspace, repository, config)
            )
    (workspace / "logs").mkdir(exist_ok=True)
    compatibility_hashes = (
        {
            relative: sha256_file(compatibility / relative)
            for relative in (
                "run_ims_pipeline.py",
                "run_stitched_pipeline.py",
                "Stitching/IMSTIFF_trim_zero_padding.py",
                "Stitching/tiff_safety.py",
            )
        }
        if compatibility is not None
        else {}
    )
    payload = {
        "schema": SCHEMA_VERSION,
        "kind": MANIFEST_KIND,
        "input_mode": mode,
        "input_root": str(input_root),
        "workspace": str(workspace),
        "repository": str(repository),
        "repository_git": git_state,
        "allow_dirty_repository": allow_dirty_repository,
        "config_path": str(config_path),
        "config_sha256": sha256_file(config_path),
        "pipeline_runner": str(SCRIPT_PATH),
        "pipeline_runner_sha256": sha256_file(SCRIPT_PATH),
        "slurm_runner": str(SBATCH_PATH),
        "slurm_runner_sha256": sha256_file(SBATCH_PATH),
        "cellprofiler_pipeline": config["_cellprofiler_pipeline"],
        "cellprofiler_pipeline_sha256": sha256_file(
            Path(config["_cellprofiler_pipeline"])
        ),
        "compat_repository": str(compatibility) if compatibility else None,
        "compatibility_sha256": compatibility_hashes,
        "scientific_sha256": scientific_hashes(repository),
        "jobs": jobs,
        "skipped_acquisitions": issues,
    }
    atomic_json(workspace / "manifest.json", payload)
    print(f"Prepared {len(jobs)} {mode} acquisition(s) in {workspace}")
    if issues:
        print(
            f"Skipped {len(issues)} incomplete acquisition(s); details are recorded in the manifest"
        )
    print(f"Manifest: {workspace / 'manifest.json'}")
    print(
        f"Preflight first acquisition: python {SCRIPT_PATH} preflight "
        f"--manifest {workspace / 'manifest.json'} --index 0"
    )
    return 0


def _manifest_job(
    manifest_path: Path, index: int
) -> tuple[dict[str, Any], dict[str, Any]]:
    payload = load_manifest(manifest_path)
    if index < 0 or index >= len(payload["jobs"]):
        fail(
            f"Invalid job index {index}; valid indices are 0..{len(payload['jobs']) - 1}"
        )
    validate_manifest(payload, check_sources=True, source_indices={index})
    return payload, payload["jobs"][index]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    for command in ("scan", "prepare"):
        action = subparsers.add_parser(command)
        action.add_argument("--mode", choices=SUPPORTED_MODES, required=True)
        action.add_argument("--input", type=Path, required=True)
        action.add_argument("--config", type=Path, required=True)
        action.add_argument("--repository", type=Path, default=REPOSITORY_ROOT)
        if command == "prepare":
            action.add_argument("--workspace", type=Path, required=True)
            action.add_argument("--allow-dirty-repository", action="store_true")
    validate = subparsers.add_parser("validate")
    validate.add_argument("--manifest", type=Path, required=True)
    for command in ("preflight", "run", "resume"):
        action = subparsers.add_parser(command)
        action.add_argument("--manifest", type=Path, required=True)
        action.add_argument("--index", type=int, required=True)
    inspect = subparsers.add_parser("inspect")
    inspect.add_argument("--manifest", type=Path, required=True)
    submit = subparsers.add_parser("submit")
    submit.add_argument("--manifest", type=Path, required=True)
    submit.add_argument("--max-parallel", type=int, default=2)
    submit.add_argument("--print-only", action="store_true")
    submit.add_argument("--allow-duplicate", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        if args.command in {"scan", "prepare"}:
            input_root = absolute(args.input)
            config_path = absolute(args.config)
            repository = absolute(args.repository)
            if args.command == "scan":
                config = load_config(config_path, args.mode, repository)
                return scan(args.mode, input_root, config, repository)
            return prepare(
                args.mode,
                input_root,
                config_path,
                absolute(args.workspace),
                repository,
                allow_dirty_repository=args.allow_dirty_repository,
            )
        manifest_path = absolute(args.manifest)
        if args.command == "validate":
            summary = validate_manifest(
                load_manifest(manifest_path), check_sources=True
            )
            print("MANIFEST VALID: " + json.dumps(summary, sort_keys=True))
            return 0
        if args.command == "inspect":
            batch_runner.inspect(load_manifest(manifest_path))
            return 0
        if args.command == "submit":
            batch_runner.submit(
                manifest_path,
                args.max_parallel,
                print_only=args.print_only,
                allow_duplicate=args.allow_duplicate,
            )
            return 0
        payload, job = _manifest_job(manifest_path, args.index)
        if args.command == "preflight":
            batch_runner.preflight_job(payload, job)
        else:
            batch_runner.run_job(payload, job, resume=args.command == "resume")
        return 0
    except (
        ValueError,
        FileNotFoundError,
        KeyError,
        IndexError,
        json.JSONDecodeError,
        RuntimeError,
        subprocess.CalledProcessError,
    ) as exc:
        parser.exit(2, f"ERROR: {exc}\n")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
