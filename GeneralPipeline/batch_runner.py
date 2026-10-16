"""Execute prepared jobs and submit reproducible GeneralPipeline Slurm arrays."""

from __future__ import annotations

import importlib.util
import json
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path
from types import ModuleType
from typing import Any

from .manifest import atomic_json, fail, load_manifest, validate_manifest


def _module(name: str, path: Path) -> ModuleType:
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load Python module from {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _command(command: list[str], label: str) -> None:
    print(
        f"\n### {label}\n"
        + " ".join(repr(item) if " " in item else item for item in command),
        flush=True,
    )
    subprocess.run(command, check=True)


def _raw_preflight(
    job: dict[str, Any], config: dict[str, Any], compatibility: Path
) -> list[str]:
    from .channel_metadata import validate_tile_metadata

    sources = [Path(item["staged"]) for item in job["sources"]]
    physical = validate_tile_metadata(
        sources, config["channel_names"], config["ims_unit"]
    )
    if physical != job["physical_metadata"]:
        fail(f"Raw IMS metadata changed after preparation for {job['name']}")
    stitch = config["stitching"]
    grid = stitch["grid"]
    preprocessing = Path(job["output_dir"]) / "preprocessing"
    command = [
        sys.executable,
        str(compatibility / "Stitching/IMSTIFF_trim_zero_padding.py"),
        str(sources[0]),
        "--fiji-sif",
        config["fiji_sif"],
        "--output-root",
        str(preprocessing),
        "--channel-names",
        *config["channel_names"],
        "--grid",
        str(grid[0]),
        str(grid[1]),
        "--overlap",
        str(stitch.get("overlap", 10.0)),
        "--min-correlation",
        str(stitch.get("min_correlation", 0.55)),
        "--basic-sample-planes",
        str(stitch.get("basic_sample_planes", 12)),
        "--ims-unit",
        config["ims_unit"],
        "--intensity-downsampling",
        str(stitch.get("intensity_downsampling", 32)),
        "--intensity-max-inliers",
        str(stitch.get("intensity_max_inliers", 10000)),
        "--intensity-offset-only",
        str(stitch.get("intensity_offset_only", 0.5)),
        "--intensity-unmodified",
        str(stitch.get("intensity_unmodified", 0.5)),
        "--keep-previous",
        "--cleanup-stitching-intermediates",
    ]
    if not config["trim_zero_padding"]:
        command.append("--no-trim-zero-padding")
    return command


def preflight_job(payload: dict[str, Any], job: dict[str, Any]) -> None:
    config = json.loads(Path(job["config_path"]).read_text(encoding="utf-8"))
    if job["input_mode"] == "raw_ims_tiles":
        command = _raw_preflight(job, config, Path(payload["compat_repository"]))
        print(
            "Raw IMS metadata, calibration, physical channels, and biological mapping: OK"
        )
        print("Planned preprocessing command: " + " ".join(command))
        print(f"Shared downstream runner: {job['downstream_runner']}")
        return
    downstream = _module(
        "general_pipeline_stitched_preflight", Path(job["downstream_runner"])
    )
    parsed = downstream.read_config(Path(job["config_path"]))
    descriptions = downstream.validate_inputs(parsed)
    print("Stitched TIFF mapping and shared downstream configuration: OK")
    for name, description in descriptions.items():
        print(
            f"  {name}: shape={description['shape']}, axes={description['axes']}, "
            f"projection_required={description['projection_required']}"
        )
    print("IMS export, BaSiC, and BigStitcher: BYPASSED")
    print(f"Shared downstream runner: {job['downstream_runner']}")


def _run_raw(payload: dict[str, Any], job: dict[str, Any], *, resume: bool) -> None:
    import tifffile

    compatibility = Path(payload["compat_repository"])
    repository = Path(payload["repository"])
    config_path = Path(job["config_path"])
    config = json.loads(config_path.read_text(encoding="utf-8"))
    command = _raw_preflight(job, config, compatibility)
    _command(command, "IMS export, BaSiC correction, BigStitcher, and projections")

    ims = _module(
        "general_pipeline_ims_controller", compatibility / "run_ims_pipeline.py"
    )
    tiff_safety = _module(
        "general_pipeline_tiff_safety", repository / "Stitching/tiff_safety.py"
    )
    sample = config["stage_sample_name"]
    preprocessing = Path(job["output_dir"]) / "preprocessing"
    full_projections = preprocessing / sample / "04_max_projections"
    analysis_projections = preprocessing / sample / "05_analysis_projections"
    factor = config["analysis_downsampling"]["factor"]
    outputs: dict[str, Path] = {}
    downsampling_geometry: dict[str, Any] = {}
    for channel in config["channel_names"]:
        source = full_projections / f"{sample}_{channel}_stitched_max_projection.tif"
        if not source.is_file():
            raise FileNotFoundError(f"Missing expected BigStitcher output: {source}")
        if factor == 1:
            outputs[channel] = source
            continue
        image = tifffile.imread(source)
        averaged, details = ims.area_mean_downsample(image, factor)
        target = (
            analysis_projections
            / f"{sample}_{channel}_stitched_max_projection_downsample{factor}.tif"
        )
        target.parent.mkdir(parents=True, exist_ok=True)
        if not tiff_safety.is_complete_tiff(
            target, tuple(averaged.shape), averaged.dtype
        ):
            tiff_safety.atomic_tiff_write(
                target,
                averaged,
                bigtiff=True,
                compression="zlib",
                metadata={"axes": "YX"},
            )
        outputs[channel] = target
        downsampling_geometry[channel] = details

    pixel = job["physical_metadata"]["voxel_size_um"]
    downstream_config = dict(config["analysis"])
    downstream_config.update(
        {
            "nuclear_channel": config["nuclear_channel"],
            "channels": {name: str(path) for name, path in outputs.items()},
            "output_dir": str(Path(job["output_dir"]) / "analysis"),
            "pixel_size_um": {
                "x": float(pixel["x"]) * factor,
                "y": float(pixel["y"]) * factor,
            },
        }
    )
    generated = config_path.with_suffix(".stitched.json")
    if (
        generated.exists()
        and json.loads(generated.read_text(encoding="utf-8")) != downstream_config
    ):
        fail(f"Generated downstream configuration changed during resume: {generated}")
    atomic_json(generated, downstream_config)
    provenance = {
        "input_mode": "raw_ims_tiles",
        "source_acquisition": job["source_acquisition"],
        "sources": job["sources"],
        "channel_names": config["channel_names"],
        "nuclear_channel": config["nuclear_channel"],
        "physical_metadata": job["physical_metadata"],
        "analysis_downsampling": config["analysis_downsampling"],
        "downsampling_geometry": downsampling_geometry,
        "generated_stitched_config": str(generated),
        "downstream_runner": job["downstream_runner"],
    }
    atomic_json(
        Path(job["output_dir"]) / "general_pipeline_provenance.json", provenance
    )
    downstream_command = [
        sys.executable,
        job["downstream_runner"],
        "--config",
        str(generated),
    ]
    if resume:
        downstream_command.append("--resume")
    _command(downstream_command, "Shared stitched-image downstream analysis")


def _run_stitched(job: dict[str, Any], *, resume: bool) -> None:
    atomic_json(
        Path(job["output_dir"]) / "general_pipeline_provenance.json",
        {
            "input_mode": "stitched_tiffs",
            "source_acquisition": job["source_acquisition"],
            "sources": job["sources"],
            "channel_names": job["channel_names"],
            "nuclear_channel": job["nuclear_channel"],
            "pixel_size_um": job["pixel_size_um"],
            "generated_stitched_config": job["config_path"],
            "downstream_runner": job["downstream_runner"],
            "bypassed_stages": job["bypassed_stages"],
        },
    )
    command = [sys.executable, job["downstream_runner"], "--config", job["config_path"]]
    if resume:
        command.append("--resume")
    _command(command, "Shared stitched-image downstream analysis")


def run_job(payload: dict[str, Any], job: dict[str, Any], *, resume: bool) -> None:
    output = Path(job["output_dir"])
    if not resume and output.exists() and any(output.iterdir()):
        fail(
            f"Output already contains files for {job['name']}; use the resume command "
            "only when inputs and parameters are unchanged"
        )
    if job["input_mode"] == "raw_ims_tiles":
        _run_raw(payload, job, resume=resume)
    else:
        _run_stitched(job, resume=resume)


def inspect(payload: dict[str, Any]) -> None:
    print(f"Workspace: {payload['workspace']}  mode={payload['input_mode']}")
    for job in payload["jobs"]:
        analysis_manifest = Path(job["output_dir"]) / "analysis/run_manifest.json"
        state = "COMPLETE/EXISTS" if analysis_manifest.is_file() else "NOT COMPLETE"
        print(f"[{job['index']:03d}] {job['name']}: {state} ({analysis_manifest})")


def submit(
    manifest_path: Path,
    maximum_parallel: int,
    *,
    print_only: bool = False,
    allow_duplicate: bool = False,
    sbatch: str = "sbatch",
) -> str | None:
    manifest_path = manifest_path.expanduser().resolve()
    payload = load_manifest(manifest_path)
    validate_manifest(payload, check_sources=True)
    workspace = Path(payload["workspace"])
    if manifest_path != workspace / "manifest.json":
        fail("Manifest must be located at <workspace>/manifest.json")
    if maximum_parallel < 1:
        fail("--max-parallel must be at least 1")
    jobs = payload["jobs"]
    throttle = min(len(jobs), maximum_parallel)
    array = f"0-{len(jobs) - 1}%{throttle}"
    recorded = workspace / "slurm_submission.json"
    if recorded.exists() and not (print_only or allow_duplicate):
        previous = json.loads(recorded.read_text(encoding="utf-8"))
        fail(
            f"Workspace was already submitted as job {previous.get('slurm_job_id')}; "
            "inspect squeue/sacct before using --allow-duplicate"
        )
    command = [
        sbatch,
        "--parsable",
        f"--array={array}",
        f"--chdir={workspace}",
        payload["slurm_runner"],
        str(manifest_path),
        payload["pipeline_runner"],
    ]
    print(f"Selected {len(jobs)} acquisitions; maximum simultaneous tasks: {throttle}")
    print("Submit command: " + " ".join(command))
    if print_only:
        print("PRINT ONLY: no Slurm job submitted")
        return None
    if not shutil.which(sbatch):
        fail(f"Could not find Slurm command {sbatch!r}")
    result = subprocess.run(
        command, cwd=workspace, text=True, capture_output=True, check=True
    )
    value = result.stdout.strip()
    match = re.fullmatch(r"(\d+)(?:;\S+)?", value)
    if match is None:
        fail(
            f"Slurm returned unexpected output {value!r}; inspect the queue before retrying"
        )
    job_id = match.group(1)
    atomic_json(
        recorded,
        {
            "slurm_job_id": job_id,
            "array": array,
            "manifest": str(manifest_path),
            "workspace": str(workspace),
            "input_mode": payload["input_mode"],
            "command": command,
            "submitted_unix_time": int(time.time()),
        },
    )
    print(f"Submitted Slurm array {job_id}")
    return job_id
