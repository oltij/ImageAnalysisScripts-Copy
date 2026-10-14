"""Crash-safe manifests, provenance, and source/workspace safety checks."""

from __future__ import annotations

import hashlib
import json
import os
import subprocess
from pathlib import Path
from typing import Any


SCHEMA_VERSION = 1
MANIFEST_KIND = "general_microscopy_pipeline"
TURBO_IMS_ROOT = Path("/nfs/turbo/umms-parent/andor_microscope_data_dump")
SCRATCH_ROOT = Path("/scratch")

SCIENTIFIC_FILES = (
    "run_ims_pipeline.py",
    "run_stitched_pipeline.py",
    "Stitching/IMSTIFF_trim_zero_padding.py",
    "Stitching/tiff_safety.py",
    "Alignment/CASTalign_two_channel_registration.py",
    "CellProfiler/cellprofilerdriver.py",
    "CellProfiler/MGEOPVFinal.cppipe",
    "ExtractingROIs/extract_rois.py",
    "ROIFiltering/Intensity/intensityfilter.py",
    "ROIFiltering/FinalFilter/filter.py",
    "ROIFiltering/Area/areafilter.py",
    "ROIFiltering/Circularity/circularityfilter.py",
    "ROIFiltering/Eccentricity/eccentricityfilter.py",
    "ROIFiltering/Solidity/solidityfilter.py",
    "QC/organoid_mask_qc.py",
    "Colocalization/colocalizationdapiscript1.py",
    "Colocalization/colocalizationdapiscript2.py",
)


def fail(message: str) -> None:
    raise ValueError(message)


def absolute(value: str | Path) -> Path:
    return Path(value).expanduser().resolve()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def file_stats(path: Path) -> dict[str, Any]:
    stat = path.stat()
    return {"path": str(path), "size": stat.st_size, "mtime_ns": stat.st_mtime_ns}


def atomic_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    try:
        temporary.write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def repository_git_state(repository: Path) -> dict[str, Any]:
    try:
        commit = subprocess.run(
            ["git", "-C", str(repository), "rev-parse", "HEAD"],
            check=True,
            text=True,
            capture_output=True,
        ).stdout.strip()
        status = subprocess.run(
            ["git", "-C", str(repository), "status", "--porcelain"],
            check=True,
            text=True,
            capture_output=True,
        ).stdout.splitlines()
    except (FileNotFoundError, subprocess.CalledProcessError) as exc:
        return {"commit": None, "dirty": None, "status": [], "error": str(exc)}
    return {"commit": commit, "dirty": bool(status), "status": status}


def require_reproducible_checkout(state: dict[str, Any], allow_dirty: bool) -> None:
    if state.get("commit") is None:
        fail(
            "Preparation requires a Git checkout so the repository revision can be recorded"
        )
    if state.get("dirty") and not allow_dirty:
        fail(
            "Preparation refuses a dirty repository; commit/stash changes or use --allow-dirty-repository for development"
        )


def scientific_hashes(repository: Path) -> dict[str, str]:
    result: dict[str, str] = {}
    for relative in SCIENTIFIC_FILES:
        path = repository / relative
        if not path.is_file():
            raise FileNotFoundError(
                f"Required scientific pipeline file is missing: {path}"
            )
        result[relative] = sha256_file(path)
    return result


def validate_storage_locations(source_root: Path, workspace: Path) -> None:
    source = source_root.resolve(strict=False)
    destination = workspace.resolve(strict=False)
    turbo = TURBO_IMS_ROOT.resolve(strict=False)
    scratch = SCRATCH_ROOT.resolve(strict=False)
    if (
        destination == source
        or source in destination.parents
        or destination in source.parents
    ):
        fail("Workspace and source acquisition trees must be separate")
    protected = source == turbo or turbo in source.parents
    on_scratch = destination == scratch or scratch in destination.parents
    if protected and not on_scratch:
        fail(
            f"Inputs below {TURBO_IMS_ROOT} require a workspace below {SCRATCH_ROOT}; got {workspace}"
        )


def fresh_workspace(path: Path) -> None:
    if path.exists() and (not path.is_dir() or any(path.iterdir())):
        fail(f"Workspace exists and is not empty: {path}; choose a new directory")
    path.mkdir(parents=True, exist_ok=True)


def create_controller_compatibility(repository: Path, workspace: Path) -> Path:
    """Preserve staged F00-F03 symlink names in the existing IMS controllers."""
    destination = workspace / "controller_compat"
    (destination / "Stitching").mkdir(parents=True, exist_ok=True)
    replacements = (
        (
            "run_ims_pipeline.py",
            "            path.resolve()\n            for path in directory.iterdir()",
            "            path.absolute()\n            for path in directory.iterdir()",
        ),
        (
            "Stitching/IMSTIFF_trim_zero_padding.py",
            "    input_path = args.path.expanduser().resolve()\n",
            "    input_path = args.path.expanduser().absolute()\n",
        ),
    )
    for relative, old, new in replacements:
        source = repository / relative
        text = source.read_text(encoding="utf-8")
        if text.count(old) != 1:
            fail(
                f"Cannot generate staged-symlink compatibility controller from {source}"
            )
        target = destination / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(text.replace(old, new, 1), encoding="utf-8")
    for relative in ("run_stitched_pipeline.py", "Stitching/tiff_safety.py"):
        target = destination / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        target.symlink_to((repository / relative).resolve())
    return destination


def load_manifest(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if (
        not isinstance(payload, dict)
        or payload.get("schema") != SCHEMA_VERSION
        or payload.get("kind") != MANIFEST_KIND
        or payload.get("input_mode") not in {"raw_ims_tiles", "stitched_tiffs"}
        or not isinstance(payload.get("jobs"), list)
    ):
        fail(f"Not a supported GeneralPipeline manifest: {path}")
    return payload


def _within(path: Path, root: Path) -> bool:
    resolved = path.resolve(strict=False)
    return resolved == root or root in resolved.parents


def validate_manifest(
    payload: dict[str, Any],
    *,
    check_sources: bool = True,
    source_indices: set[int] | None = None,
) -> dict[str, Any]:
    jobs = payload["jobs"]
    if not jobs or [job.get("index") for job in jobs] != list(range(len(jobs))):
        fail("Manifest jobs must use deterministic sequential indices 0..N-1")
    workspace = Path(payload["workspace"]).resolve(strict=False)
    source_root = Path(payload["input_root"]).resolve(strict=False)
    validate_storage_locations(source_root, workspace)
    for job in jobs:
        for key in ("staging_dir", "output_dir", "config_path"):
            if not _within(Path(job[key]), workspace):
                fail(f"Job {job['index']} {key} is outside its workspace: {job[key]}")
        generated_config = Path(job["config_path"])
        if not generated_config.is_file() or sha256_file(generated_config) != job.get(
            "config_sha256"
        ):
            fail(
                f"Generated job configuration changed after preparation: {generated_config}"
            )
        if job.get("downstream_runner") != str(
            Path(payload["repository"]) / "run_stitched_pipeline.py"
        ):
            fail(f"Job {job['index']} does not use the shared downstream runner")
        if check_sources and (source_indices is None or job["index"] in source_indices):
            for source in job["sources"]:
                original = Path(source["source"])
                staged = Path(source["staged"])
                validate_storage_locations(original, workspace)
                if not staged.is_symlink() or staged.resolve() != original:
                    fail(f"Staged input symlink is missing or changed: {staged}")
                if file_stats(original) != source["source_stats"]:
                    fail(f"Source input changed after preparation: {original}")

    tracked = {
        Path(payload["config_path"]): payload["config_sha256"],
        Path(payload["pipeline_runner"]): payload["pipeline_runner_sha256"],
        Path(payload["slurm_runner"]): payload["slurm_runner_sha256"],
        Path(payload["cellprofiler_pipeline"]): payload["cellprofiler_pipeline_sha256"],
    }
    for path, expected in tracked.items():
        if not path.is_file() or sha256_file(path) != expected:
            fail(f"Prepared file is missing or changed: {path}")
    compatibility = payload.get("compat_repository")
    if compatibility is not None:
        compatibility_root = Path(compatibility)
        if not _within(compatibility_root, workspace):
            fail(
                f"Compatibility repository is outside its workspace: {compatibility_root}"
            )
        for relative, expected in payload.get("compatibility_sha256", {}).items():
            path = compatibility_root / relative
            if not path.is_file() or sha256_file(path) != expected:
                fail(
                    f"Generated compatibility controller changed after preparation: {path}"
                )
    repository = Path(payload["repository"])
    for relative, expected in payload["scientific_sha256"].items():
        path = repository / relative
        if not path.is_file() or sha256_file(path) != expected:
            fail(f"Scientific pipeline file changed after preparation: {path}")
    current = repository_git_state(repository)
    recorded = payload["repository_git"]
    if current.get("commit") != recorded.get("commit"):
        fail(
            "Repository commit changed after preparation; create a new workspace manifest"
        )
    return {
        "input_mode": payload["input_mode"],
        "jobs": len(jobs),
        "sources": sum(len(job["sources"]) for job in jobs),
        "sources_checked": (
            sum(len(job["sources"]) for job in jobs)
            if check_sources and source_indices is None
            else sum(
                len(job["sources"])
                for job in jobs
                if check_sources and job["index"] in (source_indices or set())
            )
        ),
        "git_commit": recorded.get("commit"),
    }


def copy_or_link(source: Path, destination: Path) -> None:
    """Create a source-preserving symlink and fail instead of replacing data."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() or destination.is_symlink():
        fail(f"Refusing to replace an existing staged input: {destination}")
    destination.symlink_to(source.resolve())
