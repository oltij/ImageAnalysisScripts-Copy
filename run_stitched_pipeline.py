#!/usr/bin/env python3
"""Run the ImageAnalysisScripts stages from channel-specific stitched TIFFs.

Initial CellProfiler segmentation uses the configured pipeline unchanged. For final
marker segmentation, a derived copy loads the aligned Hoechst organoid mask as a
second input while preserving marker-specific cell detection and other settings.
See STITCHED_WORKFLOW.md before using this on experimental data.
"""

from __future__ import annotations

import argparse
import ast
import hashlib
import importlib.util
import io
import itertools
import json
import math
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import tifffile

ROOT = Path(__file__).resolve().parent

INTENSITY_METRICS = (
    "Mean intensity", "Median intensity", "Maximum intensity",
    "Integrated intensity", "Intensity SD",
)
SHAPE_METRICS = (
    "Area_um2", "Perimeter_um", "Solidity", "Eccentricity",
    "Circularity", "ConvexHullArea_um2", "MajorAxisLength_um",
    "MinorAxisLength_um", "AspectRatio", "EquivalentDiameter_um", "Extent",
)
SHAPE_REPORTS = {
    "area": "ROIFiltering/Area/areafilter.py",
    "circularity": "ROIFiltering/Circularity/circularityfilter.py",
    "eccentricity": "ROIFiltering/Eccentricity/eccentricityfilter.py",
    "solidity": "ROIFiltering/Solidity/solidityfilter.py",
}

ANALYSIS_FILES = (
    "Alignment/CASTalign_two_channel_registration.py",
    "CellProfiler/cellprofilerdriver.py",
    "ExtractingROIs/extract_rois.py",
    "ROIFiltering/Intensity/intensityfilter.py",
    "ROIFiltering/FinalFilter/filter.py",
    *SHAPE_REPORTS.values(),
    "QC/organoid_mask_qc.py",
    "Colocalization/colocalizationdapiscript1.py",
    "Colocalization/colocalizationdapiscript2.py",
)

EXPECTED_CP_PIPELINE_SETTINGS = (
    "Name the output objects:FilterObjects",
    "Name the output objects:FilterObjects2",
    "Enter single file name:OrganoidMask",
    "Enter single file name:CellMask",
    "Filename prefix:MyExpt_",
)

MAIN_RUNTIME_MODULES = (
    "h5py", "matplotlib", "numpy", "pandas", "PIL", "scipy", "skimage", "tifffile",
)

SHARED_MASK_IMAGE_NAME = "HoechstOrganoidMask"
SHARED_FLUORESCENCE_TOKEN = "__analysis_fluorescence__"
SHARED_MASK_TOKEN = "__hoechst_organoid_mask__"


def path_value(value: object) -> str:
    """A Python expression for a pathlib.Path (or literal None)."""
    return "None" if value is None else f"Path({str(value)!r})"


def patched_source(path: Path, replacements: dict[str, str]) -> bytes:
    """Replace only top-level assignment RHS values, preserving all other bytes.

    ast coordinates are UTF-8 byte offsets, so patch bytes, not characters.
    Values are caller-supplied Python expressions, never raw shell commands.
    """
    src = path.read_bytes()
    tree = ast.parse(src, filename=str(path))
    line_offsets = [0]
    for line in src.splitlines(keepends=True):
        line_offsets.append(line_offsets[-1] + len(line))

    edits: list[tuple[int, int, bytes]] = []
    found: set[str] = set()
    for statement in tree.body:
        if isinstance(statement, ast.Assign):
            names = [t.id for t in statement.targets if isinstance(t, ast.Name)]
            value = statement.value
        elif isinstance(statement, ast.AnnAssign) and isinstance(statement.target, ast.Name):
            names = [statement.target.id]
            value = statement.value
        else:
            continue
        chosen = [name for name in names if name in replacements]
        if not chosen or value is None:
            continue
        if len(chosen) != 1 or len(names) != 1:
            raise ValueError(f"Ambiguous assignment in {path}: {chosen}")
        name = chosen[0]
        if name in found:
            raise ValueError(f"Repeated top-level assignment to {name} in {path}")
        expression = replacements[name]
        ast.parse(expression, mode="eval")
        start = line_offsets[value.lineno - 1] + value.col_offset
        end = line_offsets[value.end_lineno - 1] + value.end_col_offset
        edits.append((start, end, expression.encode("utf-8")))
        found.add(name)
    missing = set(replacements) - found
    if missing:
        raise ValueError(f"Missing configurable assignments in {path}: {sorted(missing)}")
    for start, end, replacement in sorted(edits, reverse=True):
        src = src[:start] + replacement + src[end:]
    ast.parse(src, filename=str(path))
    return src


def resolve_path(value: str, base: Path) -> Path:
    p = Path(value).expanduser()
    return (p if p.is_absolute() else base / p).resolve()


def validate_cellprofiler_pipeline(path: Path) -> str:
    """Require a real text .cppipe with the outputs used by this workflow."""
    if not path.is_file():
        raise FileNotFoundError(f"CellProfiler pipeline missing: {path}")
    prefix = path.read_bytes()[:8]
    if prefix == b"\x89HDF\r\n\x1a\n":
        raise ValueError(
            f"{path} is an HDF5 CellProfiler project, not a text .cppipe. "
            "Export its pipeline in CellProfiler or use the text MGEOPVFinal.cppipe supplied here."
        )
    try:
        pipeline_text = path.read_text(encoding="utf-8")
    except UnicodeDecodeError as exc:
        raise ValueError(f"CellProfiler pipeline is not UTF-8 text: {path}") from exc
    if not pipeline_text.startswith("CellProfiler Pipeline:"):
        raise ValueError(f"Not a recognizable CellProfiler text pipeline: {path}")
    missing = [token for token in EXPECTED_CP_PIPELINE_SETTINGS if token not in pipeline_text]
    if missing:
        raise ValueError(
            "CellProfiler pipeline is incompatible with the automated handoff; "
            f"missing settings for expected outputs: {missing}"
        )
    return pipeline_text


def cellprofiler_module(text: str, name: str) -> tuple[re.Match[str], str]:
    match = re.search(
        rf"(?ms)^{re.escape(name)}:\[.*?(?=^[A-Za-z][A-Za-z0-9]*:\[|\Z)",
        text,
    )
    if match is None:
        raise ValueError(f"CellProfiler pipeline is missing required module {name}")
    return match, match.group(0).rstrip()


def write_shared_mask_pipeline(source: Path, destination: Path) -> Path:
    """Derive a final-marker pipeline that uses an aligned Hoechst tissue mask."""
    text = validate_cellprofiler_pipeline(source)
    names_match, names_block = cellprofiler_module(text, "NamesAndTypes")
    required_names = (
        "Assign a name to:All images",
        "Assignments count:1",
        "Name to assign these images:DNA",
    )
    missing = [setting for setting in required_names if setting not in names_block]
    if missing:
        raise ValueError(
            "Cannot derive the shared-mask pipeline from this NamesAndTypes "
            f"configuration; missing {missing}"
        )
    names_header = names_block.splitlines()[0]
    shared_names = (
        f"{names_header}\n"
        "    Assign a name to:Images matching rules\n"
        "    Select the image type:Grayscale image\n"
        "    Name to assign these images:DNA\n"
        "    Match metadata:[]\n"
        "    Image set matching method:Order\n"
        "    Set intensity range from:Image metadata\n"
        "    Assignments count:2\n"
        "    Single images count:0\n"
        "    Maximum intensity:255.0\n"
        "    Process as 3D?:No\n"
        "    Relative pixel spacing in X:1.0\n"
        "    Relative pixel spacing in Y:1.0\n"
        "    Relative pixel spacing in Z:1.0\n"
        f'    Select the rule criteria:and (file does contain "{SHARED_FLUORESCENCE_TOKEN}")\n'
        "    Name to assign these images:DNA\n"
        "    Name to assign these objects:Cell\n"
        "    Select the image type:Grayscale image\n"
        "    Set intensity range from:Image metadata\n"
        "    Maximum intensity:255.0\n"
        f'    Select the rule criteria:and (file does contain "{SHARED_MASK_TOKEN}")\n'
        f"    Name to assign these images:{SHARED_MASK_IMAGE_NAME}\n"
        "    Name to assign these objects:Cell\n"
        "    Select the image type:Grayscale image\n"
        "    Set intensity range from:Image metadata\n"
        "    Maximum intensity:255.0"
    )
    text = text[: names_match.start()] + shared_names + "\n\n" + text[names_match.end() :]

    convert_match, convert_block = cellprofiler_module(text, "ConvertImageToObjects")
    required_convert = (
        "Select the input image:Threshold",
        "Name the output object:ConvertImageToObjects",
    )
    missing = [
        setting for setting in required_convert if setting not in convert_block
    ]
    if missing:
        raise ValueError(
            "Cannot derive the shared-mask pipeline from ConvertImageToObjects; "
            f"missing {missing}"
        )
    shared_convert = convert_block.replace(
        "Select the input image:Threshold",
        f"Select the input image:{SHARED_MASK_IMAGE_NAME}",
        1,
    )
    text = (
        text[: convert_match.start()]
        + shared_convert
        + "\n\n"
        + text[convert_match.end() :]
    )
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(".cppipe.tmp")
    temporary.write_text(text.rstrip() + "\n", encoding="utf-8")
    temporary.replace(destination)
    validate_cellprofiler_pipeline(destination)
    return destination


def read_binary_mask(path: Path, label: str) -> np.ndarray:
    array = np.squeeze(tifffile.imread(path))
    if array.ndim != 2:
        raise ValueError(f"{label} must be a 2-D TIFF, got {array.shape}: {path}")
    return array > 0


def validate_shared_mask_canvas(image: Path, organoid_mask: Path) -> None:
    image_info = inspect_tiff(image)
    if image_info["projection_required"]:
        raise ValueError(f"Shared-mask fluorescence must already be 2-D: {image}")
    mask = read_binary_mask(organoid_mask, "Hoechst organoid mask")
    if not np.any(mask):
        raise ValueError(
            f"Hoechst organoid mask contains no foreground pixels: {organoid_mask}"
        )
    if tuple(image_info["yx_shape"]) != mask.shape:
        raise ValueError(
            "Aligned marker and Hoechst organoid mask have different XY shapes: "
            f"{image_info['yx_shape']} versus {mask.shape}"
        )


def verify_shared_mask_outputs(
    reference_path: Path, organoid_path: Path, cell_path: Path
) -> None:
    reference = read_binary_mask(reference_path, "Hoechst organoid mask")
    organoid = read_binary_mask(organoid_path, "marker organoid mask")
    cells = read_binary_mask(cell_path, "marker cell mask")
    if not np.any(reference):
        raise RuntimeError("Hoechst reference organoid mask is empty")
    if not np.any(organoid):
        raise RuntimeError(
            "Marker organoid mask is empty; shared-mask segmentation did not succeed"
        )
    if reference.shape != organoid.shape or reference.shape != cells.shape:
        raise RuntimeError(
            "Shared-mask CellProfiler outputs do not match the Hoechst mask canvas"
        )
    disagreement = int(np.count_nonzero(organoid ^ reference))
    outside = int(np.count_nonzero(cells & ~reference))
    if disagreement or outside:
        raise RuntimeError(
            "Shared Hoechst-mask verification failed: "
            f"{disagreement} organoid-mask pixels disagree and "
            f"{outside} segmented cell pixels are outside the Hoechst mask"
        )


def inspect_tiff(path: Path) -> dict[str, object]:
    """Return the supported single-channel TIFF layout without loading all pixels."""
    import tifffile

    try:
        with tifffile.TiffFile(path) as tf:
            series = tf.series[0]
            shape = tuple(int(v) for v in series.shape)
            axes = str(series.axes)
            dtype = str(series.dtype)
    except Exception as exc:
        raise ValueError(f"Could not read TIFF metadata from {path}: {exc}") from exc

    if len(shape) == 2 and axes == "YX":
        projection_required = False
    elif len(shape) == 3 and axes in ("ZYX", "QYX", "IYX"):
        projection_required = True
    else:
        raise ValueError(
            f"Expected a single-channel YX image or ZYX page stack, got "
            f"shape={shape}, axes={axes!r}: {path}"
        )
    return {
        "path": str(path),
        "shape": shape,
        "axes": axes,
        "dtype": dtype,
        "yx_shape": shape[-2:],
        "projection_required": projection_required,
    }


def validate_inputs(cfg: dict) -> dict[str, dict[str, object]]:
    """Fail before creating outputs if stitched TIFFs are missing or incompatible."""
    descriptions: dict[str, dict[str, object]] = {}
    for name, source in cfg["_channels"].items():
        if not source.is_file():
            raise FileNotFoundError(f"Stitched input for {name} does not exist: {source}")
        if source.stat().st_size == 0:
            raise ValueError(f"Stitched input for {name} is empty: {source}")
        descriptions[name] = inspect_tiff(source)
    dimensions = {tuple(info["yx_shape"]) for info in descriptions.values()}
    if len(dimensions) != 1:
        detail = ", ".join(
            f"{name}={tuple(info['yx_shape'])}" for name, info in descriptions.items()
        )
        raise ValueError(
            "All channels must share one XY pixel canvas before registration/colocalization; "
            + detail
        )
    validate_cellprofiler_pipeline(cfg["_pipeline"])
    return descriptions


def validate_runtime(cfg: dict) -> None:
    """Check both Python environments before a long run creates partial results."""
    modules = list(MAIN_RUNTIME_MODULES)
    if cfg.get("alignment", {}).get("enabled", True):
        modules.append("castalign")
    missing = [name for name in modules if importlib.util.find_spec(name) is None]
    if missing:
        raise RuntimeError(
            "The active analysis environment is missing Python package(s): "
            + ", ".join(missing)
            + ". Activate the ims-mosaic environment described in STITCHED_WORKFLOW.md."
        )

    cp = cfg.get("cellprofiler", {})
    conda = str(cp.get("conda_executable", "conda"))
    if shutil.which(conda) is None:
        raise RuntimeError(f"Conda executable not found: {conda}")
    environment = str(cp.get("conda_env", "cellprofiler-native"))
    check = subprocess.run(
        [conda, "run", "-n", environment, "python", "-c", "import cellprofiler"],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    if check.returncode:
        raise RuntimeError(
            f"CellProfiler environment {environment!r} is unavailable or cannot import "
            f"CellProfiler. Conda said:\n{check.stdout.strip()}"
        )


def analysis_code_sha256(pipeline: Path) -> str:
    """Fingerprint the controller, pipeline, and scientific scripts used by --resume."""
    digest = hashlib.sha256()
    for path in [Path(__file__).resolve(), pipeline, *(ROOT / p for p in ANALYSIS_FILES)]:
        if not path.is_file():
            raise FileNotFoundError(f"Required analysis file missing: {path}")
        digest.update(str(path.relative_to(ROOT) if path.is_relative_to(ROOT) else path).encode())
        digest.update(b"\0")
        digest.update(path.read_bytes())
        digest.update(b"\0")
    return digest.hexdigest()


def filters_for_pass(settings: dict, pass_name: str) -> dict:
    if pass_name not in ("intensity", "shape"):
        raise ValueError(pass_name)
    allowed = INTENSITY_METRICS if pass_name == "intensity" else SHAPE_METRICS
    defaults = {
        name: {"percentile": 0, "keep": "below" if name in ("Eccentricity", "AspectRatio") else "above"}
        for name in (*INTENSITY_METRICS, *SHAPE_METRICS)
    }
    user_settings = settings.get(pass_name, {})
    if not isinstance(user_settings, dict):
        raise ValueError(f"filters.{pass_name} must be an object")
    unknown = set(user_settings) - set(allowed)
    if unknown:
        raise ValueError(f"{pass_name} pass includes wrong/unknown metrics: {sorted(unknown)}")
    for name, rule in user_settings.items():
        if not isinstance(rule, dict) or set(rule) != {"percentile", "keep"}:
            raise ValueError(f"{name} must define percentile and keep")
        pct = rule["percentile"]
        if isinstance(pct, bool) or not isinstance(pct, (int, float)) or not 0 <= pct <= 100:
            raise ValueError(f"Invalid percentile for {name}: {pct}")
        if rule["keep"] not in ("above", "below"):
            raise ValueError(f"Invalid keep direction for {name}")
        defaults[name] = rule
    return defaults


def read_config(config_file: Path) -> dict:
    raw = json.loads(config_file.read_text(encoding="utf-8"))
    if not isinstance(raw, dict):
        raise ValueError("The top level of the configuration must be a JSON object")
    nuclear = raw["nuclear_channel"]
    raw_channels = raw["channels"]
    if not isinstance(nuclear, str) or not isinstance(raw_channels, dict) or not raw_channels or nuclear not in raw_channels:
        raise ValueError("channels must contain the named nuclear_channel")
    for name, value in raw_channels.items():
        if not re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*", name):
            raise ValueError(f"Use simple ASCII alphanumeric/underscore channel names: {name!r}")
        if not isinstance(value, str) or not value.strip():
            raise ValueError(f"Channel {name!r} must have a nonempty TIFF path")
    channels = {name: resolve_path(path, config_file.parent) for name, path in raw_channels.items()}
    marker_names = [name for name in channels if name != nuclear]
    pairs = raw.get("pairs", list(itertools.combinations(marker_names, 2)))
    if not isinstance(pairs, list):
        raise ValueError("pairs must be a JSON list")
    normalized_pairs: list[tuple[str, str]] = []
    for pair in pairs:
        if not isinstance(pair, (list, tuple)) or len(pair) != 2:
            raise ValueError(f"Each pair must contain two channel names: {pair!r}")
        if pair[0] == pair[1] or any(x not in marker_names for x in pair):
            raise ValueError(f"Pair must contain two distinct, non-nuclear channel names: {pair}")
        normalized_pairs.append((pair[0], pair[1]))
    if len(set(normalized_pairs)) != len(normalized_pairs):
        raise ValueError("pairs contains a duplicate directional pair")
    thresholds = raw.get("thresholds", {})
    if not isinstance(thresholds, dict):
        raise ValueError("thresholds must be a JSON object")
    if marker_names and "marker_nuclear" not in thresholds:
        raise ValueError("Explicit thresholds.marker_nuclear is required")
    if normalized_pairs and "marker_marker" not in thresholds:
        raise ValueError("Explicit thresholds.marker_marker is required")
    unknown_thresholds = set(thresholds) - {"marker_nuclear", "marker_marker"}
    if unknown_thresholds:
        raise ValueError(f"Unknown threshold settings: {sorted(unknown_thresholds)}")
    for label, v in thresholds.items():
        if isinstance(v, bool) or not isinstance(v, (int, float)) or not 0 <= v <= 1:
            raise ValueError(f"{label} threshold must lie in [0, 1]")
    pixels = raw["pixel_size_um"]
    if not isinstance(pixels, dict) or set(pixels) != {"x", "y"}:
        raise ValueError("pixel_size_um must contain exactly x and y")
    try:
        pixel_x = float(pixels["x"])
        pixel_y = float(pixels["y"])
    except (TypeError, ValueError) as exc:
        raise ValueError("Pixel sizes must be finite positive numbers") from exc
    if (
        isinstance(pixels["x"], bool)
        or isinstance(pixels["y"], bool)
        or not math.isfinite(pixel_x)
        or not math.isfinite(pixel_y)
        or pixel_x <= 0
        or pixel_y <= 0
    ):
        raise ValueError("Pixel sizes must be positive")
    filters = raw.get("filters", {})
    if not isinstance(filters, dict):
        raise ValueError("filters must be a JSON object")
    unknown_filter_groups = set(filters) - {"intensity", "shape"}
    if unknown_filter_groups:
        raise ValueError(f"Unknown filter groups: {sorted(unknown_filter_groups)}")
    filters_for_pass(filters, "intensity")
    filters_for_pass(filters, "shape")
    cp = raw.get("cellprofiler", {})
    if not isinstance(cp, dict):
        raise ValueError("cellprofiler must be a JSON object")
    unknown_cp = set(cp) - {"pipeline", "conda_env", "conda_executable"}
    if unknown_cp:
        raise ValueError(f"Unknown cellprofiler settings: {sorted(unknown_cp)}")
    pipeline_value = cp.get("pipeline", "CellProfiler/MGEOPVFinal.cppipe")
    if not isinstance(pipeline_value, str) or not pipeline_value.strip():
        raise ValueError("cellprofiler.pipeline must be a nonempty path")
    for setting in ("conda_env", "conda_executable"):
        if setting in cp and (not isinstance(cp[setting], str) or not cp[setting].strip()):
            raise ValueError(f"cellprofiler.{setting} must be a nonempty string")
    pipeline = resolve_path(pipeline_value, ROOT)
    validate_cellprofiler_pipeline(pipeline)
    alignment = raw.get("alignment", {})
    if not isinstance(alignment, dict) or not isinstance(alignment.get("enabled", True), bool):
        raise ValueError("alignment must be an object with a boolean enabled setting")
    extra_args = alignment.get("extra_args", [])
    if not isinstance(extra_args, list) or not all(isinstance(v, str) for v in extra_args):
        raise ValueError("alignment.extra_args must be a list of command-line strings")
    protected_args = {
        "--fixed-image", "--moving-image", "--fixed-csv", "--moving-csv", "--output",
        "--pixel-size-x", "--pixel-size-y", "--fixed-x-col", "--fixed-y-col",
        "--moving-x-col", "--moving-y-col",
    }
    if any(
        value in protected_args or any(value.startswith(flag + "=") for flag in protected_args)
        for value in extra_args
    ):
        raise ValueError("alignment.extra_args cannot override runner-managed file/calibration arguments")
    reports = raw.get("reports", {})
    allowed_reports = {"organoid_mask", "intensity", "shape", "colocalization"}
    if not isinstance(reports, dict) or set(reports) - allowed_reports:
        raise ValueError(
            "reports may contain only organoid_mask, intensity, shape, "
            "and colocalization"
        )
    if not all(isinstance(value, bool) for value in reports.values()):
        raise ValueError("Every reports setting must be true or false")
    output_value = raw.get("output_dir")
    if not isinstance(output_value, str) or not output_value.strip():
        raise ValueError("output_dir must be a nonempty path")
    result = dict(raw)
    result["_channels"] = channels
    result["_markers"] = marker_names
    result["_pairs"] = normalized_pairs
    result["_nuclear"] = nuclear
    result["_pipeline"] = pipeline
    result["_output"] = resolve_path(output_value, config_file.parent)
    if result["_output"].exists() and not result["_output"].is_dir():
        raise ValueError(f"output_dir exists but is not a directory: {result['_output']}")
    return result


class Runner:
    def __init__(self, cfg: dict, resume: bool):
        self.cfg = cfg
        self.out: Path = cfg["_output"]
        self.out.mkdir(parents=True, exist_ok=True)
        self.manifest = self.out / "run_manifest.json"
        self.log = self.out / "run.log"
        # Include input file metadata so --resume cannot silently reuse changed inputs.
        file_info = {
            name: {"path": str(path), "size": path.stat().st_size, "mtime_ns": path.stat().st_mtime_ns}
            for name, path in cfg["_channels"].items()
        }
        digest_data = {k: v for k, v in cfg.items() if not k.startswith("_")}
        self.code_digest = analysis_code_sha256(cfg["_pipeline"])
        fingerprint = {
            "config": digest_data,
            "inputs": file_info,
            "analysis_code_sha256": self.code_digest,
        }
        self.digest = hashlib.sha256(json.dumps(fingerprint, sort_keys=True).encode()).hexdigest()
        if self.manifest.exists():
            old = json.loads(self.manifest.read_text())
            if not resume:
                raise RuntimeError(f"Output already has a manifest; choose a new output_dir or --resume: {self.out}")
            if old.get("run_fingerprint_sha256") != self.digest:
                raise RuntimeError(
                    "--resume refused: configuration, input metadata, pipeline, controller, "
                    "or analysis scripts changed. Use a fresh output_dir."
                )
            self.completed = old["completed"]
            self.started_at = old["started_at"]
        else:
            if resume and any(self.out.iterdir()):
                raise RuntimeError("Cannot resume an untracked output directory")
            if any(self.out.iterdir()):
                raise RuntimeError(f"Refusing to mix with existing untracked files: {self.out}")
            self.completed = {}
            self.started_at = datetime.now(timezone.utc).isoformat()
            self.save()
        self.shared_mask_pipeline = write_shared_mask_pipeline(
            cfg["_pipeline"],
            self.out
            / "generated_cellprofiler"
            / "final_marker_with_hoechst_mask.cppipe",
        )

    def save(self) -> None:
        payload = {
            "manifest_version": 1,
            "started_at": self.started_at,
            "updated_at": datetime.now(timezone.utc).isoformat(),
            "run_fingerprint_sha256": self.digest,
            "analysis_code_sha256": self.code_digest,
            "pipeline": str(self.cfg["_pipeline"]),
            "completed": self.completed,
        }
        temp = self.manifest.with_suffix(".json.tmp")
        temp.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        temp.replace(self.manifest)

    def execute(self, label: str, cmd: list[str], expected: list[Path]) -> None:
        if label in self.completed and all(p.is_file() and p.stat().st_size for p in expected):
            print(f"[resume] {label}", flush=True)
            return
        print(f"\n{'=' * 70}\nRUNNING {label}\n{'=' * 70}", flush=True)
        command_text = shlex.join(cmd)
        print(command_text, flush=True)
        with self.log.open("a", encoding="utf-8") as log_file:
            log_file.write(
                f"\n[{datetime.now(timezone.utc).isoformat()}] {label}\n{command_text}\n"
            )
            log_file.flush()
            process = subprocess.Popen(
                cmd,
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
                log_file.write(line)
            returncode = process.wait()
        if returncode:
            raise subprocess.CalledProcessError(returncode, cmd)
        missing = [str(p) for p in expected if not p.is_file() or not p.stat().st_size]
        if missing:
            raise RuntimeError(f"{label} finished without required output(s): {missing}")
        self.completed[label] = [str(p) for p in expected]
        self.save()

    def patched(self, label: str, relative_script: str, replacements: dict[str, str], expected: list[Path]) -> None:
        script = ROOT / relative_script
        if not script.is_file():
            raise FileNotFoundError(script)
        if label in self.completed and all(p.is_file() and p.stat().st_size for p in expected):
            print(f"[resume] {label}", flush=True)
            return
        with tempfile.TemporaryDirectory(prefix="runner-script-", dir=self.out) as temp:
            target = Path(temp) / script.name
            target.write_bytes(patched_source(script, replacements))
            self.execute(label, [sys.executable, str(target)], expected)

    def cp(
        self,
        label: str,
        image: Path,
        out: Path,
        organoid_mask: Path | None = None,
    ) -> None:
        cp = self.cfg.get("cellprofiler", {})
        environment = cp.get("conda_env", "cellprofiler-native")
        conda = cp.get("conda_executable", "conda")
        pipeline = self.cfg["_pipeline"]
        if organoid_mask is not None:
            validate_shared_mask_canvas(image, organoid_mask)
            pipeline = self.shared_mask_pipeline
        cmd = [conda, "run", "--no-capture-output", "-n", environment, "python",
               str(ROOT / "CellProfiler/cellprofilerdriver.py"),
               "--pipeline", str(pipeline),
               "--input", str(image), "--output", str(out)]
        if organoid_mask is not None:
            cmd.extend(["--organoid-mask", str(organoid_mask)])
        expected = [
            out / "OrganoidMask.tiff",
            out / "CellMask.tiff",
            out / "MyExpt_FilterObjects2.csv",
            out / "MyExpt_FilterObjects.csv",
        ]
        self.execute(
            label,
            cmd,
            expected,
        )
        if organoid_mask is not None:
            try:
                verify_shared_mask_outputs(
                    organoid_mask,
                    out / "OrganoidMask.tiff",
                    out / "CellMask.tiff",
                )
            except Exception:
                self.completed.pop(label, None)
                self.save()
                raise


def ensure_projection(source: Path, destination: Path) -> Path:
    """Accept 2D channel TIFFs or project a single-channel ZYX stitched TIFF."""
    import numpy as np
    import tifffile
    description = inspect_tiff(source)
    if not description["projection_required"]:
        return source
    if destination.is_file() and destination.stat().st_size:
        projected = inspect_tiff(destination)
        if projected["projection_required"] or projected["yx_shape"] != description["yx_shape"]:
            raise ValueError(f"Existing projection is incompatible with its source: {destination}")
        return destination
    with tifffile.TiffFile(source) as tf:
        series = tf.series[0]
        shape, axes = series.shape, series.axes
        destination.parent.mkdir(parents=True, exist_ok=True)
        print(f"Projecting {source.name}, Z={shape[0]} (max across axis 0)", flush=True)
        if len(tf.pages) == shape[0] and all(len(p.shape) == 2 for p in tf.pages):
            projection = tf.pages[0].asarray().copy()
            for page in tf.pages[1:]:
                np.maximum(projection, page.asarray(), out=projection)
        else:
            projection = np.max(series.asarray(), axis=0)
        tifffile.imwrite(
            destination,
            projection,
            bigtiff=True,
            compression="zlib",
            metadata={"axes": "YX"},
        )
    return destination


def export_nuclear_positive(h5path: Path, csv_path: Path) -> None:
    import h5py
    import pandas as pd
    with h5py.File(h5path, "r") as h5:
        df = pd.read_csv(io.BytesIO(h5["tables/filtered_object_a_pixels"][()].tobytes()))
    if not {"X", "Y"}.issubset(df.columns) or not ({"ROI", "CellID"} & set(df.columns)):
        raise ValueError(f"Unexpected nuclear-positive HDF5 columns: {list(df.columns)}")
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(csv_path, index=False)


def input_rows(names: list[str], csvs: dict[str, Path], images: dict[str, Path]) -> list[dict]:
    return [{"name": n, "roi_csv": str(csvs[n]), "original": str(images[n])} for n in names]


def ensure_named_alignment(registered: Path, named: Path) -> Path:
    """Expose CASTalign's fixed output name under the actual marker name."""
    if named == registered:
        return registered
    if named.is_symlink():
        if named.resolve() == registered.resolve():
            return named
        named.unlink()
    elif named.exists():
        raise RuntimeError(f"Refusing to replace unexpected alignment output: {named}")
    named.symlink_to(registered.name)
    return named


def run(cfg: dict, resume: bool) -> None:
    r = Runner(cfg, resume)
    names = list(cfg["_channels"])
    markers = cfg["_markers"]
    nuc = cfg["_nuclear"]
    x = float(cfg["pixel_size_um"]["x"])
    y = float(cfg["pixel_size_um"]["y"])
    reports = cfg.get("reports", {})
    images: dict[str, Path] = {}

    # External stitching is already complete. Only channel max projections are needed.
    for name, input_path in cfg["_channels"].items():
        images[name] = ensure_projection(input_path, r.out / "00_max_projections" / f"{name}.tif")

    aligned = {nuc: images[nuc]}
    if cfg.get("alignment", {}).get("enabled", True):
        first_cp: dict[str, Path] = {}
        for name in names:
            first_cp[name] = r.out / "01_initial_segmentation" / name
            r.cp(f"02_initial_cp_{name}", images[name], first_cp[name])
        for name in markers:
            target = r.out / "02_registration" / name
            command = [sys.executable, str(ROOT / "Alignment/CASTalign_two_channel_registration.py"),
                       "--fixed-image", str(images[nuc]), "--moving-image", str(images[name]),
                       "--fixed-csv", str(first_cp[nuc] / "MyExpt_FilterObjects2.csv"),
                       "--moving-csv", str(first_cp[name] / "MyExpt_FilterObjects2.csv"),
                       "--fixed-x-col", "AreaShape_Center_X", "--fixed-y-col", "AreaShape_Center_Y",
                       "--moving-x-col", "AreaShape_Center_X", "--moving-y-col", "AreaShape_Center_Y",
                       "--pixel-size-x", str(x), "--pixel-size-y", str(y), "--output", str(target)]
            command += list(cfg.get("alignment", {}).get("extra_args", []))
            registered = target / "aligned_PV_max.tif"  # Hard-coded by existing CASTalign script.
            r.execute(f"03_align_{name}", command, [registered])
            named = target / f"aligned_{name}_max.tif"
            aligned[name] = ensure_named_alignment(registered, named)
    else:
        aligned.update({name: images[name] for name in markers})

    aligned_cp: dict[str, Path] = {}
    aligned_cp[nuc] = r.out / "03_aligned_segmentation" / nuc
    r.cp(f"04_aligned_cp_{nuc}", aligned[nuc], aligned_cp[nuc])
    reference_mask = aligned_cp[nuc] / "OrganoidMask.tiff"
    if not np.any(read_binary_mask(reference_mask, "Hoechst organoid mask")):
        raise RuntimeError(
            f"Final Hoechst organoid mask contains no foreground pixels: {reference_mask}"
        )
    for name in markers:
        aligned_cp[name] = r.out / "03_aligned_segmentation" / name
        r.cp(
            f"04_aligned_cp_{name}",
            aligned[name],
            aligned_cp[name],
            organoid_mask=reference_mask,
        )

    if reports.get("organoid_mask", True):
        out = r.out / "04_organoid_mask_qc"
        command = [
            sys.executable,
            str(ROOT / "QC/organoid_mask_qc.py"),
            "--output-dir",
            str(out),
            "--reference-mask",
            str(reference_mask),
        ]
        for name in names:
            command.extend(
                [
                    "--entry",
                    name,
                    str(aligned[name]),
                    str(aligned_cp[name] / "OrganoidMask.tiff"),
                    str(aligned_cp[name] / "CellMask.tiff"),
                ]
            )
        r.execute(
            "04_organoid_mask_qc",
            command,
            [
                out / "organoid_mask_qc.png",
                out / "organoid_mask_statistics.csv",
            ],
        )

    extracted: dict[str, Path] = {}
    for name in names:
        out = r.out / "04_roi_extraction" / name
        extracted[name] = out / f"{name}_ROI_pixels.csv"
        r.patched(f"04_roi_{name}", "ExtractingROIs/extract_rois.py", {
            "tiff_path": path_value(aligned_cp[name] / "CellMask.tiff"),
            "channel_csv_path": path_value(aligned_cp[name] / "MyExpt_FilterObjects2.csv"),
            "organoid_csv_path": path_value(aligned_cp[name] / "MyExpt_FilterObjects.csv"),
            "CHANNEL_NAME": repr(name), "output_dir": path_value(out),
        }, [extracted[name]])

    if reports.get("intensity", True):
        out = r.out / "05_intensity_report"
        r.patched("05_intensity_report", "ROIFiltering/Intensity/intensityfilter.py", {
            "output_dir": path_value(out), "channels": repr(input_rows(names, extracted, aligned)),
        }, [out / "All_Channel_ROI_Intensity_Report.pdf"])

    filter_settings = cfg.get("filters", {})
    intensity_csv: dict[str, Path] = {}
    intensity_out = r.out / "05_intensity_filtered"
    for name in names:
        intensity_csv[name] = intensity_out / name / f"{name}_ROI_pixels_FILTERED.csv"
    r.patched("05_intensity_gate", "ROIFiltering/FinalFilter/filter.py", {
        "output_dir": path_value(intensity_out),
        "channels": repr(input_rows(names, extracted, aligned)),
        "FILTERS": repr(filters_for_pass(filter_settings, "intensity")),
        "PIXEL_SIZE_X_UM": repr(x), "PIXEL_SIZE_Y_UM": repr(y),
    }, list(intensity_csv.values()))

    if reports.get("shape", True):
        for metric, script in SHAPE_REPORTS.items():
            out = r.out / "06_shape_reports" / metric
            r.patched(f"06_shape_report_{metric}", script, {
                "output_dir": path_value(out),
                "channels": repr(input_rows(names, intensity_csv, aligned)),
                "PIXEL_SIZE_X_UM": repr(x), "PIXEL_SIZE_Y_UM": repr(y),
            }, [out / "All_Channel_ROI_Shape_Report.pdf"])

    final_out = r.out / "06_final_shape_filtered"
    final_csv: dict[str, Path] = {
        name: final_out / name / f"{name}_ROI_pixels_FILTERED.csv" for name in names
    }
    r.patched("06_shape_gate", "ROIFiltering/FinalFilter/filter.py", {
        "output_dir": path_value(final_out),
        "channels": repr(input_rows(names, intensity_csv, aligned)),
        "FILTERS": repr(filters_for_pass(filter_settings, "shape")),
        "PIXEL_SIZE_X_UM": repr(x), "PIXEL_SIZE_Y_UM": repr(y),
    }, list(final_csv.values()))

    nuclear_h5: dict[str, Path] = {}
    associated_csv: dict[str, Path] = {}
    organoid = aligned_cp[nuc] / "MyExpt_FilterObjects.csv"

    def colocalize(label: str, a: str, b: str, a_csv: Path, b_csv: Path,
                   directory: Path, threshold: float, nuclear_sources=None) -> Path:
        h5path = directory / "analysis.h5"
        settings = {
            "OBJECT_A_NAME": repr(a), "OBJECT_B_NAME": repr(b),
            "OBJECT_A_CSV": path_value(a_csv), "OBJECT_B_CSV": path_value(b_csv),
            "OBJECT_A_TIFF": path_value(aligned[a]), "OBJECT_B_TIFF": path_value(aligned[b]),
            "ORGANOID_CSV": path_value(organoid), "RESULTS_FILE": path_value(h5path),
            "THRESHOLD": repr(float(threshold)),
            "PIXEL_SIZE_X_UM": repr(x), "PIXEL_SIZE_Y_UM": repr(y),
        }
        r.patched(label, "Colocalization/colocalizationdapiscript1.py", settings, [h5path])
        if reports.get("colocalization", True):
            viz_settings = {
                "RESULTS_FILE": path_value(h5path), "OUTPUT_DIR": path_value(directory / "visualizations"),
                "ENABLE_ORGANOID_AREA_SUMMARY": "True",
                "ORGANOID_MEASUREMENTS_CSV": path_value(organoid),
                "ENABLE_NUCLEAR_ALL_CELL_SUMMARY": "True" if nuclear_sources else "False",
                "OBJECT_A_NUCLEAR_RESULTS_FILE": path_value(nuclear_sources[0]) if nuclear_sources else "None",
                "OBJECT_B_NUCLEAR_RESULTS_FILE": path_value(nuclear_sources[1]) if nuclear_sources else "None",
            }
            r.patched(label + "_viz", "Colocalization/colocalizationdapiscript2.py", viz_settings,
                      [directory / "visualizations" / f"{a}_{b}_filtering_stage_counts.csv"])
        return h5path

    for marker in markers:
        directory = r.out / "07_nuclear_colocalization" / f"{marker}_vs_{nuc}"
        nuclear_h5[marker] = colocalize(f"07_match_{marker}", marker, nuc,
                                        final_csv[marker], final_csv[nuc], directory,
                                        cfg["thresholds"]["marker_nuclear"])
        associated_csv[marker] = directory / f"{marker}_nuclear_positive_cells.csv"
        label = f"07_export_{marker}"
        if (
            label not in r.completed
            or not associated_csv[marker].is_file()
            or not associated_csv[marker].stat().st_size
        ):
            export_nuclear_positive(nuclear_h5[marker], associated_csv[marker])
            r.completed[label] = [str(associated_csv[marker])]
            r.save()

    for a, b in cfg["_pairs"]:
        directory = r.out / "08_marker_colocalization" / f"{a}_vs_{b}"
        colocalize(f"08_pair_{a}_{b}", a, b, associated_csv[a], associated_csv[b],
                   directory, cfg["thresholds"]["marker_marker"],
                   (nuclear_h5[a], nuclear_h5[b]))
    print(f"\nFinished. Outputs and provenance: {r.out}\nManifest: {r.manifest}", flush=True)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", type=Path, required=True, help="JSON configuration; see stitched_pipeline.example.json")
    ap.add_argument("--resume", action="store_true", help="Resume only if configuration and input metadata match")
    ap.add_argument("--dry-run", action="store_true", help="Validate config and print the ordered stages without executing")
    args = ap.parse_args(argv)
    cfg = read_config(args.config.resolve())
    descriptions = validate_inputs(cfg)
    if args.dry_run:
        print("Configuration and stitched TIFF preflight: OK")
        for name, info in descriptions.items():
            action = "max projection required" if info["projection_required"] else "already 2D"
            print(
                f"  {name}: shape={info['shape']}, axes={info['axes']}, "
                f"dtype={info['dtype']} ({action})"
            )
        print("Ordered stages:")
        print("  00  maximum projection for the inputs marked above")
        if cfg.get("alignment", {}).get("enabled", True):
            print("  01  initial CellProfiler segmentation for registration")
            print("  02  CASTalign each marker to the nuclear reference")
        print("  03a final-frame Hoechst segmentation establishes tissue mask")
        if cfg["_markers"]:
            print("  03b marker-specific detection inside the shared Hoechst mask")
        if cfg.get("reports", {}).get("organoid_mask", True):
            print("  04a Hoechst-reference mask equality and containment QC")
        print("  04b ROI reconstruction and pixel-coordinate export")
        if cfg.get("reports", {}).get("intensity", True):
            print("  05a intensity QC report")
        print("  05b intensity-only filter pass")
        if cfg.get("reports", {}).get("shape", True):
            print("  06a shape QC reports on intensity survivors")
        print("  06b shape-only filter pass")
        if cfg["_markers"]:
            print("  07  marker/nuclear matching and nuclear-positive exports")
        if cfg["_pairs"]:
            print("  08  requested marker/marker comparisons")
        if cfg.get("reports", {}).get("colocalization", True) and cfg["_markers"]:
            print("      colocalization visualizations/QC enabled")
        print(f"Nuclear: {cfg['_nuclear']}; channels: {list(cfg['_channels'])}")
        print(f"Pairs (A-inside-B direction): {cfg['_pairs']}")
        print(f"CellProfiler pipeline: {cfg['_pipeline']}")
        print(f"Output: {cfg['_output']}")
        print("Runtime environments are checked when execution starts (without --dry-run).")
        return 0
    validate_runtime(cfg)
    run(cfg, resume=args.resume)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
