#!/usr/bin/env python3
"""Run the existing ImageAnalysisScripts stages from channel-specific stitched TIFFs.

No scientific analysis code is modified. Scripts with file-level USER SETTINGS are
copied temporarily and ONLY their top-level assignment values are replaced.
See STITCHED_WORKFLOW.md before using this on experimental data.
"""

from __future__ import annotations

import argparse
import ast
import hashlib
import io
import itertools
import json
import re
import subprocess
import sys
import tempfile
from pathlib import Path

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
        if set(rule) != {"percentile", "keep"}:
            raise ValueError(f"{name} must define percentile and keep")
        pct = rule["percentile"]
        if not isinstance(pct, (int, float)) or not 0 <= pct <= 100:
            raise ValueError(f"Invalid percentile for {name}: {pct}")
        if rule["keep"] not in ("above", "below"):
            raise ValueError(f"Invalid keep direction for {name}")
        defaults[name] = rule
    return defaults


def read_config(config_file: Path) -> dict:
    raw = json.loads(config_file.read_text(encoding="utf-8"))
    nuclear = raw["nuclear_channel"]
    raw_channels = raw["channels"]
    if not isinstance(raw_channels, dict) or not raw_channels or nuclear not in raw_channels:
        raise ValueError("channels must contain the named nuclear_channel")
    for name in raw_channels:
        if not re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*", name):
            raise ValueError(f"Use simple ASCII alphanumeric/underscore channel names: {name!r}")
    channels = {name: resolve_path(path, config_file.parent) for name, path in raw_channels.items()}
    marker_names = [name for name in channels if name != nuclear]
    pairs = raw.get("pairs", list(itertools.combinations(marker_names, 2)))
    for pair in pairs:
        if len(pair) != 2 or pair[0] == pair[1] or any(x not in marker_names for x in pair):
            raise ValueError(f"Pair must contain two distinct, non-nuclear channel names: {pair}")
    thresholds = raw.get("thresholds", {})
    if marker_names and "marker_nuclear" not in thresholds:
        raise ValueError("Explicit thresholds.marker_nuclear is required")
    if pairs and "marker_marker" not in thresholds:
        raise ValueError("Explicit thresholds.marker_marker is required")
    for label, v in thresholds.items():
        if not isinstance(v, (int, float)) or not 0 <= v <= 1:
            raise ValueError(f"{label} threshold must lie in [0, 1]")
    pixels = raw["pixel_size_um"]
    if float(pixels["x"]) <= 0 or float(pixels["y"]) <= 0:
        raise ValueError("Pixel sizes must be positive")
    filters = raw.get("filters", {})
    filters_for_pass(filters, "intensity")
    filters_for_pass(filters, "shape")
    cp = raw.get("cellprofiler", {})
    pipeline = resolve_path(cp.get("pipeline", "CellProfiler/MGEOPVFinal.cppipe"), ROOT)
    if not pipeline.exists():
        raise FileNotFoundError(f"CellProfiler pipeline missing: {pipeline}")
    result = dict(raw)
    result["_channels"] = channels
    result["_markers"] = marker_names
    result["_pairs"] = pairs
    result["_nuclear"] = nuclear
    result["_pipeline"] = pipeline
    result["_output"] = resolve_path(raw["output_dir"], config_file.parent)
    return result


class Runner:
    def __init__(self, cfg: dict, resume: bool):
        self.cfg = cfg
        self.out: Path = cfg["_output"]
        self.out.mkdir(parents=True, exist_ok=True)
        self.manifest = self.out / "run_manifest.json"
        # Include input file metadata so --resume cannot silently reuse changed inputs.
        file_info = {
            name: {"path": str(path), "size": path.stat().st_size, "mtime_ns": path.stat().st_mtime_ns}
            for name, path in cfg["_channels"].items()
        }
        digest_data = {k: v for k, v in cfg.items() if not k.startswith("_")}
        self.digest = hashlib.sha256(json.dumps({"config": digest_data, "inputs": file_info}, sort_keys=True).encode()).hexdigest()
        if self.manifest.exists():
            old = json.loads(self.manifest.read_text())
            if not resume:
                raise RuntimeError(f"Output already has a manifest; choose a new output_dir or --resume: {self.out}")
            if old["config_input_sha256"] != self.digest:
                raise RuntimeError("--resume refused: config or input file metadata changed. Use a fresh output_dir.")
            self.completed = old["completed"]
        else:
            if resume and any(self.out.iterdir()):
                raise RuntimeError("Cannot resume an untracked output directory")
            if any(self.out.iterdir()):
                raise RuntimeError(f"Refusing to mix with existing untracked files: {self.out}")
            self.completed = {}
            self.save()

    def save(self) -> None:
        payload = {"config_input_sha256": self.digest, "completed": self.completed}
        temp = self.manifest.with_suffix(".json.tmp")
        temp.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        temp.replace(self.manifest)

    def execute(self, label: str, cmd: list[str], expected: list[Path]) -> None:
        if label in self.completed and all(p.is_file() and p.stat().st_size for p in expected):
            print(f"[resume] {label}", flush=True)
            return
        print(f"\n{'=' * 70}\nRUNNING {label}\n{'=' * 70}", flush=True)
        subprocess.run(cmd, cwd=ROOT, check=True)
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

    def cp(self, label: str, image: Path, out: Path) -> None:
        cp = self.cfg.get("cellprofiler", {})
        environment = cp.get("conda_env", "cellprofiler-native")
        conda = cp.get("conda_executable", "conda")
        cmd = [conda, "run", "--no-capture-output", "-n", environment, "python",
               str(ROOT / "CellProfiler/cellprofilerdriver.py"),
               "--pipeline", str(self.cfg["_pipeline"]),
               "--input", str(image), "--output", str(out)]
        self.execute(label, cmd, [out / "CellMask.tiff", out / "MyExpt_FilterObjects2.csv", out / "MyExpt_FilterObjects.csv"])


def ensure_projection(source: Path, destination: Path) -> Path:
    """Accept 2D channel TIFFs or project a single-channel ZYX stitched TIFF."""
    import numpy as np
    import tifffile
    with tifffile.TiffFile(source) as tf:
        series = tf.series[0]
        shape, axes = series.shape, series.axes
        if len(shape) == 2:
            return source
        if len(shape) != 3 or axes not in ("ZYX", "QYX", "IYX") :
            raise ValueError(f"Expected one-channel YX/ZYX TIFF, got {shape}, axes={axes} in {source}")
        if destination.is_file():
            return destination
        destination.parent.mkdir(parents=True, exist_ok=True)
        print(f"Projecting {source.name}, Z={shape[0]} (max across axis 0)", flush=True)
        if len(tf.pages) == shape[0] and all(len(p.shape) == 2 for p in tf.pages):
            projection = tf.pages[0].asarray()
            for page in tf.pages[1:]:
                np.maximum(projection, page.asarray(), out=projection)
        else:
            projection = np.max(series.asarray(), axis=0)
        tifffile.imwrite(destination, projection, bigtiff=True, compression="zlib")
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


def run(cfg: dict, resume: bool) -> None:
    for n, source in cfg["_channels"].items():
        if not source.is_file():
            raise FileNotFoundError(f"Stitched input for {n} does not exist: {source}")
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

    first_cp: dict[str, Path] = {}
    for name in names:
        first_cp[name] = r.out / "01_initial_segmentation" / name
        r.cp(f"02_initial_cp_{name}", images[name], first_cp[name])

    aligned = {nuc: images[nuc]}
    if cfg.get("alignment", {}).get("enabled", True):
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
            if name != "PV" and not named.exists():
                named.symlink_to(registered.name)  # No second multi-GB TIFF copy.
            aligned[name] = named
    else:
        aligned.update({name: images[name] for name in markers})

    aligned_cp: dict[str, Path] = {}
    for name in names:
        aligned_cp[name] = r.out / "03_aligned_segmentation" / name
        r.cp(f"04_aligned_cp_{name}", aligned[name], aligned_cp[name])

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
        if label not in r.completed or not associated_csv[marker].is_file():
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
    if args.dry_run:
        print("Stitched TIFF -> max projection (when ZYX) -> initial CP -> CASTalign -> aligned CP")
        print(" -> ROI extraction -> intensity report/gate -> shape reports/gate -> marker/nuclear")
        print(" -> nuclear-positive CSV export -> marker/marker -> visualization/QC")
        print(f"Nuclear: {cfg['_nuclear']}; channels: {list(cfg['_channels'])}")
        print(f"Pairs (A-inside-B direction): {cfg['_pairs']}")
        print(f"Output: {cfg['_output']}")
        return 0
    run(cfg, resume=args.resume)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
