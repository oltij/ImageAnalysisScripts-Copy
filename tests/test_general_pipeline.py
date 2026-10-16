import contextlib
import hashlib
import io
import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import h5py
import numpy as np
import tifffile

from GeneralPipeline import batch_runner
from GeneralPipeline import channel_metadata
from GeneralPipeline import discovery
from GeneralPipeline import manifest as manifest_module
from GeneralPipeline import run_pipeline


ROOT = Path(__file__).resolve().parents[1]


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def analysis_settings(alignment=True):
    return {
        "cellprofiler": {
            "pipeline": "CellProfiler/MGEOPVFinal.cppipe",
            "conda_env": "cellprofiler-native",
            "conda_executable": "conda",
        },
        "alignment": {"enabled": alignment, "extra_args": []},
        "pairs": "auto",
        "thresholds": {"marker_nuclear": 0.25, "marker_marker": 0.25},
        "filters": {"intensity": {}, "shape": {}},
        "reports": {
            "organoid_mask": True,
            "intensity": False,
            "shape": False,
            "colocalization": False,
        },
    }


def write_ims(path, biological_channels, *, stored_names=None, shape=(2, 8, 9)):
    path.parent.mkdir(parents=True, exist_ok=True)
    stored_names = stored_names or [
        f"Physical channel {index}" for index in range(len(biological_channels))
    ]
    with h5py.File(path, "w") as handle:
        image = handle.create_group("DataSetInfo/Image")
        image.attrs["Unit"] = np.bytes_("um")
        for axis, (pixels, spacing) in enumerate(
            zip((shape[2], shape[1], shape[0]), (0.4, 0.5, 1.2))
        ):
            image.attrs[f"ExtMin{axis}"] = np.bytes_("0")
            image.attrs[f"ExtMax{axis}"] = np.bytes_(str(pixels * spacing))
            image.attrs[("X", "Y", "Z")[axis]] = np.bytes_(str(pixels))
        timepoint = handle.create_group("DataSet/ResolutionLevel 0/TimePoint 0")
        for index, name in enumerate(stored_names):
            group = timepoint.create_group(f"Channel {index}")
            group.create_dataset("Data", data=np.zeros(shape, dtype=np.uint16))
            info = handle.create_group(f"DataSetInfo/Channel {index}")
            info.attrs["Name"] = np.bytes_(name)
            info.attrs["EmissionWavelength"] = np.bytes_(str(405 + index * 70))


def create_field(
    root,
    channels,
    *,
    group_prefix="Sample",
    tile_ids=(2, 4, 6, 8),
    session="Experiment_A",
    nested=False,
    set_number=None,
    acquisition_identifier=None,
    stored_names=None,
):
    markers = "_".join(channels)
    set_token = f"_Set{set_number}" if set_number is not None else ""
    identifier = (
        f"_{acquisition_identifier}" if acquisition_identifier is not None else ""
    )
    group = f"{group_prefix}{set_token}_{markers}{identifier}"
    parent = root / session / "2026-06-18" if session is not None else root
    field = parent / f"{group}_Field_1"
    sources = []
    for tile_id in tile_ids:
        filename = f"{group}_F{tile_id:02d}.ims"
        path = field / f"tile_F{tile_id:02d}" / filename if nested else field / filename
        write_ims(path, channels, stored_names=stored_names)
        sources.append(path)
    return field, sources


def raw_config(path, *, channel_names="auto", nuclear_channel="auto"):
    payload = {
        "input_mode": "raw_ims_tiles",
        "channel_names": channel_names,
        "nuclear_channel": nuclear_channel,
        "filename_markers": [],
        "raw_ims_tiles": {
            "grid": [2, 2],
            "max_depth": 5,
            "require_contiguous_tiles": False,
            "skip_incomplete": True,
            "sessions": [],
            "group_pattern": None,
        },
        "preprocessing": {
            "fiji_sif": "/software/fiji_bigstitcher.sif",
            "ims_unit": "um",
            "trim_zero_padding": True,
            "analysis_downsampling": {"factor": 1, "method": "area_mean"},
            "stitching": {"grid": [2, 2]},
        },
        "analysis": analysis_settings(),
    }
    path.write_text(json.dumps(payload), encoding="utf-8")
    return path


def write_tiff(path, shape, axes):
    path.parent.mkdir(parents=True, exist_ok=True)
    tifffile.imwrite(
        path,
        np.zeros(shape, dtype=np.uint16),
        metadata={"axes": axes},
        photometric="rgb" if axes == "YXS" else "minisblack",
    )


def stitched_config(path, samples=None, *, mapping=None, alignment=True, channels=None):
    channels = channels or ["DAPI", "GFP", "LHX6"]
    payload = {
        "input_mode": "stitched_tiffs",
        "channel_names": channels,
        "nuclear_channel": "DAPI",
        "pixel_size_um": {"x": 0.4, "y": 0.5},
        "stitched_tiffs": {"mapping": mapping or {"samples": samples}},
        "analysis": analysis_settings(alignment=alignment),
    }
    path.write_text(json.dumps(payload), encoding="utf-8")
    return path


class RawIMSDiscoveryTests(unittest.TestCase):
    def settings(self, names="auto"):
        return {
            "grid": [2, 2],
            "max_depth": 5,
            "require_contiguous_tiles": False,
            "skip_incomplete": True,
            "sessions": [],
            "group_pattern": None,
            "channel_names": names,
            "nuclear_channel": "auto",
            "filename_markers": channel_metadata.marker_vocabulary(
                names if isinstance(names, list) else []
            ),
            "ims_unit": "um",
        }

    def test_nested_discovery_preserves_original_tile_order_and_optional_session(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            field, sources = create_field(
                root, ["DAPI", "GFP", "LHX6"], nested=True, set_number=3
            )
            acquisitions = discovery.discover_raw_acquisitions(root, self.settings())
            self.assertEqual(len(acquisitions), 1)
            acquisition = acquisitions[0]
            self.assertEqual(acquisition.path, field.resolve())
            self.assertEqual(acquisition.session, "Experiment_A")
            self.assertEqual(acquisition.experimental_set_number, 3)
            self.assertEqual(list(acquisition.channels), ["DAPI", "GFP", "LHX6"])
            self.assertEqual(
                [tile.source_index for tile in acquisition.tiles], [2, 4, 6, 8]
            )
            self.assertEqual(
                [tile.path for tile in acquisition.tiles],
                [path.resolve() for path in sources],
            )

    def test_no_session_directory_is_required(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            create_field(root, ["DAPI", "GFP"], session=None)
            acquisition = discovery.discover_raw_acquisitions(root, self.settings())[0]
            self.assertIsNone(acquisition.session)

    def test_two_three_and_four_channel_suffixes_are_supported(self):
        panels = (
            ["DAPI", "GFP"],
            ["DAPI", "ASCL1", "NKX2_1"],
            ["DAPI", "GFP", "BiVe3", "PV"],
        )
        for panel in panels:
            with self.subTest(panel=panel), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                create_field(root, panel)
                acquisition = discovery.discover_raw_acquisitions(
                    root, self.settings()
                )[0]
                self.assertEqual(list(acquisition.channels), panel)

    def test_set_identifier_is_metadata_not_a_channel(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            create_field(root, ["DAPI", "GFP"], set_number=3)
            acquisition = discovery.discover_raw_acquisitions(root, self.settings())[0]
            self.assertEqual(list(acquisition.channels), ["DAPI", "GFP"])
            self.assertEqual(acquisition.experimental_set_number, 3)

    def test_standalone_acquisition_identifier_must_match(self):
        markers = channel_metadata.marker_vocabulary()
        self.assertIsNone(
            channel_metadata.optional_acquisition_identifier(
                "Sample_DAPI_NKX2_1", markers, filename=False
            )
        )
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            field, sources = create_field(
                root, ["DAPI", "GFP"], acquisition_identifier=7
            )
            bad = sources[0]
            replacement = bad.with_name(bad.name.replace("_7_F", "_8_F"))
            bad.rename(replacement)
            with self.assertRaisesRegex(ValueError, "Acquisition identifier mismatch"):
                discovery.discover_raw_acquisitions(root, self.settings())

    def test_incomplete_fields_are_recorded_but_duplicate_tiles_are_fatal(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            create_field(root, ["DAPI", "GFP"], tile_ids=(0, 1, 2))
            create_field(root, ["DAPI", "GFP"], group_prefix="Complete")
            issues = []
            acquisitions = discovery.discover_raw_acquisitions(
                root, self.settings(), issues=issues
            )
            self.assertEqual(len(acquisitions), 1)
            self.assertEqual(issues[0]["kind"], "incomplete")

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            field, _ = create_field(root, ["DAPI", "GFP"], tile_ids=(0, 1, 2, 3))
            write_ims(field / "duplicate_DAPI_GFP_F00.ims", ["DAPI", "GFP"])
            with self.assertRaisesRegex(ValueError, "Duplicate tile index"):
                discovery.discover_raw_acquisitions(root, self.settings())

    def test_fusionstitcher_outputs_are_not_rediscovered_as_raw_tiles(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            field, _ = create_field(root, ["DAPI", "GFP"], tile_ids=(0, 1, 2, 3))
            (field / "Sample_FusionStitcher_F00.ims").write_bytes(
                b"derived stitched output, not an original tile"
            )
            acquisitions = discovery.discover_raw_acquisitions(root, self.settings())
            self.assertEqual(len(acquisitions), 1)
            self.assertEqual(
                [tile.source_index for tile in acquisitions[0].tiles], [0, 1, 2, 3]
            )

    def test_phase_channel_is_rejected_without_index_shifting(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            create_field(
                root,
                ["DAPI", "GFP", "LHX6", "PV"],
                stored_names=["Phase Contrast", "Blue", "Green", "Red"],
            )
            with self.assertRaisesRegex(ValueError, "unsupported physical phase"):
                discovery.discover_raw_acquisitions(root, self.settings())

    def test_explicit_mapping_is_required_when_filenames_do_not_name_channels(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            field = root / "Experiment_A/RawSample_Field_1"
            for tile_id in (0, 1, 2, 3):
                write_ims(
                    field / f"RawSample_F{tile_id:02d}.ims",
                    ["DAPI", "GFP"],
                    stored_names=["Confocal - Blue", "Confocal - Green"],
                )
            with self.assertRaisesRegex(ValueError, "biological-channel suffix"):
                discovery.discover_raw_acquisitions(root, self.settings())
            acquisition = discovery.discover_raw_acquisitions(
                root, self.settings(["DAPI", "GFP"])
            )[0]
            self.assertEqual(acquisition.channel_naming_basis, "explicit")
            self.assertEqual(list(acquisition.channels), ["DAPI", "GFP"])

    def test_cross_tile_channel_metadata_mismatch_is_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            _field, sources = create_field(root, ["DAPI", "GFP"])
            write_ims(sources[-1], ["DAPI", "GFP"], stored_names=["Blue", "Changed"])
            with self.assertRaisesRegex(ValueError, "differ across IMS tiles"):
                discovery.discover_raw_acquisitions(root, self.settings())


class StitchedTIFFDiscoveryTests(unittest.TestCase):
    def test_explicit_mapping_accepts_yx_and_zyx_and_different_registered_canvases(
        self,
    ):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            write_tiff(root / "s_DAPI.tif", (12, 13), "YX")
            write_tiff(root / "s_GFP.tif", (3, 10, 11), "ZYX")
            write_tiff(root / "s_LHX6.tif", (9, 12), "YX")
            settings = {
                "channel_names": ["DAPI", "GFP", "LHX6"],
                "nuclear_channel": "DAPI",
                "pixel_size_um": {"x": 0.4, "y": 0.5},
                "mapping": {
                    "samples": [
                        {
                            "name": "s",
                            "channels": {
                                "DAPI": "s_DAPI.tif",
                                "GFP": "s_GFP.tif",
                                "LHX6": "s_LHX6.tif",
                            },
                        }
                    ]
                },
                "analysis": analysis_settings(alignment=True),
            }
            downstream = run_pipeline._downstream_module(ROOT)
            acquisition = discovery.discover_stitched_acquisitions(
                root, settings, downstream.inspect_tiff
            )[0]
            self.assertFalse(acquisition.tiff_metadata["DAPI"]["projection_required"])
            self.assertTrue(acquisition.tiff_metadata["GFP"]["projection_required"])

    def test_rgb_and_unsupported_hyperstack_are_rejected(self):
        downstream = run_pipeline._downstream_module(ROOT)
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "rgb.tif"
            write_tiff(path, (8, 9, 3), "YXS")
            with self.assertRaisesRegex(ValueError, "single-channel YX image or ZYX"):
                downstream.inspect_tiff(path)

    def test_missing_channel_and_alignment_disabled_dimension_mismatch_fail(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            write_tiff(root / "s_DAPI.tif", (12, 13), "YX")
            write_tiff(root / "s_GFP.tif", (10, 11), "YX")
            downstream = run_pipeline._downstream_module(ROOT)
            base = {
                "channel_names": ["DAPI", "GFP"],
                "nuclear_channel": "DAPI",
                "pixel_size_um": {"x": 0.4, "y": 0.5},
                "analysis": analysis_settings(alignment=False),
            }
            with self.assertRaisesRegex(ValueError, "missing=.*GFP"):
                discovery.discover_stitched_acquisitions(
                    root,
                    {
                        **base,
                        "mapping": {
                            "samples": [
                                {"name": "s", "channels": {"DAPI": "s_DAPI.tif"}}
                            ]
                        },
                    },
                    downstream.inspect_tiff,
                )
            with self.assertRaisesRegex(ValueError, "different XY dimensions"):
                discovery.discover_stitched_acquisitions(
                    root,
                    {
                        **base,
                        "mapping": {
                            "samples": [
                                {
                                    "name": "s",
                                    "channels": {
                                        "DAPI": "s_DAPI.tif",
                                        "GFP": "s_GFP.tif",
                                    },
                                }
                            ]
                        },
                    },
                    downstream.inspect_tiff,
                )

    def test_named_regex_groups_multiple_samples_and_rejects_duplicates(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for sample in ("a", "b"):
                for channel in ("DAPI", "GFP"):
                    write_tiff(root / f"{sample}__{channel}.tif", (8, 9), "YX")
            downstream = run_pipeline._downstream_module(ROOT)
            settings = {
                "channel_names": ["DAPI", "GFP"],
                "nuclear_channel": "DAPI",
                "pixel_size_um": {"x": 0.4, "y": 0.5},
                "mapping": {
                    "filename_regex": r"(?P<sample>.+)__(?P<channel>DAPI|GFP)\.tif"
                },
                "analysis": analysis_settings(),
            }
            acquisitions = discovery.discover_stitched_acquisitions(
                root, settings, downstream.inspect_tiff
            )
            self.assertEqual([item.name for item in acquisitions], ["a", "b"])

            write_tiff(root / "a__DAPI_copy.tif", (8, 9), "YX")
            settings["mapping"] = {
                "filename_regex": r"(?P<sample>.+)__(?P<channel>DAPI|GFP)(?:_copy)?\.tif"
            }
            with self.assertRaisesRegex(ValueError, "Duplicate stitched TIFF"):
                discovery.discover_stitched_acquisitions(
                    root, settings, downstream.inspect_tiff
                )

    def test_explicit_sample_manifest_is_supported(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            inputs = root / "inputs"
            write_tiff(inputs / "s_DAPI.tif", (8, 9), "YX")
            write_tiff(inputs / "s_GFP.tif", (8, 9), "YX")
            sample_manifest = root / "samples.json"
            sample_manifest.write_text(
                json.dumps(
                    {
                        "schema": 1,
                        "samples": [
                            {
                                "name": "s",
                                "channels": {
                                    "DAPI": "s_DAPI.tif",
                                    "GFP": "s_GFP.tif",
                                },
                            }
                        ],
                    }
                )
            )
            downstream = run_pipeline._downstream_module(ROOT)
            acquisitions = discovery.discover_stitched_acquisitions(
                inputs,
                {
                    "channel_names": ["DAPI", "GFP"],
                    "nuclear_channel": "DAPI",
                    "pixel_size_um": {"x": 0.4, "y": 0.5},
                    "mapping": {"sample_manifest": "samples.json"},
                    "analysis": analysis_settings(),
                },
                downstream.inspect_tiff,
                manifest_base=root,
            )
            self.assertEqual([item.name for item in acquisitions], ["s"])

    def test_calibration_is_mandatory_and_never_defaulted(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "config.json"
            stitched_config(path, samples=[])
            payload = json.loads(path.read_text())
            payload.pop("pixel_size_um")
            path.write_text(json.dumps(payload))
            with self.assertRaisesRegex(ValueError, "pixel_size_um"):
                run_pipeline.load_config(path, "stitched_tiffs", ROOT)


class GeneralPipelinePreparationTests(unittest.TestCase):
    def test_raw_prepare_is_source_immutable_and_records_shared_handoff(self):
        with tempfile.TemporaryDirectory() as temporary:
            temporary = Path(temporary)
            source = temporary / "source"
            _field, paths = create_field(source, ["DAPI", "GFP"])
            before = {path: sha256(path) for path in paths}
            config = raw_config(temporary / "raw.json")
            workspace = temporary / "workspace"
            with contextlib.redirect_stdout(io.StringIO()):
                run_pipeline.prepare(
                    "raw_ims_tiles",
                    source.resolve(),
                    config.resolve(),
                    workspace.resolve(),
                    ROOT,
                    allow_dirty_repository=True,
                )
            payload = manifest_module.load_manifest(workspace / "manifest.json")
            job = payload["jobs"][0]
            self.assertEqual(
                job["downstream_runner"], str(ROOT / "run_stitched_pipeline.py")
            )
            self.assertEqual(job["bypassed_stages"], [])
            self.assertEqual(
                [item["source_tile_index"] for item in job["sources"]], [2, 4, 6, 8]
            )
            self.assertEqual(
                [Path(item["staged"]).name[-7:] for item in job["sources"]],
                ["F00.ims", "F01.ims", "F02.ims", "F03.ims"],
            )
            self.assertTrue(
                (workspace / "controller_compat/run_ims_pipeline.py").is_file()
            )
            raw_command = batch_runner._raw_preflight(
                job,
                json.loads(Path(job["config_path"]).read_text(encoding="utf-8")),
                Path(payload["compat_repository"]),
            )
            self.assertIn("--cleanup-stitching-intermediates", raw_command)
            self.assertEqual({path: sha256(path) for path in paths}, before)
            summary = manifest_module.validate_manifest(payload)
            self.assertEqual(summary["sources_checked"], 4)
            generated = Path(job["config_path"])
            generated.write_text(generated.read_text() + "\n")
            with self.assertRaisesRegex(ValueError, "configuration changed"):
                manifest_module.validate_manifest(payload)

    def test_stitched_prepare_bypasses_preprocessing_and_uses_same_downstream_runner(
        self,
    ):
        with tempfile.TemporaryDirectory() as temporary:
            temporary = Path(temporary)
            source = temporary / "tiffs"
            channels = ["DAPI", "GFP", "LHX6"]
            for channel in channels:
                write_tiff(source / f"s_{channel}.tif", (8, 9), "YX")
            samples = [
                {
                    "name": "s",
                    "channels": {channel: f"s_{channel}.tif" for channel in channels},
                }
            ]
            config = stitched_config(temporary / "stitched.json", samples)
            workspace = temporary / "workspace"
            with contextlib.redirect_stdout(io.StringIO()):
                run_pipeline.prepare(
                    "stitched_tiffs",
                    source.resolve(),
                    config.resolve(),
                    workspace.resolve(),
                    ROOT,
                    allow_dirty_repository=True,
                )
            payload = manifest_module.load_manifest(workspace / "manifest.json")
            job = payload["jobs"][0]
            self.assertEqual(
                job["bypassed_stages"], ["IMS export", "BaSiC", "BigStitcher"]
            )
            self.assertEqual(job["preprocessing_stages"], [])
            self.assertEqual(
                job["downstream_runner"], str(ROOT / "run_stitched_pipeline.py")
            )
            generated = json.loads(Path(job["config_path"]).read_text())
            self.assertEqual(generated["pixel_size_um"], {"x": 0.4, "y": 0.5})
            self.assertEqual(generated["nuclear_channel"], "DAPI")

    def test_print_only_submission_never_calls_slurm(self):
        with tempfile.TemporaryDirectory() as temporary:
            temporary = Path(temporary)
            source = temporary / "tiffs"
            write_tiff(source / "s_DAPI.tif", (8, 9), "YX")
            write_tiff(source / "s_GFP.tif", (8, 9), "YX")
            samples = [
                {
                    "name": "s",
                    "channels": {"DAPI": "s_DAPI.tif", "GFP": "s_GFP.tif"},
                }
            ]
            config = stitched_config(
                temporary / "stitched.json", samples, channels=["DAPI", "GFP"]
            )
            workspace = temporary / "workspace"
            with contextlib.redirect_stdout(io.StringIO()):
                run_pipeline.prepare(
                    "stitched_tiffs",
                    source.resolve(),
                    config.resolve(),
                    workspace.resolve(),
                    ROOT,
                    allow_dirty_repository=True,
                )
            with mock.patch.object(batch_runner.shutil, "which") as which_mock:
                with contextlib.redirect_stdout(io.StringIO()):
                    result = batch_runner.submit(
                        workspace / "manifest.json", 2, print_only=True
                    )
            self.assertIsNone(result)
            which_mock.assert_not_called()
            self.assertFalse((workspace / "slurm_submission.json").exists())

    def test_nonempty_workspace_is_never_overwritten(self):
        with tempfile.TemporaryDirectory() as temporary:
            workspace = Path(temporary) / "workspace"
            workspace.mkdir()
            (workspace / "keep.txt").write_text("user data")
            with self.assertRaisesRegex(ValueError, "not empty"):
                manifest_module.fresh_workspace(workspace)
            self.assertEqual((workspace / "keep.txt").read_text(), "user data")

    def test_absolute_turbo_tiff_mapping_still_requires_scratch_workspace(self):
        with tempfile.TemporaryDirectory() as temporary:
            temporary = Path(temporary)
            turbo = temporary / "turbo"
            input_root = temporary / "mapping_root"
            input_root.mkdir()
            write_tiff(turbo / "s_DAPI.tif", (8, 9), "YX")
            write_tiff(turbo / "s_GFP.tif", (8, 9), "YX")
            samples = [
                {
                    "name": "s",
                    "channels": {
                        "DAPI": str(turbo / "s_DAPI.tif"),
                        "GFP": str(turbo / "s_GFP.tif"),
                    },
                }
            ]
            config = stitched_config(
                temporary / "stitched.json", samples, channels=["DAPI", "GFP"]
            )
            with (
                mock.patch.object(manifest_module, "TURBO_IMS_ROOT", turbo),
                mock.patch.object(
                    manifest_module, "SCRATCH_ROOT", temporary / "scratch"
                ),
                self.assertRaisesRegex(ValueError, "require a workspace below"),
            ):
                run_pipeline.prepare(
                    "stitched_tiffs",
                    input_root.resolve(),
                    config.resolve(),
                    (temporary / "home/workspace").resolve(),
                    ROOT,
                    allow_dirty_repository=True,
                )


class BackwardCompatibilityTests(unittest.TestCase):
    EXPECTED_HASHES = {
        "IGF Organoids/run_folder_batch.py": "35660bcdc2317cbff4a34363f478351ad585e23ce1ef085af66299e0275a56be",
        "IGF Organoids/submit_folder_batch.py": "b44c10db069a490dea979306d15b5ee94ef94601481ec091d2697645110ba99a",
        "IGF Organoids/submit_production_batches.py": "d1e3868f3733529441f97328c109f1d8e7ed38521cdc0efc01271302ea359823",
        "IGF Organoids/run_folder_batch.sbatch": "ee4cbd4916bed6af06ff23a49649aef461b15047eccd08884d7b7e4eebd2a184",
        "IGF Organoids/configs/ims_mge_generic.example.json": "55b2633bb6b4d99cbf86d1f02d3aaa65ce0e530741532146aac9c1cd236937db",
        "IGF Organoids/configs/production_94_fields.json": "dff3f19706aeb11de6914d1a689c260ba9917da8c27a46d5a4644e84b584e22b",
        "Alignment/CASTalign_two_channel_registration.py": "0b6100725ee98f3cc58596f6a80d2d751031c15c2732018edf4649a00808c795",
        "run_stitched_pipeline.py": "53346f12c06d69d6895d69e8e76ffaec7308524cb718ad75ffc13a6ec21026d5",
        "run_ims_pipeline.py": "7c99fe219eefe944efd818a81f07743bf10f59aea7dc676415aeaf6f178496c7",
    }

    def test_igf_entrypoints_configs_and_scientific_controllers_are_byte_unchanged(
        self,
    ):
        actual = {
            relative: sha256(ROOT / relative) for relative in self.EXPECTED_HASHES
        }
        self.assertEqual(actual, self.EXPECTED_HASHES)

    def test_igf_94_field_plan_and_paths_remain_available(self):
        plan = json.loads(
            (ROOT / "IGF Organoids/configs/production_94_fields.json").read_text()
        )
        self.assertEqual(plan["expected"]["complete_fields"], 94)
        self.assertEqual(plan["expected"]["original_tiles"], 376)
        self.assertEqual(plan["expected"]["incomplete_fields"], 19)
        self.assertEqual(len(plan["sessions"]), 12)
        self.assertTrue((ROOT / "IGF Organoids/run_folder_batch.py").is_file())
        self.assertTrue((ROOT / "IGF Organoids/run_folder_batch.sbatch").is_file())


if __name__ == "__main__":
    unittest.main()
