import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import h5py
import numpy as np
import tifffile

import run_ims_pipeline as ims


def create_synthetic_ims(path, channel_count=4, shape=(3, 6, 7)):
    wavelengths = [405.0, 488.0, 561.0, 640.0]
    names = [
        "Confocal - Blue",
        "Confocal - Green",
        "Confocal - Red",
        "Confocal - Far Red",
    ]
    with h5py.File(path, "w") as handle:
        image = handle.create_group("DataSetInfo/Image")
        image.attrs["Unit"] = np.bytes_("um")
        for axis, (pixels, spacing) in enumerate(
            zip((shape[2], shape[1], shape[0]), (0.5, 0.6, 1.25))
        ):
            image.attrs[f"ExtMin{axis}"] = np.bytes_("0")
            image.attrs[f"ExtMax{axis}"] = np.bytes_(str(pixels * spacing))
            image.attrs[("X", "Y", "Z")[axis]] = np.bytes_(str(pixels))

        timepoint = handle.create_group("DataSet/ResolutionLevel 0/TimePoint 0")
        for index in range(channel_count):
            data = np.zeros(shape, dtype=np.uint16)
            data[:, :-1, :-1] = (index + 1) * 100 + np.arange(
                shape[0], dtype=np.uint16
            )[:, None, None]
            channel = timepoint.create_group(f"Channel {index}")
            channel.create_dataset("Data", data=data)
            info = handle.create_group(f"DataSetInfo/Channel {index}")
            info.attrs["Name"] = np.bytes_(
                names[index] if index < len(names) else f"Extra {index}"
            )
            info.attrs["EmissionWavelength"] = np.bytes_(
                str(wavelengths[index] if index < len(wavelengths) else 700 + index)
            )


def load_mosaic_module(repo_root):
    path = repo_root / "Stitching/IMSTIFF_trim_zero_padding.py"
    spec = importlib.util.spec_from_file_location("ims_mosaic_script", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


class IMSPipelineTests(unittest.TestCase):
    def _analysis_settings(self):
        return {
            "cellprofiler": {
                "pipeline": "CellProfiler/MGEOPVFinal.cppipe",
                "conda_env": "cellprofiler-native",
                "conda_executable": "conda",
            },
            "alignment": {"enabled": False, "extra_args": []},
            "pairs": [["mNeonGreen", "BiVe3"]],
            "thresholds": {"marker_nuclear": 0.25, "marker_marker": 0.25},
            "filters": {"intensity": {}, "shape": {}},
            "reports": {"intensity": False, "shape": False, "colocalization": False},
        }

    def _write_config(self, directory, source, output="output", mode="auto"):
        config = {
            "ims_input": str(source),
            "input_mode": mode,
            "output_dir": output,
            "fiji_sif": None,
            "channel_order": list(ims.CHANNEL_ORDER),
            "ims_unit": "um",
            "trim_zero_padding": True,
            "stitching": {"grid": [2, 2]},
            "analysis": self._analysis_settings(),
        }
        path = Path(directory) / "ims_config.json"
        path.write_text(json.dumps(config), encoding="utf-8")
        return path

    def test_inspection_extracts_calibration_and_fixed_channel_mapping(self):
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / "already_stitched.ims"
            create_synthetic_ims(source)
            info = ims.inspect_ims(source, "um")
            self.assertEqual(
                [channel["canonical_name"] for channel in info["channels"]],
                list(ims.CHANNEL_ORDER),
            )
            self.assertEqual(
                [channel["ims_name"] for channel in info["channels"]],
                [
                    "Confocal - Blue",
                    "Confocal - Green",
                    "Confocal - Red",
                    "Confocal - Far Red",
                ],
            )
            self.assertEqual(info["voxel_size_um"], {"x": 0.5, "y": 0.6, "z": 1.25})
            self.assertTrue(info["wavelength_order_verified"])

    def test_wrong_channel_count_fails_before_export(self):
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / "three_channels.ims"
            create_synthetic_ims(source, channel_count=3)
            with self.assertRaisesRegex(ValueError, "requires exactly 4"):
                ims.inspect_ims(source, "um")

    def test_single_mosaic_field_discovers_all_siblings(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            for field in range(4):
                create_synthetic_ims(tmp / f"Sample_F{field:02d}.ims")
            acquisitions = ims.discover_acquisitions(tmp / "Sample_F00.ims", "auto")
            self.assertEqual(len(acquisitions), 1)
            self.assertEqual(acquisitions[0].name, "Sample")
            self.assertEqual(acquisitions[0].mode, "mosaic")
            self.assertEqual(
                [path.name for path in acquisitions[0].files],
                [
                    "Sample_F00.ims",
                    "Sample_F01.ims",
                    "Sample_F02.ims",
                    "Sample_F03.ims",
                ],
            )
            mosaic = load_mosaic_module(ims.ROOT)
            groups = mosaic.discover_groups(tmp / "Sample_F00.ims")
            self.assertEqual([tile.field for tile in groups["Sample"]], [0, 1, 2, 3])

    def test_direct_stitched_export_writes_four_stacks_and_projections(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            source = tmp / "StitchedSample.ims"
            create_synthetic_ims(source)
            acquisition = ims.Acquisition("StitchedSample", "stitched", (source,))
            outputs = ims.export_stitched_ims(acquisition, tmp / "preprocessing", True)
            self.assertEqual(set(outputs), set(ims.CHANNEL_ORDER))
            for index, channel in enumerate(ims.CHANNEL_ORDER):
                stack_path = (
                    tmp / "preprocessing/StitchedSample/01_exported" / f"{channel}.tif"
                )
                self.assertTrue(stack_path.is_file())
                self.assertEqual(tifffile.imread(stack_path).shape, (3, 5, 6))
                projection = tifffile.imread(outputs[channel])
                self.assertEqual(projection.shape, (5, 6))
                self.assertTrue(np.all(projection == (index + 1) * 100 + 2))

    def test_existing_mosaic_export_accepts_canonical_names(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            source = tmp / "Sample_F00.ims"
            create_synthetic_ims(source)
            mosaic = load_mosaic_module(ims.ROOT)
            names, calibration = mosaic.export_acquisition(
                [mosaic.Tile(source, 0)],
                tmp / "exported",
                overwrite=False,
                fallback_unit="um",
                channel_names=list(ims.CHANNEL_ORDER),
            )
            self.assertEqual(names, list(ims.CHANNEL_ORDER))
            self.assertEqual(calibration, (0.5, 0.6, 1.25))
            for channel in ims.CHANNEL_ORDER:
                self.assertTrue((tmp / "exported" / f"F00_{channel}.tif").is_file())

    def test_generated_downstream_config_uses_ims_xy_calibration(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            source = tmp / "Sample.ims"
            create_synthetic_ims(source)
            config_path = self._write_config(tmp, source)
            cfg = ims.read_config(config_path)
            acquisition = ims.discover_acquisitions(source, "auto")[0]
            metadata = ims.inspect_acquisition(acquisition, "um")
            outputs = ims.projection_paths(tmp / "preprocessing", acquisition)
            generated = ims.downstream_config(cfg, acquisition, metadata, outputs)
            self.assertEqual(generated["nuclear_channel"], "Hoechst")
            self.assertEqual(generated["pixel_size_um"], {"x": 0.5, "y": 0.6})
            self.assertEqual(list(generated["channels"]), list(ims.CHANNEL_ORDER))

    def test_dry_run_writes_no_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            source = tmp / "Sample.ims"
            create_synthetic_ims(source)
            config_path = self._write_config(tmp, source)
            self.assertEqual(ims.main(["--config", str(config_path), "--dry-run"]), 0)
            self.assertFalse((tmp / "output").exists())

    def test_mosaic_command_passes_fixed_channel_order(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            source = tmp / "Sample_F00.ims"
            create_synthetic_ims(source)
            config_path = self._write_config(tmp, source)
            raw = json.loads(config_path.read_text(encoding="utf-8"))
            raw["fiji_sif"] = str(tmp / "fiji.sif")
            config_path.write_text(json.dumps(raw), encoding="utf-8")
            cfg = ims.read_config(config_path)
            command = ims.mosaic_command(cfg, tmp / "preprocessing")
            index = command.index("--channel-names")
            self.assertEqual(command[index + 1 : index + 5], list(ims.CHANNEL_ORDER))
            self.assertIn("--output-root", command)

    def test_mosaic_requires_contiguous_field_numbers(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            for field in (0, 1, 2, 4):
                create_synthetic_ims(tmp / f"Sample_F{field:02d}.ims")
            config_path = self._write_config(tmp, tmp)
            with self.assertRaisesRegex(ValueError, "F00 through F03"):
                ims.main(["--config", str(config_path), "--dry-run"])

    def test_direct_ims_orchestration_builds_downstream_handoff(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            source = tmp / "Sample.ims"
            create_synthetic_ims(source)
            config_path = self._write_config(tmp, source)
            cfg = ims.read_config(config_path)
            acquisition = ims.discover_acquisitions(source, "auto")[0]
            metadata = {"Sample": ims.inspect_acquisition(acquisition, "um")}

            def fake_execute(state, label, command):
                self.assertEqual(label, "downstream analysis Sample")
                generated_path = Path(command[command.index("--config") + 1])
                generated = json.loads(generated_path.read_text(encoding="utf-8"))
                manifest = Path(generated["output_dir"]) / "run_manifest.json"
                manifest.parent.mkdir(parents=True, exist_ok=True)
                manifest.write_text("{}\n", encoding="utf-8")

            with (
                mock.patch.object(ims, "validate_runtime"),
                mock.patch.object(ims.IMSRun, "execute", new=fake_execute),
            ):
                ims.run(cfg, [acquisition], metadata, resume=False)

            generated_path = (
                tmp / "output/generated_configs/Sample.stitched_pipeline.json"
            )
            generated = json.loads(generated_path.read_text(encoding="utf-8"))
            self.assertEqual(generated["pixel_size_um"], {"x": 0.5, "y": 0.6})
            for channel in ims.CHANNEL_ORDER:
                self.assertTrue(Path(generated["channels"][channel]).is_file())
            metadata_path = tmp / "output/preprocessing/Sample/ims_metadata.json"
            self.assertTrue(metadata_path.is_file())
            top_manifest = json.loads(
                (tmp / "output/ims_run_manifest.json").read_text()
            )
            self.assertIn("analysis_Sample", top_manifest["completed"])

    def test_resume_refuses_changed_ims_input(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            source = tmp / "Sample.ims"
            create_synthetic_ims(source)
            cfg = ims.read_config(self._write_config(tmp, source))
            acquisitions = ims.discover_acquisitions(source, "auto")
            ims.IMSRun(cfg, acquisitions, resume=False)
            with h5py.File(source, "a") as handle:
                handle.attrs["Changed"] = "yes"
            with self.assertRaisesRegex(RuntimeError, "IMS inputs"):
                ims.IMSRun(cfg, acquisitions, resume=True)


if __name__ == "__main__":
    unittest.main()
