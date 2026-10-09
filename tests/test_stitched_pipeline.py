import importlib.util
import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np
import tifffile


SRC = Path(__file__).resolve().parents[1] / "run_stitched_pipeline.py"
SPEC = importlib.util.spec_from_file_location("stitched_runner", SRC)
runner = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(runner)


class StitchedPipelineTests(unittest.TestCase):
    def test_runner_resolves_repository_root(self):
        self.assertEqual(runner.ROOT, SRC.parent)
        self.assertTrue((runner.ROOT / "Alignment/CASTalign_two_channel_registration.py").is_file())

    def test_replace_top_level_assignments_only(self):
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / "script.py"
            source.write_text(
                'from pathlib import Path\nA = Path(\n  "old"\n)\nB = {"x": 1}\n'
                'def f():\n    A = "keep"\n',
                encoding="utf-8",
            )
            changed = runner.patched_source(
                source, {"A": runner.path_value("/new"), "B": repr({"z": 3})}
            ).decode()
            self.assertIn("A = Path('/new')", changed)
            self.assertIn("B = {'z': 3}", changed)
            self.assertIn('A = "keep"', changed)

    def test_missing_patch_assignment_fails(self):
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / "source.py"
            source.write_text("A=1\n", encoding="utf-8")
            with self.assertRaises(ValueError):
                runner.patched_source(source, {"B": "2"})

    def test_every_runner_patch_contract_matches_real_scripts(self):
        path = runner.path_value("/tmp/example")
        channels = repr([{"name": "A", "roi_csv": "/tmp/a.csv", "original": "/tmp/a.tif"}])
        contracts = {
            "ExtractingROIs/extract_rois.py": {
                "tiff_path": path,
                "channel_csv_path": path,
                "organoid_csv_path": path,
                "CHANNEL_NAME": repr("A"),
                "output_dir": path,
            },
            "ROIFiltering/Intensity/intensityfilter.py": {
                "output_dir": path,
                "channels": channels,
            },
            "ROIFiltering/FinalFilter/filter.py": {
                "output_dir": path,
                "channels": channels,
                "FILTERS": repr(runner.filters_for_pass({}, "shape")),
                "PIXEL_SIZE_X_UM": "1.0",
                "PIXEL_SIZE_Y_UM": "1.0",
            },
            "Colocalization/colocalizationdapiscript1.py": {
                "OBJECT_A_NAME": repr("A"),
                "OBJECT_B_NAME": repr("B"),
                "OBJECT_A_CSV": path,
                "OBJECT_B_CSV": path,
                "OBJECT_A_TIFF": path,
                "OBJECT_B_TIFF": path,
                "ORGANOID_CSV": path,
                "RESULTS_FILE": path,
                "THRESHOLD": "0.25",
                "PIXEL_SIZE_X_UM": "1.0",
                "PIXEL_SIZE_Y_UM": "1.0",
            },
            "Colocalization/colocalizationdapiscript2.py": {
                "RESULTS_FILE": path,
                "OUTPUT_DIR": path,
                "ENABLE_ORGANOID_AREA_SUMMARY": "True",
                "ORGANOID_MEASUREMENTS_CSV": path,
                "ENABLE_NUCLEAR_ALL_CELL_SUMMARY": "False",
                "OBJECT_A_NUCLEAR_RESULTS_FILE": "None",
                "OBJECT_B_NUCLEAR_RESULTS_FILE": "None",
            },
        }
        for relative in runner.SHAPE_REPORTS.values():
            contracts[relative] = {
                "output_dir": path,
                "channels": channels,
                "PIXEL_SIZE_X_UM": "1.0",
                "PIXEL_SIZE_Y_UM": "1.0",
            }
        for relative, replacements in contracts.items():
            with self.subTest(script=relative):
                patched = runner.patched_source(runner.ROOT / relative, replacements)
                compile(patched, relative, "exec")

    def test_filter_passes_remain_separate(self):
        preferences = {
            "intensity": {"Mean intensity": {"percentile": 15, "keep": "above"}},
            "shape": {"Area_um2": {"percentile": 20, "keep": "above"}},
        }
        intensity = runner.filters_for_pass(preferences, "intensity")
        shape = runner.filters_for_pass(preferences, "shape")
        self.assertEqual(intensity["Mean intensity"]["percentile"], 15)
        self.assertEqual(intensity["Area_um2"]["percentile"], 0)
        self.assertEqual(shape["Mean intensity"]["percentile"], 0)
        self.assertEqual(shape["Area_um2"]["percentile"], 20)

    def test_wrong_filter_metric_and_boolean_percentile_fail(self):
        with self.assertRaises(ValueError):
            runner.filters_for_pass(
                {"intensity": {"Area_um2": {"percentile": 10, "keep": "above"}}},
                "intensity",
            )
        with self.assertRaises(ValueError):
            runner.filters_for_pass(
                {"intensity": {"Mean intensity": {"percentile": True, "keep": "above"}}},
                "intensity",
            )

    def test_checked_in_cellprofiler_pipeline_is_text_and_compatible(self):
        pipeline = runner.ROOT / "CellProfiler/MGEOPVFinal.cppipe"
        text = runner.validate_cellprofiler_pipeline(pipeline)
        self.assertTrue(text.startswith("CellProfiler Pipeline:"))
        self.assertIn("ModuleCount:24", text)

    def test_binary_project_renamed_cppipe_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            fake = Path(tmp) / "bad.cppipe"
            fake.write_bytes(b"\x89HDF\r\n\x1a\nmore")
            with self.assertRaisesRegex(ValueError, "HDF5 CellProfiler project"):
                runner.validate_cellprofiler_pipeline(fake)

    def _write_config(self, directory, channels, output="output", pairs=None):
        config = {
            "nuclear_channel": "DNA",
            "channels": channels,
            "output_dir": output,
            "pixel_size_um": {"x": 1.0, "y": 1.0},
            "cellprofiler": {
                "pipeline": str(runner.ROOT / "CellProfiler/MGEOPVFinal.cppipe")
            },
            "alignment": {"enabled": False, "extra_args": []},
            "thresholds": {"marker_nuclear": 0.25, "marker_marker": 0.25},
            "filters": {"intensity": {}, "shape": {}},
            "reports": {"intensity": False, "shape": False, "colocalization": False},
        }
        if pairs is not None:
            config["pairs"] = pairs
        path = Path(directory) / "config.json"
        path.write_text(json.dumps(config), encoding="utf-8")
        return path

    def test_duplicate_directional_pair_fails_config(self):
        with tempfile.TemporaryDirectory() as tmp:
            config = self._write_config(
                tmp,
                {"DNA": "dna.tif", "A": "a.tif", "B": "b.tif"},
                pairs=[["A", "B"], ["A", "B"]],
            )
            with self.assertRaisesRegex(ValueError, "duplicate"):
                runner.read_config(config)

    def test_tiff_preflight_and_projection(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            dna = tmp / "dna.tif"
            marker = tmp / "marker.tif"
            stack = np.arange(3 * 4 * 5, dtype=np.uint16).reshape(3, 4, 5)
            tifffile.imwrite(dna, np.ones((4, 5), dtype=np.uint16), metadata={"axes": "YX"})
            tifffile.imwrite(
                marker, stack, metadata={"axes": "ZYX"}, photometric="minisblack"
            )
            config = self._write_config(tmp, {"DNA": str(dna), "Marker": str(marker)}, pairs=[])
            cfg = runner.read_config(config)
            descriptions = runner.validate_inputs(cfg)
            self.assertFalse(descriptions["DNA"]["projection_required"])
            self.assertTrue(descriptions["Marker"]["projection_required"])
            projected_path = tmp / "projection.tif"
            result = runner.ensure_projection(marker, projected_path)
            np.testing.assert_array_equal(tifffile.imread(result), stack.max(axis=0))
            self.assertEqual(runner.inspect_tiff(result)["axes"], "YX")

    def test_mismatched_xy_dimensions_fail_preflight(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            dna = tmp / "dna.tif"
            marker = tmp / "marker.tif"
            tifffile.imwrite(dna, np.ones((4, 5), dtype=np.uint8), metadata={"axes": "YX"})
            tifffile.imwrite(marker, np.ones((5, 5), dtype=np.uint8), metadata={"axes": "YX"})
            config = self._write_config(tmp, {"DNA": str(dna), "Marker": str(marker)}, pairs=[])
            with self.assertRaisesRegex(ValueError, "one XY pixel canvas"):
                runner.validate_inputs(runner.read_config(config))

    def test_resume_refuses_changed_pipeline(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            dna = tmp / "dna.tif"
            tifffile.imwrite(dna, np.ones((2, 2), dtype=np.uint8), metadata={"axes": "YX"})
            pipeline = tmp / "pipeline.cppipe"
            pipeline.write_text(
                (runner.ROOT / "CellProfiler/MGEOPVFinal.cppipe").read_text(encoding="utf-8"),
                encoding="utf-8",
            )
            config = json.loads(
                self._write_config(tmp, {"DNA": str(dna)}, output="analysis", pairs=[]).read_text()
            )
            config["cellprofiler"]["pipeline"] = str(pipeline)
            config_path = tmp / "config.json"
            config_path.write_text(json.dumps(config), encoding="utf-8")
            cfg = runner.read_config(config_path)
            runner.Runner(cfg, resume=False)
            pipeline.write_text(pipeline.read_text(encoding="utf-8") + "\n# changed\n", encoding="utf-8")
            cfg = runner.read_config(config_path)
            with self.assertRaisesRegex(RuntimeError, "analysis scripts changed"):
                runner.Runner(cfg, resume=True)

    def test_full_orchestration_routes_each_handoff_in_order(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            channels = {}
            for name in ("DNA", "LHX6", "PV"):
                image = tmp / f"{name}.tif"
                tifffile.imwrite(
                    image, np.ones((4, 5), dtype=np.uint8), metadata={"axes": "YX"}
                )
                channels[name] = str(image)
            config_path = self._write_config(
                tmp, channels, output="analysis", pairs=[["LHX6", "PV"]]
            )
            raw = json.loads(config_path.read_text(encoding="utf-8"))
            raw["alignment"]["enabled"] = True
            raw["reports"] = {
                "organoid_mask": True,
                "intensity": True,
                "shape": True,
                "colocalization": True,
            }
            config_path.write_text(json.dumps(raw), encoding="utf-8")
            cfg = runner.read_config(config_path)

            class FakeRunner:
                def __init__(self):
                    self.out = cfg["_output"]
                    self.out.mkdir()
                    self.manifest = self.out / "run_manifest.json"
                    self.completed = {}
                    self.calls = []

                @staticmethod
                def _create(paths):
                    for path in paths:
                        path.parent.mkdir(parents=True, exist_ok=True)
                        path.write_bytes(b"test-output")

                def cp(self, label, image, out):
                    self.calls.append(("cp", label, image, out))
                    expected = [
                        out / "OrganoidMask.tiff",
                        out / "CellMask.tiff",
                        out / "MyExpt_FilterObjects2.csv",
                        out / "MyExpt_FilterObjects.csv",
                    ]
                    self._create(expected)
                    self.completed[label] = [str(path) for path in expected]

                def execute(self, label, command, expected):
                    self.calls.append(("execute", label, command, expected))
                    self._create(expected)
                    self.completed[label] = [str(path) for path in expected]

                def patched(self, label, script, replacements, expected):
                    self.calls.append(("patched", label, script, replacements, expected))
                    self._create(expected)
                    self.completed[label] = [str(path) for path in expected]

                def save(self):
                    pass

            fake = FakeRunner()

            def fake_export(_h5path, csv_path):
                csv_path.write_text("ROI,X,Y\n1,1,1\n", encoding="utf-8")

            with mock.patch.object(runner, "Runner", return_value=fake), mock.patch.object(
                runner, "export_nuclear_positive", side_effect=fake_export
            ):
                runner.run(cfg, resume=False)

            labels = [call[1] for call in fake.calls]
            expected_labels = [
                "02_initial_cp_DNA",
                "02_initial_cp_LHX6",
                "02_initial_cp_PV",
                "03_align_LHX6",
                "03_align_PV",
                "04_aligned_cp_DNA",
                "04_aligned_cp_LHX6",
                "04_aligned_cp_PV",
                "04_organoid_mask_qc",
                "04_roi_DNA",
                "04_roi_LHX6",
                "04_roi_PV",
                "05_intensity_report",
                "05_intensity_gate",
                "06_shape_report_area",
                "06_shape_report_circularity",
                "06_shape_report_eccentricity",
                "06_shape_report_solidity",
                "06_shape_gate",
                "07_match_LHX6",
                "07_match_LHX6_viz",
                "07_match_PV",
                "07_match_PV_viz",
                "08_pair_LHX6_PV",
                "08_pair_LHX6_PV_viz",
            ]
            self.assertEqual(labels, expected_labels)

            calls_by_label = {call[1]: call for call in fake.calls}
            shape_gate = calls_by_label["06_shape_gate"]
            shape_channels = shape_gate[3]["channels"]
            self.assertIn("05_intensity_filtered/LHX6/LHX6_ROI_pixels_FILTERED.csv", shape_channels)
            nuclear_match = calls_by_label["07_match_LHX6"]
            self.assertIn("06_final_shape_filtered/LHX6", nuclear_match[3]["OBJECT_A_CSV"])
            marker_pair = calls_by_label["08_pair_LHX6_PV"]
            self.assertIn("LHX6_nuclear_positive_cells.csv", marker_pair[3]["OBJECT_A_CSV"])
            self.assertIn("PV_nuclear_positive_cells.csv", marker_pair[3]["OBJECT_B_CSV"])


if __name__ == "__main__":
    unittest.main()
