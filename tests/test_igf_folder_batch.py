import argparse
import contextlib
import importlib.util
import io
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import h5py


ROOT = Path(__file__).resolve().parents[1]
IGF_DIR = ROOT / "IGF Organoids"
CONFIG_PATH = IGF_DIR / "configs/ims_mge_generic.example.json"
PRODUCTION_PLAN_PATH = IGF_DIR / "configs/production_94_fields.json"


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


batch = load_module("run_folder_batch", IGF_DIR / "run_folder_batch.py")
submitter = load_module("submit_folder_batch", IGF_DIR / "submit_folder_batch.py")


class IGFFolderBatchTests(unittest.TestCase):
    def _write_ims(self, path, channel_names):
        path.parent.mkdir(parents=True, exist_ok=True)
        with h5py.File(path, "w") as handle:
            timepoint = handle.create_group("DataSet/ResolutionLevel 0/TimePoint 0")
            info = handle.create_group("DataSetInfo")
            for index, name in enumerate(channel_names):
                timepoint.create_group(f"Channel {index}")
                channel = info.create_group(f"Channel {index}")
                channel.attrs["Name"] = name

    def _named_field(
        self, root, session, group, field_number=1, tile_count=4,
        channel_names=None,
    ):
        field = root / session / "2026-09-28" / f"{group}_Field_{field_number}"
        paths = []
        for tile_id in range(tile_count):
            path = field / f"{group}_F{tile_id:02d}.ims"
            if channel_names is None:
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_bytes(f"source tile {tile_id}".encode())
            else:
                self._write_ims(path, channel_names)
            paths.append(path)
        return field, paths

    def _field(self, root, *, directory_id=7, filename_id=7, tile_ids=(2, 4, 6, 8)):
        session = root / "Olti_session" / "2026-09-28"
        field = session / f"Sample_DAPI_GFP_LHX6_PV_{directory_id}_Field_1"
        field.mkdir(parents=True)
        paths = []
        for tile_id in tile_ids:
            path = field / f"Sample_DAPI_GFP_LHX6_PV_{filename_id}_F{tile_id:02d}.ims"
            path.write_bytes(f"source tile {tile_id}".encode())
            paths.append(path)
        return field, paths

    def _prepare(self, source_root, workspace, **overrides):
        args = argparse.Namespace(
            command="prepare",
            input=str(source_root),
            group=None,
            field=None,
            max_depth=5,
            grid=[2, 2],
            require_contiguous=False,
            strict_incomplete=False,
            skip_incomplete=False,
            spec=str(CONFIG_PATH),
            workspace=str(workspace),
            repository=str(ROOT),
            session=None,
            run_plan=None,
            allow_dirty_repository=False,
        )
        for key, value in overrides.items():
            setattr(args, key, value)
        with contextlib.redirect_stdout(io.StringIO()):
            batch.do_prepare(args)
        return json.loads((workspace / "manifest.json").read_text())

    def test_recursive_prepare_preserves_sources_and_renumbers_staged_symlinks(self):
        with tempfile.TemporaryDirectory() as temp:
            temp = Path(temp)
            field, sources = self._field(temp / "source")
            before = {path: path.read_bytes() for path in sources}
            workspace = temp / "scratch_workspace"

            manifest = self._prepare(temp / "source", workspace)

            self.assertEqual(len(manifest["jobs"]), 1)
            job = manifest["jobs"][0]
            self.assertEqual(job["source_field"], str(field.resolve()))
            self.assertEqual(job["channel_names"], ["DAPI", "GFP", "LHX6", "PV"])
            self.assertEqual(job["nuclear_channel"], "DAPI")
            self.assertEqual(
                [tile["source_index"] for tile in job["tile_map"]],
                [2, 4, 6, 8],
            )
            self.assertEqual(
                [Path(tile["staged"]).name[-7:] for tile in job["tile_map"]],
                ["F00.ims", "F01.ims", "F02.ims", "F03.ims"],
            )
            for path in sources:
                self.assertEqual(path.read_bytes(), before[path])
            for tile in job["tile_map"]:
                staged = Path(tile["staged"])
                self.assertTrue(staged.is_symlink())
                self.assertEqual(staged.resolve(), Path(tile["source"]))
            self.assertEqual(Path(manifest["repository"]), ROOT)
            self.assertEqual(Path(manifest["batch_runner"]), batch.SCRIPT_PATH)
            self.assertEqual(Path(manifest["slurm_runner"]), batch.SBATCH_PATH)
            self.assertTrue((workspace / "controller_compat/run_ims_pipeline.py").is_file())

    def test_optional_acquisition_identifier_must_agree_between_directory_and_files(self):
        with tempfile.TemporaryDirectory() as temp:
            temp = Path(temp)
            field_path, _ = self._field(temp, directory_id=7, filename_id=8)
            fields = batch.discover_fields(
                temp, max_depth=5, tile_count=4, require_contiguous=False,
                group_pattern=None, field_numbers=set(),
            )
            self.assertEqual(fields[0].path, field_path.resolve())
            with self.assertRaisesRegex(
                ValueError,
                r"directory=7, IMS filenames=8",
            ):
                batch.validate_acquisition_identifiers(
                    fields[0], batch.DEFAULT_FILENAME_MARKERS
                )

    def test_explicit_channel_names_cannot_hide_acquisition_mismatch(self):
        with tempfile.TemporaryDirectory() as temp:
            temp = Path(temp)
            self._field(temp, directory_id=7, filename_id=8)
            field = batch.discover_fields(
                temp, max_depth=5, tile_count=4, require_contiguous=False,
                group_pattern=None, field_numbers=set(),
            )[0]
            spec = batch.load_spec(CONFIG_PATH, require_fiji=False)
            spec["channel_names"] = ["DAPI", "GFP", "LHX6", "PV"]
            with self.assertRaisesRegex(ValueError, "Acquisition identifier mismatch"):
                batch.resolve_field_channels(field, spec)

    def test_numeric_marker_suffix_is_not_treated_as_acquisition_identifier(self):
        markers = batch.DEFAULT_FILENAME_MARKERS
        self.assertIsNone(
            batch.optional_acquisition_identifier(
                "Sample_DAPI_NKX2_1", markers, filename=False
            )
        )
        self.assertEqual(
            batch.optional_acquisition_identifier(
                "Sample_DAPI_NKX2_1_7_F00.ims", markers, filename=True
            ),
            7,
        )

    def test_configuration_uses_exact_successful_alignment_arguments(self):
        config = json.loads(CONFIG_PATH.read_text())
        self.assertEqual(config["analysis_downsampling"], {
            "factor": 1, "method": "area_mean"
        })
        self.assertEqual(
            config["analysis"]["alignment"]["extra_args"],
            [
                "--skip-phase-init",
                "--match-gate-um",
                "25",
                "--iterations",
                "5",
                "--max-points",
                "1500",
            ],
        )

    def test_submit_command_passes_absolute_python_runner_to_spooled_sbatch(self):
        with tempfile.TemporaryDirectory() as temp:
            temp = Path(temp)
            self._field(temp / "source")
            workspace = temp / "scratch_workspace"
            self._prepare(temp / "source", workspace)
            output = io.StringIO()
            with contextlib.redirect_stdout(output):
                submitter.submit(
                    workspace / "manifest.json", 2, print_only=True
                )
            command = output.getvalue()
            self.assertIn(str(batch.SBATCH_PATH), command)
            self.assertIn(str(batch.SCRIPT_PATH), command)
            self.assertNotIn("/home/oltij/ImageAnalysisScripts-Copy", command)

    def test_slurm_profile_preserves_great_lakes_runtime(self):
        sbatch = (IGF_DIR / "run_folder_batch.sbatch").read_text()
        self.assertIn("#SBATCH --account=parent0", sbatch)
        self.assertIn("#SBATCH --partition=standard", sbatch)
        self.assertIn("#SBATCH --cpus-per-task=4", sbatch)
        self.assertIn("#SBATCH --mem=64G", sbatch)
        self.assertIn("#SBATCH --time=1-00:00:00", sbatch)
        self.assertNotIn("#SBATCH --gres", sbatch)
        self.assertIn("module load singularity/4.4.1", sbatch)
        self.assertIn('SINGULARITY_BINDPATH="/gpfs:/gpfs,/scratch:/scratch"', sbatch)
        self.assertIn("conda activate ims-mosaic", sbatch)
        self.assertIn('python "$BATCH_RUNNER"', sbatch)
        self.assertIn("--resume", sbatch)
        self.assertNotIn('SCRIPT_DIR="/home/', sbatch)

    def test_turbo_inputs_require_scratch_workspace(self):
        source = batch.TURBO_IMS_ROOT / "session"
        batch.validate_storage_locations(source, Path("/scratch/parent_root/run"))
        with self.assertRaisesRegex(ValueError, "require a workspace below /scratch"):
            batch.validate_storage_locations(source, Path("/home/user/run"))

    def test_scratch_symlink_is_compared_canonically(self):
        with tempfile.TemporaryDirectory() as temp:
            temp = Path(temp)
            turbo = temp / "nfs/turbo/umms-parent/andor_microscope_data_dump"
            turbo.mkdir(parents=True)
            canonical_scratch = temp / "gpfs/accounts"
            canonical_scratch.mkdir(parents=True)
            scratch_alias = temp / "scratch"
            scratch_alias.symlink_to(canonical_scratch, target_is_directory=True)
            with mock.patch.object(batch, "TURBO_IMS_ROOT", turbo), mock.patch.object(
                batch, "SCRATCH_ROOT", scratch_alias
            ):
                batch.validate_storage_locations(
                    turbo / "Olti_session",
                    scratch_alias / "parent_root/parent0/oltij/run",
                )
                with self.assertRaisesRegex(ValueError, "require a workspace below"):
                    batch.validate_storage_locations(
                        turbo / "Olti_session", temp / "home/oltij/run"
                    )

    def test_repeatable_sessions_limit_recursive_discovery_and_order(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            self._named_field(root, "Olti_B", "Sample_DAPI_GFP")
            self._named_field(root, "Olti_A", "Sample_DAPI_GFP")
            self._named_field(root, "Olti_unrelated", "Sample_DAPI_GFP")
            fields = batch.discover_fields(
                root, max_depth=5, tile_count=4, require_contiguous=False,
                group_pattern=None, field_numbers=set(),
                sessions=["Olti_B", "Olti_A"],
            )
            self.assertEqual([field.session for field in fields], ["Olti_A", "Olti_B"])
            self.assertNotIn("Olti_unrelated", {field.session for field in fields})
            self.assertEqual([field.path for field in fields], sorted(
                [field.path for field in fields], key=lambda path: str(path).casefold()
            ))

    def test_session_selectors_reject_missing_duplicate_and_traversal(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root / "Olti_A").mkdir()
            for sessions, message in (
                (["Olti_missing"], "does not exist"),
                (["Olti_A", "olti_a"], "Duplicate"),
                (["../Olti_A"], "Invalid"),
                (["/absolute"], "Invalid"),
            ):
                with self.subTest(sessions=sessions), self.assertRaisesRegex(
                    ValueError, message
                ):
                    batch.selected_session_roots(root, sessions)

    def test_selected_session_incomplete_fields_are_reported_only_from_selection(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            self._named_field(root, "Olti_A", "Sample_DAPI_GFP", tile_count=3)
            self._named_field(root, "Olti_B", "Sample_DAPI_GFP")
            issues = []
            fields = batch.discover_fields(
                root, max_depth=5, tile_count=4, require_contiguous=False,
                group_pattern=None, field_numbers=set(), sessions=["Olti_A"],
                issues=issues, scan_only=True,
            )
            self.assertEqual(fields, [])
            self.assertEqual(len(issues), 1)
            self.assertEqual(issues[0]["session"], "Olti_A")
            self.assertNotIn("Olti_B", issues[0]["path"])

    def test_set_identifier_is_metadata_and_five_channel_layout_is_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            group = "Exp17_157-2_Set3_DAPI_GFP_LHX6_PV"
            field_path, paths = self._named_field(
                root, "Olti_06_18_2026", group,
                channel_names=[
                    "White Light - Phase", "Confocal - Blue", "Confocal - Green",
                    "Confocal - Red", "Confocal - Far Red",
                ],
            )
            field = batch.Field(
                field_path.resolve(), group, 1,
                tuple(batch.Tile(i, path.resolve()) for i, path in enumerate(paths)),
                "Olti_06_18_2026",
            )
            self.assertEqual(
                batch.channel_suffix(paths[0].name, batch.DEFAULT_FILENAME_MARKERS),
                ["DAPI", "GFP", "LHX6", "PV"],
            )
            self.assertEqual(batch.validate_experimental_set_identifiers(field), 3)
            with self.assertRaisesRegex(ValueError, "unsupported physical phase"):
                batch.detect_field_channels(field, batch.DEFAULT_FILENAME_MARKERS)

    def test_filename_markers_never_absorb_set_identifier_for_five_channels(self):
        filename = "Exp17_157-2_Set3_DAPI_GFP_LHX6_PV_F00.ims"
        with self.assertRaisesRegex(ValueError, "5 physical channels"):
            batch.channels_from_filename(
                filename, batch.DEFAULT_FILENAME_MARKERS, count=5
            )

    def test_valid_two_three_four_channel_suffixes_and_multitoken_marker(self):
        cases = (
            ("Sample_DAPI_GFP_F00.ims", 2, ["DAPI", "GFP"]),
            ("Sample_DAPI_GFP_PV_F00.ims", 3, ["DAPI", "GFP", "PV"]),
            ("Sample_DAPI_GFP_LHX6_PV_F00.ims", 4, ["DAPI", "GFP", "LHX6", "PV"]),
            ("Sample_DAPI_NKX2_1_F00.ims", 2, ["DAPI", "NKX2_1"]),
        )
        for filename, count, expected in cases:
            with self.subTest(filename=filename):
                self.assertEqual(
                    batch.channels_from_filename(
                        filename, batch.DEFAULT_FILENAME_MARKERS, count
                    ),
                    expected,
                )

    def test_submit_interface_forwards_repeatable_sessions(self):
        with tempfile.TemporaryDirectory() as temp, mock.patch.object(
            submitter.batch, "do_prepare"
        ) as prepare, mock.patch.object(submitter, "submit"):
            result = submitter.main([
                "--input", str(Path(temp) / "input"),
                "--spec", str(CONFIG_PATH),
                "--workspace", str(Path(temp) / "workspace"),
                "--session", "Olti_A", "--session", "Olti_B",
                "--print-only",
            ])
            self.assertEqual(result, 0)
            forwarded = prepare.call_args.args[0]
            self.assertEqual(forwarded.session, ["Olti_A", "Olti_B"])

    def test_submit_help_formats_literal_percent(self):
        result = subprocess.run(
            [sys.executable, str(IGF_DIR / "submit_folder_batch.py"), "--help"],
            text=True, capture_output=True, check=True,
        )
        help_text = result.stdout
        self.assertIn("Slurm %", help_text)
        self.assertIn("throttle", help_text)

    def test_committed_production_plan_has_exact_audited_scope(self):
        plan = batch.load_run_plan(PRODUCTION_PLAN_PATH)
        self.assertEqual(len(plan["sessions"]), 12)
        self.assertEqual(sum(x["complete_fields"] for x in plan["sessions"]), 94)
        self.assertEqual(plan["expected"]["original_tiles"], 376)
        self.assertEqual(plan["expected"]["incomplete_fields"], 19)
        self.assertNotIn("Olti_07_16_2026", [x["name"] for x in plan["sessions"]])
        self.assertEqual(
            plan["expected"]["biological_channels"],
            ["DAPI", "GFP", "LHX6", "PV"],
        )
        self.assertEqual(
            plan["configuration"]["sha256"], batch.sha256_file(CONFIG_PATH)
        )

    def test_production_git_state_requires_commit_and_clean_checkout(self):
        with self.assertRaisesRegex(ValueError, "requires a Git checkout"):
            batch.require_production_git_state(
                {"commit": None, "dirty": None}, allow_dirty=True
            )
        with self.assertRaisesRegex(ValueError, "refuses a dirty"):
            batch.require_production_git_state(
                {"commit": "a" * 40, "dirty": True}, allow_dirty=False
            )
        batch.require_production_git_state(
            {"commit": "a" * 40, "dirty": True}, allow_dirty=True
        )

    def test_run_plan_builds_and_validates_deterministic_manifest(self):
        with tempfile.TemporaryDirectory() as temp:
            temp = Path(temp)
            source = temp / "source"
            _, first_tiles = self._named_field(
                source, "Olti_A", "Sample_Set1_DAPI_GFP",
                channel_names=["Confocal - Blue", "Confocal - Green"],
            )
            self._named_field(
                source, "Olti_A", "Incomplete_Set1_DAPI_GFP", field_number=2,
                tile_count=3,
                channel_names=["Confocal - Blue", "Confocal - Green"],
            )
            _, second_tiles = self._named_field(
                source, "Olti_B", "Sample_Set2_DAPI_GFP",
                channel_names=["Confocal - Blue", "Confocal - Green"],
            )
            plan_path = temp / "small_plan.json"
            plan_path.write_text(json.dumps({
                "schema": 1,
                "name": "synthetic production audit",
                "input_root": str(source.resolve()),
                "sessions": [
                    {"name": "Olti_A", "complete_fields": 1},
                    {"name": "Olti_B", "complete_fields": 1},
                ],
                "expected": {
                    "complete_fields": 2,
                    "incomplete_fields": 1,
                    "original_tiles": 8,
                    "physical_channels_per_tile": 2,
                    "biological_channels": ["DAPI", "GFP"],
                },
                "configuration": {
                    "filename": CONFIG_PATH.name,
                    "sha256": batch.sha256_file(CONFIG_PATH),
                },
            }))
            source_bytes = {
                path: path.read_bytes() for path in first_tiles + second_tiles
            }
            clean_git = {"commit": "a" * 40, "dirty": False, "status": []}
            workspace = temp / "workspace"
            with mock.patch.object(
                batch, "repository_git_state", return_value=clean_git
            ):
                manifest = self._prepare(
                    source, workspace, run_plan=str(plan_path),
                    allow_dirty_repository=False,
                )
                summary = batch.validate_manifest_for_submission(manifest)

            self.assertEqual(summary["fields"], 2)
            self.assertEqual(summary["tiles"], 8)
            self.assertEqual([job["index"] for job in manifest["jobs"]], [0, 1])
            self.assertEqual(
                [job["session"] for job in manifest["jobs"]],
                ["Olti_A", "Olti_B"],
            )
            self.assertEqual(manifest["selected_sessions"], ["Olti_A", "Olti_B"])
            self.assertTrue(manifest["run_plan"]["validation"]["passed"])
            self.assertEqual(manifest["run_plan"]["validation"]["incomplete_fields"], 1)
            for job in manifest["jobs"]:
                self.assertEqual(job["physical_channel_count"], 2)
                self.assertEqual(job["physical_channel_indices"], [0, 1])
                self.assertEqual(
                    job["physical_channel_names"],
                    ["Confocal - Blue", "Confocal - Green"],
                )
                self.assertIn(job["experimental_set_number"], (1, 2))
            for path, before in source_bytes.items():
                self.assertEqual(path.read_bytes(), before)


if __name__ == "__main__":
    unittest.main()
