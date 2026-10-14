import argparse
import contextlib
import importlib.util
import io
import json
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
IGF_DIR = ROOT / "IGF Organoids"
CONFIG_PATH = IGF_DIR / "configs/ims_mge_generic.example.json"


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


batch = load_module("run_folder_batch", IGF_DIR / "run_folder_batch.py")
submitter = load_module("submit_folder_batch", IGF_DIR / "submit_folder_batch.py")


class IGFFolderBatchTests(unittest.TestCase):
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

    def _prepare(self, source_root, workspace):
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
        )
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


if __name__ == "__main__":
    unittest.main()
