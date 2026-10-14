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


ROOT = Path(__file__).resolve().parents[1]
IGF_DIR = ROOT / "IGF Organoids"


def load_module(name, path):
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


batch = load_module("run_folder_batch", IGF_DIR / "run_folder_batch.py")
production = load_module(
    "submit_production_batches", IGF_DIR / "submit_production_batches.py"
)


class IGFProductionSubmissionTests(unittest.TestCase):
    def _context(self, workspace):
        manifest_path = workspace / "manifest.json"
        manifest_path.write_text("{}")
        manifest = {
            "workspace": str(workspace),
            "jobs": [{"index": index} for index in range(94)],
            "run_plan": {"validation": {"passed": True}},
            "slurm_runner": str(IGF_DIR / "run_folder_batch.sbatch"),
            "batch_runner": str(IGF_DIR / "run_folder_batch.py"),
        }
        return manifest_path, manifest

    def _validation_mocks(self, manifest):
        return (
            mock.patch.object(production.batch, "load_manifest", return_value=manifest),
            mock.patch.object(
                production.batch,
                "validate_manifest_for_submission",
                return_value={"fields": 94, "tiles": 376},
            ),
        )

    def test_print_only_shows_exact_arrays_and_dependency_without_submission(self):
        with tempfile.TemporaryDirectory() as temp:
            workspace = Path(temp)
            manifest_path, manifest = self._context(workspace)
            load_mock, validate_mock = self._validation_mocks(manifest)
            output = io.StringIO()
            with load_mock, validate_mock, mock.patch.object(
                production.subprocess, "run"
            ) as run, contextlib.redirect_stdout(output):
                result = production.submit_two_batches(
                    manifest_path, print_only=True
                )
            self.assertIsNone(result)
            run.assert_not_called()
            text = output.getvalue()
            self.assertIn("--array=0-46%47", text)
            self.assertIn("--array=47-93%47", text)
            self.assertIn("--dependency=afterok:<BATCH1_JOB_ID>", text)
            self.assertFalse((workspace / production.RECORD_NAME).exists())

    def test_submission_uses_afterok_and_records_both_job_ids(self):
        with tempfile.TemporaryDirectory() as temp:
            workspace = Path(temp)
            manifest_path, manifest = self._context(workspace)
            load_mock, validate_mock = self._validation_mocks(manifest)
            results = [
                subprocess.CompletedProcess([], 0, stdout="111;cluster\n", stderr=""),
                subprocess.CompletedProcess([], 0, stdout="222\n", stderr=""),
            ]
            with load_mock, validate_mock, mock.patch.object(
                production.shutil, "which", return_value="/usr/bin/sbatch"
            ), mock.patch.object(
                production.subprocess, "run", side_effect=results
            ) as run, contextlib.redirect_stdout(io.StringIO()):
                ids = production.submit_two_batches(manifest_path)

            self.assertEqual(ids, ("111", "222"))
            first = run.call_args_list[0].args[0]
            second = run.call_args_list[1].args[0]
            self.assertIn("--array=0-46%47", first)
            self.assertNotIn("--dependency=afterok:111", first)
            self.assertIn("--array=47-93%47", second)
            self.assertIn("--dependency=afterok:111", second)
            self.assertTrue(Path(first[-3]).is_absolute())
            self.assertTrue(Path(first[-2]).is_absolute())
            self.assertTrue(Path(first[-1]).is_absolute())
            record = json.loads((workspace / production.RECORD_NAME).read_text())
            self.assertEqual(record["status"], "both_submitted")
            self.assertEqual(record["batch1_job_id"], "111")
            self.assertEqual(record["batch2_job_id"], "222")
            self.assertEqual(record["batch2_dependency"], "afterok:111")

    def test_duplicate_submission_is_rejected_before_sbatch(self):
        with tempfile.TemporaryDirectory() as temp:
            workspace = Path(temp)
            manifest_path, manifest = self._context(workspace)
            (workspace / production.RECORD_NAME).write_text(json.dumps({
                "batch1_job_id": "111", "batch2_job_id": "222"
            }))
            load_mock, validate_mock = self._validation_mocks(manifest)
            with load_mock, validate_mock, mock.patch.object(
                production.subprocess, "run"
            ) as run, self.assertRaisesRegex(ValueError, "already recorded"):
                production.submit_two_batches(manifest_path)
            run.assert_not_called()

    def test_partial_submission_can_resume_only_second_batch(self):
        with tempfile.TemporaryDirectory() as temp:
            workspace = Path(temp)
            manifest_path, manifest = self._context(workspace)
            (workspace / production.RECORD_NAME).write_text(json.dumps({
                "status": "batch1_submitted",
                "batch1_job_id": "111",
                "batch2_job_id": None,
                "batch1_command": ["sbatch", "--array=0-46%47"],
                "submitted_unix_time": 1,
            }))
            load_mock, validate_mock = self._validation_mocks(manifest)
            result = subprocess.CompletedProcess([], 0, stdout="222\n", stderr="")
            with load_mock, validate_mock, mock.patch.object(
                production.shutil, "which", return_value="/usr/bin/sbatch"
            ), mock.patch.object(
                production.subprocess, "run", return_value=result
            ) as run, contextlib.redirect_stdout(io.StringIO()):
                ids = production.submit_two_batches(
                    manifest_path, resume_second_batch=True
                )
            self.assertEqual(ids, ("111", "222"))
            self.assertEqual(run.call_count, 1)
            self.assertIn("--dependency=afterok:111", run.call_args.args[0])


if __name__ == "__main__":
    unittest.main()
