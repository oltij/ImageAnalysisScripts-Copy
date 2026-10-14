import importlib.util
import json
import os
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
LAUNCHER_PATH = ROOT / "IGF Organoids/run_igf_organoids.py"
EXAMPLE_PATH = ROOT / "IGF Organoids/configs/igf_stitched.example.json"
SPEC = importlib.util.spec_from_file_location("igf_launcher", LAUNCHER_PATH)
launcher = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(launcher)


class IGFWorkflowTests(unittest.TestCase):
    def test_launcher_resolves_root_and_forwards_supported_flags_from_any_cwd(self):
        with tempfile.TemporaryDirectory() as tmp:
            config = Path(tmp) / "igf.json"
            raw = json.loads(EXAMPLE_PATH.read_text(encoding="utf-8"))
            config.write_text(json.dumps(raw), encoding="utf-8")
            previous = Path.cwd()
            try:
                os.chdir(tmp)
                command = launcher.build_command(
                    Path("igf.json"), dry_run=True, resume=True
                )
            finally:
                os.chdir(previous)

        self.assertEqual(launcher.REPOSITORY_ROOT, ROOT)
        self.assertEqual(launcher.PIPELINE_RUNNER, ROOT / "run_stitched_pipeline.py")
        self.assertEqual(command[0], sys.executable)
        self.assertEqual(command[1], str(ROOT / "run_stitched_pipeline.py"))
        self.assertIn("--dry-run", command)
        self.assertIn("--resume", command)
        self.assertEqual(
            Path(command[command.index("--config") + 1]), config.resolve()
        )

    def test_launcher_executes_shared_runner_with_same_python(self):
        with tempfile.TemporaryDirectory() as tmp:
            config = Path(tmp) / "igf.json"
            config.write_text(EXAMPLE_PATH.read_text(encoding="utf-8"), encoding="utf-8")
            with mock.patch.object(launcher.os, "execv") as execute:
                launcher.main(["--config", str(config), "--dry-run"])
            execute.assert_called_once()
            executable, command = execute.call_args.args
            self.assertEqual(executable, sys.executable)
            self.assertEqual(command, launcher.build_command(config, dry_run=True))

    def test_example_has_exact_validated_igf_registration_profile(self):
        raw = json.loads(EXAMPLE_PATH.read_text(encoding="utf-8"))
        self.assertEqual(raw["nuclear_channel"], "DAPI")
        self.assertEqual(list(raw["channels"]), ["DAPI", "GFP", "LHX6", "PV"])
        self.assertEqual(
            raw["pairs"],
            [["GFP", "LHX6"], ["GFP", "PV"], ["LHX6", "PV"]],
        )
        self.assertTrue(raw["alignment"]["enabled"])
        self.assertEqual(
            raw["alignment"]["extra_args"], list(launcher.IGF_ALIGNMENT_ARGS)
        )
        self.assertNotIn("--try-affine", raw["alignment"]["extra_args"])

    def test_launcher_rejects_changed_registration_profile(self):
        with tempfile.TemporaryDirectory() as tmp:
            config = Path(tmp) / "igf.json"
            raw = json.loads(EXAMPLE_PATH.read_text(encoding="utf-8"))
            raw["alignment"]["extra_args"].remove("--skip-phase-init")
            config.write_text(json.dumps(raw), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "must be exactly"):
                launcher.validate_igf_alignment_profile(config)


if __name__ == "__main__":
    unittest.main()
