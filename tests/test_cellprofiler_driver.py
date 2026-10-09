import importlib.util
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock


SRC = Path(__file__).resolve().parents[1] / "CellProfiler/cellprofilerdriver.py"
SPEC = importlib.util.spec_from_file_location("cellprofiler_driver", SRC)
driver = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(driver)


class CellProfilerDriverTests(unittest.TestCase):
    def test_optional_mask_is_linked_as_a_second_named_input(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            pipeline = tmp / "pipeline.cppipe"
            fluorescence = tmp / "marker.tif"
            mask = tmp / "OrganoidMask.tiff"
            output = tmp / "output"
            pipeline.write_text("pipeline", encoding="utf-8")
            fluorescence.write_bytes(b"fluorescence")
            mask.write_bytes(b"mask")

            def fake_run(command):
                input_dir = Path(command[command.index("-i") + 1])
                fluorescence_link = input_dir / "__analysis_fluorescence__.tif"
                mask_link = input_dir / "__hoechst_organoid_mask__.tif"
                self.assertTrue(fluorescence_link.is_symlink())
                self.assertTrue(mask_link.is_symlink())
                self.assertEqual(fluorescence_link.resolve(), fluorescence.resolve())
                self.assertEqual(mask_link.resolve(), mask.resolve())
                return SimpleNamespace(returncode=0)

            with mock.patch.object(driver.subprocess, "run", side_effect=fake_run):
                driver.run_cellprofiler(
                    pipeline,
                    fluorescence,
                    output,
                    organoid_mask=mask,
                )


if __name__ == "__main__":
    unittest.main()
