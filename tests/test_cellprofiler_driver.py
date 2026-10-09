import importlib.util
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

import numpy as np
import tifffile


SRC = Path(__file__).resolve().parents[1] / "CellProfiler/cellprofilerdriver.py"
SPEC = importlib.util.spec_from_file_location("cellprofiler_driver", SRC)
driver = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(driver)


class CellProfilerDriverTests(unittest.TestCase):
    def test_optional_mask_is_normalized_as_uint8_0_255_second_input(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            pipeline = tmp / "pipeline.cppipe"
            fluorescence = tmp / "marker.tif"
            mask = tmp / "OrganoidMask.tiff"
            output = tmp / "output"
            pipeline.write_text("pipeline", encoding="utf-8")
            fluorescence.write_bytes(b"fluorescence")
            source_mask = np.zeros((5, 6), dtype=np.uint16)
            source_mask[1:4, 2:5] = 1
            tifffile.imwrite(mask, source_mask)

            def fake_run(command):
                input_dir = Path(command[command.index("-i") + 1])
                fluorescence_link = input_dir / "__analysis_fluorescence__.tif"
                mask_link = input_dir / "__hoechst_organoid_mask__.tif"
                self.assertTrue(fluorescence_link.is_symlink())
                self.assertTrue(mask_link.is_file())
                self.assertFalse(mask_link.is_symlink())
                self.assertEqual(fluorescence_link.resolve(), fluorescence.resolve())
                normalized = tifffile.imread(mask_link)
                self.assertEqual(normalized.dtype, np.uint8)
                self.assertEqual(set(np.unique(normalized)), {0, 255})
                np.testing.assert_array_equal(normalized > 0, source_mask > 0)
                return SimpleNamespace(returncode=0)

            with mock.patch.object(driver.subprocess, "run", side_effect=fake_run):
                driver.run_cellprofiler(
                    pipeline,
                    fluorescence,
                    output,
                    organoid_mask=mask,
                )

            source_after = tifffile.imread(mask)
            self.assertEqual(source_after.dtype, np.uint16)
            np.testing.assert_array_equal(source_after, source_mask)

    def test_empty_optional_mask_is_rejected_before_cellprofiler(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            pipeline = tmp / "pipeline.cppipe"
            fluorescence = tmp / "marker.tif"
            mask = tmp / "OrganoidMask.tiff"
            pipeline.write_text("pipeline", encoding="utf-8")
            fluorescence.write_bytes(b"fluorescence")
            tifffile.imwrite(mask, np.zeros((5, 6), dtype=np.uint16))

            with (
                mock.patch.object(driver.subprocess, "run") as run,
                self.assertRaisesRegex(ValueError, "no foreground pixels"),
            ):
                driver.run_cellprofiler(
                    pipeline,
                    fluorescence,
                    tmp / "output",
                    organoid_mask=mask,
                )
            run.assert_not_called()


if __name__ == "__main__":
    unittest.main()
