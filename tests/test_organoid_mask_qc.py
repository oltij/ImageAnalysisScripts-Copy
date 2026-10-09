import csv
import importlib.util
import tempfile
import unittest
from pathlib import Path

import numpy as np
import tifffile


SRC = Path(__file__).resolve().parents[1] / "QC/organoid_mask_qc.py"
SPEC = importlib.util.spec_from_file_location("organoid_mask_qc", SRC)
qc = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(qc)


class OrganoidMaskQCTests(unittest.TestCase):
    def test_qc_uses_actual_masks_and_reports_outside_cell_pixels(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            fluorescence = np.arange(36, dtype=np.uint16).reshape(6, 6)
            organoid = np.zeros((6, 6), dtype=np.uint16)
            organoid[0:3, 0:3] = 1
            organoid[4:6, 4:6] = 2
            cells = np.zeros((6, 6), dtype=np.uint16)
            cells[1:3, 1:3] = 7
            cells[3, 3] = 8

            fluorescence_path = tmp / "fluorescence.tif"
            organoid_path = tmp / "OrganoidMask.tiff"
            reference_path = tmp / "HoechstOrganoidMask.tiff"
            cell_path = tmp / "CellMask.tiff"
            reference = np.zeros((6, 6), dtype=np.uint16)
            reference[0:3, 0:3] = 1
            tifffile.imwrite(fluorescence_path, fluorescence)
            tifffile.imwrite(organoid_path, organoid)
            tifffile.imwrite(reference_path, reference)
            tifffile.imwrite(cell_path, cells)

            figure_path, csv_path = qc.run_qc(
                [("Hoechst", fluorescence_path, organoid_path, cell_path)],
                tmp / "qc",
                reference_path,
            )

            self.assertTrue(figure_path.is_file())
            self.assertGreater(figure_path.stat().st_size, 0)
            with csv_path.open(newline="", encoding="utf-8") as handle:
                row = next(csv.DictReader(handle))
            self.assertEqual(int(row["organoid_mask_area_pixels"]), 13)
            self.assertEqual(int(row["connected_components"]), 2)
            self.assertEqual(int(row["largest_component_area_pixels"]), 9)
            self.assertEqual(row["touches_image_border"], "True")
            self.assertEqual(int(row["segmented_cell_pixels"]), 5)
            self.assertEqual(int(row["cell_pixels_outside_organoid"]), 1)
            self.assertAlmostEqual(
                float(row["cell_pixels_outside_organoid_percent"]), 20.0
            )
            self.assertEqual(int(row["hoechst_reference_mask_area_pixels"]), 9)
            self.assertEqual(int(row["organoid_mask_disagreement_pixels"]), 4)
            self.assertEqual(int(row["organoid_mask_pixels_outside_reference"]), 4)
            self.assertEqual(int(row["reference_mask_pixels_missing_from_organoid"]), 0)
            self.assertEqual(int(row["cell_pixels_outside_hoechst_reference"]), 1)
            self.assertAlmostEqual(
                float(row["cell_pixels_outside_hoechst_reference_percent"]),
                20.0,
            )

    def test_qc_rejects_shape_mismatch(self):
        with self.assertRaisesRegex(ValueError, "different shapes"):
            qc.mask_statistics(
                "PV",
                np.ones((3, 3), dtype=bool),
                np.ones((4, 3), dtype=bool),
                np.ones((3, 3), dtype=bool),
            )


if __name__ == "__main__":
    unittest.main()
