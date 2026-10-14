import importlib.util
import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest import mock

import numpy as np
import tifffile


SRC = (
    Path(__file__).resolve().parents[1]
    / "Alignment/CASTalign_two_channel_registration.py"
)

# The local test environment does not need CASTalign for the canvas helpers.
# Supply a placeholder module so the script can be imported; production runtime
# validation still requires the real package before registration starts.
ADDED_FAKE_CASTALIGN = "castalign" not in sys.modules
if ADDED_FAKE_CASTALIGN:
    sys.modules["castalign"] = types.ModuleType("castalign")

SPEC = importlib.util.spec_from_file_location("castalign_registration", SRC)
registration = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(registration)
if ADDED_FAKE_CASTALIGN:
    del sys.modules["castalign"]


class CASTalignCanvasTests(unittest.TestCase):
    def test_equal_sized_phase_inputs_are_not_copied_or_padded(self):
        fixed = np.arange(20, dtype=np.uint16).reshape(4, 5)
        moving = fixed + 1
        fixed_result, moving_result = registration.phase_correlation_canvases(
            fixed, moving
        )
        self.assertIs(fixed_result, fixed)
        self.assertIs(moving_result, moving)

    def test_unequal_phase_inputs_are_temporarily_padded_without_mutation(self):
        fixed = np.arange(20, dtype=np.uint16).reshape(4, 5)
        moving = (100 + np.arange(35, dtype=np.uint16)).reshape(5, 7)
        fixed_before = fixed.copy()
        moving_before = moving.copy()

        with mock.patch.object(
            registration,
            "phase_cross_correlation",
            return_value=(np.array([1.0, -2.0]), 0.1, 0.2),
        ) as phase:
            shift, error, difference = registration.phase_correlation_shift(
                fixed, moving, upsample_factor=7
            )

        padded_fixed, padded_moving = phase.call_args.args[:2]
        self.assertEqual(padded_fixed.shape, (5, 7))
        self.assertEqual(padded_moving.shape, (5, 7))
        np.testing.assert_array_equal(padded_fixed[:4, :5], fixed)
        np.testing.assert_array_equal(padded_moving, moving)
        self.assertTrue(np.all(padded_fixed[4:, :] == 0))
        self.assertTrue(np.all(padded_fixed[:, 5:] == 0))
        np.testing.assert_array_equal(fixed, fixed_before)
        np.testing.assert_array_equal(moving, moving_before)
        np.testing.assert_array_equal(shift, [1.0, -2.0])
        self.assertEqual((error, difference), (0.1, 0.2))
        self.assertEqual(phase.call_args.kwargs["upsample_factor"], 7)

    def test_skip_phase_initializer_never_calls_phase_correlation(self):
        fixed = np.ones((4, 5), dtype=np.uint8)
        moving = np.ones((5, 7), dtype=np.uint8)
        with mock.patch.object(registration, "phase_correlation_shift") as phase:
            shift, error, difference = registration.phase_initializer(
                fixed, moving, skip=True, upsample_factor=10
            )
        phase.assert_not_called()
        np.testing.assert_array_equal(shift, [0.0, 0.0])
        self.assertTrue(np.isnan(error))
        self.assertTrue(np.isnan(difference))

    def test_transform_writes_fixed_canvas_with_one_resampling_call(self):
        moving = np.arange(35, dtype=np.uint16).reshape(5, 7)
        moving_before = moving.copy()

        class FakeTransform:
            def __init__(self):
                self.calls = []

            def transform_image(self, image, **kwargs):
                self.calls.append((image, kwargs))
                return np.ones(kwargs["output_size"], dtype=np.float32)

        transform = FakeTransform()
        registered = registration.transform_moving_to_fixed_canvas(
            transform, moving, (4, 5)
        )
        self.assertEqual(registered.shape, (4, 5))
        self.assertEqual(len(transform.calls), 1)
        self.assertIs(transform.calls[0][0], moving)
        self.assertEqual(transform.calls[0][1]["output_size"], (1, 4, 5))
        self.assertTrue(transform.calls[0][1]["force_size"])
        self.assertFalse(transform.calls[0][1]["labels"])
        np.testing.assert_array_equal(moving, moving_before)

    def test_canvas_helpers_do_not_modify_original_tiff_files(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            fixed_path = tmp / "dapi.tif"
            moving_path = tmp / "pv.tif"
            tifffile.imwrite(fixed_path, np.ones((4, 5), dtype=np.uint16))
            tifffile.imwrite(moving_path, np.ones((5, 7), dtype=np.uint16) * 2)
            fixed_bytes = fixed_path.read_bytes()
            moving_bytes = moving_path.read_bytes()

            fixed = registration.load_tiff_2d(fixed_path)
            moving = registration.load_tiff_2d(moving_path)
            registration.phase_initializer(
                fixed, moving, skip=True, upsample_factor=10
            )

            self.assertEqual(fixed_path.read_bytes(), fixed_bytes)
            self.assertEqual(moving_path.read_bytes(), moving_bytes)


if __name__ == "__main__":
    unittest.main()
