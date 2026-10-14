import importlib.util
import sys
import types
import unittest
from pathlib import Path
from unittest import mock

import numpy as np


SRC = (
    Path(__file__).resolve().parents[1]
    / "Alignment/CASTalign_two_channel_registration.py"
)

# These tests exercise the image-canvas compatibility path, not CASTalign
# itself. Production dependency checks still require the real package.
added_fake_castalign = "castalign" not in sys.modules
if added_fake_castalign:
    sys.modules["castalign"] = types.ModuleType("castalign")
spec = importlib.util.spec_from_file_location("castalign_registration", SRC)
registration = importlib.util.module_from_spec(spec)
spec.loader.exec_module(registration)
if added_fake_castalign:
    del sys.modules["castalign"]


class CASTalignCanvasTests(unittest.TestCase):
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

    def test_equal_phase_inputs_retain_original_behavior(self):
        fixed = np.arange(20, dtype=np.uint16).reshape(4, 5)
        moving = fixed + 1
        with mock.patch.object(
            registration,
            "phase_cross_correlation",
            return_value=(np.array([0.0, 0.0]), 0.0, 0.0),
        ) as phase:
            registration.phase_correlation_shift(fixed, moving)
        phase_fixed, phase_moving = phase.call_args.args[:2]
        np.testing.assert_array_equal(phase_fixed, fixed.astype(np.float32))
        np.testing.assert_array_equal(phase_moving, moving.astype(np.float32))
        self.assertEqual(phase_fixed.shape, fixed.shape)
        self.assertEqual(phase_moving.shape, moving.shape)


if __name__ == "__main__":
    unittest.main()
