import json
import tempfile
import unittest
from pathlib import Path
import importlib.util

SRC = Path(__file__).resolve().parents[1] / 'run_stitched_pipeline.py'
spec = importlib.util.spec_from_file_location('runner', SRC)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

class TestPatching(unittest.TestCase):
    def test_replace_top_level_only(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)/'script.py'
            p.write_text('from pathlib import Path\nA = Path(\n  "old"\n)\nB = {"x": 1}\ndef f():\n    A = "keep"\n')
            changed=mod.patched_source(p, {'A':mod.path_value('/new'), 'B': repr({'z':3})}).decode()
            self.assertIn('A = Path(\'/new\')', changed)
            self.assertIn('B = {\'z\': 3}',changed)
            self.assertIn('A = "keep"',changed)
            self.assertIn('def f():', changed)
    def test_missing_assignment(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'source.py';p.write_text('A=1\n')
            with self.assertRaises(ValueError): mod.patched_source(p,{'B':'2'})
    def test_separate_filter_phases(self):
        prefs={'intensity':{'Mean intensity':{'percentile':15,'keep':'above'}},'shape':{'Area_um2':{'percentile':20,'keep':'above'}}}
        i=mod.filters_for_pass(prefs,'intensity');s=mod.filters_for_pass(prefs,'shape')
        self.assertEqual(i['Mean intensity']['percentile'],15)
        self.assertEqual(i['Area_um2']['percentile'],0)
        self.assertEqual(s['Mean intensity']['percentile'],0)
        self.assertEqual(s['Area_um2']['percentile'],20)
    def test_wrong_metric_fails(self):
        with self.assertRaises(ValueError):mod.filters_for_pass({'intensity':{'Area_um2':{'percentile':10,'keep':'above'}}},'intensity')

if __name__=='__main__': unittest.main()
