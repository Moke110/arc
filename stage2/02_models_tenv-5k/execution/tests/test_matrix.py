"""Fail-closed matrix planning/resume seams, without model imports."""
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
SCRIPT=Path(__file__).resolve().parents[1]/'run_matrix.py'
spec=importlib.util.spec_from_file_location('matrix231',SCRIPT)
matrix=importlib.util.module_from_spec(spec);spec.loader.exec_module(matrix)

class MatrixTests(unittest.TestCase):
    def test_plan_is_complete_and_unique(self):
        runs=matrix.plan()
        self.assertEqual(len(runs),90)
        self.assertEqual(len({tuple(r.values()) for r in runs}),90)
        self.assertEqual({r['seed'] for r in runs},{17,42,2026})
    def test_failed_smoke_missing_digest_resume_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary)
            for status,smoke in [('failed',False),('complete',True),('complete',False)]:
                (root/'terminal.json').write_text(json.dumps(dict(status=status,smoke=smoke)))
                with self.assertRaises(ValueError): matrix.verify_run(root,dict(adapter='mlp',scenario='cold-drug',fold=0,seed=17),'input')
    def test_matrix_requires_every_registered_identity(self):
        with self.assertRaises(ValueError): matrix.require_complete(plan=matrix.plan(),found=matrix.plan()[:-1])
        with self.assertRaises(ValueError): matrix.require_complete(plan=matrix.plan(),found=matrix.plan()[:-1]+matrix.plan()[:1])

if __name__=='__main__':unittest.main()
