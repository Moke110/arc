"""Approved plan/execution/publication boundaries for optimization round1."""
import importlib.util
import subprocess
import sys
import unittest
from pathlib import Path

STEP=Path(__file__).resolve().parents[2]
SPEC=importlib.util.spec_from_file_location('round1_matrix',STEP/'execution/run_matrix.py')
MATRIX=importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MATRIX)

class Round1(unittest.TestCase):
    def test_round1_plan_contains_only_126_outer_test_origins(self):
        rows=MATRIX.plan()
        self.assertEqual(len(rows),126)
        self.assertEqual({r['scenario'] for r in rows},{'cold-drug','held-dose'})
        self.assertEqual(len(MATRIX.plan('cold-drug')),63)
        self.assertEqual(len(MATRIX.plan('held-dose')),63)
        with self.assertRaisesRegex(ValueError,'final|round1'):
            MATRIX.plan('final')

    def test_cli_rejects_final_before_creating_candidate(self):
        for entry in ['run_matrix.py','transigen.py']:
            result=subprocess.run([sys.executable,str(STEP/'execution'/entry),*(['run'] if entry=='run_matrix.py' else []),'--inputs','/missing','--candidate','/missing','--source','/missing','--scenario','final'],capture_output=True,text=True)
            self.assertNotEqual(result.returncode,0)
            self.assertIn('final',result.stderr)

    def test_deployment_publication_is_retired(self):
        result=subprocess.run([sys.executable,str(STEP/'execution/publish_release.py'),'publish'],capture_output=True,text=True)
        self.assertNotEqual(result.returncode,0)
        self.assertIn('retired',result.stderr.lower())

    def test_frozen_execution_identities_resolve_to_preserved_bytes(self):
        import hashlib
        expected={'execution/run_matrix.py':'f979849ca66f33f0eb208178082d6bfe9061875c867e64cb2c4408c9384938ce',
                  'execution/transigen.py':'690324d3e38050e5f3d9b2bf18bd450179eb85c3eaf68a06357f76c3e60cbfe8'}
        for original,preserved in MATRIX.historical_code_mapping().items():
            self.assertEqual(hashlib.sha256((STEP/preserved).read_bytes()).hexdigest(),expected[original])

if __name__=='__main__':unittest.main()
