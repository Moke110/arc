"""Public policy CLI tests use worked native-iteration examples."""
import json
from pathlib import Path
import subprocess
import sys
import unittest

CLI=Path(__file__).resolve().parents[1]/'script/cmonge_v2_protocol.py'
class PolicyTest(unittest.TestCase):
    def call(self, payload):
        return subprocess.run([sys.executable,str(CLI)],input=json.dumps(payload),text=True,capture_output=True)
    def test_plateau_earliest_1000_and_strict_best_independent(self):
        events=[{'iteration':i,'loss':100.-i/100000} for i in range(50,1001,50)]
        result=self.call({'events':events,'prior_training_seconds':13975.509756658023})
        self.assertEqual(result.returncode,0,result.stderr)
        value=json.loads(result.stdout)
        self.assertEqual(value['stop_iteration'],1000)
        self.assertEqual(value['best_iteration'],1000)
        self.assertAlmostEqual(value['remaining_seconds'],424.490243341977)
    def test_missing_debit_resume_and_exhaustion_rejected(self):
        for payload in ({'events':[]},{'events':[],'prior_training_seconds':14400},
                        {'events':[],'prior_training_seconds':0,'resume':False}):
            self.assertNotEqual(self.call(payload).returncode,0)
    def test_non_cadence_rejected_and_relative_improvement_resets(self):
        self.assertNotEqual(self.call({'events':[{'iteration':51,'loss':1}],
                                     'prior_training_seconds':0}).returncode,0)
        events=[{'iteration':i,'loss':100 if i<950 else 99} for i in range(50,1401,50)]
        result=self.call({'events':events,'prior_training_seconds':0})
        self.assertEqual(result.returncode,0,result.stderr)
        self.assertIsNone(json.loads(result.stdout)['stop_iteration'])
if __name__=='__main__':unittest.main()
