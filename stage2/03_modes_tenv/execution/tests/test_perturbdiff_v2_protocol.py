"""Offline artifact CLI policy evidence; no native model imports."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

SCRIPT=Path(__file__).resolve().parents[1]/'script/perturbdiff_v2_protocol.py'
class ProtocolCLI(unittest.TestCase):
    def run_policy(self,events):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'history.json'; p.write_text(json.dumps(events))
            return subprocess.run([sys.executable,str(SCRIPT),'--history',str(p)],capture_output=True,text=True)
    def test_warmup_and_ten_subsequent_checks_best_is_independent(self):
        events=[{'step':s,'validation_loss':1-s/10000000} for s in range(1000,31000,1000)]
        result=self.run_policy(events)
        self.assertEqual(result.returncode,0,result.stderr)
        audit=json.loads(result.stdout)
        self.assertEqual(audit['stop_step'],30000)
        self.assertEqual(audit['best_step'],30000)
        self.assertEqual(audit['counter'],10)
    def test_legacy_or_missing_debit_request_is_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'request.json'
            for request in ({},{'prior_training_seconds':0,'resume_checkpoint':'v1.pt'},{'prior_training_seconds':14400}):
                p.write_text(json.dumps(request))
                result=subprocess.run([sys.executable,str(SCRIPT),'--request',str(p)],capture_output=True,text=True)
                self.assertNotEqual(result.returncode,0)
            p.write_text(json.dumps({'prior_training_seconds':12475.943569269002}))
            result=subprocess.run([sys.executable,str(SCRIPT),'--request',str(p)],capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stderr)
            self.assertAlmostEqual(json.loads(result.stdout)['remaining_training_seconds'],1924.056430730998)

    def test_partial_validation_batch_has_group_weight(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'batches.json'; p.write_text(json.dumps([{'loss':1,'groups':8},{'loss':9,'groups':2}]))
            result=subprocess.run([sys.executable,str(SCRIPT),'--batches',str(p)],capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stderr)
            self.assertEqual(json.loads(result.stdout)['validation_loss'],2.6)

    def test_noncheckpoint_history_is_rejected(self):
        self.assertNotEqual(self.run_policy([{'step':98,'validation_loss':1}]).returncode,0)

if __name__=='__main__': unittest.main()
