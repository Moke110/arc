import sys,json,tempfile,hashlib
from pathlib import Path
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from run_matrix import plan,require_complete,lock_selection,verify_test_sequence,write,event,execution_registration,verify_run

class Matrix(unittest.TestCase):
    def test_complete_registered_task_graph(self):
        rows=plan()
        self.assertEqual(len(rows),126)
        self.assertEqual([len(plan(s)) for s in ['cold-drug','held-dose']],[63,63])
        require_complete(plan=rows,found=rows)
        with self.assertRaises(ValueError):require_complete(plan=rows,found=rows[:-1])
        with self.assertRaises(ValueError):require_complete(plan=rows,found=rows+[rows[0]])

    def test_locked_choice_reads_only_validation_and_rejects_modified_scores(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);input_digest='input'
            for run in plan('cold-drug'):
                if run['fold']!=0:continue
                path=root/'runs'/run['scenario']/run['method']/str(run['fold'])/str(run['seed']);path.mkdir(parents=True)
                config=dict(run,smoke_rows=None)
                (path/'config.json').write_text(json.dumps(config))
                (path/'terminal.json').write_text(json.dumps(dict(status='validation-complete',parameter_count=10)))
                (path/'validation_metrics.json').write_text(json.dumps(dict(pcc=.5 if run['method']=='scaling' else .1)))
                (path/'test_prediction.npy').write_text('unreadable outer data')
                hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in path.iterdir() if p.name not in ['terminal.json','test_prediction.npy']}
                (path/'validation_manifest.json').write_text(json.dumps(dict(schema='tsg-dose-validation-v1',input_manifest_sha256=input_digest,parameter_count=10,outputs=hashes)))
            receipt=lock_selection(root,'cold-drug',0,input_digest)
            self.assertEqual(json.loads(receipt.read_text())['method'],'scaling')
            terminals=[root/'runs/cold-drug/no-dose/0'/str(seed)/'terminal.json' for seed in [17,42,2026]]
            for terminal in terminals:terminal.write_text(json.dumps(dict(status='validation-complete',parameter_count=1)))
            with self.assertRaisesRegex(ValueError,'parameter count changed'):
                lock_selection(root,'cold-drug',0,input_digest)
            for terminal in terminals:terminal.write_text(json.dumps(dict(status='validation-complete',parameter_count=10)))
            changed=root/'runs/cold-drug/scaling/0/17/validation_metrics.json'
            changed.write_text('{"pcc":0}')
            with self.assertRaises(ValueError):lock_selection(root,'cold-drug',0,input_digest)

    def test_outer_access_requires_matching_prior_choice(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);run=plan('cold-drug')[0]
            receipt=root/'selections/cold-drug/0/receipt.json';write(receipt,dict(method='no-dose'))
            digest=hashlib.sha256(receipt.read_bytes()).hexdigest()
            event(root,'test-access-start',run=run,receipt_sha256=digest)
            with self.assertRaises(ValueError):verify_test_sequence(root,'cold-drug',0,[run])
            (root/'test-access-sequence.jsonl').unlink()
            event(root,'choice-locked',scenario='cold-drug',fold=0,receipt_sha256=digest)
            event(root,'test-access-start',run=run,receipt_sha256=digest)
            path=root/'runs/cold-drug/no-dose/0/17'
            write(path/'selection_receipt.json',dict(receipt=dict(method='no-dose'),sha256=digest))
            verify_test_sequence(root,'cold-drug',0,[run])
            receipt.write_text('{"method":"film"}')
            with self.assertRaises(ValueError):verify_test_sequence(root,'cold-drug',0,[run])

    def test_acceptance_rejects_a_foreign_formal_origin(self):
        run=plan('cold-drug')[0]
        frozen=dict(device='cuda',input_manifest_sha256='one')
        registration=execution_registration(run,frozen,'a'*64)
        with tempfile.TemporaryDirectory() as tmp:
            changed=dict(frozen,input_manifest_sha256='other')
            with self.assertRaisesRegex(ValueError,'registration changed'):
                verify_run(tmp,run,changed,registration)

if __name__=='__main__':unittest.main()
