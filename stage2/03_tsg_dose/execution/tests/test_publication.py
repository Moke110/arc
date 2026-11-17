"""Release acceptance and filesystem installation boundaries, small fixtures only."""
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
import hashlib
from unittest.mock import patch

STEP=Path(__file__).resolve().parents[2]


def publication():
    spec=importlib.util.spec_from_file_location('dose_publication',STEP/'execution/publish_release.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


def complete_matrix(root):
    methods=['no-dose','scaling','input-numeric','latent-numeric','input-embedding','latent-embedding','film']
    runs=[dict(scenario=s,fold=f,method=m,seed=t) for s in ['cold-drug','held-dose','final'] for f in range(1 if s=='final' else 3) for m in methods for t in [17,42,2026]]
    (root/'plan.json').write_text(json.dumps(dict(runs=runs)))
    plan_digest=hashlib.sha256((root/'plan.json').read_bytes()).hexdigest()
    accepted=dict(schema='tsg-dose-matrix-acceptance-v1',scenarios=['cold-drug','held-dose','final'],n_formal_fits=147,runs=runs,plan_sha256=plan_digest)
    (root/'acceptance-cold-drug-held-dose-final.json').write_text(json.dumps(accepted))
    for run in runs:
        path=root/'runs'/run['scenario']/run['method']/str(run['fold'])/str(run['seed']);path.mkdir(parents=True)
        (path/'terminal.json').write_text(json.dumps(dict(status='validation-complete' if run['scenario']=='final' else 'complete')))
        digest=hashlib.sha256((path/'terminal.json').read_bytes()).hexdigest()
        (path/'manifest.json').write_text(json.dumps(dict(schema='tsg-dose-fit-v1',outputs={'terminal.json':digest})))
    return runs


class PublicationTests(unittest.TestCase):
    def test_incomplete_matrix_cannot_become_a_release_candidate(self):
        module=publication()
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            (root/'acceptance-cold-drug-held-dose-final.json').write_text(json.dumps(dict(n_formal_fits=146)))
            with self.assertRaisesRegex(ValueError,'147'):module.verify_matrix(root)

    def test_complete_unique_origins_are_required_and_failed_attempt_blocks_publication(self):
        module=publication()
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);runs=complete_matrix(root)
            self.assertEqual(module.verify_matrix(root)['n_formal_fits'],147)
            path=root/'attempts/interrupted-origin';path.mkdir(parents=True)
            (path/'terminal.json').write_text(json.dumps(dict(status='failed',error='nonfinite')))
            with self.assertRaisesRegex(ValueError,'failed'):module.verify_matrix(root)

    def test_corrupted_origin_payload_and_duplicate_coverage_are_rejected(self):
        module=publication()
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);runs=complete_matrix(root)
            path=root/'runs/cold-drug/no-dose/0/17/terminal.json';original=path.read_bytes()
            path.write_text(json.dumps(dict(status='complete',changed=True)))
            with self.assertRaisesRegex(ValueError,'digest'):module.verify_matrix(root)
            path.write_bytes(original)
            receipt=root/'acceptance-cold-drug-held-dose-final.json';accepted=json.loads(receipt.read_text())
            accepted['runs'][-1]=accepted['runs'][0];receipt.write_text(json.dumps(accepted))
            with self.assertRaisesRegex(ValueError,'coverage'):module.verify_matrix(root)

    def test_failed_install_restores_old_accepted_bytes_and_preserves_candidate(self):
        module=publication()
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);candidate=root/'candidate';stable=root/'model'
            candidate.mkdir();stable.mkdir()
            (candidate/'manifest.json').write_text('new accepted manifest')
            (stable/'manifest.json').write_text('old accepted manifest')
            receipt=dict(status='passed',n_formal_fits=147,downstream_exact_replay=True,csv_only_figures=True,candidate_manifest_sha256=hashlib.sha256((candidate/'manifest.json').read_bytes()).hexdigest())
            original_replace=module.os.replace
            def fail_new_install(source,destination):
                if Path(source)==candidate:raise OSError('injected filesystem replacement failure')
                return original_replace(source,destination)
            with patch.object(module.os,'replace',side_effect=fail_new_install):
                with self.assertRaisesRegex(OSError,'injected'):module.install(candidate,stable,receipt)
            self.assertEqual((stable/'manifest.json').read_text(),'old accepted manifest')
            self.assertEqual((candidate/'manifest.json').read_text(),'new accepted manifest')
            module.install(candidate,stable,receipt)
            self.assertEqual((stable/'manifest.json').read_text(),'new accepted manifest')
            self.assertFalse(candidate.exists())

    def test_install_requires_successful_downstream_and_csv_only_acceptance(self):
        module=publication()
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);candidate=root/'candidate';candidate.mkdir()
            (candidate/'manifest.json').write_text('manifest')
            with self.assertRaisesRegex(ValueError,'acceptance'):module.install(candidate,root/'model',dict(status='pending'))

    def test_preservation_snapshot_detects_protected_byte_changes(self):
        module=publication()
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);(root/'accepted.json').write_text('prior accepted bytes')
            snapshot=module.file_snapshot(root,['accepted.json'])
            module.verify_file_snapshot(root,snapshot)
            (root/'accepted.json').write_text('changed bytes')
            with self.assertRaisesRegex(ValueError,'preservation'):module.verify_file_snapshot(root,snapshot)


if __name__=='__main__':unittest.main()
