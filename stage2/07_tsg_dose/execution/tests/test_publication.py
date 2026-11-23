"""Release acceptance and filesystem installation boundaries, small fixtures only."""
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
import hashlib
import csv
import sys
from unittest.mock import patch

STEP=Path(__file__).resolve().parents[2]


def publication():
    spec=importlib.util.spec_from_file_location('dose_publication',STEP/'execution/historical/publish_release_v1.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    module.STEP=STEP
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
    def test_publication_cli_finds_checkout_marker_above_relocated_step(self):
        module=publication()
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);(root/'.git').write_text('gitdir: elsewhere\n')
            step=root/'arbitrary'/'extra'/'stage2'/'07_tsg_dose';step.mkdir(parents=True)
            argv=['publish_release.py','publish']
            for name in ['inputs','matrix','candidate','source','log','destination']:
                argv.extend(['--'+name,str(root/name)])
            with patch.object(module,'STEP',step),patch.object(sys,'argv',argv),patch.object(module,'publish') as publish:
                module.main()
            self.assertEqual(publish.call_args.args[0].root,root)

    def test_held_dose_intervals_resample_parents_carrying_all_seed_metrics(self):
        module=publication();rows=[]
        for parent,scores in [('A',[.1,.3,.5]),('B',[.5,.7,.9])]:
            for seed,value in zip([17,42,2026],scores):
                rows.append(dict(method='film',dose_um=1,drug_parent=parent,seed=seed,pcc=value,spearman=value,rmse=value,r2='' if parent=='A' else value))
        intervals=module.held_dose_uncertainty(rows)
        pcc=next(row for row in intervals if row['metric']=='pcc')
        self.assertAlmostEqual(pcc['estimate'],.5)
        self.assertAlmostEqual(pcc['low'],.3);self.assertAlmostEqual(pcc['high'],.7)
        self.assertEqual((pcc['n_parents'],pcc['n_seeds'],pcc['repetitions']),(2,3,10000))
        r2=next(row for row in intervals if row['metric']=='r2')
        self.assertAlmostEqual(r2['low'],.7);self.assertAlmostEqual(r2['high'],.7)
        self.assertEqual((r2['n_defined_parents'],r2['n_excluded_undefined']),(1,1))
        self.assertEqual(intervals,module.held_dose_uncertainty(list(reversed(rows))))
        with self.assertRaisesRegex(ValueError,'seed'):
            module.held_dose_uncertainty(rows[:-1])

    def test_export_keeps_curve_coordinates_and_auxiliary_selection_references(self):
        module=publication()
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);matrix=root/'matrix';matrix.mkdir()
            for scenario in ['cold-drug','held-dose','final']:
                report=matrix/'reports'/scenario;report.mkdir(parents=True)
                (report/'learning_curves.csv').write_text('scenario,method,fold,seed,epoch,train_loss,validation_mse\n'+scenario+',film,0,17,1,2.5,1.25\n')
                (report/'fit_resources.csv').write_text('scenario,method,fold,seed,parameter_count,elapsed_seconds,peak_rss_kib\n'+scenario+',film,0,17,123,4.5,1000\n')
                if scenario!='final':
                    (report/'summary.csv').write_text('method,pcc,spearman,rmse,r2,n_parents,undefined_pcc,undefined_r2\nfilm,0.4,0.3,0.2,-0.1,1,0,0\n')
                    (report/'per_parent.csv').write_text('method,drug_parent,n_seeds,pcc,spearman,rmse,r2,truth_constant_count,prediction_constant_count\nfilm,A,3,0.4,0.3,0.2,-0.1,0,0\n')
                    (report/'per_dose.csv').write_text('method,dose_um,pcc\nfilm,1.0,0.3\n')
                for fold in range(1 if scenario=='final' else 3):
                    selection=matrix/'selections'/scenario/str(fold);selection.mkdir(parents=True)
                    (selection/'receipt.json').write_text(json.dumps(dict(method='film' if scenario=='final' else 'scaling')))
                    (selection/'validation_table.csv').write_text('method,seed,pcc\nfilm,17,0.5\n')
            (matrix/'reports/cold-drug/primary_paired.csv').write_text('metric,estimate,low,high\npcc,0.1,-0.1,0.2\n')
            (matrix/'reports/cold-drug/descriptive_paired.csv').write_text('metric,estimate,low,high\npcc,0.1,-0.1,0.2\n')
            (matrix/'reports/held-dose/per_parent_dose_seed.csv').write_text('method,dose_um,drug_parent,seed,pcc,spearman,rmse,r2\n'+''.join(f'film,1.0,A,{seed},0.3,0.2,1.0,0.1\n' for seed in [17,42,2026]))
            destination=root/'data';module.export_plot_data(matrix,destination)
            rows=module.csv_rows(destination/'learning_curves_all.csv')
            self.assertEqual(len(rows),3);self.assertEqual(rows[2]['validation_mse'],'1.25');self.assertEqual(rows[2]['train_loss'],'2.5')
            held=module.csv_rows(destination/'held_dose_values.csv')[0]
            self.assertEqual(held['final_selected_method'],'film');self.assertEqual(held['primary_fold_2_selected_method'],'scaling')
            self.assertEqual(held['evaluation'],'interpolation')
            self.assertEqual(float(held['pcc_low']),.3);self.assertEqual(float(held['pcc_high']),.3)
            self.assertEqual(len(module.csv_rows(destination/'held_dose_uncertainty.csv')),4)
            historical_step=root/'historical-step';(historical_step/'figure/script').mkdir(parents=True)
            import shutil
            shutil.copy2(STEP/'execution/historical/render_v1.py',historical_step/'figure/script/render.py')
            with patch.object(module,'STEP',historical_step):
                rendered=module.render_isolated(destination,root/'figures',sys.executable)
            self.assertEqual(len(rendered),12)
            metrics=module.csv_rows(destination/'metrics_values.csv')
            self.assertEqual(len(metrics),8)
            self.assertEqual({r['metric'] for r in metrics},{'pcc','spearman','rmse','r2'})
            self.assertEqual(next(r['value'] for r in metrics if r['scenario']=='cold-drug' and r['metric']=='r2'),'-0.1')
            self.assertEqual({r['scenario'] for r in metrics},{'cold-drug','held-dose'})
            self.assertIn('learning_curve.png',rendered)
            self.assertIn('metrics.png',rendered)
            self.assertTrue(all((root/'figures'/name).read_bytes().startswith(b'\x89PNG') for name in rendered))

    def test_assembly_publishes_only_selected_three_checkpoint_origins(self):
        module=publication()
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);matrix=root/'matrix';matrix.mkdir();complete_matrix(matrix)
            inputs=root/'inputs';inputs.mkdir()
            (inputs/'manifest.json').write_text('{}');input_digest=module.digest(inputs/'manifest.json')
            (inputs/'config.json').write_text(json.dumps(dict(axis_size=5000)))
            (inputs/'genes.csv').write_text('gene_index\nG1\n')
            selection=matrix/'selections/final/0';selection.mkdir(parents=True)
            (selection/'receipt.json').write_text(json.dumps(dict(method='film',input_manifest_sha256=input_digest)))
            for seed in [17,42,2026]:
                path=matrix/f'runs/final/film/0/{seed}'
                (path/'checkpoint.pt').write_bytes(f'checkpoint-{seed}'.encode())
                (path/'environment.json').write_text(json.dumps(dict(device='cuda')))
                (path/'terminal.json').write_text(json.dumps(dict(status='validation-complete',parameter_count=123)))
                (path/'manifest.json').write_text(json.dumps(dict(schema='tsg-dose-fit-v1',outputs={p.name:module.digest(p) for p in path.iterdir() if p.name!='manifest.json'})))
            candidate=root/'release';receipt=module.assemble(matrix,inputs,candidate)
            self.assertEqual(receipt['method'],'film')
            self.assertEqual(set(receipt['outputs']),{'genes.csv','input_contract.json','selection_receipt.json','checkpoints/17.pt','checkpoints/42.pt','checkpoints/2026.pt'})
            self.assertEqual((candidate/'checkpoints/42.pt').read_bytes(),b'checkpoint-42')
            self.assertEqual(module.read(candidate/'input_contract.json')['input_manifest_sha256'],input_digest)
            with self.assertRaisesRegex(ValueError,'empty'):module.assemble(matrix,inputs,candidate)

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

    def test_failed_figure_install_restores_accepted_model_csv_and_png(self):
        module=publication()
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);candidate=root/'candidate';stable=root/'model'
            staged=root/'staged-figure';figures=root/'figure'
            for path in [candidate,stable,staged,figures]:path.mkdir()
            (candidate/'manifest.json').write_text('new model')
            (stable/'manifest.json').write_text('old model')
            for path,label in [(staged,'new'),(figures,'old')]:
                (path/'data').mkdir();(path/'data/held.csv').write_text(label+' csv')
                (path/'held.png').write_bytes((label+' png').encode())
            receipt=dict(status='passed',n_formal_fits=147,downstream_exact_replay=True,csv_only_figures=True,candidate_manifest_sha256=module.digest(candidate/'manifest.json'))
            original_replace=module.os.replace
            def fail_figures(source,destination):
                if Path(source)==staged:raise OSError('injected figure replacement failure')
                return original_replace(source,destination)
            with patch.object(module.os,'replace',side_effect=fail_figures):
                with self.assertRaisesRegex(OSError,'injected'):
                    module.install(candidate,stable,receipt,figures_candidate=staged,figures_destination=figures)
            self.assertEqual((stable/'manifest.json').read_text(),'old model')
            self.assertEqual((figures/'data/held.csv').read_text(),'old csv')
            self.assertEqual((figures/'held.png').read_bytes(),b'old png')
            self.assertEqual((candidate/'manifest.json').read_text(),'new model')
            module.install(candidate,stable,receipt,figures_candidate=staged,figures_destination=figures)
            self.assertEqual((stable/'manifest.json').read_text(),'new model')
            self.assertEqual((figures/'data/held.csv').read_text(),'new csv')
            self.assertEqual((figures/'held.png').read_bytes(),b'new png')

    def test_preservation_snapshot_detects_protected_byte_changes(self):
        module=publication()
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);(root/'accepted.json').write_text('prior accepted bytes')
            snapshot=module.file_snapshot(root,['accepted.json'])
            module.verify_file_snapshot(root,snapshot)
            (root/'accepted.json').write_text('changed bytes')
            with self.assertRaisesRegex(ValueError,'preservation'):module.verify_file_snapshot(root,snapshot)


if __name__=='__main__':unittest.main()
