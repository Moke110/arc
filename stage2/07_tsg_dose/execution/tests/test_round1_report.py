"""Approved saved-report export, CSV renderer and publication boundaries."""
import csv
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

STEP=Path(__file__).resolve().parents[2]

def report():
    spec=importlib.util.spec_from_file_location('round1_report',STEP/'execution/publish_report.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module

def write_csv(path,rows):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)

def fixture(matrix):
    methods=['no-dose','scaling','input-numeric','latent-numeric','input-embedding','latent-embedding','film']
    for scenario in ['cold-drug','held-dose']:
        for fold in range(3):
            selected=matrix/f'selections/{scenario}/{fold}/receipt.json';selected.parent.mkdir(parents=True)
            selected.write_text(json.dumps(dict(method=['scaling','film','input-numeric'][fold])))
            for method in methods:
                for index,seed in enumerate([17,42,2026]):
                    value=(fold*3+index+1)/10
                    rows=[dict(drug_parent=parent,dose_um=dose,n_wells=2,pcc=score,spearman=score,rmse=score,r2=score,constant_truth=0,constant_prediction=0) for parent,dose,score in [('A',.2,value-.1),('A',1,value+.1),('B',10,value)]]
                    write_csv(matrix/f'runs/{scenario}/{method}/{fold}/{seed}/test_per_group.csv',rows)

def scientific_fixture(module,matrix,scientific):
    fixture(matrix)
    (matrix/'plan.json').write_text('{}')
    runs=[dict(scenario=scenario,fold=fold,method=method,seed=seed)
          for scenario in module.SCENARIOS for fold in range(3)
          for method in module.METHODS for seed in module.SEEDS]
    manifests={}
    for run in runs:
        origin=matrix/f"runs/{run['scenario']}/{run['method']}/{run['fold']}/{run['seed']}"
        manifest=origin/'manifest.json'
        manifest.write_text(json.dumps(dict(outputs={'test_per_group.csv':module.digest(origin/'test_per_group.csv')})))
        manifests[str(manifest.relative_to(matrix))]=module.digest(manifest)
    sources={}
    for name in ['reports/cold-drug/per_dose.csv','reports/cold-drug/learning_curves.csv','reports/held-dose/learning_curves.csv']:
        source=matrix/name;source.parent.mkdir(parents=True,exist_ok=True)
        source.write_text('saved scientific report\n')
        sources[name]=module.digest(source)
    scientific.write_text(json.dumps(dict(schema='tsg-dose-round1-acceptance-v1',
        scenarios=module.SCENARIOS,n_formal_fits=126,runs=runs,
        plan_sha256=module.digest(matrix/'plan.json'),origin_manifest_sha256=manifests,
        report_source_sha256=sources,selection_sha256={str(p.relative_to(matrix)):module.digest(p)
                     for p in matrix.glob('selections/*/*/receipt.json')})))

class Round1ReportTests(unittest.TestCase):
    def test_export_uses_nine_parent_equal_run_metrics_and_sample_sd(self):
        module=report()
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);fixture(root/'matrix')
            module.export_run_metrics(root/'matrix',root/'data')
            runs=module.rows(root/'data/round1_runs.csv');summaries=module.rows(root/'data/round1_summary.csv')
            part=[r for r in runs if r['scenario']=='cold-drug' and r['method']=='selected' and r['metric']=='pcc']
            self.assertEqual(len(part),9)
            self.assertEqual([r['actual_method'] for r in part],['scaling']*3+['film']*3+['input-numeric']*3)
            summary=next(r for r in summaries if r['scenario']=='cold-drug' and r['method']=='selected' and r['metric']=='pcc')
            self.assertAlmostEqual(float(summary['mean']),.5)
            self.assertAlmostEqual(float(summary['sd']),.27386127875258304)
            self.assertEqual((summary['n_total'],summary['n_defined'],summary['n_undefined']),('9','9','0'))
            self.assertFalse(any(r['scenario']=='held-dose' and r['method']=='selected' for r in runs))

    def test_export_preserves_four_dose_metrics_and_every_curve_coordinate(self):
        module=report()
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);matrix=root/'matrix'
            doses=[dict(method=m,dose_um=d,pcc=.1,spearman=.2,rmse=.3,r2=-.4,n_parents=2,n_seeds=3) for m in module.METHODS+['selected'] for d in [.2,1,10]]
            write_csv(matrix/'reports/cold-drug/per_dose.csv',doses)
            for scenario in module.SCENARIOS:
                curves=[dict(scenario=scenario,method=m,fold=f,seed=s,epoch=e,train_loss=12.5+e,validation_mse=1.5/e) for m in module.METHODS for f in range(3) for s in module.SEEDS for e in [1,2]]
                write_csv(matrix/f'reports/{scenario}/learning_curves.csv',curves)
            module.export_coordinates(matrix,root/'data')
            actual=module.rows(root/'data/round1_dose.csv')
            self.assertEqual(len(actual),96)
            self.assertEqual({r['metric'] for r in actual},{'pcc','spearman','rmse','r2'})
            self.assertEqual(next(r['value'] for r in actual if r['metric']=='r2'),'-0.4')
            curves=module.rows(root/'data/round1_curves.csv')
            self.assertEqual(len(curves),252)
            self.assertEqual({r['scenario'] for r in curves},{'cold-drug','held-dose'})
            self.assertEqual(curves[-1]['train_loss'],'14.5')
            self.assertEqual(len(curves[-1]['source_sha256']),64)

    def test_csv_only_renderer_emits_exactly_three_round1_figures(self):
        import subprocess,sys
        module=report()
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);matrix=root/'matrix';fixture(matrix);data=root/'data'
            module.export_run_metrics(matrix,data)
            write_csv(matrix/'reports/cold-drug/per_dose.csv',[dict(method=m,dose_um=d,pcc=.1,spearman=.2,rmse=.3,r2=-.4,n_parents=2,n_seeds=3) for m in module.METHODS+['selected'] for d in [.2,1,10]])
            for scenario in module.SCENARIOS:
                write_csv(matrix/f'reports/{scenario}/learning_curves.csv',[dict(scenario=scenario,method=m,fold=f,seed=s,epoch=e,train_loss=12.5/e,validation_mse=1.5/e) for m in module.METHODS for f in range(3) for s in module.SEEDS for e in [1,2]])
            module.export_coordinates(matrix,data)
            import shutil
            shutil.copy2(module.historical_report_root()/'data/plot_style.csv',data/'plot_style.csv')
            rendered=module.render_isolated(data,root/'figures',sys.executable)
            self.assertEqual(set(rendered),module.EXPECTED)
            self.assertTrue(all((root/'figures'/n).read_bytes().startswith(b'\x89PNG') for n in rendered))

    def test_nonfinite_parent_metrics_are_counted_without_zero_substitution(self):
        module=report()
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);fixture(root/'matrix')
            source=root/'matrix/runs/cold-drug/no-dose/0/17/test_per_group.csv'
            part=module.rows(source)
            for r in part:r.update(pcc='nan',constant_truth='False',constant_prediction='True')
            write_csv(source,part)
            module.export_run_metrics(root/'matrix',root/'data')
            summary=next(r for r in module.rows(root/'data/round1_summary.csv') if (r['scenario'],r['method'],r['metric'])==('cold-drug','no-dose','pcc'))
            self.assertAlmostEqual(float(summary['mean']),.55)
            self.assertEqual((summary['n_defined'],summary['n_undefined'],summary['n_nonfinite_parent_dose'],summary['n_undefined_parents']),('8','1','3','2'))

    def test_failed_report_install_restores_prior_report_and_rejects_changed_candidate(self):
        from unittest.mock import patch
        module=report()
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);candidate=root/'candidate';stable=root/'figure'
            candidate.mkdir();stable.mkdir();(stable/'old.png').write_bytes(b'old')
            for name in module.EXPECTED:(candidate/name).write_bytes(b'new')
            receipt=dict(status='passed',n_formal_fits=126,csv_only_figures=True,visual_inspection=True,files=module.snapshot(candidate))
            original=module.os.replace
            def fail_new(source,destination):
                if Path(source)==candidate:raise OSError('injected replacement failure')
                return original(source,destination)
            with patch.object(module.os,'replace',side_effect=fail_new):
                with self.assertRaisesRegex(OSError,'injected'):module.install(candidate,stable,receipt)
            self.assertEqual((stable/'old.png').read_bytes(),b'old');self.assertTrue(candidate.exists())
            (candidate/'round1_metrics.png').write_bytes(b'corrupted')
            with self.assertRaisesRegex(ValueError,'candidate'):module.install(candidate,stable,receipt)
            self.assertEqual((stable/'old.png').read_bytes(),b'old')

    def test_changed_sealed_choice_rejects_preparation_and_preserves_prior_report(self):
        module=report()
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);matrix=root/'matrix';scientific=root/'scientific.json'
            scientific_fixture(module,matrix,scientific)
            module.verify_scientific_receipt(matrix,scientific)
            stable=root/'figure';stable.mkdir();(stable/'old.png').write_bytes(b'previous accepted report')
            choice=matrix/'selections/cold-drug/0/receipt.json'
            choice.write_text(json.dumps(dict(method='film')))
            candidate=root/'candidate'
            with self.assertRaisesRegex(ValueError,'selection|choice'):
                module.prepare(matrix,candidate,stable,scientific,'unused-python')
            self.assertFalse(candidate.exists())
            self.assertEqual((stable/'old.png').read_bytes(),b'previous accepted report')

    def test_changed_origin_payload_and_manifest_reject_scientific_evidence(self):
        module=report()
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);matrix=root/'matrix';scientific=root/'scientific.json'
            scientific_fixture(module,matrix,scientific)
            module.verify_scientific_receipt(matrix,scientific)
            origin=matrix/'runs/cold-drug/no-dose/0/17'
            changed=module.rows(origin/'test_per_group.csv');changed[0]['pcc']='.99'
            write_csv(origin/'test_per_group.csv',changed)
            (origin/'manifest.json').write_text(json.dumps(dict(outputs={'test_per_group.csv':module.digest(origin/'test_per_group.csv')})))
            with self.assertRaisesRegex(ValueError,'manifest'):
                module.verify_scientific_receipt(matrix,scientific)

    def test_changed_saved_curve_or_dose_report_rejects_scientific_evidence(self):
        module=report()
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);matrix=root/'matrix';scientific=root/'scientific.json'
            scientific_fixture(module,matrix,scientific)
            module.verify_scientific_receipt(matrix,scientific)
            for name in ['reports/cold-drug/per_dose.csv','reports/cold-drug/learning_curves.csv','reports/held-dose/learning_curves.csv']:
                with self.subTest(source=name):
                    source=matrix/name;original=source.read_bytes()
                    source.write_bytes(b'changed scientific values\n')
                    with self.assertRaisesRegex(ValueError,'scientific report'):
                        module.verify_scientific_receipt(matrix,scientific)
                    source.write_bytes(original)
