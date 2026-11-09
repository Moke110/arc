import tempfile
import unittest
from pathlib import Path

from artifacts import sha256, write_json
from paired_publication import build_figure, main, plot_points, resources


class PairedPublication(unittest.TestCase):
    def test_scheduler_rejects_completed_prediction_for_an_old_plan(self):
        with tempfile.TemporaryDirectory() as temp:
            runtime = Path(temp)
            step = runtime / 'stage2/06_pdiff_transfer_optimize'
            root = step / 'output/paired-transfer'
            cell = root / 'reference/seed-17'
            write_json(cell / 'plan.json', {'version': 'current'})
            prediction = cell / 'Transfer/prediction.json'
            write_json(prediction, {'version': 'current'})
            slot = root / 'schedule/reference-predict'
            scheduler = slot / 'scheduler'
            scheduler.mkdir(parents=True)
            job = 'reference-17-predict'
            queued = {
                'id': job, 'cwd': str(runtime), 'kind': 'gpu', 'memory_gib': 10,
                'disk_gib': 2, 'timeout_seconds': 1800,
                'command': [str(runtime / 'stage2/_envs/perturbdiff/bin/python'),
                            str(step / 'execution/script/paired_reference.py'), 'predict', str(cell),
                            '--runtime', str(runtime), '--seed', '17', '--plan-sha256', '0' * 64],
            }
            write_json(slot / 'queue.json', {'record_dir': 'scheduler', 'jobs': [queued]})
            log = scheduler / f'{job}.attempt1.log'
            log.write_text(f'{prediction} SHA-256 {sha256(prediction)}\n')
            write_json(scheduler / 'status.json', {'jobs': {job: {
                'identity': queued, 'status': 'completed',
                'attempts': [{'log': str(log), 'outcome': 'completed', 'returncode': 0,
                              'started_at': 1, 'finished_at': 2}],
            }}})
            (scheduler / 'telemetry.jsonl').touch()
            with self.assertRaisesRegex(ValueError, f'{job}.*sealed plan'):
                resources(root, runtime)

    def test_figure_replays_absolute_means_with_the_declared_panels_and_undefined_support(self):
        runs = {}
        for i, arm in enumerate(('Transfer', 'random', 'log', 'all-half', 'main-half', 'adaLN-only')):
            runs[arm] = [
                {'seed': seed, 'per_dose': [
                    {'dose_um': dose, 'mean_ap100': 0.1 * i + value,
                     'pcc': None if arm == 'log' and seed == 2026 and dose == '10.0' else value}
                    for dose in ('0.2', '1.0', '10.0')
                ]}
                for seed, value in [(17, 0.1), (42, 0.2), (2026, 0.6)]
            ]
        points = plot_points(runs)
        fig = build_figure(points)
        top = ['Original Transfer', 'random', 'log']
        bottom = ['Original Transfer', 'all-half', 'main-half', 'adaLN-only']
        for ax, names, ylabel in zip(fig.axes, (top, top, bottom, bottom),
                                     ('AUPRC (mean AP@100)', 'PCC') * 2, strict=True):
            self.assertEqual([line.get_label() for line in ax.lines], names)
            self.assertEqual(ax.get_ylabel(), ylabel)
            self.assertEqual([label.get_text() for label in ax.get_xticklabels()], ['0.2', '1', '10'])
            self.assertAlmostEqual(ax.lines[0].get_ydata()[0], 0.3)
        self.assertAlmostEqual(fig.axes[0].lines[1].get_ydata()[0], 0.4)
        self.assertIsNone(fig.axes[1].lines[2].get_ydata()[-1])
        self.assertEqual(points['log'][-1]['pcc_undefined_seeds'], [2026])
        self.assertEqual({ax.lines[0].get_color() for ax in fig.axes}, {fig.axes[0].lines[0].get_color()})

        with tempfile.TemporaryDirectory() as temp:
            temp = Path(temp)
            original, record, replay = (temp / name for name in ('original.png', 'record.json', 'replay.png'))
            fig.savefig(original)
            write_json(record, {'schema': 'step06-paired-publication-v1', 'runs': runs,
                                'plotted': points, 'figure_sha256': sha256(original)})
            main(['replay', str(record), '--figure', str(replay)])
            self.assertEqual(original.read_bytes(), replay.read_bytes())

        runs['random'].pop()
        with self.assertRaisesRegex(ValueError, 'three paired seeds'):
            plot_points(runs)
