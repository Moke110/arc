import tempfile
import unittest
from pathlib import Path

from artifacts import sha256, write_json
from paired_publication import build_figure, main, plot_points


class PairedPublication(unittest.TestCase):
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
