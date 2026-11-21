"""Coordinates and CSV I/O shared by report preparation and rendering."""

import csv

METHODS = ('no_dose', 'latent_numeric', 'latent_embedding', 'FiLM')
METHOD_LABELS = ('No dose', 'Latent numeric', 'Latent embedding', 'FiLM')
METHOD_COLORS = ('#386095', '#ad7745', '#8970b4', '#a65170')
WIDTHS = ('original', '1.5x', '2x', '2.5x')
FOLDS = (0, 1, 2)
FOLD_COLORS = ('#23669d', '#c16041', '#43866a')
METRICS = {
    'pcc': 'Δ Pearson r',
    'spearman': 'Δ Spearman ρ',
    'rmse': 'Δ RMSE',
    'r2': 'Δ R²',
}
OVERALL = 'metrics_balloon'
BAR = 'metrics_bar'
PER_DOSE = 'metrics_per_dose_ballone'
CURVES = 'learning_curves'


def read_rows(directory, name):
    with (directory / f'{name}.csv').open(newline='', encoding='utf-8') as stream:
        return list(csv.DictReader(stream))


def write_rows(directory, name, rows):
    with (directory / f'{name}.csv').open('w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
