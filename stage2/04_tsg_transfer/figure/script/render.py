"""Render Stage 2 Step 04's paired Delta PCC comparison from plotted CSVs."""
import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def render(data, output):
    with (data / 'transfer_scores.csv').open(newline='') as stream:
        scores = list(csv.DictReader(stream))
    with (data / 'transfer_differences.csv').open(newline='') as stream:
        differences = list(csv.DictReader(stream))
    folds = [row['fold'] for row in differences]
    labels = [f"Fold {row['fold']}\n{float(row['held_dose_um']):g} µM"
              if row['fold'] != 'mean' else 'Equal-fold\nmean' for row in differences]
    x = list(range(len(folds)))
    fig, (absolute, contrasted) = plt.subplots(ncols=2, figsize=(10, 4.3), layout='constrained')
    colors = ('#566575', '#B76942', '#206D73')
    for offset, model, label, color in zip(
        (-0.23, 0, 0.23), ('tenv-only', 'tahoe-only', 'transfer'),
        ('tenv-only', 'Tahoe-only', 'Transfer'), colors, strict=True,
    ):
        values = {row['fold']: float(row['pcc']) for row in scores if row['path'] == model}
        absolute.bar([index + offset for index in x], [values[fold] for fold in folds],
                     width=0.20, color=color, label=label)
    for offset, key, label, color in (
        (-0.12, 'transfer_minus_tenv_only', 'Transfer − tenv-only', colors[0]),
        (0.12, 'transfer_minus_tahoe_only', 'Transfer − Tahoe-only', colors[1]),
    ):
        contrasted.bar([index + offset for index in x], [float(row[key]) for row in differences],
                       width=0.20, color=color, label=label)
    for axis in (absolute, contrasted):
        axis.set_xticks(x, labels)
        axis.axhline(0, linewidth=0.8, color='#1F2932')
        axis.axvline(2.5, linewidth=0.7, color='#B5BDC6')
        axis.spines[['top', 'right']].set_visible(False)
    absolute.legend(loc='upper center', bbox_to_anchor=(0.5, -0.22), ncol=3, frameon=False, fontsize=8)
    contrasted.legend(loc='upper center', bbox_to_anchor=(0.5, -0.22), ncol=2, frameon=False, fontsize=8)
    absolute.set_title('Held-out scores', loc='left')
    absolute.set_ylabel('Sample-equal Delta PCC')
    contrasted.set_title('Paired differences', loc='left')
    contrasted.set_ylabel('Difference in Delta PCC')
    fig.savefig(output, dpi=180, facecolor='white')
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    render(args.data, args.output)


if __name__ == '__main__':
    main()
