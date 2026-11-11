"""Preserve the non-retrieval Tahoe source-exposure figure."""
import argparse
import csv
from pathlib import Path
REPO = next(parent for parent in Path(__file__).resolve().parents if (parent / '.git').exists())
STEP = REPO / 'stage2/06_pdiff_transfer_optimize'


def render_exposure(exposure, destination):
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure
    destination.mkdir(parents=True, exist_ok=True)
    figure = Figure(figsize=(8.5, 4.1), layout='constrained')
    FigureCanvasAgg(figure)
    ax = figure.subplots()
    roles = ('train', 'validation', 'test')
    cohorts = {(row['source_release'], row['tahoe_role']) for row in exposure}
    counts = {(row['source_release'], row['tahoe_role'], row['tenv_role']): int(row['shared_parent_count'])
              for row in exposure}
    for index, (release, role) in enumerate(sorted(cohorts, reverse=True)):
        offset = (index - (len(cohorts) - 1) / 2) * .19
        ax.bar([i + offset for i in range(3)], [counts[release, role, key] for key in roles], width=.18,
               label=f'{release.replace("_", " ").capitalize()} Tahoe {role}')
    ax.set_xticks(range(3), [key.capitalize() for key in roles])
    ax.set(xlabel='Tenv fold 0 role', ylabel='Shared parent-connectivity blocks',
           title='Tahoe fitting exposure before and after global exclusion')
    ax.legend(frameon=False, fontsize=8)
    figure.savefig(destination / "02-tahoe-exposure.png", dpi=144)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=STEP / 'source_exposure.csv')
    parser.add_argument('--output', type=Path, default=STEP / 'output')
    args = parser.parse_args()
    with args.source.open() as stream:
        render_exposure(list(csv.DictReader(stream)), args.output)


if __name__ == '__main__':
    main()
