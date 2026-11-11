"""Remove superseded retrieval statistics and figures; Git retains history."""
import argparse
import json
from pathlib import Path
import re

REPO = next(parent for parent in Path(__file__).resolve().parents if (parent / '.git').exists())
ROOT = REPO / 'stage2'
LEGACY = re.compile(r'(?:AUPRC_(?:positive|negative)@\d+|auprc_(?:full|g)|(?:positive|negative|mean)_ap100|AP@(?:100|200)|auprc@(?:full|G))', re.I)


def old_values(value):
    if isinstance(value, dict):
        if LEGACY.search(str(value.get('metric', ''))):
            return True
        for key, item in value.items():
            if (LEGACY.search(key) or 'auprc' in key.lower()) and isinstance(item, (int, float, list, dict)):
                return True
            if old_values(item):
                return True
    elif isinstance(value, list):
        return any(old_values(x) for x in value)
    return False


def candidates():
    paths = set()
    # Fixed workflow areas include ignored old scoring releases as well as tracked ones.
    for step in ROOT.iterdir():
        if not step.is_dir() or step.name.startswith(('_envs', 'evaluation', 'presentation')):
            continue
        if not (step.name[:2].isdigit() or step.name in ('_archive', 'audit')):
            continue
        for path in step.rglob('*'):
            if not path.is_file() or any(p in ('source', 'local', '__pycache__', '.git') for p in path.parts):
                continue
            if 'retrieval' in path.parts or 'labels' in path.parts or 'input' in path.parts or 'corpus' in path.parts:
                continue
            if path.suffix not in ('.json', '.csv', '.tsv', '.npz', '.png', '.svg'):
                continue
            name = path.name.lower()
            if path.suffix in ('.png', '.svg'):
                old = (re.search(r'(^|[_-])(pr|prc|auprc)([_-]|\.)', name) or name.startswith('best_run')
                       or (step.name.startswith(('03_', '04_', '05_', '09_')) and 'metrics' in name)
                       or (step.name.startswith('06_') and name.startswith(('01-paired', '03-paired'))))
                if old:
                    paths.add(path)
                continue
            if path.suffix == '.npz':
                if name.startswith('pr_'):
                    paths.add(path)
                continue
            if name.startswith(('pr_curves', 'pr_points', 'pr_bands', 'paired_ap100')):
                paths.add(path)
                continue
            if step.name.startswith('09_') and any(p in ('test_scores', 'validation_scores', 'validation_scores_directional', 'legacy_unsigned') for p in path.parts):
                paths.add(path)
                continue
            if path.suffix in ('.csv', '.tsv'):
                with path.open(errors='replace') as handle:
                    header = handle.readline()
                if LEGACY.search(header) or re.search(r'auprc', header, re.I):
                    paths.add(path)
                continue
            try:
                data = json.loads(path.read_text())
            except (ValueError, UnicodeDecodeError):
                continue
            if old_values(data):
                paths.add(path)
            elif name in ('evaluation.json', 'publication.json', 'directional-pr.json') and step.name.startswith(('03_', '04_', '05_', '06_', '09_')):
                paths.add(path)
    return sorted(paths)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--delete', action='store_true')
    args = parser.parse_args()
    paths = candidates()
    for path in paths:
        print(path.relative_to(REPO))
        if args.delete:
            path.unlink()
    print(f"{'Deleted' if args.delete else 'Found'} {len(paths)} obsolete retrieval artifacts")


if __name__ == '__main__':
    main()
