"""Read removed metadata through Git without restoring historical score files."""
import hashlib
import json
from pathlib import Path
import subprocess

REPO = next(parent for parent in Path(__file__).resolve().parents if (parent / '.git').exists())
COMMIT = '301f70d93ce17518fc842009bfbca727aed2ef0c'
# Only original admission/exposure metadata may use this compatibility path.
PATHS = {
    'stage2/03_modes_tenv/execution/final_record/evaluation.json',
    'stage2/04_pdiff_hyperparameter/execution/formal_record/tenv/evaluation.json',
    'stage2/04_pdiff_hyperparameter/execution/formal_record/tenv/per_run.csv',
    'stage2/05_pdiff_transfer/execution/formal_record/publication.json',
    'stage2/05_pdiff_transfer/execution/formal_record/paired_ap100.json',
    'stage2/05_pdiff_transfer/execution/tahoe_zero_shot_record/qualification.json',
    'stage2/05_pdiff_transfer/execution/transfer_record/qualification.json',
}


def original_bytes(path, expected=None):
    path = Path(path).absolute()
    relative = str(path.relative_to(REPO))
    if relative not in PATHS:
        raise ValueError(f'not an authorized historical metadata path: {relative}')
    subprocess.run(['git', '-C', str(REPO), 'merge-base', '--is-ancestor', COMMIT, 'dev'], check=True)
    value = subprocess.check_output(['git', '-C', str(REPO), 'show', f'{COMMIT}:{relative}'])
    if expected is not None and hashlib.sha256(value).hexdigest() != expected:
        raise ValueError(f'historical metadata digest differs: {relative}')
    return value


def original_json(path, expected=None):
    return json.loads(original_bytes(path, expected))


def metadata_digest(path):
    return hashlib.sha256(original_bytes(path)).hexdigest()
