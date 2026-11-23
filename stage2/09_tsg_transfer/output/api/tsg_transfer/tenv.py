"""Freeze the Step 03 tenv cold-drug-dose split and evaluate its test profiles."""
import json
from pathlib import Path
import secrets
import tempfile

import numpy as np

from tsg_round1.inputs import sha256
from tsg_round1.protocol import make_splits, score_prediction, summary

DOSES = (0.2, 1.0, 10.0)
FOLDS = range(len(DOSES))
PARTITIONS = ('train', 'validation', 'test', 'excluded-held-parent', 'excluded-held-dose')


def _split(data):
    if data.source != 'tenv' or tuple(sorted(data.obs.dose_um.unique())) != DOSES:
        raise ValueError('expected tenv with the three registered physical doses')
    splits = make_splits(data.obs)
    return splits.loc[splits.scenario == 'cold-drug-dose']


def _write_json(path, value):
    path.write_text(json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + '\n')


def _split_receipt(data, directory, splits):
    return dict(
        schema='tsg-transfer-tenv-split-v1',
        input_manifest_sha256=sha256(data.directory / 'manifest.json'),
        splits_sha256=sha256(directory / 'splits.csv'),
        folds=[dict(
            fold=fold,
            held_dose_um=dose,
            counts={role: int(((splits.fold == fold) & (splits.partition == role)).sum())
                    for role in PARTITIONS},
        ) for fold, dose in enumerate(DOSES)],
    )


def freeze_split(data, directory):
    """Seal a complete observation-by-fold role table under the input identity."""
    directory = Path(directory)
    splits = _split(data)
    if directory.exists():
        _verify_split(data, directory, splits)
        return directory / 'manifest.json'
    directory.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=directory.parent) as temporary:
        staged = Path(temporary)
        (staged / 'splits.csv').write_text(splits.to_csv(index=False))
        _write_json(staged / 'manifest.json', _split_receipt(data, staged, splits))
        staged.rename(directory)
    return directory / 'manifest.json'


def _verify_split(data, directory, splits):
    receipt = json.loads((directory / 'manifest.json').read_text())
    if (receipt != _split_receipt(data, directory, splits)
            or (directory / 'splits.csv').read_text() != splits.to_csv(index=False)):
        raise ValueError('tenv split differs from sealed Step 03 identities')


def load_split(data, directory, fold):
    """Verify the full frozen role table, then return ordered row IDs by role."""
    if fold not in FOLDS:
        raise ValueError('unknown tenv fold')
    splits = _split(data)
    _verify_split(data, Path(directory), splits)
    selected = splits.loc[splits.fold == fold]
    return {role: selected.loc[selected.partition == role, 'row'].to_numpy(dtype=np.int64)
            for role in PARTITIONS}


EVALUATION_PATHS = ('tenv-only', 'tahoe-only', 'transfer')
TEST_FILES = ('identity.csv', 'prediction.npy', 'scores.csv')
PROCESS_TOKEN = secrets.token_hex(32)


def _test_identity(data, rows):
    return data.obs.iloc[rows][['profile_id', 'drug_parent_inchikey', 'dose_um', 'cell_line']]


def _test_values(data, rows, prediction, prediction_ids):
    if prediction.dtype != np.float32 or not np.isfinite(prediction).all():
        raise ValueError('expected finite float32 test predictions')
    return score_prediction(data.obs.iloc[rows], data.expression[rows], prediction,
                            data.control[rows], prediction_ids=prediction_ids)


def _test_binding(data, split, fold, path, lineage):
    if path not in EVALUATION_PATHS or not lineage:
        raise ValueError('unregistered test path or missing checkpoint lineage')
    return dict(
        schema='tsg-transfer-tenv-test-v1', fold=fold, path=path, lineage=lineage,
        input_manifest_sha256=sha256(data.directory / 'manifest.json'),
        split_manifest_sha256=sha256(Path(split) / 'manifest.json'),
    )


def export_test(data, split, fold, output, prediction, *, path, lineage, prediction_ids):
    """Write the ordered model output, sample scores and content-bound receipt."""
    output = Path(output)
    rows = load_split(data, split, fold)['test']
    binding = _test_binding(data, split, fold, path, lineage)
    if output.exists():
        raise FileExistsError(output)
    prediction = np.asarray(prediction)
    scores = _test_values(data, rows, prediction, prediction_ids)
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=output.parent) as temporary:
        staged = Path(temporary)
        _test_identity(data, rows).to_csv(staged / 'identity.csv', index=False)
        scores.to_csv(staged / 'scores.csv', index=False)
        np.save(staged / 'prediction.npy', prediction, allow_pickle=False)
        _write_json(staged / 'manifest.json', dict(
            binding, metrics=summary(scores), evaluation_process=PROCESS_TOKEN,
            files={name: sha256(staged / name) for name in TEST_FILES},
        ))
        staged.rename(output)
    return output / 'manifest.json'


def verify_test(data, split, fold, output, *, path, lineage, require_replay=True):
    """Validate lineage, exact ordered test identities and sample-equal scoring."""
    output = Path(output)
    receipt = json.loads((output / 'manifest.json').read_text())
    binding = _test_binding(data, split, fold, path, lineage)
    if (any(receipt.get(key) != value for key, value in binding.items())
            or not receipt.get('evaluation_process')
            or receipt.get('files') != {name: sha256(output / name) for name in TEST_FILES}):
        raise ValueError('test receipt lineage or artifact digests differ')

    rows = load_split(data, split, fold)['test']
    identity = _test_identity(data, rows)
    scores = _test_values(data, rows, np.load(output / 'prediction.npy', mmap_mode='r'),
                          identity.profile_id)
    if (identity.to_csv(index=False) != (output / 'identity.csv').read_text()
            or scores.to_csv(index=False) != (output / 'scores.csv').read_text()
            or receipt.get('metrics') != summary(scores)):
        raise ValueError('test identity or Delta score evidence differs')
    if require_replay:
        proof = output / 'replay.json'
        if not proof.is_file():
            raise ValueError('test requires fresh-process prediction replay')
        replay = json.loads(proof.read_text())
        if (replay.get('test_manifest_sha256') != sha256(output / 'manifest.json')
                or not replay.get('process') or replay['process'] == receipt['evaluation_process']):
            raise ValueError('test replay does not bind a separate process and the test receipt')
    return receipt


def replay_test(data, split, fold, output, prediction, *, path, lineage, prediction_ids):
    """Compare predictions generated in a separate process against the sealed array."""
    output = Path(output)
    receipt = verify_test(data, split, fold, output, path=path, lineage=lineage,
                          require_replay=False)
    if receipt['evaluation_process'] == PROCESS_TOKEN:
        raise ValueError('test replay must run in a fresh process')
    rows = load_split(data, split, fold)['test']
    if list(prediction_ids) != list(data.obs.index[rows]):
        raise ValueError('replayed test prediction identity order differs')
    stored = np.load(output / 'prediction.npy', mmap_mode='r')
    if prediction.dtype != stored.dtype or not np.array_equal(stored, prediction):
        raise ValueError('exact fresh-process test prediction replay differs')
    with (output / 'replay.json').open('x') as proof:
        json.dump(dict(test_manifest_sha256=sha256(output / 'manifest.json'),
                       process=PROCESS_TOKEN), proof, sort_keys=True, indent=2)
        proof.write('\n')
    return verify_test(data, split, fold, output, path=path, lineage=lineage)


def control_lineage(fit_dir):
    """Bind an evaluation to its verified tenv-only fit and replay proof."""
    fit_dir = Path(fit_dir)
    return dict(fit_receipt_sha256=sha256(fit_dir / 'fit.json'),
                fit_replay_sha256=sha256(fit_dir / 'validation-replay.json'))
