"""Bounded CPU training and independently verifiable saved Morgan runs."""
import hashlib
import json
from pathlib import Path
import platform
import random
import numpy as np
from metrics import METRIC_CONVENTION, evaluate_predictions, hard_predictions
from workflow import digest, encoded

CONFIG = dict(input_features=2049, hidden_units=128, classes=[-1, 0, 1],
              learning_rate=0.001, weight_decay=0.0001, batch_size=64,
              max_epochs=300, patience=30)


def _torch():
    import torch
    torch.set_num_threads(2)
    torch.use_deterministic_algorithms(True)
    return torch


def _model(torch):
    return torch.nn.Sequential(torch.nn.Linear(2049, 128), torch.nn.ReLU(), torch.nn.Linear(128, 3))


def _features(rows, statistics):
    from rdkit import Chem
    from rdkit.Chem import rdFingerprintGenerator
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048, includeChirality=True,
        countSimulation=False, useBondTypes=True, onlyNonzeroInvariants=False,
        includeRingMembership=True, includeRedundantEnvironments=False)
    values = np.empty((len(rows), 2049), dtype=np.float32)
    for i, row in enumerate(rows):
        molecule = Chem.MolFromSmiles(row['canonical_smiles'])
        dose = float(row['dose_um'])
        if molecule is None or not np.isfinite(dose) or dose <= 0:
            raise ValueError('Invalid structure or nonpositive/nonfinite dose')
        values[i, :2048] = generator.GetFingerprintAsNumPy(molecule)
        values[i, 2048] = (np.log10(dose)-statistics['mean'])/statistics['std']
    return values


def _sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _feature_configuration():
    import rdkit
    return dict(radius=2, bits=2048, chirality=True, rdkit_version=rdkit.__version__,
                effective_configuration=dict(countSimulation=False, useBondTypes=True,
                    onlyNonzeroInvariants=False, includeRingMembership=True,
                    includeRedundantEnvironments=False))


def _prepared(path):
    from workflow import SPLIT_SEED, TRAINING_SEEDS, CLASSES, reconstruct_splits
    value = json.loads(Path(path).read_text())
    if digest(value['rows']) != value['input_fingerprint'] or digest(value['splits']) != value['split_fingerprint']:
        raise ValueError('Prepared input or split fingerprint mismatch')
    if value['features'] != _feature_configuration():
        raise ValueError('Incompatible fingerprint configuration or RDKit version')
    if (value.get('schema_version') != 1 or value.get('scope') != '+'.join(value.get('settings', ['cold_drug']))
            or value.get('split_seed') != SPLIT_SEED or value.get('training_seeds') != TRAINING_SEEDS
            or value.get('class_order') != CLASSES):
        raise ValueError('Incompatible prepared configuration')
    rows = value['rows']
    if not rows or len({r['measurement_id'] for r in rows}) != len(rows):
        raise ValueError('Missing or duplicate measurement identity')
    from rdkit import Chem
    for row in rows:
        molecule = Chem.MolFromSmiles(row['canonical_smiles'])
        if (molecule is None or Chem.MolToSmiles(molecule, isomericSmiles=True) != row['canonical_smiles']
                or Chem.MolToInchiKey(molecule) != row['parent_id']
                or row['direction'] not in CLASSES or not np.isfinite(float(row['dose_um']))
                or float(row['dose_um']) <= 0):
            raise ValueError('Invalid prepared measurement')
    # Reconstruct the entire canonical partition, including counts and mappings.
    if value['splits'] != reconstruct_splits(rows, value.get('settings', ['cold_drug'])):
        raise ValueError('Prepared split membership/counts differ from canonical partitions')
    expected = [dict(condition=c, phenotype=p, setting=s['setting'], fold=s['fold'], seed=seed)
                for c, p in sorted({(r['condition'], r['phenotype']) for r in rows})
                for s in value['splits'] for seed in TRAINING_SEEDS]
    if value['expected_runs'] != expected:
        raise ValueError('Prepared expected run identities differ from canonical configuration')
    return value


def _statistics(rows):
    doses = np.log10([float(r['dose_um']) for r in rows])
    return dict(mean=float(doses.mean()), std=float(doses.std()) or 1.0, n=len(doses), ddof=0)


def predict(checkpoint_path, rows):
    """Infer with saved best weights and training-only dose statistics."""
    torch = _torch()
    checkpoint = torch.load(Path(checkpoint_path), map_location='cpu', weights_only=True)
    if checkpoint['config'] != CONFIG or checkpoint['features'] != _feature_configuration():
        raise ValueError('Incompatible checkpoint configuration/software')
    model = _model(torch)
    model.load_state_dict(checkpoint['state_dict'])
    model.eval()
    with torch.no_grad():
        return torch.softmax(model(torch.from_numpy(_features(rows, checkpoint['dose_statistics']))), dim=1).tolist()


def train(prepared_path, output_path, condition, phenotype, setting, fold, seed, smoke_epochs=None):
    """Write one immutable run directory; smoke runs cannot count as formal fits."""
    available = int(next(line.split()[1] for line in Path('/proc/meminfo').read_text().splitlines()
                         if line.startswith('MemAvailable:'))) * 1024
    if available < 12 * 1024**3:
        raise ValueError('Training requires at least 12 GiB available memory')
    prepared = _prepared(prepared_path)
    identity = dict(condition=condition, phenotype=phenotype, setting=setting, fold=fold, seed=seed)
    if seed not in [17, 42, 2026] or identity not in prepared['expected_runs']:
        raise ValueError('Unexpected run identity')
    if smoke_epochs is not None and (not isinstance(smoke_epochs, int) or not 1 <= smoke_epochs <= 300):
        raise ValueError('smoke_epochs must be between 1 and 300')
    output = Path(output_path)
    if output.exists():
        raise ValueError('Run output already exists; choose a fresh directory')
    split = next(s for s in prepared['splits'] if s['setting'] == setting and s['fold'] == fold)
    if len(split['roles']) != len(prepared['rows']):
        raise ValueError('Split roles do not align with rows')
    groups = {role: [r for r, assigned in zip(prepared['rows'], split['roles'])
                     if assigned == role and (r['condition'], r['phenotype']) == (condition, phenotype)]
              for role in ['train', 'validation', 'test']}
    if any(not rows for rows in groups.values()):
        raise ValueError('Every target needs nonempty train/validation/test partitions')
    doses = np.log10([float(r['dose_um']) for r in groups['train']])
    if not np.isfinite(doses).all():
        raise ValueError('Training dose must be positive and finite')
    statistics = _statistics(groups['train'])
    torch = _torch()
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    model = _model(torch)
    optimizer = torch.optim.AdamW(model.parameters(), lr=0.001, weight_decay=0.0001)
    x = torch.from_numpy(_features(groups['train'], statistics))
    y = torch.tensor([r['direction'] + 1 for r in groups['train']], dtype=torch.long)
    validation = torch.from_numpy(_features(groups['validation'], statistics))
    truth = [r['direction'] for r in groups['validation']]
    history, best_score, best_epoch, best_state, stale = [], -1.0, 0, None, 0
    generator = torch.Generator().manual_seed(seed)
    for epoch in range(1, (smoke_epochs or 300)+1):
        model.train()
        total = 0.0
        order = torch.randperm(len(x), generator=generator)
        for indices in order.split(64):
            optimizer.zero_grad()
            loss = torch.nn.functional.cross_entropy(model(x[indices]), y[indices])
            loss.backward()
            optimizer.step()
            total += float(loss.detach()) * len(indices)
        model.eval()
        with torch.no_grad():
            probabilities = torch.softmax(model(validation), dim=1).tolist()
        score = evaluate_predictions(truth, probabilities)['macro_f1']
        history.append(dict(epoch=epoch, train_loss=total/len(x), validation_macro_f1=score))
        if score > best_score:
            best_score, best_epoch, stale = score, epoch, 0
            best_state = {k: v.detach().clone() for k, v in model.state_dict().items()}
        else:
            stale += 1
        if stale >= 30:
            break
    output.mkdir(parents=True)
    checkpoint = dict(state_dict=best_state, config=CONFIG, features=prepared['features'],
                      dose_statistics=statistics, identity=identity, selected_epoch=best_epoch,
                      input_fingerprint=prepared['input_fingerprint'], split_fingerprint=prepared['split_fingerprint'])
    torch.save(checkpoint, output/'checkpoint.pt')
    probabilities = predict(output/'checkpoint.pt', groups['test'])
    labels = hard_predictions(probabilities).tolist()
    records = [dict(row, run_identity=identity, probabilities=p, prediction=label)
               for row, p, label in zip(groups['test'], probabilities, labels)]
    (output/'predictions.json').write_bytes(encoded(records)+b'\n')
    (output/'history.json').write_bytes(encoded(history)+b'\n')
    metrics = evaluate_predictions([r['direction'] for r in records], probabilities, labels)
    (output/'metrics.json').write_bytes(encoded(metrics)+b'\n')
    manifest = dict(identity=identity, status='completed', mode='smoke' if smoke_epochs else 'formal',
                    smoke_epochs=smoke_epochs, config=CONFIG, features=prepared['features'],
                    dose_statistics=statistics, selected_epoch=best_epoch, best_validation_macro_f1=best_score,
                    input_fingerprint=prepared['input_fingerprint'], split_fingerprint=prepared['split_fingerprint'],
                    metric_convention=METRIC_CONVENTION,
                    execution=dict(device='cpu', torch_version=str(torch.__version__), numpy_version=np.__version__,
                                   python_version=platform.python_version(), platform=platform.platform(),
                                   threads=2, deterministic_algorithms=True, available_memory_bytes=available),
                    artifact_hashes={name: _sha(output/name) for name in ['checkpoint.pt', 'predictions.json', 'history.json', 'metrics.json']})
    (output/'run.json').write_bytes(encoded(manifest)+b'\n')
    return manifest


def evaluate(run_path):
    """Regenerate metrics solely from saved prediction rows, without fitting."""
    records = json.loads((Path(run_path)/'predictions.json').read_text())
    return evaluate_predictions([r['direction'] for r in records], [r['probabilities'] for r in records],
                                [r['prediction'] for r in records])


def verify(run_path, prepared_path=None):
    """Reject changed bytes, mismatched test identities or inconsistent reload inference."""
    run = Path(run_path)
    manifest = json.loads((run/'run.json').read_text())
    for name in ['checkpoint.pt', 'predictions.json', 'history.json', 'metrics.json']:
        if _sha(run/name) != manifest['artifact_hashes'][name]:
            raise ValueError(f'Artifact hash mismatch: {name}')
    torch = _torch()
    checkpoint = torch.load(run/'checkpoint.pt', map_location='cpu', weights_only=True)
    for key in ['identity', 'features', 'config', 'dose_statistics', 'selected_epoch', 'input_fingerprint', 'split_fingerprint']:
        if checkpoint[key] != manifest[key]:
            raise ValueError(f'Checkpoint metadata mismatch: {key}')
    if manifest['config'] != CONFIG or manifest['features'] != _feature_configuration():
        raise ValueError('Incompatible saved configuration')
    smoke = manifest['smoke_epochs']
    if (manifest['status'] != 'completed' or manifest['mode'] != ('formal' if smoke is None else 'smoke')
            or (smoke is not None and (type(smoke) is not int or not 1 <= smoke <= CONFIG['max_epochs']))):
        raise ValueError('Inconsistent formal/smoke configuration')
    history = json.loads((run/'history.json').read_text())
    limit = smoke or CONFIG['max_epochs']
    if (not history or len(history) > limit
            or [h['epoch'] for h in history] != list(range(1, len(history)+1))):
        raise ValueError('Invalid training history schedule')
    best = max(history, key=lambda item: item['validation_macro_f1'])
    if best['epoch'] != manifest['selected_epoch'] or best['validation_macro_f1'] != manifest['best_validation_macro_f1']:
        raise ValueError('Checkpoint selection inconsistent with validation history')
    records = json.loads((run/'predictions.json').read_text())
    if len({r['measurement_id'] for r in records}) != len(records):
        raise ValueError('Duplicate prediction identity')
    if any(r['run_identity'] != manifest['identity'] for r in records):
        raise ValueError('Prediction run identity mismatch')
    if evaluate(run) != json.loads((run/'metrics.json').read_text()):
        raise ValueError('Saved metrics inconsistent with predictions')
    if predict(run/'checkpoint.pt', records) != [r['probabilities'] for r in records]:
        raise ValueError('Checkpoint reload predictions inconsistent')
    if prepared_path is not None:
        prepared = _prepared(prepared_path)
        if any(prepared[k] != manifest[k] for k in ['input_fingerprint', 'split_fingerprint', 'features']):
            raise ValueError('Incompatible prepared fingerprint')
        identity = manifest['identity']
        if identity not in prepared['expected_runs']:
            raise ValueError('Unexpected saved run identity')
        split = next(s for s in prepared['splits'] if s['setting'] == identity['setting'] and s['fold'] == identity['fold'])
        training_rows = [r for r, role in zip(prepared['rows'], split['roles']) if role == 'train'
                         and (r['condition'], r['phenotype']) == (identity['condition'], identity['phenotype'])]
        if manifest['dose_statistics'] != _statistics(training_rows):
            raise ValueError('Dose statistics differ from prepared training rows')
        expected = [r for r, role in zip(prepared['rows'], split['roles']) if role == 'test'
                    and (r['condition'], r['phenotype']) == (identity['condition'], identity['phenotype'])]
        actual = [{k: v for k, v in r.items() if k not in ['probabilities', 'prediction', 'run_identity']} for r in records]
        if actual != expected:
            raise ValueError('Prediction identities/provenance differ from prepared test partition')
    return True
