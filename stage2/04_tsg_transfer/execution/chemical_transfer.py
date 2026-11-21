"""Compare complete and chemical-only Tahoe handoffs in tenv development."""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

STEP = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(STEP / 'output/api'))
from tsg_transfer import chemical, development, fit, tenv, transfer

FOLD = 0
EPOCHS = 20
LR = 1e-4
SEEDS = (42, 17)
BINS = (('low', 0., .4), ('medium', .4, .6), ('high', .6, 1.01))
PARENT = 'drug_parent_inchikey'


def nearest_training(data, train, assessment):
    """Nearest nonoverlapping fit-parent similarity on frozen 2048-bit ECFP4."""
    parents = data.obs[PARENT].to_numpy()
    fit_parents = sorted(set(parents[train]))
    query = sorted(set(parents[assessment]))
    if set(fit_parents).intersection(query):
        raise ValueError('assessment and training parents overlap')
    lookup = data.parent_lookup
    reference = data.morgan[[lookup[key] for key in fit_parents]].astype(np.int16)
    targets = data.morgan[[lookup[key] for key in query]].astype(np.int16)
    intersection = targets @ reference.T
    union = targets.sum(axis=1)[:, None] + reference.sum(axis=1) - intersection
    if (union <= 0).any():
        raise ValueError('Tanimoto similarity requires nonempty fingerprint union')
    values = (intersection / union).max(axis=1)
    return dict(zip(query, values, strict=True))


def summarize_bins(scores, similarity):
    matched = scores[PARENT].map(similarity)
    result = {}
    for name, lower, upper in BINS:
        subset = scores[(matched >= lower) & (matched < upper)]
        metrics = development.summarize(subset)
        result[name] = {key: (value if np.isfinite(value) else None) for key, value in metrics.items()}
    return result


def run(args):
    data = fit.load_inputs(args.inputs, 'tenv')
    formal = tenv.load_split(data, args.experiment / 'tenv-split', FOLD)
    roles = development.split_validation(data, formal['train'], formal['validation'])
    tahoe, winner, source_doses = development.load_tahoe_winner(
        args.tahoe_inputs, args.experiment / 'tahoe', args.source,
    )
    held = tenv.DOSES[FOLD]
    chemical_keys = sorted('feat_embeddings.' + name for name in tahoe.feat_embeddings.state_dict())
    ownership = {name: dict(shape=list(tensor.shape), dtype=str(tensor.dtype))
                 for name, tensor in tahoe.state_dict().items() if name in chemical_keys}
    nearest = nearest_training(data, roles['train'], roles['assessment'])
    identities = {}
    for role, indices in roles.items():
        rows = data.obs.iloc[indices]
        identities[role] = dict(
            n_profiles=len(rows), n_parents=int(rows[PARENT].nunique()),
            ordered_profile_sha256=hashlib.sha256(
                json.dumps(list(rows.index), separators=(',', ':')).encode()).hexdigest(),
            doses_um=sorted(rows.dose_um.unique()),
        )

    def initialize_full(recipient, transform):
        return transfer.initialize_from_tahoe(tahoe, recipient, winner['transform'], transform,
                                             source_doses + tenv.DOSES, held_dose=held)

    def initialize_chemical(recipient, _):
        copied = chemical.initialize_chemical(tahoe, recipient)
        if sorted(copied) != chemical_keys:
            raise ValueError('chemical-only transfer differs from recorded ownership')

    results = []
    for seed in SEEDS:
        for arm, initializer in (('random', None), ('full', initialize_full),
                                 ('chemical-only', initialize_chemical)):
            print(f'seed={seed} arm={arm}', file=sys.stderr, flush=True)
            run = development.run_fit(data, args.source, roles, learning_rate=LR, seed=seed,
                                      epochs=EPOCHS, initialize=initializer, verify_reload=True)
            chosen = run['selected']['pcc']
            results.append(dict(seed=seed, arm=arm, selected=chosen['selected'],
                                assessment=chosen['assessment'], history=run['history'],
                                exact_reload_verified=chosen['exact_reload_verified'],
                                by_similarity=summarize_bins(chosen['scores'], nearest)))
    return dict(fold=FOLD, excluded_dose_um=held, epochs=EPOCHS, learning_rate=LR,
                seeds=SEEDS, identities=identities, chemical_state=ownership,
                nonchemical_keys=sorted(set(tahoe.state_dict()) - set(chemical_keys)),
                winner_sha256=winner['checkpoint_sha256'],
                split_sha256=fit.sha256(args.experiment / 'tenv-split' / 'manifest.json'),
                source_tenv_inputs_sha256=fit.sha256(args.inputs / 'manifest.json'),
                source_tahoe_inputs_sha256=fit.sha256(args.tahoe_inputs / 'manifest.json'),
                similarity_bins=BINS, assessment_parent_tanimoto=nearest, results=results)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('inputs', 'tahoe-inputs', 'experiment', 'source', 'output'):
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('--output must be new')
    fit.write_json(args.output, run(args))


if __name__ == '__main__':
    main()
