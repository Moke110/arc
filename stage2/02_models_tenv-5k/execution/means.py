"""Run deterministic fold Means from verified shared prepared inputs."""
import argparse
import json
from pathlib import Path
import sys
import numpy as np
import pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'output/api'))
from tenv5k.inputs import Inputs, sha256
from tenv5k.protocol import mean_prediction, score_prediction


def run(inputs, candidate):
    data = Inputs(inputs)
    candidate = Path(candidate)
    if candidate.exists() and any(candidate.iterdir()):
        raise ValueError('candidate must be empty')
    candidate.mkdir(parents=True, exist_ok=True)
    statuses = []
    for scenario in ('cold-drug', 'warm-start'):
        for fold in range(5):
            train = data.rows(scenario, fold, 'train')
            test = data.rows(scenario, fold, 'test')
            for method in ('Mean', 'Mean_cell', 'Mean_drug'):
                prediction = mean_prediction(data.expression, data.obs, data.controls, train, test, method, scenario)
                identity = {'scenario': scenario, 'fold': fold, 'method': method}
                if prediction is None:
                    statuses.append({**identity, 'status': 'not-applicable', 'reason': 'unseen test drug'})
                    continue
                target = candidate / scenario / str(fold) / method
                target.mkdir(parents=True)
                np.save(target / 'prediction.npy', prediction)
                data.obs.iloc[test][['profile_id', 'drug_parent_inchikey']].to_csv(target / 'identity.csv', index=False)
                wells, parents = score_prediction(data.obs.iloc[test], data.expression[test], prediction, data.control[test], prediction_ids=data.obs.index[test])
                wells.to_csv(target / 'per_well.csv', index=False)
                parents.to_csv(target / 'per_drug.csv', index=False)
                statuses.append({**identity, 'status': 'complete', 'n_test': len(test)})
    pd.DataFrame(statuses).to_csv(candidate / 'status.csv', index=False)
    outputs = {str(p.relative_to(candidate)): sha256(p) for p in sorted(candidate.rglob('*')) if p.is_file()}
    (candidate / 'manifest.json').write_text(json.dumps({'schema': 'tenv5k-means-v1', 'input_manifest_sha256': sha256(Path(inputs) / 'manifest.json'), 'outputs': outputs}, indent=2) + '\n')
    return statuses

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs', type=Path, required=True)
    parser.add_argument('--candidate', type=Path, required=True)
    args = parser.parse_args()
    run(args.inputs, args.candidate)
