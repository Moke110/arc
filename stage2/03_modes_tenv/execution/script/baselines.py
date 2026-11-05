"""Lightweight subprocess adapters for the two deterministic references."""
import argparse
import time
from pathlib import Path
import numpy as np
from artifacts import read_json, write_json, rows


def fit_predict(request):
    directory = Path(request['output_directory'])
    files = request['predictor_files']
    cases = rows(files['test_cases'])
    training = rows(files['train_cases'])
    n_genes = request['n_genes']
    if request['identity']['model'] == 'no_change':
        delta = np.zeros(n_genes, dtype=np.float32)
    elif request['identity']['model'] == 'delta_mean':
        responses = np.load(files['train_expression'], mmap_mode='r')
        groups = sorted({r['group_id'] for r in training})
        if not groups:
            raise ValueError('delta_mean requires training groups')
        delta = np.zeros(n_genes, dtype=np.float64)
        for group in groups:
            selected = [i for i,r in enumerate(training) if r['group_id'] == group]
            delta += np.mean(responses[selected],axis=0,dtype=np.float64) / len(groups)
        delta = delta.astype(np.float32)
    else:
        raise ValueError('unknown baseline')
    write_json(directory / 'progress.json',{'phase':'inference','training_seconds':0.})
    inference_start=time.monotonic()
    prediction = np.lib.format.open_memmap(directory / 'prediction.npy',mode='w+',dtype='float32',shape=(len(cases),n_genes))
    for start in range(0,len(cases),32):
        prediction[start:start+32] = delta
    prediction.flush()
    write_json(directory / 'case_ids.json',[r['case_id'] for r in cases])
    write_json(directory / 'checkpoint.json',{'identity':request['identity'], 'deterministic_reference':True,
               'training_seconds':0., 'inference_seconds':time.monotonic()-inference_start,
               'author_objective':None, 'selection':'no fitted checkpoint',
               'trainable_parameters':0, 'group_weighting':'equal training drug-dose groups'})

if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--request',required=True)
    fit_predict(read_json(parser.parse_args().request))
