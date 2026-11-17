"""Fresh process, frozen-runtime replay that never writes original evidence."""
import argparse
import contextlib
import hashlib
import json
import os
import tempfile
from pathlib import Path
import platform
import sys


def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1048576),b''):h.update(block)
    return h.hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('runtime','inputs','run','source'):parser.add_argument('--'+name,type=Path,required=True)
    args=parser.parse_args()
    # Native third-party imports need writable caches even when the frozen
    # runtime and environment installation are read-only. No source is changed.
    caches=tempfile.TemporaryDirectory(prefix='tenv5k-import-cache-')
    os.environ['NUMBA_CACHE_DIR']=str(Path(caches.name)/'numba')
    os.environ['MPLCONFIGDIR']=str(Path(caches.name)/'matplotlib')
    sys.dont_write_bytecode=True
    sys.path.insert(0,str(args.runtime))
    import numpy as np
    import torch
    from tenv5k.learner import replay
    torch.set_num_threads(2)
    config=json.loads((args.run/'config.json').read_text())
    saved=torch.load(args.run/'checkpoint.pt',map_location='cpu',weights_only=False)
    if saved['config']!=config:raise ValueError('original checkpoint/config mismatch')
    environment=json.loads((args.run/'environment.json').read_text())
    actual=dict(python=platform.python_version(),torch=torch.__version__,numpy=np.__version__,
                device=config['device'],cuda=torch.version.cuda)
    if environment!=actual:raise ValueError('original replay environment changed')
    with contextlib.redirect_stdout(sys.stderr):
        prediction=replay(args.inputs,args.run/'checkpoint.pt',args.source,adapter='prnet')
    recorded=np.load(args.run/'prediction.npy')
    if prediction.dtype!=recorded.dtype or not np.array_equal(prediction,recorded):
        raise ValueError('original frozen-runtime sampled replay differs')
    print(json.dumps(dict(status='exact',prediction_sha256=digest(args.run/'prediction.npy'),
        checkpoint_sha256=digest(args.run/'checkpoint.pt'),input_manifest_sha256=digest(args.inputs/'manifest.json'),
        config_sha256=digest(args.run/'config.json'),environment_sha256=digest(args.run/'environment.json'),
        environment=actual,loaded_learner=str(Path(sys.modules['tenv5k.learner'].__file__).resolve()))))


if __name__=='__main__':main()
