"""Run the unchanged CMonge adapter with its isolated pinned CUDA plugin overlay."""
import os
from pathlib import Path
import runpy
import sys
sys.dont_write_bytecode = True

ROOT = next(p for p in Path(__file__).resolve().parents if (p / '.git').exists())
sys.path.insert(0, str(ROOT / 'pert_model/CMonge/local/tenv_cuda'))
os.environ['JAX_PLATFORM_NAME'] = 'cuda'
os.environ['XLA_PYTHON_CLIENT_PREALLOCATE'] = 'false'
runpy.run_path(str(Path(__file__).with_name('cmonge_adapter.py')), run_name='__main__')
