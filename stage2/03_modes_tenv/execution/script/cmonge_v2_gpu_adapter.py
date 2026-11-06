"""Run the isolated CMonge v2 adapter with its isolated pinned CUDA plugin overlay."""
import os
from pathlib import Path
import runpy
import sys
sys.dont_write_bytecode = True

ROOT = next(p for p in Path(__file__).resolve().parents if (p / '.git').exists())
sys.path.insert(0, str(ROOT / 'pert_model/CMonge/local/tenv_cuda'))
cuda_tools = ROOT / 'pert_model/CMonge/local/tenv_cuda/nvidia/cuda_nvcc'
os.environ['PATH'] = str(cuda_tools / 'bin') + os.pathsep + os.environ.get('PATH', '')
os.environ['XLA_FLAGS'] = os.environ.get('XLA_FLAGS', '') + ' --xla_gpu_cuda_data_dir=' + str(cuda_tools)
os.environ['JAX_PLATFORM_NAME'] = 'cuda'
os.environ['XLA_PYTHON_CLIENT_PREALLOCATE'] = 'false'
runpy.run_path(str(Path(__file__).with_name('cmonge_v2_adapter.py')), run_name='__main__')
