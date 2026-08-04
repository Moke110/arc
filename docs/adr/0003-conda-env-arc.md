# Conda environment renamed from xpert to arc

All subprojects now assume a single conda env named `arc` (Python 3.13), created from a root `environment.yml` and shared by every component and workflow. The env is renamed because the old name (`xpert`) referred to one component, which no longer reflects a shared environment; a root-owned spec also makes the environment reproducible, which it previously was not. Every `run.sh`, notebook, and config that referenced `conda activate xpert` is rewritten to `arc`.
