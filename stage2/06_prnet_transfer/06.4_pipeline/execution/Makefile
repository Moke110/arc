# Stage 2 · 06.4 — the PRnet formal-transfer pipeline (ticket #4).
#
# Run from this directory (`<substep>/execution/`):
#   make test        offline unittest suite (arc env, CPU, synthetic fixture)
#   make smoke-cpu   the CPU end-to-end smoke test alone
#   make gpu-smoke   Phase-B real-data smoke on the pinned GPU env
#   make fixture     write the synthetic frozen input (default /tmp location)
#   make tahoe-cache build/verify the derived 4.9 GB log1p memmap for #5
#   make check       dry-run the resolved config against the frozen input

PY ?= /home/chang/miniconda3/envs/arc/bin/python
GPU_PY ?= /home/chang/miniconda3/envs/arc-stage2-bench-prnet/bin/python
ROOT := $(shell git -C $(CURDIR) rev-parse --show-toplevel)
DATA_DIR ?= $(ROOT)/stage2/06_prnet_transfer/06.3_input/output
OUT ?= $(ROOT)/stage2/06_prnet_transfer/06.4_pipeline/output
CACHE ?= $(OUT)/cache/tahoe_log1p_cp10k_v1.npy
FIXTURE ?= /tmp/ptlib-fixture
PYTHONDONTWRITEBYTECODE = 1
export PYTHONDONTWRITEBYTECODE

.PHONY: test smoke-cpu gpu-smoke fixture frozen-tables tahoe-cache check

test:
	$(PY) -m unittest discover -s tests -v

smoke-cpu:
	$(PY) -m unittest tests.test_end_to_end -v

gpu-smoke:
	$(GPU_PY) smoke_real.py --data-dir "$(DATA_DIR)" --tahoe-cache "$(CACHE)" --build-cache

fixture:
	$(PY) tools/make_fixture.py --out "$(FIXTURE)"

frozen-tables:
	$(PY) tools/convert_frozen_tables.py --data-dir "$(DATA_DIR)" --cache-dir "$(OUT)/input_cache"

tahoe-cache:
	$(GPU_PY) tools/build_cache.py --data-dir "$(DATA_DIR)" --cache "$(CACHE)"

check:
	$(GPU_PY) train.py --stage pretrain --check --out "$(OUT)/_check" --data-dir "$(DATA_DIR)"
