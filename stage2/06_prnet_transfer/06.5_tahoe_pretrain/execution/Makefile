# Stage 2 · 06.5 — seeded formal Tahoe pretrains (ticket #5).
#
# Run from this directory (`<substep>/execution/`):
#   make test           offline record/checker tests (arc env, no GPU, no data)
#   make plan           print the formal plan (commands + config digests), no run
#   make profile        short real-data profile (bounded smoke + one full epoch), 12G cap
#   make formal         the three frozen seeds, locked budget, sequential, 18G cap
#   make check-formal   validate output/formal/tahoe_pretrain_record.json
#   make check-profile  validate output/profile/profile_record.json
#
# The wrapper imports the 06.4 pipeline library; it does not copy it.

PY ?= /home/chang/miniconda3/envs/arc/bin/python
GPU_PY ?= /home/chang/miniconda3/envs/arc-stage2-bench-prnet/bin/python
ROOT := $(shell git -C $(CURDIR) rev-parse --show-toplevel)
DATA_DIR ?= $(ROOT)/stage2/06_prnet_transfer/06.3_input/output
OUT ?= $(ROOT)/stage2/06_prnet_transfer/06.5_tahoe_pretrain/output
CACHE ?= $(ROOT)/stage2/06_prnet_transfer/06.4_pipeline/output/cache/tahoe_log1p_cp10k_v1.npy
TABLE_CACHE ?= $(ROOT)/stage2/06_prnet_transfer/06.4_pipeline/output/input_cache
PYTHONDONTWRITEBYTECODE = 1
export PYTHONDONTWRITEBYTECODE

.PHONY: test plan profile formal check-formal check-profile

test:
	$(PY) -m unittest discover -s tests -v

plan:
	$(GPU_PY) run_pretrain.py --mode formal --dry-run --data-dir "$(DATA_DIR)" \
		--cache "$(CACHE)" --table-cache "$(TABLE_CACHE)"

profile:
	PY=$(GPU_PY) ./run.sh profile --data-dir "$(DATA_DIR)" --cache "$(CACHE)" \
		--table-cache "$(TABLE_CACHE)"

formal:
	PY=$(GPU_PY) ./run.sh formal --data-dir "$(DATA_DIR)" --cache "$(CACHE)" \
		--table-cache "$(TABLE_CACHE)"

check-formal:
	$(PY) check_record.py --record "$(OUT)/formal/tahoe_pretrain_record.json" --expect-mode formal

check-profile:
	$(PY) check_record.py --record "$(OUT)/profile/profile_record.json" --expect-mode profile
