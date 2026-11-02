# Stage 2 · 06.6 — 45 random-init CM candidates and the per-fold LR selection (ticket #6).
#
# Run from this directory (`stage2/06_prnet_transfer/06.6_cm_random/execution/`):
#   make test        offline unittest suite (arc env, CPU, synthetic run dirs)
#   make tables      derived CSV table cache for the pinned env (arc env, ~1 s)
#   make queue       write output/units.json (45 units) from the frozen locked config
#   make plan        ask the ticket-#4 queue what it would do (resume-aware dry run)
#   make matrix      run all 45 candidates under the systemd-run 12 GiB scope, sampled
#   make record      derive output/record/* from output/runs/*
#   make check       re-derive + gate the record (one-to-one, selection, no test read)
#   make resume      prove resume: 45/45 skipped when complete, interrupted unit re-runs
#
# The canonical heavy command is `make matrix` (execution/run.sh scoped); every
# path is overridable.  Never run two of these at once.

PY ?= /home/chang/miniconda3/envs/arc/bin/python
GPU_PY ?= /home/chang/miniconda3/envs/arc-stage2-bench-prnet/bin/python
ROOT := $(shell git -C $(CURDIR) rev-parse --show-toplevel)
# Ticket #3's frozen input is 16.4 GB and gitignored: it exists once, in the
# integration tree.  A ticket worktree carries only the small tracked members, so
# fall back to the canonical frozen input when this checkout has no CM view.
INTEGRATION_DIR ?= /home/chang/arc
DATA_DIR ?= $(if $(wildcard $(ROOT)/stage2/06_prnet_transfer/06.3_input/output/cm_transfer.h5ad),$(ROOT)/stage2/06_prnet_transfer/06.3_input/output,$(INTEGRATION_DIR)/stage2/06_prnet_transfer/06.3_input/output)
OUT ?= $(ROOT)/stage2/06_prnet_transfer/06.6_cm_random/output
PIPELINE := $(ROOT)/stage2/06_prnet_transfer/06.4_pipeline/execution
TABLE_CACHE ?= $(OUT)/input_cache
WORKERS ?= 2
OMP_THREADS ?= 2
SCOPE := systemd-run --user --scope -p MemoryMax=12G -p MemorySwapMax=12G
PYTHONDONTWRITEBYTECODE = 1
export PYTHONDONTWRITEBYTECODE

.PHONY: test tables queue plan matrix record check resume smoke check-input guard-data

guard-data:
	@test -f "$(DATA_DIR)/cm_transfer.h5ad" || { \
		echo "error: no frozen CM view under DATA_DIR=$(DATA_DIR)"; \
		echo "       pass DATA_DIR=<integration>/stage2/06_prnet_transfer/06.3_input/output"; \
		exit 1; }

test:
	$(PY) -m unittest discover -s tests -v

tables: guard-data
	$(PY) $(PIPELINE)/tools/convert_frozen_tables.py --data-dir "$(DATA_DIR)" --cache-dir "$(TABLE_CACHE)"

queue: guard-data
	$(PY) build_queue.py --data-dir "$(DATA_DIR)" --out-dir "$(OUT)"

plan: guard-data
	$(GPU_PY) run_matrix.py --data-dir "$(DATA_DIR)" --out-dir "$(OUT)" \
		--table-cache "$(TABLE_CACHE)" --workers "$(WORKERS)" --dry-run

matrix: guard-data
	$(SCOPE) -- $(GPU_PY) run_matrix.py --data-dir "$(DATA_DIR)" --out-dir "$(OUT)" \
		--table-cache "$(TABLE_CACHE)" --workers "$(WORKERS)" \
		--omp-threads "$(OMP_THREADS)" --scope "$(SCOPE)"

record:
	$(PY) write_record.py --data-dir "$(DATA_DIR)" --out-dir "$(OUT)"

check:
	$(PY) write_record.py --check --data-dir "$(DATA_DIR)" --out-dir "$(OUT)"

resume: guard-data
	$(SCOPE) -- $(GPU_PY) run_matrix.py --verify-resume --data-dir "$(DATA_DIR)" \
		--out-dir "$(OUT)" --table-cache "$(TABLE_CACHE)"

# frozen-input + resolved-config contract check through the ticket-#4 CLI (no training)
check-input: guard-data
	$(GPU_PY) $(PIPELINE)/train.py --stage finetune --fold 0 --check \
		--out /tmp/cm_random_check --data-dir "$(DATA_DIR)" --table-cache "$(TABLE_CACHE)"
