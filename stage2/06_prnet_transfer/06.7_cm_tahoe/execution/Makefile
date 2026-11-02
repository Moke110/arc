# Stage 2 · 06.7 — 45 Tahoe-init CM candidates and the per-fold LR selection (ticket #7).
#
# Run from this directory (`stage2/06_prnet_transfer/06.7_cm_tahoe/execution/`):
#   make test        offline unittest suite (arc env, CPU, synthetic checkpoints/runs)
#   make tables      derived CSV table cache for the pinned env (arc env, ~1 s)
#   make init        verify the three #5 formal pretrain checkpoints (warm-start proof)
#   make queue       write output/units.json (45 units) from the frozen locked config
#   make plan        ask the ticket-#4 queue what it would do (resume-aware dry run)
#   make matrix      run all 45 candidates under the systemd-run 12 GiB scope, sampled
#   make profile     bounded one-unit wrapper profile (2 epochs) against the first
#                    available #5 checkpoint; never the matrix
#   make record      derive output/record/* from output/runs/*
#   make check       re-derive + gate the record (warm start, one-to-one, no test read)
#   make resume      prove resume: 45/45 skipped when complete, interrupted unit re-runs
#
# The canonical heavy command is `make matrix` (execution/run.sh scoped); every
# path is overridable.  Never run two heavy jobs at once: the matrix must not start
# until ticket #5's formal sweep has finished and the Lead has signalled it.

PY ?= /home/chang/miniconda3/envs/arc/bin/python
GPU_PY ?= /home/chang/miniconda3/envs/arc-stage2-bench-prnet/bin/python
ROOT := $(shell git -C $(CURDIR) rev-parse --show-toplevel)
# Ticket #3's frozen input is 16.4 GB and gitignored: it exists once, in the
# integration tree.  A ticket worktree carries only the small tracked members, so
# fall back to the canonical frozen input when this checkout has no CM view.
INTEGRATION_DIR ?= /home/chang/arc
DATA_DIR ?= $(if $(wildcard $(ROOT)/stage2/06_prnet_transfer/06.3_input/output/cm_transfer.h5ad),$(ROOT)/stage2/06_prnet_transfer/06.3_input/output,$(INTEGRATION_DIR)/stage2/06_prnet_transfer/06.3_input/output)
# Ticket #5's formal pretrain checkpoints (gitignored *.pt) live in the Lead-run
# integration tree until they are copied; a local formal/runs tree wins.
PRETRAIN_DIR ?= $(if $(wildcard $(ROOT)/stage2/06_prnet_transfer/06.5_tahoe_pretrain/output/formal/runs/tahoe_seed20261021/checkpoint/model.pt),$(ROOT)/stage2/06_prnet_transfer/06.5_tahoe_pretrain/output/formal/runs,$(INTEGRATION_DIR)/stage2/06_prnet_transfer/06.5_tahoe_pretrain/output/formal/runs)
OUT ?= $(ROOT)/stage2/06_prnet_transfer/06.7_cm_tahoe/output
PIPELINE := $(ROOT)/stage2/06_prnet_transfer/06.4_pipeline/execution
TABLE_CACHE ?= $(OUT)/input_cache
WORKERS ?= 2
OMP_THREADS ?= 2
PROFILE_EPOCHS ?= 2
SCOPE := systemd-run --user --scope -p MemoryMax=12G -p MemorySwapMax=12G
PYTHONDONTWRITEBYTECODE = 1
export PYTHONDONTWRITEBYTECODE

.PHONY: test tables init queue plan matrix profile record check resume check-input guard-data guard-init

guard-data:
	@test -f "$(DATA_DIR)/cm_transfer.h5ad" || { \
		echo "error: no frozen CM view under DATA_DIR=$(DATA_DIR)"; \
		echo "       pass DATA_DIR=<integration>/stage2/06_prnet_transfer/06.3_input/output"; \
		exit 1; }

guard-init:
	@test -f "$(PRETRAIN_DIR)/tahoe_seed20261021/checkpoint/model.pt" || { \
		echo "error: no ticket-#5 formal checkpoint under PRETRAIN_DIR=$(PRETRAIN_DIR)"; \
		echo "       the formal sweep must finish before the ticket-#7 matrix"; \
		exit 1; }

test:
	$(PY) -m unittest discover -s tests -v

tables: guard-data
	$(PY) $(PIPELINE)/tools/convert_frozen_tables.py --data-dir "$(DATA_DIR)" --cache-dir "$(TABLE_CACHE)"

init: guard-data guard-init
	$(GPU_PY) check_init.py --data-dir "$(DATA_DIR)" --checkpoint-root "$(PRETRAIN_DIR)" \
		--out "$(OUT)/init/tahoe_init.json"

queue: guard-data guard-init
	$(GPU_PY) build_queue.py --data-dir "$(DATA_DIR)" --checkpoint-root "$(PRETRAIN_DIR)" \
		--out-dir "$(OUT)"

plan: guard-data guard-init
	$(GPU_PY) run_matrix.py --data-dir "$(DATA_DIR)" --checkpoint-root "$(PRETRAIN_DIR)" \
		--out-dir "$(OUT)" --table-cache "$(TABLE_CACHE)" --workers "$(WORKERS)" --dry-run

matrix: guard-data guard-init
	$(SCOPE) -- $(GPU_PY) run_matrix.py --data-dir "$(DATA_DIR)" \
		--checkpoint-root "$(PRETRAIN_DIR)" --out-dir "$(OUT)" \
		--table-cache "$(TABLE_CACHE)" --workers "$(WORKERS)" \
		--omp-threads "$(OMP_THREADS)" --scope "$(SCOPE)"

profile: guard-data guard-init
	$(SCOPE) -- $(GPU_PY) run_matrix.py --profile --profile-epochs "$(PROFILE_EPOCHS)" \
		--data-dir "$(DATA_DIR)" --checkpoint-root "$(PRETRAIN_DIR)" --out-dir "$(OUT)" \
		--table-cache "$(TABLE_CACHE)" --hash-cache "$(OUT)/hash_cache.json" \
		--workers "$(WORKERS)" --omp-threads "$(OMP_THREADS)" --scope "$(SCOPE)"

record:
	$(PY) write_record.py --data-dir "$(DATA_DIR)" --out-dir "$(OUT)"

check:
	$(PY) write_record.py --check --data-dir "$(DATA_DIR)" --out-dir "$(OUT)"

resume: guard-data guard-init
	$(SCOPE) -- $(GPU_PY) run_matrix.py --verify-resume --data-dir "$(DATA_DIR)" \
		--checkpoint-root "$(PRETRAIN_DIR)" --out-dir "$(OUT)" --table-cache "$(TABLE_CACHE)"

# frozen-input + resolved-config contract check through the ticket-#4 CLI (no training)
check-input: guard-data
	$(GPU_PY) $(PIPELINE)/train.py --stage finetune --fold 0 --check \
		--out /tmp/cm_tahoe_check --data-dir "$(DATA_DIR)" --table-cache "$(TABLE_CACHE)"
