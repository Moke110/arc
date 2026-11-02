# Stage 2 · 06.8 — ticket #8 five-fold baselines.
#
# Data/frozen input and compute outputs live in the main checkout, not the
# worktree, so DATA_DIR / OUT_ROOT are explicit (06.2 precedent).  The tracked
# record is copied back into this step with `make record`; the five record files
# are tracked through a narrow `.gitignore` exception on `06.8_baselines/output/`.
#
#   make test         offline unittest suite (arc env, synthetic fixture)
#   make check        contract check against the real frozen input, no writes
#   make run          all five folds + the tracked record (pinned GPU env)
#   make record       copy the small tracked record into ./output/
#   make fixture      write the synthetic frozen input (test aid)
#   make dry-run      full record on the synthetic fixture (pipeline check only)
#   make mem-check    dry-run on the row-replicated fixture, under /usr/bin/time -v

PY ?= /home/chang/miniconda3/envs/arc/bin/python
GPU_PY ?= /home/chang/miniconda3/envs/arc-stage2-bench-prnet/bin/python
ROOT := $(shell git -C $(CURDIR) rev-parse --show-toplevel)
HERE := $(patsubst %/,%,$(dir $(abspath $(lastword $(MAKEFILE_LIST)))))
DATA_DIR ?= $(ROOT)/stage2/06_prnet_transfer/06.3_input/output
OUT_ROOT ?= $(ROOT)/stage2/06_prnet_transfer
TABLE_CACHE ?= $(ROOT)/stage2/06_prnet_transfer/06.4_pipeline/output/input_cache
OUT := $(OUT_ROOT)/06.8_baselines/output
RECORD := $(HERE)/../output
FIXTURE ?= /tmp/t08-baselines-fixture
STRESS ?= /tmp/t08-baselines-stress
STRESS_REPLICATES ?= 18
RECORD_FILES := summary_table.csv metrics.csv per_drug_dose.csv fold_coverage.csv \
                manifest.json
PYTHONDONTWRITEBYTECODE = 1
export PYTHONDONTWRITEBYTECODE

.PHONY: test check run record fixture dry-run stress-fixture mem-check tables clean

test:
	cd "$(HERE)" && $(PY) -m unittest discover -s tests -v

check:
	cd "$(HERE)" && $(PY) baselines.py --check --data-dir "$(DATA_DIR)" \
		--table-cache "$(TABLE_CACHE)"

run:
	cd "$(HERE)" && $(GPU_PY) baselines.py --all --data-dir "$(DATA_DIR)" \
		--out-root "$(OUT_ROOT)" --table-cache "$(TABLE_CACHE)"

tables:
	$(PY) $(ROOT)/stage2/06_prnet_transfer/06.4_pipeline/execution/tools/convert_frozen_tables.py \
		--data-dir "$(DATA_DIR)" --cache-dir "$(TABLE_CACHE)"

record:
	mkdir -p "$(RECORD)"
	for f in $(RECORD_FILES); do cp "$(OUT)/$$f" "$(RECORD)/$$f"; done

fixture:
	$(PY) $(ROOT)/stage2/06_prnet_transfer/06.4_pipeline/execution/tools/make_fixture.py \
		--out "$(FIXTURE)"

dry-run:
	cd "$(HERE)" && $(PY) baselines.py --all --data-dir "$(FIXTURE)" \
		--out-root /tmp/t08-baselines-dry-run --force

# Memory gate (ticket #8): the standard fixture is too small to show the
# row-level allocations, so replicate its rows to a real-data row count and
# measure the peak.  Wrap it in the WSL2 memory cap:
#   systemd-run --user --scope -p MemoryMax=4G -p MemorySwapMax=4G -- make mem-check
stress-fixture:
	$(PY) "$(HERE)/tools/replicate_fixture.py" --src "$(FIXTURE)" --out "$(STRESS)" \
		--replicates "$(STRESS_REPLICATES)"

mem-check: stress-fixture
	cd "$(HERE)" && /usr/bin/time -v $(PY) baselines.py --all --data-dir "$(STRESS)" \
		--out-root /tmp/t08-baselines-mem-check --force
