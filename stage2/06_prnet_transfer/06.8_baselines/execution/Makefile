# Stage 2 · 06.8 — ticket #8 five-fold baselines.
#
# Data/frozen input and compute outputs live in the main checkout, not the
# worktree, so DATA_DIR / OUT_ROOT are explicit (06.2 precedent).  The tracked
# record is copied back into this step with `make record` (then `git add -f`).
#
#   make test         offline unittest suite (arc env, synthetic fixture)
#   make check        contract check against the real frozen input, no writes
#   make run          all five folds + the tracked record (pinned GPU env)
#   make record       copy the small tracked record into ./output/
#   make fixture      write the synthetic frozen input (test aid)
#   make dry-run      full record on the synthetic fixture (pipeline check only)

PY ?= /home/chang/miniconda3/envs/arc/bin/python
GPU_PY ?= /home/chang/miniconda3/envs/arc-stage2-bench-prnet/bin/python
ROOT := $(shell git -C $(CURDIR) rev-parse --show-toplevel)
HERE := $(patsubst %/,%,$(dir $(abspath $(lastword $(MAKEFILE_LIST)))))
DATA_DIR ?= /home/chang/arc/stage2/06_prnet_transfer/06.3_input/output
OUT_ROOT ?= /home/chang/arc/stage2/06_prnet_transfer
TABLE_CACHE ?= /home/chang/arc/stage2/06_prnet_transfer/06.4_pipeline/output/input_cache
OUT := $(OUT_ROOT)/06.8_baselines/output
RECORD := $(HERE)/../output
FIXTURE ?= /tmp/t08-baselines-fixture
RECORD_FILES := summary_table.csv metrics.csv per_drug_dose.csv fold_coverage.csv \
                manifest.json
PYTHONDONTWRITEBYTECODE = 1
export PYTHONDONTWRITEBYTECODE

.PHONY: test check run record fixture dry-run tables clean

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
