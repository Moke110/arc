# Stage 2 · 06.9 — ticket #9: evaluate the 30 selected checkpoints and write the
# first-round report.
#
# Frozen input, the #6/#7 run trees and #8's record live in the main checkout
# (worktrees are removed after merge), so DATA_DIR / TABLE_CACHE are explicit
# (06.2/06.8 precedent); OUT is this substep's output in the current tree, and
# the tracked record is `output/record/` (force-added).
#
#   make test              offline unittest suite (arc env, no data, no GPU)
#   make units             verify + write the 30-unit manifest (fails closed until #7)
#   make units-random      harness dry run: random branch only
#   make plan              print the 30-unit scoring plan, score nothing
#   make score-valid       harness smoke: predict+evaluate on the VALIDATION role only
#   make score-test        THE test run: requires CONFIRM_FROZEN_SELECTION=1
#   make score-test-scoped the same under a systemd-run 8 GiB memory cap (preferred)
#   make record            aggregate -> output/record (tables, summary, gates)
#   make record-valid      aggregate the validation smoke -> output/dry_run (labelled)
#   make cards             30 model cards -> output/model_cards (test record)
#   make report            results section + README block from the record
#   make check             re-derive the record and byte-compare (no test read)
#   make all               units -> score-test -> record -> cards -> report -> check

PY ?= /home/chang/miniconda3/envs/arc/bin/python
GPU_PY ?= /home/chang/miniconda3/envs/arc-stage2-bench-prnet/bin/python
ROOT := $(shell git -C $(CURDIR) rev-parse --show-toplevel)
HERE := $(patsubst %/,%,$(dir $(abspath $(lastword $(MAKEFILE_LIST)))))
STEP := $(ROOT)/stage2/06_prnet_transfer/06.9_evaluation
DATA_DIR ?= /home/chang/arc/stage2/06_prnet_transfer/06.3_input/output
TABLE_CACHE ?= /home/chang/arc/stage2/06_prnet_transfer/06.4_pipeline/output/input_cache
BASELINES_DIR ?= $(ROOT)/stage2/06_prnet_transfer/06.8_baselines/output
RANDOM_RECORD ?= $(ROOT)/stage2/06_prnet_transfer/06.6_cm_random/output/record/selected_checkpoints.csv
TAHOE_RECORD ?= $(ROOT)/stage2/06_prnet_transfer/06.7_cm_tahoe/output/record/selected_checkpoints.csv
OUT ?= $(STEP)/output
UNITS ?= $(OUT)/record/units.json
HASH_CACHE ?= $(OUT)/hash_cache.json
SCOPE ?= systemd-run --user --scope -p MemoryMax=8G -p MemorySwapMax=8G --
PYTHONDONTWRITEBYTECODE = 1
export PYTHONDONTWRITEBYTECODE

SCORE_ARGS := --units "$(UNITS)" --out "$(OUT)" --data-dir "$(DATA_DIR)" \
	--table-cache "$(TABLE_CACHE)" --hash-cache "$(HASH_CACHE)" --gpu-py "$(GPU_PY)"
UNITS_ARGS := --data-dir "$(DATA_DIR)" --baselines-dir "$(BASELINES_DIR)" \
	--random-record "$(RANDOM_RECORD)" --tahoe-record "$(TAHOE_RECORD)" --out "$(UNITS)"
RECORD_ARGS := --units "$(UNITS)" --out "$(OUT)" --data-dir "$(DATA_DIR)" \
	--baselines-dir "$(BASELINES_DIR)"

.PHONY: test units units-random plan score-valid score-test score-test-scoped \
	record record-valid cards report check all tables

test:
	cd "$(HERE)" && $(PY) -m unittest discover -s tests -v

units:
	cd "$(HERE)" && $(PY) build_units.py $(UNITS_ARGS)

units-random:
	cd "$(HERE)" && $(PY) build_units.py $(UNITS_ARGS) --allow-incomplete

plan:
	cd "$(HERE)" && $(PY) score.py $(SCORE_ARGS) --phase test --plan --confirm-frozen-selection

score-valid:
	cd "$(HERE)" && $(PY) score.py $(SCORE_ARGS) --phase valid --limit 2

score-test:
	@test "$(CONFIRM_FROZEN_SELECTION)" = "1" || \
		(echo "refusing: run 'make score-test CONFIRM_FROZEN_SELECTION=1' only after the" \
		 "Lead confirms #6/#7 selections are final and #5's sweep is done" >&2; exit 4)
	cd "$(HERE)" && $(PY) score.py $(SCORE_ARGS) --phase test --confirm-frozen-selection

score-test-scoped:
	$(SCOPE) $(MAKE) -C "$(HERE)" score-test CONFIRM_FROZEN_SELECTION=1

record:
	cd "$(HERE)" && $(PY) write_record.py $(RECORD_ARGS) --phase test

record-valid:
	cd "$(HERE)" && $(PY) write_record.py $(RECORD_ARGS) --phase valid --allow-subset --no-figures

cards:
	cd "$(HERE)" && $(PY) make_cards.py --units "$(UNITS)" --out "$(OUT)" --phase test

report:
	cd "$(HERE)" && $(PY) report.py --record-dir "$(OUT)/record" --write-readme

check:
	cd "$(HERE)" && $(PY) write_record.py $(RECORD_ARGS) --phase test --check

all: units score-test-scoped record cards report check

tables:
	$(PY) $(ROOT)/stage2/06_prnet_transfer/06.4_pipeline/execution/tools/convert_frozen_tables.py \
		--data-dir "$(DATA_DIR)" --cache-dir "$(TABLE_CACHE)"
