PY ?= /home/chang/miniconda3/envs/arc/bin/python

# Run from this directory (`<substep>/execution/`).  `test` is the offline
# suite; `check` reads the frozen sources only; `build`/`verify` act on the
# main checkout's artifact tree by default (override OUT=/path/to/06_prnet_transfer).
# Set ARC_DATA_ROOT=/home/chang/arc when running from a tracked-files-only worktree.
OUT ?= $(shell cd ../.. && pwd)

.PHONY: test check build verify

test:
	PYTHONDONTWRITEBYTECODE=1 $(PY) -m unittest discover -s tests -v

check:
	$(PY) build_input.py --check

build:
	$(PY) build_input.py --out-root $(OUT)

stub:
	$(PY) build_input.py --stub-handoffs --skip-tahoe --out-root /tmp/06_3_stub

verify:
	$(PY) build_input.py --verify --out-root $(OUT)
