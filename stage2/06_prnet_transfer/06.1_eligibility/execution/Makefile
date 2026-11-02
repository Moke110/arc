PY ?= /home/chang/miniconda3/envs/arc/bin/python
ROOT := $(shell git rev-parse --show-toplevel)
EXEC := stage2/06_prnet_transfer/06.1_eligibility/execution

.PHONY: all build offline check test

all: build check test

build:
	cd $(ROOT) && $(PY) $(EXEC)/build_eligibility.py

offline:
	cd $(ROOT) && $(PY) $(EXEC)/build_eligibility.py --offline

check:
	cd $(ROOT) && $(PY) $(EXEC)/build_eligibility.py --check

# The tests import the execution modules directly (``import rules``), so run
# unittest from ``execution/`` and discover the ``tests/`` subdirectory; running
# from the repo root only worked by accident through an import side effect.
test:
	cd $(ROOT)/$(EXEC) && PYTHONDONTWRITEBYTECODE=1 $(PY) -m unittest discover -s tests -v
