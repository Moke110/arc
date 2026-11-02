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

test:
	cd $(ROOT) && PYTHONDONTWRITEBYTECODE=1 $(PY) -m unittest discover -s $(EXEC)/tests -v
