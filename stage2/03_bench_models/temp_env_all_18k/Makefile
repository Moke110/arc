PY ?= /home/chang/miniconda3/envs/arc/bin/python

.PHONY: dry-run test build-data fidelity

dry-run:
	bash scripts/run.sh --dry-run

test:
	$(PY) -m unittest discover -s scripts/tests -v

build-data:
	$(PY) scripts/build_data.py

fidelity:
	$(PY) scripts/aggregate_fidelity.py
