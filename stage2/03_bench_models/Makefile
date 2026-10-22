PY ?= /home/chang/miniconda3/envs/arc/bin/python

.PHONY: dry-run test build-data fidelity

dry-run:
	bash run.sh --dry-run

test:
	$(PY) -m unittest discover -s tests -v

build-data:
	$(PY) build_data.py

fidelity:
	$(PY) aggregate_fidelity.py
