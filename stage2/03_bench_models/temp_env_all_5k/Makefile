PY ?= /home/chang/miniconda3/envs/arc/bin/python

.PHONY: dry-run test build-data gene-set fidelity compare

dry-run:
	bash scripts/run.sh --dry-run

test:
	$(PY) -m unittest discover -s scripts/tests -v

gene-set:
	$(PY) scripts/build_gene_set.py

build-data:
	$(PY) scripts/build_data.py

fidelity:
	$(PY) scripts/aggregate_fidelity.py

compare:
	$(PY) scripts/compare_18k.py
