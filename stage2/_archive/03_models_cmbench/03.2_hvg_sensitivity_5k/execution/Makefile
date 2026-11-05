PY ?= /home/chang/miniconda3/envs/arc/bin/python

# Run from this directory (`<substep>/execution/`).  `test` is the offline
# suite; `gene-set` + `build-data` regenerate the frozen input under
# ../output/input/.
.PHONY: dry-run test gene-set build-data fidelity compare

dry-run:
	bash script/run.sh --dry-run

test:
	$(PY) -m unittest discover -s tests -v

gene-set:
	$(PY) script/build_gene_set.py

build-data:
	$(PY) script/build_data.py

fidelity:
	$(PY) script/aggregate_fidelity.py

compare:
	$(PY) script/compare_18k.py
