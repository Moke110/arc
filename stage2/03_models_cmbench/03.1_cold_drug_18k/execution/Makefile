PY ?= /home/chang/miniconda3/envs/arc/bin/python

# Run from this directory (`<substep>/execution/`).  `test` is the offline
# suite; `build-data` regenerates the frozen input under ../output/input/.
.PHONY: dry-run test build-data fidelity

dry-run:
	bash script/run.sh --dry-run

test:
	$(PY) -m unittest discover -s tests -v

build-data:
	$(PY) script/build_data.py

fidelity:
	$(PY) script/aggregate_fidelity.py
