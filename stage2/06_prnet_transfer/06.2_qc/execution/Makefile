# stage2/06_prnet_transfer/06.2_qc — ticket #2 harness
#
# Data lives in the main checkout (pert_data/, stage2/02|04 output/ are gitignored),
# so DATA_ROOT is explicit. The tracked record is copied back with `make record`.

PYTHON   ?= /home/chang/miniconda3/envs/arc/bin/python
DATA_ROOT ?= /home/chang/arc
OUT_ROOT  ?= $(DATA_ROOT)
HERE      := $(patsubst %/,%,$(dir $(abspath $(lastword $(MAKEFILE_LIST)))))
RECORD    := $(abspath $(HERE)/../output)
OUT       := $(OUT_ROOT)/stage2/06_prnet_transfer/06.2_qc/output
TRACKED   := qc_rules.json qc_summary.json handoff.json manifest.json \
             counts_by_source_role.csv counts_by_plate_dose.csv \
             counts_tahoe_by_cell_line.csv exclusions_by_rule.csv \
             primary_reason_counts.csv plate_drug_confounding.csv \
             matched_control_availability.csv temp_env_dmso_controls.csv \
             tahoe_dmso_controls.csv

export PYTHONDONTWRITEBYTECODE = 1

.PHONY: all test lock build verify record status clean

all: lock build verify

test:
	cd "$(HERE)" && "$(PYTHON)" -m unittest discover -s tests -v

lock:
	cd "$(HERE)" && "$(PYTHON)" build_qc.py --lock --data-root "$(DATA_ROOT)" --out-root "$(OUT_ROOT)"

build:
	cd "$(HERE)" && "$(PYTHON)" build_qc.py --build --data-root "$(DATA_ROOT)" --out-root "$(OUT_ROOT)"

verify:
	cd "$(HERE)" && "$(PYTHON)" build_qc.py --verify --data-root "$(DATA_ROOT)" --out-root "$(OUT_ROOT)"

# copy the small tracked record into the step's own output/ (git add -f is manual)
record:
	mkdir -p "$(RECORD)"
	for f in $(TRACKED); do cp "$(OUT)/$$f" "$(RECORD)/$$f"; done
	@echo "copied $(words $(TRACKED)) files to $(RECORD); now: git add -f stage2/06_prnet_transfer/06.2_qc/output"

status:
	@echo "out:    $(OUT)"
	@echo "record: $(RECORD)"
	@ls -la "$(OUT)" 2>/dev/null || true

# removes only the QC output directory of this step
clean:
	rm -rf "$(OUT)"
