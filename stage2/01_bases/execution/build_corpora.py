"""Build/replay the indexed training release without changing native sources."""

import argparse
import json
from pathlib import Path
from bases.corpora import publish
from bases import raw

root = next(p for p in Path(__file__).resolve().parents if (p / ".git").exists())
p = argparse.ArgumentParser(description=__doc__)
raw.add_binding_arguments(p)
p.add_argument(
    "--shared",
    type=Path,
    default=root / "stage2/01_bases/output/shared-v2",
)
p.add_argument(
    "--drugs",
    type=Path,
    default=root / "pert_data/Tahoe-100M/source/metadata/drug_metadata.parquet",
)
p.add_argument(
    "--output", type=Path, default=root / "stage2/01_bases/output/training-v2"
)
p.add_argument("--verify", action="store_true")
p.add_argument("--record", type=Path)
a = p.parse_args()
result = publish(raw.cli_bindings(a,root), a.shared, a.drugs, a.output, verify=a.verify)
serialized = json.dumps(result, indent=2, sort_keys=True) + "\n"
if a.record:
    a.record.parent.mkdir(parents=True, exist_ok=True)
    a.record.write_text(serialized)
print(serialized)
