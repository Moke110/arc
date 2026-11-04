#!/usr/bin/env bash
# Read-only offline source validation via the shared archive contract.
set -euo pipefail
item_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
python3 "$item_dir/script/download.py" --verify
