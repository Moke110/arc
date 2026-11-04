#!/usr/bin/env bash
set -euo pipefail
item_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
python3 "$item_dir/script/download.py" --verify
