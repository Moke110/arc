#!/usr/bin/env bash
set -euo pipefail
item_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
cd -- "$item_dir/source"
sha256sum -c ../manifest.sha256
