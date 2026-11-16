#!/usr/bin/env bash
set -euo pipefail
ITEM_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
cd -- "$ITEM_DIR/source"
sha256sum -c ../manifest.sha256
