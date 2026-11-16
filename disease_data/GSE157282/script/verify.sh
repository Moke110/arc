#!/usr/bin/env bash
set -euo pipefail
ITEM="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ITEM/source"
sha256sum -c ../manifest.sha256
