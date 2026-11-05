#!/usr/bin/env bash
set -euo pipefail
python "$(dirname -- "${BASH_SOURCE[0]}")/download.py" --verify
