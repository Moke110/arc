#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
(cd source && sha256sum --check ../manifest.sha256)
