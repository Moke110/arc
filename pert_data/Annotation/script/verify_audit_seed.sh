#!/usr/bin/env bash
# Verify the separately sealed local audit seed after installation.
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd -- "$script_dir/../local"
sha256sum -c "$script_dir/audit_seed.sha256"
