#!/usr/bin/env bash
set -euo pipefail
ITEM_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
FILE=44161_2022_28_MOESM4_ESM.xlsx
EXPECTED=4dc5a8197a855bd2394423c72b86076518233c510a480c2f98c2681bc9b554b1
mkdir -p "$ITEM_DIR/source"
if [[ -f "$ITEM_DIR/source/$FILE" ]]; then
  printf '%s  %s\n' "$EXPECTED" "$ITEM_DIR/source/$FILE" | sha256sum -c -
  exit 0
fi
curl --fail --silent --show-error --location --retry 3 \
  'https://static-content.springer.com/esm/art%3A10.1038%2Fs44161-022-00028-6/MediaObjects/44161_2022_28_MOESM4_ESM.xlsx' \
  --output "$ITEM_DIR/source/$FILE.part"
printf '%s  %s\n' "$EXPECTED" "$ITEM_DIR/source/$FILE.part" | sha256sum -c -
mv -- "$ITEM_DIR/source/$FILE.part" "$ITEM_DIR/source/$FILE"
