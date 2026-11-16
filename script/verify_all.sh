#!/usr/bin/env bash
# verify_all.sh — run every external item's offline verifier.
#
# Why this exists
#   Under the provenance contract (ADR 0030) each external item owns
#   `script/verify.sh`, which checks the downloaded bytes against the item's
#   `manifest.sha256` without the network and without the original archive. This
#   wrapper gives later tickets one red/green answer across all five external
#   categories.
#
# Usage
#   script/verify_all.sh        run the check
#   script/verify_all.sh --help print this text
#
# What it does
#   1. Checks the core external items plus active ADR 0038 phenotype items. A missing
#      category/item/verifier, or an unknown extra item directory, is a failure.
#   2. An item without `manifest.sha256` is a failure, not a skip: every item is
#      expected to carry one (initial27/27 hashed2026-10-30; additions hashed
#      on admission), so a deleted manifest
#      must turn this check red. The message says how to regenerate it.
#   3. For every expected item it runs its verifier from the item directory
#      and collects the result. A verifier failure makes the whole check fail.
#   4. Prints a PASS/FAIL summary and the exit code.
#
# Exit status
#   0  all expected active items exist and passed
#   1  the inventory is incomplete/unexpected, or at least one verifier failed
#   2  the script could not locate the repository root
#
# The repository root is found by walking up from this script until a directory
# containing .git is found — never by counting `../..` segments.
#
# No network. Read-only; idempotent: re-running it changes nothing.

set -uo pipefail

usage() {
  sed -n '2,/^set -uo pipefail$/p' "${BASH_SOURCE[0]}" | sed '$d' | sed 's/^# \{0,1\}//'
}

find_repo_root() {
  local dir
  dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)" || return 1
  while [ -n "$dir" ] && [ "$dir" != "/" ]; do
    if [ -e "$dir/.git" ]; then
      printf '%s\n' "$dir"
      return 0
    fi
    dir="$(dirname -- "$dir")"
  done
  return 1
}

main() {
  case "${1:-}" in
    -h | --help)
      usage
      exit 0
      ;;
    '') ;;
    *)
      echo "ERROR: unknown argument: $1 (try --help)" >&2
      exit 2
      ;;
  esac

  local repo_root
  repo_root="$(find_repo_root)" || {
    echo "ERROR: cannot find the repository root (no .git directory above this script)" >&2
    echo "exit: 2" >&2
    exit 2
  }
  cd -- "$repo_root" || exit 2

  local tmpdir
  tmpdir="$(mktemp -d)" || exit 2
  trap 'rm -rf -- "${tmpdir:-}"' EXIT

  local categories=(pert_data pert_model disease_data drug_emb_model pheno_data)
  local expected_items=(
    pert_data/Annotation
    pert_data/DToxS-iPSC
    pert_data/GSE244740
    pert_data/GSE247541
    pert_data/GSE262419
    pert_data/LINCS-LDS-1593
    pert_data/Tahoe-100M
    pert_data/XPert-L1000
    pert_model/PRnet
    pert_model/CPA
    pert_model/chemCPA
    pert_model/biolord
    pert_model/PerturbNet
    pert_model/scVIDR
    pert_model/Squidiff
    pert_model/CMonge
    pert_model/DECANT
    pert_model/WAVE
    pert_model/ExPO
    pert_model/State
    pert_model/TranSiGen
    pert_model/XPert
    pert_model/PerturbDiff
    disease_data/GSE106297
    disease_data/GSE157282
    disease_data/GSE232331
    disease_data/GSE235470
    disease_data/GSE242692
    disease_data/GSE255505
    disease_data/GSE313746
    disease_data/GSE108940
    disease_data/GSE116250
    disease_data/GSE128188
    disease_data/GSE137293
    disease_data/GSE138252
    disease_data/GSE141910
    disease_data/GSE160516
    disease_data/GSE165303
    disease_data/GSE183852
    disease_data/GSE221945
    disease_data/GSE248187
    disease_data/GSE298023
    disease_data/GSE329393
    disease_data/GSE194297
    disease_data/GSE224129
    disease_data/GSE234141
    disease_data/GSE238242
    disease_data/GSE243674
    disease_data/GSE253223
    disease_data/GSE255612
    disease_data/GSE269882
    disease_data/GSE270788
    disease_data/GSE282618
    disease_data/GSE304575
    disease_data/GSE305648
    disease_data/GSE35108
    disease_data/Kuppe-MI
    disease_data/Reichart-CM
    drug_emb_model/HG
    drug_emb_model/KPGT
  )
  # ADR 0038 extends the legacy inventory with a declared phenotype registry.
  # Invalid or absent registry is fatal, including on otherwise intact clones.
  local phenotype_items
  phenotype_items="$(python3 - "$repo_root/pheno_data/registry.json" <<'REGISTRY'
import json, pathlib, re, sys
try:
    data = json.loads(pathlib.Path(sys.argv[1]).read_text())
    assert data["schema_version"] == 1 and isinstance(data["items"], list)
    names, ids = set(), set()
    for row in data["items"]:
        name, key, state = row["item"], row["dataset_id"], row["status"]
        assert re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]*", name)
        assert name not in names and key not in ids and not key.startswith("test")
        assert state in {"active", "candidate", "unavailable"}
        names.add(name); ids.add(key)
        if state == "active":
            print("pheno_data/" + name)
        elif (pathlib.Path(sys.argv[1]).parent / name).exists():
            raise ValueError("non-active item has a directory: " + name)
except (OSError, ValueError, KeyError, TypeError, AssertionError) as exc:
    print("invalid phenotype inventory: " + str(exc), file=sys.stderr)
    sys.exit(1)
REGISTRY
)" || { echo "FAIL: phenotype category missing or invalid inventory" >&2; return 1; }
  while IFS= read -r item; do
    [ -z "$item" ] || expected_items+=("$item")
  done <<<"$phenotype_items"
  local category item verifier manifest
  local total=0 passed=0 failed=0
  local failed_items="$tmpdir/failed.txt"
  : >"$failed_items"

  echo "== verify_all: offline manifests for every external item =="
  echo "repo root: $repo_root"
  echo "categories: ${categories[*]}"

  # Reject directories that are not part of ADR 0030's exact inventory. This
  # catches a missing item being accidentally replaced by a differently named
  # directory while preserving the total count.
  for category in "${categories[@]}"; do
    if [ ! -d "$category" ]; then
      failed=$((failed + 1))
      printf '%s\n' "$category" >>"$failed_items"
      echo "  FAIL  $category  (category missing)"
      continue
    fi
    for item in "$category"/*/; do
      item="${item%/}"
      [ -d "$item" ] || continue
      if ! printf '%s\n' "${expected_items[@]}" | grep -Fxq -- "$item"; then
        failed=$((failed + 1))
        printf '%s\n' "$item" >>"$failed_items"
        echo "  FAIL  $item  (unknown item directory)"
      fi
    done
  done

  for item in "${expected_items[@]}"; do
    total=$((total + 1))
    if [ ! -d "$item" ]; then
      failed=$((failed + 1))
      printf '%s\n' "$item" >>"$failed_items"
      echo "  FAIL  $item  (item missing)"
      continue
    fi
    verifier="$item/script/verify.sh"
    if [ ! -f "$verifier" ]; then
      failed=$((failed + 1))
      printf '%s\n' "$item" >>"$failed_items"
      echo "  FAIL  $item  (script/verify.sh missing)"
      continue
    fi
    manifest="$item/manifest.sha256"
    if [ ! -s "$manifest" ]; then
      failed=$((failed + 1))
      printf '%s\n' "$item" >>"$failed_items"
      echo "  FAIL  $item  (manifest.sha256 missing)"
      continue
    fi
    local outfile="$tmpdir/verify.out"
    local rc=0
    (cd -- "$item" && bash script/verify.sh) >"$outfile" 2>&1 || rc=$?
    if [ "$rc" -eq 0 ]; then
      passed=$((passed + 1))
      echo "  PASS  $item"
    else
      failed=$((failed + 1))
      printf '%s\n' "$item" >>"$failed_items"
      echo "  FAIL  $item  (verify.sh exit $rc)"
      sed -n '1,20p' "$outfile" | sed 's/^/        /'
    fi
  done

  echo "expected items: $total"

  if [ "$failed" -eq 0 ]; then
    echo
    echo "RESULT: PASS — verified $passed/$total item(s)"
    echo "exit: 0"
    return 0
  fi

  echo
  echo "RESULT: FAIL — $failed of $total item(s) failed verification"
  sed 's/^/  failed  /' "$failed_items"
  echo "exit: 1"
  return 1
}

main "$@"
