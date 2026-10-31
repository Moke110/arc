#!/usr/bin/env bash
# verify_all.sh — run every external item's offline verifier.
#
# Why this exists
#   Under the provenance contract (ADR 0030) each external item owns
#   `script/verify.sh`, which checks the downloaded bytes against the item's
#   `manifest.sha256` without the network and without the original archive. This
#   wrapper gives later tickets one red/green answer across all four external
#   categories.
#
# Usage
#   script/verify_all.sh        run the check
#   script/verify_all.sh --help print this text
#
# What it does
#   1. Checks the exact 27-item inventory accepted by ADR 0030. A missing
#      category/item/verifier, or an unknown extra item directory, is a failure.
#   2. An item without `manifest.sha256` is a failure, not a skip: every item is
#      expected to carry one (27/27 hashed 2026-10-30), so a deleted manifest
#      must turn this check red. The message says how to regenerate it.
#   3. For every expected item it runs its verifier from the item directory
#      and collects the result. A verifier failure makes the whole check fail.
#   4. Prints a PASS/FAIL summary and the exit code.
#
# Exit status
#   0  all 27 expected items exist and passed
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

  local categories=(pert_data pert_model disease_data drug_emb_model)
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
    disease_data/GSE108940
    disease_data/GSE137293
    disease_data/GSE160516
    drug_emb_model/HG
    drug_emb_model/KPGT
  )
  local category item verifier manifest
  local total=0 passed=0 failed=0
  local failed_items="$tmpdir/failed.txt"
  : >"$failed_items"

  echo "== verify_all: offline manifests for every external item =="
  echo "repo root: $repo_root"
  echo "categories: ${categories[*]}"

  # Reject directories that are not part of ADR 0030's exact inventory. This
  # catches a missing item being accidentally replaced by a differently named
  # directory while preserving the total count of 27.
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
    if [ ! -f "$manifest" ]; then
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
