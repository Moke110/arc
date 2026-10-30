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
#   1. Walks the four external categories: pert_data/, pert_model/,
#      disease_data/, drug_emb_model/ — one directory per item, one level deep.
#   2. An item participates when it has `script/verify.sh`.
#   3. An item without `manifest.sha256` is skipped and counted; the manifests are
#      produced by ticket 09, so this is expected until then.
#   4. For every participating item it runs its verifier from the item directory
#      and collects the result. A verifier failure makes the whole check fail.
#   5. Prints a PASS/FAIL summary and the exit code.
#
# Exit status
#   0  every item that has a manifest verified, or there was nothing to verify
#      yet (manifests pending); a clear message says which
#   1  at least one item's verifier failed
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
  local category item verifier manifest
  local total=0 passed=0 failed=0 pending=0
  local pending_items="$tmpdir/pending.txt"
  local failed_items="$tmpdir/failed.txt"
  : >"$pending_items"
  : >"$failed_items"

  echo "== verify_all: offline manifests for every external item =="
  echo "repo root: $repo_root"
  echo "categories: ${categories[*]}"

  for category in "${categories[@]}"; do
    [ -d "$category" ] || continue
    for item in "$category"/*/; do
      item="${item%/}"
      [ -d "$item" ] || continue
      verifier="$item/script/verify.sh"
      [ -f "$verifier" ] || continue
      total=$((total + 1))
      manifest="$item/manifest.sha256"
      if [ ! -f "$manifest" ]; then
        pending=$((pending + 1))
        printf '%s\n' "$item" >>"$pending_items"
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
  done

  echo "items with script/verify.sh: $total"
  if [ "$pending" -gt 0 ]; then
    echo "manifests pending: $pending 项 (skipped — manifest.sha256 is produced by ticket 09)"
    sed 's/^/  pending  /' "$pending_items"
  else
    echo "manifests pending: 0 项"
  fi

  if [ "$total" -eq 0 ]; then
    echo "note: no item has script/verify.sh yet — the item layout lands in tickets 03–05,"
    echo "      and the manifests in ticket 09."
    echo
    echo "RESULT: PASS — nothing to verify yet"
    echo "exit: 0"
    return 0
  fi

  if [ "$failed" -eq 0 ]; then
    echo
    echo "RESULT: PASS — verified $passed/$total participating item(s), $pending pending"
    echo "exit: 0"
    return 0
  fi

  echo
  echo "RESULT: FAIL — $failed of $total participating item(s) failed verification"
  sed 's/^/  failed  /' "$failed_items"
  echo "exit: 1"
  return 1
}

main "$@"
