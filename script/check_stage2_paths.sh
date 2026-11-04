#!/usr/bin/env bash
# check_stage2_paths.sh — assert stage2/ executable code has no retired step paths.
#
# Why this exists
#   Ticket #15/#14 folded `stage2/02_cm_data/` and `stage2/04_tahoe_data/` into
#   `stage2/02_dataset_corpus/` (ADR 0033).  A live Python/bash/R file that still
#   points at either retired step is a dead path waiting to fail at runtime.  The
#   `datasets-1` bundle name is retired too, but only as a *non-frozen* path
#   component: the frozen v1 evidence under
#   `stage2/_frozen/cm-data-v1/datasets-1/` keeps its historical name on purpose.
#
# Usage
#   script/check_stage2_paths.sh          scan and report
#   script/check_stage2_paths.sh --help   print this text
#
# What it does
#   1. Takes the tracked file list from `git ls-files stage2` and keeps the
#      executable-code suffixes `*.py`, `*.sh`, `*.R`, excluding the archival
#      trees `stage2/_frozen/` and `stage2/_archive/`.
#   2. Flags `stage2/02_cm_data/...` and `stage2/04_tahoe_data/...` literals.
#   3. Flags `datasets-1` except when it is bound to the frozen root: either the
#      literal contains `_frozen/cm-data-v1/datasets-1` or the `datasets-1`
#      component sits on a `FROZEN...`-rooted expression (e.g.
#      `FROZEN / "datasets-1/..."`).
#   4. Exempts recorded-evidence payload strings via
#      `script/stage2_path_whitelist.txt`: one `<file-glob>:<literal-ERE>` pair
#      (or a bare `<file-glob>`) per line, `#` comments.  An entry matches only
#      the named file and the exact literal; it is not a blanket file skip.
#   5. Reports counts and the location of every stale literal.
#
# Exit status
#   0  no stale literal
#   1  at least one stale literal outside the whitelist
#   2  the script could not locate the repository root
#
# The repository root is found by walking up from this script until a directory
# containing .git is found — never by counting `../..` segments (ADR 0030 §4).
# No network. Read-only; idempotent: repeated runs give identical results.

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

# A whitelist entry is `<file-glob>:<literal-ERE>` (or a bare `<file-glob>`).
# Both halves are matched against the repo-relative source path and the exact
# extracted literal, so an entry can never exempt an unrelated future string.
is_whitelisted() {
  local file="$1" literal="$2" line
  [ -f "$WHITELIST" ] || return 1
  while IFS= read -r line || [ -n "$line" ]; do
    line="${line%$'\r'}"
    case "$line" in
      '' | '#'*) continue ;;
    esac
    if printf '%s' "$line" | grep -q ':'; then
      local file_pat="${line%%:*}" literal_pat="${line#*:}"
      # shellcheck disable=SC2254  # glob match is intended here
      case "$file" in
        $file_pat) ;;
        *) continue ;;
      esac
      if printf '%s' "$literal" | grep -qE -e "$literal_pat" 2>/dev/null; then
        return 0
      fi
    else
      # shellcheck disable=SC2254
      case "$file" in
        $line) return 0 ;;
      esac
    fi
  done <"$WHITELIST"
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
  export WHITELIST="$repo_root/script/stage2_path_whitelist.txt"

  local tmpdir
  tmpdir="$(mktemp -d)" || exit 2
  trap 'rm -rf -- "${tmpdir:-}"' EXIT

  local stale_rows="$tmpdir/stale.tsv"       # literal<TAB>file:line
  : >"$stale_rows"
  local files_scanned=0 occurrences=0 whitelisted=0 stale=0

  # The two retired step roots, then the retired bundle name handled per line.
  local step_pattern='stage2/(02_cm_data|04_tahoe_data)[[:alnum:]_./+@-]*'

  # Production scans the tracked tree.  `STAGE2_PATH_AUDIT_FILES` is a
  # newline-separated file list used only by script/tests/ to exercise the
  # detection logic on fixtures without dirtying the tracked tree.
  local files_list="$tmpdir/files.z"
  if [ -n "${STAGE2_PATH_AUDIT_FILES:-}" ]; then
    tr '\n' '\0' <"$STAGE2_PATH_AUDIT_FILES" >"$files_list" || exit 2
  else
    git ls-files -z stage2 >"$files_list" || exit 2
  fi

  local f
  while IFS= read -r -d '' f; do
    case "$f" in
      stage2/_frozen/* | stage2/_archive/*) continue ;;
      *.py | *.sh | *.R) ;;
      *) continue ;;
    esac
    [ -f "$f" ] || continue
    files_scanned=$((files_scanned + 1))

    # 1. retired step path literals (one or more per line)
    local hit lineno literal
    while IFS= read -r hit; do
      [ -n "$hit" ] || continue
      lineno="${hit%%:*}"
      literal="${hit#*:}"
      occurrences=$((occurrences + 1))
      if is_whitelisted "$f" "$literal"; then
        whitelisted=$((whitelisted + 1))
      else
        stale=$((stale + 1))
        printf '%s\t%s:%s\n' "$literal" "$f" "$lineno" >>"$stale_rows"
      fi
    done < <(grep -onE -- "$step_pattern" "$f" 2>/dev/null)

    # 2. `datasets-1`, allowed only when bound to the frozen v1 root.
    local line
    while IFS= read -r line; do
      [ -n "$line" ] || continue
      lineno="${line%%:*}"
      line="${line#*:}"
      local stripped
      stripped="$(printf '%s' "$line" | sed -E 's#_frozen/cm-data-v1/datasets-1##g')"
      printf '%s' "$stripped" | grep -q 'datasets-1' || continue
      # `<FROZEN...> / "datasets-1..."` is a legitimate frozen-root binding.
      if printf '%s' "$stripped" | grep -qE 'FROZEN[A-Za-z0-9_]*[[:space:]]*/[[:space:]]*"[^"]*datasets-1'; then
        continue
      fi
      literal="$(printf '%s' "$line" | grep -oE '[[:alnum:]_./+@-]*datasets-1[[:alnum:]_./+@-]*' | head -n 1)"
      [ -n "$literal" ] || literal="datasets-1"
      occurrences=$((occurrences + 1))
      if is_whitelisted "$f" "$literal"; then
        whitelisted=$((whitelisted + 1))
      else
        stale=$((stale + 1))
        printf '%s\t%s:%s\n' "$literal" "$f" "$lineno" >>"$stale_rows"
      fi
    done < <(grep -nE -- 'datasets-1' "$f" 2>/dev/null)
  done <"$files_list"

  local distinct_stale
  distinct_stale="$(cut -f1 "$stale_rows" | LC_ALL=C sort -u | grep -c '' || true)"
  local whitelist_entries=0
  if [ -f "$WHITELIST" ]; then
    whitelist_entries="$(grep -cvE '^[[:space:]]*(#|$)' "$WHITELIST" || true)"
  fi

  echo "== check_stage2_paths: no retired step paths in tracked stage2 executable code =="
  echo "repo root: $repo_root"
  echo "tracked stage2 executable files scanned: $files_scanned"
  echo "stage2 path literals: $occurrences occurrences (retired 02/04 step roots + datasets-1)"
  echo "whitelisted: $whitelisted occurrences (entries in script/stage2_path_whitelist.txt: $whitelist_entries)"
  echo "stale: $stale occurrences, ${distinct_stale:-0} distinct"
  if [ "$stale" -gt 0 ]; then
    echo
    echo "STALE LITERALS (literal: location)"
    LC_ALL=C sort "$stale_rows" | awk -F'\t' '{
      if ($1 != prev) { printf "  %s\n", $1; prev = $1 }
      printf "      %s\n", $2
    }'
    echo
    echo "RESULT: FAIL — $stale stale occurrence(s) of ${distinct_stale:-0} distinct retired stage2 path literal(s)"
    echo "hint: repoint the literal at the current step, or add an entry to"
    echo "      script/stage2_path_whitelist.txt only for frozen recorded evidence."
    echo "exit: 1"
    return 1
  fi
  echo
  echo "RESULT: PASS — 0 stale stage2 path literals"
  echo "exit: 0"
  return 0
}

main "$@"
