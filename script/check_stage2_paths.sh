#!/usr/bin/env bash
# check_stage2_paths.sh — assert stage2/ executable code has no retired step paths.
#
# Stage 2's active roots are 01_bases (shared foundation) and 02_models_tenv
# (first model-comparison step, ADR 0063).  Scan tracked active executable code
# for references to any other Stage 2 root.
# Historical payload literals require exact file/literal whitelist entries.
# Every invocation rejects any other immediate on-disk Stage 2 entry (directories
# and regular files alike; the routing document stage2/AGENTS.md is the only
# allowed root file) and checks each active root's owning AGENTS.md/README.md.
# The layout check
# runs in fixture mode too; `STAGE2_PATH_AUDIT_ROOT` overrides the scanned root
# for script/tests/.
# No network or dataset access. Exit 0 clean, 1 stale/layout failure, 2 usage error.

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

  # The retired step roots, then the retired bundle name handled per line.
  local step_pattern='stage2/[[:alnum:]_./+@-]+'

  # Accepted steps and the approved successor under construction (ADR 0074).
  # Allowing its paths does not accept its candidate or retire the old step.
  local active_roots=(01_bases 02_models_tenv 02_models_tenv-5k)

  # Layout guard: fail closed on retired step entries.  Only the active step
  # directories and the Stage 2 routing document stage2/AGENTS.md (present at the
  # review base a53c9ff) may sit at the Stage 2 root; any other regular file,
  # symlink or directory is a retired-step entry.  The earlier form replaced the
  # existence test with a bare `-d` test, which silently accepted any stray file
  # (#157 review S2).  `STAGE2_PATH_AUDIT_ROOT` is a script/tests/ hook;
  # production uses `stage2`.
  local stage2_root="${STAGE2_PATH_AUDIT_ROOT:-stage2}"
  local entry root active
  for entry in "$stage2_root"/* "$stage2_root"/.[!.]* "$stage2_root"/..?*; do
    [ -e "$entry" ] || [ -L "$entry" ] || continue
    active=0
    for root in "${active_roots[@]}"; do
      if [ "$entry" = "$stage2_root/$root" ] && [ -d "$entry" ] && [ ! -L "$entry" ]; then
        active=1
        break
      fi
    done
    [ "$active" -eq 1 ] && continue
    if [ "$entry" = "$stage2_root/AGENTS.md" ] && [ -f "$entry" ] && [ ! -L "$entry" ]; then
      continue
    fi
    echo "RESULT: FAIL — unexpected Stage 2 entry: $entry" >&2
    return 1
  done
  for root in "${active_roots[@]}"; do
    for entry in "$stage2_root/$root/AGENTS.md" "$stage2_root/$root/README.md"; do
      if [ ! -f "$entry" ]; then
        echo "RESULT: FAIL — missing step routing: $entry" >&2
        return 1
      fi
    done
  done

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
      local active_literal=0
      for root in "${active_roots[@]}"; do
        case "$literal" in "stage2/$root" | "stage2/$root"/*) active_literal=1 ;; esac
      done
      [ "$active_literal" -eq 1 ] && continue
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
      stripped="$line"
      printf '%s' "$stripped" | grep -q 'datasets-1' || continue
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
  echo "stage2 path literals: $occurrences occurrences (retired roots + datasets-1)"
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
