#!/usr/bin/env bash
# check_paths.sh — assert that every external-data path literal in tracked files exists.
#
# Why this exists
#   The provenance restructure (ADR 0030) renames and re-homes four external-source
#   trees (`pert_data/`, `pert_model/`, `disease_data/`, `drug_emb_model/`). Every
#   later ticket needs a mechanical red/green answer to "does the path literal I
#   just wrote or moved still point at something that exists on disk?".
#
# Usage
#   script/check_paths.sh                    scan and report
#   script/check_paths.sh --baseline [FILE]  scan, report, and write the counts
#                                            (default FILE:
#                                            .scratch/data-provenance-restructure/ops/path-baseline.txt)
#   script/check_paths.sh --help             print this text
#
# What it does
#   1. Takes the tracked file list from `git ls-files` (untracked-but-not-ignored
#      files are optional and are off by default, so scratch is never scanned).
#   2. Keeps text files only (binary files are skipped).
#   3. Extracts path literals beginning with `pert_data/`, `pert_model/`,
#      `disease_data/`, or `drug_emb_model/`, trims trailing prose punctuation, and
#      tests each one with `test -e` against the repository root.
#   4. Reports total / resolved / whitelisted / unresolved, and lists every
#      unresolved literal with a location.
#   5. `script/path_whitelist.txt` exempts frozen-evidence files that legitimately
#      keep historical strings: one glob-or-ERE per line, `#` comments. An entry is
#      matched against the repo-relative source file path and against the literal.
#
# Exit status
#   0  every extracted literal resolves (or is whitelisted)
#   1  at least one unresolved literal outside the whitelist
#   2  the script could not locate the repository root
#
# The repository root is found by walking up from this script until a directory
# containing .git is found — never by counting `../..` segments.
#
# No network; only writes the baseline file when --baseline is given.
# Idempotent: repeated runs give identical results.

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

# Strip prose punctuation that commonly follows a path in Markdown prose, e.g.
# "`pert_data/Tahoe-100M/source/`," or "(see pert_model/CPA/source/)." — never part of a path.
trim_literal() {
  printf '%s' "$1" | sed -E 's/[].,:;(){}<>|!?"'"'"'`]+$//'
}

# Whitelist entries are glob-or-ERE; match against the source file path or the
# literal. A match means "this occurrence is frozen evidence, keep the old string".
is_whitelisted() {
  local file="$1" literal="$2" line pat
  [ -f "$WHITELIST" ] || return 1
  while IFS= read -r line || [ -n "$line" ]; do
    line="${line%$'\r'}"
    if printf '%s' "$line" | grep -qE '^[[:space:]]*(#|$)'; then
      continue
    fi
    pat="$(printf '%s' "$line" | sed -E 's/^[[:space:]]+//; s/[[:space:]]+$//')"
    [ -n "$pat" ] || continue
    # shellcheck disable=SC2254  # glob match is intended here
    case "$file" in
      $pat) return 0 ;;
    esac
    if printf '%s\n%s\n' "$file" "$literal" | grep -qE -e "$pat" 2>/dev/null; then
      return 0
    fi
  done <"$WHITELIST"
  return 1
}

main() {
  local baseline_path=""
  case "${1:-}" in
    -h | --help)
      usage
      exit 0
      ;;
    --baseline)
      baseline_path="${2:-.scratch/data-provenance-restructure/ops/path-baseline.txt}"
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

  export WHITELIST="$repo_root/script/path_whitelist.txt"

  local tmpdir
  tmpdir="$(mktemp -d)" || exit 2
  trap 'rm -rf -- "${tmpdir:-}"' EXIT

  local files_list="$tmpdir/files.z"
  local occurrences="$tmpdir/occurrences.tsv"    # file<TAB>line<TAB>literal
  local unresolved_rows="$tmpdir/unresolved.tsv" # file:line<TAB>literal (failing)
  : >"$occurrences"
  : >"$unresolved_rows"

  # `(^|not-a-path-char)` prevents matching the tail of a longer identifier.
  local pattern='(^|[^[:alnum:]_./@-])(pert_data|pert_model|disease_data|drug_emb_model|pheno_data)/[[:alnum:]_./+@-]*'

  local files_scanned=0
  local f
  git ls-files -z | LC_ALL=C sort -z >"$files_list"
  while IFS= read -r -d '' f; do
    [ -f "$f" ] || continue
    # text files only: `grep -I` treats binary as non-matching; empty files skip.
    grep -Iq . -- "$f" 2>/dev/null || continue
    files_scanned=$((files_scanned + 1))
    grep -onE -- "$pattern" "$f" 2>/dev/null | while IFS= read -r hit; do
      local lineno="${hit%%:*}"
      local match="${hit#*:}"
      local literal
      case "$match" in
        pert_data/* | pert_model/* | disease_data/* | drug_emb_model/* | pheno_data/*) literal="$match" ;;
        *) literal="${match:1}" ;; # drop the one boundary character grep captured
      esac
      literal="$(trim_literal "$literal")"
      [ -n "$literal" ] || continue
      printf '%s\t%s\t%s\n' "$f" "$lineno" "$literal" >>"$occurrences"
    done
  done <"$files_list"

  local occ_total=0 occ_resolved=0 occ_whitelisted=0 occ_unresolved=0
  local src lineno literal
  while IFS=$'\t' read -r src lineno literal; do
    [ -n "$literal" ] || continue
    occ_total=$((occ_total + 1))
    if [ -e "$literal" ]; then
      occ_resolved=$((occ_resolved + 1))
    elif is_whitelisted "$src" "$literal"; then
      occ_whitelisted=$((occ_whitelisted + 1))
    else
      occ_unresolved=$((occ_unresolved + 1))
      printf '%s:%s\t%s\n' "$src" "$lineno" "$literal" >>"$unresolved_rows"
    fi
  done <"$occurrences"

  local distinct_total distinct_resolved distinct_whitelisted distinct_unresolved
  local all_literals="$tmpdir/all_literals.txt" failing_literals="$tmpdir/failing_literals.txt"
  cut -f3 "$occurrences" | LC_ALL=C sort -u >"$all_literals"
  cut -f2 "$unresolved_rows" | LC_ALL=C sort -u >"$failing_literals"
  distinct_total="$(grep -c '' "$all_literals" || true)"
  distinct_resolved=0
  distinct_whitelisted=0
  distinct_unresolved="$(grep -c '' "$failing_literals" || true)"
  local lit
  while IFS= read -r lit; do
    [ -n "$lit" ] || continue
    if [ -e "$lit" ]; then
      distinct_resolved=$((distinct_resolved + 1))
    elif ! grep -Fxq -- "$lit" "$failing_literals"; then
      distinct_whitelisted=$((distinct_whitelisted + 1))
    fi
  done <"$all_literals"

  local whitelist_entries=0
  if [ -f "$WHITELIST" ]; then
    whitelist_entries="$(grep -cvE '^[[:space:]]*(#|$)' "$WHITELIST" || true)"
  fi

  echo "== check_paths: every external-data path literal must resolve =="
  echo "repo root: $repo_root"
  echo "tracked text files scanned: $files_scanned"
  echo "path literals scanned: $occ_total occurrences, ${distinct_total:-0} distinct" \
    "(pert_data, pert_model, disease_data, drug_emb_model, pheno_data)"
  echo "resolved:    $occ_resolved occurrences, ${distinct_resolved:-0} distinct"
  echo "whitelisted: $occ_whitelisted occurrences, $distinct_whitelisted distinct" \
    "(entries in script/path_whitelist.txt: $whitelist_entries)"
  echo "unresolved:  $occ_unresolved occurrences, ${distinct_unresolved:-0} distinct"
  if [ "$occ_unresolved" -gt 0 ]; then
    echo
    echo "UNRESOLVED LITERALS (location: file:line)"
    LC_ALL=C sort "$unresolved_rows" | awk -F'\t' '{
      if ($2 != prev) { printf "  %s\n", $2; prev = $2 }
      printf "      %s\n", $1
    }'
  fi

  if [ -n "$baseline_path" ]; then
    case "$baseline_path" in
      /*) : ;;
      *) baseline_path="$repo_root/$baseline_path" ;;
    esac
    mkdir -p -- "$(dirname -- "$baseline_path")" || true
    local branch commit observed
    branch="$(git rev-parse --abbrev-ref HEAD 2>/dev/null || echo unknown)"
    commit="$(git rev-parse --short HEAD 2>/dev/null || echo unknown)"
    observed="$(date -Iseconds 2>/dev/null || date)"
    {
      echo "# check_paths.sh baseline record — pre-restructure external-data path literals"
      echo "# recorded_at: $observed"
      echo "# repo: $repo_root  branch: $branch  commit: $commit"
      echo "# categories scanned: pert_data/ pert_model/ disease_data/ drug_emb_model/"
      echo "tracked_text_files: $files_scanned"
      echo "literal_occurrences_total: $occ_total"
      echo "literal_occurrences_resolved: $occ_resolved"
      echo "literal_occurrences_whitelisted: $occ_whitelisted"
      echo "literal_occurrences_unresolved: $occ_unresolved"
      echo "literal_distinct_total: ${distinct_total:-0}"
      echo "literal_distinct_resolved: ${distinct_resolved:-0}"
      echo "literal_distinct_whitelisted: $distinct_whitelisted"
      echo "literal_distinct_unresolved: ${distinct_unresolved:-0}"
      echo "# unresolved literals (occurrence<tab>literal):"
      LC_ALL=C sort "$unresolved_rows"
    } >"$baseline_path"
    echo
    echo "baseline written: $baseline_path"
  fi

  # Ticket #17: the stage2 retired-step-path audit is part of this health check.
  # It is a separate script so the stage2 owners can run it alone; a red audit
  # makes this combined check red too.
  local stage2_status=0
  echo
  if [ -f "$repo_root/script/check_stage2_paths.sh" ]; then
    bash "$repo_root/script/check_stage2_paths.sh" || stage2_status=$?
  else
    echo "ERROR: script/check_stage2_paths.sh is missing (ticket #17)" >&2
    stage2_status=1
  fi

  if [ "$occ_unresolved" -eq 0 ] && [ "$stage2_status" -eq 0 ]; then
    echo
    echo "RESULT: PASS — every extracted path literal resolves"
    echo "stage2 retired-step-path audit: PASS"
    echo "exit: 0"
    return 0
  fi
  echo
  if [ "$occ_unresolved" -gt 0 ]; then
    echo "RESULT: FAIL — $occ_unresolved occurrence(s) of ${distinct_unresolved:-0} unresolved literal(s) outside the whitelist"
  fi
  if [ "$stage2_status" -ne 0 ]; then
    echo "RESULT: FAIL — the stage2 retired-step-path audit is red (see above)"
  fi
  echo "hint: add the owning files to script/path_whitelist.txt only for frozen evidence;"
  echo "      otherwise fix the literal or the path."
  echo "exit: 1"
  return 1
}

main "$@"
