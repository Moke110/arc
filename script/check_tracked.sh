#!/usr/bin/env bash
# check_tracked.sh — assert that every executable file in the worktree is tracked by git.
#
# Why this exists
#   The repository's .gitignore once ignored whole directory trees, so a code file
#   could be silently swallowed by a rule meant for data. The historical example is
#   the vendored KPGT library: its `src/data/` tree was untracked because
#   the global `data/` rule matched any directory named `data`, not just `<repo>/data/`.
#   Ticket 05 narrowed that rule to `/data/` and the KPGT code (now tracked at
#   `drug_emb_model/KPGT/source/src/data/`) is in git; the detection below remains as
#   a regression guard. A file that exists on disk but not in git cannot be
#   reproduced from a clone.
#
# Usage
#   script/check_tracked.sh             run the check
#   script/check_tracked.sh --help      print this text
#
# What it does
#   1. Finds every *.py / *.sh / *.R in the worktree, pruning:
#        - .git/ and __pycache__/
#        - virtualenv trees (any directory named _envs, .venv, or venv)
#        - external item trees (any directory named source/ or local/) — the
#          download and the declared local additions are not repository content
#          (ADR 0030)
#        - workspace/archive trees that are not repository content:
#          .scratch/, .worktrees/ (separate Git checkouts), tmp/, reference/, wet_data/, .uploads/
#   2. Compares that list with `git ls-files`.
#   3. Groups every untracked file by the .gitignore rule that hides it, flags any
#      file hidden by an over-broad `data/` rule explicitly, and reports trees that
#      are untracked *by documented design* (EXPECTED_IGNORED below) separately, so
#      nothing is silent.
#   4. Prints a PASS/FAIL summary and the exit code.
#
# Exit status
#   0  every scanned executable file is tracked, or lives in a documented
#      expected-ignored tree
#   1  at least one unexpected untracked executable file exists (a hole)
#   2  the script could not locate the repository root
#
# The repository root is found by walking up from this script until a directory
# containing .git is found — never by counting `../..` segments.
#
# No network, no writes, idempotent.

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

# Trees whose code is untracked on purpose. Each case must name the documented
# reason; these files are counted and reported (never silent) but do not fail the
# check. Remove entries as the restructure retires the tree.
is_expected_ignored() {
  case "$1" in
    # Frozen pre-refactor MIRI snapshot, deliberately kept out of git
    # (.gitignore: "Frozen pre-refactor MIRI pipeline snapshots ...").
    stage1/xpert-fda-miri-v0/*)
      return 0 ;;
  esac
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
  local disk_list="$tmpdir/disk.z" tracked_list="$tmpdir/tracked.z"
  local holes_file="$tmpdir/holes.z"
  local fail_rows="$tmpdir/fail.tsv" exp_rows="$tmpdir/expected.txt"
  : >"$fail_rows"
  : >"$exp_rows"

  # 1. Code files on disk (find prunes, NUL-safe for odd path names).
  find . \
    \( -path './.git' -o -path './.worktrees' -o -path './.scratch' -o -path './tmp' -o -path './reference' \
       -o -path './wet_data' -o -path './.uploads' \
       -o -name '_envs' -o -name '.venv' -o -name 'venv' -o -name 'source' \
       -o -name 'local' -o -name '__pycache__' \) -prune \
    -o -type f \( -name '*.py' -o -name '*.sh' -o -name '*.R' \) -print0 \
    | sed -z 's|^\./||' | LC_ALL=C sort -z >"$disk_list"

  # 2. Code files tracked by git (NUL-safe; -z avoids the C-quoting of non-ASCII
  #    path names that would otherwise break the comparison).
  git ls-files -z -- '*.py' '*.sh' '*.R' | LC_ALL=C sort -z >"$tracked_list"

  comm -z -23 "$disk_list" "$tracked_list" >"$holes_file" 2>/dev/null || true

  local on_disk tracked
  on_disk="$(grep -zc '' "$disk_list" 2>/dev/null || true)"
  tracked="$(grep -zc '' "$tracked_list" 2>/dev/null || true)"

  # 3. Classify every untracked file.
  local f rule exp_count=0 fail_count=0
  while IFS= read -r -d '' f; do
    if is_expected_ignored "$f"; then
      exp_count=$((exp_count + 1))
      printf '%s\n' "$f" >>"$exp_rows"
    else
      rule="$(git check-ignore -v -- "$f" 2>/dev/null | cut -f1)"
      [ -n "$rule" ] || rule="(not ignored by .gitignore)"
      fail_count=$((fail_count + 1))
      printf '%s\t%s\n' "$rule" "$f" >>"$fail_rows"
    fi
  done <"$holes_file"

  echo "== check_tracked: every executable file must be tracked by git =="
  echo "repo root: $repo_root"
  echo "scan: $(printf '%s' "$on_disk") code files on disk (excluding virtualenvs, source/, __pycache__/, .git/, workspace archives); $(printf '%s' "$tracked") tracked by git"
  echo "untracked: $((exp_count + fail_count))  (documented expected-ignored: $exp_count, holes: $fail_count)"
  echo

  # 4a. Regression guard: an over-broad `data/` rule swallowing code again.
  local data_rule_rows="$tmpdir/data_rule.tsv"
  awk -F'\t' '$1 ~ /:data\/$/ { print $2 }' "$fail_rows" | LC_ALL=C sort >"$data_rule_rows"
  local data_rule_count
  data_rule_count="$(grep -c '' "$data_rule_rows" 2>/dev/null || true)"
  if [ "${data_rule_count:-0}" -gt 0 ]; then
    echo "HOLE — an over-broad \`data/\` ignore rule swallows code ($data_rule_count files)"
    echo "  rule: a `data/` pattern matches ANY directory named data, not just <repo>/data/"
    awk -F/ '{
      for (i = 1; i < NF; i++) if ($i == "data") {
        d = "";
        for (j = 1; j < i; j++) d = d (j > 1 ? "/" : "") $j;
        print d "/data/";
        next
      }
    }' "$data_rule_rows" | LC_ALL=C sort | uniq -c | awk '{ printf "  %6d  %s\n", $1, $2 }'
    echo "  fix: track those files, or narrow the rule to \`/data/\` (root-only)"
    echo
  fi

  # 4b. Other unexpected holes, grouped by suppressing rule.
  local other_rows="$tmpdir/other.tsv"
  awk -F'\t' '$1 !~ /:data\/$/ { print }' "$fail_rows" | LC_ALL=C sort >"$other_rows"
  local other_count
  other_count="$(grep -c '' "$other_rows" 2>/dev/null || true)"
  if [ "${other_count:-0}" -gt 0 ]; then
    echo "UNEXPECTED UNTRACKED CODE ($other_count files, grouped by suppressing rule)"
    awk -F'\t' '{
      if ($1 != prev) { if (prev != "") print ""; printf "  [%s]\n", $1; prev=$1 }
      printf "    %s\n", $2
    }' "$other_rows"
    echo
  fi

  # 4c. Documented expected-ignored trees (reported, not failed).
  if [ "$exp_count" -gt 0 ]; then
    echo "EXPECTED-IGNORED TREES (documented; reported, not failed)"
    sed -E 's|^(pert_model/[^/]+)/.*|\1|; s|^(stage1/[^/]+)/.*|\1|' "$exp_rows" \
      | LC_ALL=C sort | uniq -c | awk '{ printf "  %6d  %s\n", $1, $2 }'
    echo "  rationale: frozen v0 and ADR 0042 old formal recovery copies stay out of git. External items"
    echo "  (source/ and local/) are pruned by the scan itself, not listed here. Remove an"
    echo "  entry from EXPECTED_IGNORED when its tree retires."
    echo
  fi

  if [ "$fail_count" -eq 0 ]; then
    echo "RESULT: PASS — all scanned executable files are tracked (or documented expected-ignored)"
    echo "exit: 0"
    return 0
  fi
  echo "RESULT: FAIL — $fail_count untracked executable file(s) are not accounted for"
  echo "exit: 1"
  return 1
}

main "$@"
