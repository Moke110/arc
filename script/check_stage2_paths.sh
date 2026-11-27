#!/usr/bin/env bash
# check_stage2_paths.sh — assert stage2_CarPS-VAE/ executable code has no retired step paths.
#
# Scan tracked Stage 2 executable code for retired or archived paths.
# The layout guard allows Steps 01–08, util, routing docs, and the ignored local
# history directory. STAGE2_PATH_AUDIT_ROOT supplies a fixture layout.
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

  local old_alias='03[.][1-7]_(tenv|tahoe)_(hvg|denoise)-deg-[345]k'
  local step_pattern="(stage2_CarPS-VAE|stage2)/[[:alnum:]_./+@-]+|(06_tsg_drug_emb|07_tsg_dose|08_tsg_width_dose|09_tsg_transfer|03[.]1_hvg-deg-5k|$old_alias|log/316-runs)(/[[:alnum:]_./+@-]*)?"
  local active_roots=(util 01_data_raw 02_gene_axis 03_data_train 04_tsg_axis 05_tsg_loss 06_tsg_shape 07_cpsv_alpha 08_cpsv_chem-emb 09_cpsv_axis 10_cpsv_dose)
  local shape_rounds=(06_tsg_shape/06.1_round1 06_tsg_shape/06.2_round2 06_tsg_shape/06.3_round3 06_tsg_shape/06.4_round4)
  local stage2_root="${STAGE2_PATH_AUDIT_ROOT:-stage2_CarPS-VAE}"
  if [ -z "${STAGE2_PATH_AUDIT_ROOT:-}" ]; then
    if [ -e stage2 ] || [ -L stage2 ]; then
      echo "RESULT: FAIL — previous Stage 2 root remains: stage2" >&2
      return 1
    fi
    local archived_files
    archived_files="$(git ls-files -- stage2_CarPS-VAE/history)" || return 2
    if [ -n "$archived_files" ] || ! git check-ignore -q -- stage2_CarPS-VAE/history/.archive-probe; then
      echo "RESULT: FAIL — Stage 2 history must be ignored and untracked" >&2
      return 1
    fi
  fi
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
    case "$entry" in
      "$stage2_root/AGENTS.md" | "$stage2_root/README.md")
        [ -f "$entry" ] && [ ! -L "$entry" ] && continue ;;
      "$stage2_root/history")
        [ -d "$entry" ] && [ ! -L "$entry" ] && continue ;;
    esac
    echo "RESULT: FAIL — unexpected Stage 2 entry: $entry" >&2
    return 1
  done
  for entry in "$stage2_root"/06_tsg_shape/06.*; do
    [ -e "$entry" ] || [ -L "$entry" ] || continue
    for root in "${shape_rounds[@]}"; do
      if [ "$entry" = "$stage2_root/$root" ] && [ -d "$entry" ] && [ ! -L "$entry" ]; then
        continue 2
      fi
    done
    echo "RESULT: FAIL — unexpected Stage 2 shape round: $entry" >&2
    return 1
  done
  for root in "${active_roots[@]}" 05_tsg_loss/{05.1_loss_metrics,05.2_mse_l1,05.3_mse-l1_mse-x,05.4_mse-l1_pcc,05.5_mse-l1_cos} "${shape_rounds[@]}"; do
    for entry in "$stage2_root/$root/AGENTS.md" "$stage2_root/$root/README.md"; do
      if [ ! -f "$entry" ]; then
        echo "RESULT: FAIL — missing step routing: $entry" >&2
        return 1
      fi
    done
  done

  for entry in "$stage2_root/AGENTS.md" "$stage2_root/README.md"; do
    if [ ! -f "$entry" ]; then
      echo "RESULT: FAIL — missing Stage 2 routing: $entry" >&2
      return 1
    fi
  done
  for entry in "$stage2_root"/03_data_train/{03.1_hvg-deg-5k,output/tenv,output/tahoe-bulk,output/manifest.json} \
               "$stage2_root"/03_data_train/03.[1-7]_{tenv,tahoe}_{hvg,denoise}-deg-* \
               "$stage2_root"/04_tsg_axis/{output/formal,output/report,output/.report.publication.lock,log/316-runs}; do
    if [ -e "$entry" ] || [ -L "$entry" ]; then
      echo "RESULT: FAIL — superseded Stage 2 release in active tree: $entry" >&2
      return 1
    fi
  done
  local history_outside_archive
  history_outside_archive="$(find "${active_roots[@]/#/$stage2_root/}" \
    \( -name historical -o -name history \) -print -quit)" || return 2
  if [ -n "$history_outside_archive" ]; then
    echo "RESULT: FAIL — step-local Stage 2 history: $history_outside_archive" >&2
    return 1
  fi

  # Production scans the tracked tree.  `STAGE2_PATH_AUDIT_FILES` is a
  # newline-separated file list used only by script/tests/ to exercise the
  # detection logic on fixtures without dirtying the tracked tree.
  local files_list="$tmpdir/files.z"
  if [ -n "${STAGE2_PATH_AUDIT_FILES:-}" ]; then
    tr '\n' '\0' <"$STAGE2_PATH_AUDIT_FILES" >"$files_list" || exit 2
  else
    git ls-files -z stage2_CarPS-VAE >"$files_list" || exit 2
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
        case "$literal" in "stage2_CarPS-VAE/$root" | "stage2_CarPS-VAE/$root"/*) active_literal=1 ;; esac
      done
      case "$literal" in
        stage2_CarPS-VAE/06_tsg_shape/AGENTS.md | stage2_CarPS-VAE/06_tsg_shape/README.md) ;;
        stage2_CarPS-VAE/06_tsg_shape/?*)
          active_literal=0
          for root in "${shape_rounds[@]}"; do
            case "$literal" in "stage2_CarPS-VAE/$root" | "stage2_CarPS-VAE/$root"/*) active_literal=1 ;; esac
          done ;;
      esac
      case "$literal" in
        stage2_CarPS-VAE/05_tsg_loss/execution | stage2_CarPS-VAE/05_tsg_loss/execution/* | \
        stage2_CarPS-VAE/05_tsg_loss/output | stage2_CarPS-VAE/05_tsg_loss/output/* | \
        stage2_CarPS-VAE/05_tsg_loss/log | stage2_CarPS-VAE/05_tsg_loss/log/* | \
        stage2_CarPS-VAE/05_tsg_loss/figure | stage2_CarPS-VAE/05_tsg_loss/figure/*)
          active_literal=0 ;;
        */historical | */historical/* | */03.1_hvg-deg-5k | */03.1_hvg-deg-5k/*)
          active_literal=0 ;;
      esac
      if [[ "$literal" =~ (^|/)($old_alias|log/316-runs)(/|$) ]]; then
        active_literal=0
      fi
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
