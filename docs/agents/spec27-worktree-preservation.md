# Spec 27 local worktree preservation

Ticket #408 emptied the repository's root `.worktrees/` on 2026-11-28.
The `dev` checkout moved to `tmp/spec27-dev/` at `b6f588bc`. The remaining
entries moved intact into the ignored local archive
`archived/spec27-retained-worktrees/`. Git repaired 17 registrations that
still pointed to the former `/home/chang/arc` checkout. The serial cleanup
worktree also moved with its Git registration.

The relocation audit compared 44,117 paths and 19 worktree statuses before
and after the move. It checked filesystem identity, file sizes and
modification times, symlink targets, the Git index status, and HEAD.
`stage2-deg1000-tenv-pooling` and `stage2-reorder-plan` retained 37 and two
status entries respectively. Ignored payloads moved alongside the tracked
and untracked files. All existing branch tips were preserved. The detached
historical tip `986f8ea` is an ancestor of `dev`.

From the root checkout, run the saved local audit to verify the preserved
archive and its Git registrations. The original full comparison also
includes the pre-integration `dev` snapshot, which will differ once new work
has been merged.

```bash
(
  cd archived/spec27-retained-worktrees/.relocation-audit
  sha256sum -c SHA256SUMS
)
python archived/spec27-retained-worktrees/.relocation-audit/snapshot.py archive
(
  cd archived/spec27-retained-worktrees/spec23-retained-provenance
  sha256sum -c SHA256SUMS
)
```

The archive is local evidence. Keep it outside disposable checkout scratch
when running cleanup. The archive contains the complete saved inventory and
the pre-integration verification output under `.relocation-audit/`.
