"""Serialize Stage 2 release installation and preserve interrupted backups."""
from contextlib import contextmanager
import fcntl
import os
from pathlib import Path
import shutil
from tempfile import mkdtemp


def _previous(output):
    return output.with_name(f".{output.name}.previous")


def _recover(output):
    previous = _previous(output)
    if output.is_symlink() or previous.is_symlink():
        raise ValueError(f"release and backup must be actual directories: {output}")
    if not output.exists() and previous.exists():
        previous.rename(output)


@contextmanager
def publication_lock(destination):
    """Serialize all operations against a release, independent of its staging path."""
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    output = destination.parent.resolve() / destination.name
    with output.with_name(f".{output.name}.publication.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        _recover(output)
        yield output


def replace_release(candidate, output, replace, *, verify_existing=None):
    """Install a candidate under its publication lock and recover on a failed swap."""
    previous = _previous(output)
    if previous.exists():
        if verify_existing is None:
            raise ValueError(f"cannot discard an unverified release backup: {previous}")
        verify_existing(output)
        shutil.rmtree(previous)
    try:
        return replace(candidate, output)
    except BaseException:
        _recover(output)
        raise


def publish_candidate(release, output, build, verify, *, staging=None, verify_existing=None):
    """Build and verify a candidate, then install it under an already-held lock."""
    candidate = release.candidate(output, staging)
    try:
        build(candidate)
        verify(candidate)
        return replace_release(candidate, output, release.replace,
                               verify_existing=verify_existing or verify)
    finally:
        shutil.rmtree(candidate, ignore_errors=True)


def publish_files(updates, *, staging):
    """Publish a preflighted file set under its destination lock, rolling back on error."""
    if not updates:
        return
    staging.mkdir(parents=True, exist_ok=True)
    work = Path(mkdtemp(prefix="cpst-publication-", dir=staging))
    candidates = []
    published = []
    preserve_backup = False
    try:
        for index, (target, contents) in enumerate(updates.items()):
            # Keep the old inode. A copied payload can share an inode with its sealed source.
            backup = work / f"{index}.previous" if target.exists() else None
            if backup is not None:
                os.link(target, backup)
            elif target.parent.stat().st_dev != work.stat().st_dev:
                raise ValueError(f"publication staging and destination use different filesystems: {target}")
            candidate = work / str(index)
            candidate.write_bytes(contents)
            candidate.chmod(backup.stat().st_mode if backup is not None else 0o644)
            if candidate.read_bytes() != contents:
                raise OSError(f"staged file differs from verified publication: {target}")
            candidates.append((candidate, target, backup))
        for candidate, target, backup in candidates:
            published.append((candidate, target, backup))
            candidate.replace(target)
    except BaseException as error:
        failures = []
        for candidate, target, backup in reversed(published):
            try:
                if candidate.exists():
                    continue
                if backup is None:
                    target.unlink()
                else:
                    backup.replace(target)
            except OSError as failure:
                failure.add_note(f"Restore {target} from {backup}")
                failures.append(failure)
        if failures:
            preserve_backup = True
            raise BaseExceptionGroup(f"publication rollback incomplete; backups retained at {work}",
                                     [error, *failures])
        raise
    finally:
        if not preserve_backup:
            shutil.rmtree(work, ignore_errors=True)
