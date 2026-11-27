"""Serialize Stage 2 release installation and preserve interrupted backups."""
from contextlib import contextmanager
import fcntl
from pathlib import Path
import shutil


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
