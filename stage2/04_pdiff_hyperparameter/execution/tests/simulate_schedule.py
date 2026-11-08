"""Test-only scheduler entry point with fixture telemetry and a private lock."""
from contextlib import contextmanager
import fcntl
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'script'))
import schedule


def main():
    fixture = Path(sys.argv.pop(1))
    snapshot = schedule.snapshot

    @contextmanager
    def isolated_lock():
        with fixture.with_suffix('.lock').open('a') as lock:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            yield

    schedule.snapshot = lambda record, pids: snapshot(record, fixture, pids)
    schedule.workstation_lock = isolated_lock
    return schedule.main()


if __name__ == '__main__':
    sys.exit(main())
