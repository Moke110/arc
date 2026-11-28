"""Resource envelope inherited by the downstream CPST studies."""

from pathlib import Path

from matrix import GIB

MEMORY_CAP, SWAP_CAP = 12 * GIB, 8 * GIB


def scope_directory():
    group = Path("/proc/self/cgroup").read_text().strip().split(":", 2)[2]
    scope = Path("/sys/fs/cgroup") / Path(group).relative_to("/")
    for filename, limit in (("memory.max", MEMORY_CAP), ("memory.swap.max", SWAP_CAP)):
        if (scope / filename).read_text().strip() != str(limit):
            raise ValueError("CPST comparison matrix requires a 12 GiB RAM and 8 GiB swap scope")
    return scope
