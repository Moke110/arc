"""Read-only native source loading, including transitive sibling imports."""
from contextlib import contextmanager
import sys


@contextmanager
def no_source_bytecode():
    """Python must not add cache bytes anywhere in immutable source trees."""
    previous = sys.dont_write_bytecode
    sys.dont_write_bytecode = True
    try:
        yield
    finally:
        sys.dont_write_bytecode = previous
