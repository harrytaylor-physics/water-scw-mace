"""Atomic output writes (A86).

A producer that fails partway leaves a *partial* file with a *new* mtime, which
looks fresher than its inputs and passes every freshness check A85 added. That
happened during A85: diagnosing the S9 failure left `S9_case_evolution.csv` with
three of its four cases, and nothing could tell.

The fix is the usual one: write to a temporary file in the same directory, then
``os.replace`` it over the target, which is atomic on every platform this runs
on. A run that dies leaves the previous good file untouched and the temporary
behind, rather than a half-written output that reads as current.

Long, streaming producers are where this actually bites -- the four scan and
evolution CSVs, which write row by row over minutes. Single-call JSON dumps are
very nearly atomic already, but they go through the same helper so there is one
rule rather than two.
"""

from __future__ import annotations

import os
from contextlib import contextmanager
from pathlib import Path


def write_text(path, text: str, encoding: str = "utf-8") -> None:
    """Replace ``path`` with ``text``, atomically."""
    path = Path(path)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(text, encoding=encoding)
    os.replace(tmp, path)


@contextmanager
def open_write(path, mode: str = "w", **kwargs):
    """Open a temporary file and move it over ``path`` on clean exit only.

    On an exception the temporary is removed and ``path`` keeps whatever it had,
    so a producer that dies mid-write cannot leave a plausible-looking partial.
    """
    path = Path(path)
    tmp = path.with_name(path.name + ".tmp")
    fh = tmp.open(mode, **kwargs)
    try:
        yield fh
    except BaseException:
        fh.close()
        tmp.unlink(missing_ok=True)
        raise
    fh.close()
    os.replace(tmp, path)
