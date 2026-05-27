"""
Nextflow execution-trace parser.

Reads all ``execution_trace_*.txt`` files from a Nextflow ``pipeline_info``
directory and returns a list of task records as plain dicts.  For every task
the ``work_dir`` key is resolved to the actual directory under
``<pipeline_info>/../../work`` by expanding the short hash (e.g. ``02/5852bd``)
to the matching full subdirectory on disk.

Usage
-----
    from gosh_cli.cli.parse_trace import parse_traces

    records = parse_traces("/path/to/pipeline_info")
    for r in records:
        print(r["name"], r["work_dir"])
"""

from __future__ import annotations

import glob
import os
from pathlib import Path
from typing import Optional


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _resolve_work_dir(hash_field: str, work_root: Path) -> Optional[str]:
    """Expand a Nextflow short hash (``XX/YYYYYY``) to the full work directory.

    Nextflow stores task outputs in ``work/<prefix>/<full_hash>/``.
    The trace file only records the first few characters of the full hash as
    the suffix part (e.g. ``02/5852bd``).  We glob for
    ``work/02/5852bd*/`` and return the first match.

    Returns ``None`` when the directory cannot be found on disk (e.g. the
    work directory has been cleaned).
    """
    if not hash_field or hash_field == "-":
        return None

    parts = hash_field.split("/", 1)
    if len(parts) != 2:
        return None

    prefix, suffix = parts
    pattern = str(work_root / prefix / f"{suffix}*")
    matches = glob.glob(pattern)

    if not matches:
        return None

    # Prefer an exact directory match; fall back to the first glob hit.
    dirs = [m for m in matches if os.path.isdir(m)]
    return dirs[0] if dirs else None


def _parse_value(value: str) -> str:
    """Strip surrounding whitespace; keep the raw string representation.

    Nextflow uses ``-`` for missing / not-applicable values.
    """
    return value.strip()


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def parse_trace_file(trace_path: str | Path, work_root: Path) -> list[dict]:
    """Parse a single Nextflow execution-trace TSV file.

    Parameters
    ----------
    trace_path:
        Absolute (or resolvable) path to an ``execution_trace_*.txt`` file.
    work_root:
        Absolute path to the ``work/`` directory that corresponds to the run.

    Returns
    -------
    list[dict]
        One dict per task row.  All original trace columns are preserved as
        string values.  An additional ``work_dir`` key is added with the
        resolved absolute path (or ``None`` when the directory is not found).
        A ``trace_file`` key records which source file the row came from.
    """
    trace_path = Path(trace_path)
    records: list[dict] = []

    with trace_path.open(encoding="utf-8") as fh:
        header_line = fh.readline()
        if not header_line:
            return records

        columns = [c.strip() for c in header_line.rstrip("\n").split("\t")]
        if "hash" not in columns:
            raise ValueError(
                f"Expected a 'hash' column in {trace_path}; "
                f"found columns: {columns}"
            )

        for lineno, line in enumerate(fh, start=2):
            line = line.rstrip("\n")
            if not line:
                continue

            fields = line.split("\t")
            # Guard against ragged rows (shouldn't happen, but be safe).
            if len(fields) < len(columns):
                fields.extend([""] * (len(columns) - len(fields)))

            row: dict = {col: _parse_value(fields[i]) for i, col in enumerate(columns)}
            row["work_dir"] = _resolve_work_dir(row.get("hash", ""), work_root)
            row["trace_file"] = trace_path.name
            records.append(row)

    return records


def parse_traces(
    pipeline_info_dir: str | Path,
    *,
    work_dir: str | Path | None = None,
    pattern: str = "execution_trace_*.txt",
) -> list[dict]:
    """Parse all execution-trace files found in *pipeline_info_dir*.

    Parameters
    ----------
    pipeline_info_dir:
        Path to the Nextflow ``pipeline_info`` output directory.
    work_dir:
        Override the default work directory (``../../work`` relative to
        *pipeline_info_dir*).  Pass an absolute path to be explicit.
    pattern:
        Glob pattern used to discover trace files inside *pipeline_info_dir*.
        Defaults to ``execution_trace_*.txt``.

    Returns
    -------
    list[dict]
        Combined list of task records from all matching trace files, sorted by
        ``submit`` timestamp then ``task_id``.  Duplicate task entries across
        multiple trace files are **not** deduplicated — the caller can filter
        by ``trace_file`` or ``hash`` as needed.
    """
    pipeline_info_dir = Path(pipeline_info_dir).resolve()

    if work_dir is None:
        work_root = (pipeline_info_dir / "../../work").resolve()
    else:
        work_root = Path(work_dir).resolve()

    trace_files = sorted(pipeline_info_dir.glob(pattern))
    if not trace_files:
        raise FileNotFoundError(
            f"No files matching '{pattern}' found in {pipeline_info_dir}"
        )

    all_records: list[dict] = []
    for tf in trace_files:
        all_records.extend(parse_trace_file(tf, work_root))

    # Sort by submit time (string sort works for the ISO-like NF format),
    # then by task_id as a tie-breaker.
    all_records.sort(
        key=lambda r: (r.get("submit", ""), r.get("task_id", ""))
    )

    return all_records


# ---------------------------------------------------------------------------
# CLI convenience (``python -m gosh_cli.cli.parse_trace <pipeline_info_dir>``)
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import sys
    import json as _json

    if len(sys.argv) < 2:
        print(
            "Usage: python parse_trace.py <pipeline_info_dir> [work_dir]",
            file=sys.stderr,
        )
        sys.exit(1)

    _pipeline_info = sys.argv[1]
    _work_dir = sys.argv[2] if len(sys.argv) > 2 else None

    _records = parse_traces(_pipeline_info, work_dir=_work_dir)
    _json.dump(_records, sys.stdout, indent=2)
    print()  # trailing newline
