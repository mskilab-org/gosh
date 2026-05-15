"""
gosh cache — build a unified work-dir cache for protect and purge.

Scans the results directory for symlinks into the work directory (BFS,
transitively following staged-input symlinks), then does a shallow two-level
enumeration of every hash dir on disk.  The difference gives two disjoint
sets:

  protect   — hash dirs reachable from results (completed processes)
  purge     — hash dirs on disk but *not* reachable (failed / orphaned)

The cache is written as JSON to <work_dir>/.gosh_cache.json (overrideable
with --output).  Both `gosh protect` and `gosh purge` accept --cache-file to
consume this file instead of re-running the scan themselves, making repeated
invocations essentially free.
"""

import json
import os
import sys
from datetime import datetime
from pathlib import Path

import click

CACHE_FILENAME = ".gosh_cache.json"


def _load_cache(cache_path: Path) -> dict | None:
    """Load and validate a previously written cache file.

    Returns the parsed dict on success, or None if the file is missing,
    malformed, or has an unexpected schema.  Emits a yellow warning so
    callers can fall back gracefully.
    """
    if not cache_path.exists():
        return None
    try:
        with open(cache_path) as f:
            data = json.load(f)
    except (json.JSONDecodeError, OSError) as exc:
        click.secho(f"Warning: could not read cache {cache_path}: {exc}", fg="yellow")
        return None
    if not isinstance(data, dict) or "protect" not in data or "purge" not in data:
        click.secho(
            f"Warning: cache file {cache_path} has unexpected schema; ignoring.",
            fg="yellow",
        )
        return None
    return data


def _save_cache(cache_path: Path, protect: set[Path], purge: set[Path]) -> None:
    """Atomically write the cache to *cache_path*.

    Paths are serialised as absolute POSIX strings and sorted for stable diffs.
    A write-then-rename ensures readers never see a partial file.
    """
    payload = {
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "protect": sorted(str(p) for p in protect),
        "purge": sorted(str(p) for p in purge),
    }
    tmp = cache_path.with_name(cache_path.name + ".tmp")
    with open(tmp, "w") as f:
        json.dump(payload, f, indent=2)
    os.replace(tmp, cache_path)


@click.command(name="cache")
@click.option(
    "-p",
    "--pipeline-output-dir",
    "results_dir",
    default="./results/",
    show_default=True,
    type=click.Path(file_okay=False, dir_okay=True),
    help="Results directory (Nextflow publishDir root) to scan.",
)
@click.option(
    "-w",
    "--work-dir",
    default="./work/",
    show_default=True,
    type=click.Path(file_okay=False, dir_okay=True),
    help="Nextflow work directory to enumerate.",
)
@click.option(
    "-o",
    "--output",
    "output_path",
    default=None,
    type=click.Path(dir_okay=False),
    help=(
        "Where to write the cache JSON.  "
        f"Defaults to <work-dir>/{CACHE_FILENAME}."
    ),
)
@click.option(
    "-c",
    "--cache-file",
    "cache_file",
    default=None,
    type=click.Path(exists=True, dir_okay=False),
    help=(
        "Verify (or display) an existing cache file instead of building a new one.  "
        "Requires --verify or --show-purge / --show-protect; skips writing."
    ),
)
@click.option(
    "-n",
    "--dry-run",
    is_flag=True,
    help="Print the cache contents without writing to disk.",
)
@click.option(
    "--show-purge",
    is_flag=True,
    help="Also print the list of purgeable work hash dirs.",
)
@click.option(
    "--show-protect",
    is_flag=True,
    help="Also print the list of work hash dirs to protect.",
)
@click.option(
    "-V",
    "--verify",
    "do_verify",
    is_flag=True,
    help=(
        "After building (or loading) the cache, re-scan results independently "
        "and assert that no symlink there points into the purge set.  "
        "Exits non-zero if any violations are found."
    ),
)
def cache_cli(results_dir, work_dir, output_path, cache_file, dry_run,
              show_purge, show_protect, do_verify):
    """
    Build a unified protect/purge cache for the work directory.

    The cache records which work hash directories are referenced from results
    (and should be protected) and which are not (and may be purged).  Both
    `gosh protect` and `gosh purge` can consume this file with --cache-file to
    avoid rescanning the filesystem.

    \b
    Scan strategy (fast even with thousands of failures):
      1. BFS from results symlinks → reachable (completed) hash dirs.
         Only the top level of each hash dir is scanned — no os.walk.
      2. Shallow 2-level os.scandir of work/ → every hash dir on disk.
      3. purge = all_on_disk - reachable  (pure set subtraction, no I/O).

    \b
    Verification (--verify):
      Re-walks results with a plain os.walk — entirely separate from the BFS
      used to build the cache — and checks that no symlink there resolves into
      the purge set.  Exits 1 if violations are found.
    """
    from ..utils.find_work_deps import build_work_dir_cache, verify_cache

    work_path = Path(work_dir).resolve()
    if not work_path.is_dir():
        click.secho(f"Error: work dir not found: {work_path}", fg="red")
        sys.exit(1)

    results_path = Path(results_dir).resolve()
    if not results_path.is_dir():
        click.secho(f"Error: results dir not found: {results_path}", fg="red")
        sys.exit(1)

    # ------------------------------------------------------------------ #
    # Resolve protect / purge sets — from a pre-built file or live scan.  #
    # ------------------------------------------------------------------ #
    if cache_file:
        loaded = _load_cache(Path(cache_file))
        if loaded is None:
            click.secho(
                f"Error: cache file {cache_file} could not be loaded.", fg="red"
            )
            sys.exit(1)
        protect_dirs: set[Path] = {Path(p) for p in loaded["protect"]}
        purge_dirs: set[Path] = {Path(p) for p in loaded["purge"]}
        click.secho(
            f"Loaded cache from {cache_file}: "
            f"{len(protect_dirs)} protect, {len(purge_dirs)} purge.",
            fg="blue",
        )
    else:
        cache_path = (
            Path(output_path).resolve() if output_path else work_path / CACHE_FILENAME
        )
        click.secho(f"Scanning {results_path} -> {work_path}...", fg="blue")
        built = build_work_dir_cache(results_path, work_path)
        protect_dirs = built["protect"]
        purge_dirs = built["purge"]
        click.secho(
            f"Found {len(protect_dirs)} work dirs to protect, "
            f"{len(purge_dirs)} to purge.",
            fg="blue",
        )

    # ------------------------------------------------------------------ #
    # Optional display.                                                    #
    # ------------------------------------------------------------------ #
    if show_protect:
        click.secho("\n# protect:", fg="cyan")
        for d in sorted(protect_dirs):
            click.echo(f"  {d.relative_to(work_path)}")

    if show_purge:
        click.secho("\n# purge:", fg="cyan")
        for d in sorted(purge_dirs):
            click.echo(f"  {d.relative_to(work_path)}")

    # ------------------------------------------------------------------ #
    # Write (skipped when loading from an existing file or dry-running).  #
    # ------------------------------------------------------------------ #
    if not cache_file:
        if dry_run:
            click.secho(
                f"\nDry run — cache not written to {cache_path}.", fg="yellow"
            )
        else:
            try:
                _save_cache(cache_path, protect_dirs, purge_dirs)
            except OSError as exc:
                click.secho(
                    f"Error: could not write cache to {cache_path}: {exc}", fg="red"
                )
                sys.exit(1)
            click.secho(f"\nCache written to {cache_path}.", fg="green")

    # ------------------------------------------------------------------ #
    # Orthogonal verification — independent os.walk, separate code path.  #
    # ------------------------------------------------------------------ #
    if not do_verify:
        return

    click.secho(
        f"\nVerifying: re-scanning {results_path} for purge-set collisions...",
        fg="blue",
    )
    violations = verify_cache(results_path, purge_dirs, work_path)

    if not violations:
        click.secho(
            f"OK — no results symlinks point into the purge set "
            f"({len(purge_dirs)} purge dirs checked).",
            fg="green",
        )
        return

    # Violations found — report every symlink that triggered each one.
    click.secho(
        f"\nERROR: {len(violations)} purge-set collision(s) detected!\n",
        fg="red",
        bold=True,
    )
    for wdir in sorted(violations):
        click.secho(f"  {wdir.relative_to(work_path)}", fg="red")
        for symlink in sorted(violations[wdir]):
            try:
                rel_sym = symlink.relative_to(results_path)
            except ValueError:
                rel_sym = symlink
            click.echo(f"    <- {rel_sym}")

    click.secho(
        "\nThe purge set contains work dirs still referenced from results.  "
        "Do not purge until this is resolved.",
        fg="red",
    )
    sys.exit(1)
