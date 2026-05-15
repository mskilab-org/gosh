"""
gosh purge — remove Nextflow work directories not referenced from a results directory.

Scans every 2-level hash dir under the work directory, subtracts the set that
is referenced (directly or transitively) by symlinks in the results directory,
and deletes the remainder.

Pass --cache-file (or -c) to point at a cache produced by `gosh cache` so the
filesystem scan is skipped entirely and the command runs in near-zero time.

Safe-guards:
  • --dry-run  shows exactly what would be removed without touching the filesystem.
  • Confirmation prompt before any deletion (bypass with --yes / -y).
  • Any hash dir recorded in the protect state file is skipped with a warning,
    even if it is no longer reachable from the results directory.
"""

import os
import shutil
import stat
import sys
from pathlib import Path

import click

from .protect import STATE_FILENAME, _load_state
from .cache import CACHE_FILENAME, _load_cache


def _dir_size(path: Path) -> int:
    """Return the sum of st_size for every non-symlink file under path."""
    total = 0
    try:
        for dirpath, _, filenames in os.walk(path, followlinks=False):
            for name in filenames:
                try:
                    s = os.stat(os.path.join(dirpath, name), follow_symlinks=False)
                    if not stat.S_ISLNK(s.st_mode):
                        total += s.st_size
                except OSError:
                    pass
    except OSError:
        pass
    return total


def _human_size(n: int) -> str:
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if n < 1024.0:
            return f"{n:.1f} {unit}"
        n /= 1024.0
    return f"{n:.1f} PB"


def _rmtree_writable(path: Path) -> None:
    """Remove a directory tree, adding owner-write to blocked directories first.

    Works bottom-up so that by the time shutil.rmtree tries to rmdir a
    directory, its children are already gone and its own write bit is set."""
    for dirpath, dirnames, _ in os.walk(path, topdown=False, followlinks=False):
        for name in dirnames:
            d = Path(dirpath) / name
            if d.is_symlink():
                continue
            try:
                m = stat.S_IMODE(os.stat(d, follow_symlinks=False).st_mode)
                if not (m & stat.S_IWUSR):
                    os.chmod(d, m | stat.S_IWUSR)
            except OSError:
                pass
    shutil.rmtree(path)


@click.command(name="purge")
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
    help="Nextflow work directory to purge unreferenced hash dirs from.",
)
@click.option(
    "-c",
    "--cache-file",
    "cache_file",
    default=None,
    type=click.Path(exists=True, dir_okay=False),
    help=(
        "Path to a cache file produced by `gosh cache`.  "
        "When supplied the filesystem scan is skipped entirely."
    ),
)
@click.option(
    "-n",
    "--dry-run",
    is_flag=True,
    help="Show what would be removed without deleting anything.",
)
@click.option(
    "-y",
    "--yes",
    is_flag=True,
    help="Skip the confirmation prompt.",
)
@click.option(
    "-s",
    "--show-size",
    is_flag=True,
    help="Calculate and display disk usage for each purgeable directory (slower).",
)
def purge_cli(results_dir, work_dir, cache_file, dry_run, yes, show_size):
    """
    Remove work directories not referenced from the results directory.

    Referenced directories (those reachable via symlinks in the results
    directory, transitively) are kept.  Any directory recorded in the protect
    state file is also kept, even if it is no longer referenced — run
    `gosh protect --unprotect` first if you want those removed too.

    Run `gosh cache` first and pass the result via --cache-file to avoid
    rescanning the filesystem on repeated invocations.
    """
    from ..utils.find_work_deps import build_work_dir_cache

    work_path = Path(work_dir).resolve()
    if not work_path.is_dir():
        click.secho(f"Error: work dir not found: {work_path}", fg="red")
        sys.exit(1)

    # Load the protect state so we never accidentally purge a protected dir.
    state_path = work_path / STATE_FILENAME
    protected = set(_load_state(state_path).keys())

    # --- Resolve the purge/protect split ---------------------------------
    # Prefer a pre-built cache file when one is supplied; fall back to a live
    # scan only when needed.
    if cache_file:
        loaded = _load_cache(Path(cache_file))
        if loaded is None:
            click.secho(
                f"Error: cache file {cache_file} could not be loaded.",
                fg="red",
            )
            sys.exit(1)
        purge_set = {Path(p) for p in loaded["purge"]}
        n_protect = len(loaded["protect"])
        click.secho(
            f"Loaded cache from {cache_file}: "
            f"{n_protect} protect, {len(purge_set)} purge.",
            fg="blue",
        )
    else:
        results_path = Path(results_dir).resolve()
        if not results_path.is_dir():
            click.secho(f"Error: results dir not found: {results_path}", fg="red")
            sys.exit(1)
        click.secho(f"Scanning {results_path} -> {work_path}...", fg="blue")
        cache = build_work_dir_cache(results_path, work_path)
        purge_set = cache["purge"]
        click.secho(
            f"Found {len(cache['protect'])} referenced, "
            f"{len(purge_set)} purgeable work hash directories.",
            fg="blue",
        )

    # Partition the purge candidates: skip any that are write-protected.
    purgeable: list[Path] = []
    guarded: list[Path] = []
    for hd in sorted(purge_set):
        if str(hd) in protected:
            guarded.append(hd)
        else:
            purgeable.append(hd)

    if guarded:
        click.secho(
            f"Skipping {len(guarded)} protected-but-unreferenced "
            f"director{'y' if len(guarded) == 1 else 'ies'} "
            f"(run `gosh protect --unprotect` first to allow removal).",
            fg="yellow",
        )

    if not purgeable:
        click.secho("Nothing to purge.", fg="green")
        return

    # Display the list of dirs to be removed.
    prefix = "[dry-run] " if dry_run else ""
    total_size = 0
    for hd in purgeable:
        rel = hd.relative_to(work_path)
        if show_size:
            sz = _dir_size(hd)
            total_size += sz
            click.echo(f"  {prefix}{rel}  ({_human_size(sz)})")
        else:
            click.echo(f"  {prefix}{rel}")

    n = len(purgeable)
    dirs_word = "directory" if n == 1 else "directories"
    size_str = f"  (~{_human_size(total_size)} on disk)" if show_size else ""
    click.secho(f"\n{n} {dirs_word} to purge{size_str}.", fg="yellow")

    if dry_run:
        click.secho("Dry run — nothing deleted.", fg="blue")
        return

    if not yes:
        click.confirm("Proceed with deletion?", abort=True)

    removed = 0
    errors = 0
    for hd in purgeable:
        try:
            _rmtree_writable(hd)
            removed += 1
        except OSError as e:
            click.secho(f"Warning: could not remove {hd}: {e}", fg="yellow")
            errors += 1

    # Remove any 2-char prefix dirs that are now empty.
    for prefix_dir in sorted(work_path.iterdir()):
        if prefix_dir.is_dir() and not prefix_dir.is_symlink():
            try:
                prefix_dir.rmdir()  # no-op unless the directory is empty
            except OSError:
                pass

    removed_word = "directory" if removed == 1 else "directories"
    click.secho(
        f"Removed {removed} {removed_word}; {errors} error(s).",
        fg="green" if errors == 0 else "yellow",
    )
