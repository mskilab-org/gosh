"""
gosh protect — write-protect Nextflow work dirs referenced from a results directory.

Scans a results directory for symlinks into a work directory (transitively,
following work-dir-to-work-dir references), then removes write permission
from each referenced work-hash directory.

Reversible: original modes are recorded in <work_dir>/.gosh_protected.json
so `gosh protect --unprotect` restores them exactly.
"""

import json
import os
import stat
import sys
from datetime import datetime
from pathlib import Path

import click


STATE_FILENAME = ".gosh_protected.json"


def _load_state(state_path: Path) -> dict:
    if not state_path.exists():
        return {}
    try:
        with open(state_path) as f:
            data = json.load(f)
        if not isinstance(data, dict):
            click.secho(
                f"Warning: {state_path} is not a JSON object; ignoring.", fg="yellow"
            )
            return {}
        return data
    except (json.JSONDecodeError, OSError) as e:
        click.secho(f"Warning: could not read {state_path}: {e}", fg="yellow")
        return {}


def _save_state(state_path: Path, state: dict) -> None:
    tmp = state_path.with_name(state_path.name + ".tmp")
    with open(tmp, "w") as f:
        json.dump(state, f, indent=2, sort_keys=True)
    os.replace(tmp, state_path)


def _remove_write_bits(mode: int) -> int:
    return mode & ~(stat.S_IWUSR | stat.S_IWGRP | stat.S_IWOTH)


def _iter_targets(hash_dir: Path, recursive: bool):
    """Yield the hash dir itself, plus (if recursive) all non-symlink descendants."""
    yield hash_dir
    if not recursive:
        return
    for dirpath, dirnames, filenames in os.walk(hash_dir, followlinks=False):
        dp = Path(dirpath)
        for name in filenames + dirnames:
            p = dp / name
            if not p.is_symlink():
                yield p


@click.command(name="protect")
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
    help="Nextflow work directory whose hash dirs to protect.",
)
@click.option(
    "-u",
    "--unprotect",
    is_flag=True,
    help="Restore recorded permissions instead of removing write access.",
)
@click.option(
    "-R",
    "--recursive",
    is_flag=True,
    default=True,
    show_default=True,
    help="Also chmod every non-symlink file/dir inside each hash dir.",
)
@click.option(
    "-n",
    "--dry-run",
    is_flag=True,
    help="Show what would change without touching the filesystem.",
)
def protect_cli(results_dir, work_dir, unprotect, recursive, dry_run):
    """
    Write-protect work directories referenced from the results directory.

    With --unprotect, restores the modes recorded when protect was run.
    """
    from ..utils.find_work_deps import find_work_dirs

    work_path = Path(work_dir).resolve()
    if not work_path.is_dir():
        click.secho(f"Error: work dir not found: {work_path}", fg="red")
        sys.exit(1)

    state_path = work_path / STATE_FILENAME
    state = _load_state(state_path)

    if unprotect:
        targets = sorted(state.keys())
        if not targets:
            click.secho(f"No protected paths recorded in {state_path}.", fg="yellow")
            return
        click.secho(f"Restoring {len(targets)} paths from {state_path}...", fg="blue")
    else:
        results_path = Path(results_dir).resolve()
        if not results_path.is_dir():
            click.secho(f"Error: results dir not found: {results_path}", fg="red")
            sys.exit(1)
        click.secho(f"Scanning {results_path} -> {work_path}...", fg="blue")
        all_dirs, _ = find_work_dirs(results_path, work_path)
        click.secho(
            f"Found {len(all_dirs)} referenced work hash directories.", fg="blue"
        )
        targets = [
            str(p) for hd in sorted(all_dirs) for p in _iter_targets(hd, recursive)
        ]

    changed = 0
    skipped = 0
    missing = 0
    errors = 0

    for key in targets:
        target = Path(key)

        if unprotect:
            entry = state.get(key)
            if entry is None:
                skipped += 1
                continue
            orig_mode = entry["mode"]
            if dry_run:
                click.echo(f"[dry-run] restore {oct(orig_mode)} {target}")
                changed += 1
                continue
            try:
                os.chmod(target, orig_mode)
                del state[key]
                changed += 1
            except FileNotFoundError:
                del state[key]
                missing += 1
            except OSError as e:
                click.secho(f"Warning: chmod failed on {target}: {e}", fg="yellow")
                errors += 1
            continue

        # Protect path
        try:
            st = os.stat(target, follow_symlinks=False)
        except FileNotFoundError:
            missing += 1
            continue
        except OSError as e:
            click.secho(f"Warning: stat failed on {target}: {e}", fg="yellow")
            errors += 1
            continue

        cur_mode = stat.S_IMODE(st.st_mode)
        new_mode = _remove_write_bits(cur_mode)
        if new_mode == cur_mode:
            skipped += 1
            continue

        if key not in state:
            state[key] = {
                "mode": cur_mode,
                "protected_at": datetime.now().isoformat(timespec="seconds"),
            }

        if dry_run:
            click.echo(f"[dry-run] chmod {oct(cur_mode)}->{oct(new_mode)} {target}")
            changed += 1
            continue
        try:
            os.chmod(target, new_mode)
            changed += 1
        except OSError as e:
            click.secho(f"Warning: chmod failed on {target}: {e}", fg="yellow")
            errors += 1

    if not dry_run:
        try:
            if state:
                _save_state(state_path, state)
            elif state_path.exists():
                state_path.unlink()
        except OSError as e:
            click.secho(
                f"Warning: could not update state file {state_path}: {e}", fg="yellow"
            )

    verb = "would change" if dry_run else ("restored" if unprotect else "protected")
    click.secho(
        f"{verb} {changed}; skipped {skipped}; missing {missing}; errors {errors}. "
        f"State: {state_path}",
        fg="green" if errors == 0 else "yellow",
    )
