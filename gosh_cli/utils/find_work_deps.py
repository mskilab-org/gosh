#!/usr/bin/env python3
"""
Find all work directories referenced (directly or transitively) by symlinks
in a Nextflow results directory.
"""

import os
import sys
from pathlib import Path
from collections import defaultdict


def resolve_to_work_hash_dir(path: Path, work_dir: Path) -> Path | None:
    """
    Given a resolved symlink target, return the work hash directory
    (the 2-level hash dir, e.g. work/ab/cdef1234...) if it lives under work_dir.
    """
    try:
        rel = path.relative_to(work_dir)
        parts = rel.parts
        if len(parts) >= 2:
            return work_dir / parts[0] / parts[1]
    except ValueError:
        pass
    return None


def read_symlink_target(entry: Path) -> Path | None:
    """
    Return the symlink's target as an absolute, normalized path, without
    requiring the target to exist. Robust to broken relative symlinks.

    Unlike Path.resolve(), this does not chase chains or touch the filesystem
    beyond a single readlink() — so dangling links don't raise and can still
    be classified against the work dir.
    """
    try:
        raw = os.readlink(entry)
    except OSError:
        return None
    target = Path(raw)
    if not target.is_absolute():
        target = entry.parent / target
    return Path(os.path.normpath(str(target)))


def enumerate_all_work_hash_dirs(work_dir: Path) -> set[Path]:
    """
    Return the set of every two-level hash directory that physically exists
    under work_dir (e.g. work/ab/cdef1234...).

    This is a shallow two-level scan — no recursion into hash dir contents —
    so it is fast regardless of how many files live inside each hash dir.
    Failed, cached, and completed dirs are all included; callers subtract the
    reachable set to find candidates for purging.
    """
    result: set[Path] = set()
    try:
        with os.scandir(work_dir) as top_it:
            for prefix_entry in top_it:
                if not prefix_entry.is_dir(follow_symlinks=False):
                    continue
                try:
                    with os.scandir(prefix_entry.path) as hash_it:
                        for hash_entry in hash_it:
                            if hash_entry.is_dir(follow_symlinks=False):
                                result.add(Path(hash_entry.path))
                except OSError:
                    pass
    except OSError:
        pass
    return result


def verify_cache(
    results_dir: Path,
    purge_set: set[Path],
    work_dir: Path,
) -> dict[Path, list[Path]]:
    """
    Orthogonally verify that no symlink in *results_dir* points into a work
    hash directory that is marked for purging.

    This function is intentionally **independent** of ``find_work_dirs`` and
    ``build_work_dir_cache``.  It performs a plain ``os.walk`` of the results
    directory — a completely separate code path — so it can surface bugs in
    the BFS logic itself, not just confirm its own output.

    Algorithm
    ---------
    For every symlink encountered anywhere under *results_dir*:

    1. Read the raw link target with a single ``os.readlink`` (no chain
       following, no existence check).
    2. Normalize to an absolute path relative to the symlink's parent.
    3. Resolve to the two-level work hash dir (``work/ab/cdef…``) if the
       target lives under *work_dir*.
    4. If that hash dir is present in *purge_set*, record the violation.

    Parameters
    ----------
    results_dir:
        Root of the Nextflow publishDir to scan.
    purge_set:
        The set of work hash dirs classified as purgeable by the cache.
        Typically ``cache["purge"]`` from ``build_work_dir_cache`` or a
        loaded ``.gosh_cache.json``.
    work_dir:
        Resolved path to the Nextflow ``work/`` directory.

    Returns
    -------
    A dict mapping each *violating work hash dir* to the list of result
    symlinks that reference it.  An empty dict means the cache is clean.
    """
    results_dir = results_dir.resolve()
    work_dir = work_dir.resolve()

    # violations: work_hash_dir -> [symlink_path, ...]
    violations: dict[Path, list[Path]] = {}

    for dirpath, dirnames, filenames in os.walk(results_dir, followlinks=False):
        dp = Path(dirpath)
        for name in filenames + dirnames:
            entry = dp / name
            if not entry.is_symlink():
                continue
            target = read_symlink_target(entry)
            if target is None:
                continue
            wdir = resolve_to_work_hash_dir(target, work_dir)
            if wdir is None:
                continue
            if wdir in purge_set:
                violations.setdefault(wdir, []).append(entry)

    return violations


def build_work_dir_cache(
    results_dir: Path,
    work_dir: Path,
) -> dict[str, set[Path]]:
    """
    Build a unified cache that partitions every work hash directory into two
    buckets in a single coordinated pass:

    ``"protect"``
        Work hash dirs that are reachable (directly or transitively) from
        symlinks in *results_dir*.  These correspond to processes that
        completed and whose outputs were published.  They should be
        write-protected so accidental re-runs or ``nextflow clean`` cannot
        delete them.

    ``"purge"``
        Work hash dirs that exist on disk but are *not* reachable from
        results.  These are failed attempts, cached-but-unpublished
        intermediates, or orphaned retries.  They are safe to delete once
        the pipeline is complete.

    Why this is faster than two separate scans
    ------------------------------------------
    Previously, finding the "purge" candidates would require walking the
    entire work directory tree — including the contents of every hash dir —
    in addition to the BFS that builds the reachable set.  This function
    avoids that by:

    1. Running the results-dir BFS (``find_work_dirs``) to get the reachable
       set.  The BFS only ever calls ``os.scandir`` at the *top level* of
       each work hash dir, not ``os.walk`` into its contents.
    2. Running a two-level shallow scan of *work_dir* itself
       (``enumerate_all_work_hash_dirs``) to learn every hash dir that exists
       on disk.  This is two levels of ``os.scandir`` with no recursion into
       hash dir contents, so it is O(number_of_hash_dirs) in syscalls rather
       than O(total_files_in_work_dir).
    3. Computing ``purge = all_existing - reachable`` with a plain set
       difference — no extra I/O at all.

    Parameters
    ----------
    results_dir:
        Resolved path to the Nextflow ``publishDir`` / results root.
    work_dir:
        Resolved path to the Nextflow ``work/`` directory.

    Returns
    -------
    dict with keys ``"protect"`` and ``"purge"``, each mapping to a
    ``set[Path]`` of absolute two-level work hash directory paths.
    """
    results_dir = results_dir.resolve()
    work_dir = work_dir.resolve()

    # Step 1: BFS from results → reachable (completed) work hash dirs.
    reachable, _ = find_work_dirs(results_dir, work_dir)

    # Step 2: Shallow two-level enumeration of every hash dir on disk.
    all_existing = enumerate_all_work_hash_dirs(work_dir)

    # Step 3: Anything on disk but not reachable from results is a purge candidate.
    return {
        "protect": reachable,
        "purge": all_existing - reachable,
    }


def find_work_dirs(results_dir: Path, work_dir: Path) -> tuple[set[Path], dict[Path, set[Path]]]:
    """
    Recursively find all work hash directories referenced by symlinks,
    starting from results_dir and following transitive dependencies.

    Returns a dict mapping each work hash dir to the set of work hash dirs
    it depends on (its own staged-input symlinks).
    """
    results_dir = results_dir.resolve()
    work_dir = work_dir.resolve()

    visited_work_dirs: set[Path] = set()
    # maps work_hash_dir -> set of work_hash_dirs it pulls from
    dependency_map: dict[Path, set[Path]] = defaultdict(set)

    queue: list[Path] = [results_dir]

    def scan_dir_for_symlinks(scan_root: Path) -> set[Path]:
        """Recursively walk scan_root, collect work hash dirs referenced by any symlink.
        Used for the results directory where outputs may be nested arbitrarily deep."""
        found = set()
        for dirpath, dirnames, filenames in os.walk(scan_root, followlinks=False):
            for name in filenames + dirnames:
                entry = Path(dirpath) / name
                if not entry.is_symlink():
                    continue
                target = read_symlink_target(entry)
                if target is None:
                    continue
                wdir = resolve_to_work_hash_dir(target, work_dir)
                if wdir is not None:
                    found.add(wdir)
        return found

    def scan_work_dir_for_symlinks(scan_root: Path) -> set[Path]:
        """Scan only the top level of a work hash dir for symlinks into work_dir.
        Nextflow always stages inputs as top-level symlinks in the hash dir, so
        there is no need to recurse — avoiding an expensive os.walk() over logs,
        scripts, and output files for every (including failed) work directory."""
        found = set()
        try:
            with os.scandir(scan_root) as it:
                for entry in it:
                    if not entry.is_symlink():
                        continue
                    target = read_symlink_target(Path(entry.path))
                    if target is None:
                        continue
                    wdir = resolve_to_work_hash_dir(target, work_dir)
                    if wdir is not None:
                        found.add(wdir)
        except OSError:
            pass
        return found

    # Seed: scan the results directory
    initial = scan_dir_for_symlinks(results_dir)
    for w in initial:
        if w not in visited_work_dirs:
            visited_work_dirs.add(w)
            queue.append(w)

    # BFS over work dirs, following their internal symlinks transitively
    processed: set[Path] = {results_dir}
    while queue:
        current = queue.pop()
        if current in processed:
            continue
        processed.add(current)

        if not current.exists():
            continue

        deps = scan_work_dir_for_symlinks(current)
        for dep in deps:
            dependency_map[current].add(dep)
            if dep not in visited_work_dirs:
                visited_work_dirs.add(dep)
                queue.append(dep)

    return visited_work_dirs, dict(dependency_map)


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Find all Nextflow work dirs referenced from a results directory."
    )
    parser.add_argument("results_dir", help="Path to publishDir / results directory")
    parser.add_argument("work_dir", help="Path to Nextflow work/ directory")
    parser.add_argument(
        "--show-deps", action="store_true", help="Show per-work-dir dependency graph"
    )
    parser.add_argument(
        "--show-purge",
        action="store_true",
        help="Also list work dirs not reachable from results (purge candidates).",
    )
    args = parser.parse_args()

    results_dir = Path(args.results_dir).resolve()
    work_dir = Path(args.work_dir).resolve()

    if not results_dir.exists():
        print(f"ERROR: results dir not found: {results_dir}", file=sys.stderr)
        sys.exit(1)
    if not work_dir.exists():
        print(f"ERROR: work dir not found: {work_dir}", file=sys.stderr)
        sys.exit(1)

    print(f"Scanning: {results_dir}", file=sys.stderr)
    print(f"Work dir: {work_dir}", file=sys.stderr)

    if args.show_purge:
        cache = build_work_dir_cache(results_dir, work_dir)
        protect_dirs = cache["protect"]
        purge_dirs = cache["purge"]
        print(
            f"\nFound {len(protect_dirs)} referenced (protect) + "
            f"{len(purge_dirs)} unreachable (purge) work directories.\n",
            file=sys.stderr,
        )
        print("# protect", file=sys.stderr)
        for d in sorted(protect_dirs):
            print(d.relative_to(work_dir))
        print("\n# purge", file=sys.stderr)
        for d in sorted(purge_dirs):
            print(d.relative_to(work_dir))
        return

    all_work_dirs, dep_map = find_work_dirs(results_dir, work_dir)

    print(
        f"\nFound {len(all_work_dirs)} referenced work directories:\n", file=sys.stderr
    )
    for d in sorted(all_work_dirs):
        print(d.relative_to(work_dir))

    if args.show_deps and dep_map:
        print(
            "\n--- Dependency graph (work dir -> deps it pulls from) ---",
            file=sys.stderr,
        )
        for src, deps in sorted(dep_map.items()):
            for dep in sorted(deps):
                print(
                    f"  {src.relative_to(work_dir)} -> {dep.relative_to(work_dir)}",
                    file=sys.stderr,
                )


if __name__ == "__main__":
    main()
