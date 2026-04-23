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


def find_work_dirs(results_dir: Path, work_dir: Path) -> dict[Path, set[Path]]:
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
        """Walk scan_root, collect work hash dirs referenced by any symlink."""
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

        deps = scan_dir_for_symlinks(current)
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
