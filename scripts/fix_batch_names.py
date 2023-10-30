#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import os
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, List, TypeVar, Iterator


T = TypeVar("T")


def eprint(*args: object, **kwargs: object) -> None:
    kwargs["file"] = sys.stderr
    print(*args, **kwargs)


def fragment(size: int, items: list[T]) -> Iterator[list[T]]:
    """Faster alternative to grouper for lists/strings."""
    return (items[i : i + size] for i in range(0, len(items), size))


def sha256(items: list[str]) -> str:
    """Calculates the SHA256 hash for a set of values"""
    hasher = hashlib.sha256()
    for it in items:
        if isinstance(it, Path):
            it = str(it)

        hasher.update(it.encode("utf-8"))
    return hasher.hexdigest()


def _collect_fast5s(source, groups=None):
    # fast5s are grouped by dir, to ensure that adding/removing whole folders
    # does not result in the re-processing of existing folders
    if groups is None:
        groups = defaultdict(list)

    for it in source.iterdir():
        if it.is_dir():
            _collect_fast5s(it, groups)
        elif it.suffix.lower() in (".fast5", ".pod5"):
            groups[it.parent].append(str(it))

    return groups.values()


def remove_roots(root: str, filenames: list[str]) -> list[str]:
    # Exclude the source folder from the hashed paths, to
    # prevent minor changes in folder structure from causing
    # files to be re-run.
    result = []
    for filename in filenames:
        filepath = Path(filename)
        assert filepath.is_relative_to(root)

        result.append(filepath.relative_to(root))

    return result


def _collect_batches(
    destination: str, source: str, batch_size: int = 90
) -> Iterator[tuple[str, list[str]]]:
    for it in sorted(Path(source).iterdir()):
        for group in _collect_fast5s(it):
            group.sort()  # Ensure stable batches even filesystem order changes

            if batch_size is None:
                yield it.name, group
            else:
                for batch in fragment(batch_size, group):
                    yield it.name, batch


class HelpFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self, *args: Any, **kwargs: Any) -> None:
        kwargs.setdefault("width", 79)

        super().__init__(*args, **kwargs)


def parse_args(argv: List[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(formatter_class=HelpFormatter)
    parser.add_argument("source", type=Path)
    parser.add_argument("destination", type=Path)
    parser.add_argument("--commit", action="store_true")
    parser.add_argument("--batch-size", default=90, type=int)
    parser.add_argument("--actual-source")

    return parser.parse_args(argv)


def main(argv: List[str]) -> int:
    args = parse_args(argv)
    source: str = args.source
    destination: str = args.destination
    batch_size: int = args.batch_size
    commit: bool | None = args.commit
    actual_source: str = source
    if args.actual_source is not None:
        actual_source = args.actual_source

    tasks: list[tuple[str, str]] = []
    for sample, batch in _collect_batches(
        destination=destination,
        source=source,
        batch_size=batch_size,
    ):
        normalized_batch = remove_roots(source, batch)
        normalized_hash = sha256(normalized_batch)

        original_batch = [os.path.join(actual_source, it) for it in normalized_batch]
        original_hash = sha256(original_batch)

        cache = os.path.join(destination, f"{sample}.cache")
        for ext in (".fast5s", ".model.txt", ".fq.gz"):
            source_filename = os.path.join(cache, f"{original_hash}{ext}")
            source_exists = os.path.exists(source_filename)

            destination_filename = os.path.join(cache, f"{normalized_hash}{ext}")
            destination_exists = os.path.exists(destination_filename)

            if source_exists and not destination_exists:
                tasks.append((source_filename, destination_filename))
            elif source_exists and destination_exists:
                eprint(f"WARNING: Duplicate source/destination: {source_filename!r}")
            elif not (source_exists or destination_exists):
                eprint(f"WARNING: Source not found: {source_filename!r}")

    for source_filename, destination_filename in tasks:
        eprint(f"Renaming {source_filename!r}")
        eprint(f"      to {destination_filename!r}")
        if commit:
            Path(source_filename).rename(destination_filename)

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
