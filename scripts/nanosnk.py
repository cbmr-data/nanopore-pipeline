from __future__ import annotations

import hashlib
import os
import sys
from collections import defaultdict
from collections.abc import Iterable, Iterator
from dataclasses import dataclass, field
from pathlib import Path
from typing import TypeVar

T = TypeVar("T")

__all__ = [
    "abort",
    "collect_samples",
    "collect_samples_for_genotyping",
]


@dataclass
class Chunk:
    sample: str
    batch: str
    name: str

    # FIXME: Sort? Should be sorted by default
    files: list[str] = field(default_factory=list)


@dataclass
class Batch:
    sample: str
    name: str
    chunks: dict[str, Chunk] = field(default_factory=dict)


@dataclass
class Sample:
    name: str
    batches: dict[str, Batch] = field(default_factory=dict)

    def fastq_chunks(self, destination: str, batch: str) -> list[str]:
        filenames: list[str] = []
        for _, chunk in sorted(self.batches[batch].chunks.items()):
            filenames.append(
                os.path.join(
                    destination,
                    "reads",
                    chunk.sample,
                    batch,
                    f"{chunk.name}.fq.gz",
                )
            )

        return filenames

    def bam_batches(self, destination: str, genome: str, kind: str) -> list[str]:
        return [
            os.path.join(
                destination,
                "alignments",
                batch.sample,
                f"{batch.name}.{genome}.{kind}.bam",
            )
            for _, batch in sorted(self.batches.items())
        ]


def abort(*args: object):
    print("ERROR:", *args, file=sys.stderr)
    sys.exit(1)


def _fragment(size: int | None, items: list[T]) -> Iterator[list[T]]:
    """Faster alternative to grouper for lists/strings."""
    if size is None:
        yield items
    else:
        for i in range(0, len(items), size):
            yield items[i : i + size]


def _sha256(items: Iterable[str]) -> str:
    """Calculates the SHA256 hash for a set of values"""
    hasher = hashlib.sha256()
    for it in items:
        hasher.update(it.encode("utf-8"))
    return hasher.hexdigest()


def _collect_fast5s(
    source: Path,
    groups: dict[Path, list[str]] | None = None,
) -> list[list[str]]:
    # fast5s are grouped by dir, to ensure that adding/removing whole folders
    # does not result in the re-processing of existing folders
    if groups is None:
        groups = defaultdict(list)

    for it in source.iterdir():
        if it.is_dir():
            _collect_fast5s(it, groups)
        elif it.suffix.lower() in (".fast5", ".pod5"):
            groups[it.parent].append(str(it))

    return [value for _, value in sorted(groups.items())]


def _hash_source_files(source: str, filenames: Iterable[str]) -> str:
    # Exclude the source folder from the hashed paths, to prevent global
    # changes in folder structure (e.g. Esrum v. Computerome) from
    # causing files to be re-run.
    result: list[str] = []
    for filename in filenames:
        filepath = Path(filename)
        assert filepath.is_relative_to(source)

        result.append(str(filepath.relative_to(source)))

    return _sha256(sorted(result))


def collect_samples(source: str, batch_size: int | None = 90) -> dict[str, Sample]:
    samples: dict[str, Sample] = {}
    for it in sorted(Path(source).iterdir()):
        samples[it.name] = sample = Sample(name=it.name)

        for batch_dir in it.iterdir():
            if not batch_dir.is_dir():
                abort("Batch is not a folder:", batch_dir)

            batch = sample.batches[batch_dir.name] = Batch(
                sample=it.name,
                name=batch_dir.name,
            )

            for group in _collect_fast5s(batch_dir):
                group.sort()  # Ensure stable batches even if filesystem order changes

                for chunk in _fragment(batch_size, group):
                    namehash = _hash_source_files(source, chunk)
                    batch.chunks[namehash] = Chunk(
                        sample=it.name,
                        batch=batch_dir.name,
                        name=namehash,
                        files=chunk,
                    )

    return samples


def collect_samples_for_genotyping(
    samples: Iterable[str],
    blacklist: Iterable[str],
) -> tuple[str, ...]:
    unknown_samples = set(blacklist).difference(samples)
    if unknown_samples:
        raise AssertionError(f"Unexpected samples in blacklist: {unknown_samples}")

    return tuple(sorted(set(samples).difference(blacklist)))
