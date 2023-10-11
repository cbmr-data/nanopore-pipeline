#!/usr/bin/env python3
# -*- coding: utf8 -*-
# pyright: strict
from __future__ import annotations

import argparse
import hashlib
import random
import sys
from pathlib import Path
from collections import defaultdict
from typing import (
    Any,
    List,
    Optional,
    Dict,
    TypeVar,
    Sequence,
    Iterable,
    Iterator,
    Tuple,
    NoReturn,
)
from dataclasses import dataclass

import pysam
from tqdm import tqdm


@dataclass
class BAMStats:
    sampled_batches: int = 0
    completed_batches: int = 0
    genome_size: int = 0
    short_reads_excluded: int = 0
    unmapped_reads_excluded: int = 0
    long_reads_included: int = 0
    total_query_length: int = 0
    bases_mapped: int = 0
    bases_deleted: int = 0
    bases_inserted: int = 0
    bases_clipped: int = 0

    def estimate(self, num_batches: int) -> BAMStats:
        if not self.sampled_batches:
            return self

        m = num_batches / self.sampled_batches
        return BAMStats(
            sampled_batches=self.sampled_batches,
            completed_batches=num_batches,
            genome_size=self.genome_size,
            short_reads_excluded=int(self.short_reads_excluded * m),
            unmapped_reads_excluded=int(self.unmapped_reads_excluded * m),
            long_reads_included=int(self.long_reads_included * m),
            total_query_length=int(self.total_query_length * m),
            bases_mapped=int(self.bases_mapped * m),
            bases_deleted=int(self.bases_deleted * m),
            bases_inserted=int(self.bases_inserted * m),
            bases_clipped=int(self.bases_clipped * m),
        )


def eprint(*args: object) -> None:
    print(*args, file=sys.stderr)


def warning(*args: object) -> None:
    eprint("WARNING:", *args)


def collect_result_stats(
    root: Path,
    batches: Iterable[str],
    min_length: int,
    sample_bams: int = -1,
    sample_name: str = "unknown sample",
) -> BAMStats:
    stats = BAMStats()
    batches = frozenset(batches)
    if not root.exists():
        return stats

    available_batches: List[Path] = []
    for filepath in root.iterdir():
        batch_key, suffix = filepath.name.lower().split(".", 1)
        if suffix not in ("bam", "fq.gz"):
            continue
        elif batch_key not in batches:
            warning("unknown batch", filepath)
            continue

        if suffix == "bam":
            available_batches.append(filepath)
        elif suffix == "fq.gz":
            stats.completed_batches += 1

    if 0 <= sample_bams < len(available_batches):
        available_batches = random.sample(available_batches, k=sample_bams)

    for idx, filepath in enumerate(sorted(available_batches), start=1):
        with pysam.AlignmentFile(str(filepath), threads=4) as handle:
            if not stats.genome_size:
                stats.genome_size = sum(handle.lengths)

            cigars = [0] * 9
            desc = f"{sample_name} [{idx}/{len(available_batches)}]"
            for record in tqdm(handle, desc=desc, unit_scale=True):
                # Skip supplementary and alternative alignments
                if record.flag & 0x900:
                    continue

                stats.total_query_length += record.query_length
                if record.query_length < min_length:
                    stats.short_reads_excluded += 1
                    continue
                elif record.is_unmapped:
                    stats.unmapped_reads_excluded += 1
                    continue
                else:
                    stats.long_reads_included += 1

                for op, length in record.cigartuples:
                    cigars[op] += length

            stats.bases_mapped += cigars[0] + cigars[7] + cigars[8]
            stats.bases_deleted += cigars[2]
            stats.bases_inserted += cigars[1]
            stats.bases_clipped += cigars[4] + cigars[5]
            stats.sampled_batches += 1

    return stats


def collect_batch_sizes(batches: Dict[str, List[Path]]) -> Dict[str, int]:
    sizes: Dict[str, int] = {}
    for key, filepaths in batches.items():
        raw_size = 0
        for it in filepaths:
            raw_size += it.stat().st_size

        sizes[key] = raw_size

    return sizes


T = TypeVar("T")


def fragment(size: int, items: Sequence[T]) -> Iterator[Sequence[T]]:
    """Faster alternative to grouper for lists/strings."""
    return (items[i : i + size] for i in range(0, len(items), size))


def sha256(items: Iterable[Path]) -> str:
    """Calculates the SHA256 hash for a set of values"""
    hasher = hashlib.sha256()
    for it in items:
        hasher.update(str(it).encode("utf-8"))
    return hasher.hexdigest()


def _collect_fast5s(
    source: Path,
    groups: Optional[Dict[Path, List[Path]]] = None,
    tqdm: Optional[tqdm[NoReturn]] = None,
) -> List[List[Path]]:
    # fast5s are grouped by dir, to ensure that adding/removing whole folders
    # does not result in the re-processing of existing folders
    if groups is None:
        groups = defaultdict(list)

    for it in source.iterdir():
        if tqdm is not None:
            tqdm.update(1)
        if it.is_dir():
            _collect_fast5s(it, groups, tqdm)
        elif it.suffix.lower() in (".fast5", ".pod5"):
            groups[it.parent].append(it)

    return list(groups.values())


def collect_batches(
    source: str,
    batch_size: int = 25,
) -> Dict[str, Dict[str, List[Path]]]:
    samples: Dict[str, Dict[str, List[Path]]] = {}

    with tqdm() as progress:
        for it in Path(source).iterdir():
            samples[it.name] = sample = {}

            for group in _collect_fast5s(it, tqdm=progress):
                group.sort()  # Ensure stable batches even filesystem order changes

                for batch in fragment(batch_size, group):
                    sample[sha256(batch)] = batch

    return samples


class HelpFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self, *args: Any, **kwargs: Any) -> None:
        kwargs.setdefault("width", 79)

        super().__init__(*args, **kwargs)


def parse_args(argv: List[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(formatter_class=HelpFormatter)
    parser.add_argument("batch", type=Path)
    parser.add_argument("--batches-root", type=Path, default=Path("batches"))
    parser.add_argument("--results-root", type=Path, default=Path("results"))
    parser.add_argument("--min-length", type=int, default=500)
    parser.add_argument("--batch-size", type=int, default=25)
    parser.add_argument("--sample-bams", type=int, default=5)

    return parser.parse_args(argv)


def main(argv: List[str]) -> int:
    pysam.set_verbosity(0)
    args = parse_args(argv)

    batches_root = args.batches_root / args.batch
    eprint("collecting batches")
    batches = collect_batches(batches_root, batch_size=args.batch_size)
    num_batches = sum(map(len, batches.values()))
    eprint("Found", len(batches), "samples and", num_batches, "batches")

    print(
        "Sample",
        "RawSize",
        "Batches",
        "BatchesDonePct",
        "BatchesSamples",
        "EstCoverage",
        "MaxCoverage",
        "ShortReadsExcluded",
        "UnmappedReadsExcluded",
        "LongReadsIncluded",
        "QueryLengthGbp",
        "MappedGbp",
        "MappedPct",
        "DeletedGbp",
        "DeletedPct",
        "InsertedGbp",
        "InsertedPct",
        "ClippedGbp",
        "ClippedPct",
        sep="\t",
    )

    for sample_name, sample_batches in sorted(batches.items()):
        eprint("processing batches for", sample_name)
        results_root = args.results_root / args.batch / f"{sample_name}.cache"
        batch_sizes = collect_batch_sizes(sample_batches)

        stats = collect_result_stats(
            root=results_root,
            batches=sample_batches,
            min_length=args.min_length,
            sample_bams=args.sample_bams,
            sample_name=sample_name,
        )
        est_stats = stats.estimate(len(sample_batches))

        def _fmt(n: int) -> Tuple[str, str]:
            return (
                f"{n / 1e9:.1f}",
                f"{(n * 100) / est_stats.total_query_length:.1f}"
                if est_stats.total_query_length
                else "NA",
            )

        print(
            sample_name,
            sum(batch_sizes.values()),
            len(sample_batches),
            f"{(100 * stats.completed_batches) / len(sample_batches):.1f}",
            stats.sampled_batches,
            f"{est_stats.bases_mapped / stats.genome_size:.1f}"
            if stats.genome_size
            else "NA",
            f"{est_stats.total_query_length / stats.genome_size:.1f}"
            if stats.genome_size
            else "NA",
            est_stats.short_reads_excluded,
            est_stats.unmapped_reads_excluded,
            est_stats.long_reads_included,
            f"{est_stats.total_query_length / 1e9:.1f}",
            *_fmt(est_stats.bases_mapped),
            *_fmt(est_stats.bases_deleted),
            *_fmt(est_stats.bases_inserted),
            *_fmt(est_stats.bases_clipped),
            sep="\t",
        )

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
