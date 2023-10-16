#!/usr/bin/env python3
# -*- coding: utf8 -*-
# pyright: strict
import argparse
import json
import os
import random
import sys
from collections import Counter, defaultdict
from itertools import product
from pathlib import Path
from typing import Any, List, Dict, Tuple, Optional, Iterable, TypeVar
from typing_extensions import TypedDict

import pysam

QueryLens = TypedDict(
    "QueryLens",
    {
        "mapped": Dict[int, int],
        "unmapped": Dict[int, int],
        "supplementary": Dict[int, int],
        "secondary": Dict[int, int],
    },
)

Statistics = TypedDict(
    "Statistics",
    {
        "query_len": QueryLens,
        "mapq": Dict[Optional[int], int],
        "qs": Dict[int, int],
        "n": Dict[int, int],
        "cigar": Dict[str, Dict[int, int]],
        "alignment": Dict[str, int],
        "pct_mapped": Dict[int, int],
        "pct_mapped_mm": Dict[int, int],
    },
)

T = TypeVar("T")


def progress(iterable: Iterable[T], desc: Optional[str] = None) -> Iterable[T]:
    try:
        if sys.stderr.isatty():
            import tqdm

            return tqdm.tqdm(iterable, unit_scale=True, desc=desc)
    except ImportError:
        pass

    return iterable


class ReservoirSample:
    def __init__(self, downsample_to: int) -> None:
        self._downsample_to = downsample_to
        self._rng = random.Random()
        self.items: List[Any] = []
        self._index = 0

        if downsample_to < 0:
            raise ValueError(
                "Negative value for 'downsample_to': %i" % (downsample_to,)
            )

    def add(self, item: Any) -> None:
        if self._index >= self._downsample_to:
            index = self._rng.randint(0, self._index)
            if index < self._downsample_to:
                self.items[index] = item
        else:
            self.items.append(item)

        self._index += 1


def new_sample() -> Statistics:
    return {
        "query_len": {
            "mapped": Counter(),
            "unmapped": Counter(),
            "supplementary": Counter(),
            "secondary": Counter(),
        },
        "mapq": Counter(),
        "qs": Counter(),
        "n": Counter(),
        "cigar": {key: Counter() for key in "MIDNSHP"},
        "alignment": {
            f"{key1}/{key2}": 0
            for (key1, key2) in product("NACGT-", repeat=2)
            if key1 != "-" or key2 != "-"
        },
        "pct_mapped": Counter(),
        "pct_mapped_mm": Counter(),
    }


def flatten_counts(value: Any) -> Any:
    if isinstance(value, Counter):
        if all(isinstance(key, int) for key in value):
            return {"x": list(value), "n": list(value.values())}
    elif isinstance(value, dict):
        return {key: flatten_counts(value) for key, value in value.items()}
    else:
        return value


def collect_basic_stats(
    sample: str,
    stats: Statistics,
    subsample: ReservoirSample,
    filename: Path,
    nthreads: int = 4,
) -> None:
    with pysam.AlignmentFile(os.fspath(filename), threads=nthreads) as handle:
        for record in progress(handle, desc=f"{sample} [{filename}]"):
            alignment_type = "mapped"
            if record.is_unmapped:
                alignment_type = "unmapped"
            elif record.is_supplementary:
                alignment_type = "supplementary"
            elif record.is_secondary:
                alignment_type = "secondary"

            # Query lengths
            stats["query_len"][alignment_type][record.query_length] += 1

            # Alternative/etc. alignments are ignored for the remaining statistics
            if record.is_supplementary or record.is_secondary:
                continue

            # Mean phred encoded base error probability
            mapq = None if record.is_unmapped else record.mapping_quality
            stats["mapq"][mapq] += 1

            mean_quality_score: int = record.get_tag("qs")  # type: ignore
            stats["qs"][mean_quality_score] += 1

            # Number of uncalled bases in the query
            stats["n"][record.query_sequence.count("N")] += 1

            # matches, mismatches, etc.
            if not record.is_unmapped:
                subsample.add(record)


def collect_detailed_stats(
    sample: str,
    stats: Statistics,
    subsample: ReservoirSample,
) -> None:
    for record in progress(subsample.items, desc=sample):
        query_sequence = record.query_sequence
        assert query_sequence is not None

        # Cigar string elements
        assert record.cigartuples is not None

        cigar_sums: Dict[str, int] = Counter()
        cigar_counts: Dict[Tuple[int, int], int] = Counter(record.cigartuples)
        for (op, length), count in cigar_counts.items():
            # The ops are 'MIDNSHP=X' but we merge match (=)/mismatch (X) into M
            op = "MIDNSHPMM"[op]
            cigar_sums[op] += length * count
            stats["cigar"][op][length] += count

        query_sequence = record.query_sequence
        assert query_sequence is not None, record.query_name
        query_length = len(query_sequence)

        alignment_counts: Dict[Any, int] = Counter()
        for qpos, _, ref in record.get_aligned_pairs(with_seq=True):  # type: ignore
            ref = "-" if ref is None else ref.upper()
            query = "-" if qpos is None else query_sequence[qpos]
            alignment_counts[f"{ref}/{query}"] += 1

        for key, count in alignment_counts.items():
            stats["alignment"][key] += count

        # Percent mapped (excludes indels and clipping, etc.)
        frac_mapped = int(round(cigar_sums["M"] * 100 / query_length))
        stats["pct_mapped"][frac_mapped] += 1

        # Percent mapped bases that are mismatches to reference
        # The NM tag as calculated by `samtools calmd` counts indels as mismatches:
        # https://github.com/samtools/samtools/blob/75a780628e8b7a4e69bd3b2112979f635ed82320/bam_md.c#L92
        n_mismatches: int = record.get_tag("NM") - cigar_sums["D"] - cigar_sums["I"]  # type: ignore
        pct_mismatches = int(round(n_mismatches * 100 / cigar_sums["M"]))
        stats["pct_mapped_mm"][pct_mismatches] += 1


def main_stats(args: argparse.Namespace) -> int:
    samples: Dict[str, List[Path]] = defaultdict(list)
    for filename in args.files:
        sample, _ = Path(filename).stem.split(".", 1)
        samples[sample].append(filename)

    stats: Dict[str, Statistics] = {}
    for sample, filenames in samples.items():
        sample_stats = stats[sample] = new_sample()
        subsample = ReservoirSample(args.nsample)
        for filename in filenames:
            collect_basic_stats(
                sample=sample,
                stats=sample_stats,
                subsample=subsample,
                filename=filename,
                nthreads=args.nthreads,
            )

        collect_detailed_stats(
            sample=sample,
            stats=sample_stats,
            subsample=subsample,
        )

    json.dump(flatten_counts(stats), sys.stdout)

    return 0


def main_join(args: argparse.Namespace) -> int:
    stats = {}
    for filename in args.files:
        with filename.open() as handle:
            data = json.load(handle)
            for key, values in data.items():
                assert key not in stats, key
                stats[key] = values

    json.dump(stats, sys.stdout)

    return 0


class HelpFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self, *args: Any, **kwargs: Any) -> None:
        kwargs.setdefault("width", 79)

        super().__init__(*args, **kwargs)


def parse_args(argv: List[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(formatter_class=HelpFormatter)
    parser.set_defaults(main=lambda _: parser.print_usage())

    subparsers = parser.add_subparsers(title="command")

    stats = subparsers.add_parser("stats")
    stats.set_defaults(main=main_stats)
    stats.add_argument("files", nargs="+", type=Path)
    stats.add_argument("--nthreads", default=4, type=int)
    stats.add_argument("--nsample", default=1_000_000, type=int)

    join = subparsers.add_parser("join")
    join.set_defaults(main=main_join)
    join.add_argument("files", nargs="+", type=Path)

    return parser.parse_args(argv)


def main(argv: List[str]) -> int:
    args = parse_args(argv)

    return args.main(args)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
