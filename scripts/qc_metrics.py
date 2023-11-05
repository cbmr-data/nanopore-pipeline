#!/usr/bin/env python3
from __future__ import annotations

import dataclasses
import json
import os
import random
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import (
    Annotated,
    Any,
    Iterable,
    NoReturn,
    Optional,
    TypeAlias,
    TypedDict,
    TypeVar,
)

import pysam
import typed_argparse as tap
from koda import Just, Maybe, nothing
from koda_validate import IntValidator, coercer

T = TypeVar("T")


class CustomJSONEncoder(json.JSONEncoder):
    def default(self, o: object) -> object:
        if dataclasses.is_dataclass(o):
            return CustomJSONEncoder._as_dict(o)
        return super().default(o)

    @staticmethod
    def _as_dict(obj: object) -> object:
        if dataclasses.is_dataclass(obj):
            obj = {
                field.name: CustomJSONEncoder._as_dict(getattr(obj, field.name))
                for field in dataclasses.fields(obj)
            }

        return obj


@coercer(int, str)
def allow_coerce_str_to_int(val: Any) -> Maybe[int]:
    if isinstance(val, int):
        return Just(val)
    elif isinstance(val, str) and val.isdigit():
        return Just(int(val))
    else:
        return nothing


IntKey: TypeAlias = Annotated[int, IntValidator(coerce=allow_coerce_str_to_int)]


class SAMFlags:
    paired = 0x0001  # the read is paired in sequencing
    proper_pair = 0x0002  # the read is mapped in a proper pair
    unmapped = 0x0004  # the query sequence itself is unmapped
    mate_unmapped = 0x0008  # the mate is unmapped
    reverse = 0x0010  # strand of the query (1 for reverse)
    mate_reverse = 0x0020  # strand of the mate
    first_in_pair = 0x0040  # the read is the first read in a pair
    second_in_pair = 0x0080  # the read is the second read in a pair
    secondary = 0x0100  # the alignment is not primary
    qc_failed = 0x0200  # QC failure
    duplicate = 0x0400  # optical or PCR duplicate
    supplementary = 0x0800  # supplementary alignment


@dataclass
class MetaData:
    name: str

    sample_size: int
    min_query_length: int
    min_qscore: int
    exclude_flags: int
    base_composition_length: int
    genome_size: int

    filenames: list[str]


@dataclass
class BasicStatistics:
    reads: int = 0
    excluded: int = 0


def _cigar_counts() -> dict[str, dict[IntKey, int]]:
    return defaultdict(Counter)


@dataclass
class UnmappedHistograms:
    query_lengths: dict[IntKey, int] = field(default_factory=Counter)
    qscore_counts: dict[IntKey, int] = field(default_factory=Counter)


@dataclass
class FilteredHistograms:
    length: UnmappedHistograms = field(default_factory=UnmappedHistograms)
    qscore: UnmappedHistograms = field(default_factory=UnmappedHistograms)


@dataclass
class MappedHistograms:
    query_lengths: dict[IntKey, int] = field(default_factory=Counter)
    qscore_counts: dict[IntKey, int] = field(default_factory=Counter)
    mapq_counts: dict[IntKey, int] = field(default_factory=Counter)
    cigar_counts: dict[str, dict[IntKey, int]] = field(default_factory=_cigar_counts)
    mapped_pct: list[int] = field(default_factory=lambda: [0] * 101)
    clipped_pct: list[int] = field(default_factory=lambda: [0] * 101)
    inserted_pct: list[int] = field(default_factory=lambda: [0] * 101)
    deleted_pct: list[int] = field(default_factory=lambda: [0] * 101)
    mismatches: dict[str, int] = field(default_factory=Counter)


class ACGTCounts(TypedDict):
    A: list[int]
    C: list[int]
    G: list[int]
    T: list[int]


def new_acgt_dict(length: int) -> ACGTCounts:
    return ACGTCounts(
        A=[0] * length,
        C=[0] * length,
        G=[0] * length,
        T=[0] * length,
    )


@dataclass
class BaseComposition:
    head: ACGTCounts
    tail: ACGTCounts


@dataclass
class Statistics:
    metadata: MetaData
    compositions: BaseComposition
    basics: BasicStatistics = field(default_factory=BasicStatistics)
    filtered: FilteredHistograms = field(default_factory=FilteredHistograms)
    unmapped: UnmappedHistograms = field(default_factory=UnmappedHistograms)
    mapped: MappedHistograms = field(default_factory=MappedHistograms)


########################################################################################


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
        self.items: list[Any] = []
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


########################################################################################


class StatsArgs(tap.TypedArgs):
    files: list[Path] = tap.arg(
        positional=True,
        help="One or more JSON files generated by the 'stats' command",
    )

    threads: int = tap.arg(
        default=8,
        help="Number of threads for BAM decompression",
    )
    sample_size: int = tap.arg(
        default=1_000_000,
        help="Number of reads sampled for detailed statistics; basic statistics are "
        "collected using all reads",
    )
    min_query_length: int = tap.arg(
        default=500,
        help="Exclude reads with a query-length shorter than this number of base-pairs",
    )
    min_qscore: int = tap.arg(
        default=10,
        help="Minimum average base quality score (QS tag)",
    )
    exclude_flags: int = tap.arg(
        default=SAMFlags.secondary | SAMFlags.supplementary,
        help="Minimum average base quality score (QS tag)",
    )
    base_composition_length: int = tap.arg(
        default=100,
        help="Collect base composition for the first/last N base-pairs",
    )


def main_stats(args: StatsArgs) -> NoReturn:
    samples: dict[str, list[Path]] = defaultdict(list)
    for filename in args.files:
        sample, _ = Path(filename).stem.split(".", 1)
        samples[sample].append(filename)

    sample_stats: list[Statistics] = []
    for sample, filenames in samples.items():
        stats = Statistics(
            metadata=MetaData(
                name=sample,
                sample_size=args.sample_size,
                min_query_length=args.min_query_length,
                min_qscore=args.min_qscore,
                exclude_flags=args.exclude_flags,
                base_composition_length=args.base_composition_length,
                genome_size=-1,
                filenames=[str(it) for it in filenames],
            ),
            compositions=BaseComposition(
                head=new_acgt_dict(args.base_composition_length),
                tail=new_acgt_dict(args.base_composition_length),
            ),
        )

        subsample = ReservoirSample(args.sample_size)
        for filename in filenames:
            collect_basic_stats(
                sample=sample,
                stats=stats,
                subsample=subsample,
                filename=filename,
                args=args,
            )

        collect_detailed_stats(
            sample=sample,
            stats=stats,
            subsample=subsample,
            args=args,
        )

        sample_stats.append(stats)

    json.dump(sample_stats, sys.stdout, cls=CustomJSONEncoder)

    sys.exit(0)


def collect_basic_stats(
    sample: str,
    stats: Statistics,
    subsample: ReservoirSample,
    filename: Path,
    args: StatsArgs,
) -> None:
    with pysam.AlignmentFile(os.fspath(filename), threads=args.threads) as handle:
        if stats.metadata.genome_size == -1:
            stats.metadata.genome_size = sum(handle.lengths)

        for record in progress(handle, desc=f"{sample} [{filename}]"):
            if record.flag & args.exclude_flags:
                stats.basics.excluded += 1
                continue

            stats.basics.reads += 1
            mean_quality_score: int = record.get_tag("qs")  # type: ignore

            if record.query_length < args.min_query_length:
                histograms = stats.filtered.length
            elif mean_quality_score < args.min_qscore:
                histograms = stats.filtered.qscore
            elif record.is_unmapped:
                histograms = stats.unmapped
            else:
                histograms = stats.mapped
                stats.mapped.mapq_counts[record.mapping_quality] += 1

                # matches, mismatches, etc.
                subsample.add(record)

            histograms.query_lengths[record.query_length] += 1
            histograms.qscore_counts[mean_quality_score] += 1


def collect_detailed_stats(
    sample: str,
    stats: Statistics,
    subsample: ReservoirSample,
    args: StatsArgs,
) -> None:
    for record in progress(subsample.items, desc=sample):
        query_sequence = record.query_sequence
        assert query_sequence is not None

        # Cigar string elements
        cigar_tuples: list[tuple[int, int]] | None = record.cigartuples
        assert cigar_tuples is not None

        cigar_sums: dict[str, int] = Counter()
        cigar_counts = Counter(cigar_tuples)
        for (op, length), count in cigar_counts.items():
            # The ops are 'MIDNSHP=X' but we merge match (=) and mismatch (X) into M,
            # and soft (S) and hard (H) clipping into pseudo-op C.
            op = "MIDNCCPMM"[op]
            cigar_sums[op] += length * count
            stats.mapped.cigar_counts[op][length] += count

        def _increment_pct(pct: list[int], key: str, length: int) -> None:
            value = round(cigar_sums[key] * 100 / length)
            assert value < len(pct), (value, key, length, cigar_sums)

            pct[value] += 1

        query_length: int = record.query_length
        # All indels are included, so that the lenght represents the full alignment
        alignment_length: int = sum(cigar_sums[key] for key in "MIDS")

        _increment_pct(stats.mapped.mapped_pct, "M", query_length)
        _increment_pct(stats.mapped.clipped_pct, "C", query_length)
        _increment_pct(stats.mapped.inserted_pct, "I", alignment_length)
        _increment_pct(stats.mapped.deleted_pct, "D", alignment_length)

        if len(query_sequence) >= 2 * args.base_composition_length:
            head = stats.compositions.head
            for idx, nuc in enumerate(query_sequence[: args.base_composition_length]):
                if nuc in "ACGT":
                    head[nuc][idx] += 1

            tail = stats.compositions.tail
            for idx, nuc in enumerate(query_sequence[-args.base_composition_length :]):
                if nuc in "ACGT":
                    tail[nuc][idx] += 1

        mismatches = stats.mapped.mismatches
        for qpos, _, ref in record.get_aligned_pairs(with_seq=True):  # type: ignore
            ref = "-" if ref is None else ref.upper()
            query = "-" if qpos is None else query_sequence[qpos]
            if ref in "ACGT-" and query in "ACGT-":
                mismatches[f"{ref}/{query}"] += 1


########################################################################################


class JoinArgs(tap.TypedArgs):
    files: list[Path] = tap.arg(
        positional=True,
        help="One or more JSON files generated by the 'stats' command",
    )


def main_join(args: JoinArgs) -> NoReturn:
    stats = {}
    for filename in args.files:
        with filename.open() as handle:
            data = json.load(handle)
            for key, values in data.items():
                assert key not in stats, key
                stats[key] = values

    json.dump(stats, sys.stdout)

    sys.exit(0)


def main(argv: list[str]) -> None:
    tap.Parser(
        tap.SubParserGroup(
            tap.SubParser("stats", StatsArgs),
            tap.SubParser("join", JoinArgs),
        ),
    ).bind(main_stats, main_join).run(argv)


if __name__ == "__main__":
    main(sys.argv[1:])
