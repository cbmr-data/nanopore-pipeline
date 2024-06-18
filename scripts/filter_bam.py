#!/usr/bin/env python3
from __future__ import annotations

import argparse
import functools
import sys
from itertools import groupby
from pathlib import Path

import pysam


def queryname(record: pysam.AlignedSegment) -> str | None:
    return record.query_name


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        formatter_class=functools.partial(
            argparse.ArgumentDefaultsHelpFormatter,
            width=79,
        ),
        allow_abbrev=False,
    )

    parser.add_argument(
        "--input",
        metavar="BAM",
        type=Path,
        required=True,
        help="Name ordered BAM",
    )
    parser.add_argument(
        "--output-pass",
        metavar="BAM",
        type=Path,
        required=True,
        help="Uncompressed BAM in input order",
    )
    parser.add_argument(
        "--output-fail",
        metavar="BAM",
        type=Path,
        required=True,
        help="Compressed BAM in input order",
    )
    parser.add_argument(
        "--min-query-length",
        type=int,
        default=500,
        help="Minimum length of query sequences. This minimum does not apply to "
        "supplementary alignments",
    )
    parser.add_argument(
        "--min-base-quality",
        type=int,
        default=10,
        help="Minimum average base quality based on the 'qs' tag",
    )
    parser.add_argument(
        "--exclude-flags",
        default="0x704",
        help="Flags to exclude from output",
    )

    # FIXME: Triggers deadlock in pysam==0.21.0
    # parser.add_argument("--threads", type=int, default=2)

    return parser.parse_args(argv)


def main(argv: list[str]) -> int:
    args = parse_args(argv)
    if args.exclude_flags.lower().startswith("0x"):
        exclude_flags = int(args.exclude_flags, 16)
    else:
        exclude_flags = int(args.exclude_flags)

    if sys.stdin.isatty():
        print("ERROR: STDIN is terminal", file=sys.stderr)
        return 1
    elif sys.stdout.isatty():
        print("ERROR: STDOUT is terminal", file=sys.stderr)
        return 1

    with (
        pysam.AlignmentFile(args.input) as in_file,
        pysam.AlignmentFile(args.output_pass, "wb", template=in_file) as out_pass,
        pysam.AlignmentFile(args.output_fail, "wb", template=in_file) as out_fail,
    ):
        sorting = in_file.header.to_dict().get("HD", {}).get("SO")
        # Unsorted is assumed to be direct output from mapper, and hence grouped by name
        if sorting not in ("queryname", "unsorted"):
            print(f"ERROR: Input is sorted by {sorting}, not by name", file=sys.stderr)
            return 1

        for _, group in groupby(in_file, key=queryname):
            reads: list[pysam.AlignedSegment] = list(group)
            reads_pass: list[pysam.AlignedSegment] = []
            reads_fail: list[pysam.AlignedSegment] = []

            main_passes = False
            for it in reads:
                passes = not it.flag & exclude_flags

                # Fail all reads if the "main" read fails QC checks
                if passes and not (it.is_secondary or it.is_supplementary):
                    main_passes = passes = passes and (
                        it.get_tag("qs") >= args.min_base_quality
                        and it.query_length >= args.min_query_length
                    )

                if passes:
                    reads_pass.append(it)
                else:
                    reads_fail.append(it)

            if main_passes:
                for it in reads_pass:
                    out_pass.write(it)
            else:
                # Preserve input order also when we fail all reads
                reads_fail = reads

            for it in reads_fail:
                out_fail.write(it)

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
