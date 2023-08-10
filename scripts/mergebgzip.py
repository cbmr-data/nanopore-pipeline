#!/usr/bin/env python3
# -*- coding: utf8 -*-
# pyright: strict
import argparse
import sys
from pathlib import Path
from typing import Any, List, BinaryIO

# BGZIP end-of-file marker
EOF_MARKER = b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42\x43\x02\x00\x1b\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00"


def copy_bgzip(src: BinaryIO, dst: BinaryIO, length: int) -> bool:
    chunk_1 = src.read(length)
    chunk_2 = src.read(length)
    while True:
        chunk_3 = src.read(length)
        if not chunk_3:
            break

        dst.write(chunk_1)
        chunk_1, chunk_2 = chunk_2, chunk_3

    chunk_1 = chunk_1 + chunk_2
    if not chunk_1.endswith(EOF_MARKER):
        return False

    chunk_1 = chunk_1[: -len(EOF_MARKER)]
    dst.write(chunk_1)

    return True


class HelpFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self, *args: Any, **kwargs: Any) -> None:
        kwargs.setdefault("width", 79)

        super().__init__(*args, **kwargs)


def main(argv: List[str]) -> int:
    parser = argparse.ArgumentParser(
        formatter_class=HelpFormatter,
        usage="%(prog)s 1.gz 2.gz .... > out.gz",
    )
    parser.add_argument("files", nargs="+", type=Path)
    parser.add_argument(
        "--buffer",
        type=lambda it: max(len(EOF_MARKER), int(it)),
        default=64 * 1024,
    )

    if sys.stdout.isatty():
        parser.print_usage()

    args = parser.parse_args(argv)

    buffer = sys.stdout.buffer

    for filepath in args.files:
        with filepath.open("rb") as handle:
            if not copy_bgzip(handle, buffer, length=args.buffer):
                print(f"ERROR: Not a BGZIP file: {filepath}", file=sys.stderr)
                return 1

    buffer.write(EOF_MARKER)

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
