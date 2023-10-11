#!/usr/bin/env python3
# -*- coding: utf8 -*-
# pyright: strict
import argparse
import sys
from pathlib import Path
from typing import Any, List, Tuple, Optional

from multiprocessing import Pool

import h5py
import tqdm


def walk(items: List[Path]) -> List[Path]:
    queue: List[Path] = list(items)
    result: List[Path] = []

    while queue:
        it = queue.pop()
        if it.is_file():
            if it.suffix.lower() == ".fast5":
                result.append(it)
        elif it.is_dir():
            queue.extend(it.iterdir())
        elif not it.exists():
            raise KeyError(it)

    result.sort()
    return result


def test(filename: Path) -> List[Tuple[Path, str]]:
    errors: List[Tuple[Path, str]] = []
    try:
        with h5py.File(filename, "r") as handle:
            for name, item in handle.items():
                if "channel_id" not in item.keys():
                    errors.append((filename, f"{name}: {list(item.keys())[:5]}"))
    except Exception as error:
        errors.append((filename, f"Unhandled exception: {error}"))

    return errors


class HelpFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self, *args: Any, **kwargs: Any) -> None:
        kwargs.setdefault("width", 79)

        super().__init__(*args, **kwargs)


def parse_args(argv: List[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(formatter_class=HelpFormatter)
    parser.add_argument("files", nargs="+", type=Path)
    parser.add_argument("--threads", type=int, default=4)

    return parser.parse_args(argv)


def main(argv: List[str]) -> int:
    args = parse_args(argv)
    filenames = walk(args.files)
    last_filename: Optional[Path] = None

    pool = Pool(args.threads)
    for errors in tqdm.tqdm(pool.imap(test, filenames), total=len(filenames)):
        for filename, error in errors:
            if filename != last_filename:
                print(filename)
                last_filename = filename

            print("  -", error)

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
