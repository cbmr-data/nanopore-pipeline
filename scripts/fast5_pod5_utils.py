#!/usr/bin/env python3
# -*- coding: utf8 -*-
# pyright: strict
import argparse
import sys
from dataclasses import dataclass, field
from multiprocessing import Pool
from pathlib import Path
from typing import (
    Any,
    Callable,
    Dict,
    Iterable,
    Iterator,
    List,
    Optional,
    Set,
    Tuple,
    TypeVar,
    Union,
)

import h5py
import pod5
import tqdm
from typing_extensions import Literal

T = TypeVar("T")


InfoDict = Dict[str, Union[int, float, str, bytes, "InfoDict"]]


def read_raw_info(filename: Path) -> InfoDict:
    config: InfoDict = {}
    suffix = filename.suffix.lower()

    if suffix == ".fast5":
        with h5py.File(filename, "r") as handle:
            for name, item in handle.items():
                for group in ("channel_id", "context_tags"):
                    config[group] = dict(item[group].attrs.items())
                break
    elif suffix == ".pod5":
        with pod5.Reader(filename) as handle:
            for record in handle:
                run_info = record.run_info
                for key in dir(run_info):
                    if not key.startswith("_"):
                        config[key] = getattr(run_info, key)

                break

    return config


def collect_raw_files(
    items: List[Path],
    filetypes: Tuple[Literal[".fast5", ".pod5"], ...] = (".fast5", ".pod5"),
) -> Iterator[Path]:
    queue: List[Path] = list(items)
    while queue:
        it = queue.pop()
        if it.is_file():
            if it.suffix.lower() in filetypes:
                yield it
        elif it.is_dir():
            queue.extend(it.iterdir())
        elif not it.exists():
            raise KeyError(it)


def process(
    func: Callable[[Path], T],
    filenames: Iterable[Path],
    *,
    threads: int = 1,
    quiet: bool = False,
) -> Iterator[T]:
    pool = Pool(threads)
    imap = pool.imap(func, filenames)

    if not quiet:
        try:
            total = len(filenames)
        except TypeError:
            total = None

        yield from tqdm.tqdm(imap, total=total)
    else:
        yield from imap


#######################################################################################


def validate_fast5(filename: Path) -> List[Tuple[Path, str]]:
    errors: List[Tuple[Path, str]] = []
    try:
        with h5py.File(filename, "r") as handle:
            for name, item in handle.items():
                if "channel_id" not in item.keys():
                    errors.append((filename, f"{name}: {list(item.keys())[:5]}"))
                    break
    except Exception as error:
        errors.append((filename, f"Unhandled exception: {error}"))

    return errors


def main_validate_fast5(args: argparse.Namespace) -> int:
    filenames = collect_raw_files(args.roots, filetypes=(".fast5",))
    last_filename: Optional[Path] = None
    for errors in process(validate_fast5, filenames, threads=args.threads):
        for filename, error in errors:
            if filename != last_filename:
                print(filename)
                last_filename = filename

            print("  -", error)

    return 0


#######################################################################################


def worker_metadata(filename: Path) -> Tuple[Path, InfoDict]:
    return filename, read_raw_info(filename)


def padding_metadata(info: InfoDict, step: int = 2) -> int:
    padding = 0
    for key, value in info.items():
        if len(key) + step > padding:
            padding = len(key) + step

        if isinstance(value, dict):
            padding = max(padding, 2 + padding_metadata(value, step))

    return padding


def print_metadata(info: InfoDict, padding: int, indent: int = 2) -> None:
    indent_s = " " * indent
    for key, value in sorted(info.items()):
        if isinstance(value, dict):
            print(f"{indent_s}{key}/")
            print_metadata(value, padding - 2, indent=indent + 2)
        else:
            if isinstance(value, bytes):
                value = value.decode("utf-8")

            print(f"{indent_s}{key.ljust(padding)} = {value}")


def main_metadata(args: argparse.Namespace) -> int:
    filenames = collect_raw_files(args.roots)
    for filepath, info in process(
        worker_metadata,
        filenames,
        threads=args.threads,
        quiet=True,
    ):
        print(f"{filepath}/")
        print_metadata(info, padding=padding_metadata(info))

    return 0


#######################################################################################


@dataclass
class BatchRates:
    path: Path
    sample_rate: Set[int] = field(default_factory=set)
    bases_per_second: Set[int] = field(default_factory=set)

    def update(self, info: InfoDict) -> None:
        for key in ("sample_rate", "sample_frequency", "sampling_rate"):
            value = info.get(key)
            if value is not None:
                assert isinstance(value, (str, bytes, float, int)), info
                self.sample_rate.add(int(value))

        if "selected_speed_bases_per_second" in info:
            value = info["selected_speed_bases_per_second"]
            assert isinstance(value, (str, bytes, int, float)), info
            self.bases_per_second.add(int(value))

        for child in info.values():
            if isinstance(child, dict):
                self.update(child)


def collect_rates(it: Path) -> BatchRates:
    rates = BatchRates(path=it.parent)
    rates.update(read_raw_info(it))

    return rates


def main_params(args: argparse.Namespace) -> int:
    def _fmt(it: Set[Any]) -> str:
        return "/".join(map(str, it)) or "NA"

    results: Dict[Path, BatchRates] = {}
    filenames = collect_raw_files(args.roots)
    for it in process(collect_rates, filenames, threads=args.threads):
        try:
            dest = results[it.path]
        except KeyError:
            results[it.path] = it
            continue

        dest.sample_rate.update(it.sample_rate)
        dest.bases_per_second.update(it.bases_per_second)

    print(
        "path",
        "sample_rate",
        "bases_per_second",
        sep="\t",
    )

    for it in results.values():
        print(
            it.path,
            _fmt(it.sample_rate),
            _fmt(it.bases_per_second),
            sep="\t",
        )

    return 0


#######################################################################################


class HelpFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self, *args: Any, **kwargs: Any) -> None:
        kwargs.setdefault("width", 79)

        super().__init__(*args, **kwargs)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        formatter_class=HelpFormatter,
        description="Utility for performing (read) operations on fast5 and pod5 files",
    )
    parser.set_defaults(main=None)

    default_threads = 8

    subparsers = parser.add_subparsers()
    subparser = subparsers.add_parser(
        "validate",
        description="Checks that all records have a `channel_id` entry, the absence "
        "which was observed for some malformed files, causing dorado to crash.",
    )
    subparser.set_defaults(main=main_validate_fast5)
    subparser.add_argument("roots", nargs="+", type=Path)
    subparser.add_argument("--threads", type=int, default=default_threads)

    subparser = subparsers.add_parser(
        "metadata",
        description="Prints metadata from the first record in fast5 or pod5 files",
    )
    subparser.set_defaults(main=main_metadata)
    subparser.add_argument("roots", nargs="+", type=Path)
    subparser.add_argument("--threads", type=int, default=default_threads)

    subparser = subparsers.add_parser(
        "params",
        description="Prints read length and sample frequency paramters from the first "
        "record in fast5 or pod5 files. Results are aggregated by folder.",
    )
    subparser.set_defaults(main=main_params)
    subparser.add_argument("roots", nargs="+", type=Path)
    subparser.add_argument("--threads", type=int, default=default_threads)

    return parser


def main(argv: List[str]) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.main is None:
        parser.print_usage()
        return 1

    return args.main(args)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
