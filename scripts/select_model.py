#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterator, List, NoReturn, TypeVar, Union

import h5py
import pod5
from typing_extensions import Literal

T = TypeVar("T")


InfoDict = Dict[str, Union[int, float, str, bytes, "InfoDict"]]


def eprint(*args: object) -> None:
    print(*args, file=sys.stderr)


def abort(*args: object) -> NoReturn:
    eprint("ERROR:", *args)
    sys.exit(1)


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
    filetypes: tuple[Literal[".fast5", ".pod5"], ...] = (".fast5", ".pod5"),
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


#######################################################################################


@dataclass
class BatchRates:
    path: Path
    sample_rate: set[int] = field(default_factory=set)
    bases_per_second: set[int] = field(default_factory=set)

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


#######################################################################################


@dataclass
class Model:
    path: str
    sample_rate: int | None
    bases_per_second: int | None


def parse_integer_param(value: str) -> int | None:
    return None if value == "." else int(value)


def read_model_table(filepath: Path) -> list[Model]:
    required_columns = {"sample_rate", "bases_per_second", "model"}
    models: list[Model] = []

    if not filepath.is_file():
        abort(f"model table is not a file: {filepath!r}")

    with filepath.open() as handle:
        header = handle.readline().rstrip("\r\n")
        columns = header.split("\t")

        missing_columns = required_columns.difference(columns)
        if missing_columns:
            abort(f"Missing columns in {filepath}: {missing_columns}")

        for linenum, line in enumerate(handle, start=2):
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != len(columns):
                abort(f"malformed line {linenum} in {filepath}")

            row = dict(zip(columns, fields))
            sample_rate = row["sample_rate"]
            bases_per_second = row["bases_per_second"]
            model = row["model"]

            if sample_rate not in ("5000", "4000", "."):
                abort(
                    f"invalid sample_rate at line {linenum} in {filepath}: {sample_rate!r}"
                )
            elif bases_per_second not in ("400", "260", "."):
                abort(
                    f"invalid bases_per_second at line {linenum} in {filepath}: {bases_per_second!r}"
                )
            elif not (Path(model) / "config.toml").is_file():
                abort(f"model at line {linenum} in {filepath} not a model: {model!r}")

            models.append(
                Model(
                    path=model,
                    sample_rate=parse_integer_param(sample_rate),
                    bases_per_second=parse_integer_param(bases_per_second),
                )
            )

    return models


#######################################################################################


class HelpFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self, *args: Any, **kwargs: Any) -> None:
        kwargs.setdefault("width", 79)

        super().__init__(*args, **kwargs)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        formatter_class=HelpFormatter,
        description="Prints read length and sample frequency paramters from the first "
        "record in fast5 or pod5 files. Results are aggregated by folder.",
    )

    parser.add_argument(
        "models", type=Path, help="table specifying what models to use for what "
    )
    parser.add_argument("root", type=Path, nargs="+")

    return parser


def main(argv: list[str]) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    models = read_model_table(args.models)
    rates: BatchRates | None = None

    for root in args.root:
        for filename in collect_raw_files([root]):
            it = collect_rates(filename)
            if rates is None:
                rates = it
            else:
                rates.sample_rate.update(it.sample_rate)
                rates.bases_per_second.update(it.bases_per_second)

    sample_rate: int | None = None
    bases_per_second: int | None = None

    if rates is None:
        abort("Rates could not be determined")
    elif len(rates.sample_rate) > 1:
        abort(f"Multiple sample-rates: {rates.sample_rate}")
    elif len(rates.bases_per_second) > 1:
        abort(f"Multiple bases-per-second: {rates.bases_per_second}")

    if rates.sample_rate:
        (sample_rate,) = rates.sample_rate

    if rates.bases_per_second:
        (bases_per_second,) = rates.bases_per_second

    candidates: list[Model] = [
        model
        for model in models
        if (
            model.sample_rate == sample_rate
            and model.bases_per_second == bases_per_second
        )
    ]

    if not candidates:
        abort(f"No candidate model for {sample_rate=} and {bases_per_second=}")
    elif len(candidates) > 1:
        abort(f"Multiple candidate models for {sample_rate=} and {bases_per_second=}")

    (candidate,) = candidates
    print(candidate.path)

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
