#!/usr/bin/env python3
# -*- coding: utf8 -*-
# pyright: strict
from __future__ import annotations

import argparse
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, NamedTuple, NoReturn, cast

_REGEX_FOLDER = re.compile(
    r"^\d+_[0-9]+_[\d+][A-Z]_[A-Z]{3}[0-9]+_[0-9a-f]+_barcode[0-9]+$"
)

_REGEX_SAMPLE_NAME = re.compile("^[a-z0-9_]+$", re.I)

_REQUIRED_HEADER = frozenset(
    (
        "ExperimentID",
        "FlowCellID",
        "FlowCell",
        "SampleID",
        "BarcodeID",
    )
)


class RunKey(NamedTuple):
    # Date derived from ExperimentID or filename
    date: str
    # Flowcell ID in the form PAK[0-9]+
    flowcell_id: str
    # Barcode ID of run
    barcode_id: int


@dataclass
class RunMetadata:
    key: RunKey
    # ID/Name of sample
    sample_id: str
    # Complete meta-data for the expected run
    meta: Dict[str, str]

    @classmethod
    def parse(cls, row: Dict[str, str]) -> Iterable[RunMetadata]:
        date, _ = row["ExperimentID"].split("_", 1)
        if not _REGEX_SAMPLE_NAME.match(row["SampleID"]):
            raise ValueError(f"Invalid sample name {row['SampleID']!r}")

        for barcode_id in row["BarcodeID"].split("_"):
            yield RunMetadata(
                key=RunKey(
                    date=date,
                    flowcell_id=row["FlowCellID"],
                    barcode_id=int(barcode_id),
                ),
                sample_id=row["SampleID"],
                meta=dict(row),
            )


@dataclass
class RunFolder:
    key: RunKey
    # Filepath to fast5 files
    path: Path

    @classmethod
    def parse(cls, filepath: Path) -> RunFolder:
        date, _, _, flowcell_id, _, barcode = filepath.name.split("_")
        if not barcode.lower().startswith("barcode"):
            raise ValueError(f"invalid barcode {barcode!r}")

        return RunFolder(
            key=RunKey(
                date=date,
                flowcell_id=flowcell_id,
                barcode_id=int(barcode[7:]),
            ),
            path=filepath,
        )


@dataclass
class Run:
    meta: RunMetadata
    folder: RunFolder

    @property
    def source(self) -> Path:
        return self.folder.path.resolve()

    @property
    def destination(self) -> Path:
        return Path(self.meta.sample_id) / self.folder.path.name


def abort(*args: object, **kwargs: object) -> NoReturn:
    eprint("ERROR:", *args, **kwargs)
    sys.exit(1)


def eprint(*args: object, **kwargs: Any) -> None:
    kwargs.setdefault("file", sys.stderr)
    print(*args, **kwargs)


def read_expected_runs(filepath: Path) -> Dict[RunKey, RunMetadata]:
    result: Dict[RunKey, RunMetadata] = {}
    with filepath.open() as handle:
        header = handle.readline().rstrip("\r\n").split("\t")
        if missing := _REQUIRED_HEADER.difference(header):
            abort(f"missing columns in table header: {missing}")

        for idx, line in enumerate(handle, start=2):
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != len(header):
                abort(f"line {idx} is malformed; wrong number of columns")

            for expected_run in RunMetadata.parse(dict(zip(header, fields))):
                if expected_run.key in result:
                    abort(f"duplicate run {expected_run!r} in {filepath}")

                result[expected_run.key] = expected_run

    eprint(f"found {len(result)} samples in {filepath}")

    return result


def find_run_folders(root: Path) -> Dict[RunKey, RunFolder]:
    queue = [root]
    result: Dict[RunKey, RunFolder] = {}

    while queue:
        current = queue.pop()
        if current.is_dir():
            if _REGEX_FOLDER.match(current.name):
                folder = RunFolder.parse(current)
                if folder.key in result:
                    abort(f"duplicate run {folder!r}")

                result[folder.key] = folder
            else:
                queue.extend(current.iterdir())
        else:
            eprint(f"skipping file {current}")

    eprint(f"found {len(result)} runs in {root}")

    return result


def group_runs(
    meta: Dict[RunKey, RunMetadata],
    folders: Dict[RunKey, RunFolder],
) -> List[Run]:
    any_errors = False
    if missing_folders := meta.keys() - folders:
        eprint("ERROR: missing run folders for")
        for it in sorted(missing_folders):
            eprint(f"    {it}")
        any_errors = True

    if missing_meta := folders.keys() - meta:
        eprint("ERROR: missing meta data for")
        for it in sorted(missing_meta):
            eprint(f"    {it}")
        any_errors = True

    if any_errors:
        abort("cannot proceed")

    result: List[Run] = []
    for key in sorted(meta.keys() | folders):
        result.append(Run(meta=meta[key], folder=folders[key]))

    return result


def create_mapping_folder(root: Path, runs: List[Run]) -> Dict[Path, Path]:
    sources: set[Path] = set()
    result: Dict[Path, Path] = {}
    for it in runs:
        source = it.source
        destination = root / it.destination
        if destination in result:
            abort(f"path is clobbered at {destination}")
        elif source in sources:
            abort(f"source is duplicated: {source}")

        result[destination] = source
        sources.add(source)

    return result


def walk_mapping_folder(root: Path) -> Dict[Path, Path]:
    result: Dict[Path, Path] = {}
    queue: List[Path] = [root]
    while queue:
        it = queue.pop()
        if it.is_symlink():
            result[it] = it.resolve()
        elif it.is_dir():
            queue.extend(it.iterdir())
        else:
            abort(f"unexpected file at {it}")

    return result


def create_link(destination: Path, source: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)

    if destination.is_symlink():
        destination.unlink()
    elif destination.exists():
        abort(f"file exists at {destination}")

    destination.symlink_to(source)


class HelpFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self, *args: Any, **kwargs: Any) -> None:
        kwargs.setdefault("width", 79)

        super().__init__(*args, **kwargs)


class Args(argparse.Namespace):
    name: str
    batches: Path
    commit: bool


def parse_args(argv: List[str]) -> Args:
    parser = argparse.ArgumentParser(formatter_class=HelpFormatter)
    parser.add_argument("name")
    parser.add_argument(
        "--batches",
        default=Path("batches"),
        help="Folder containing batches",
    )
    parser.add_argument(
        "--commit",
        action="store_true",
        help="Apply changes to destination dir",
    )

    return cast(Args, parser.parse_args(argv))


def main(argv: List[str]) -> int:
    args = parse_args(argv)

    table_file = args.batches / (args.name + ".tsv")
    if not table_file.is_file():
        abort("table file does not exist/is not a file:", table_file)

    source_dir = args.batches / (args.name + ".source")
    if not source_dir.is_dir():
        abort("source dir does not exist/is not a dir:", source_dir)

    destination_dir = args.batches / args.name
    if destination_dir.exists() and not destination_dir.is_dir():
        abort("destination dir does exists but is not a dir:", destination_dir)

    runs = group_runs(
        meta=read_expected_runs(table_file),
        folders=find_run_folders(source_dir),
    )

    current_setup = walk_mapping_folder(destination_dir)
    planned_setup = create_mapping_folder(destination_dir, runs)

    for extra in current_setup.keys() - planned_setup:
        eprint("removing extranous symlink", extra)
        if args.commit:
            extra.unlink()

    for destination, source in planned_setup.items():
        current_source = current_setup.get(destination)
        if current_source is None:
            eprint("creating symlink", destination)
            if args.commit:
                create_link(destination=destination, source=source)
        elif current_source != source:
            eprint("replacing symlink", destination)
            eprint("             with", source)
            if args.commit:
                create_link(destination=destination, source=source)

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
