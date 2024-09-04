#!/usr/bin/env python3
# ruff: noqa: UP006, UP035, FA100
import argparse
import getpass
import os
import shlex
import subprocess
import sys
from pathlib import Path
from typing import Any, List, NoReturn, TextIO, Union


def eprint(*args: Any, file: TextIO = sys.stderr, **kwargs: Any):
    print(*args, file=file, **kwargs)


def abort(*args: Any, **kwargs: Any) -> NoReturn:
    eprint("ERROR:", *args, **kwargs)
    sys.exit(1)


def quote(value: Union[str, Path]) -> str:
    return shlex.quote(str(value))


class HelpFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self, *args: Any, **kwargs: Any) -> None:
        kwargs.setdefault("width", 79)

        super().__init__(*args, **kwargs)


def parse_args(argv: List[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(formatter_class=HelpFormatter)
    parser.add_argument("project")
    parser.add_argument("args", nargs="*")
    parser.add_argument(
        "--projects",
        type=Path,
        default=Path("projects"),
    )
    parser.add_argument(
        "--profile",
        type=Path,
        default=Path("profile"),
    )
    parser.add_argument(
        "--config",
        action="append",
        default=[],
    )
    parser.add_argument(
        "--shadow-prefix",
        default=Path("/scratch") / getpass.getuser(),
        help="Shadow prefix for temporary files used for certain tasks to improve "
        "performance",
    )
    parser.add_argument(
        "--bin",
        default=Path("./venv/bin").resolve(),
        help="Bin folder added (prefixed) to PATH",
    )

    return parser.parse_args(argv)


def main(argv: List[str]) -> int:
    args = parse_args(argv)

    project: Path = args.projects / args.project
    data: Path = project / "data"
    results: Path = project / "results"
    models: Path = project / "models.tsv"

    profile: Path = args.profile

    if not os.path.lexists(project):
        abort(f"no project directory at {quote(project)}")
    elif not project.is_dir():
        abort(f"invalid project directory at {quote(project)}")
    elif not os.path.lexists(data):
        abort(f"no data directory at {quote(data)}")
    elif not data.is_dir():
        abort(f"invalid data directory at {quote(data)}")
    elif not models.is_file():
        abort(f"models.tsv table not found at {quote(models)}")
    elif not args.bin.is_dir():
        abort(f"--bin is not a directory: {quote(args.bin)}")

    if not profile.is_dir():
        abort(f"invalid profile directory at {quote(profile)}")

    command: List[Union[str, Path]] = [
        "snakemake",
        # excluding defaults 'code' and 'software-env'
        "--rerun-triggers",
        "mtime",
        "params",
        "input",
        "--profile",
        profile,
        "--shadow-prefix",
        args.shadow_prefix,
    ]

    configfile = project / "config.yaml"
    if configfile.exists():
        command += ["--configfile", configfile]
    else:
        eprint(f"WARNING: Using default config: {str(configfile)!r} not found")

    # Set *after* the configfile; not clear if this is required
    command += [
        "--config",
        f"input_dir={data}",
        f"results_dir={results}",
        f"dorado_models={models}",
    ]

    # Extra config values are guaranteed to be added before other args
    command.extend(args.config)
    # Additional command-line arguments; any format
    command.extend(args.args)

    env = dict(os.environ)
    env["PATH"] = "{}:{}".format(args.bin, env["PATH"])

    eprint("Command =", " ".join(quote(it) for it in command))
    return subprocess.call(
        command,
        stdin=subprocess.DEVNULL,
        env=env,
    )


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
