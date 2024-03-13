#!/usr/bin/env python3
# pyright: basic
from __future__ import annotations

import dataclasses
import json
import re
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import (
    Any,
    Callable,
    Iterable,
    Iterator,
    Literal,
    NoReturn,
    Optional,
    TypeVar,
)

import altair as alt
import pandas as pd
import simplereport
import typed_argparse as tap
from koda_validate import DataclassValidator, Invalid, ListValidator, Validator
from koda_validate.typehints import get_typehint_validator
from qc_metrics import Statistics
from simplereport import ImageFormat, Row, cell, unspecified

T = TypeVar("T")


class Args(tap.TypedArgs):
    metrics: Path = tap.arg(positional=True)
    output_prefix: Path = tap.arg(positional=True)

    head: Optional[int] = tap.arg(
        help="Process only the first N records; defaults to processing every record",
    )

    image_format: Literal["auto", "vega", "svg", "png", "jpg"] = tap.arg(
        help="Output format for images; if auto 'vega' will be used for 10 or fewer "
        "samples and jpg if there are more samples",
        default="auto",
    )

    jpeg_quality: int = tap.arg(
        default=80,
        help="Quality of JPEG images; a value in the rnage 0 to 100",
    )

    sample_info: Optional[Path] = tap.arg(
        default=None,
        help="Path to TSV table file containing sample information. An 'ExperimentID' "
        "column and a 'SampleID' column is required, all other columns are ignored.",
    )


def eprint(*args: object) -> None:
    print(*args, file=sys.stderr)


def abort(*args: object) -> NoReturn:
    eprint(*args)
    sys.exit(1)


def custom_resolver(annotations: Any) -> Validator[Any]:
    if dataclasses.is_dataclass(annotations):
        return DataclassValidator(
            annotations,
            typehint_resolver=custom_resolver,
            fail_on_unknown_keys=True,
        )

    return get_typehint_validator(annotations)


class timer:
    def __init__(self, desc: str) -> None:
        self._desc = desc
        self._start = 0

    def __enter__(self) -> "timer":
        eprint(self._desc)
        self._start = time.time()
        return self

    def __exit__(self, typ: object, value: object, tb: object) -> None:
        eprint(f" .. finished in {time.time() - self._start:.1f}s")


@dataclasses.dataclass
class Sample:
    name: str
    experiments: list[Experiment]


@dataclasses.dataclass
class Experiment:
    name: str
    timestamp: int | None = None
    samples: list[Sample] = dataclasses.field(default_factory=list)

    @property
    def key(self) -> str:
        if self.timestamp is None:
            return self.name

        return f"{self.timestamp}_{self.name}"

    @staticmethod
    def parse_id(label: str) -> tuple[str, int | None]:
        match = re.match(r"([0-9]{8})_(.*)", label)
        if match is not None:
            timestamp, name = match.groups()
            return name, int(timestamp)

        return label, None


@dataclasses.dataclass
class Metadata:
    samples: dict[str, Sample] = dataclasses.field(default_factory=dict)
    experiments: dict[str, Experiment] = dataclasses.field(default_factory=dict)


def read_sample_information(filepath: Path) -> Metadata:
    samples: dict[str, Sample] = {}
    experiments: dict[str, Experiment] = {}
    with filepath.open() as handle:
        header = handle.readline().rstrip("\r\n").split("\t")
        if {"SampleID", "ExperimentID"} - set(header):
            abort("Required columns SampleID and ExperimentID not found in", filepath)

        for linenum, line in enumerate(handle, start=2):
            values = line.rstrip("\r\n").split("\t")
            if len(values) != len(header):
                abort(f"Malformed line {linenum} in {filepath}")

            row = dict(zip(header, values))
            sample_id = row["SampleID"].strip()
            experiment_id = row["ExperimentID"].strip()
            if not (sample_id and experiment_id):
                abort(f"Missing SampleID or ExperimentID at {filepath}:{linenum}")

            try:
                sample = samples[sample_id]
            except KeyError:
                sample = samples[sample_id] = Sample(name=sample_id, experiments=[])

            try:
                experiment = experiments[experiment_id]
            except KeyError:
                name, timestamp = Experiment.parse_id(experiment_id)
                experiment = experiments[experiment_id] = Experiment(
                    name=name, timestamp=timestamp, samples=[]
                )

            sample.experiments.append(experiment)
            experiment.samples.append(sample)

    return Metadata(samples=samples, experiments=experiments)


def prune_metadata(metadata: Metadata, samples: list[Statistics]) -> None:
    whitelist = set(it.metadata.name for it in samples)
    filtered_samples: dict[str, Sample] = {}
    for key, sample in metadata.samples.items():
        if sample.name not in whitelist:
            for experiment in sample.experiments:
                experiment.samples.remove(sample)
        else:
            filtered_samples[key] = sample

    for sample in samples:
        name = sample.metadata.name
        if name not in filtered_samples:
            filtered_samples[name] = Sample(name=name, experiments=[])

    metadata.samples = filtered_samples
    metadata.experiments = {
        key: value for key, value in metadata.experiments.items() if value.samples
    }


def generate_report(args: Args) -> None:
    metadata = Metadata()
    if args.sample_info is not None:
        with timer("reading sample/experiment data"):
            metadata = read_sample_information(args.sample_info)

    with timer("reading data"):
        try:
            data: list[object] = []
            with args.metrics.open() as handle:
                for line in handle:
                    if args.head is not None and len(data) >= args.head:
                        break

                    data.append(json.loads(line))
        except OSError as error:
            abort(f"Error reading metrics file: {error}")

    with timer("validating data"):
        validator = ListValidator(
            DataclassValidator(
                Statistics,
                typehint_resolver=custom_resolver,
            )
        )

        result = validator(data)

    if isinstance(result, Invalid):
        abort(result)

    samples = sorted(result.val, key=lambda it: it.metadata.name)
    if args.image_format == "auto":
        if len(samples) > 10:
            eprint("WARNING: More than 10 samples; rendering images as static JPGs.")
            eprint("         Use --image-format to choose output format.")
            image_format = "jpg"
        else:
            image_format = "vega"
    else:
        image_format = args.image_format

    prune_metadata(metadata=metadata, samples=samples)

    report = simplereport.report(
        title="Nanopore quality metrics",
        image_format=image_format,
    )

    def timed(
        label: str,
        func: Callable[[simplereport.report, list[Statistics], Metadata], T],
    ) -> T:
        with timer(f"building {label}"):
            return func(report, samples, metadata)

    experiments_table = timed("experiment table", add_experiments_table)
    sequencing_table = timed("sequencing table", add_sequencing_table)
    mapping_table = timed("mapping table", add_mapping_table)

    alt.data_transformers.disable_max_rows()

    with timer("building query length plot"):
        plot_query_lengths(report, samples, metadata, image_format=image_format)

    timed("cigar percentage plot", plot_cigar_percentages)
    timed("mismatch rates plot", plot_overall_mismatch_rates)
    timed("indel length plot", plot_indel_lengths)
    timed("base composition plot", plot_base_composition)

    prefix = args.output_prefix
    report_filepath = prefix.parent / f"{prefix.name}.{image_format}.html"
    experiments_filepath = prefix.parent / f"{prefix.name}.experiments.tsv"
    sequencing_filepath = prefix.parent / f"{prefix.name}.sequencing.tsv"
    mapping_filepath = prefix.parent / f"{prefix.name}.mapping.tsv"

    eprint("Writing HTML report to", report_filepath)
    report_filepath.write_text(report.render())
    eprint("Writing experiments table to", experiments_filepath)
    experiments_table.to_csv(experiments_filepath, sep="\t", index=False)
    eprint("Writing sequencing table to", sequencing_filepath)
    sequencing_table.to_csv(sequencing_filepath, sep="\t", index=False)
    eprint("Writing mapping table to", mapping_filepath)
    mapping_table.to_csv(mapping_filepath, sep="\t", index=False)


########################################################################################


# based on Python's 'inclusive' quantile algorithm
def quantiles_from_counts(
    counts: dict[int, int],
    quantiles: Iterable[float],
) -> Iterator[float]:
    m = sum(counts.values()) - 1
    if m < 1:
        raise ValueError("at least two values required to calculate quantile")

    include_0 = False
    include_1 = False

    points = []
    for value in quantiles:
        if value == 0:
            include_0 = True
            continue
        elif value == 1:
            include_1 = True
            continue
        elif not 0 < value < 1:
            raise ValueError(f"invalid quantile {value!r}")

        n = 1.0 / value
        idx, mod = divmod(m, n)
        points.append((idx, mod, n))

    points.sort(reverse=True)
    if not points:
        raise ValueError("no quantiles specified")

    if include_0:
        yield 0.0

    start = 0
    idx, mod, n = points.pop()
    sorted_counts = sorted(counts.items())
    for nth, (value, count) in enumerate(sorted_counts):
        end = start + count

        while start <= idx < end:
            value_2 = value
            if end <= idx + 1:
                value_2 = sorted_counts[nth + 1][0]

            yield (value * (n - mod) + value_2 * mod) / n
            if not points:
                break

            idx, mod, n = points.pop()

        start = end

    if include_1:
        yield sorted_counts[-1][0]


def mean(counts: list[int]) -> float:
    return sum(idx * count for idx, count in enumerate(counts)) / sum(counts)


########################################################################################


def table_to_pandas(
    rows: list[Row] | list[list[str | int | None | float | cell]],
    columns: list[str | None],
) -> pd.DataFrame:
    data: dict[str, list[object]] = {key: [] for key in columns if key is not None}

    for row in rows:
        for key, value in dict(zip(columns, row)).items():
            if key is not None:
                if isinstance(value, cell):
                    if not isinstance(value.raw_data, unspecified):
                        value = value.raw_data

                data[key].append(value)

    return pd.DataFrame(data, columns=columns, index=None)


########################################################################################


class checkbox(cell):
    def __init__(self, **data: str) -> None:
        super().__init__(True)
        self._data = data

    def __str__(self) -> str:
        checked = "checked " if self.value else ""
        data: list[str] = []
        for key, value in self._data.items():
            data.append(f"data-{key}={value!r}")

        return f'<input type="checkbox" {" ".join(data)} onchange="on_click_checkbox(this)" {checked} />'


class bp_cell(cell):
    def __init__(self, value: int, *, shading: float | Literal["DYNAMIC"] = 0) -> None:
        super().__init__(value, sort_data=value, raw_data=value, shading=shading)

    def __str__(self) -> str:
        assert isinstance(self.value, int), self.value
        return f"{self.value / 1e9:.2f}"


class count_cell(cell):
    def __init__(self, value: int, *, shading: float | Literal["DYNAMIC"] = 0) -> None:
        super().__init__(value, sort_data=value, raw_data=value, shading=shading)

    def __str__(self) -> str:
        assert isinstance(self.value, int), self.value
        return f"{self.value:,}"


def sum_qlengths(counts: dict[int, int]) -> tuple[int, int]:
    reads = sum(counts.values())
    bp = sum(length * count for length, count in counts.items())
    return reads, bp


def add_experiments_table(
    doc: simplereport.report,
    samples: list[Statistics],
    metadata: Metadata,
) -> pd.DataFrame:
    rows: list[Row] = []
    for experiment in sorted(
        metadata.experiments.values(), key=lambda it: (it.timestamp or 0, it.name)
    ):
        repeat_samples = 0
        for sample in experiment.samples:
            if len(sample.experiments) > 1:
                repeat_samples += 1

        rows.append(
            Row(
                checkbox(experiment=experiment.key),
                cell(experiment.timestamp),
                experiment.name,
                len(experiment.samples),
                repeat_samples,
                cls="experiment",
                data={"experiment": experiment.key},
            )
        )

    section = doc.add()
    section.set_title("Experiments")
    section.add_table(
        rows,
        columns=[
            None,
            "Date",
            "Experiment",
            "Samples",
            "Repeat Samples*",
        ],
    )

    section.add_paragraph(
        "<sup>*</sup> The number of samples sequenced in more than one experiment"
    )

    return table_to_pandas(
        rows,
        columns=[
            None,
            "Experiment",
            "Samples",
            "RepeatSamples*",
        ],
    )


def add_sequencing_table(
    doc: simplereport.report,
    samples: list[Statistics],
    metadata: Metadata,
) -> pd.DataFrame:
    max_query_bp = 1
    for it in samples:
        current_query_bp = 0
        for it in (it.filtered.length, it.filtered.qscore, it.unmapped, it.mapped):
            _, current_bp = sum_qlengths(it.query_lengths)
            current_query_bp += current_bp
        max_query_bp = max(max_query_bp, current_query_bp)

    rows: list[Row] = []
    for sample in samples:
        lowq_reads, lowq_bp = sum_qlengths(sample.filtered.qscore.query_lengths)
        short_reads, short_bp = sum_qlengths(sample.filtered.length.query_lengths)
        mapped_reads, mapped_bp = sum_qlengths(sample.mapped.query_lengths)
        unmapped_reads, unmapped_bp = sum_qlengths(sample.unmapped.query_lengths)

        query_reads = sample.basics.reads
        query_bp = lowq_bp + short_bp + mapped_bp + unmapped_bp

        passed_reads = mapped_reads + unmapped_reads
        passed_bp = mapped_bp + unmapped_bp

        query_lengths: dict[int, int] = Counter()
        for it in (sample.unmapped, sample.mapped):
            query_lengths.update(it.query_lengths)

        query_lengths_raw: dict[int, int] = Counter(query_lengths)
        for it in (sample.filtered.length, sample.filtered.qscore):
            query_lengths_raw.update(it.query_lengths)

        median_length = round(next(quantiles_from_counts(query_lengths, [0.5])))
        median_length_raw = round(next(quantiles_from_counts(query_lengths_raw, [0.5])))

        sample_meta = metadata.samples[sample.metadata.name]
        experiment_keys = [experiment.key for experiment in sample_meta.experiments]
        experiments_names = [experiment.name for experiment in sample_meta.experiments]

        rows.append(
            Row(
                checkbox(sample=sample.metadata.name),
                " ".join(sorted(experiments_names)),
                sample.metadata.name,
                count_cell(query_reads, shading="DYNAMIC"),
                count_cell(median_length_raw, shading="DYNAMIC"),
                bp_cell(query_bp, shading="DYNAMIC"),
                f"{query_bp / sample.metadata.genome_size:.1f}",
                None,
                count_cell(
                    short_reads,
                    shading=short_reads / query_reads,
                ),
                count_cell(
                    lowq_reads,
                    shading=lowq_reads / query_reads,
                ),
                count_cell(
                    passed_reads,
                    shading=passed_reads / query_reads,
                ),
                None,
                count_cell(median_length, shading="DYNAMIC"),
                bp_cell(
                    passed_bp,
                    shading=passed_bp / max_query_bp,
                ),
                cell(
                    f"{passed_bp / sample.metadata.genome_size:.1f}",
                    shading=passed_bp / max_query_bp,
                ),
                f"{(passed_bp - query_bp) / sample.metadata.genome_size:.1f}",
                cls="sample",
                data={
                    "sample": sample.metadata.name,
                    "experiments": " ".join(experiment_keys),
                },
            )
        )
    min_lengths = {f"{it.metadata.min_query_length}" for it in samples}
    min_length = next(iter(min_lengths)) if len(min_lengths) == 1 else "N"
    min_qscores = {f"{it.metadata.min_qscore}" for it in samples}
    min_qscore = next(iter(min_qscores)) if len(min_qscores) == 1 else "N"

    section = doc.add()
    section.set_title("Sequencing statistics")
    section.add_table(
        rows,
        columns=[
            None,
            "Experiments",
            "Sample",
            "Reads<sup>*</sup>",
            "Length<sup>†</sup>",
            "Gbp",
            "×",
            None,
            "Short<sup>‡</sup>",
            "LowQ<sup>‡</sup>",
            "Filtered",
            None,
            "Length<sup>†</sup>",
            "Gbp",
            "×",
            "Δ",
        ],
        id="tbl_sequencing",
    )

    notes: list[str] = [
        "<sup>*</sup> This number should correspond to the number of query sequences "
        "in the input FASTQ(s), but is estimated from the BAM file by excluding "
        "supplementary/alternative alignments",
        "<sup>†</sup> Median length of query sequences before/after filtering",
        f"<sup>‡</sup> Quality filtering of short reads (< {min_length} bp) and reads "
        f"with an average base Quality score less than {min_qscore}",
    ]

    if samples:
        notes.append(
            f"× Coverage assuming a {samples[0].metadata.genome_size / 1e9:.2f}Gbp "
            "genome"
        )

    notes.append("Δ Change in measure due to filtering")

    section.add_paragraph("; ".join(notes) + ".")

    return table_to_pandas(
        rows,
        columns=[
            None,
            "Groups",
            "Sample",
            "Reads",
            "Length",
            "bp",
            "X",
            None,
            "ShortReads",
            "LowQReads",
            "FilteredReads",
            None,
            "FilteredLength",
            "Filteredbp",
            "FilteredX",
            None,
        ],
    )


########################################################################################


def sum_cigars(cigars: dict[str, dict[int, int]]) -> dict[str, int]:
    totals: dict[str, int] = Counter()
    for key, counts in cigars.items():
        totals[key] = sum(length * count for length, count in counts.items())
    return totals


def median_of_pct(counts: list[int]) -> cell | str:
    if len(counts) < 2:
        return ""

    median = next(quantiles_from_counts(dict(enumerate(counts)), [0.5]))
    return cell(f"{median:.1f}", shading=median / 100)


def mismatch_rate(counts: dict[str, int]) -> float:
    sites = 0
    diff = 0
    for (ref, _, query), count in counts.items():
        if ref in "ACGT" and query in "ACGT":
            sites += count
            if ref != query:
                diff += count
        elif ref not in "ACGTN-" or query not in "ACGTN-":
            raise AssertionError((ref, query))

    return diff / sites


@dataclasses.dataclass
class mapping_summary:
    mapped_reads: int = 0
    mapped_bp: int = 0
    unmapped_reads: int = 0
    unmapped_bp: int = 0
    mismatch_rate: float = 0

    @property
    def query_reads(self):
        return self.mapped_reads + self.unmapped_reads

    @property
    def query_bp(self):
        return self.mapped_bp + self.unmapped_bp


def add_mapping_table(
    doc: simplereport.report,
    samples: list[Statistics],
    metadata: Metadata,
) -> pd.DataFrame:
    max_query_gbp = 0
    max_mismatch_rate = 0

    summary_statistics: list[tuple[Statistics, mapping_summary]] = []
    for sample in samples:
        it = mapping_summary()
        it.mapped_reads, it.mapped_bp = sum_qlengths(sample.mapped.query_lengths)
        it.unmapped_reads, it.unmapped_bp = sum_qlengths(sample.unmapped.query_lengths)
        it.mismatch_rate = mismatch_rate(sample.mapped.mismatches)

        max_query_gbp = max(max_query_gbp, it.query_bp)
        max_mismatch_rate = max(max_mismatch_rate, it.mismatch_rate)

        summary_statistics.append((sample, it))

    rows: list[Row] = [
        Row(
            sample.metadata.name,
            count_cell(it.query_reads, shading="DYNAMIC"),
            None,
            count_cell(it.unmapped_reads),
            f"{100 * it.unmapped_reads / it.query_reads:.2f}",
            bp_cell(
                it.unmapped_bp,
                shading=it.unmapped_bp / max_query_gbp,
            ),
            None,
            count_cell(it.mapped_reads),
            cell(
                f"{100 * it.mapped_reads / it.query_reads:.2f}",
                shading=it.mapped_reads / it.query_reads,
            ),
            bp_cell(it.mapped_bp),
            cell(
                f"{it.mapped_bp / sample.metadata.genome_size:.1f}",
                shading=it.mapped_bp / max_query_gbp,
            ),
            None,
            median_of_pct(sample.mapped.mapped_pct),
            median_of_pct(sample.mapped.deleted_pct),
            median_of_pct(sample.mapped.inserted_pct),
            median_of_pct(sample.mapped.clipped_pct),
            None,
            cell(
                f"{it.mismatch_rate * 100:.1f}",
                shading=it.mismatch_rate / max(0.05, max_mismatch_rate),
            ),
            cls="mapping",
            data={
                "sample": sample.metadata.name,
                "experiments": " ".join(
                    experiment.key
                    for experiment in metadata.samples[sample.metadata.name].experiments
                ),
            },
        )
        for sample, it in summary_statistics
    ]

    section = doc.add()
    section.set_title("Mapping statistics")
    section.add_table(
        rows,
        columns=[
            "Sample",
            "Query reads<sup>*</sup>",
            None,
            "Unmapped",
            "%",
            "Gbp",
            None,
            "Mapped",
            "%",
            "Gbp",
            "×",
            None,
            "M%<sup>†</sup>",
            "D%<sup>†</sup>",
            "I%<sup>†</sup>",
            "S%<sup>†</sup>",
            None,
            "MM%<sup>‡</sup>",
        ],
        id="tbl_mapping",
    )

    notes: list[str] = [
        "<sup>*</sup> The number of length/Qscore filtered reads for which mapping was "
        "attempted",
    ]

    if samples:
        notes.append(
            f"× Coverage assuming a {samples[0].metadata.genome_size / 1e9:.2f}Gbp "
            "genome"
        )

    notes.append(
        "<sup>†</sup> The median percentage of alignments consisting of a given "
        "alignment type: <u>M</u>atch, <u>D</u>eletion, <u>I</u>nsertion, and "
        "<u>S</u>oft clipped"
    )
    notes.append(
        "<sup>‡</sup> Average mismatch rate (%) for aligned (M) bases. Default scale "
        "for bars is 0% to 5%"
    )

    section.add_paragraph("; ".join(notes) + ".")

    return table_to_pandas(
        rows,
        columns=[
            "Sample",
            "QueryReads",
            None,
            "UnmappedReads",
            None,
            "Unmappedbp",
            None,
            "MappedReads",
            None,
            "Mappedbp",
            "MappedX",
            None,
            "MedianPctM",
            "MedianPctD",
            "MedianPctI",
            "MedianPctS",
            None,
            "MismatchPct",
        ],
    )


########################################################################################


def create_length_plot(
    data: pd.DataFrame,
    groups: list[str],
    start_x_axis_from: int = 1,
    preselected_group: str | None = None,
) -> alt.FacetChart:
    # toggle based on legend selection
    selection = alt.selection_point(
        fields=["group"],
        bind="legend",
        value=preselected_group,
    )

    bars = (
        alt.Chart(data)
        .mark_line(interpolate="step-after")
        .encode(
            x=alt.X("length:Q")
            .title(None)
            .scale(type="log", domainMin=start_x_axis_from),
            y=alt.Y("quantile:Q")
            .title("Quantile over bp")
            # tickCount is required or values has no effect
            .axis(values=[0, 25, 50, 75, 100], tickCount=5)
            .scale(domainMin=0, domainMax=120),
            color=alt.Color(
                "group",
                scale=alt.Scale(domain=groups),
                legend=alt.Legend(
                    title=None,
                    orient="none",
                    legendX=350,
                    legendY=-60,
                    direction="horizontal",
                    titleAnchor="middle",
                ),
            ),
            opacity=alt.condition(selection, alt.value(1), alt.value(0.1)),
        )
        .add_params(
            selection,
        )
    )

    # Select the nearest x-value
    nearest = alt.selection_point(
        nearest=True,
        on="mouseover",
        fields=["quantile"],
        empty=False,
    )

    selectors = (
        alt.Chart(data)
        .mark_point()
        .encode(y="quantile:Q", opacity=alt.value(0))
        .add_params(nearest)
    )

    text = bars.mark_text(align="right", dx=10, dy=-10, fontSize=14).encode(
        text=alt.condition(nearest, "length:Q", alt.value(" ")),
    )

    ruler = (
        alt.Chart(data)
        .mark_rule(color="gray")
        .encode(y="quantile:Q")
        .transform_filter(nearest)
    )

    return (
        alt.layer(bars, selectors, text, ruler)
        .properties(width=525, height=100)
        .facet(facet="sample:N", columns=2)
        # Repeat x-axis
        .resolve_axis(x="independent")
        # Disable facet header showing column name
        .configure_header(title=None)
    )


def length_pct_quantiles(lengths: dict[int, int]) -> list[int]:
    return [
        round(length)
        for length in quantiles_from_counts(
            lengths,
            [x / 100 for x in range(0, 101)],
        )
    ]


def plot_query_lengths(
    doc: simplereport.report,
    samples: list[Statistics],
    metadata: Metadata,
    image_format: ImageFormat,
) -> None:
    data = pd.concat(
        pd.DataFrame(
            {
                "sample": it.metadata.name,
                "group": group,
                "quantile": range(0, 101),
                "length": length_pct_quantiles(counts.query_lengths),
            }
        )
        for it in samples
        for (group, counts) in (
            ("mapped", it.mapped),
            ("unmapped", it.unmapped),
            ("filtered (quality)", it.filtered.qscore),
            ("filtered (length)", it.filtered.length),
        )
        if sum(counts.query_lengths.values()) > 1
    )

    section = doc.add()
    section.set_title("Query lengths")
    section.add_chart(
        create_length_plot(
            data=data,
            groups=[
                "mapped",
                "unmapped",
                "filtered (quality)",
                "filtered (length)",
            ],
            start_x_axis_from=100,
            preselected_group="mapped" if image_format == "vega" else None,
        )
    )

    section.add_paragraph(
        "X axis shows the length of query sequences in bp (prior to mapping/clipping)"
    )


########################################################################################


def plot_cigar_percentages(
    doc: simplereport.report,
    samples: list[Statistics],
    metadata: Metadata,
) -> None:
    data = pd.concat(
        pd.DataFrame(
            {
                "sample": it.metadata.name,
                "group": group,
                "percent": range(0, 101),
                "count": counts,
            }
        )
        for it in samples
        for (group, counts) in (
            ("mapped", it.mapped.mapped_pct),
            ("deleted", it.mapped.deleted_pct),
            ("inserted", it.mapped.inserted_pct),
            ("clipped", it.mapped.clipped_pct),
        )
    )

    # toggle based on legend selection
    selection = alt.selection_point(fields=["group"], bind="legend")

    bars = (
        alt.Chart(data)
        .mark_line(interpolate="step-after")
        .encode(
            x=alt.X("percent:Q").title(None).scale(domainMin=0, domainMax=100),
            y=alt.Y("count:Q").title("Count over %"),
            color=alt.Color(
                "group",
                scale=alt.Scale(
                    domain=[
                        "mapped",
                        "deleted",
                        "inserted",
                        "clipped",
                    ]
                ),
                legend=alt.Legend(
                    title=None,
                    orient="none",
                    legendX=350,
                    legendY=-60,
                    direction="horizontal",
                    titleAnchor="middle",
                ),
            ),
            opacity=alt.condition(selection, alt.value(1), alt.value(0.1)),
        )
        .add_params(
            selection,
        )
    )

    # Select the nearest x-value
    nearest = alt.selection_point(
        nearest=True,
        on="mouseover",
        fields=["percent"],
        empty=False,
    )

    selectors = (
        alt.Chart(data)
        .mark_point()
        .encode(x="percent:Q", opacity=alt.value(0))
        .add_params(nearest)
    )

    # Draw points on the line, and highlight based on selection
    points = bars.mark_point().encode(
        opacity=alt.condition(
            nearest,
            alt.value(1),
            alt.value(0),
        ),
        color=alt.value("white"),
    )

    text = bars.mark_text(align="left", dx=10, dy=-10, fontSize=14).encode(
        text=alt.condition(nearest, "count:Q", alt.value(" ")),
    )

    ruler = (
        alt.Chart(data)
        .mark_rule(color="gray")
        .encode(x="percent:Q")
        .transform_filter(nearest)
    )

    chart = (
        alt.layer(bars, selectors, points, text, ruler)
        .properties(width=525, height=100)
        .facet(facet="sample:N", columns=2)
        # Repeat x-axis
        .resolve_axis(x="independent")
        # Disable facet header showing column name
        .configure_header(title=None)
        .interactive()
    )

    section = doc.add()
    section.set_title("CIGAR Percentages")
    section.add_chart(chart)


########################################################################################


def plot_overall_mismatch_rates(
    doc: simplereport.report,
    samples: list[Statistics],
    metadata: Metadata,
) -> None:
    data: list[list[str | int | float]] = []
    for it in samples:
        totals: dict[str, int] = Counter()
        for (ref, _, query), count in it.mapped.mismatches.items():
            if ref in "ACGT" and query in "ACGT":
                totals[ref] += count
            elif ref not in "ACGTN-" or query not in "ACGTN-":
                raise AssertionError((ref, query))

        for (ref, _, query), count in it.mapped.mismatches.items():
            if ref in "ACGT" and query in "ACGT" and ref != query:
                data.append(
                    [
                        it.metadata.name,
                        ref,
                        query,
                        f"{ref}to{query}",
                        count / totals[ref],
                    ]
                )

    columns = ["sample", "reference", "query", "event", "rate"]
    chart = (
        alt.Chart(pd.DataFrame(data, columns=columns))
        .mark_bar(width=15)
        .encode(
            x=alt.X("event:N", axis=alt.Axis(labelAngle=-45)).title(None),
            y=alt.Y("rate:Q").title("Rate").scale(domainMax=0.01, clamp=True),
            color=alt.Color(
                "event:N",
                legend=None,
            ).scale(scheme="paired"),
            xOffset="reference:N",
            tooltip=columns,
        )
        .properties(
            width=250,
            height=100,
        )
        .facet(
            facet="sample:N",
            columns=4,
        )
        # Disable facet header showing column name
        .configure_header(title=None)
        # Repeat x-axis
        .resolve_axis(x="independent")
    )
    section = doc.add()
    section.set_title("Mismatch rates")
    section.add_chart(chart)


########################################################################################


def plot_indel_lengths(
    doc: simplereport.report,
    samples: list[Statistics],
    metadata: Metadata,
) -> None:
    data = pd.concat(
        pd.DataFrame(
            {
                "sample": it.metadata.name,
                "group": group,
                "quantile": range(0, 101),
                "length": length_pct_quantiles(it.mapped.cigar_counts[key]),
            }
        )
        for it in samples
        for (key, group) in (
            ("I", "Insertions"),
            ("D", "Deletions"),
            ("S", "Soft Clipped"),
        )
    )

    section = doc.add()
    section.set_title("CIGAR lengths")
    section.add_chart(
        create_length_plot(
            data=data,
            groups=[
                "Insertions",
                "Deletions",
                "Soft Clipped",
            ],
        )
    )


########################################################################################


def plot_base_composition(
    doc: simplereport.report,
    samples: list[Statistics],
    metadata: Metadata,
) -> None:
    dataframes: list[pd.DataFrame] = []
    for it in samples:
        for key, counts in (
            ("head", it.compositions.head),
            ("tail", it.compositions.tail),
        ):
            totals = [
                sum(row)
                for row in zip(counts["A"], counts["C"], counts["G"], counts["T"])
            ]
            for group in ("A", "C", "G", "T"):
                data = pd.DataFrame(
                    {
                        "sample": f"{it.metadata.name} ({key})",
                        "group": group,
                        "position": range(it.metadata.base_composition_length),
                        "frequency": [
                            round(value / totals[idx], 3)
                            for idx, value in enumerate(counts[group])
                        ],
                    }
                )

                dataframes.append(data)

    data = pd.concat(dataframes)
    # toggle based on legend selection
    selection = alt.selection_point(fields=["group"], bind="legend")

    bars = (
        alt.Chart(data)
        .mark_line(interpolate="step-after")
        .encode(
            x=alt.X("position:Q").title(None),
            y=alt.Y("frequency:Q"),
            color=alt.Color(
                "group",
                scale=alt.Scale(domain=list("ACGT")),
                legend=alt.Legend(
                    title=None,
                    orient="none",
                    legendX=350,
                    legendY=-60,
                    direction="horizontal",
                    titleAnchor="middle",
                ),
            ),
            opacity=alt.condition(selection, alt.value(1), alt.value(0.1)),
        )
        .add_params(
            selection,
        )
    )

    # Select the nearest x-value
    nearest = alt.selection_point(
        nearest=True,
        on="mouseover",
        fields=["position"],
        empty=False,
    )

    selectors = (
        alt.Chart(data)
        .mark_point()
        .encode(x="position:Q", opacity=alt.value(0))
        .add_params(nearest)
    )

    text = bars.mark_text(align="left", dx=10, dy=-10, fontSize=14).encode(
        text=alt.condition(nearest, "frequency:Q", alt.value(" ")),
    )

    ruler = (
        alt.Chart(data)
        .mark_rule(color="gray")
        .encode(x="position:Q")
        .transform_filter(nearest)
    )

    chart = (
        alt.layer(bars, selectors, text, ruler)
        .properties(width=550, height=100)
        .facet(facet="sample:N", columns=2)
        # Repeat x-axis
        .resolve_axis(x="independent")
        # Disable facet header showing column name
        .configure_header(title=None)
    )

    section = doc.add()
    section.set_title("Base compositon")
    section.add_chart(chart)


########################################################################################


def main(argv: list[str]) -> None:
    return tap.Parser(Args).bind(generate_report).run(argv)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
