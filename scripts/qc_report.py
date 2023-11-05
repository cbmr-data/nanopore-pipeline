#!/usr/bin/env python3
# -*- coding: utf8 -*-
# pyright: basic
import dataclasses
import json
import sys
from collections import Counter
from pathlib import Path
from typing import Any, Iterable, Iterator, Literal

import altair as alt
import pandas as pd
import simplereport
import typed_argparse as tap
from koda_validate import DataclassValidator, Invalid, ListValidator, Validator
from koda_validate.typehints import get_typehint_validator
from qc_metrics import Statistics
from simplereport import cell


class Args(tap.TypedArgs):
    metrics: Path = tap.arg("metrics", positional=True)


def custom_resolver(annotations: Any) -> Validator[Any]:
    if dataclasses.is_dataclass(annotations):
        return DataclassValidator(
            annotations,
            typehint_resolver=custom_resolver,
            fail_on_unknown_keys=True,
        )

    return get_typehint_validator(annotations)


def generate_report(args: Args) -> None:
    try:
        with args.metrics.open() as handle:
            data = json.load(handle)
    except OSError as error:
        sys.exit(f"Error reading metrics file: {error}")

    validator = ListValidator(
        DataclassValidator(
            Statistics,
            typehint_resolver=custom_resolver,
        )
    )

    result = validator(data)
    if isinstance(result, Invalid):
        print(result, file=sys.stderr)
        sys.exit(1)

    samples = sorted(result.val, key=lambda it: it.metadata.name)
    report = simplereport.report("Nanopore quality metrics")

    add_sequencing_table(report, samples)
    add_mapping_table(report, samples)

    plot_query_lengths(report, samples)
    plot_cigar_percentages(report, samples)
    plot_overall_mismatch_rates(report, samples)

    plot_indel_lengths(report, samples)

    print(report.render())


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


class bp_cell(cell):
    def __init__(self, value: int, *, shading: float | Literal["DYNAMIC"] = 0) -> None:
        super().__init__(value, sort_data=value, shading=shading)

    def __str__(self) -> str:
        assert isinstance(self.value, int), self.value
        return f"{self.value / 1e9:.2f}"


class count_cell(cell):
    def __init__(self, value: int, *, shading: float | Literal["DYNAMIC"] = 0) -> None:
        super().__init__(value, sort_data=value, shading=shading)

    def __str__(self) -> str:
        assert isinstance(self.value, int), self.value
        return f"{self.value:,}"


def sum_qlengths(counts: dict[int, int]) -> tuple[int, int]:
    reads = sum(counts.values())
    bp = sum(length * count for length, count in counts.items())
    return reads, bp


def add_sequencing_table(doc: simplereport.report, samples: list[Statistics]) -> None:
    max_query_bp = 1
    for it in samples:
        current_query_bp = 0
        for it in (it.filtered.length, it.filtered.qscore, it.unmapped, it.mapped):
            _, current_bp = sum_qlengths(it.query_lengths)
            current_query_bp += current_bp
        max_query_bp = max(max_query_bp, current_query_bp)

    rows: list[list[int | str | float | None | simplereport.cell]] = []
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

        rows.append(
            [
                sample.metadata.name,
                count_cell(query_reads, shading="DYNAMIC"),
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
                None,
                count_cell(passed_reads),
                count_cell(passed_reads - query_reads),
                bp_cell(
                    passed_bp,
                    shading=passed_bp / max_query_bp,
                ),
                bp_cell(passed_bp - query_bp),
                cell(
                    f"{passed_bp / sample.metadata.genome_size:.1f}",
                    shading=passed_bp / max_query_bp,
                ),
                f"{(passed_bp - query_bp) / sample.metadata.genome_size:.1f}",
                None,
                count_cell(median_length),
                count_cell(median_length - median_length_raw),
            ]
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
            "Sample",
            "Reads<sup>*</sup>",
            "Gbp",
            "×Cov",
            "",
            f"Length < {min_length}",
            f"Q < {min_qscore}",
            "",
            "Filtered",
            "Δ",
            "Gbp",
            "Δ",
            "×Cov",
            "Δ",
            "",
            "Length<sup>†</sup>",
            "Δ",
        ],
    )

    section.add_paragraph(
        "<sup>*</sup> This number should correspond to the number of query sequences "
        "in the input FASTQ(s), but is estimated from the BAM file by excluding "
        "filtered reads."
    )
    section.add_paragraph("<sup>Δ</sup> Change in measure due to filtering.")
    section.add_paragraph(
        "<sup>†</sup> Median length of query sequences after filtering."
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


def add_mapping_table(doc: simplereport.report, samples: list[Statistics]) -> None:
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

    rows: list[list[int | str | float | None | simplereport.cell]] = [
        [
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
        ]
        for sample, it in summary_statistics
    ]

    section = doc.add()
    section.set_title("Mapping statistics")
    section.add_table(
        rows,
        columns=[
            "Sample",
            "Query reads<sup>*</sup>",
            "",
            "Unmapped",
            "%",
            "Gbp",
            "",
            "Mapped",
            "%",
            "Gbp",
            "×Cov",
            "",
            "M%<sup>†</sup>",
            "D%<sup>†</sup>",
            "I%<sup>†</sup>",
            "S%<sup>†</sup>",
            "",
            "MM%<sup>‡</sup>",
        ],
    )

    section.add_paragraph(
        "<sup>*</sup> The number of length/Qscore filtered reads for which mapping was "
        "attempted. "
    )
    section.add_paragraph(
        "<sup>†</sup> The median percentage of alignments consisting of a given "
        "alignment type: [M]atch, [D]eletion, [I]nsertion, and [S]oft Clipped."
    )
    section.add_paragraph(
        "<sup>‡</sup> Average mismatch rate (%) for aligned (M) bases. Default scale "
        "for bars is 0% to 5%"
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
        text=alt.condition(nearest, "length:Q", alt.value(" ")),
    )

    ruler = (
        alt.Chart(data)
        .mark_rule(color="gray")
        .encode(y="quantile:Q")
        .transform_filter(nearest)
    )

    return (
        alt.layer(bars, selectors, points, text, ruler)
        .properties(width=525, height=100)
        .facet(facet="sample:N", columns=2)
        # Repeat x-axis
        .resolve_axis(x="independent")
        # Disable facet header showing column name
        .configure_header(title=None)
        .interactive()
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
            preselected_group="mapped",
        )
    )

    section.add_paragraph(
        "X axis shows the length of query sequences in bp (prior to mapping/clipping)"
    )


########################################################################################


def plot_cigar_percentages(
    doc: simplereport.report,
    samples: list[Statistics],
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
            y=alt.Y("rate:Q").title("Rate"),
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


def main(argv: list[str]) -> None:
    return tap.Parser(Args).bind(generate_report).run(argv)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
