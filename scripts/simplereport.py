from __future__ import annotations

import itertools
from typing import Any, List, Optional

import altair as alt
import pandas as pd
from typing_extensions import Literal


class unspecified:
    pass


class cell:
    def __init__(
        self,
        value: None | int | float | str,
        *,
        sort_data: unspecified | int | float | str | None = unspecified(),
        shading: float | Literal["DYNAMIC"] = 0,
    ) -> None:
        self.value = value
        self.shading = (
            max(0, round(shading, 3)) if isinstance(shading, float) else shading
        )

        self.sort_data = sort_data

    @property
    def has_sort_data(self) -> bool:
        return not isinstance(self.sort_data, unspecified)

    def __str__(self) -> str:
        return str(self.value)


class section:
    _title: Optional[str]
    _paragraphs: List[str]
    _running_id: int = 0

    def __init__(self) -> None:
        self._title = None
        self._paragraphs = []

    def set_title(self, title: str) -> "section":
        self._title = title
        return self

    def add_subtitle(self, title: str) -> "section":
        self._paragraphs.append("        <h2>{}</h2>".format(title))

        return self

    def add_chart(
        self,
        chart: alt.Chart | alt.FacetChart,
        caption: Optional[str] = None,
    ) -> "section":
        section._running_id += 1
        self._paragraphs.append(
            _TEMPLATE_ALTAIR.format(
                id=section._running_id,
                spec=chart.to_json(indent=0),  # type: ignore
            )
        )

        if caption is not None:
            self._paragraphs.append(f"      <figcaption>{caption}</figcaption>")

        return self

    def add_table(
        self,
        rows: list[list[str | int | None | float | cell]],
        columns: List[str],
    ) -> "section":
        self._update_shading(rows)
        lines: List[str] = []
        add = lines.append

        add('      <table class="sortable pure-table io-table pure-table-striped">')
        if columns:
            add("        <thead>")
            add("          <tr>")
            for name in columns:
                add(f"            <th>{name}</th>")
            add("          </tr>")
            add("        </thead>")
        add("        <tbody>")
        for rowidx, row in enumerate(rows):
            add("          <tr>")

            for colidx, value in enumerate(row):
                attrs = ""
                if (
                    isinstance(value, cell)
                    and isinstance(value.shading, float)
                    and value.shading > 0
                ):
                    pct = round(100 * value.shading, 1)
                    attrs = f" class='percent' style='background-size: {pct}% 100%'"
                elif value is None:
                    value = ""

                if isinstance(value, cell) and value.has_sort_data:
                    attrs = f'{attrs} data-sort="{value.sort_data}" '

                add(f"            <td{attrs}>{value}</td>")
            add("          </tr>")
        add("        </tbody>")
        add("      </table>")

        self._paragraphs.append("\n".join(lines))

        return self

    def add_dataframe(self, table: pd.DataFrame) -> "section":
        columns: List[Any] = list(table.columns)  # type: ignore
        if any(isinstance(value, int) for value in columns):
            columns = []

        rows: List[Any] = list(list(row) for _, row in table.iterrows())  # type: ignore

        return self.add_table(rows, columns)

    def add_paragraph(self, *text: str) -> "section":
        self._paragraphs.append("      <p>{}</p>".format(" ".join(text)))

        return self

    def render(self) -> str:
        elements: List[str] = ['<div class="section">']
        if self._title is not None:
            elements.append(_TEMPLATE_TITLE.format(title=self._title))

        elements.extend(self._paragraphs)
        elements.append("</div>")

        return "\n".join(elements)

    def _update_shading(
        self,
        rows: list[list[str | int | float | None | cell]],
    ) -> None:
        for column in itertools.zip_longest(*rows):
            values: list[int | float] = []
            for it in column:
                if isinstance(it, cell):
                    if it.shading == "DYNAMIC":
                        if isinstance(it.value, (int, float)):
                            values.append(it.value)
                        else:
                            break
                elif it is not None:
                    break
            else:
                max_value = max(values, default=1)
                for it in column:
                    if isinstance(it, cell) and it.shading == "DYNAMIC":
                        assert isinstance(it.value, (int, float))
                        it.shading = it.value / max_value


class report:
    _title: str
    _sections: List[section]

    def __init__(self, title: str) -> None:
        self._title = title
        self._sections = []

    def add(self) -> section:
        s = section()
        self._sections.append(s)
        return s

    def render(self) -> str:
        return (
            _TEMPLATE_DOC.strip("\r\n")
            .replace("{title}", self._title)
            .replace(
                "{body}",
                "\n\n".join(s.render().rstrip("\r\n") for s in self._sections),
            )
        )


_TEMPLATE_DOC = """
<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <meta http-equiv="X-UA-Compatible" content="ie=edge">
  <title>{title}</title>
  <link
    rel="stylesheet"
    href="https://cdn.jsdelivr.net/npm/purecss@3.0.0/build/pure-min.css"
    integrity="sha384-X38yfunGUhNzHpBaEBsWLO+A0HDYOQi8ufWDkZ0k9e0eXz/tH3II7uKZ9msv++Ls"
    crossorigin="anonymous">
  <script
    src="https://unpkg.com/vega@5.21.0/build/vega.min.js"
    integrity="sha384-s2nYi9D0FfKNopEKsfINeS1Ffhcf+5uvwIrb7Zqso2II+HPhzBTWvXClt+NdUwFc"
    crossorigin="anonymous"></script>
   <script
    src="https://unpkg.com/vega-lite@5.2.0/build/vega-lite.min.js"
    integrity="sha384-tU6fj0fI2gxrcWwC7uBMp70QvipC9ukjcXyOs85VMmdCq33CrA7xQ3nJkJu0SmDm"
    crossorigin="anonymous"></script>
   <script
    src="https://unpkg.com/vega-embed@6.20.2/build/vega-embed.min.js"
    integrity="sha384-oP1rwLY7weRZ5jvAVzfnJsAn+sYA69rQC4geH82Y9oMvr8ruA1oeE9Jkft2noCHR"
    crossorigin="anonymous"></script>
   <link
    rel="stylesheet"
    href="https://cdn.jsdelivr.net/gh/tofsjonas/sortable@3.0.0/sortable-base.min.css"
    integrity="sha384-RMPvgKdhV7JWj5RH7yq5bQgwAt02lpUEAdzaUPj4RLA9BdifNNiI1gtIobFLfeuO"
    crossorigin="anonymous">
   <script
    src="https://cdn.jsdelivr.net/gh/tofsjonas/sortable@latest/sortable.min.js"
    integrity="sha384-Ui7TCZUUp8xuvYhwek30kUmzgl+cbbS0TWhUdE7J3/dF/O7kuF0wd6ZlLZ3KhLVn"
    crossorigin="anonymous"></script>
    <style type='text/css'>
      body {
        background-color: #E3E2DE;
      }

      div#layout {
          max-width: 1280px;
          margin-left: auto;
          margin-right: auto;
          font-size: smaller;
      }

      div.title {
          background-color: #8C9CC0;
          margin: -10px !important;
          text-align: center;
          border-radius: 5px;
      }

      div.title>h1,
      div.title>h2 {
          padding: 5px;
      }

      h5 {
          margin-bottom: 2px;
          margin-left: 6px;
      }

      .pure-table {
          margin-left: 1em;
          text-align: right;
      }

      .pure-table thead>tr>th {
          font-weight: bold;
          background-color: #C4CCDB;
      }

      .pure-table tbody>tr:nth-child(even)>td.percent {
        background-image: linear-gradient(to right, #C4CCDB, #C4CCDB);
        background-repeat: no-repeat;
        background-position: 100% 100%; /* right aligned */
      }

      .pure-table tbody>tr:nth-child(odd)>td.percent {
        background-image: linear-gradient(to right, #D4DCEB, #D4DCEB);
        background-repeat: no-repeat;
        background-position: 100% 100%; /* right aligned */
      }

      .section {
          background-color: #FFF;
          border-radius: 5px;
          margin-bottom: 10px;
          padding: 10px;
          padding-top: 0px;
      }

      .epilogue,
      .note {
          color: #777;
          font-size: small;
          padding-top: 10px;
      }
    </style>
</head>
<body>
  <main>
    <div id='layout'>
      <div class="title">
        <h1>{title}</h1>
      </div>
{body}
    </div>
  </main>
</body>
</html>
"""

_TEMPLATE_TITLE = """
      <div class="title">
        <h2>{title}</h2>
      </div>
"""

_TEMPLATE_ALTAIR = """
      <div id="vega{id}"></div>
      <script type="text/javascript">
          var spec = {spec};
          vegaEmbed("#vega{id}", spec);
      </script>
"""
