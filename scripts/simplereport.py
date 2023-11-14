from __future__ import annotations

import base64
import itertools
from typing import Any, List, Optional

import altair as alt
import pandas as pd
import vl_convert as vlc
from typing_extensions import Literal, TypeAlias

ImageFormat: TypeAlias = Literal["vega", "png", "jpg", "svg"]


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

    def __init__(
        self,
        image_format: ImageFormat,
        jpeg_quality: int = 90,
    ) -> None:
        self._title = None
        self._paragraphs = []
        self._image_format = image_format
        self._jpeg_quality = jpeg_quality

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
        spec: str = chart.to_json(indent=0)

        if self._image_format == "vega":
            section._running_id += 1
            self._paragraphs.append(
                _TEMPLATE_ALTAIR.format(
                    id=section._running_id,
                    spec=spec,  # type: ignore
                )
            )
        else:
            fmt = self._image_format
            if self._image_format == "svg":
                fmt = "svg+xml"
                data = vlc.vegalite_to_svg(spec).encode()  # type: ignore
            elif self._image_format == "png":
                data = vlc.vegalite_to_png(spec)  # type: ignore
            elif self._image_format == "jpg":
                data = vlc.vegalite_to_jpeg(spec, quality=self._jpeg_quality)  # type: ignore
            else:
                raise NotImplementedError(self._image_format)

            data = base64.b64encode(data).decode()
            self._paragraphs.append(f'<img src="data:image/{fmt};base64, {data}">')

        if caption is not None:
            self._paragraphs.append(f"      <figcaption>{caption}</figcaption>")

        return self

    def add_table(
        self,
        rows: list[list[str | int | None | float | cell]],
        columns: List[str | None],
    ) -> "section":
        self._update_shading(rows)
        lines: List[str] = []
        add = lines.append

        add('      <table class="sortable pure-table io-table pure-table-striped">')
        if columns:
            add("        <thead>")
            add("          <tr>")
            for name in columns:
                attrs = ""
                if name is None:
                    name = ""
                    attrs = ' class="no-sort"'

                add(f"            <th{attrs}>{name}</th>")
            add("          </tr>")
            add("        </thead>")
        add("        <tbody>")
        for row in rows:
            add("          <tr>")

            for value in row:
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
    _image_format: ImageFormat
    _jpeg_quality: int

    def __init__(
        self,
        title: str,
        image_format: ImageFormat = "jpg",
        jpeg_quality: int = 90,
    ) -> None:
        self._title = title
        self._sections = []
        self._image_format = image_format
        self._jpeg_quality = jpeg_quality

    def add(self) -> section:
        s = section(
            image_format=self._image_format,
            jpeg_quality=self._jpeg_quality,
        )
        self._sections.append(s)
        return s

    def render(self) -> str:
        return (
            _TEMPLATE_DOC.strip("\r\n")
            .replace("{title}", self._title)
            .replace("{external}", self._render_external())
            .replace(
                "{body}",
                "\n\n".join(s.render().rstrip("\r\n") for s in self._sections),
            )
        )

    def _render_external(self) -> str:
        external: list[tuple[str, str]] = [
            (
                "https://cdn.jsdelivr.net/npm/purecss@3.0.0/build/pure-min.css",
                "X38yfunGUhNzHpBaEBsWLO+A0HDYOQi8ufWDkZ0k9e0eXz/tH3II7uKZ9msv++Ls",
            ),
            (
                "https://cdn.jsdelivr.net/gh/tofsjonas/sortable@3.0.0/sortable-base.min.css",
                "RMPvgKdhV7JWj5RH7yq5bQgwAt02lpUEAdzaUPj4RLA9BdifNNiI1gtIobFLfeuO",
            ),
            (
                "https://cdn.jsdelivr.net/gh/tofsjonas/sortable@latest/sortable.min.js",
                "Ui7TCZUUp8xuvYhwek30kUmzgl+cbbS0TWhUdE7J3/dF/O7kuF0wd6ZlLZ3KhLVn",
            ),
        ]

        if self._image_format == "vega":
            external += [
                (
                    "https://cdn.jsdelivr.net/npm/vega@5.25.0/build/vega.min.js",
                    "iY3zZAtrtgjJoD8rliThCLEeLUYo8aSNWYQkL+Jaa3KQEAACPnaw/lQIRrFbPCsj",
                ),
                (
                    "https://cdn.jsdelivr.net/npm/vega-lite@5.16.3/build/vega-lite.min.js",
                    "OpYOZH0bO1dKgCRBxral2WygJs8r9nrCCA73wdSb8UQlxpXf362P3+v78uKEKiz2",
                ),
                (
                    "https://cdn.jsdelivr.net/npm/vega-embed@6.22.2/build/vega-embed.min.js",
                    "EA8k5FkiwPXfiSQeH8xlNaljrtD6qj7T49n8VoweOD7Tlm/DHHaoKLDbtJ+8ly5+",
                ),
            ]

        lines: list[str] = []
        for url, sha384 in external:
            if url.endswith(".css"):
                lines.append(_TEMPLATE_CSS.format(url=url, hash=sha384))
            elif url.endswith(".js"):
                lines.append(_TEMPLATE_JS.format(url=url, hash=sha384))
            else:
                raise NotImplementedError(url)

        return "".join(lines)


_TEMPLATE_DOC = """
<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <meta http-equiv="X-UA-Compatible" content="ie=edge">
  <title>{title}</title>
  {external}
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

# sha384 hashes calculated as
#   $ openssl dgst -sha384 -binary ${filename} | openssl base64

_TEMPLATE_CSS = """   <link
    rel="stylesheet"
    href="{url}"
    integrity="sha384-{hash}"
    crossorigin="anonymous">
"""

_TEMPLATE_JS = """   <script
    src="{url}"
    integrity="sha384-{hash}"
    crossorigin="anonymous"></script>
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
