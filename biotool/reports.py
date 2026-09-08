"""Accessible HTML companions for the interactive Plotly figures."""

from collections import Counter
from html import escape
from pathlib import Path
from typing import TYPE_CHECKING
import webbrowser

from .theme import PALETTES
from .web_theme import STYLE, THEME_SCRIPT

if TYPE_CHECKING:
    from plotly.graph_objs import Figure
    from .app import ProteinData


def document(title: str, pdb_id: str, protein_title: str, chart: str, body: str,
             *, structure: bool = False, appearance: str = "dark") -> str:
    """Wrap trusted report markup with escaped metadata and semantic navigation."""
    overview_current = '' if structure else ' aria-current="page"'
    structure_current = ' aria-current="page"' if structure else ''
    if appearance not in PALETTES:
        raise ValueError(f"Apariencia desconocida: {appearance}")
    choices = "".join(
        f'<option value="{value}"{" selected" if appearance == value else ""}>{label}</option>'
        for value, label in (("light", "Claro"), ("dark", "Oscuro"))
    )
    return f"""<!doctype html>
<html lang="es" data-theme="{appearance}">
<head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>{escape(pdb_id)} · {escape(title)} | BioTool</title><style>{STYLE}</style></head>
<body><main>
<header><div class="toolbar"><p class="brand">BioTool / {escape(pdb_id)}</p>
<label class="theme-picker" for="appearance">Apariencia
<select id="appearance">{choices}</select></label></div>
<h1>{escape(title)}</h1><p class="muted">{escape(protein_title)}</p></header>
<p id="theme-status" role="status"></p>
<nav aria-label="Informes del análisis">
<a href="simple_plot.html"{overview_current}>Composición y estructura 3D</a>
<a href="basic_pie_chart.html"{structure_current}>Clasificación de tripletes</a>
</nav>
<section class="panel chart" aria-label="Gráfica interactiva">{chart}</section>
{body}
<footer>Fuente: RCSB Protein Data Bank. Residuos observados del primer modelo;
no representan necesariamente la secuencia biológica completa.</footer>
</main>{THEME_SCRIPT}</body></html>"""


def write_reports(protein: "ProteinData", pdb_id: str, figure: "Figure",
                  pie: "Figure", output_dir: Path,
                  *, auto_open: bool, appearance: str = "dark") -> None:
    """Save self-contained reports, then open them only after both writes succeed."""
    config = {"responsive": True, "displaylogo": False, "scrollZoom": False}
    sequence = "".join(protein.chains.values())
    counts = Counter(sequence)
    fasta = []
    for chain, residues in protein.chains.items():
        fasta.append(f">{pdb_id}|cadena {chain.strip() or '(sin identificador)'}")
        fasta.extend(residues[i:i + 80] for i in range(0, len(residues), 80))
    fasta_text = escape("\n".join(fasta))
    chain_label = "cadena" if len(protein.chains) == 1 else "cadenas"
    rows = "".join(
        f"<tr><th scope=\"row\">{escape(aa)}</th><td>{count}</td>"
        f"<td>{100 * count / len(sequence):.1f}%</td></tr>"
        for aa, count in sorted(counts.items())
    )
    overview_body = f"""
<section class="panel"><h2>Resumen del análisis</h2>
<p>{len(sequence)} aminoácidos · {len(protein.chains)} {chain_label} ·
{len(protein.coordinates)} coordenadas atómicas</p>
<p class="muted">La vista 3D muestra átomos, no enlaces. Arrastra para girar;
usa los controles de la gráfica para ampliar o descargar una imagen.</p>
<details><summary>Ver composición en tabla</summary>
<table><caption>Aminoácidos observados</caption>
<thead><tr><th scope="col">Aminoácido</th><th scope="col">Cantidad</th>
<th scope="col">Porcentaje</th></tr></thead><tbody>{rows}</tbody></table></details>
<details open><summary>Secuencia FASTA por cadena</summary>
<pre>{fasta_text}</pre></details></section>"""
    triplet_rows = ""
    if pie.data:
        triplet_rows = "".join(
            f'<tr><th scope="row">{escape(label)}</th><td>{count}</td></tr>'
            for label, count in zip(pie.data[0].labels, pie.data[0].values)
        )
    structure_body = """
<section class="panel"><h2>Cómo interpretar estos resultados</h2>
<p>Clasificación heurística de tripletes, no experimental.
No es una asignación estructural ni una predicción validada.</p>
<p class="muted">Se cuentan tripletes no solapados dentro de cada cadena.
Los residuos finales incompletos se omiten aquí, pero sí se incluyen en la composición.</p>
"""
    if triplet_rows:
        structure_body += f"""<table><caption>Tripletes por clasificación</caption>
<thead><tr><th scope="col">Clasificación</th><th scope="col">Tripletes</th></tr></thead>
<tbody>{triplet_rows}</tbody></table>"""
        structure_body += "<details><summary>Ver tripletes por clasificación</summary>"
        for label, patterns in zip(pie.data[0].labels, pie.data[0].customdata):
            structure_body += (
                f"<h3>{escape(label)}</h3><pre>{escape(patterns) or 'Sin tripletes.'}</pre>"
            )
        structure_body += "</details>"
    else:
        structure_body += "<p>No hay tripletes completos en las cadenas observadas.</p>"
    structure_body += "</section>"
    paths = [output_dir / "simple_plot.html", output_dir / "basic_pie_chart.html"]
    for path, title, chart, body, structure in (
        (paths[0], "Composición y estructura 3D", figure, overview_body, False),
        (paths[1], "Clasificación de tripletes", pie, structure_body, True),
    ):
        path.write_text(document(
            title, pdb_id, protein.title,
            chart.to_html(full_html=False, include_plotlyjs=True, config=config),
            body, structure=structure, appearance=appearance,
        ), encoding="utf-8")
    if auto_open:
        for path in paths:
            webbrowser.open(path.resolve().as_uri())
