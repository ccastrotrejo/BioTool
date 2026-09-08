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
        raise ValueError(f"Unknown appearance: {appearance}")
    choices = "".join(
        f'<option value="{value}"{" selected" if appearance == value else ""}>{label}</option>'
        for value, label in (("light", "Light"), ("dark", "Dark"))
    )
    return f"""<!doctype html>
<html lang="en" data-theme="{appearance}">
<head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>{escape(pdb_id)} · {escape(title)} | BioTool</title><style>{STYLE}</style></head>
<body><main>
<header><div class="toolbar"><p class="brand">BioTool / {escape(pdb_id)}</p>
<label class="theme-picker" for="appearance">Appearance
<select id="appearance">{choices}</select></label></div>
<h1>{escape(title)}</h1><p class="muted">{escape(protein_title)}</p></header>
<p id="theme-status" role="status"></p>
<nav aria-label="Analysis reports">
<a href="simple_plot.html"{overview_current}>Composition and 3D structure</a>
<a href="basic_pie_chart.html"{structure_current}>Triplet classification</a>
</nav>
<section class="panel chart" aria-label="Interactive chart">{chart}</section>
{body}
<footer>Source: RCSB Protein Data Bank. Observed residues from the first model;
these do not necessarily represent the complete biological sequence.</footer>
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
        fasta.append(f">{pdb_id}|chain {chain.strip() or '(no identifier)'}")
        fasta.extend(residues[i:i + 80] for i in range(0, len(residues), 80))
    fasta_text = escape("\n".join(fasta))
    chain_label = "chain" if len(protein.chains) == 1 else "chains"
    rows = "".join(
        f"<tr><th scope=\"row\">{escape(aa)}</th><td>{count}</td>"
        f"<td>{100 * count / len(sequence):.1f}%</td></tr>"
        for aa, count in sorted(counts.items())
    )
    overview_body = f"""
<section class="panel"><h2>Analysis summary</h2>
<p>{len(sequence)} amino acids · {len(protein.chains)} {chain_label} ·
{len(protein.coordinates)} atomic coordinates</p>
<p class="muted">The 3D view shows atoms, not bonds. Drag to rotate;
use the chart controls to zoom or download an image.</p>
<details><summary>View composition table</summary>
<table><caption>Observed amino acids</caption>
<thead><tr><th scope="col">Amino acid</th><th scope="col">Count</th>
<th scope="col">Percentage</th></tr></thead><tbody>{rows}</tbody></table></details>
<details open><summary>FASTA sequence by chain</summary>
<pre>{fasta_text}</pre></details></section>"""
    triplet_rows = ""
    if pie.data:
        triplet_rows = "".join(
            f'<tr><th scope="row">{escape(label)}</th><td>{count}</td></tr>'
            for label, count in zip(pie.data[0].labels, pie.data[0].values)
        )
    structure_body = """
<section class="panel"><h2>How to interpret these results</h2>
<p>Heuristic triplet classification, not experimental.
This is neither a structural assignment nor a validated prediction.</p>
<p class="muted">Non-overlapping triplets are counted within each chain.
Trailing residues that do not form a complete triplet are omitted here,
but are included in the composition.</p>
"""
    if triplet_rows:
        structure_body += f"""<table><caption>Triplets by classification</caption>
<thead><tr><th scope="col">Classification</th><th scope="col">Triplets</th></tr></thead>
<tbody>{triplet_rows}</tbody></table>"""
        structure_body += "<details><summary>View triplets by classification</summary>"
        for label, patterns in zip(pie.data[0].labels, pie.data[0].customdata):
            structure_body += (
                f"<h3>{escape(label)}</h3><pre>{escape(patterns) or 'No triplets.'}</pre>"
            )
        structure_body += "</details>"
    else:
        structure_body += "<p>No complete triplets in the observed chains.</p>"
    structure_body += "</section>"
    paths = [output_dir / "simple_plot.html", output_dir / "basic_pie_chart.html"]
    for path, title, chart, body, structure in (
        (paths[0], "Composition and 3D structure", figure, overview_body, False),
        (paths[1], "Triplet classification", pie, structure_body, True),
    ):
        path.write_text(document(
            title, pdb_id, protein.title,
            chart.to_html(full_html=False, include_plotlyjs=True, config=config),
            body, structure=structure, appearance=appearance,
        ), encoding="utf-8")
    if auto_open:
        for path in paths:
            webbrowser.open(path.resolve().as_uri())
