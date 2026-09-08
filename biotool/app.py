#!/usr/bin/env python3
"""Analyze observed protein residues from the first model of an RCSB PDB file.

Install with ``pip install .`` and launch with ``biotool`` or ``python gui.py``.
The desktop interface requires a Python installation with Tk support.
"""

from collections import Counter
from dataclasses import dataclass
from math import isfinite
from pathlib import Path
from typing import Optional
import re
import urllib.error
import urllib.request

from .theme import CHART_FONT, GROUP_PALETTES, PALETTES


Estructura = {
    "A": [1.25, 0.89, 0.78], "R": [0.99, 1.02, 0.88],
    "N": [0.87, 0.86, 1.28], "D": [1.03, 0.74, 1.41],
    "C": [1.12, 0.85, 0.8], "E": [1.45, 0.65, 1],
    "Q": [1.24, 0.82, 0.97], "G": [0.57, 0.93, 1.64],
    "H": [1.25, 1.04, 0.69], "I": [0.94, 1.41, 0.51],
    "L": [1.32, 1.03, 0.59], "K": [1.24, 0.81, 0.96],
    "M": [1.43, 0.99, 0.39], "F": [1.08, 1.22, 0.58],
    "P": [0.6, 0.71, 1.91], "S": [0.82, 0.96, 1.33],
    "T": [0.81, 1.13, 1.03], "W": [1.03, 1.15, 0.75],
    "Y": [0.75, 1.25, 1.05], "V": [0.88, 1.48, 0.47],
}

AA3_TO_1 = {
    "ALA": "A", "VAL": "V", "PHE": "F", "PRO": "P", "MET": "M",
    "ILE": "I", "LEU": "L", "ASP": "D", "GLU": "E", "LYS": "K",
    "ARG": "R", "SER": "S", "THR": "T", "TYR": "Y", "HIS": "H",
    "CYS": "C", "ASN": "N", "GLN": "Q", "TRP": "W", "GLY": "G",
    "MSE": "M",
}

AMINO_ACID_GROUPS = (
    ("Group A: Nonpolar Amino Acids (Hydrophobic)", (
        ("G", "GLY", "Glycine", "Gly"), ("A", "ALA", "Alanine", "Ala"),
        ("V", "VAL", "Valine", "Val"), ("L", "LEU", "Leucine", "Leu"),
        ("I", "ILE", "Isoleucine", "Ile"), ("M", "MET", "Methionine", "Met"),
        ("F", "PHE", "Phenylalanine", "Phe"),
        ("W", "TRP", "Tryptophan", "Trp"), ("P", "PRO", "Proline", "Pro"),
    )),
    ("Group B: Polar, Uncharged Amino Acids (Hydrophilic)", (
        ("S", "SER", "Serine", "Ser"), ("T", "THR", "Threonine", "Thr"),
        ("C", "CYS", "Cysteine", "Cysteine"), ("Y", "TYR", "Tyrosine", "Tyr"),
        ("N", "ASN", "Asparagine", "Asn"), ("Q", "GLN", "Glutamine", "Gln"),
    )),
    ("Group C: Polar, Negatively Charged Amino Acids (Hydrophilic)", (
        ("D", "ASP", "Aspartic Acid", "Asp"), ("E", "GLU", "Glutamic Acid", "Glu"),
    )),
    ("Group D: Polar, Positively Charged Amino Acids (Hydrophilic)", (
        ("K", "LYS", "Lysine", "Lys"), ("R", "ARG", "Arginine", "Arg"),
        ("H", "HIS", "Histidine", "His"),
    )),
)

@dataclass
class ProteinData:
    title: str
    chains: dict[str, str]
    coordinates: list[tuple[float, float, float]]


def probAA(text: str) -> dict[str, list[float]]:
    """Score complete non-overlapping triplets using the original propensities."""
    invalid = set(text) - Estructura.keys()
    if invalid:
        raise ValueError(f"Aminoácidos no reconocidos: {', '.join(sorted(invalid))}")
    return {
        text[i:i + 3]: [
            round(sum(Estructura[aa][axis] for aa in text[i:i + 3]) / 3, 4)
            for axis in range(3)
        ]
        for i in range(0, len(text) - 2, 3)
    }


def estrAA(text: str) -> dict[str, str]:
    """Classify triplets with the legacy heuristic, not a structural assignment."""
    structures = {}
    for pattern, (alpha, beta, turn) in probAA(text).items():
        if alpha > 1.1 and beta < 1.2 and turn < 1.3:
            structures[pattern] = "alfa"
        elif alpha < 1.1 and beta > 1 and turn < 1.3:
            structures[pattern] = "beta"
        elif alpha < 1.25 and beta < 1 and turn > 1.15:
            structures[pattern] = "giro beta"
        else:
            structures[pattern] = "azar"
    return structures


def countEstr(text: str) -> list:
    """Return each group's triplets and count, retaining repeated occurrences."""
    groups = {"alfa": [], "beta": [], "giro beta": [], "azar": []}
    structures = estrAA(text)
    for i in range(0, len(text) - 2, 3):
        pattern = text[i:i + 3]
        groups[structures[pattern]].append(pattern)
    return [item for patterns in groups.values() for item in (patterns, len(patterns))]


def normalize_pdb_id(value: str) -> str:
    """Validate a classic four-character PDB identifier before network or file I/O."""
    pdb_id = value.strip().upper()
    if not re.fullmatch(r"[0-9][A-Z0-9]{3}", pdb_id):
        raise ValueError("Introduce un PDB ID de cuatro caracteres, por ejemplo 1CRN.")
    return pdb_id


def download_pdb(pdb_id: str) -> str:
    """Download a PDB file with a timeout and release the HTTP connection."""
    pdb_id = normalize_pdb_id(pdb_id)
    with urllib.request.urlopen(
        f"https://files.rcsb.org/download/{pdb_id}.pdb", timeout=20
    ) as response:
        return response.read().decode("utf-8")


def parse_pdb(text: str, pdb_id: str) -> ProteinData:
    """Read all chains and full-width coordinates once, using only the first model."""
    chains: dict[str, str] = {}
    atoms = {}
    titles = []
    compounds = []
    model_started = False

    for line_number, line in enumerate(text.splitlines(), start=1):
        record = line[:6].strip()
        if record == "TITLE":
            titles.append(line[10:80].strip())
        elif record == "COMPND":
            compounds.append(line[10:80].strip())
        elif record == "MODEL":
            if model_started:
                break
            model_started = True
        elif record in ("ENDMDL", "END"):
            break
        elif record in ("ATOM", "HETATM"):
            residue_name = line[17:20].strip()
            if record == "HETATM" and residue_name != "MSE":
                continue
            if len(line) < 54:
                raise ValueError(f"Registro de coordenadas incompleto en la línea {line_number}.")
            chain = line[21]
            atom_name = line[12:16].strip()
            atom = (chain, line[22:26], line[26], atom_name)
            # Prefer blank/A conformers; if a site only has B, keep it once.
            alternate = line[16]
            priority = 0 if alternate == " " else 1 if alternate == "A" else 2
            if atom in atoms and atoms[atom][0] <= priority:
                continue
            try:
                point = tuple(float(line[start:start + 8]) for start in (30, 38, 46))
            except ValueError as error:
                raise ValueError(f"Coordenadas inválidas en la línea {line_number}.") from error
            if not all(isfinite(value) for value in point):
                raise ValueError(f"Coordenadas no finitas en la línea {line_number}.")
            atoms[atom] = (priority, residue_name, point)

    for (chain, position, insertion, atom_name), (_, residue_name, _) in atoms.items():
        if atom_name == "CA":
            if residue_name not in AA3_TO_1:
                raise ValueError(
                    f"Residuo no compatible: {residue_name}, cadena {chain}, "
                    f"posición {position.strip()}{insertion.strip()}."
                )
            chains[chain] = chains.get(chain, "") + AA3_TO_1[residue_name]

    if not chains:
        raise ValueError("El archivo PDB no contiene aminoácidos con átomos CA.")
    title = " ".join(titles)
    if not title:
        molecules = re.findall(r"MOLECULE:\s*([^;]+)", " ".join(compounds))
        title = "; ".join(molecules) or pdb_id
    return ProteinData(title, chains, [point for _, _, point in atoms.values()])


def create_figures(protein: ProteinData, pdb_id: str, *, appearance: str = "dark"):
    """Build composition, atom and heuristic charts without opening a browser."""
    import plotly.graph_objs as go
    from plotly.subplots import make_subplots

    colors = PALETTES[appearance]
    group_colors = GROUP_PALETTES[appearance]
    sequence = "".join(protein.chains.values())
    counts = Counter(sequence)
    figure = make_subplots(
        rows=2, cols=1, specs=[[{"type": "xy"}], [{"type": "scene"}]],
        row_heights=[0.4, 0.6], vertical_spacing=0.16,
        subplot_titles=("Composición de aminoácidos", "Coordenadas atómicas · vista 3D"),
    )
    group_names = ("No polares", "Polares sin carga", "Carga negativa", "Carga positiva")
    for index, (_, amino_acids) in enumerate(AMINO_ACID_GROUPS):
        figure.add_trace(go.Bar(
            x=[
                f'<a href="https://www.aminoacidsguide.com/{page}.html">'
                f"{three}</a>"
                for one, three, full, page in amino_acids
            ],
            y=[100 * counts[one] / len(sequence) for one, *_ in amino_acids],
            customdata=[full for _, _, full, _ in amino_acids],
            hovertemplate="%{customdata}<br>%{y:.1f}%<extra>%{fullData.name}</extra>",
            marker_color=group_colors[index],
            name=group_names[index],
        ), row=1, col=1)
    x, y, z = zip(*protein.coordinates)
    figure.add_trace(go.Scatter3d(
        x=x, y=y, z=z, mode="markers", showlegend=False,
        marker=dict(size=3, color=z,
                    colorscale=[[0, colors["muted"]], [0.5, colors["accent"]],
                                [1, colors["gold"]]],
                    opacity=0.9),
        name="Átomos",
        hovertemplate="x: %{x:.3f} Å<br>y: %{y:.3f} Å<br>z: %{z:.3f} Å<extra>Átomo</extra>",
    ), row=2, col=1)
    groups = [[], [], [], []]
    for chain_sequence in protein.chains.values():
        result = countEstr(chain_sequence)
        for index in range(4):
            groups[index].extend(result[index * 2])

    theme = dict(
        template="plotly_dark" if appearance == "dark" else "plotly_white",
        paper_bgcolor=colors["surface"], plot_bgcolor=colors["surface"],
        font=dict(color=colors["text"], family=CHART_FONT, size=13),
        hoverlabel=dict(bgcolor=colors["elevated"],
                        font=dict(color=colors["text"], size=14)),
        modebar=dict(bgcolor=colors["surface"], color=colors["muted"],
                     activecolor=colors["accent"]),
    )
    figure.update_layout(
        **theme, height=940, autosize=True,
        margin=dict(r=24, t=120, b=24, l=56),
        scene=dict(
            bgcolor=colors["surface"],
            xaxis=dict(title="x (Å)", gridcolor=colors["border"], showbackground=False),
            yaxis=dict(title="y (Å)", gridcolor=colors["border"], showbackground=False),
            zaxis=dict(title="z (Å)", gridcolor=colors["border"], showbackground=False),
            aspectmode="data",
        ),
        legend=dict(orientation="h", x=0, y=1.12, yanchor="bottom"),
        xaxis=dict(tickangle=-45, automargin=True),
        yaxis=dict(title="Porcentaje (%)", ticksuffix="%", rangemode="tozero",
                   gridcolor=colors["border"], zerolinecolor=colors["border"]),
        bargap=0.25,
    )
    pie = go.Figure()
    pie.update_layout(
        **theme, height=440, margin=dict(l=24, r=24, t=24, b=72),
        legend=dict(orientation="h", x=0.5, xanchor="center", y=-0.12),
    )
    if any(groups):
        pie.add_trace(go.Pie(
            labels=["Alfa", "Beta", "Giro beta", "Azar"],
            values=[len(patterns) for patterns in groups],
            customdata=[" ".join(patterns) for patterns in groups],
            hole=0.6, sort=False, marker=dict(colors=group_colors),
            textinfo="percent", textposition="inside",
            insidetextfont=dict(color=colors["on_accent"]),
            hovertemplate="%{label}<br>%{value} tripletes · %{percent}<extra></extra>",
        ))
        pie.add_annotation(
            text=f"<b>{sum(map(len, groups))}</b><br>tripletes",
            x=0.5, y=0.5, showarrow=False, font_size=20,
        )
    else:
        pie.add_annotation(text="No hay tripletes completos.", showarrow=False)
    return figure, pie


def build_guipro(pdb_id: str, output_dir: Path = Path("."), *, auto_open: bool = True,
                 library_dir: Optional[Path] = None, refresh: bool = False,
                 appearance: str = "dark") -> None:
    """Download, analyze and save the PDB and both self-contained HTML charts."""
    pdb_id = normalize_pdb_id(pdb_id)
    if library_dir is None:
        text = download_pdb(pdb_id)
        protein = parse_pdb(text, pdb_id)
    else:
        from .library import load_structure
        text, protein = load_structure(pdb_id, library_dir, refresh=refresh)
    figure, pie = create_figures(protein, pdb_id, appearance=appearance)
    from .reports import write_reports

    (output_dir / f"{pdb_id}.pdb").write_text(text, encoding="utf-8")
    write_reports(protein, pdb_id, figure, pie, output_dir, auto_open=auto_open,
                  appearance=appearance)


def main() -> None:
    """Launch the Spanish-language desktop interface only when executed directly."""
    from .desktop import run
    run()


if __name__ == "__main__":
    main()
