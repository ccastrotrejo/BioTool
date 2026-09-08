#!/usr/bin/env python3
"""Analyze observed protein residues from the first model of an RCSB PDB file.

Install with ``pip install .`` and launch with ``biotool`` or ``python gui.py``.
The desktop interface requires a Python installation with Tk support.
"""

from collections import Counter
from dataclasses import dataclass
from html import escape
from math import isfinite
from pathlib import Path
import re
import urllib.error
import urllib.request


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


def create_figures(protein: ProteinData, pdb_id: str):
    """Build composition, atom and heuristic charts without opening a browser."""
    import plotly.graph_objs as go

    sequence = "".join(protein.chains.values())
    counts = Counter(sequence)
    traces = []
    for name, amino_acids in AMINO_ACID_GROUPS:
        traces.append(go.Bar(
            x=[
                f'<a href="https://www.aminoacidsguide.com/{page}.html">'
                f"{three}({one}): {full}</a>"
                for one, three, full, page in amino_acids
            ],
            y=[100 * counts[one] / len(sequence) for one, *_ in amino_acids],
            name=name,
        ))
    x, y, z = zip(*protein.coordinates)
    traces.append(go.Scatter3d(
        x=x, y=y, z=z, mode="markers", showlegend=False,
        marker=dict(size=4, color=z, colorscale="Viridis"),
        name="Átomos",
    ))
    fasta = []
    groups = [[], [], [], []]
    for chain, chain_sequence in protein.chains.items():
        fasta.append(escape(f">{pdb_id}|cadena {chain.strip() or '(sin identificador)'}"))
        fasta.extend(chain_sequence[i:i + 80] for i in range(0, len(chain_sequence), 80))
        result = countEstr(chain_sequence)
        for index in range(4):
            groups[index].extend(result[index * 2])

    figure = go.Figure(data=traces, layout=dict(
        plot_bgcolor="black", paper_bgcolor="black", font=dict(color="white"),
        title=dict(
            text=f"{escape(protein.title)}<br>Número total de aminoácidos: {len(sequence)}",
            font=dict(size=15, family="Raleway"),
        ),
        margin=dict(r=10, t=130, b=20, l=50),
        scene=dict(
            domain=dict(x=[0.52, 0.97], y=[0.3, 1]),
            xaxis=dict(gridcolor="white"), yaxis=dict(gridcolor="white"),
            zaxis=dict(gridcolor="white"), aspectmode="data",
        ),
        showlegend=True, legend=dict(x=0, y=1.2),
        xaxis=dict(anchor="y", domain=[0.01, 0.45]),
        yaxis=dict(anchor="x", domain=[0.26, 0.95], showgrid=False, title="Porcentaje"),
        annotations=[dict(
            text="FASTA (residuos observados, primer modelo):<br>" + "<br>".join(fasta)
            + '<br><a href="https://www.bachem.com/fileadmin/user_upload/pdf/Flyers/'
            'Periodic_Chart_Amino_Acids.pdf"><b>Tablas de Aminoácidos</b></a>',
            showarrow=False, xref="paper", yref="paper", x=1, y=0,
        )],
    ))
    labels = [
        f"{name}: {' '.join(patterns)}"
        for name, patterns in zip(("Alfa", "Beta", "Beta giro", "Azar"), groups)
    ]
    pie = go.Figure()
    pie.update_layout(title=dict(
        text="Clasificación heurística de tripletes (no experimental)"
        "<br>Tripletes no solapados por cadena; se omiten residuos finales incompletos."
    ))
    if any(groups):
        pie.add_trace(go.Pie(labels=labels, values=[len(patterns) for patterns in groups]))
    else:
        pie.add_annotation(text="No hay tripletes completos.", showarrow=False)
    return figure, pie


def build_guipro(pdb_id: str, output_dir: Path = Path("."), *, auto_open: bool = True) -> None:
    """Download, analyze and save the PDB and both self-contained HTML charts."""
    pdb_id = normalize_pdb_id(pdb_id)
    text = download_pdb(pdb_id)
    protein = parse_pdb(text, pdb_id)
    figure, pie = create_figures(protein, pdb_id)
    (output_dir / f"{pdb_id}.pdb").write_text(text, encoding="utf-8")
    figure.write_html(str(output_dir / "simple_plot.html"), auto_open=auto_open)
    pie.write_html(str(output_dir / "basic_pie_chart.html"), auto_open=auto_open)


def main() -> None:
    """Launch the Spanish-language desktop interface only when executed directly."""
    import tkinter as tk
    from tkinter import messagebox
    from PIL import Image, ImageTk

    root = tk.Tk()
    root.title("BioTool")
    with Image.open(Path(__file__).with_name("3.png")) as image:
        logo = ImageTk.PhotoImage(image)
    tk.Label(root, image=logo).pack(side="top", fill="both", expand=True)
    tk.Label(root, text="Herramienta bioinformática para analizar proteínas").pack()
    tk.Label(root, text="PDB ID:").pack()
    entry_text = tk.StringVar(root)
    entry = tk.Entry(root, width=10, textvariable=entry_text)
    entry.pack()
    label_text = tk.StringVar(root)

    def analyze() -> None:
        button.configure(state="disabled")
        label_text.set("Descargando y analizando...")
        root.update_idletasks()
        try:
            build_guipro(entry_text.get())
        except urllib.error.HTTPError as error:
            label_text.set("No se pudo descargar el archivo PDB.")
            messagebox.showerror("BioTool", f"RCSB respondió con HTTP {error.code}.", parent=root)
        except (OSError, ValueError) as error:
            label_text.set("No se pudo completar el análisis.")
            messagebox.showerror("BioTool", str(error), parent=root)
        else:
            label_text.set("Análisis completado. Gráficas HTML guardadas en la carpeta actual.")
        finally:
            button.configure(state="normal")

    button = tk.Button(root, text="Desplegar Información", command=analyze, height=3, width=35)
    button.pack()
    tk.Label(root, textvariable=label_text).pack()
    entry.focus_set()
    root.mainloop()


if __name__ == "__main__":
    main()
