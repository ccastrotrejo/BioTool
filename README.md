# BioTool

A Spanish-language Python desktop application for exploring protein structures
from the RCSB Protein Data Bank.

## Install and run

Use Python 3.9 or newer with Tk support. On Debian/Ubuntu, install the
`python3-tk` system package if your Python installation does not include Tk.
Other Python distributions may also require installing their Tk component.

On recent macOS versions, Apple's bundled Python/Tk 8.5 can display an empty
window. Use a current Python distribution with Tk instead, for example
`brew install python-tk@3.14`, and create the environment with `python3.14`.

Create and activate a virtual environment:

```sh
python -m venv .venv
```

On macOS/Linux, run `source .venv/bin/activate`. On Windows PowerShell, run
`.venv\Scripts\Activate.ps1`.

Install and launch from a source checkout:

```sh
python -m pip install .
biotool
```

`python gui.py` also launches the application from the checkout. Enter a
four-character PDB ID, such as `1CRN`, and select **Analizar estructura** or
press Enter. Alternatively, open **Biblioteca local** and select a starter
protein. An internet connection is required only for the first download or
an explicit refresh; saved structures can be analyzed offline.

BioTool saves the downloaded `<PDB_ID>.pdb`, `simple_plot.html`, and
`basic_pie_chart.html` in the current working directory. The HTML charts include
Plotly and open in your browser. Running another analysis replaces the two chart
files.

The desktop window includes input guidance, analysis status, and a reminder of
the output behavior. The linked, Spanish-language HTML reports include stacked
composition/3D charts, a readable per-chain FASTA section, composition and triplet
count tables, and an explanation of the heuristic's limitations. Reports adapt
to the browser width; the interactive charts also include image export controls.
Use **Claro / Oscuro** to switch the desktop appearance. Newly generated reports
start in that appearance and also have their own light/dark selector. The palettes
live in `biotool/theme.py`. Reports use Liquid Glass-inspired translucent surfaces
and backdrop blur, with opaque fallbacks for unsupported browsers and reduced
transparency preferences. Tkinter uses a matching opaque, glass-inspired theme:
it does not support Apple's native Liquid Glass refraction or background blur.
Both desktop themes explicitly style controls to avoid mixed system colors.

## Local protein library

The library is a small educational catalog, not a download of the entire PDB.
Catalog names and descriptions are bundled with the application; coordinates
are downloaded from RCSB only when you choose a structure. **Ficha RCSB** opens
the original entry and citation information.

The starter catalog includes [crambin (1CRN)](https://www.rcsb.org/structure/1CRN),
[ubiquitin (1UBQ)](https://www.rcsb.org/structure/1UBQ),
[porcine insulin (4INS)](https://www.rcsb.org/structure/4INS),
[lysozyme (1LYZ)](https://www.rcsb.org/structure/1LYZ),
[myoglobin (1MBN)](https://www.rcsb.org/structure/1MBN),
[GFP (1GFL)](https://www.rcsb.org/structure/1GFL), and
[hemoglobin (2HHB)](https://www.rcsb.org/structure/2HHB).
These classic-PDB downloads were compatible with BioTool's parser when selected.
Catalog residue counts describe observed residues across the file's first-model
chains, not complete biological sequences or generated biological assemblies.
Heme, zinc and other non-MSE HETATM records are excluded from the visualization.

Every structure analyzed through the desktop, including a manually entered ID,
is saved in the per-user library after it passes BioTool's parser:

- macOS: `~/Library/Application Support/BioTool/pdb`
- Windows: `%LOCALAPPDATA%/BioTool/pdb`
- Linux: `$XDG_DATA_HOME/BioTool/pdb`, or `~/.local/share/BioTool/pdb`

**En este equipo** identifies saved files; **Por descargar** needs internet.
Saved files are reused without network access. **Actualizar desde RCSB**
explicitly downloads a fresh copy. Failed or incompatible downloads never replace
the previous file. A corrupt local copy produces an error with refresh guidance
rather than silently falling back to the network.

Downloads and report generation run in a worker so the desktop remains usable.
Only one analysis runs at a time. Closing the window stops UI updates; an active
download may finish in the worker before the process exits.

The Python `build_guipro` function preserves its uncached default for existing
callers. Pass `library_dir=Path(...)` to opt into caching and `refresh=True` to
replace a cached structure. `appearance="light"` generates light-theme reports.

### PDB sources and limitations

Downloads follow the [RCSB file-download services](https://www.rcsb.org/docs/programmatic-access/file-download-services).
PDB archive data is [CC0](https://www.wwpdb.org/about/usage-policies); please
[cite the structure authors and PDB entry](https://www.wwpdb.org/about/cite-us)
when using these structures in your work. This does not imply that website
illustrations or articles share the archive's license.

PDBx/mmCIF is the primary archive format. BioTool currently supports only
[entries available in legacy PDB format](https://www.rcsb.org/docs/general-help/structures-without-legacy-pdb-format-files).
It does not yet accept mmCIF or
[extended PDB identifiers](https://www.wwpdb.org/documentation/pdb-id-extension-faq).
The starter catalog is therefore deliberately small and compatible, rather
than a claim of support for every structure in the archive.

## Analysis scope

- Reads observed residues from all chains in the first PDB model; chain labels
  and artificial starting residues are not included in the sequence.
- Uses one coordinate per atom site, preferring blank/A alternate locations.
  Selenomethionine (`MSE`) is treated as methionine.
- Reports amino-acid composition and displays atomic coordinates without
  implying bonds between consecutive records.
- Preserves the original secondary-structure propensity rules. This is a
  heuristic classification, not an experimentally determined structure or a
  validated secondary-structure predictor.
- Counts every complete, non-overlapping triplet independently within each
  chain. One or two trailing residues are omitted from this classification,
  but remain in the composition calculation.

Sequences reflect resolved CA atoms, not the full biological sequence.
Unsupported residues with CA atoms are reported as errors rather than silently
changing the analysis. This version accepts classic PDB files, not mmCIF.

## Tests

After installing the application dependencies:

```sh
python -m unittest discover -v
```

Tests cover parsing, classification, composition, downloads, chart generation,
and desktop callbacks. They use synthetic PDB records and mocked network/browser
access, and do not require a display server.

## Build

```sh
python -m pip install build
python -m build
```

This creates a source archive and a platform-independent wheel in `dist/`.
The package includes the application logo and the `biotool` launcher. Install
the wheel with `python -m pip install dist/biotool-0.1.0-py3-none-any.whl`.
Python, Tk, and the declared Python dependencies are still required; the wheel
is not a standalone desktop executable.

GitHub Actions runs the tests on Linux, macOS, and Windows for pushes and pull
requests. After the tests pass, it builds both distributions, checks the
installed wheel outside the checkout, and uploads the
`biotool-distributions` artifact. The workflow can also be started manually.
