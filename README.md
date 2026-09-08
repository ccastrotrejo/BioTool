# BioTool

A Spanish-language Python desktop application for exploring protein structures
from the RCSB Protein Data Bank.

## Install and run

Use Python 3.9 or newer with Tk support. On Debian/Ubuntu, install the
`python3-tk` system package if your Python installation does not include Tk.
Other Python distributions may also require installing their Tk component.

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

`python gui.py` also launches the application from the checkout. An internet
connection is required to download structures. Enter a four-character PDB ID,
such as `1CRN`, and select **Desplegar Información**.

BioTool saves the downloaded `<PDB_ID>.pdb`, `simple_plot.html`, and
`basic_pie_chart.html` in the current working directory. The HTML charts include
Plotly and open in your browser. Running another analysis replaces the two chart
files.

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
