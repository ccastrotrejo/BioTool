"""A small, on-demand local PDB library; never a mirror of the PDB archive."""

import os
from pathlib import Path
import re
import sys
import tempfile


def default_library_dir() -> Path:
    """Return the platform's per-user application-data folder without creating it."""
    if sys.platform == "darwin":
        base = Path.home() / "Library" / "Application Support"
    elif sys.platform == "win32":
        base = Path(os.environ.get("LOCALAPPDATA", Path.home() / "AppData" / "Local"))
    else:
        base = Path(os.environ.get("XDG_DATA_HOME", Path.home() / ".local" / "share"))
    return base / "BioTool" / "pdb"


def saved_ids(directory: Path) -> list[str]:
    """List only canonical PDB filenames owned by the local library."""
    if not directory.exists():
        return []
    return sorted(
        path.stem for path in directory.iterdir()
        if path.is_file() and re.fullmatch(r"[0-9][A-Z0-9]{3}\.pdb", path.name)
    )


def load_structure(pdb_id: str, directory: Path, *, refresh: bool = False):
    """Read offline first, or validate a download before atomically caching it."""
    from .app import download_pdb, normalize_pdb_id, parse_pdb

    pdb_id = normalize_pdb_id(pdb_id)
    destination = directory / f"{pdb_id}.pdb"
    if destination.exists() and not refresh:
        text = destination.read_text(encoding="utf-8")
        try:
            protein = parse_pdb(text, pdb_id)
        except ValueError as error:
            raise ValueError(
                f"The local copy of {pdb_id} is incompatible. "
                'Select "Refresh from RCSB" in the library to replace it.'
            ) from error
        return text, protein

    text = download_pdb(pdb_id)
    protein = parse_pdb(text, pdb_id)
    directory.mkdir(parents=True, exist_ok=True)
    temporary_path = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w", encoding="utf-8", dir=directory, prefix=f".{pdb_id}-",
            suffix=".tmp", delete=False,
        ) as temporary:
            temporary_path = Path(temporary.name)
            temporary.write(text)
        temporary_path.replace(destination)
    finally:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)
    return text, protein
