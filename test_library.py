"""Offline regression coverage for the on-demand PDB library."""

from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import urllib.error

from biotool.library import default_library_dir, load_structure, saved_ids
from test_gui import pdb_for_chains


class LocalLibraryTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.directory = Path(temporary.name) / "library"
        self.text = pdb_for_chains({"A": ["ALA", "GLY", "VAL"]})

    @patch("biotool.app.download_pdb")
    def test_download_once_then_reuse_offline(self, download):
        download.return_value = self.text
        first = load_structure(" 1abc ", self.directory)
        download.side_effect = urllib.error.URLError("offline")
        second = load_structure("1ABC", self.directory)
        self.assertEqual(first, second)
        self.assertEqual(first[1].chains, {"A": "AGV"})
        download.assert_called_once_with("1ABC")
        self.assertEqual(saved_ids(self.directory), ["1ABC"])

    @patch("biotool.app.download_pdb")
    def test_invalid_input_never_touches_network_or_disk(self, download):
        with self.assertRaises(ValueError):
            load_structure("../path", self.directory)
        download.assert_not_called()
        self.assertFalse(self.directory.exists())

    @patch("biotool.app.download_pdb", return_value="<html>not a PDB</html>")
    def test_invalid_download_is_not_cached(self, download):
        with self.assertRaises(ValueError):
            load_structure("1ABC", self.directory)
        self.assertFalse(self.directory.exists())

    @patch("biotool.app.download_pdb")
    def test_refresh_keeps_old_copy_on_network_or_parse_failure(self, download):
        download.return_value = self.text
        load_structure("1ABC", self.directory)
        download.side_effect = urllib.error.URLError("offline")
        with self.assertRaises(urllib.error.URLError):
            load_structure("1ABC", self.directory, refresh=True)
        download.side_effect = None
        download.return_value = "not a PDB"
        with self.assertRaises(ValueError):
            load_structure("1ABC", self.directory, refresh=True)
        self.assertEqual((self.directory / "1ABC.pdb").read_text(), self.text)

    @patch("biotool.app.download_pdb")
    def test_refresh_replaces_valid_structure_atomically(self, download):
        download.return_value = self.text
        load_structure("1ABC", self.directory)
        download.return_value = pdb_for_chains({"B": ["GLY"]})
        _, protein = load_structure("1ABC", self.directory, refresh=True)
        self.assertEqual(protein.chains, {"B": "G"})
        self.assertEqual(len(list(self.directory.iterdir())), 1)

    @patch("biotool.app.download_pdb")
    def test_corrupt_local_copy_requires_explicit_refresh(self, download):
        self.directory.mkdir()
        (self.directory / "1ABC.pdb").write_text("broken", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "Actualizar desde RCSB"):
            load_structure("1ABC", self.directory)
        download.assert_not_called()
        download.return_value = self.text
        self.assertEqual(load_structure("1ABC", self.directory, refresh=True)[0], self.text)

    @patch("biotool.app.download_pdb")
    def test_failed_atomic_replace_removes_temporary_file(self, download):
        download.return_value = self.text
        with patch.object(Path, "replace", side_effect=PermissionError("read-only")):
            with self.assertRaises(PermissionError):
                load_structure("1ABC", self.directory)
        self.assertEqual(list(self.directory.iterdir()), [])

    def test_listing_ignores_unrelated_files(self):
        self.assertEqual(saved_ids(self.directory), [])
        self.directory.mkdir()
        for filename in ("1ABC.pdb", "notes.txt", ".1ABC.tmp", "ABC.pdb", "1!BC.pdb"):
            (self.directory / filename).touch()
        self.assertEqual(saved_ids(self.directory), ["1ABC"])

    def test_platform_data_locations(self):
        with patch("biotool.library.sys.platform", "darwin"), patch.object(
            Path, "home", return_value=Path("/home/test")
        ):
            self.assertEqual(default_library_dir(),
                             Path("/home/test/Library/Application Support/BioTool/pdb"))
        with patch("biotool.library.sys.platform", "linux"), patch.dict(
            "os.environ", {"XDG_DATA_HOME": "/custom/data"}
        ):
            self.assertEqual(default_library_dir(), Path("/custom/data/BioTool/pdb"))
        with patch("biotool.library.sys.platform", "win32"), patch.dict(
            "os.environ", {"LOCALAPPDATA": "/custom/local"}
        ):
            self.assertEqual(default_library_dir(), Path("/custom/local/BioTool/pdb"))


if __name__ == "__main__":
    unittest.main()
