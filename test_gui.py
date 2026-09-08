"""Offline regression coverage: python -m unittest -v."""

from contextlib import ExitStack
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import Mock, patch
import urllib.error

from biotool import app as gui


def atom_line(
    serial=1, name="CA", residue="ALA", chain="A", position=1,
    point=(-123.456, 1234.567, -987.654), record="ATOM",
    alternate=" ", insertion=" ",
):
    return (
        f"{record:<6}{serial:5d} {name:^4}{alternate}{residue:>3} "
        f"{chain}{position:4d}{insertion}   "
        f"{point[0]:8.3f}{point[1]:8.3f}{point[2]:8.3f}"
        "  1.00 20.00           C"
    )


def pdb_for_chains(chains):
    lines = ["TITLE     EXAMPLE PROTEIN"]
    for chain, residues in chains.items():
        for position, residue in enumerate(residues, start=1):
            lines.append(atom_line(
                serial=len(lines), residue=residue, chain=chain, position=position,
            ))
    return "\n".join(lines + ["END"])


class StructureTests(unittest.TestCase):
    def test_original_propensities_and_thresholds(self):
        self.assertEqual(gui.probAA("AAA"), {"AAA": [1.25, 0.89, 0.78]})
        self.assertEqual(gui.estrAA("AAAVVVGGGAGM"), {
            "AAA": "alfa", "VVV": "beta", "GGG": "giro beta", "AGM": "azar",
        })

    def test_repeated_triplets_count_every_occurrence(self):
        result = gui.countEstr("AAAAAAVVV")
        self.assertEqual(result, [["AAA", "AAA"], 2, ["VVV"], 1, [], 0, [], 0])

    def test_counts_do_not_depend_on_a_random_structure(self):
        self.assertEqual(gui.countEstr("AGMAAA")[-1], 1)
        self.assertEqual(gui.countEstr("AGMAAA")[1], 1)
        self.assertEqual(gui.countEstr("AAA")[1], 1)

    def test_empty_short_and_incomplete_triplets(self):
        for sequence in ("", "A", "AA"):
            with self.subTest(sequence=sequence):
                self.assertEqual(gui.countEstr(sequence), [[], 0, [], 0, [], 0, [], 0])
        self.assertEqual(gui.countEstr("AAAAV")[1], 1)

    def test_invalid_residues_are_reported_even_in_a_trailing_fragment(self):
        for sequence in ("AAZ", "AAA?", "aaa"):
            with self.subTest(sequence=sequence), self.assertRaisesRegex(ValueError, "no reconocidos"):
                gui.countEstr(sequence)


class ParserTests(unittest.TestCase):
    def test_all_chains_without_synthetic_residues(self):
        protein = gui.parse_pdb(pdb_for_chains({
            "A": ["ALA", "GLY"], "B": ["VAL"], " ": ["SER"],
        }), "1ABC")
        self.assertEqual(protein.chains, {"A": "AG", "B": "V", " ": "S"})
        self.assertEqual(len(protein.coordinates), 4)
        self.assertEqual(protein.title, "EXAMPLE PROTEIN")

    def test_full_width_signed_coordinates(self):
        protein = gui.parse_pdb(atom_line(), "1ABC")
        self.assertEqual(protein.coordinates, [(-123.456, 1234.567, -987.654)])

    def test_only_first_model_is_read(self):
        text = "\n".join([
            "MODEL        1", atom_line(), "ENDMDL", "MODEL        2",
            atom_line(residue="GLY"), "ENDMDL",
        ])
        protein = gui.parse_pdb(text, "1ABC")
        self.assertEqual(protein.chains, {"A": "A"})
        self.assertEqual(len(protein.coordinates), 1)

    def test_insertion_codes_are_distinct_residues(self):
        text = "\n".join([
            atom_line(residue="ALA"), atom_line(residue="GLY", insertion="A"),
        ])
        self.assertEqual(gui.parse_pdb(text, "1ABC").chains, {"A": "AG"})

    def test_alternate_locations_prefer_a_regardless_of_record_order(self):
        preferred = atom_line(alternate="A", point=(1, 2, 3))
        alternate = atom_line(alternate="B", point=(4, 5, 6))
        for records in ((alternate, preferred), (preferred, alternate)):
            with self.subTest(records=records):
                protein = gui.parse_pdb("\n".join(records), "1ABC")
                self.assertEqual(protein.chains, {"A": "A"})
                self.assertEqual(protein.coordinates, [(1, 2, 3)])

    def test_blank_and_b_only_conformers(self):
        records = [
            atom_line(alternate="A", point=(1, 2, 3)),
            atom_line(point=(4, 5, 6)),
            atom_line(position=2, alternate="B", point=(7, 8, 9)),
        ]
        protein = gui.parse_pdb("\n".join(records), "1ABC")
        self.assertEqual(protein.chains, {"A": "AA"})
        self.assertEqual(protein.coordinates, [(4, 5, 6), (7, 8, 9)])

    def test_mse_is_methionine_and_water_is_excluded(self):
        text = "\n".join([
            atom_line(record="HETATM", residue="MSE"),
            atom_line(record="HETATM", residue="HOH", name="O", position=2),
        ])
        protein = gui.parse_pdb(text, "1ABC")
        self.assertEqual(protein.chains, {"A": "M"})
        self.assertEqual(len(protein.coordinates), 1)

    def test_non_ca_atoms_are_plotted_but_not_counted_as_residues(self):
        text = "\n".join([atom_line(), atom_line(name="N")])
        protein = gui.parse_pdb(text, "1ABC")
        self.assertEqual(protein.chains, {"A": "A"})
        self.assertEqual(len(protein.coordinates), 2)

    def test_missing_protein_and_unknown_residues_are_reported(self):
        for text in ("", "<html>Not a PDB</html>", atom_line(residue="UNK")):
            with self.subTest(text=text), self.assertRaises(ValueError):
                gui.parse_pdb(text, "1ABC")

    def test_bad_coordinates_are_reported(self):
        for text in (
            "ATOM      1", atom_line(point=(float("nan"), 1, 2)),
            atom_line()[:30] + "invalid " + atom_line()[38:],
        ):
            with self.subTest(text=text), self.assertRaises(ValueError):
                gui.parse_pdb(text, "1ABC")

    def test_multiline_title_and_compound_fallback(self):
        text = "\n".join(["TITLE     FIRST", "TITLE    2 SECOND", atom_line()])
        self.assertEqual(gui.parse_pdb(text, "1ABC").title, "FIRST SECOND")
        text = "\n".join([
            "COMPND    MOL_ID: 1;", "COMPND   2 MOLECULE: TEST",
            "COMPND   3 PROTEIN;", atom_line(),
        ])
        self.assertEqual(gui.parse_pdb(text, "1ABC").title, "TEST PROTEIN")
        self.assertEqual(gui.parse_pdb(atom_line(), "1ABC").title, "1ABC")


class DownloadTests(unittest.TestCase):
    def test_normalizes_id(self):
        self.assertEqual(gui.normalize_pdb_id(" 1crn "), "1CRN")

    @patch("urllib.request.urlopen")
    def test_invalid_id_never_reaches_network(self, urlopen):
        for value in ("", "../file", "1ABCD", "ABCD", "1A/B", "1 A1"):
            with self.subTest(value=value), self.assertRaises(ValueError):
                gui.download_pdb(value)
        urlopen.assert_not_called()

    @patch("urllib.request.urlopen")
    def test_timeout_and_connection_cleanup(self, urlopen):
        response = urlopen.return_value.__enter__.return_value
        response.read.return_value = b"HEADER example"
        self.assertEqual(gui.download_pdb("1crn"), "HEADER example")
        urlopen.assert_called_once_with(
            "https://files.rcsb.org/download/1CRN.pdb", timeout=20,
        )
        urlopen.return_value.__exit__.assert_called_once()

    @patch("urllib.request.urlopen", side_effect=urllib.error.URLError("offline"))
    def test_network_errors_propagate(self, urlopen):
        with self.assertRaises(urllib.error.URLError):
            gui.download_pdb("1ABC")


class ChartTests(unittest.TestCase):
    def test_composition_and_structure_include_all_chains(self):
        protein = gui.parse_pdb(pdb_for_chains({
            "A": ["ALA"] * 6, "B": ["VAL"] * 3,
        }), "1ABC")
        figure, pie = gui.create_figures(protein, "1ABC")
        percentages = [value for trace in figure.data[:4] for value in trace.y]
        self.assertAlmostEqual(sum(percentages), 100)
        self.assertAlmostEqual(figure.data[0].y[1], 100 * 6 / 9)
        self.assertAlmostEqual(figure.data[0].y[2], 100 * 3 / 9)
        self.assertEqual(tuple(pie.data[0].values), (2, 1, 0, 0))
        self.assertEqual(len(figure.data[4].x), 9)
        self.assertEqual(figure.data[4].mode, "markers")
        self.assertIn("AAAAAA", figure.layout.annotations[0].text)
        self.assertIn("VVV", figure.layout.annotations[0].text)
        self.assertIn("cadena B", figure.layout.annotations[0].text)

    def test_structure_triplets_never_cross_chain_boundaries(self):
        protein = gui.parse_pdb(pdb_for_chains({
            "A": ["ALA", "ALA"], "B": ["VAL"],
        }), "1ABC")
        _, pie = gui.create_figures(protein, "1ABC")
        self.assertEqual(len(pie.data), 0)
        self.assertEqual(pie.layout.annotations[0].text, "No hay tripletes completos.")

    def test_metadata_is_escaped_in_chart_markup(self):
        protein = gui.parse_pdb(
            "TITLE     <b>NOT MARKUP</b>\n" + atom_line(chain="<"), "1ABC",
        )
        figure, _ = gui.create_figures(protein, "1ABC")
        self.assertIn("&lt;b&gt;NOT MARKUP&lt;/b&gt;", figure.layout.title.text)
        self.assertIn("cadena &lt;", figure.layout.annotations[0].text)

    @patch("urllib.request.urlopen")
    def test_end_to_end_writes_pdb_and_valid_html_without_browser(self, urlopen):
        text = pdb_for_chains({"A": ["ALA"] * 3})
        urlopen.return_value.__enter__.return_value.read.return_value = text.encode("utf-8")
        with tempfile.TemporaryDirectory() as directory, patch("webbrowser.open") as browser:
            output = Path(directory)
            gui.build_guipro("1abc", output, auto_open=False)
            self.assertEqual((output / "1ABC.pdb").read_text(encoding="utf-8"), text)
            for filename in ("simple_plot.html", "basic_pie_chart.html"):
                html = (output / filename).read_text(encoding="utf-8")
                self.assertIn("<html>", html)
                self.assertIn("Plotly.newPlot", html)
            browser.assert_not_called()


class StartupTests(unittest.TestCase):
    def test_import_does_not_load_gui_or_plot_dependencies(self):
        result = subprocess.run([
            sys.executable, "-S", "-c",
            "import gui, biotool.app, sys; assert 'tkinter' not in sys.modules; "
            "assert 'plotly' not in sys.modules; assert 'PIL' not in sys.modules",
        ], cwd=Path(__file__).parent, capture_output=True, text=True, timeout=10)
        self.assertEqual(result.returncode, 0, result.stderr)


class DesktopTests(unittest.TestCase):
    def setUp(self):
        stack = ExitStack()
        self.addCleanup(stack.close)
        directory = stack.enter_context(tempfile.TemporaryDirectory())
        previous_directory = Path.cwd()
        stack.callback(os.chdir, previous_directory)
        os.chdir(directory)
        self.widgets = {
            name: stack.enter_context(patch(f"tkinter.{name}"))
            for name in ("Tk", "Label", "Entry", "StringVar", "Button")
        }
        self.entry_text = Mock()
        self.entry_text.get.return_value = "1ABC"
        self.label_text = Mock()
        self.widgets["StringVar"].side_effect = [self.entry_text, self.label_text]
        self.open_image = stack.enter_context(patch("PIL.Image.open"))
        stack.enter_context(patch("PIL.ImageTk.PhotoImage"))
        self.showerror = stack.enter_context(patch("tkinter.messagebox.showerror"))
        self.browser = stack.enter_context(patch("webbrowser.open"))
        self.urlopen = stack.enter_context(patch("urllib.request.urlopen"))
        self.urlopen.return_value.__enter__.return_value.read.return_value = (
            pdb_for_chains({"A": ["ALA"] * 3}).encode("utf-8")
        )
        gui.main()
        self.analyze = self.widgets["Button"].call_args.kwargs["command"]

    def test_lifecycle_and_asset_path_are_independent_of_working_directory(self):
        self.widgets["Tk"].return_value.mainloop.assert_called_once()
        self.open_image.assert_called_once_with(Path(gui.__file__).with_name("3.png"))

    def test_success_generates_both_charts_and_restores_button(self):
        self.analyze()
        self.assertTrue(Path("1ABC.pdb").is_file())
        self.assertTrue(Path("simple_plot.html").is_file())
        self.assertTrue(Path("basic_pie_chart.html").is_file())
        self.assertEqual(self.browser.call_count, 2)
        self.showerror.assert_not_called()
        self.widgets["Button"].return_value.configure.assert_called_with(state="normal")
        self.assertIn("completado", self.label_text.set.call_args.args[0])

    def test_invalid_input_is_displayed_without_network_access(self):
        self.entry_text.get.return_value = "../invalid"
        self.analyze()
        self.showerror.assert_called_once()
        self.urlopen.assert_not_called()
        self.browser.assert_not_called()
        self.widgets["Button"].return_value.configure.assert_called_with(state="normal")

    def test_network_and_file_errors_are_visible_and_retryable(self):
        for error in (
            urllib.error.URLError("offline"), TimeoutError("timeout"),
            PermissionError("read-only folder"),
            urllib.error.HTTPError("https://files.rcsb.org", 404, "Not found", {}, None),
        ):
            with self.subTest(error=error):
                self.showerror.reset_mock()
                self.urlopen.side_effect = error
                self.analyze()
                self.showerror.assert_called_once()
                self.browser.assert_not_called()
                self.widgets["Button"].return_value.configure.assert_called_with(state="normal")
                self.assertIn("No se pudo", self.label_text.set.call_args.args[0])

    def test_failed_file_write_does_not_open_browser(self):
        with patch.object(Path, "write_text", side_effect=PermissionError("read-only folder")):
            self.analyze()
        self.showerror.assert_called_once()
        self.browser.assert_not_called()

    def test_unexpected_errors_are_not_swallowed(self):
        self.urlopen.side_effect = RuntimeError("unexpected failure")
        with self.assertRaisesRegex(RuntimeError, "unexpected failure"):
            self.analyze()
        self.showerror.assert_not_called()
        self.widgets["Button"].return_value.configure.assert_called_with(state="normal")


if __name__ == "__main__":
    unittest.main()
