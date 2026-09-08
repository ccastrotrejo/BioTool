"""Native desktop workflow with offline library, themes and nonblocking analysis."""

from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import tkinter as tk
from tkinter import font, messagebox, ttk
import urllib.error
import webbrowser

from PIL import Image, ImageTk

from . import app
from .library import default_library_dir
from .library_view import LibraryView
from .theme import PALETTES, SPACE, TYPE, configure_desktop_styles


class Desktop:
    """Keep all widget updates on Tk's thread while analyses run in one worker."""

    def __init__(self, root):
        self.root = root
        self.executor = ThreadPoolExecutor(max_workers=1)
        self.future = None
        self.poll_id = None
        self.library_dir = default_library_dir()
        self.appearance = tk.StringVar(root, value="dark")
        self.entry_text = tk.StringVar(root)
        self.label_text = tk.StringVar(root, value="Ready to analyze a structure.")
        self.family = font.nametofont("TkDefaultFont").actual("family")
        self.style = ttk.Style(root)
        root.title("BioTool")
        root.minsize(540, 780)
        root.geometry("720x840")
        root.columnconfigure(0, weight=1)
        root.rowconfigure(0, weight=1)
        self.change_theme()
        content = ttk.Frame(root, padding=SPACE["lg"], style="Bio.TFrame")
        content.grid(row=0, column=0, sticky="nsew")
        content.columnconfigure(0, weight=1)
        content.rowconfigure(1, weight=1)
        self.build_header(content)
        self.tabs = ttk.Notebook(content)
        self.tabs.grid(row=1, column=0, sticky="nsew", pady=(16, 0))
        self.form = ttk.Frame(self.tabs, padding=24, style="Surface.TFrame")
        self.form.columnconfigure(0, weight=1)
        self.tabs.add(self.form, text="Analyze PDB")
        self.build_form()
        self.library = LibraryView(self.tabs, self.library_dir, self.analyze_selection)
        self.tabs.add(self.library.frame, text="Local library")
        self.status = ttk.Label(content, textvariable=self.label_text,
                                style="Status.TLabel", wraplength=560)
        self.status.grid(row=2, column=0, sticky="ew", pady=(16, 8))
        self.progress = ttk.Progressbar(content, mode="indeterminate")
        self.progress.grid(row=3, column=0, sticky="ew")
        self.progress.grid_remove()
        content.bind("<Configure>", self.resize_text)
        root.protocol("WM_DELETE_WINDOW", self.close)
        self.entry.focus_set()

    def build_header(self, parent):
        header = ttk.Frame(parent, style="Bio.TFrame")
        header.grid(row=0, column=0, sticky="ew")
        header.columnconfigure(1, weight=1)
        with Image.open(Path(__file__).with_name("3.png")) as image:
            image.thumbnail((64, 64), Image.Resampling.LANCZOS)
            self.logo = ImageTk.PhotoImage(image)
        ttk.Label(header, image=self.logo, style="Bio.TLabel").grid(
            row=0, column=0, rowspan=2, padx=(0, 16))
        ttk.Label(header, text="BioTool", style="Title.TLabel").grid(
            row=0, column=1, sticky="sw")
        ttk.Label(header, text="Protein explorer", style="Muted.TLabel").grid(
            row=1, column=1, sticky="nw")
        controls = ttk.Frame(header, style="Bio.TFrame")
        controls.grid(row=2, column=0, columnspan=2, sticky="e", pady=(16, 0))
        ttk.Label(controls, text="Appearance", style="Muted.TLabel").pack(side="left", padx=8)
        for label, value in (("Light", "light"), ("Dark", "dark")):
            ttk.Radiobutton(controls, text=label, value=value, variable=self.appearance,
                            command=self.change_theme).pack(side="left")

    def build_form(self):
        ttk.Label(self.form, text="Analyze a structure", style="Heading.TLabel").grid(
            row=0, column=0, sticky="w", pady=(0, 24))
        ttk.Label(self.form, text="PDB identifier", style="Surface.TLabel").grid(
            row=1, column=0, sticky="w", pady=(0, 8))
        self.entry = ttk.Entry(self.form, textvariable=self.entry_text,
                               font=(self.family, TYPE["heading"]), style="Pdb.TEntry")
        self.entry.grid(row=2, column=0, sticky="ew")
        self.entry.bind("<Return>", self.analyze)
        ttk.Label(
            self.form, text="For example, 1CRN. Downloaded once and saved in your library.",
            style="Hint.TLabel", wraplength=420,
        ).grid(row=3, column=0, sticky="w", pady=(8, 16))
        self.button = ttk.Button(self.form, text="Analyze structure",
                                 style="Analyze.TButton", command=self.analyze)
        self.button.grid(row=4, column=0, sticky="ew")
        ttk.Separator(self.form).grid(row=5, column=0, sticky="ew", pady=24)
        ttk.Label(self.form, text="Results in your browser", style="Heading.TLabel").grid(
            row=6, column=0, sticky="w", pady=(0, 8))
        self.notes = ttk.Label(
            self.form, style="Hint.TLabel", wraplength=420, justify="left",
            text="Amino acid composition, 3D coordinates, and observed sequence.\n\n"
            "Triplet classification is a heuristic, not a validated prediction.\n\n"
            "The PDB file and two HTML reports are saved in the current folder. "
            "A new analysis replaces the previous reports.",
        )
        self.notes.grid(row=7, column=0, sticky="ew")

    def change_theme(self):
        appearance = self.appearance.get()
        configure_desktop_styles(self.style, self.family, appearance)
        self.root.configure(background=PALETTES[appearance]["background"])

    def resize_text(self, event):
        self.status.configure(wraplength=max(240, event.width - 80))
        self.notes.configure(wraplength=max(240, event.width - 96))

    def analyze_selection(self, pdb_id, *, refresh=False):
        if self.future is not None:
            self.label_text.set("An analysis is in progress. Please wait for it to finish.")
            return
        self.entry_text.set(pdb_id)
        self.tabs.select(self.form)
        self.analyze(refresh=refresh)

    def analyze(self, event=None, *, refresh=False):
        if self.future is not None:
            self.label_text.set("An analysis is in progress. Please wait for it to finish.")
            return
        try:
            pdb_id = app.normalize_pdb_id(self.entry_text.get())
        except ValueError as error:
            self.show_error(error)
            return
        self.entry_text.set(pdb_id)
        self.button.configure(state="disabled")
        self.entry.configure(state="disabled")
        self.label_text.set(f"Analyzing {pdb_id} · local library / RCSB…")
        self.progress.grid()
        self.progress.start(12)
        self.future = self.executor.submit(
            app.build_guipro, pdb_id, auto_open=False, library_dir=self.library_dir,
            refresh=refresh, appearance=self.appearance.get(),
        )
        self.poll()

    def poll(self):
        if not self.future.done():
            self.poll_id = self.root.after(75, self.poll)
            return
        self.poll_id = None
        try:
            self.future.result()
            self.library.refresh()
        except (OSError, ValueError) as error:
            self.show_error(error)
        else:
            self.label_text.set("Analysis complete. Reports saved in the current folder.")
            for filename in ("simple_plot.html", "basic_pie_chart.html"):
                webbrowser.open(Path(filename).resolve().as_uri())
        finally:
            self.future = None
            self.progress.stop()
            self.progress.grid_remove()
            self.button.configure(state="normal")
            self.entry.configure(state="normal")
            self.entry.focus_set()

    def show_error(self, error):
        self.label_text.set("Could not complete the analysis. Fix the issue and try again.")
        detail = f"RCSB returned HTTP {error.code}." if isinstance(
            error, urllib.error.HTTPError) else str(error)
        messagebox.showerror("BioTool", detail, parent=self.root)
        self.entry.focus_set()

    def close(self):
        if self.poll_id is not None:
            self.root.after_cancel(self.poll_id)
        self.executor.shutdown(wait=False, cancel_futures=True)
        self.root.destroy()


def run() -> None:
    """Launch the native application."""
    root = tk.Tk()
    Desktop(root)
    root.mainloop()
