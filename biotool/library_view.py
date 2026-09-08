"""Keyboard-accessible browsing of starter proteins and saved PDB structures."""

from pathlib import Path
import tkinter as tk
from tkinter import ttk
import webbrowser

from .catalog import EXAMPLES
from .library import saved_ids


class LibraryView:
    """Display catalog metadata without downloading structures until requested."""

    def __init__(self, parent, directory: Path, analyze):
        self.directory = directory
        self.analyze = analyze
        self.frame = ttk.Frame(parent, padding=24, style="Surface.TFrame")
        self.frame.columnconfigure(0, weight=1)
        self.frame.rowconfigure(2, weight=1)
        ttk.Label(self.frame, text="Protein library", style="Heading.TLabel").grid(
            row=0, column=0, sticky="w", pady=(0, 8))
        ttk.Label(
            self.frame, text="Choose an example. After the first download, it works offline.",
            style="Hint.TLabel", wraplength=420,
        ).grid(row=1, column=0, sticky="w", pady=(0, 16))
        listing = ttk.Frame(self.frame, style="Surface.TFrame")
        listing.grid(row=2, column=0, sticky="nsew")
        listing.columnconfigure(0, weight=1)
        listing.rowconfigure(0, weight=1)
        self.tree = ttk.Treeview(
            listing, columns=("name", "saved"), show="tree headings", selectmode="browse",
            height=6,
        )
        self.tree.heading("#0", text="PDB")
        self.tree.heading("name", text="Protein")
        self.tree.heading("saved", text="Availability")
        self.tree.column("#0", width=64, minwidth=56, stretch=False)
        self.tree.column("name", width=220, minwidth=120)
        self.tree.column("saved", width=120, minwidth=108, stretch=False)
        self.tree.grid(row=0, column=0, sticky="nsew")
        scrollbar = ttk.Scrollbar(listing, orient="vertical", command=self.tree.yview)
        scrollbar.grid(row=0, column=1, sticky="ns")
        self.tree.configure(yscrollcommand=scrollbar.set)
        self.description = tk.StringVar(value="Select a protein to view its details.")
        description = ttk.Label(
            self.frame, textvariable=self.description, style="Surface.TLabel",
            wraplength=420, justify="left",
        )
        description.grid(row=3, column=0, sticky="ew", pady=16)
        self.primary = ttk.Button(
            self.frame, text="Analyze selection", style="Analyze.TButton",
            command=self.use_selected, state="disabled",
        )
        self.primary.grid(row=4, column=0, sticky="ew", pady=(0, 8))
        actions = ttk.Frame(self.frame, style="Surface.TFrame")
        actions.grid(row=5, column=0, sticky="ew")
        self.update = ttk.Button(
            actions, text="Refresh from RCSB", command=lambda: self.use_selected(refresh=True),
            state="disabled",
        )
        self.update.pack(side="left")
        self.source = ttk.Button(actions, text="RCSB entry", command=self.open_source,
                                 state="disabled")
        self.source.pack(side="right")
        location = ttk.Label(self.frame, text=f"Local files: {directory}",
                             style="Hint.TLabel", wraplength=420)
        location.grid(row=6, column=0, sticky="ew", pady=(16, 0))
        self.frame.bind("<Configure>", lambda event: (
            description.configure(wraplength=max(240, event.width - 48)),
            location.configure(wraplength=max(240, event.width - 48)),
        ))
        self.tree.bind("<<TreeviewSelect>>", self.show_selected)
        self.tree.bind("<Return>", self.use_selected)
        self.refresh()

    def refresh(self):
        """Refresh saved-state labels without resetting the user's selection."""
        selected = self.tree.selection()
        local = set(saved_ids(self.directory))
        self.tree.delete(*self.tree.get_children())
        examples = {example.pdb_id: example for example in EXAMPLES}
        for pdb_id in list(examples) + sorted(local - examples.keys()):
            name = examples[pdb_id].name if pdb_id in examples else "Saved structure"
            self.tree.insert("", "end", iid=pdb_id, text=pdb_id,
                             values=(name, "On this device" if pdb_id in local else "Not downloaded"))
        if selected:
            self.tree.selection_set(selected[0])

    def show_selected(self, event=None):
        selection = self.tree.selection()
        state = "normal" if selection else "disabled"
        for button in (self.primary, self.update, self.source):
            button.configure(state=state)
        if not selection:
            self.description.set("Select a protein to view its details.")
            return
        example = next((item for item in EXAMPLES if item.pdb_id == selection[0]), None)
        self.description.set(
            f"{example.organism}\n{example.description}" if example
            else f"{selection[0]} · Local copy. See its RCSB entry for context."
        )

    def use_selected(self, event=None, *, refresh=False):
        selection = self.tree.selection()
        if not selection:
            self.description.set("Select a protein from the list before analyzing.")
            self.tree.focus_set()
            return
        self.analyze(selection[0], refresh=refresh)

    def open_source(self):
        selection = self.tree.selection()
        if selection:
            webbrowser.open(f"https://www.rcsb.org/structure/{selection[0]}")
