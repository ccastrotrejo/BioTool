"""Liquid Glass-inspired light/dark tokens, with opaque native control surfaces."""

PALETTES = {"dark": {
    "background": "#101b24",
    "surface": "#182833",
    "elevated": "#213744",
    "text": "#edf4f7",
    "muted": "#a9bdc9",
    "accent": "#65c5df",
    "accent_hover": "#8bd8eb",
    "on_accent": "#10212c",
    "gold": "#dbc184",
    "border": "#3c5362",
    "glass_edge": "#72909f",
    "glass": "rgba(24, 40, 51, 0.88)",
    "shadow": "rgba(0, 8, 18, 0.28)",
}, "light": {
    "background": "#e4edf3",
    "surface": "#f4f8fb",
    "elevated": "#dce9f0",
    "text": "#172d3a",
    "muted": "#425e70",
    "accent": "#006383",
    "accent_hover": "#004f6b",
    "on_accent": "#ffffff",
    "gold": "#805d1d",
    "border": "#98aebc",
    "glass_edge": "#ffffff",
    "glass": "rgba(244, 248, 251, 0.88)",
    "shadow": "rgba(30, 65, 86, 0.14)",
}}
COLORS = PALETTES["dark"]
GROUP_PALETTES = {
    "dark": (COLORS["accent"], "#86c6a4", COLORS["gold"], "#9eafe8"),
    "light": ("#006383", "#237557", "#8f651d", "#5c64a6"),
}
GROUP_COLORS = GROUP_PALETTES["dark"]
SPACE = {"xs": 4, "sm": 8, "md": 16, "lg": 24, "xl": 32}
TYPE = {"caption": 11, "body": 12, "heading": 16, "title": 28}
CHART_FONT = "Avenir Next, Segoe UI, sans-serif"


def configure_desktop_styles(style, family: str, appearance: str = "dark") -> None:
    """Use explicit ttk styles so system appearance cannot mix light/dark colors."""
    # Aqua ignores custom backgrounds and padding; clam honors our complete palette.
    style.theme_use("clam")
    colors = PALETTES[appearance]
    style.configure(
        ".", background=colors["background"], foreground=colors["text"],
        font=(family, TYPE["body"]), bordercolor=colors["border"],
        lightcolor=colors["glass_edge"], darkcolor=colors["border"],
    )
    for name, surface in (("Bio", "background"), ("Surface", "surface")):
        style.configure(f"{name}.TFrame", background=colors[surface])
        style.configure(f"{name}.TLabel", background=colors[surface],
                        foreground=colors["text"], font=(family, TYPE["body"]))
    for name, surface, color, size, weight in (
        ("Muted", "background", "muted", "caption", "normal"),
        ("Hint", "surface", "muted", "caption", "normal"),
        ("Title", "background", "text", "title", "bold"),
        ("Eyebrow", "background", "gold", "caption", "bold"),
        ("Heading", "surface", "text", "heading", "bold"),
        ("Status", "elevated", "text", "caption", "normal"),
    ):
        style.configure(f"{name}.TLabel", background=colors[surface],
                        foreground=colors[color], font=(family, TYPE[size], weight))
    style.configure("Surface.TFrame", borderwidth=1, relief="raised",
                    lightcolor=colors["glass_edge"], darkcolor=colors["border"])
    style.configure("Status.TLabel", padding=(SPACE["md"], SPACE["sm"]))
    style.configure("TSeparator", background=colors["border"])
    style.configure(
        "Analyze.TButton", font=(family, TYPE["body"], "bold"),
        padding=(SPACE["md"], SPACE["md"]), background=colors["accent"],
        foreground=colors["on_accent"], bordercolor=colors["accent"],
        lightcolor=colors["glass_edge"], darkcolor=colors["accent"],
        focuscolor=colors["on_accent"], focusthickness=2, borderwidth=1,
    )
    style.map(
        "Analyze.TButton",
        background=[("disabled", colors["elevated"]), ("active", colors["accent_hover"])],
        foreground=[("disabled", colors["muted"])],
        lightcolor=[("active", colors["glass_edge"])],
        darkcolor=[("active", colors["accent_hover"])],
    )
    style.configure(
        "Pdb.TEntry", fieldbackground=colors["background"],
        foreground=colors["text"], bordercolor=colors["border"],
        insertcolor=colors["text"], padding=SPACE["md"],
        selectbackground=colors["accent"], selectforeground=colors["on_accent"],
        lightcolor=colors["border"], darkcolor=colors["glass_edge"],
    )
    style.map(
        "Pdb.TEntry", bordercolor=[("focus", colors["accent"])],
        lightcolor=[("focus", colors["accent"])],
        darkcolor=[("focus", colors["accent"])],
    )
    for name in ("TButton", "TRadiobutton", "TNotebook.Tab"):
        style.configure(name, padding=(12, 8), background=colors["surface"],
                        foreground=colors["text"], focuscolor=colors["accent"])
        style.map(name, background=[("selected", colors["elevated"]),
                                    ("active", colors["elevated"])],
                  foreground=[("disabled", colors["muted"])])
    style.configure("TNotebook", background=colors["background"], borderwidth=0)
    style.configure("Treeview", background=colors["surface"],
                    fieldbackground=colors["surface"], foreground=colors["text"],
                    rowheight=36, borderwidth=0)
    style.map("Treeview", background=[("selected", colors["accent"])],
              foreground=[("selected", colors["on_accent"])])
    style.configure("Treeview.Heading", background=colors["elevated"],
                    foreground=colors["text"], font=(family, TYPE["caption"], "bold"))
    style.map("Treeview.Heading", background=[("active", colors["elevated"])])
    style.configure("TCombobox", fieldbackground=colors["surface"],
                    background=colors["elevated"], foreground=colors["text"],
                    arrowcolor=colors["text"], padding=8)
    style.map("TCombobox", fieldbackground=[("readonly", colors["surface"])],
              foreground=[("readonly", colors["text"])])
    style.configure("Horizontal.TProgressbar", background=colors["accent"],
                    troughcolor=colors["elevated"], borderwidth=0)
