"""Frosted report surfaces with high-contrast light/dark fallbacks."""

import json

from .theme import GROUP_PALETTES, PALETTES


STYLE = "\n".join(
    f':root[data-theme="{name}"] {{ color-scheme: {name}; ' + "; ".join(
        f"--{key.replace('_', '-')}: {value}" for key, value in colors.items()
    ) + "; }" for name, colors in PALETTES.items()
) + """
* { box-sizing: border-box; }
body { margin: 0; background: var(--background); color: var(--text);
  font: 16px/1.6 "Avenir Next", "Segoe UI", sans-serif; }
body::before { content: ""; position: fixed; inset: 0; z-index: -1;
  background: radial-gradient(ellipse at 0 10%, var(--elevated), transparent 65%),
  linear-gradient(150deg, transparent 35%, var(--border)); opacity: 0.5; }
main { max-width: 1120px; margin: auto; padding: 32px 24px; }
header { overflow-wrap: anywhere; }
.toolbar { display: flex; justify-content: space-between; align-items: center;
  gap: 16px; flex-wrap: wrap; margin-bottom: 24px; }
.brand { color: var(--accent); font-size: 18px; font-weight: 700; margin: 0; }
.theme-picker { display: flex; gap: 8px; align-items: center; color: var(--muted);
  font-size: 14px; }
select { color: var(--text); background: var(--surface); border: 1px solid var(--border);
  border-radius: 12px; padding: 8px 12px; font: inherit; min-height: 44px; }
nav { display: flex; gap: 8px 24px; flex-wrap: wrap; margin: 24px 0;
  padding: 8px 24px; }
nav a { padding: 12px 0; text-decoration: none; border-bottom: 2px solid transparent; }
a { color: var(--accent); text-underline-offset: 4px; }
a:hover { color: var(--accent-hover); }
a[aria-current] { font-weight: 700; color: var(--text); border-color: var(--accent); }
a:focus-visible, summary:focus-visible, select:focus-visible {
  outline: 3px solid var(--accent); outline-offset: 4px; }
h1 { font-size: clamp(28px, 4vw, 40px); letter-spacing: -0.03em;
  line-height: 1.2; margin: 0 0 16px; text-wrap: balance; }
h2 { font-size: 20px; margin: 0 0 16px; }
p { margin: 0 0 16px; }
.muted { color: var(--muted); max-width: 75ch; }
.panel, nav { background: var(--surface); border: 1px solid var(--border);
  border-top-color: var(--glass-edge); border-radius: 24px;
  box-shadow: 0 12px 32px var(--shadow), inset 0 1px 0 var(--glass-edge); }
.panel { padding: 24px; margin: 24px 0; min-width: 0; }
.chart { padding: 8px; overflow: hidden; }
@supports ((backdrop-filter: blur(20px)) or (-webkit-backdrop-filter: blur(20px))) {
  .panel, nav { background: var(--glass); backdrop-filter: blur(20px) saturate(130%);
    -webkit-backdrop-filter: blur(20px) saturate(130%); }
}
summary { cursor: pointer; font-weight: 600; padding: 12px 0; color: var(--accent); }
pre { white-space: pre-wrap; overflow-wrap: anywhere; font-size: 14px;
  background: var(--background); padding: 16px; border: 1px solid var(--border);
  border-radius: 12px; }
table { width: 100%; border-collapse: collapse; margin-top: 16px; }
caption { text-align: left; font-weight: 600; color: var(--gold); }
th, td { padding: 8px; text-align: left; border-bottom: 1px solid var(--border); }
td:last-child, th:last-child { text-align: right; font-variant-numeric: tabular-nums; }
footer { color: var(--muted); font-size: 14px; }
.js-plotly-plot .plotly .modebar-btn:focus-visible { outline: 2px solid var(--accent); }
#theme-status:empty { display: none; }
@media (max-width: 600px) {
  main { padding: 24px 12px; } .panel { padding: 16px; } .chart { padding: 0; }
  nav { padding: 8px 16px; }
}
@media (prefers-reduced-transparency: reduce), (prefers-contrast: more) {
  .panel, nav { background: var(--surface); backdrop-filter: none;
    -webkit-backdrop-filter: none; box-shadow: none; }
  body::before { display: none; }
}
@media print {
  .theme-picker, nav, .modebar { display: none; }
  .panel { break-inside: avoid; box-shadow: none; }
}
"""

THEME_SCRIPT = (
    "<script>\nconst palettes = " + json.dumps(PALETTES) + ";\nconst groupPalettes = "
    + json.dumps(GROUP_PALETTES) + ";\n" + """
const themePicker = document.getElementById("appearance");
themePicker.addEventListener("change", async () => {
  const theme = themePicker.value;
  const colors = palettes[theme];
  const groups = groupPalettes[theme];
  themePicker.disabled = true;
  document.documentElement.dataset.theme = theme;
  document.getElementById("theme-status").textContent = "";
  try {
    for (const plot of document.querySelectorAll(".js-plotly-plot")) {
      const layout = {
        paper_bgcolor: colors.surface, plot_bgcolor: colors.surface,
        "font.color": colors.text,
        "hoverlabel.bgcolor": colors.elevated, "hoverlabel.font.color": colors.text,
        "modebar.bgcolor": colors.surface, "modebar.color": colors.muted,
        "modebar.activecolor": colors.accent
      };
      if (plot.data.some(trace => trace.type === "scatter3d")) {
        layout["scene.bgcolor"] = colors.surface;
        layout["xaxis.color"] = colors.text;
        layout["yaxis.color"] = colors.text;
        layout["yaxis.gridcolor"] = colors.border;
        layout["yaxis.zerolinecolor"] = colors.border;
        for (const axis of ["xaxis", "yaxis", "zaxis"]) {
          layout["scene." + axis + ".color"] = colors.text;
          layout["scene." + axis + ".gridcolor"] = colors.border;
        }
      }
      // Plotly templates give annotations explicit colors; update those as well.
      for (let i = 0; i < (plot.layout.annotations || []).length; i++) {
        layout["annotations[" + i + "].font.color"] = colors.text;
      }
      await Plotly.relayout(plot, layout);
      let bar = 0;
      for (let i = 0; i < plot.data.length; i++) {
        const trace = plot.data[i];
        if (trace.type === "bar") {
          await Plotly.restyle(plot, {"marker.color": groups[bar++]}, [i]);
        } else if (trace.type === "pie") {
          await Plotly.restyle(plot, {
            "marker.colors": [groups], "insidetextfont.color": colors.on_accent
          }, [i]);
        } else if (trace.type === "scatter3d") {
          await Plotly.restyle(plot, {"marker.colorscale": [
            [[0, colors.muted], [0.5, colors.accent], [1, colors.gold]]
          ]}, [i]);
        }
      }
    }
  } catch (error) {
    document.getElementById("theme-status").textContent =
      "No se pudo actualizar la gráfica. Recarga el informe para restaurar su apariencia.";
    console.error("BioTool theme update failed", error);
  } finally {
    themePicker.disabled = false;
  }
});
</script>
"""
)
