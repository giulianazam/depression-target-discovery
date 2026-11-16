"""Vector (SVG) versions of figures 3 and 4 for the poster, rebuilt from the result tables.

Same plots as src/04_figures.py, without in-figure titles (the poster panels carry them).
Figures 1 and 2 need data/processed/cells.h5ad, so the poster uses the existing PNGs.
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = Path(__file__).resolve().parent / "assets"
OUT.mkdir(exist_ok=True)

BLUE, ORANGE, AQUA, VIOLET = "#2a78d6", "#eb6834", "#1baf7a", "#6250d6"
INK, MUTED = "#0b0b0b", "#52514e"
plt.rcParams.update({
    "font.family": "Inter", "font.size": 11, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
    "text.color": INK, "xtick.color": MUTED, "ytick.color": MUTED, "axes.spines.top": False,
    "axes.spines.right": False, "savefig.bbox": "tight", "svg.fonttype": "none",
    "figure.facecolor": "none", "axes.facecolor": "none",
})

de = pd.read_csv(ROOT / "results/tables/de_mdd_vs_control_by_celltype.csv")
sig = de[de["padj"] < 0.1]
counts = pd.DataFrame({
    "up": sig[sig["log2FoldChange"] > 0].groupby("cell_type").size(),
    "down": sig[sig["log2FoldChange"] < 0].groupby("cell_type").size(),
}).reindex(de["cell_type"].unique()).fillna(0)
counts = counts.loc[(counts.sum(axis=1)).sort_values().index]
donors = de.groupby("cell_type")["n_donors"].first()
fig, ax = plt.subplots(figsize=(7.5, 3.6))
y = np.arange(len(counts))
ax.barh(y, -counts["down"], color=BLUE, height=0.6, label="Lower in MDD")
ax.barh(y, counts["up"], color=ORANGE, height=0.6, label="Higher in MDD")
ax.set_yticks(y, [f"{c} ({donors[c]} donors{', not ranked' if donors[c] < 25 else ''})" for c in counts.index])
ax.axvline(0, color=MUTED, lw=1)
for yi, (u, d) in enumerate(zip(counts["up"], counts["down"])):
    ax.text(u + 8, yi, f"{int(u)}", va="center", fontsize=10, color=MUTED)
    ax.text(-d - 8, yi, f"{int(d)}", va="center", ha="right", fontsize=10, color=MUTED)
lim = counts.values.max() * 1.25
ax.set_xlim(-lim, lim)
ax.set_xticks(ax.get_xticks(), [f"{abs(int(t))}" for t in ax.get_xticks()])
ax.set_xlim(-lim, lim)
ax.set_xlabel("Genes changed (FDR < 0.1, corrected within each cell type)")
ax.legend(frameon=False, loc="lower right")
fig.savefig(OUT / "3_genes_changed.svg")

sl = pd.read_csv(ROOT / "results/tables/target_shortlist.csv").head(10).iloc[::-1]
parts = [("evidence", 0.35, "Evidence", BLUE), ("druggability", 0.30, "Druggability", ORANGE),
         ("brain_focus", 0.20, "Brain focus", AQUA), ("specificity", 0.15, "Cell specificity", VIOLET)]
fig, ax = plt.subplots(figsize=(8, 4.8))
left = np.zeros(len(sl))
for col, w, label, c in parts:
    width = sl[col].values * w
    ax.barh(np.arange(len(sl)), width, left=left, color=c, height=0.62, label=label,
            edgecolor="white", linewidth=1.5)
    left += width
labels = [f"{g}  ({ct.replace('Excitatory neurons', 'Exc. neurons')}, {d.split()[0]})"
          for g, ct, d in zip(sl["gene"], sl["cell_type"], sl["direction"])]
ax.set_yticks(np.arange(len(sl)), labels)
ax.get_yticklabels()[-1].set_fontweight("bold")
for yi, v in enumerate(left):
    ax.text(v + 0.01, yi, f"{v:.2f}", va="center", fontsize=10, color=MUTED)
ax.set_xlim(0, 1)
ax.set_xlabel("Priority score (0 to 1)")
ax.legend(frameon=False, ncol=4, loc="upper center", bbox_to_anchor=(0.4, -0.13))
fig.savefig(OUT / "4_target_shortlist.svg")
print("saved to", OUT)
