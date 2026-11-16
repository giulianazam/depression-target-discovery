"""Step 6: figures for the README."""
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "results/figures"
FIG.mkdir(parents=True, exist_ok=True)

# Validated categorical palette (colorblind-checked); text stays in neutral ink.
BLUE, ORANGE, AQUA, YELLOW, MAGENTA, GREEN, VIOLET, RED = (
    "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948")
GROUP = {"Control": BLUE, "MDD": ORANGE}
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e4e3df"
plt.rcParams.update({
    "font.size": 10, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": MUTED, "ytick.color": MUTED, "axes.spines.top": False,
    "axes.spines.right": False, "figure.dpi": 150, "savefig.bbox": "tight",
    "figure.facecolor": "white",
})

cells = ad.read_h5ad(ROOT / "data/processed/cells.h5ad", backed="r")
obs = cells.obs.copy()
umap = cells.obsm["X_umap"]

# 1. Cell map: what cell types are in the data (direct labels, so color is not the only cue)
order = ["Excitatory neurons", "Inhibitory neurons", "Astrocytes", "Oligodendrocytes",
         "OPCs", "Microglia/Macrophages", "Endothelial", "Mixed/unclear"]
colors = dict(zip(order, [BLUE, ORANGE, AQUA, YELLOW, MAGENTA, GREEN, VIOLET, "#b8b7b2"]))
fig, axes = plt.subplots(1, 2, figsize=(11, 5))
rng = np.random.default_rng(0)
idx = rng.permutation(len(obs))
ax = axes[0]
ax.scatter(umap[idx, 0], umap[idx, 1], s=0.3, c=obs["cell_type"].iloc[idx].map(colors).values,
           rasterized=True, linewidths=0)
for ct in order:
    m = (obs["cell_type"] == ct).values
    if m.sum() == 0:
        continue
    x, y = np.median(umap[m], axis=0)
    ax.text(x, y, ct.replace("Microglia/Macrophages", "Microglia"), fontsize=8.5, ha="center",
            va="center", color=INK, bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.8))
ax.set_title("78,886 brain nuclei, by cell type", loc="left")
ax = axes[1]
for cond in ["Control", "MDD"]:
    m = (obs["condition"].iloc[idx] == cond).values
    ax.scatter(umap[idx][m, 0], umap[idx][m, 1], s=0.3, c=GROUP[cond], label=cond,
               rasterized=True, linewidths=0, alpha=0.6)
ax.legend(markerscale=12, frameon=False, loc="upper right")
ax.set_title("Same cells, by diagnosis", loc="left")
for ax in axes:
    ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")
fig.savefig(FIG / "1_cell_map.png")

# 2. Astrocyte share per donor
# share among labelled nuclei: "Mixed/unclear" is dropped here, as in every other step
labelled = obs[obs["cell_type"] != "Mixed/unclear"]
prop = pd.crosstab(labelled["donor"], labelled["cell_type"].astype(str), normalize="index") * 100
cond = obs.drop_duplicates("donor").set_index("donor")["condition"]
fig, ax = plt.subplots(figsize=(4.2, 4.2))
vals = {g: prop.loc[cond[cond == g].index, "Astrocytes"] for g in ["Control", "MDD"]}
for i, g in enumerate(vals):
    v = vals[g]
    ax.scatter(np.full(len(v), i) + rng.uniform(-0.12, 0.12, len(v)), v, s=36, color=GROUP[g],
               edgecolor="white", linewidth=1, zorder=3)
    ax.hlines(v.median(), i - 0.25, i + 0.25, color=INK, lw=2, zorder=4)
p = mannwhitneyu(vals["Control"], vals["MDD"]).pvalue
ax.set_xticks([0, 1], [f"Control\n(n={len(vals['Control'])})", f"MDD\n(n={len(vals['MDD'])})"])
ax.set_xlim(-0.6, 1.6)
ax.set_ylabel("Astrocytes, % of each donor's labelled nuclei")
ax.yaxis.grid(True, color=GRID); ax.set_axisbelow(True)
ax.set_title(f"Fewer astrocyte nuclei in depression\n(Mann-Whitney p = {p:.3f}; line = median)", loc="left")
fig.savefig(FIG / "2_astrocyte_share.png")

# 3. How many genes change in each cell type
de = pd.read_csv(ROOT / "results/tables/de_mdd_vs_control_by_celltype.csv")
sig = de[de["padj"] < 0.1]
counts = pd.DataFrame({
    "up": sig[sig["log2FoldChange"] > 0].groupby("cell_type").size(),
    "down": sig[sig["log2FoldChange"] < 0].groupby("cell_type").size(),
}).reindex(de["cell_type"].unique()).fillna(0)
counts = counts.loc[(counts.sum(axis=1)).sort_values().index]
donors = de.groupby("cell_type")["n_donors"].first()
fig, ax = plt.subplots(figsize=(7, 3.6))
y = np.arange(len(counts))
ax.barh(y, -counts["down"], color=BLUE, height=0.6, label="Lower in MDD")
ax.barh(y, counts["up"], color=ORANGE, height=0.6, label="Higher in MDD")
ax.set_yticks(y, [f"{c} ({donors[c]} donors{', not ranked' if donors[c] < 25 else ''})" for c in counts.index])
ax.axvline(0, color=MUTED, lw=1)
for yi, (u, d) in enumerate(zip(counts["up"], counts["down"])):
    ax.text(u + 8, yi, f"{int(u)}", va="center", fontsize=8, color=MUTED)
    ax.text(-d - 8, yi, f"{int(d)}", va="center", ha="right", fontsize=8, color=MUTED)
lim = counts.values.max() * 1.25
ax.set_xlim(-lim, lim)
ax.set_xticks(ax.get_xticks(), [f"{abs(int(t))}" for t in ax.get_xticks()])
ax.set_xlabel("Genes changed (FDR < 0.1, corrected within each cell type)")
ax.legend(frameon=False, loc="lower right")
ax.set_title("Genes changed in depression, by cell type", loc="left")
fig.savefig(FIG / "3_genes_changed.png")

# 4. Target shortlist with what drives each score
sl = pd.read_csv(ROOT / "results/tables/target_shortlist.csv").head(12).iloc[::-1]
parts = [("evidence", 0.35, "Evidence", BLUE), ("druggability", 0.30, "Druggability", ORANGE),
         ("brain_focus", 0.20, "Brain focus", AQUA), ("specificity", 0.15, "Cell specificity", VIOLET)]
fig, ax = plt.subplots(figsize=(8, 5))
left = np.zeros(len(sl))
for col, w, label, c in parts:
    width = sl[col].values * w
    ax.barh(np.arange(len(sl)), width, left=left, color=c, height=0.62, label=label,
            edgecolor="white", linewidth=1.5)
    left += width
labels = [f"{g}  ({ct.replace('Excitatory neurons', 'Exc. neurons')}, {d.split()[0]})"
          for g, ct, d in zip(sl["gene"], sl["cell_type"], sl["direction"])]
ax.set_yticks(np.arange(len(sl)), labels)
for yi, v in enumerate(left):
    ax.text(v + 0.01, yi, f"{v:.2f}", va="center", fontsize=8, color=MUTED)
ax.set_xlim(0, 1)
ax.set_xlabel("Priority score (0 to 1)")
ax.legend(frameon=False, ncol=4, loc="upper center", bbox_to_anchor=(0.45, -0.12))
ax.set_title("Top candidate drug targets", loc="left")
fig.savefig(FIG / "4_target_shortlist.png")
print("figures saved to", FIG)
