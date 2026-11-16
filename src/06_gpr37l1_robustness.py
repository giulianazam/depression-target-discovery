"""Step 8: how solid is the GPR37L1 drop in astrocytes?

Two checks on the astrocyte pseudobulk test from step 3:

  1. Leave-one-donor-out: refit DESeq2 once per donor, each time without that
     person, and record GPR37L1's fold change, p-value and FDR. If one or two
     people drive the result, it will vanish when they are dropped.
  2. Astrocyte share as a covariate: add each donor's % of astrocyte nuclei to
     the model (~ batch + astro_pct + condition). If GPR37L1 is low only because
     samples with fewer astrocytes are noisier, its depression effect should
     shrink towards zero once astro_pct is in the model.
  3. Nucleus quality as a covariate: also add the median number of genes
     detected per astrocyte nucleus for each donor. If MDD astrocyte nuclei are
     simply lower quality (less RNA captured), that, not depression, could make
     GPR37L1 look low.

Reads the raw GEO files directly (no scanpy needed).
"""
import gzip
import shutil
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

warnings.filterwarnings("ignore")
RAW = Path(sys.argv[1])          # folder with the GSE144136 files
OUT = Path(sys.argv[2])          # where to write tables
OUT.mkdir(parents=True, exist_ok=True)
GENE = "GPR37L1"
CELL = "Astros"
MIN_CELLS_PER_DONOR, MIN_DONORS_PER_GROUP = 10, 3   # same as step 3

# ---- load counts ----------------------------------------------------------------
mtx = RAW / "GSE144136_GeneBarcodeMatrix_Annotated.mtx.gz"
if not mtx.exists():  # rebuild from split parts, as in step 1
    with open(mtx, "wb") as dst:
        for p in sorted(RAW.glob("part_*")):
            with open(p, "rb") as src:
                shutil.copyfileobj(src, dst)
genes = pd.read_csv(RAW / "GSE144136_GeneNames.csv.gz", index_col=0)["x"].values
cells = pd.read_csv(RAW / "GSE144136_CellNames.csv.gz", index_col=0)["x"].values
with gzip.open(mtx, "rb") as f:
    m = sp.csr_matrix(scipy.io.mmread(f))
if m.shape == (len(genes), len(cells)):
    m = m.T.tocsr()
meta = pd.Series(cells).str.extract(
    r"^(?P<cluster>.+)\.(?P<donor>\d+)_(?P<condition>[A-Za-z]+)_(?P<batch>B\d+)_")
meta["donor"] = "S" + meta["donor"].str.zfill(2)
meta["condition"] = meta["condition"].replace({"Suicide": "MDD"})
meta["type"] = meta["cluster"].str.split("_").str[0]

# % astrocyte nuclei per donor, over all nuclei (as in figure 2)
astro_pct = meta.groupby("donor")["type"].apply(lambda t: (t == CELL).mean() * 100)

# ---- astrocyte pseudobulk -----------------------------------------------------
idx = np.where(meta["type"] == CELL)[0]
sub = meta.iloc[idx]
rows, info = [], []
for donor, g in sub.groupby("donor"):
    if len(g) < MIN_CELLS_PER_DONOR:
        continue
    rows.append(np.asarray(m[g.index.values].sum(axis=0)).ravel())  # index = row in the matrix
    genes_per_nucleus = np.diff(m[g.index.values].indptr)  # non-zero genes in each nucleus
    info.append({"donor": donor, "condition": g["condition"].iloc[0], "batch": g["batch"].iloc[0],
                 "n_astro_nuclei": len(g), "median_genes_per_nucleus": float(np.median(genes_per_nucleus))})
info = pd.DataFrame(info).set_index("donor")
info["astro_pct"] = astro_pct.reindex(info.index)
counts_all = pd.DataFrame(np.vstack(rows), index=info.index, columns=pd.Index(genes).astype(str))
counts_all = counts_all.loc[:, ~counts_all.columns.duplicated()].astype(int)
print(f"astrocyte pseudobulk: {len(info)} donors", info["condition"].value_counts().to_dict())
print(info.groupby("condition")[["n_astro_nuclei", "median_genes_per_nucleus", "astro_pct"]].median())


def fit(meta_in, design):
    meta_in = meta_in.copy()
    meta_in["batch"] = meta_in["batch"].where(meta_in["batch"].map(meta_in["batch"].value_counts()) > 1, "other")
    c = counts_all.loc[meta_in.index]
    c = c.loc[:, (c >= 10).sum(axis=0) >= MIN_DONORS_PER_GROUP]
    dds = DeseqDataSet(counts=c, metadata=meta_in, design=design, quiet=True)
    dds.deseq2()
    st = DeseqStats(dds, contrast=["condition", "MDD", "Control"], quiet=True)
    st.summary()
    return st.results_df


# ---- 0. reproduce the original result ------------------------------------------
base = fit(info, "~batch + condition")
r = base.loc[GENE]
print(f"original: log2FC {r.log2FoldChange:.3f}  p {r.pvalue:.2e}  FDR {r.padj:.3f}")

# ---- 1. leave one donor out ----------------------------------------------------
loo = []
for d in info.index:
    res = fit(info.drop(index=d), "~batch + condition").loc[GENE]
    loo.append({"dropped_donor": d, "dropped_condition": info.at[d, "condition"],
                "log2FoldChange": res.log2FoldChange, "pvalue": res.pvalue, "padj": res.padj})
    print(d, f"{res.log2FoldChange:.3f} {res.pvalue:.2e} {res.padj:.3f}", flush=True)
loo = pd.DataFrame(loo)
loo.to_csv(OUT / "gpr37l1_leave_one_donor_out.csv", index=False)

# ---- 2. astrocyte share as covariate ------------------------------------------
cov = info.copy()
cov["astro_pct_z"] = (cov["astro_pct"] - cov["astro_pct"].mean()) / cov["astro_pct"].std()
adj = fit(cov, "~batch + astro_pct_z + condition")
a = adj.loc[GENE]
for col in ["median_genes_per_nucleus"]:
    cov[col + "_z"] = (cov[col] - cov[col].mean()) / cov[col].std()
q = fit(cov, "~batch + astro_pct_z + median_genes_per_nucleus_z + condition").loc[GENE]
summary = pd.DataFrame([
    {"model": "batch + condition (original)", "log2FoldChange": r.log2FoldChange, "pvalue": r.pvalue, "padj": r.padj},
    {"model": "batch + astrocyte % + condition", "log2FoldChange": a.log2FoldChange, "pvalue": a.pvalue, "padj": a.padj},
    {"model": "batch + astrocyte % + genes per nucleus + condition", "log2FoldChange": q.log2FoldChange,
     "pvalue": q.pvalue, "padj": q.padj},
])
summary.to_csv(OUT / "gpr37l1_astro_pct_covariate.csv", index=False)
info.to_csv(OUT / "astrocyte_pseudobulk_donors.csv")
print(summary.to_string(index=False))
print("LOO: padj<0.05 in", (loo.padj < 0.05).sum(), "of", len(loo), "; p<0.05 in", (loo.pvalue < 0.05).sum(),
      "; log2FC range", loo.log2FoldChange.min().round(2), "to", loo.log2FoldChange.max().round(2))
