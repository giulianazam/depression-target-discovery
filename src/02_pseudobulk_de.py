"""Step 3: depressed vs. control differential expression, one cell type at a time.

Single cells from the same person are not independent, so we sum counts per
donor and cell type ("pseudobulk") and test with DESeq2. This is the approach
recommended for multi-donor single-cell studies (Squair et al. 2021).
"""
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[1]
MIN_CELLS_PER_DONOR = 10  # skip donor/cell-type pairs with too few cells
MIN_DONORS_PER_GROUP = 3

adata = sc.read_h5ad(ROOT / "data/processed/cells.h5ad")
# Library prep batches: show how donors are spread across batch and diagnosis, and adjust for batch.
print(adata.obs.drop_duplicates("donor").groupby(["batch", "condition"], observed=True).size().unstack(fill_value=0))
out = ROOT / "results/tables"
out.mkdir(parents=True, exist_ok=True)

all_results = []
for cell_type in sorted(adata.obs["cell_type"].unique()):
    if cell_type == "Mixed/unclear":
        continue
    sub = adata[adata.obs["cell_type"] == cell_type]
    groups = sub.obs.groupby("donor", observed=True)
    counts, meta = [], []
    for donor, idx in groups.indices.items():
        if len(idx) < MIN_CELLS_PER_DONOR:
            continue
        counts.append(np.asarray(sub.layers["counts"][idx].sum(axis=0)).ravel())
        meta.append({"donor": donor, "condition": sub.obs["condition"].iloc[idx[0]],
                     "batch": sub.obs["batch"].iloc[idx[0]]})
    meta = pd.DataFrame(meta).set_index("donor")
    # a batch with a single donor left cannot be estimated; pool those into "other"
    meta["batch"] = meta["batch"].where(meta["batch"].map(meta["batch"].value_counts()) > 1, "other")
    if meta["condition"].value_counts().min() < MIN_DONORS_PER_GROUP or meta["condition"].nunique() < 2:
        print(f"{cell_type}: not enough donors, skipped")
        continue
    counts = pd.DataFrame(counts, index=meta.index, columns=sub.var_names).astype(int)
    counts = counts.loc[:, (counts >= 10).sum(axis=0) >= MIN_DONORS_PER_GROUP]  # drop barely-expressed genes

    dds = DeseqDataSet(counts=counts, metadata=meta, design="~batch + condition", quiet=True)
    dds.deseq2()
    stats = DeseqStats(dds, contrast=["condition", "MDD", "Control"], quiet=True)
    stats.summary()
    res = stats.results_df.assign(cell_type=cell_type, n_donors=len(meta))
    all_results.append(res.rename_axis("gene").reset_index())
    n_sig = (res["padj"] < 0.1).sum()
    print(f"{cell_type}: {len(meta)} donors, {counts.shape[1]} genes tested, {n_sig} significant (FDR < 0.1)")

results = pd.concat(all_results, ignore_index=True)
results.to_csv(out / "de_mdd_vs_control_by_celltype.csv", index=False)
print(f"saved {out / 'de_mdd_vs_control_by_celltype.csv'}")
