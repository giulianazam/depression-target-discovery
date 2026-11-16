"""Step 2: load the Douglas brain bank single-nucleus data and label cells.

Dataset: GSE144136 (Nagy et al., Nature Neuroscience 2020). ~79,000 nuclei from
the dorsolateral prefrontal cortex of 17 men with major depression who died by
suicide and 17 matched male controls, from the Douglas-Bell Canada Brain Bank.

Cell names look like "Astros_2.12_Control_B3_AAAGATGCATGCAATC", i.e.
<cluster>.<subject>_<condition>_<batch>_<barcode>. The authors already
clustered and labelled the cells, so we reuse their labels.
"""
import glob
import gzip
import shutil
from pathlib import Path

import pandas as pd
import scanpy as sc
import scipy.io
import scipy.sparse as sp

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data/raw"
OUT = ROOT / "data/processed"
OUT.mkdir(parents=True, exist_ok=True)

mtx_path = RAW / "GSE144136_GeneBarcodeMatrix_Annotated.mtx.gz"
if not mtx_path.exists():  # the file can be split into part_aa, part_ab, ... to move it around
    parts = sorted(glob.glob(str(RAW / "part_*")))
    with open(mtx_path, "wb") as dst:
        for p in parts:
            with open(p, "rb") as src:
                shutil.copyfileobj(src, dst)

genes = pd.read_csv(RAW / "GSE144136_GeneNames.csv.gz", index_col=0)["x"].values
cells = pd.read_csv(RAW / "GSE144136_CellNames.csv.gz", index_col=0)["x"].values
with gzip.open(mtx_path, "rb") as f:
    m = sp.csr_matrix(scipy.io.mmread(f))
if m.shape == (len(genes), len(cells)):
    m = m.T.tocsr()  # scanpy wants cells as rows
assert m.shape == (len(cells), len(genes)), m.shape

adata = sc.AnnData(X=m)
adata.obs_names = cells
adata.var_names = genes
adata.var_names_make_unique()

meta = adata.obs_names.to_series().str.extract(
    r"^(?P<cluster>.+)\.(?P<donor>\d+)_(?P<condition>[A-Za-z]+)_(?P<batch>B\d+)_"
)
meta["donor"] = "S" + meta["donor"].str.zfill(2)
# "Suicide" cases all had major depressive disorder (MDD); we call them MDD.
meta["condition"] = meta["condition"].replace({"Suicide": "MDD"})
broad = {"Ex": "Excitatory neurons", "Inhib": "Inhibitory neurons", "Astros": "Astrocytes",
         "Oligos": "Oligodendrocytes", "OPCs": "OPCs", "Micro/Macro": "Microglia/Macrophages",
         "Endo": "Endothelial", "Mix": "Mixed/unclear"}
meta["cell_type"] = meta["cluster"].str.split("_").str[0].map(broad)
adata.obs = meta.astype("category")

adata.layers["counts"] = adata.X.copy()
sc.pp.calculate_qc_metrics(adata, inplace=True, percent_top=None)
print(adata)
print(adata.obs.groupby(["cell_type", "condition"], observed=True).size().unstack())

# Normalise and embed for a cell map (UMAP) figure.
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="batch")
sc.pp.pca(adata, n_comps=30, mask_var="highly_variable")
sc.pp.neighbors(adata)
sc.tl.umap(adata)

adata.write_h5ad(OUT / "cells.h5ad", compression="gzip")
print(f"saved {OUT / 'cells.h5ad'}")
