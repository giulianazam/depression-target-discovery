"""Step 5: turn the list of changed genes into a ranked shortlist of drug targets.

Each gene that changes in depression in a given cell type gets four scores (0 to 1):

  evidence      how strong and consistent the change is (FDR and fold change)
  druggability  whether it belongs to a protein family drugs are known to work on
                (Hopkins & Groom 2002), using HGNC gene families
  brain focus   how much more it is expressed in brain than in other organs (GTEx),
                a rough proxy for fewer side effects elsewhere
  specificity   whether the gene really belongs to that cell type, rather than
                leaking in from another one (ambient RNA)

priority = 0.35*evidence + 0.30*druggability + 0.20*brain_focus + 0.15*specificity
The weights are a judgement call, and they are written here so anyone can change them.
"""
import re
import urllib.request
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
EXT = ROOT / "data/external"
EXT.mkdir(parents=True, exist_ok=True)
SOURCES = {
    "hgnc_complete_set.txt":
        "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt",
    "gtex_v8_median_tpm.gct.gz":
        "https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/"
        "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz",
}
for name, url in SOURCES.items():
    if not (EXT / name).exists():
        print(f"downloading {name}")
        urllib.request.urlretrieve(url, EXT / name)

FDR = 0.1
MIN_DONORS = 25  # cell types seen in fewer donors are too noisy to rank (drops endothelial cells)

# --- differential expression hits --------------------------------------------
de = pd.read_csv(ROOT / "results/tables/de_mdd_vs_control_by_celltype.csv")
hits = de[(de["padj"] < FDR) & (de["n_donors"] >= MIN_DONORS)].copy()

# --- druggability from HGNC gene families --------------------------------------
hgnc = pd.read_csv(EXT / "hgnc_complete_set.txt", sep="\t", low_memory=False)
families = {}
for col in ["symbol", "prev_symbol", "alias_symbol"]:  # the dataset uses some older gene names
    for syms, grp, name in hgnc[[col, "gene_group", "name"]].itertuples(index=False):
        for s in str(syms).strip('"').split("|"):
            families.setdefault(s, (str(grp), str(name)))

TIERS = [  # (score, label, pattern on HGNC gene group)
    (1.0, "GPCR", r"G protein-coupled receptor|Adrenoceptors|Dopamine receptors|"
                  r"5-hydroxytryptamine receptors, G protein|Metabotropic glutamate|Adenosine receptors"),
    (1.0, "Ion channel", r"channel|receptor subunits|Glutamate ionotropic|ionotropic"),
    (1.0, "Nuclear receptor", r"Nuclear hormone receptors"),
    (0.7, "Transporter", r"Solute carrier|ATP binding cassette"),
    (0.7, "Enzyme", r"phosphodiesterase|peptidase|protease|hydrolase|dehydrogenase|"
                    r"transferase|Cytochrome P450|oxidase|reductase|lipase|synthase|carbonic anhydrase|"
                    r"Glutamate-ammonia ligase|phosphatase"),
    (0.5, "Secreted / surface (antibody)", r"Receptor ligands|Interleukin|Growth factor|"
                    r"Fibroblast growth factor|Tumor necrosis factor|Complement|Apolipoproteins|"
                    r"Immunoglobulin like domain|CD molecules|Gap junction"),
]

# Kinases are judged on the protein's own name, because scaffold and regulatory
# proteins (e.g. "kinase interacting protein") sit in kinase-related families.
KINASE = r"kinase(?!.*(interacting|regulatory|activator|anchoring|inhibitor|associated|binding))"
NOT_KINASE_FAMILY = r"Membrane associated guanylate kinases|Pseudokinases"


def druggability(gene):
    grp, name = families.get(gene, ("", ""))
    text = f"{grp} | {name}"
    if re.search(KINASE, name, flags=re.I) and not re.search(NOT_KINASE_FAMILY, grp):
        return 1.0, "Kinase", grp
    helper = re.search(r"regulator|interacting|inhibitor|activator|binding protein", name, flags=re.I)
    for score, label, pat in TIERS:
        if helper and score > 0.5:  # e.g. "phosphatase and actin regulator" is not an enzyme
            continue
        if re.search(pat, text, flags=re.I):
            return score, label, grp
    return 0.0, "", grp

hits[["druggability", "drug_class", "gene_family"]] = [druggability(g) for g in hits["gene"]]

# --- brain focus from GTEx -------------------------------------------------------
gtex = pd.read_csv(EXT / "gtex_v8_median_tpm.gct.gz", sep="\t", skiprows=2)
gtex = gtex.drop(columns="Name").groupby("Description").max()
brain = gtex[[c for c in gtex.columns if c.startswith("Brain")]].median(axis=1)
other = gtex[[c for c in gtex.columns if not c.startswith(("Brain", "Cells"))]].median(axis=1)
hits["brain_tpm"] = hits["gene"].map(brain)
hits["other_organs_tpm"] = hits["gene"].map(other)
ratio = np.log2((hits["brain_tpm"] + 1) / (hits["other_organs_tpm"] + 1))
hits["brain_focus"] = ratio.clip(0, 5).div(5).fillna(0)

# --- cell-type specificity (guards against ambient RNA) ---------------------------
cells = ad.read_h5ad(ROOT / "data/processed/cells.h5ad")
cells = cells[cells.obs["cell_type"] != "Mixed/unclear"]
genes = hits["gene"].unique()
mean_expr = pd.DataFrame(
    {ct: np.asarray(cells[cells.obs["cell_type"] == ct, genes].X.mean(axis=0)).ravel()
     for ct in cells.obs["cell_type"].cat.categories if ct != "Mixed/unclear"},
    index=genes)
hits["home_cell_type"] = hits["gene"].map(mean_expr.idxmax(axis=1))
own = np.array([mean_expr.at[g, ct] for g, ct in zip(hits["gene"], hits["cell_type"])])
top = hits["gene"].map(mean_expr.max(axis=1)).values
hits["specificity"] = np.clip(own / np.maximum(top, 1e-9), 0, 1)
# a gene mostly made by another cell type, changing here, is likely contamination
hits["possible_contamination"] = hits["specificity"] < 0.2

# --- combine -----------------------------------------------------------------------
strength = -np.log10(hits["padj"].clip(lower=1e-20)) / 20 * 0.5 + hits["log2FoldChange"].abs().clip(0, 2) / 2 * 0.5
hits["evidence"] = strength.clip(0, 1)
hits["priority"] = (0.35 * hits["evidence"] + 0.30 * hits["druggability"]
                    + 0.20 * hits["brain_focus"] + 0.15 * hits["specificity"])
hits["direction"] = np.where(hits["log2FoldChange"] > 0, "up in MDD", "down in MDD")

cols = ["gene", "cell_type", "direction", "log2FoldChange", "padj", "drug_class", "gene_family",
        "brain_tpm", "other_organs_tpm", "home_cell_type", "possible_contamination",
        "evidence", "druggability", "brain_focus", "specificity", "priority"]
ranked = hits.sort_values("priority", ascending=False)[cols]
ranked.to_csv(ROOT / "results/tables/all_hits_scored.csv", index=False)

shortlist = ranked[(ranked["druggability"] > 0) & ~ranked["possible_contamination"]]
shortlist = shortlist.drop_duplicates("gene").head(15)
shortlist.to_csv(ROOT / "results/tables/target_shortlist.csv", index=False)
pd.set_option("display.width", 200)
print(shortlist[["gene", "cell_type", "direction", "log2FoldChange", "padj", "drug_class",
                 "brain_focus", "specificity", "priority"]].round(3).to_string(index=False))
print(f"\n{len(hits)} hits scored, {hits['possible_contamination'].sum()} flagged as likely contamination")
