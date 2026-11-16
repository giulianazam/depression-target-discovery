"""Step 7: does the ranking depend on the weights I chose?

Re-rank the candidate targets with 2,000 random weight combinations (all four
weights positive and summing to 1) and count how often each gene comes out on top.
If a gene stays near the top whatever the weights, the finding is robust.
"""
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SCORES = ["evidence", "druggability", "brain_focus", "specificity"]
N_DRAWS = 2000

hits = pd.read_csv(ROOT / "results/tables/all_hits_scored.csv")
hits = hits[(hits["druggability"] > 0) & ~hits["possible_contamination"]].reset_index(drop=True)

rng = np.random.default_rng(0)
weights = rng.dirichlet(np.ones(len(SCORES)), N_DRAWS)
top1, top5 = {}, {}
for w in weights:
    order = hits.assign(p=hits[SCORES].values @ w).sort_values("p", ascending=False)
    genes = order["gene"].drop_duplicates().tolist()
    top1[genes[0]] = top1.get(genes[0], 0) + 1
    for g in genes[:5]:
        top5[g] = top5.get(g, 0) + 1

summary = pd.DataFrame({"pct_ranked_1st": pd.Series(top1), "pct_in_top5": pd.Series(top5)})
summary = (summary.fillna(0) / N_DRAWS * 100).round(1).sort_values(
    ["pct_ranked_1st", "pct_in_top5"], ascending=False)
summary.index.name = "gene"
summary.to_csv(ROOT / "results/tables/weight_sensitivity.csv")
print(summary.head(10).to_string())
