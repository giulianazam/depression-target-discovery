# Independent audit of this analysis

*Done 12 November 2026 by a second, separate review of the inputs, the code and the
claims, without assuming the earlier work was right. Findings are ordered: real errors
first, then things that could bias the result, then wording to change. Each one says
what was checked and what to do about it.*

> **Status: every fix below has been applied.** This document is kept as the record of
> what was checked, not as a to-do list. The figures and the README now reflect all of
> it, including the rescaled-evidence re-ranking suggested in §2.1 and the
> genes-per-nucleus covariate suggested in §2.3, both of which left GPR37L1 on top.

**Headline:** nothing in the pipeline is broken, and the central result (GPR37L1 lower
in astrocytes of men with MDD) survived every check I could run, including a newly
published independent-ish analysis that found the same gene. The weak points are not in
the statistics of that one gene; they are in **how the "top target" was chosen** and in
**three sentences of wording** that claim more than the data shows.

---

## 1. Errors to fix

### 1.1 "In line with the original study" is only half true

README key finding 2 says excitatory neurons show the most changes, "in line with the
original study". Nagy et al. (2020) reported that deep-layer **excitatory neurons and
OPCs together** accounted for 47% of all expression changes, and the paper's title names
OPCs first. In this reanalysis OPCs have **8** significant genes out of 5,687 tested,
versus 962 in excitatory neurons. So the excitatory-neuron half agrees and **the OPC half
does not reproduce here**.

*Why it matters:* a reader who knows the paper will spot this immediately, and it is the
kind of thing a reviewer asks about first. It is also informative rather than
embarrassing: the likely reason is that Nagy et al. tested 26 fine-grained clusters with
their own model, while this pipeline tests 7 broad cell types adjusting only for batch.

*Fix:* change to something like "Excitatory neurons show the most gene changes (962 genes
at FDR < 0.1). Nagy et al. also found deep-layer excitatory neurons among the most
affected, but their other headline cell type, OPCs, shows only 8 changed genes here —
probably because they tested 26 fine-grained clusters and adjusted for more covariates
than this broad-cell-type pipeline does."

### 1.2 Maitra et al. is called a replication of the astrocyte drop, but it shares the same donors

README key finding 1 says "a later reanalysis of the same men plus a female cohort found
the same drop". Two problems:

- Maitra et al. (2023) reprocessed **these same 33 men** and added 38 women, so 33 of
  their 71 donors are the donors analysed here. That is not independent replication; it
  is the same tissue analysed better.
- Their effect is much smaller than this one. They report astrocyte proportions
  significantly reduced in cases (FDR 3.5e-4) alongside OPCs (5.3e-4), with excitatory
  neurons *increased* (FDR 0.048) — a modest compositional shift. Figure 2 here shows
  **7.8% → 2.3%, a 3.4-fold drop**. A 3.4-fold loss of astrocytes in the cortex is not
  biologically plausible and is far larger than anything in the literature, so most of
  that gap is almost certainly technical (how much white matter each dissection caught,
  how well astrocyte nuclei survived each prep).

*Fix:* say "partly the same donors, reprocessed with more covariates, found astrocyte
proportions significantly reduced (FDR 3.5e-4), though by much less than the 3.4-fold
difference seen here — so treat the size of the gap in Figure 2 as technical, and only
the direction as supported."

### 1.3 Figure 2's denominator includes the "Mixed/unclear" cluster

`04_figures.py` computes each donor's astrocyte share over **all** nuclei, including the
`Mixed/unclear` cluster, which every other step of the pipeline drops. If the share of
Mixed nuclei differs between groups (likely, since it tracks sample quality), part of the
astrocyte difference is just a different denominator.

*Fix:* recompute the proportions after dropping `Mixed/unclear`, and report both. It is a
two-line change and will probably shrink the gap.

### 1.4 The donors dropped from the astrocyte test are the ones that drive Figure 2

The DE step requires 10 astrocyte nuclei per donor. For astrocytes that leaves 30 of 34
donors (14 MDD, 16 control) — so **3 of the 4 dropped donors are MDD cases**, and they
are dropped precisely because they had almost no astrocyte nuclei. Figure 2 includes
those men; the GPR37L1 test excludes them.

*Why it matters:* this is not a bug (you cannot run DESeq2 on 4 nuclei), but it means the
two headline results use different samples, and the DE result is conditioned on MDD
donors that had enough astrocytes. Honest reading: this makes GPR37L1 *more* trustworthy,
not less, because it was found without the most extreme samples.

*Fix:* state the donor counts in the README (astrocytes: 14 MDD vs 16 controls of 34) and
say which direction the exclusion pushes.

### 1.5 Small code and claim hygiene

- `02_pseudobulk_de.py` carries the comment "Library prep batches are balanced across
  groups (checked)" but nothing in the code checks it. Print a `batch × condition`
  crosstab so the claim is reproducible. If a batch is nearly all one group, the batch
  term is partly absorbing the diagnosis effect.
- `01_load_and_qc.py` applies **no QC filtering of its own** — no minimum genes per
  nucleus, no mitochondrial cut-off, no doublet removal. It calls
  `calculate_qc_metrics` and never uses the result. This is defensible (the GEO matrix is
  the authors' already-filtered, already-annotated one) but the README's method table
  implies a QC step. Say plainly: "no extra QC; we use the authors' filtered matrix and
  their cell labels as published."
- The README cites "Open Targets tractability assessments" among frameworks the scoring
  follows. Open Targets is not queried anywhere; the scoring is a regex over HGNC family
  names. Keep the citation only as inspiration ("in the spirit of"), or drop it.

---

## 2. Bias risks (the important section)

### 2.1 The scoring cannot let evidence outrank category membership

This is the single biggest methodological weakness, and the weight-sensitivity analysis
cannot see it.

Three of the four scores are **saturating categories**: druggability is 1.0 for any GPCR,
ion channel or kinase; brain focus hits its 1.0 ceiling at a 32-fold brain/body ratio;
specificity is 1.0 for any gene whose top cell type is the one being tested. The fourth,
evidence, is built so it **can never approach 1**: its maximum across all 1,075 hits is
**0.56**, and GPR37L1's is **0.38**.

So the nominal "0.35 × evidence" contributes at most 0.35 × 0.56 = 0.20 to any gene's
score, while "0.30 × druggability" contributes a full 0.30 to hundreds of genes. In
practice **the ranking is a ranking of protein families, with strength of evidence as a
tie-breaker.**

Worse, only **8 of the 1,075 hits are GPCRs**, and GPR37L1 is the only gene in the entire
list that scores 1.0 on druggability, brain focus **and** specificity at once. Given the
scoring design, no other gene could have come first.

*Why it matters:* the README says the ranking "does not hinge on the chosen weights"
because GPR37L1 ranks first in 96.6% of 2,000 random weight draws. That statement is
true but close to circular: randomising the weights cannot dislodge the only gene that
maxes three of four axes. The sensitivity analysis tests the weights and not the much
more consequential choices — the 0-to-1 scaling, the ceilings, and the decision to score
druggability by family.

*Fix (choose one, cheapest first):*

1. **Reword only.** Say "GPR37L1 is the only hit that is both a GPCR and brain-specific
   and astrocyte-specific, so it comes first under any weighting. The weight sensitivity
   confirms this but does not independently validate it." This is honest and costs
   nothing.
2. **Rescale evidence** so it spans 0–1 across the observed hits (e.g. rank-percentile of
   p-value within cell type), re-run, and report whether the shortlist changes.
3. **Perturb the scores, not the weights** — e.g. drop one axis at a time, or score
   druggability 1.0 / 0.5 / 0 by Open Targets tractability bucket instead of HGNC family.

Doing (1) plus (2) would make this part of the project genuinely defensible.

### 2.2 Multiple testing is per cell type, not across cell types

FDR is computed inside each cell type, so 7 families of tests are each corrected
separately. I recomputed Benjamini-Hochberg across **all 66,942 tests at once**:

| | per-cell-type FDR (as published) | global FDR (my check) |
|---|---|---|
| Significant genes, FDR < 0.1 | 1,075 | 448 |
| GPR37L1 in astrocytes | 0.049 | **0.047** |
| GPR37L1 in OPCs | 0.026 | **0.026** |

**Good news:** GPR37L1 survives global correction essentially unchanged, in both cell
types. **Caveat:** the gene *counts* in Figure 3 are inflated roughly two-fold relative
to a global correction, so the figure caption should say "FDR within each cell type".

### 2.3 Astrocyte quality differs between groups, and that is the real rival explanation

Three things point the same way: the astrocyte UMAP island splits by diagnosis; many of
the 43 astrocyte hits are **neuronal** genes (AP1S1, TUBB4A, STMN3, TUBB, DYNLL2 are the
top 5 by p-value, all down in MDD); and MDD astrocyte samples are smaller. A generalised
loss of RNA complexity in MDD astrocyte nuclei would produce exactly this pattern.

What already argues against it: core astrocyte markers (GFAP, AQP4, SLC1A2/3, ALDH1L1,
SOX9) are unchanged (FDR 0.6–0.9); the effect holds with astrocyte % as a covariate
(-1.29, FDR 0.055); and DESeq2's size factors normalise total depth.

What would settle it, in order of effort: (a) add median genes-per-nucleus per donor as a
covariate alongside astrocyte % — a 3-line change to `06_gpr37l1_robustness.py` and the
most convincing single check you can still run; (b) show GPR37L1 as a fraction of
astrocyte nuclei expressing it, not just pseudobulk mean; (c) report the top-5 astrocyte
hits in the README so a reader can see they are neuronal and judge for themselves.

I would do (a) before the poster. It is the question an examiner will ask.

### 2.4 Confounds that cannot be fixed with this dataset

All cases died by suicide; medication, smoking, alcohol and cause of death are not in the
GEO cell names, so nothing adjusts for them. Age, pH and PMI are likewise unavailable
from the processed matrix (Maitra et al. had them and did adjust — and found GPR37L1 more
strongly, which is reassuring). The README's Limitations already say this; keep it there.

---

## 3. Literature check

I checked every substantive external claim. Results:

| Claim | Verdict |
|---|---|
| GPR37L1 is an astrocyte-enriched orphan class-A GPCR | **Correct.** HGNC lists it under "G protein-coupled receptors, Class A orphans". |
| GPCRs are the most successfully drugged family | **Correct.** ~34% of FDA-approved drugs act on GPCRs (Hauser et al., *Nat Rev Drug Discov* 2017). |
| GPR37L1 ~300x more expressed in brain than other organs | **Correct as computed** (GTEx v8 median: brain 42.6 TPM vs 0.14 TPM median non-brain = 312x). Say "median non-brain tissue" rather than "other organs". |
| Prosaposin / prosaptide as GPR37L1 ligand is debated | **Correct, and worth keeping hedged.** Meyer et al. (*Glia* 2018, PMC6492175) support prosaposin acting via GPR37L1/GPR37; Jolly/Richardson et al. (*Glia* 2017) failed to reproduce prosaptide signalling at GPR37L1. Both sides are live. |
| Nagy et al. used 10x Chromium | **Correct**, confirmed via Maitra et al. 2023, which reprocessed the male data: 10x Genomics Chromium 3' v2/v3, Cell Ranger, GRCh38. |
| Pseudobulk is the recommended approach (Squair et al. 2021) | **Correct and standard.** This was the right call and is the strongest design choice in the project. |
| Maitra et al.: astrocytes/neurons in men, microglia in women | **Correct.** Astrocytes ~60% of male broad-level DEGs; microglia ~87% of female broad-level DEGs, plus PV interneurons. |
| GABRD is the preferred target of zuranolone-type neurosteroids | **Overstated.** Neurosteroids are positive allosteric modulators at both synaptic and extrasynaptic GABA-A receptors, and δ-containing extrasynaptic receptors are *particularly* sensitive — but they are not δ-selective. Say "δ-containing extrasynaptic GABA-A receptors are especially sensitive to neurosteroids like those behind brexanolone and zuranolone", not "the target". |
| HCN2 is implicated in rodent depression models | **Correct but relocate.** The evidence is in nucleus accumbens cholinergic interneurons and VTA, not prefrontal excitatory neurons where the hit is. Say "HCN channels have been linked to depression-like behaviour in rodents, though in other brain regions". |
| Esketamine (Spravato) is an approved nasal spray; chitosan opens epithelial tight junctions | **Correct.** The delivery design document is sound as a proposal. |

### The important new paper: GPR37L1 was independently flagged four months ago

**Mitsuhashi et al., "Multi-omics characterization of astrocyte subtypes reveals
spatially coordinated astrocyte downregulation in depression", bioRxiv, 3 July 2026**
(PMID 42427711, PMC13345167; corresponding author Corina Nagy — the same group).

They combined snRNA-seq (71 donors), snATAC-seq (84 donors) and Visium spatial
transcriptomics (16 donors), and found:

- **GPR37L1 is reduced in protoplasmic astrocytes in MDD (adj. p = 0.0095), with reduced
  chromatin accessibility at the locus** — an independent data type pointing the same way.
- Reduced GPR37L1 expression in cortical **layers 5–6** specifically (spatial data).
- **GPR37L1 is the top-ranked receptor in their astrocyte communication factor, paired
  with PSAP (prosaposin)**, and PSAP–GPR37L1 signalling from deep-layer excitatory
  neurons to astrocytes is reduced in MDD.
- Male-driven: GPR37L1 DE in males (adj. p = 0.0064), not in females (adj. p = 1) —
  matching the male/female split found here exactly.

**What this means for the project, honestly stated:** this is strong external support for
the biology, and it independently picks the same ligand–receptor pair the delivery design
is built on, which is a real vindication of that choice. But it is **not** independent
replication of the data — it is the same McGill cohorts (Nagy 2020 + Maitra 2023) with
more modalities, from the group that generated them. And it means the finding is no
longer novel: as of July 2026 there is a preprint naming GPR37L1 as a depression-relevant
astrocyte receptor.

*What to do:* cite it. It strengthens every claim about GPR37L1 and the PSAP link, and
citing it is far better than having someone point out you missed it. Reframe the project's
contribution from "we found a new target" to "an independent, transparent pipeline
converges on GPR37L1 from the published male data alone, and the gene has since been
corroborated with chromatin and spatial data". That is a stronger and more defensible
claim for a student project, not a weaker one.

---

## 4. What I did not find

For completeness, things I looked for and did not find:

- No arithmetic or indexing errors in the pseudobulk summation, the donor/condition
  parsing, or the DESeq2 contrast direction (negative log2FC = lower in MDD — correct).
- The GPR37L1 numbers in the README match the tables exactly (-1.3926, FDR 0.0494).
- The leave-one-donor-out table is internally consistent: 30 refits, all p < 5e-4, log2FC
  -1.16 to -1.55 — no single donor drives the result.
- No circularity in the druggability or brain-focus scores: both come from external
  sources (HGNC, GTEx) that know nothing about this dataset.
- No sign of p-hacking in the DE step: one model, one contrast, applied uniformly to
  every cell type.

## 5. Suggested order of work

1. Reword README key findings 1 and 2, and the weight-sensitivity sentence (§1.1, §1.2,
   §2.1 fix 1). Half an hour, no re-running.
2. Add the Mitsuhashi et al. 2026 citation (§3). Fifteen minutes, high value.
3. Drop `Mixed/unclear` from Figure 2's denominator and state the astrocyte donor counts
   (§1.3, §1.4). One re-run of `04_figures.py`.
4. Add genes-per-nucleus as a covariate in `06_gpr37l1_robustness.py` (§2.3a). The one
   substantive new analysis worth doing before the poster.
5. Rescale the evidence score and re-run the ranking (§2.1 fix 2), if there is time.
