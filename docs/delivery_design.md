# From target to therapy: getting a GPR37L1 drug into the brain

*Proposed translational approach: a design concept connecting the data analysis to
biomaterial design. No experimental data; every number below is a target or a literature
value, not a measurement.*

## The problem

GPR37L1 is lower in astrocytes from people with depression. If restoring its signalling
helps, we would want a drug that **activates** it (an agonist). The leading candidate
ligands proposed in the literature are prosaposin-derived peptides (e.g. prosaptide,
TX14(A)), though whether they truly act through GPR37L1 is still debated.

Peptides have two big problems as brain drugs:
1. They barely cross the **blood-brain barrier** (BBB) from the bloodstream.
2. They are broken down quickly in blood and tissue.

Supporting this receptor-ligand pairing, a 2026 preprint from the group that produced the
data (Mitsuhashi et al.) found PSAP-GPR37L1 signalling from deep-layer neurons to
astrocytes reduced in depression, mainly in men [5].

## Proposed approach: a mucoadhesive nose-to-brain gel

The nose offers a route that partly bypasses the BBB: drugs placed high in the nasal
cavity can travel along the olfactory and trigeminal nerves into the brain. This is
already used clinically for depression: **esketamine (Spravato) is a nasal spray**,
approved by the FDA in 2019 for treatment-resistant depression [1].

The weak point of nasal delivery is that the nose clears anything within about
15–20 minutes (mucociliary clearance) [2]. That's where an adhesive material comes in.
Mucoadhesive and thiolated chitosan systems have been reported to extend nasal residence
from that 15–20 minute baseline to several hours [2].

**Design concept:**

| Component | Choice | Why |
|---|---|---|
| Carrier | Peptide loaded in chitosan nanoparticles (~100–200 nm) | Protects the peptide from enzymes; positively charged chitosan redistributes tight-junction proteins (claudin-1, ZO-1), reversibly increasing paracellular transport [2] |
| Gel | Thermo-responsive poloxamer 407 (~15% w/w, with ~1% poloxamer 188 to tune the transition) + mucoadhesive polymer at ~0.5% w/w (e.g. chitosan or hyaluronic acid) | Liquid at room temperature for spraying, gels at nasal temperature; a 15% P407 / 1% P188 nasal gel has been reported to gel at 34.5 ± 0.3 °C [3] |
| Adhesion | Positively charged chitosan binds negatively charged mucin | Holds the dose on the olfactory region longer than clearance would allow; raising a mucoadhesive polymer from 0.1% to 1% roughly doubled mucoadhesive strength in one P407 gel [3] |
| Targeting (optional) | Surface ligand for astrocytes or for nasal epithelial uptake | Increases the share of drug reaching the target cells |

## How this connects to my bioadhesive project

I worked on Kelserra, an algae-inspired surgical bioadhesive built from alginate,
xanthan gum and phloroglucinol, which sets when calcium is released from Ca-EDTA
[6]. I co-led the validation subteam, running the mechanical testing, the microscopy
and the degradation study, alongside our cytotoxicity work. Three design principles carry
straight over to a nasal gel. First, **adhesion to wet tissue is the hard part**: we
measured Kelserra on porcine skin both air-cured and kept wet in PBS at 37 °C, and
strength trended lower when wet, which is the same obstacle a gel faces on permanently
mucus-covered nasal mucosa. Second, **a material should stay liquid until a trigger sets
it in place**: Kelserra uses slow calcium release as that trigger, where the nasal gel
uses body temperature. Third, **a bioadhesive is
only as good as its stability**: tracking Kelserra over 28 days showed phase separation
by day 7 and visible darkening by day 28, so shelf life and storage belong in the design
from the start, not as an afterthought.

The validation methods also transfer, which is why the testing plan below mirrors them.
Lap shear and uniaxial tensile testing on tissue (ASTM F2258-05, converted per
ISO 4587:2003) becomes mucoadhesion testing on mucin discs; trypan blue viability
testing moves from bronchial to RPMI 2650 nasal epithelial cells; and the microscopy and
time-course stability design becomes the release and shelf-life study. Kelserra's
viability run was a single replicate on one cell line, so the nasal work should be
properly replicated from the start.

## Key design requirements

- Gels within seconds at 32–35 °C; sprayable below 25 °C (a 15% P407 / 1% P188 system
  sits in this range [3])
- Residence on the olfactory epithelium of 1 hour or more
- Releases most of the peptide over 1–4 hours
- Non-irritating to nasal tissue with repeated daily dosing. Keep nanoparticle zeta
  potential below about +30 mV: more cationic charge than that has been linked to
  membrane damage [2]
- Peptide stays intact (verify by HPLC after release)

## How I would test it

1. **In vitro:** gelation temperature (rheology), mucoadhesion strength on mucin discs,
   peptide release curve, toxicity on nasal epithelial cells (e.g. RPMI 2650).
2. **Ex vivo:** permeation across excised nasal mucosa in a diffusion cell.
3. **In vivo (later):** fluorescently labelled peptide in mice to measure how much
   reaches the prefrontal cortex vs. blood, then behavioural tests of depression-like
   behaviour after chronic stress.

## Open questions

- Is GPR37L1 loss a cause of depression, or a consequence of astrocyte loss?
  Restoring the receptor may not help if the astrocytes themselves are gone.
- Is prosaptide really a GPR37L1 agonist? Evidence is genuinely split [4], so a
  small-molecule agonist screen may be needed.
- Does GPR37L1 activation have effects outside the brain (it has been linked to blood
  pressure regulation in mice)?

## References

1. US Food and Drug Administration. *FDA approves new nasal spray medication for
   treatment-resistant depression* (esketamine, Spravato), 5 March 2019.
2. Matalqah S, Lafi Z, Atwan R, Asha SY, Mhaidat AN. *Intranasal chitosan
   nano-formulations: a novel frontier in brain cancer therapy.* Polymers 18, 2723
   (2026). [doi:10.3390/polym18212723](https://doi.org/10.3390/polym18212723)
3. Adnet T, Groo A-C, Picard C, et al. *Pharmacotechnical development of a nasal drug
   delivery composite nanosystem intended for Alzheimer's disease treatment.*
   Pharmaceutics 12, 251 (2020).
   [doi:10.3390/pharmaceutics12030251](https://doi.org/10.3390/pharmaceutics12030251)
4. Meyer RC, et al. *Prosaposin maintains lipid homeostasis in dopamine neurons* and
   related work supporting prosaposin signalling at GPR37/GPR37L1 (*Glia*, 2018);
   contrast with Coleman JLJ, Ngo T, Smith NJ, et al., who did not reproduce prosaptide
   signalling at GPR37L1 (*Glia*, 2017). The pairing remains disputed.
5. Mitsuhashi H, et al. *Multi-omics characterization of astrocyte subtypes reveals
   spatially coordinated astrocyte downregulation in depression.* bioRxiv, July 2026.
   [preprint](https://sciety.org/articles/activity/10.64898/2026.07.02.735636)
6. Fraser S, de Hepcée A, Fu A, et al. (incl. Zambito G). *Kelserra: an algae-inspired
   bioadhesive for use in minimally invasive surgeries.* Biotechnology Design 4, e8
   (2026). [doi:10.1017/S2977905726100195](https://doi.org/10.1017/S2977905726100195)
