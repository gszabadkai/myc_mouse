---
name: fgsea-takes-one-list-not-two
description: "Why a rank-STABLE Myc effect gives strongly NEGATIVE timeline NES - fGSEA never compares two rankings, and the negative timeline is what stability produces"
metadata: 
  node_type: memory
  type: project
  originSessionId: 71ffe8cc-b22d-412b-9e40-f22881e23306
  modified: 2026-08-16T22:55:51.762Z
---

**fGSEA never compares two rankings.** It takes ONE ranked gene list and asks where a set sits in
it. The author's 2026-08-16 question — *"if the ranking of the Myc effect doesn't change between 6W
and 12W, why does `outputs/fgsea/` show strongly negative NES on the timelines?"* — dissolves once
you count the lists. There are **four**, not two:

`myc_6W` · `myc_12W` · `6>12W_wt` · `6>12W_myc`

"The ranking is preserved" is about the first two (Spearman **0.933** over 866 *sets*). The Myc+
timeline is built from the **difference** between them — the part that was lost.

1. **NES is scale-free, so it cannot see the fade where you'd look for it.** Median |LFC| over the
   1,967 Myc-DE genes falls **0.678 → 0.311 (ratio 0.46)** while median lfcSE is flat
   (0.189 → 0.193, so not power) — and median NES *rises*, **+1.768 → +2.008**. The fade must
   therefore surface in a contrast that spans the ages.
2. **Rank stability FORCES the mirror.** `interaction = myc_12W − myc_6W`, so a rank-stable effect
   at half amplitude makes the interaction a near-perfect negative copy:
   Spearman(myc_6W, interaction) = **−0.878**, (myc_6W, `6>12W_myc`) = **−0.838**.
   **The negative timeline is not evidence against stability — it is what stability produces.**
   A rank-*unstable* fade would give a weak, incoherent timeline.
3. **The compression is toward ZERO, not downward**, so the fade *lifts* the sets Myc represses:
   median `6>12W_myc` NES is **−1.80** over the 579 sets Myc induces and **+1.03** over the 287 it
   represses. A fade gives a **sign flip relative to the Myc effect**, never a uniform decline.
4. **The wild-type timeline is a different question and contains no Myc**: Spearman(myc_6W,
   `6>12W_wt`) = **−0.068**, an average over two opposite-signed families that both *fall*
   (induced −0.314, repressed +0.345). Its top negatives are stem/TEB/alveolar sets.
5. **"Most pathways negative" describes the LIBRARY, not the transcriptome.** At gene level
   `6>12W_wt` is balanced (47.6% negative, median stat +0.085) yet 73.9% of NES are negative —
   605 of 866 sets are TF-target + mammary-development lanes. By category the WT timeline is *not*
   uniformly negative: `01_mitocarta` **+0.81**, `04_metabolism` +0.88, against `06_tf_targets`
   −1.46. **The mito compartment RISES in the maturing WT gland; only the OXPHOS arm falls (−2.67)**
   — the substrate half of the closing claim.

**Two sentence risks this exposes** (both on the `PANELS.md` list): "the ranking did not change" is a
ranking of gene **SETS** (0.933; at gene level only **0.724** on Myc-DE genes, **0.520** overall);
and a temporal NES and a genotype NES are **not on the same ruler** — each is normalised inside its
own list, so −2.17 and +2.86 are not comparable magnitudes.

Drawn as a worked simulation in `figures/panels/explainer_nes_sign.R` (four regions, real
`fgsea::fgsea()` calls: **+3.57, +3.58, +0.95 ns, −3.43**; filename deliberately not `fig*` so it
stays outside both runners' globs). Recorded in `PANELS.md` under Fig. 1H, whose legend now computes
the reconciliation. Cite `results/fgsea_percategory.rds`, **not** script 04's superseded 89-set
`outputs/fgsea/` — only the newer object carries the `interaction` ranking.

Related: [[abundance-ruler-dependence]], [[background-vs-myc-is-additive]],
[[intro-oxphos-restraint-tension]].
