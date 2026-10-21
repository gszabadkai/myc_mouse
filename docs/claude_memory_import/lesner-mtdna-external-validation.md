---
name: lesner-mtdna-external-validation
description: "Lesner/Simon bioRxiv 2026.07.13.738248: TCGA+mouse externally validate the MYC nuclear-up/mtDNA-down discordance (non-artifactual), reframe the death arm as protective mitophagy (not priming), and endorse our mtDNA-qPCR + BH3 + NRF1/DRP1/FUNDC1 experimental plan"
metadata:
  node_type: memory
  type: reference
  originSessionId: 71ffe8cc-b22d-412b-9e40-f22881e23306
  modified: 2026-07-20T09:38:36.782Z
---

**Lesner, Kim, ... M. Celeste Simon, bioRxiv 2026.07.13.738248** (`docs/MYC_mtDNA_paper.pdf` +
`docs/MYC_mtDNA_paper_supp.pdf`). Read 2026-07-20 at the author's request to assess relevance to our
mtDNA / mitonuclear-imbalance / priming question. **Highly relevant.**

**What it shows:** MYC-high tumors ENRICH nuclear-encoded mito transcripts but DEPLETE mtDNA-encoded
transcripts -- a mitonuclear DISCORDANCE -- across TCGA cancers and specifically LIHC (Fig 1D/E),
separated by genomic origin, with a non-mito MYC-score CONTROL, corroborated by proteomics and mouse
HCC mtDNA qPCR (mt-ND2 / nuclear Pecam1). Mechanism: MYC -> NRF1 -> DRP1/FUNDC1 mitophagy (mtDNA
depletion is TURNOVER, not failed biogenesis); mitophagy is PROTECTIVE (DRP1 KO improves survival);
apoptosis is engaged by forced FUSION + BH3 mimetics.

**Relevance to us (four points):**
1. **VALIDATES our imbalance PHENOMENON and kills the "our dissociation artifact" worry** -- TCGA is
   flash-frozen, no enzymatic digest, yet shows the same nuclear-up/mtDNA-down split we see. So the
   discordance is real biology, not our prep. (Directly relevant to [[mtpct-is-three-way-confounded]]:
   the CONTENT half is corroborated even though our mt% is confounded.)
2. **REFRAMES the death arm** -- mtDNA depletion is mitophagic turnover and PROTECTIVE, apoptosis
   needs forced fusion + BH3 mimetics. So it does NOT rescue our mt%<->priming coupling (consistent
   with our null; [[death-priming-change-untested]]).
3. **HONEST NUANCE** -- their discordance is a MYC CONTRAST. Our nuclear-UP half reproduces as a
   clean GENOTYPE effect (safe), but our mtDNA-DOWN half is genotype-INDEPENDENT (p~0.49, the
   confounded mt%/time axis). So Lesner validates the QUESTION and pinpoints where we are
   underpowered -- it is NOT a claim that we demonstrated the contrast.
4. **SCOPE caveat** -- their mechanism is partly liver-specific (ammonia/urea/HIF). Transferable
   core = MYC -> NRF1 -> mito-turnover with OXPHOS suppressed despite nuclear induction.

**Endorsed experiments/methods we can cite:** mtDNA qPCR (mt-ND2 / nuclear); non-mito MYC score as a
control; sgNRF1/sgDRP1/sgFUNDC1 perturbation -- which is exactly the ESRRA/PGC1a/NRF1 PERTURBATION
that script 38 concluded is needed to make the MYC-vs-PGC1a biogenesis/death contrast testable (bulk
n=24 cannot; [[death-priming-change-untested]]).

**How to apply:** cite Lesner as EXTERNAL VALIDATION of the discordance phenomenon and as endorsement
of the qPCR/BH3/NRF1 experimental plan -- NOT as support for a transcriptomic priming coupling.
Full read + honest-nuance framing: `docs/2026-07-20_lesner_mtDNA_external_validation.md` (new doc,
specced in the plan ADDENDUM 2, build on author approval). See [[mitonuclear-imbalance-untested]],
[[correlation-ceiling]], `docs/2026-07-19_oxphos_axis_biology_and_mtdna_priming.md`.
