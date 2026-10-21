---
name: lesner-mechanisms-tested-script39
description: "Script 39 tested Lesner's hypoxia/mitophagy mechanisms in our data: the HIF->MXI1/JUN Myc-suppression mechanisms do NOT transfer (liver-specific); nuclear OXPHOS declines Myc-specifically but NOT hypoxia-linked; Myc raises mito-turnover machinery but NOT mitophagy-specifically (content effect)"
metadata:
  node_type: memory
  type: project
  originSessionId: 71ffe8cc-b22d-412b-9e40-f22881e23306
  modified: 2026-07-20T11:20:28.916Z
---

**Testing Lesner et al. (bioRxiv 2026.07.13.738248) in the MMTV-Myc mammary model**
(`scripts/39_lesner_mechanism_test.R`, 2026-07-20, `results/lesner_mechanism_test.rds`; reproduced
read-only). Instrument = **DESeq2 contrasts** (clean genotype gaps `myc_6W/12W` + Myc-specific
`interaction`), NOT per-sample GSVA couplings (the [[correlation-ceiling]]). BATCH=TIMEPOINT so
`timepoint_pos/neg` are batch-confounded; the **interaction is the clean term**. See
[[lesner-mtdna-external-validation]].

**1. Lesner's two HYPOXIA -> MYC-suppression mechanisms do NOT reproduce (liver-specific).**
- M1 (HIF-1a -> MXI1 -> inhibits MYC) and M2 (HIF-2/EPAS1 -> inhibits JUN -> reduces MYC/TFAM).
- HYPOXIA is **LOWER** in Myc+ (myc_6W NES **-1.34**, padj 0.01) and does **not** rise over time in
  Myc+ (timepoint_pos -1.14, ns); the +1.49 interaction is because hypoxia falls LESS in Myc+ than WT,
  not because it rises. Gene level: Hif1a low in Myc+ (-0.34); Mxi1 higher in Myc+ (+0.33) but **flat
  over time** (+0.06); Tfam **flat** everywhere (+0.18). The one thread: Epas1/HIF-2 rises over time in
  Myc+ (interaction +0.45) -- BUT its predicted downstream FAILS: Jun **rises** over time in Myc+
  (+0.39, not down), JUN targets not suppressed (interaction NES +1.52). => not transferable.

**2. Nuclear OXPHOS declines over time MYC-SPECIFICALLY, but NOT hypoxia-linked.**
- Hallmark OXPHOS: Myc raises it (myc_6W +2.99, myc_12W +3.40, padj 1e-38) then it DECLINES over time
  in Myc+ (timepoint_pos **-2.07**, padj 7e-9) with WT timecourse **ns** (timepoint_neg -1.04) and a
  clean Myc-specific **interaction -2.21 (padj 7e-11)**. It falls in lockstep with MYC targets
  (interaction -1.80) while HYPOXIA does not rise => the **general Myc-programme fade** (canonical
  attenuation, [[issue4-attenuation-scoped]]), not a hypoxia switch. Nuance: nuclear OXPHOS subunits
  decline in WT too (timepoint_neg -2.07 = Issue #6 convergence); the Hallmark set is the cleaner
  Myc-specific signal. Confirms [[intro-oxphos-restraint-tension]].

**3. Myc raises the mito-turnover machinery (genotype) but NOT mitophagy-SPECIFICALLY (content effect).**
- Genotype: DYNAMICS myc_6W NES **+2.23** (padj 1e-8), FUSION +1.85, MITOPHAGY +1.56, FISSION +1.55,
  AUTOPHAGY +1.62 -- Lesner's MYC->mito-turnover DOWNSTREAM transfers. Attenuates over time (interaction
  DYNAMICS -1.81). **But** vs an **expression-matched MitoCarta background** the genotype rise is BELOW
  average for a mito gene: MITOPHAGY obs +0.24 vs bg +0.39 (**z -1.36, p 0.92**), DYNAMICS +0.27 vs
  +0.39 (**z -2.62, p 1.00**) -- the [[mito-content-bounded]] ~20% wholesale content effect, exactly
  the [[death-priming-change-untested]] trap. Gene level: DRP1/Dnm1l +0.22, Fundc1 +0.14 weakly
  Myc-induced; big movers Bnip3/Bnip3l/Prkn rise EQUALLY in WT (shared/batch). Matched null is
  anti-conservative -> the NON-specific verdict is safe.

**Net (citable cross-system contrast):** same **MYC -> mito-turnover machinery** Lesner report, reached
by a **different, non-hypoxic route** here; our OXPHOS decline is the Myc-programme attenuation, not a
hypoxia switch; and the mitophagy induction is the general mito-content rise, not a mitophagy-specific
programme. **How to apply:** cite the genotype-level mito-dynamics induction and the Myc-specific OXPHOS
decline (clean contrasts); do NOT claim a hypoxia mechanism or a mitophagy-specific programme.
Direct-test experiments (mtDNA qPCR, BH3, NRF1/DRP1/FUNDC1 perturbation) remain the way to settle
mechanism -- [[lesner-mtdna-external-validation]]. Figures: `outputs/lesner/A_mechanism_genes.pdf`,
`B_fgsea_contrasts.pdf`, `C_mitophagy_specificity_null.pdf`.
