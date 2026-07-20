# External validation: Lesner et al. on the MYC mitonuclear discordance

*2026-07-20. Read of `docs/MYC_mtDNA_paper.pdf` + `docs/MYC_mtDNA_paper_supp.pdf` (Lesner, Kim,
... M. Celeste Simon; bioRxiv 2026.07.13.738248), commissioned to assess relevance to our mtDNA /
mitonuclear-imbalance / death-priming thread. Companion to
`docs/2026-07-19_oxphos_axis_biology_and_mtdna_priming.md` and the evidence audit
`docs/2026-07-17_evidence_audit_and_narrative.md`.*

## 1. What the paper shows

Lesner et al. test the **MYC mitonuclear discordance** across cancer, separating mitochondrial genes
by **genomic origin**:

- **Fig 1 (TCGA, pan-cancer + LIHC):** a MitoCarta-based MYC score partitions tumors; **MYC-high
  tumors ENRICH nuclear-encoded mitochondrial transcripts but DEPLETE mtDNA-encoded transcripts** --
  a coordinated nuclear-up / mtDNA-down split, not a uniform mito change. A **non-mitochondrial MYC
  score** is used as a control (the effect is specific to the mito compartment, not a generic
  MYC-target readout).
- **Corroboration:** proteomics recapitulates the transcript split, and **mouse HCC mtDNA qPCR**
  (mt-ND2 normalised to nuclear Pecam1; methods p20) confirms reduced mtDNA content in MYC-high liver.
- **Mechanism:** MYC -> **NRF1** -> **DRP1 / FUNDC1** drives mitochondrial fission + **mitophagy**;
  the mtDNA depletion is **turnover**, not failed biogenesis. Mitophagy is **PROTECTIVE** (sgDRP1 /
  DRP1 loss worsens fitness / improves survival in their model); **apoptosis** is only engaged when
  mitochondria are forced to **fuse** (block fission) **plus** BH3 mimetics.

## 2. Why it is relevant to us

Our corpus arrived at a **mitonuclear imbalance** (nuclear OXPHOS subunits up while the mtDNA-encoded
fraction moves the other way) and asked whether that imbalance couples to **death priming**. Lesner
speaks to both halves.

**(a) It VALIDATES our imbalance PHENOMENON and kills the "our dissociation artifact" worry.**
Our biggest standing caveat was that our MECs are **enzymatically dissociated with no QC metadata**,
so any mt-fraction signal could be prep/leak (see
`[[mtpct-is-three-way-confounded]]` and `docs/2026-07-17_evidence_audit_and_narrative.md`). TCGA is
**flash-frozen bulk tissue, no digest**, and shows the *same* nuclear-up / mtDNA-down split. So the
discordance is real MYC biology, independently established, not an artifact of our preparation.

**(b) It REFRAMES the death arm -- and does NOT rescue our priming coupling.**
In Lesner's model the mtDNA loss is **protective mitophagy**, and apoptosis needs **forced fusion +
BH3 mimetics** -- i.e. mtDNA depletion is not itself a death signal. This is **consistent with our
null**: we found the mtDNA / imbalance <-> priming coupling collapses on mt%-adjustment (-0.71 ->
-0.20), and the only mt%-robust transcriptomic lead is **nuclear** OXPHOS <-> priming (+0.45 ->
+0.49, mid-pack) -> BH3 profiling. Lesner does not give us a transcriptomic death coupling; it points
the death question at the **apoptotic gateway** (BH3 / MOMP), exactly where we already said the answer
lives. See `[[death-priming-change-untested]]`.

**(c) HONEST NUANCE -- it validates the QUESTION, not our contrast.**
Lesner's discordance is a **MYC contrast** (MYC-high vs MYC-low tumors). In our data the two halves
behave differently:
- the **nuclear-up** half reproduces as a **clean genotype effect** (Myc raises content +21-27% and
  nuclear OXPHOS d~1.2; contrasts are safe -- `[[mito-content-bounded]]`,
  `[[mitonuclear-imbalance-untested]]`);
- the **mtDNA-down** half is **genotype-INDEPENDENT in our data** (the mt%/imbalance axis is
  time-associated but genotype p~0.49; it is the confounded axis of `[[correlation-ceiling]]`).

So Lesner corroborates the *direction* of the discordance and pinpoints **where we are underpowered**
(the mtDNA arm), rather than confirming that we demonstrated the MYC contrast on the mtDNA side.

**(d) SCOPE caveat.** Lesner's mechanism is partly **liver-specific** (ammonia / urea cycle / HIF
interplay in HCC). The transferable core for us is **MYC -> NRF1 -> mitochondrial turnover, with
respiration/mtDNA suppressed despite nuclear mito induction** -- which is the shape of our data
(nuclear OXPHOS up, mtDNA fraction not tracking it).

## 3. Experiments / methods it endorses for us

Lesner's toolkit maps almost one-to-one onto the experiments our audit already nominated:

| Our open question | Lesner's method we can cite |
|---|---|
| Is the mtDNA/imbalance real (vs prep)? | **mtDNA qPCR**: mt-ND2 / nuclear Pecam1 (physical ratio, immune to the mt% transcript confound) |
| Is the discordance MYC-specific? | **non-mitochondrial MYC score** as a control axis |
| Is the death arm mitochondrial? | forced **fusion** + **BH3 mimetics** to engage apoptosis; BH3 profiling reads the gateway |
| MYC- vs PGC1a/ESRRA-biogenesis toward death (untestable in our n=24 bulk, script 38) | **sgNRF1 / sgDRP1 / sgFUNDC1** perturbation -- exactly the biogenesis-axis perturbation script 38 concluded was required |

The last row is the key link: script 38 found MYC- and PGC1a/ESRRA-biogenesis are **not separable in
bulk** (r=0.93, the global-factor over-removal trap) and that the death contrast between them **needs
a perturbation**, not more n=24 correlation. Lesner **ran that class of perturbation** (NRF1/DRP1/
FUNDC1) -- so we can cite it as both external precedent and the template for the experiment we would
propose.

## 4. Bottom line

- **Cite Lesner as external validation** that the MYC nuclear-up / mtDNA-down discordance is real,
  cross-cancer, and non-artifactual (flash-frozen TCGA + mouse HCC + proteomics + qPCR).
- **Do NOT cite it as support for a transcriptomic priming coupling.** It reframes the mtDNA arm as
  **protective mitophagy / turnover**, consistent with our finding that the death arm has no
  transcriptomic support.
- **Use it to justify our experimental plan:** mtDNA qPCR (mt-ND2 / nuclear), non-mito MYC control,
  BH3 profiling, and NRF1/DRP1/FUNDC1 perturbation for the MYC-vs-PGC1a biogenesis/death question.

See `[[lesner-mtdna-external-validation]]`, `[[mitonuclear-imbalance-untested]]`,
`[[death-priming-change-untested]]`, `[[mtpct-is-three-way-confounded]]`,
`docs/2026-07-19_oxphos_axis_biology_and_mtdna_priming.md`.
