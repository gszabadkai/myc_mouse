---
name: mito-content-bounded
description: "Script 32: Myc raises mitochondrial content +21-27% (survives adjustment) - the one solid new result. Several claims from this script were RETRACTED; read the retraction list."
metadata: 
  node_type: memory
  type: project
  originSessionId: 71ffe8cc-b22d-412b-9e40-f22881e23306
---

**The question (author, 2026-07-16):** the manuscript says Myc drives a "robust increase in
mitochondrial biogenesis" -- is that mitochondrial CONTENT per cell, or only a program? Answered
by `scripts/32_mito_content_proxies.R` (written 2026-07-16, commit 50ad21c on `paper-figures`).

**Block A never tested it, and three objects that LOOK like they did, do not** -- do not be fooled
by their names again:
- `reframe_supp3$abund` ("AP-abund", script 22) = mitoPPS diff **DISPERSION**, not abundance.
- **mitoPPS is abundance-blind by construction** (`08:404` "cancels out total mitochondrial
  content") -> `mtnuc_index` (script 24) is a priority RATIO.
- `total_mito_score` (`08:1153`, captioned "reflects mito content") -- plotted, **never saved or
  tested**; double-counts genes over 142 overlapping pathways; **76.9% mt-driven**.
- `bio_comp` (`19:91`, the rho=0.90 "Myc footprint IS biogenesis") = **GSVA** composite = program.

**Result (genotype axis QC-clean, all 17 panels leave-one-out stable):**
- Myc RAISES the mito share of the transcriptome: mass markers **+27%** (chaperone-free, d=1.34,
  p=0.003), nuclear MitoCarta +19%, mitoribosome +35%, import +33%, nuclear OXPHOS +28%. Mass
  panel **coherent: 13/14 genes up, 11 at p<0.05**. Tspo the lone non-mover.
- Concordant with script 29 by a **different method** (nuclear OXPHOS d=1.15 vs 1.27).
- ~~Myc does NOT scale mtDNA-encoded output (+6%, p=0.83) = the imbalance in ABSOLUTE SHARE.~~
  **RETRACTED 2026-07-16** — this compared a significant p to a non-significant p instead of testing
  the DIFFERENCE. Nuclear-minus-mtDNA per sample: **p=0.57 raw, p=0.82 adjusted**. Script 32 neither
  supports nor refutes the imbalance. See [[mitonuclear-imbalance-untested]].
- **WHERE THE IMBALANCE IS *NOT* (script 32 PART 3b) -- an informative NEGATIVE.** Nuclear mito
  compartment median **+21%** (934 expressed genes); the mtDNA-dedicated core (Tfam/Polg/Polg2/Twnk/
  Ssbp1/Tfb2m/Polrmt/Tefm) rises WITH it, **unremarkably**: median +13%, **39th percentile**,
  **Wilcoxon p=0.12**. So the imbalance is **not a transcript-level deficit anywhere in the
  machinery** -- nothing in the transcriptome is the bottleneck. Cause sits **downstream of
  transcript abundance** (mtDNA copy number / transcription rate / turnover). **mtDNA qPCR is the
  discriminating experiment.**
- The rebuilt `total_mito_score` reproduces script 08's null (p=0.67) and **explains** it (mt
  dominance), rather than replacing it.
- Chaperone reconciliation: script 24's "Chaperones fall -> UPR^mt REFUTED" is a mitoPPS result;
  in absolute share chaperones are Myc-elevated and age-flat. **Not a contradiction** -- script
  24 stands.

**Three ceilings (state them; do not let the claim inflate):**
1. Bulk polyA, no spike-ins, no cell counts -> per-cell content is **not measurable**;
   median-of-ratios normalisation removes that scale. Every number is a **share of transcriptome**.
2. Myc **globally amplifies total RNA per cell**, so a flat share would already mean more
   mitochondria => **+27% is a LOWER BOUND**.
3. Transcript != protein != organelle volume. **Blot (TOMM20/VDAC/CS) / mtDNA qPCR / EM settles
   it** -- the SAME orthogonal experiment tension A (protein-level OXPHOS) already needs, so one
   blot answers both. **No orthogonal mass data exists** (author-confirmed 2026-07-16;
   transcriptome only -- the Menegollo companion does NOT carry it).

**MASS_MARKERS is the corpus's only hand roster** -- a flagged exception to CLAUDE.md's "do not
rebuild gene sets" (the library has no mass-marker set). Hspa9 (Grp75) / Hspd1 (Hsp60) were
author-requested; they are standard mass markers **and direct MYC targets** (the two largest
effects, +0.73/+1.03 log2), so they sit in their own arm and the **claim-bearing composite
excludes them**.

**RETRACTED same-day (2026-07-16): "COMMISSIONED BUT UNBUILT".** The first version of script 32
read set-level `MITOCARTA_TRANSCRIPTION` (+43%) / `MITOCARTA_MTDNA_NUCLEOID` (+22%, d=2.32) as "Myc
commissions the mtDNA machinery, output does not follow". **Set-composition artifact.** Gene by
gene the effects are carried by loosely-assigned Myc-induced members -- **Mrpl12 +116% (a
MITORIBOSOMAL protein), Atad3a +113%, Top1mt +76%, Poldip2 +58%** -- while the core barely moves
(Tfam +11% ns, Polg +3% ns, Twnk +1% ns, Polg2 -9% ns). The **opposite** claim ("Myc SPARES the
machinery") is equally unsupported (p=0.12) and was NOT substituted. Caught only after the author
had already run the script. See [[mitocarta-sets-are-membership-loose]].

**How to apply:** "biogenesis" can now be stated as a bounded CONTENT claim, not an enrichment
claim -- with the lower-bound rider. **This is the ONE solid new result from scripts 32/33; the rest
of what script 32 originally claimed was retracted.** Do NOT write the narrative as though the
transcriptome explains the mitonuclear imbalance. See [[mitonuclear-imbalance-untested]] for the
full ledger, [[timepoint-cohort-confound]], [[issue4-attenuation-scoped]],
[[intro-oxphos-restraint-tension]].
Written up in the walkthrough (script-32 section, tension-A (C), open questions, glossary).
**Author RAN it**; numbers reproduced exactly. **FIVE claims were retracted from script 32/33 over
2026-07-16/17** — see `docs/2026-07-17_evidence_audit_and_narrative.md` for the ledger. The CONTENT
claim (+21-27%, p=0.0006 adjusted for prep-stress + contamination) is what survived and is the
strongest new result. Do not re-litigate.

**OPEN (author question, 2026-07-16):** do the genes following the "mtDNA pattern" form a module?
The machinery does not. **Note the trap is NOT depth** (a share cancels depth by construction — that
warning was retracted) but **cohort structure + the 43% PC1**: everything correlates with everything,
so a scan would rank genes by that axis. Design within timepoint. Largely answered by script 33 PART
A/D2. See [[timepoint-cohort-confound]], [[sample-level-couplings-need-a-null]].
