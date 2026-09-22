# The five questions — narrative summary and how to write the hard parts

**Date:** 2026-07-18. **Status:** current reader-facing summary; supersedes the pre-correction
snapshots `2026-07-13_BlockA_revision_walkthrough_and_intro_alignment.md` and
`2026-07-17_evidence_audit_and_narrative.md` (both left unedited as dated records of how the
analysis stood before the coupling method was rebuilt). The full technical version with every
number and named figure is `2026-07-18_narrative_synthesis_five_questions.md`; this file is the
plain-language digest plus explicit guidance on the two threads that are easy to overclaim — the
**cell-death priming decline** and the **mtDNA / mitonuclear imbalance**.

**Which docs reflect the latest results (script 36):** only the technical synthesis
(`2026-07-18_narrative_synthesis_five_questions.md`) and this file. The 07-13 and 07-17 docs do
not; treat them as history.

---

## 1. The five questions, in one line each

1. **Do the attenuation mechanisms hold?** **Yes.** The Myc programme's 6W->12W shrinkage is a
   genotype x time contrast (~66% Myc-fade + ~34% convergence toward the wild-type). This is
   solid and unaffected by everything below.
2. **Is there mitochondrial reprioritisation?** **By genotype, yes; over time, no.** Myc raises
   mitochondrial **content ~21-27%** and raises absolute nuclear-encoded OXPHOS. These are
   comparisons between Myc and wild-type mice, and they are safe. The *time-course* mitochondrial
   story (the "imbalance") is **not** supported — see §4.
3. **Is there evidence for Bbc3/PUMA driving the death difference?** **No.** One gene out of 23 is
   exactly what chance gives. It is a fine hypothesis for a bench assay, not evidence from this data.
4. **Is "OXPHOS, not biogenesis, is central" still true?** **Half.** "Biogenesis is a bystander"
   holds and is strengthened. "OXPHOS is central" fails — and after the rebuilt analysis it is not
   even **testable** (see §2). The one coupling that survives is a non-Myc one: **lower antioxidant
   (redox) capacity tracks the tumour-like state**, worth a bench test.
5. **Does the Myc-on-development effect hold?** **Yes for the powered half** — Myc amplifies the
   pubertal programme and suppresses the basal/myoepithelial one (genotype comparisons, safe). The
   "endogenous Myc gates the wild-type programme" sub-claim is weak (small numbers) and should be
   softened.

The dividing line under all five: **comparisons between groups survived; correlations between
programmes within the samples did not.** The next section explains why, without jargon.

---

## 2. The "global factor" is itself a phenotype — explained plainly

This is the single most important idea for the paper, and it is easy to state without any
bioinformatics.

**The problem, by analogy.** Imagine photographing 24 mice, but the camera's exposure drifts from
shot to shot — some pictures come out brighter overall, some darker. Now ask: "across these
photos, does the amount of red track the amount of blue?" You will find that it does — but only
because the bright photos have more of *everything*. Red and blue rising together is not evidence
that red and blue are specially connected; it is the overall brightness moving them in lockstep.

**What our data does.** For each mouse we measure the activity of hundreds of biological
"programmes" (mitochondria, cell growth, metabolism, and so on). It turns out that in these 24
samples there is **one overall dial** — a single background level that lifts or lowers almost
every programme at once. We measured how strong it is: pick **any two unrelated programmes at
random**, and they track each other about as tightly as the "interesting" pairs we cared about.
So a strong correlation between, say, "mitochondria" and "growth" is, by default, **not special** —
it is the dial moving both.

**Where the dial comes from — we checked, and ruled out the boring explanations.** It is *not*
because we have few mice (random shuffles give almost nothing). It is *not* because the gene lists
overlap (unrelated lists still track together). It is *not* only the Myc-versus-normal difference
(the dial is still there even inside a single group of six mice of the same type and age). When we
mathematically subtract that one overall level, the background correlations collapse to what pure
chance would give. **It is a single, removable, sample-wide level that every reading rides on.**

**The dial is not noise — it is a real property of each sample, and that is why it deserves to be
reported.** When we ask what the dial actually corresponds to, it lines up with concrete things:

- **how pure the sample was** — these are cells teased out of tissue by enzyme digestion, *not*
  purified by sorting, so some samples carry more contaminating non-milk-duct cells than others;
- **how stressed the cells were during preparation** — the warm enzyme digest itself switches on a
  known stress-response programme;
- **how fast the cells were dividing**, and **how strongly the Myc programme was running**.

In other words the dial is **part biology (growth, Myc activity) and part sample-handling (purity,
prep stress) mixed together** — a genuine, measurable, whole-sample state. So we make it a result
in its own right: *the dominant source of variation between these samples is a broad
epithelial-state axis, and most apparent links between pathways are simply different windows onto
that one axis.*

**The consequence for the mitochondrial story.** To ask "are mitochondria *specifically* wired to
growth?" you must first subtract this shared dial. When we do, the mitochondrial programmes have
**almost nothing left** — they *are* the dial, essentially. So there is no specific mitochondria-
to-growth link to report, in either direction; the claim is not weak, it is **unanswerable with
this data**. The one link that does survive subtracting the dial is a non-Myc one — the cells'
**antioxidant (redox) capacity** versus the tumour-like state — and that is the lead worth taking
to the bench.

**One-sentence version for the paper:**
> Across samples, most pathway activities move together because they share a single broad
> whole-sample state (reflecting epithelial purity, preparation stress, and proliferative/Myc
> activity); once this shared state is accounted for, the mitochondrial pathways retain essentially
> no independent covariation, so we do not interpret their raw correlations as specific coupling.

---

## 3. How to write the cell-death priming decline

**The finding you can rely on is the phenotype, and it is external to this dataset.** Independent
experiments (immunohistochemistry, the inducible model) show Myc-driven cells **die at 6 weeks but
not at 12 weeks**. That stands, and the paper explores it experimentally. Keep it.

**What the RNA-seq does NOT support — do not write these as findings:**

- **"The transcriptome shows apoptotic priming rising and then falling" is not established.** The
  "priming" score is built **only from mitochondrial genes**, so any statement that "the
  mitochondrial state couples to death priming" is mitochondrial genes correlating with
  mitochondrial genes — circular by construction. In the rebuilt analysis (script 36) this pair is
  the **negative control**, and it behaves exactly like noise (its sign flips depending on the
  scoring method).
- **The proper statistical test was actually run, and it is not significant.** Whether priming
  changes *differently* over time in Myc versus normal mice was fitted and saved
  (interaction p = 0.18 for the pro-death set, 0.26 for the priming score). The earlier narrative
  quoted group averages ("~3x more at 6W") instead of this test. Use the test: it is null.
- **The apparent decline is not death-specific.** Everything in the Myc programme fades from 6W to
  12W together; the pro-death components fade with it, no faster than any other part. That is the
  textbook picture (as the oncogenic drive relaxes, its built-in death-priming relaxes too) — not a
  dedicated death switch that this data uncovered.

**How to phrase it (suggested):**
> Myc-expressing cells are eliminated at the early timepoint but tolerated later (Fig X, IHC). The
> transcriptome is consistent with this: pro-apoptotic components of the Myc programme are elevated
> early and decline as the programme as a whole attenuates. Because our apoptotic-priming readout is
> defined from mitochondrial genes and does not change significantly in a formal genotype-by-time
> test, we treat the priming decline as a hypothesis to be tested directly rather than as a
> transcriptomic result, and we test it by BH3 profiling (Fig Y).

**The experiment that settles it:** **BH3 profiling** of the 6-week versus 12-week substrate. It
measures death-readiness directly, sidesteps the circularity, and is independent of the
correlation ceiling.

---

## 4. How to write the mtDNA / mitonuclear imbalance

**The solid mitochondrial claims are genotype comparisons — keep these:**

- **Myc raises mitochondrial content ~21-27%** (adjustment-robust, two methods).
- **Myc raises absolute nuclear-encoded OXPHOS.** So Myc **reallocates and expands** the
  mitochondrial compartment; it does not "reduce OXPHOS."

**What is NOT supported — do not write it as a time-course finding:**

- **A mitonuclear imbalance that develops over time** (nuclear OXPHOS genes falling out of step
  with the mitochondrial-DNA-encoded genes) has **no supported statistical test**: the imbalance
  measure shows no significant change by time (p = 0.31), by genotype (p = 0.52), or in the
  interaction (p = 0.92). It looked cleaner than other measures only because it is a **ratio**, and
  a ratio automatically cancels the shared whole-sample dial from §2 — that is an artefact of the
  arithmetic, not a real signal.
- **What the mtDNA-fraction measure even represents is unresolved.** It is entangled with the same
  sample-purity / preparation-stress axis, and we cannot separate a technical cause from a
  biological one without dissociation-batch or viability information, which was never recorded.
  (It is *not* simple contamination — the arithmetic rules that out — but "not contamination" is
  not the same as "a genuine mitochondrial imbalance.")

**How to phrase it (suggested):**
> Myc increases mitochondrial content and absolute nuclear-encoded OXPHOS (Fig X), consistent with
> compartment expansion. We do not claim a time-dependent mitonuclear imbalance: the relevant
> ratio does not change significantly across timepoints or genotypes, and the mtDNA-fraction signal
> is confounded with sample composition and preparation stress in a way this cross-sectional bulk
> design cannot resolve. A direct test is mtDNA copy-number qPCR (Fig Y).

**The experiment that settles it:** **mtDNA qPCR** (copy number), which measures the mitochondrial-
to-nuclear balance directly and does not depend on the whole-sample dial.

---

## 5. Safe to claim vs. needs a bench experiment

| Claim | Verdict | Basis |
|---|---|---|
| Myc raises mitochondrial content ~21-27% | **Safe** | genotype comparison, two methods |
| Myc raises absolute nuclear OXPHOS; expands/reallocates the compartment | **Safe** | genotype comparison |
| The 6W->12W attenuation (66% fade + 34% convergence) | **Safe** | genotype x time contrast |
| Myc amplifies the pubertal programme; suppresses basal/myoepithelial | **Safe** | genotype comparison |
| Myc bends off the wild-type developmental trajectory | **Safe** | gene-level vectors, not a sample coupling |
| Biogenesis is a Myc-dose bystander | **Safe (a negative)** | tracks nothing once Myc is removed |
| The broad whole-sample state axis (purity + prep + growth) | **Report as a phenotype** | §2 |
| "OXPHOS is central to phenotype" | **Drop** | not testable — the axis *is* the dial |
| Death-priming *decline* as a transcriptomic mechanism | **Do not claim** | circular; formal test null; test by BH3 profiling |
| Mitonuclear imbalance as a *time* finding | **Do not claim** | no supported contrast; test by mtDNA qPCR |
| Bbc3/PUMA drives the death difference | **Do not claim** | 1 of 23 = chance; test by BH3 profiling |
| Myc kills at 6W, not 12W (the death phenotype) | **Safe** | external IHC / inducible model |
| Lower redox capacity tracks the tumour-like fork (Myc-independent) | **Lead** | survives the correction across 3 methods; test with a redox/glutathione assay |

**Three bench experiments close almost everything open:** mtDNA qPCR (the imbalance), BH3 profiling
(the death priming, and Bbc3), and a mitochondrial-mass measurement such as a TOMM20/VDAC/citrate-
synthase blot (the content claim). A redox/glutathione assay would additionally test the one new
lead.
