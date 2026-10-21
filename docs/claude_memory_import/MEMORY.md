# Memory index

## The paper's frame and narrative

- [The mitochondrial trade-off model](mitochondrial-tradeoff-model.md) — **the paper's frame**: respiration is ASSET + LIABILITY. Closing claim = not what the tumour NEEDS but what it can AFFORD. The gland de-respires (−0.255) without de-proliferating
- [Results order: SETTLED, Order D](results-order-three-narratives.md) — final narrative in `docs/2026-08-02_..._order_D.md`; §1 oncogene, §2 tissue. Three corrections to the author's spine inside
- [Dose vs competence](dose-vs-competence.md) — **the spine of the death story**: MMTV-Myc measures DOSE (protein −50% = ×0.55), MYC-ER measures COMPETENCE. Myc builds a death-ready mitochondrion at both ages and loses the TRIGGER
- [Foxo3/PUMA pair, non-linear de-dosing](foxo3-puma-pair-and-the-nonlinear-dedosing.md) — Foxo3 rank 1, Bbc3 rank 2 of 26 by interaction p. **Both genotype halves are mirror images and NEITHER is significant** — neither half may be labelled on a panel
- [OXPHOS/apoptosis module: the verdict](oxphos-apoptosis-module-verdict.md) — **no co-regulated module; PUMA IS A SOLO.** MYC's regulon owns the apoptotic genes, not PGC1a. What collapses is luminal lineage IDENTITY
- [Hypothesis-generating, not failed tests](hypothesis-generating-not-failed-tests.md) — a number that fails its null still nominates a testable proposition; present competing readings with the experiment that adjudicates each
- [Pre-specify from the perturbations](prespecify-from-the-perturbations.md) — name genes/sets in advance from the cell experiments; report the pre-specified items that FAIL too
- [Myc death-timing question](myc-death-timing-question.md) — why Myc kills at 6W not 12W; the deferred integrated analysis
- [Intro OXPHOS-restraint tension](intro-oxphos-restraint-tension.md) — resolved by AXIS + the SUBSTRATE frame. WT OXPHOS −0.255 vs Myc+ −0.405: **the only mito arm the normal gland withdraws from**

## Figures and instruments

- [Figure-phase architecture](figure-phase-architecture.md) — the `figures/panels/` layer, `theme_panel()`/`panel_legend()`/`legends.md`, three circulation PDFs. **Block B FROZEN at `block-b-full`; BLOCK C on `paper-final`, 28 panels, completeness glob-driven**
- [The normalisation ladder has a floor](normalisation-ladder-has-a-floor.md) — instrument choice: the message needs rungs 0–2 (counts/share/log2FC), the corpus drew 3–5 (mitoPPS/GSVA/NES). **Script 46 found the floor** — absolute OXPHOS level is 93% the global factor, so **KEEP mitoPPS for Fig. 2I** and say respiratory PRIORITY, not capacity
- [fGSEA takes ONE list, not two](fgsea-takes-one-list-not-two.md) — why a rank-stable Myc effect gives NEGATIVE timeline NES: four lists; NES is scale-free; rank stability **forces** the mirror. "The ranking" is over gene SETS (0.933), not genes (0.724/0.520)
- [Main panels must carry the argument](main-figure-panels-must-carry-the-argument.md) — *"the figures have to talk by themselves, people do not read the text."* One grammar per figure; name the drawn elements
- [The two-timeline grammar](two-timeline-grammar.md) — the plane 2F/2G/S2D alternatives share; **the drop from the diagonal IS the interaction**. 2H+I left the plane 2026-08-09
- [Use the defined scheme, don't invent one](use-the-defined-scheme-dont-invent-one.md) — declared palette/labels/idiom; a near-miss reads as an error. Declare anything new in `_panel_common.R`
- [Panels follow the narrative](panels-follow-the-narrative.md) — build a panel when its Results sentence exists, not when its data exists

## What the data does and does not support

- [Mito content is bounded, not measured](mito-content-bounded.md) — **Myc raises content +21–27%, p=0.0006 — the one solid new result**, and a LOWER bound. Several script-32 claims retracted
- [Abundance rulers disagree](abundance-ruler-dependence.md) — the same 87 OXPHOS subunits give **+0.061 / +0.201 / +0.226** by weighting; name the ruler. "Returns to baseline" is compartment-relative; five lowest expressers are tissue-restricted paralogs
- [Background vs Myc is additive](background-vs-myc-is-additive.md) — **never flip the shared-vector axes.** Myc RESCALES ×0.55 without reshaping; Issue #6's 34% WT-convergence is a shared-baseline ARTIFACT
- [Issues #4/#5/#6 attenuation done](issue4-attenuation-scoped.md) — magnitude not rank shrinks; no Myc-independent axis absorbs it; attenuation = ~66% Myc-fade + ~34% WT-convergence (biosynthetic arms only)
- [Attenuation: four competing hypotheses](attenuation-competing-hypotheses.md) — RESOLVED by the blot; kept for the negative scans (every E-box binder flat; 111-gene NR/corepressor/DREAM scan null)
- [The PGC1a axis did NOT change](pgc1a-axis-did-not-change.md) — the developmental OXPHOS decline is arm-selective inside a flat regulon (0.05th vs 99.55th pct). **Never name a TF from a `_GRAY_`/`_MITO` lane** — context beats factor 0.472 to 0.267. Carries the measured licence for the PGC1a/NRF1 rescue
- [Mitonuclear imbalance is UNTESTED](mitonuclear-imbalance-untested.md) — no supported contrast; the mito→death bridge is retracted as a mechanism. What SURVIVES is listed
- [The MYC x OXPHOS gate model](myc-oxphos-priming-gate-model.md) — **the model for the human arm, written as an equation and verified in the mouse (script 48)**. Pure gate, not a crossover; use `ox_rel`; the guardian half carries it — `MCL1:BCL2L1` beats `PUMA:BCL-XL`
- [Death priming: no transcriptomic support](death-priming-change-untested.md) — `pro_comp` IS a MitoCarta set (circular); the priming interaction was fitted, saved and NEVER READ; the attenuation is GLOBAL
- [The correlation ceiling](correlation-ceiling.md) — every composite correlates with every other at ~0.8. **ONE global per-sample axis** explains it. Coupling instrument = linear z-score + design→PC1→residual, **not GSVA**. Trap: mito axes ARE the global factor
- [Composition adjustment needs a positive control](composition-adjustment-needs-a-positive-control.md) — a loading belongs to a SUBSET not a ruler; and an adjustment that removes the confound may remove the biology (Mki67 +0.31 -> -0.10). **Flagged for the human arm's trap list**
- [Sample-level couplings need a null](sample-level-couplings-need-a-null.md) — at n=12 everything correlates at 0.6–0.8; script 25's death couplings fail an empirical null. MECs are enzymatically dissociated, NOT FACS
- [Gene PC1 = composition, PC2 = Myc-mito](gene-pc1-is-composition-pc2-is-myc-mito.md) — gene-PC1 (44%) is a composition nuisance, correctly dropped; the Myc-mito programme rides at PC2 and the pathway lens foregrounds it as the 76% axis
- [Figure B loaders are curated, not top](figure-b-loaders-are-curated-not-top.md) — axis is mito-LED not mito-SPECIFIC; **loading = covariation NOT primacy.** The very top is construction-inflated
- [mt% is three-way-confounded](mtpct-is-three-way-confounded.md) — content × proliferation denominator × prep-leak, unresolvable. Only nuclear-OXPHOS↔priming is mt%-robust → BH3
- [MitoCarta sets are membership-loose](mitocarta-sets-are-membership-loose.md) — never make a set-level MECHANISM claim without resolving gene by gene; small set + big effect = one gene, not a module
- [LE = low-expressing, not low-oestrogen](le-is-low-expressing.md) — Gray 2023 p.7; our provenance table is WRONG. 25 of 87 OXPHOS subunits are AP LE-markers → always report the mito-removed version

## Design, data traps and method

- [Batch = timepoint](batch-equals-timepoint.md) — batch confounded with timepoint; genotype balanced within each. **Genotype/interaction CLEAN; pure-time reading confounded and not separable**
- [Timepoint is batch-aligned — but depth is NOT a confound](timepoint-cohort-confound.md) — shares and DESeq2 both cancel depth; 12W is not degraded. mt temporal claims are IMPRECISE, not uninterpretable
- [Gene-symbol vintage reconciliation](gene-symbol-vintage-reconciliation.md) — old GMT symbols silently drop renamed genes during SET aggregation. **Any new set-based analysis MUST use `functions/reconcile_gene_symbols.R`**
- [Don't assert a mechanism without checking it can work](dont-assert-a-mechanism-without-checking-it-can-work.md) — four failures in one session, all caught by the author's plain questions. **When the author asks a simple question about a result I'm confident in — TEST, don't explain**
- [MEC consensus nomenclature](mec-consensus-nomenclature.md) — Gray/Kessenbrock/Khaled 2025 3-type (BMYO/LASP/LHS), adopted project-wide
- [Block A Day-1 gates](block-a-day1-gates.md) — gates, the GSVA engine, trajectory validation, biogenesis = developmental convergence

## External literature

- [Lesner mtDNA external validation](lesner-mtdna-external-validation.md) — bioRxiv 2026.07.13.738248: TCGA + mouse HCC validate the nuclear-up/mtDNA-down discordance (kills the dissociation-artifact worry) and endorse the qPCR/BH3 plan. Validates the QUESTION, not our contrast
- [Lesner mechanisms tested (script 39)](lesner-mechanisms-tested-script39.md) — hypoxia→MXI1 route does NOT transfer (liver-specific); Myc raises mito-turnover but not mitophagy-specifically
