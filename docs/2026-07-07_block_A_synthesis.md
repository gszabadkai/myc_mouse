# Block A — full synthesis (what was tested, the outcome, where it lives)

**Date:** 2026-07-07. **Branch:** `analysis-exploratory`, HEAD `163a53a`.
**Purpose:** one transparent picture of the whole Block A exploratory phase before
we choose figures. Supersedes the terse inventory in
`docs/2026-07-06_block_A_review.md` (that doc's tiered claim-spine and narrative
sections are still valid; this doc is the complete hypothesis -> outcome -> file map).

Read order if you want the argument, not just the table: **§1** (the question) ->
**§2** (the lenses — the part that has caused confusion) -> **§3** (the master
table) -> **§5** (the model) -> **§6** (open figure decisions).

---

## 1. The experiment and the core question

MMTV-Myc transgenic mice, bulk RNA-seq. 4 groups x 6 reps: `6W_neg`, `6W_pos`,
`12W_neg`, `12W_pos` (24 samples). DESeq2 `~ timepoint * myc_status`. Myc dose is
**constant** 6W->12W (mRNA + blot + literature) — a constant driver on a moving
substrate.

Two questions drive Block A:
1. **What is the Myc transcriptional program, and why does it centre on
   mitochondria?** (the "integrate oncogenic + metabolic programs" thesis)
2. **Why does Myc kill more at 6W than 12W?** — the phenotype comes from an
   EXTERNAL Myc-ER tamoxifen-inducible model (IHC, no RNA-seq here): an acute Myc
   pulse on a normally-developed background kills more at 6W. So the death timing
   is a property of the **developmental SUBSTRATE**, not of chronic-Myc adaptation.
   Our chronic RNA-seq characterises that 6W-vs-12W substrate.

**The honest ceiling (applies everywhere):** n=6/group, bulk RNA (not single-cell),
one timepoint per age, and SURVIVOR BIAS (we sequence the cells that did NOT die).
=> we can characterise the death-permissive STATE (association) but not prove
causation. Powered claims are LEVELS / main effects; the interaction (n=6) is
directional support, not significance.

---

## 2. The method lenses — what each one actually answers

This is the single biggest source of confusion, so it is spelled out. The same
mitochondrial signal reads differently through each lens because each normalises
differently. **They are complementary, not interchangeable.**

| Lens | Scale / normalisation | Question it answers | Where |
|---|---|---|---|
| **fGSEA** | Wald `stat` ranks vs the WHOLE transcriptome (shrinkage-independent) | How enriched / important is a pathway **relative to everything else changing**? | `fgsea_*`, `fgsea_percategory` |
| **MitoPathway abundance** | absolute pathway score WITHIN mito (geomean of normalised counts) | How much of the mito budget (by level) sits in a pathway? | `mitopps_scores$raw_*` |
| **MitoPPS** | pairwise ratio WITHIN the mito compartment (Monzel 2025) | How is priority **reallocated among mito pathways** (fixed budget)? | `mitopps_scores$mitopps_*` |
| **GSVA** | per-sample, cohort-relative (log VST, Gaussian) | Where does each SAMPLE sit on a program? (reframing, adds no power) | `gsva_scores`, `gsva_overview` |

**Concepts that are easy to conflate (they are DIFFERENT axes):**
- **Attenuation** (Gate 1) = the Myc effect MAGNITUDE shrinks 6W->12W. A size
  effect, LFC-independent.
- **Reprioritisation** (script 22 AP-abund) = the SPREAD of mitoPPS reallocation
  across pathways narrows 6W->12W (6W selective -> 12W uniform). fGSEA magnitude is
  flat here — a distinct axis from attenuation.
- **Imbalance** (the mitonuclear one) = a WITHIN-mito quantity, read on **mitoPPS
  ALONE**: nuclear-encoded OXPHOS vs mtDNA-encoded OXPHOS priority. NOT an fGSEA
  quantity.
- **Dissociation** = where fGSEA (vs transcriptome) disagrees with mitoPPS (within
  mito). A secondary annotation. **Not** the imbalance (this was a mislabel we
  corrected in script 24).

**The two temporal contrasts (do not mix them up):**
- `timepoint_neg` = **WT** 6W->12W = the developmental substrate the acute pulse
  hits. Also written `Temporal_Myc-` / `WT_6W->12W`.
- `timepoint_pos` = **Myc+** 6W->12W = the chronic-Myc window.
- `myc_6W` / `myc_12W` = Myc effect (pos vs neg) AT each age.
- `interaction` = does the Myc effect differ by age (n=6 -> directional only).

**Shrunken vs raw LFC:** fGSEA is unaffected (ranks on the Wald stat). ΔLFC /
averaged-LFC methods use RAW (unshrunken) LFCs (`combined_df_annotated_raw`).
Shrinkage is for MA-plot QC / visualisation only.

---

## 3. Master table — hypothesis / question -> outcome -> files

Grouped by theme. "Verdict": **CONFIRMED** (powered), **DIRECTIONAL** (method-
consistent, not per-set significant at n=6), **REFUTED**, **INFRA** (engine).

### A. Is the Myc effect real, and does it change with age?

| Question | Script | Verdict | Key result | results/.rds | outputs/ |
|---|---|---|---|---|---|
| Does the genotype effect attenuate 6W->12W? (**Gate 1**) | 13 | **CONFIRMED** | Genotype effect SHRINKS ~2-fold 6W->12W; biological (effect-size), not a power artifact | `gate1_divergence_timing` | `gates/` |
| Data load + DESeq2 (interaction + group) + QC | (01-12 trunk) | INFRA | 24 samples; interaction + group models | `dds_int*`, `dds_group*`, `interaction_results`, `group_results`, `combined_df_annotated*` | `deseq_qc/`, `qc/` |

### B. Why mitochondria? (is it the preferential target)

| Question | Script | Verdict | Key result | results/.rds | outputs/ |
|---|---|---|---|---|---|
| Is mito the most preferentially altered compartment vs an expression-matched null? (**AP6.2**) | 21 | **CONFIRMED** | mito/OXPHOS most preferentially altered vs expression x dispersion-matched null; proliferation does NOT survive | `ap6_permutation_null` | `ap6_null/` |
| Per-category fGSEA enrichment (**AP6.1**) | 20 | **CONFIRMED** | MitoCarta NES 2.18 (padj 1e-55); attenuation-via-convergence on temporal contrasts | `fgsea_percategory` | `fgsea_percategory/`, `fgsea/`, `fgsea_cross_sectional/` |
| Does Myc drive the METABRIC biogenesis fork? (**AP7, centrepiece**) | 18 | **CONFIRMED** | Genotype main effect p=0.005; MB2-specific Wilcoxon p=0.008 | `ap7_mb_fork` | `ap7/` |

### C. What IS the Myc mito program? (biogenesis, and which kind)

| Question | Script | Verdict | Key result | results/.rds | outputs/ |
|---|---|---|---|---|---|
| Myc = biogenesis? developmental relationship? (**AP1-4**) | 19 | **CONFIRMED** (biogenesis) | Myc~biogenesis (Felsher rho 0.90, geno p=3.7e-6); suppresses basal (p=0.016); anti-parallel to lineage development | `dev_composition` | `dev_composition/` |
| mtDNA vs nuclear OXPHOS split; fGSEA-vs-mitoPPS lens split; retention (**AP-mtDNA/abund/retention**) | 22 | **CONFIRMED / DIRECTIONAL** | WT mtDNA-OXPHOS up (0.89->1.42)/nuclear down; abund lens-split TIME-DEPENDENT (6W selective SD 0.14 -> 12W uniform 0.10, fGSEA NES flat 2.2); biogenesis retention convergence | `reframe_supp3` | `reframe_supp3/`, `mitopps/`, `mitopps_fgsea/` |
| Per-pathway within-mito imbalance; Myc-vs-ER biogenesis; TF drivers (**review pts 1/2/4d**) | 24 | **CONFIRMED / DIRECTIONAL** | Imbalance=mitoPPS alone; mtnuc peaks 6W_pos (+0.56), inverts 12W; WT nuclear-assembly(6W)->mtDNA-running(12W), Chaperones fall -> **UPR^mt REFUTED**; Myc drives MYC=ER/PGC1a biogenesis (co-opts ESRRA/NRF1/GABPA; ESR1 antagonised) | `biogenesis_discrimination` | `biogenesis_discrimination/` |
| mitoPPS engine (Monzel) | (08 trunk) | INFRA | 142 mito pathways x 24; mtDNA-encoded separated (13 mt-*) | `mitopps_scores`, `mitopps_fgsea_comparison`, `interaction_fgsea_mitopps` | `mitopps/` |

### D. The cell-death question (why kill at 6W not 12W)

| Question | Script | Verdict | Key result | results/.rds | outputs/ |
|---|---|---|---|---|---|
| Does the Apoptosis-PRO module shift? (**Gate 2**) | 14 | **CONFIRMED** | Apoptosis-PRO module coordinately shifts (t p=0.023); death stays main-figure candidate (DECISION_POINT) | `gate2_apoptosis_readout` | `gates/` |
| Directional binomial on raw dLFC (**branch 1**) | 16 | **REFUTED/null** | 52.7%, p=0.27 (raw binomial null) | `cell_death_de_raw`, `combined_df_annotated_raw` | `cell_death_raw/` |
| Tang 15-RCD fGSEA per contrast (**branch 2**) | (12 trunk) | **DIRECTIONAL** | apoptosis modality temporal_pos NES -1.38 (padj 0.016) — death program declines in Myc+ with age | `cell_death_fgsea`, `cell_death_de` | `cell_death/` |
| **H1** BH3:BCL2 / p53 rheostat | 23 | **CONFIRMED** (coupling) | Myc raises priming (PRO-ANTI geno p=0.038; myc_6W PRO p=0.041) ~3x at 6W; p53/ARF NULL (rheostat = BCL2-family) | `death_timing_substrate` | `death_timing/` |
| **H2** mitonuclear imbalance = mechanism | 23 | **CONFIRMED/DIRECTIONAL** | imbalance peaks 6W_pos; coupling to priming r 0.75->0.07 (6W->12W) | `death_timing_substrate` | `death_timing/` |
| **H3** proliferation-death coupling | 23 | **REFUTED as cause** | prolif~death INVARIANT (0.68->0.65) -> timing is mito-specific, not generic proliferation | `death_timing_substrate` | `death_timing/` |
| **H4** selection / survivor culling | 23, 25 | **DIRECTIONAL** (weakest) | CV/dispersion of priming + biogenesis narrows 6W->12W, more in Myc+ (sd_pro ratio 0.42 vs 0.60) | `death_timing_substrate`, `developmental_substrate_death` | `death_timing/`, `dev_substrate_death/` |
| "WT substrate de-primes with age" | 23 | **REFUTED** | WT baseline priming FLAT (-0.24 -> -0.08). The signal is Myc's COUPLING, not the WT baseline | `death_timing_substrate` | — |

### E. The developmental "why" (the substrate the acute pulse hits)

| Question | Script | Verdict | Key result | results/.rds | outputs/ |
|---|---|---|---|---|---|
| WT developmental substrate shift; which axis tracks death; mtDNA resolution; culling (**review pts 3/4a/4b/4c/4e**) | 25 | **CONFIRMED (mito) / DIRECTIONAL (composition)** | Substrate matures on 2 coupled axes: MITO (imbalance resolution ~88% mtDNA-driven, genotype-shared) + CELL-STATE (lineage SHIFT underpowered p>0.18, BUT per-sample stem/MASC content death-coupled 0.78 at 6W). Part B: bio_comp 0.83 > MASC 0.78 > imbalance 0.67 couple to death, all decouple by 12W | `developmental_substrate_death` | `dev_substrate_death/` |
| Where "ER/PGC1a kills" (external) lands (**review pt 4d**) | 24 | **DIRECTIONAL** | Balanced CORE/MITO_NU biogenesis-death intersections most death-coupled at 6W (r 0.77/0.83 > MYC_SPECIFIC 0.66) | `biogenesis_discrimination` | `biogenesis_discrimination/` |

### F. Infrastructure / engines

| Script | Role | results/.rds | outputs/ |
|---|---|---|---|
| 15 | GSVA engine: 884 library sets x 24 (VST, cohort-relative) | `gsva_scores` | — |
| 17 | GSVA overview + 4-method validation (camera/roast/perm/mechanism) | `gsva_overview` | `gsva_overview/` |

---

## 4. What is powered vs what is directional (state this in the paper)

- **Powered / significant (levels + main effects):** Gate 1 attenuation; AP6.2
  preferential-mito null; AP6.1 MitoCarta NES; AP7 MB-fork (p=0.005); Myc=biogenesis
  (p=3.7e-6); basal suppression (p=0.016); Gate 2 apoptosis-PRO shift (p=0.023);
  Myc priming (p=0.038/0.041); mtDNA developmental resolution (genotype-shared).
- **Directional (n=6 interaction / within-timepoint correlation / survivor-biased):**
  all the death-COUPLING correlations (r 0.75->0.07 etc.); H4 culling; branch-2
  temporal NES; the lineage-composition SHIFT; the Cat-9 coupling ranks.
- **Refuted:** WT de-priming; p53/ARF rheostat; H3 as the timing cause; branch-1
  binomial; the UPR^mt/ISR frame (Chaperones fall, not retained).

---

## 5. The model this assembles (narrative)

**Lead = A** (Myc integrates/preferentially amplifies mitochondrial biogenesis),
**with D as the mito-death bridge**.

1. Myc drives a broad mitochondrial biogenesis program (translation + OXPHOS +
   import all up on every lens at a fixed age), preferentially over the rest of the
   transcriptome (AP6/7). It co-opts the ER/PGC1a-axis machinery (ESRRA/NRF1/GABPA)
   rather than running a separable Myc-only program; it antagonises the ESR1/ER
   program.
2. This biogenesis is UNBALANCED against the developmental substrate. At **6W** the
   WT mitochondrion is in a nuclear-led ASSEMBLY state (nuclear respiratory
   complexes / import / chaperones prioritised, mtDNA-encoded OXPHOS lagging). Myc
   pushes nuclear OXPHOS higher while holding mtDNA low -> the **mitonuclear
   imbalance peaks at 6W_pos**.
3. That imbalanced, biogenesis-primed, stem/progenitor-rich 6W state is
   **death-permissive**: biogenesis, stem content, and the imbalance are all
   coupled to pro-apoptotic priming at 6W (oncogene-induced apoptosis), and Myc
   engages the death program (priming p=0.038, ~3x stronger at 6W).
4. By **12W** the substrate has MATURED — mtDNA-encoded OXPHOS catches up (~88% of
   the imbalance resolution, genotype-shared = developmental), the tissue is more
   differentiated, and priming/biogenesis have converged. The death coupling
   collapses (r 0.75->0.07) and Myc's killing fades (Gate 2). **The 12W decoupling
   IS the loss of killing.**
5. This is why the acute Myc-ER pulse kills more at 6W: it hits a young, imbalanced,
   stem-rich substrate; by 12W that substrate no longer supports the mismatch.

---

## 6. Open decisions before figures

1. **Framing** — confirm A-lead + D-bridge (this doc assumes it).
2. **Which main figure hosts the mitonuclear-imbalance panel** (the mito->death
   bridge). Candidate source: `dev_substrate_death/A1_mito_maturation_imbalance` +
   `biogenesis_discrimination/A2_mitonuclear_imbalance`.
3. **Cell death: main vs Supp.** Review doc recommends HYBRID — imbalance panel
   main, full death analysis (H1-H4 + branches + Gate 2) a strong Supp.
4. **Figure 1 / Figure 2 / Supp content split** — see the draft split in
   `docs/2026-07-06_block_A_review.md` §"Proposed figure split", to be updated for
   the 24/25 results (per-pathway imbalance map; Myc-vs-ER lanes; stem-death
   coupling).
5. **Then** create `paper-figures` off `163a53a` and build scripts 26/27/28.

---

## 7. Housekeeping notes

- **Stale outputs removed (2026-07-07):** the superseded first-run Part-A PDFs of
  script 24 (`A1_tier_three_lens_overview`, `A2_lens_divergence_scatter`,
  `A3_pathway_detail_mycpos_window`) were deleted; the current Part A produces
  `A1_imbalance_group_map`, `A2_mitonuclear_imbalance`, `A3_wt_substrate_shift`.
- `results/*.rds` and `outputs/*` are gitignored (regenerated by sourcing the
  scripts in Positron). The scripts + docs are the version-controlled record.
- Other `results/.rds` not in the tables above are trunk/DESeq intermediates
  (`coldata`, `count_matrix`, `mito_de`, `mitocarta_*`, `ortholog_table`,
  `lfc_comparison_timepoint`, `heatmap_paths_by_variant`, `extended_qc_summary`,
  `interaction_gene_characterisation`, `gene_sets_list`, `fgsea_results`,
  `fgsea_xs_results`) — inputs/engines feeding the above, not standalone results.
