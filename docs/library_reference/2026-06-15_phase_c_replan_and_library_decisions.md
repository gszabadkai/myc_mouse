---
date: 2026-06-15
tags: [project/mammary_geneset_library, planning, phase-c, method-design, nomenclature, decisions]
supersedes:
  - "2026-06-15_phase_c_dataset_review.md (the input-for-replanning note)"
  - "Phase C of PLAN.md (Pal/Bach/Giraddi/Sun shortlist; script 06 MG_*_STAGE_UP/DOWN scheme)"
  - "PLAN.md Categories 6 and 7 (TF_targets; Biogenesis_discrimination), partial"
status: decided - ready to implement in Claude Code
implements_into: PLAN.md v3 (bump from v2 when folded in, mirroring the v1 -> v2 lineage)
---

# 2026-06-15 - Phase C re-plan + library method/TF/nomenclature decisions

This note records decisions taken in a planning session that reviewed the
completed literature table
(`gene_set_resources_table_with_developmental_stages_and_supplementary_information.md`)
against the 2026-06-15 dataset-review note and the current PLAN.md.

It is the DECISIONS layer. Nothing has been written to the library yet. No set
has been named yet. Claude Code implements from this note: it first folds these
decisions into PLAN.md (-> v3), and that edited PLAN.md is the first thing to be
reviewed, before any set is built.

Scripts 00-05 (setup + Phase B: msigdb, mitocarta, ortholog, user_curated,
myc_mouse_assets) are NOT changed by this note. Their code stands. This note
affects script 06 onward, plus the QC script (12), plus the Category 6/7
descriptions and the project scope wording.

---

## Decision 1 - Phase C dataset roster (mouse)

The 2026-06-15 review note proposed two pivots that the completed table reverses:

- It proposed **Pal 2021 (BCR)** as the transcriptome backbone. The table shows
  Pal 2021 has **no directly usable gene set in its supplementary** (Supplementary
  Table S1 is QC statistics only; the TEB-vs-ductal heatmap sets exist only with
  further analysis). So Pal 2021 cannot be a backbone as-is.
- It hedged **Gray 2023** as a likely voter only ("DE tables may be thin"). The
  table shows Gray ships clean, fGSEA-ready DE tables: Suppl2 (HE vs LE subtype
  DE, logFC, padj) and **Suppl4 (TEB vs ductal DE, logFC, padj)** - the direct
  pubertal stage contrast Pal 2021 was wanted for. The note's blocking open
  question is therefore resolved: Gray has usable DE.

Final roster:

- **Backbone (clean DE, puberty + adult):** Gray 2023, Garcia Sola 2021, Scheele 2017.
- **Cross-species bridge:** Saeki 2021.
- **Consensus voters:** Pal 2017 (Suppl3 basal/luminal), Henry 2021.
- **Embryonic / stem slice:** Giraddi 2018, Chung 2019 (ATAC-derived).
- **Pal 2021 (BCR): PARKED.** No re-analysis within Phase C. Its signal still
  enters the library indirectly, because Gray Suppl2 is built on Pal 2021 scRNA
  data and carries it through Gray's DE lists. Cheng 2023 remains the documented
  workflow if its TEB sets are ever derived later (out of current scope).
- **Bach 2017:** no own supplement reviewed; reachable via the Garcia Sola and
  Saeki integrations, both of which reanalyse it.
- **Sun 2018: DROPPED** (mislabelled as pubertal; profiles adult-virgin + pregnant).
- **Garcia Sola 2021 == the "Gutierrez 2021" of the review note** (same DOI
  10.1007/s10911-021-09488-1). It is the integration of Pal 2017 + Giraddi 2018 +
  Bach 2017, so it may make those three redundant as individual downloads - check
  during the inspection pass.

Human-side datasets (Reed 2024, Kumar 2023, Pal 2021 EMBO, Chen 2022) remain
deferred to a separate Phase C-human sub-section for the downstream Hannon analysis.

---

## Decision 2 - file staging and readers (Phase C mechanics)

- **Stage all atlas supplement files now**, into a dedicated subfolder
  `data/raw/mammary_dev_atlases/` (parallel to `data/raw/user_curated/`), using the
  table's filename convention (`Gray-et-al_2023_Suppl1.xlsx`,
  `Saeki-et-al_2021_Suppl4.xlsx`, etc.). Staging now is a prerequisite Claude Code
  needs regardless; it does not help the chat-side planning (the conceptual plan
  needs no column layouts).
- **README for the folder:** the full resources table .md serves as the master
  reference README. It is a superset of folder contents (includes human-only
  Chen 2022, code-only Cheng 2023, parked Pal 2021), so Claude Code prepends a
  short provenance header listing exactly which files are present, with
  download_date, source DOI/accession, and sheet-of-interest per file. DOI,
  accession, and sheet come from the table; the **download_date must be supplied
  by the user** (one line at session start, e.g. "downloaded today, 2026-06-15");
  some exact file URLs may be left as the landing page where the direct link is
  not derivable.
- **Readers are bespoke per file, not one generic function.** The structures
  genuinely differ (Gray: per-lineage worksheets, HE vs LE; Saeki: ranked
  correlation lists with fixed sizes; Garcia Sola: clusters with logFC-ranked
  DEGs across 3 trajectory sheets; Pal 2017: clusters as columns; Giraddi: NMF
  groups). Workflow: Claude Code writes an inspection block (per file: sheet
  names, dims, head); the user runs it in Positron; Claude Code then writes one
  tailored reader per file. All readers emit one tidy schema:
  `{set_name, gene, logFC, padj, source, stage, method}`.
- The sheet-to-biology mapping (which sheet = which contrast; which clusters are
  pubertal; where to truncate) is the interpretive work already done in the
  resources table. Claude Code leans on the table for that mapping rather than
  guessing.

---

## Decision 3 - fGSEA vs ssGSEA/GSVA: scope and method model

Both methods are kept. They are complementary, not redundant.

- **fGSEA -> NES.** Contrast-level, competitive (set vs rest of genome),
  directional. Runs on a ranked list (signed DESeq2 Wald `stat`); the contrast is
  defined upfront in the design. One signed NES per set per contrast. **No
  per-sample value can be extracted.** Answers: "is this program coordinately
  shifted in the MYC contrast?"
- **ssGSEA / GSVA -> per-sample score.** Scores each sample from the expression
  matrix alone; produces a set x sample matrix. Sample-level scores drop straight
  into the **same genotype x timepoint interaction design** used at the gene
  level. Answers: "how does this program's abundance move across the timeline, and
  is the MYC effect stage-dependent?"

**Concordance is the prize.** Where a program exists in both a tight and a large
form, run both; agreement is the strong evidence, disagreement is informative.
This extends the project's existing fGSEA-vs-mitoPPS concordance logic to the
developmental sets.

**Interpretive ceiling (record in manuscript-relevant terms).** A signature
score on BULK data conflates two things: a change in the proportion of cells
running the program (composition) and a change in per-cell expression at fixed
composition (regulation). A GSVA score cannot separate them. The existing
"clonal selection against mito-high apoptosis-primed cells" reframing is a
composition claim; a GSVA score moving the right way is consistent with it but
does not prove it over per-cell regulation. Separating the two needs explicit
deconvolution (CIBERSORTx / MuSiC / Bisque against a reference) or orthogonal
evidence (flow, IHC). GSVA also does not add statistical power: the interaction
test on scores has the same small-n limits as the gene-level interaction.

**Size and DE availability both point the same way** (and the table already
pre-sorted on this):

- fGSEA needs 15 <= size <= 500. GSVA tolerates large sets (built for them). The
  asymmetry is at the top end. <15 is unusable for either.
- Sets tagged in the table "too many genes for fGSEA, quantify for average
  expression" (Giraddi ~477-937, Saeki scGSVA 160/240/500/200, Gray major-type)
  -> GSVA. Sets tagged "fGSEA -> NES" (Gray Suppl2/4, Garcia Sola clusters,
  Pal Suppl3) -> in-range, fGSEA.
- A set used in fGSEA does NOT need DE stats on itself (only a membership list +
  the ranked data). But source DE stats matter twice: for size control (logFC +
  padj let you threshold at FDR<0.05 and truncate top-N -> tight, in-range) and
  for sign (signed UP/DOWN sets are far more informative under a signed ranking).
  Source DE present -> tight, signed -> fGSEA. No source stats -> used whole -> GSVA.

**Technical cautions (carry into the GSVA-scoring script):**

1. **GSVA input scale is OPPOSITE to mitoPPS.** mitoPPS uses DESeq2-normalised
   counts on the LINEAR scale (not VST). GSVA wants the LOG scale: log-CPM or
   VST/normalised-log with `kcdf="Gaussian"`, OR raw integer counts with
   `kcdf="Poisson"`. Do NOT carry the mitoPPS scale choice over.
2. **GSVA scores are cohort-relative.** The nonparametric step uses the
   distribution across the samples in the call. Score ALL samples together in one
   run; scores are not portable to a separately-scored cohort. Fine for the
   single-matrix interaction model; never score in two batches and compare.
3. GSVA is preferred over raw mean-expression for the large sets (mean-of-scaled-
   genes is dominated by a few high-expressors); keep mean-expression only as a
   sanity check.

**Where each sits in the existing stack:** mitoPPS is already a per-sample scorer,
specialised for mitochondrial pathways. GSVA is the general per-sample scorer for
the developmental/lineage signatures. fGSEA is the contrast-level competitive
layer over both.

**Plan consequence:** every set carries a **method tag** (`fgsea` / `gsva` /
`both`). Script 12's 15-500 size filter needs a **per-set exemption flag** so the
GSVA-tagged large sets survive the QC step.

### Method assignment (which set, which method, what it buys)

| Set group | Method | What it achieves |
|---|---|---|
| MitoCarta sub-pathways, MYC sigs, Metabolism, Proliferation, generic TF targets, Apoptosis, intersection categories | fgsea | Directional "is program shifted in the MYC contrast"; continuity with myc_mouse |
| Gray Suppl2 (HE/LE), Gray Suppl4 (TEB-vs-ductal UP/DN), Garcia Sola clusters (truncated), Pal 2017 Suppl3, Scheele clusters, Henry mEC/mNEC | fgsea | Tests whether MYC coordinately shifts a named developmental program in a contrast |
| Giraddi fMaSC/balancers, Saeki scGSVA Stem/Basal/Alv/Hor, Gray major-type (BA/HS/AP) | gsva | Per-sample developmental-state abundance -> genotype x timepoint interaction |
| Any program present in both a tight and a large form | both | Concordance test (agreement = strong; disagreement = composition-vs-regulation flag) |

---

## Decision 4 - TF lanes, biogenesis framing, Category 7 rescope

**Overarching aim:** distinguish developmental from MYC-driven mitochondrial
transcriptional regulation, at hypothesis-generation level. GRN analysis
(SCENIC/pySCENIC, decoupleR) is the eventual confirmation and is **explicitly out
of scope** for this paper. Gray's CHEA output is TF->target inference filtered by
a mammary cell-state signature (not mammary ChIP), which keeps it hypothesis-level
and consistent with that scope.

**Corrections / framing fixed in this session:**

- **ESR1 != ESRRA.** The original TF roster (NRF1, GABPA, ESRRA, PPARGC1A, MYC)
  had no ligand-activated estrogen receptor. ESRRA is the orphan
  estrogen-RELATED receptor, PGC-1a's partner, on the biogenesis side. ESR1 is the
  estradiol-activated developmental signal. Add ESR1.
- **PGC1a is a coactivator, not DNA-binding.** It does not appear in CHEA as
  itself; its regulon is operationalised through its partners. "PGC1a axis" =
  `union(ESRRA, NRF1, GABPA)` at the DNA level. This is why the earlier
  "PGC1a vs MYC" framing and the core-TF framing are the same object.

**Developmental arm - ESR1-only in pass 1.** The mouse model is 6-12 weeks
(young adult; active estrous cycle; NO pregnancy). ESR1 is the dominant hormonal
input in this window. PGR (alveologenesis) and STAT5 (lactation) are
pregnancy-driven and are **deferred to a second pass**, with GATA3/FOXA1/ELF5
considered then if Gray surfaces them for HS/AP.

**USER DECISION - keep ESR1 inside the mito algebra (not developmental-GSVA only).**
The literature assumption that ESR1 acts on the mitochondrial transcriptome
exclusively via PGC1a is not accepted here. ESR1 has a wide cistrome; non-canonical
routes may exist and are worth discovering. So `DEVELOPMENTAL_MITO = ESR1 ∩
MitoCarta` is retained as its own arm, with `ESR1_NOT_CORE_MITO` as the discovery
set: ER-reachable mito genes the PGC1a axis does not explain. A near-empty result
is itself evidence for the exclusivity model (derived, not assumed).

  Caveat to record: ER's broad cistrome makes `DEVELOPMENTAL_MITO` large and
  noisier than the core sets, with non-functional occupancy included. In fGSEA it
  is a blunter instrument; read its enrichment as "ER-associated", not
  "ER-driven". Mitigated by the Gray-CHEA ESR1 context sets (sharper, mammary-
  filtered) used for GSVA - the same TF seen both genome-wide (fGSEA, generic
  regulon) and in-context (GSVA, Gray).

**Category 7 is RESCOPED, not added to.** "MYC vs PGC1a" -> "oncogenic vs
developmental biogenesis". No new top-level category; the category count stays at 9.

### Pass-1 set algebra (script 08), all gated size >= 15

```
MITO = MitoCarta (full ~1140 genes - see Decision 5)
generic *_targets supplied by the generic regulon lane (DoRothEA A/B + ChIP-Atlas)

DEVELOPMENTAL_MITO          = intersect(ESR1_targets, MITO)
CORE_MITO                   = intersect(union(ESRRA, NRF1, GABPA)_targets, MITO)   # PGC1a axis
MYC_MITO                    = intersect(MYC_targets, MITO)

MYC_SPECIFIC_MITO           = setdiff(MYC_MITO, union(DEVELOPMENTAL_MITO, CORE_MITO))   # key set
ESR1_NOT_CORE_MITO          = setdiff(DEVELOPMENTAL_MITO, CORE_MITO)                    # discovery set
ESR1_AND_CORE_MITO          = intersect(DEVELOPMENTAL_MITO, CORE_MITO)                  # canonical route
MYC_AND_DEVELOPMENTAL_MITO  = intersect(MYC_MITO, DEVELOPMENTAL_MITO)
MYC_AND_CORE_MITO           = intersect(MYC_MITO, CORE_MITO)                            # hijack set
MYC_DEV_CORE_TRIPLE         = intersect(intersect(MYC_MITO, DEVELOPMENTAL_MITO), CORE_MITO)
```

### Two TF lanes (Category 6) - side by side, do NOT merge

Same object type (TF->target lists), different methods:

| Lane | Source | Used for | Method |
|---|---|---|---|
| Generic regulon | DoRothEA A/B, ChIP-Atlas | Category 7 set algebra (∩ MitoCarta) | fgsea |
| Mammary-context | Gray Suppl3 (adult HE/LE), Suppl5 (TEB/ductal) | per-sample TF-program scoring | gsva |

Generic regulons are complete -> correct for set operations. Gray-CHEA is mammary-
and stage-resolved -> correct granularity for abundance scoring.

TF-set naming carries a provenance token so both coexist:
`TFT_<TF>_DOROTHEA_AB`, `TFT_<TF>_CHIPATLAS` (generic);
`TFT_<TF>_GRAY_<CONTEXT>` (Gray), CONTEXT in
{HS_HE, HS_LE, AP_HE, AP_LE, BA_HE, BA_LE, TEB, DUCTAL}.

### Hypotheses this enables (both tested on the genotype x timepoint interaction model)

- **H1 (convergent):** MYC hijacks the developmental TF program. GSVA - the MYC
  timeline scores high on Gray's proliferative-expansion / TEB program; the
  ribosome-biogenesis signature rises. (The parallel the dataset-review note
  already spotted, now testable.)
- **H2 (divergent):** MYC drives mito biogenesis the ER/PGC1a route does not.
  fGSEA - `MYC_SPECIFIC_MITO` enriched; GSVA - developmental (ESR1) programs flat
  or suppressed while the MYC-specific mito program rises.

Because developmental TF programs have known temporal windows (ER/branching at
puberty), scoring them across the MYC timeline asks the sharp question: does MYC's
mito effect TRACK a developmental window's TF program, or is it
timepoint-independent and cell-autonomous? That is an interaction on the GSVA
scores.

---

## Decision 5 - Gray-CHEA mito-overlap TF search (script 07 add-on)

For each Gray CHEA TF, intersect its listed overlapping signature genes with
MitoCarta, to find TFs that regulate mitochondrial genes present in each
cell-state program. Numbers are small -> treat as qualitative shortlisting, not
ranking.

**USER DECISION - keep the MitoCarta list LONG (full ~1140 genes), not the
biogenesis subset.** The project wants to see how apoptosis (and metabolism) are
regulated, so TFs acting on the apoptotic and metabolic mito aspects are wanted,
not only biogenesis regulators. The search therefore finds "TFs regulating any
mito program present in the cell-state signature" - directly relevant to the
apoptosis-evasion thread.

**Refinement - annotate by MitoCarta sub-pathway, do not just count.** A raw count
collapses a TF that hits 8 OXPHOS genes and one that hits 8 apoptosis genes into
the same number. The shortlist artifact carries a breakdown so the list can be
read functionally:

```
TF | n_mito_total | n_oxphos | n_mitoribo | n_apoptosis_pro | n_apoptosis_anti | n_metab_other | frac | hyperg_p
```

Sub-pathway membership comes from script 02's existing MitoCarta splits; the
search just tags each overlap gene with its pathway. A per-aspect cap (top ~5-8
by sub-pathway) is used rather than one global threshold, or the list grows
unwieldy with the long MitoCarta.

- **Validation:** the search should recover ESRRA / NRF1 / GABPA from Gray's own
  CHEA output (recovers the known core in-context). Novel mito-overlap TFs are
  genuine hypothesis candidates -> add to the GSVA panel; optionally emit small
  `TFT_<TF>_GRAY_<CONTEXT>_MITO` sets (GSVA-only / descriptive, given size).
- **Artifact:** `outputs/qc/gray_chea_mito_tf_shortlist.csv` - a stop-and-check the
  user reviews on the real Suppl3/Suppl5 before the TF panel locks.

**TF panel to score (selection rule: present in Gray CHEA AND assignable to one arm):**

- Developmental: ESR1 (pass 1); PGR, STAT5, +/- GATA3/FOXA1/ELF5 (pass 2).
- Biogenesis core: PPARGC1A (as concept), ESRRA, NRF1, GABPA.
- Oncogenic: MYC, E2F family, plus Gray's proliferative-expansion CHEA hits.

---

## Decision 6 - nomenclature

Cell-state names are forced: the published supplements give cell-state / lineage
signatures, not clean per-developmental-stage DE contrasts (the only genuine
contrast is Gray's TEB-vs-ductal). The `MG_*_STAGE_UP/DOWN` scheme in the old
script-06 spec is therefore REPLACED.

**Table names are collection-level, not set-level** ("Gray-Brugge major type
signature set" = 3 sets; "Saeki-Chan SC signature set" = 6). They cannot be the
set name. What makes a row unique is **cell-state + source**. Consensus filtering
needs the source token (e.g. `MG_BASAL_GRAY` vs `MG_BASAL_PAL2017` vs
`MG_BASAL_SAEKI` is the whole point of voting).

### Grammar

```
set name:    <PREFIX>_<CELLSTATE-or-CONTRAST>_<SOURCE>[_<UP|DN>][_<SUBTOKEN>]
description: "<full table collection name> | <Suppl#> | <stage annotation> | <method>"
method tag:  fgsea | gsva | both   (separate provenance column)

PREFIX:   MG_  (mammary development / cell-state),  TFT_ (TF target)
SOURCE:   short token - GRAY, GARCIASOLA, SCHEELE, SAEKI, PAL2017, HENRY, GIRADDI, CHUNG
SUBTOKEN: disambiguators only - _GSVA or _C (same state via two analyses in one source),
                                _ATAC (Chung, accessibility-derived)
```

- The **full table name stays in the GMT description field + provenance CSV**
  (verbatim), never in the token.
- **Stage goes in the description, not the name.** The same cell-state spans
  stages and the source-to-stage mapping is often fuzzy ("pubertal ~ clusters
  18,19 (?)"). Forcing a stage token would bake uncertain calls into the
  identifier. Keep stage revisable. Only the genuine contrasts carry stage in the
  name (Gray `TEB_VS_DUCTAL`).

### Examples (covering the messy cases)

```
Gray Suppl1   -> MG_BASAL_GRAY, MG_HS_GRAY, MG_AP_GRAY                 (gsva)
Gray Suppl4   -> MG_TEB_VS_DUCTAL_GRAY_UP / _DN                        (fgsea)
Pal 2017 S3   -> MG_BASAL_PAL2017, MG_LUMINAL_PAL2017                  (fgsea)
Saeki Suppl4  -> MG_BPRO_SAEKI ... MG_LHOR_SAEKI                       (fgsea)
Saeki Suppl6  -> MG_STEM_SAEKI_GSVA, MG_BASAL_SAEKI_GSVA, ...          (gsva)  # _GSVA disambiguates vs Suppl4
Chung         -> MG_FETAL_CHUNG_ATAC, ...                              (varies) # _ATAC, not pooled with RNA
Gray CHEA     -> TFT_MYC_GRAY_TEB, TFT_ESR1_GRAY_HS_HE                 (gsva)
```

### Wrinkles to expect on the generated draft (resolve in the readers)

- Within-source collision when one paper gives the same cell-state via two
  analyses (Saeki Suppl4 cluster-markers + Suppl6 scGSVA both "basal") -> add a
  `_GSVA` / `_C` sub-token.
- Chung is ATAC-derived, not RNA -> `_ATAC` token so it is not silently pooled.

### Inspection tiers (how much the user needs to look at the real files)

The convention needs no inspection. The NAMES need a graded peek - and that peek
is the SAME inspection pass used for the readers (Decision 2) and the mito-TF
search (Decision 5). Claude Code generates names programmatically from the grammar
+ the resources table; the user reviews a generated draft, not hand-authored rows.

- **Generate straight, no inspection:** Gray Suppl1/4, Pal 2017 Suppl3, Saeki Suppl6.
- **Light confirm (worksheet split + counts):** Gray Suppl2, Saeki Suppl4, Henry,
  Giraddi.
- **Real look (numbered clusters, partial identity, fuzzy stage):** Garcia Sola
  (~39 numbered clusters across 3 sheets) and Scheele (9 mixed clusters). These
  are where cluster naming and stage annotations actually get decided.

---

## Knock-on edits for Claude Code to make in PLAN.md (-> v3)

1. **Section 1 file list / tree:** remove `pal2017_supplementary.xlsx`,
   `bach2017_supplementary.xlsx`, `giraddi2018_supplementary.xlsx` from
   `data/raw/`; add `data/raw/mammary_dev_atlases/` holding the new roster files
   under the table's filename convention, plus the resources-table README.
2. **Section 1 + Section 0 scope:** add GSVA to scope; add the per-set method tag;
   add the size-filter exemption concept.
3. **Section 0:** add the "oncogenic vs developmental biogenesis" framing and the
   "GRN is out of scope" boundary.
4. **Section 3 category table:** rescope Category 6 (two TF lanes) and Category 7
   ("oncogenic vs developmental biogenesis"). Keep 9 categories.
5. **Step 8 (06_build_mammary_dev_sets.R):** replace the `MG_*_STAGE_UP/DOWN`
   spec with the cell-state/lineage spec + the Decision 6 grammar + the genuine
   contrasts (Gray `TEB_VS_DUCTAL`). Tidy schema
   `{set_name, gene, logFC, padj, source, stage, method}`. Consensus sets
   (`MG_<STATE>_CONSENSUS`, intersection across sources, size >= 15) now vote
   across the new roster.
6. **Step 9 (07_build_tf_target_sets.R):** add ESR1; add the generic vs Gray-CHEA
   two-lane build; add the mito-overlap shortlist with sub-pathway breakdown
   (Decision 5) and its `outputs/qc/gray_chea_mito_tf_shortlist.csv` artifact.
7. **Step 10 (08_build_intersection_sets.R):** widen to the Decision-4 pass-1
   algebra (ESR1 developmental arm + PGC1a core + MYC, with the discovery/hijack
   complements).
8. **Step 16 (12_qc_filter_and_dedup.R):** add the per-set size-exemption flag so
   `gsva`-tagged large sets are not dropped by the 15-500 filter.
9. **Pre-flight checklist:** update the supplementary filenames to the new roster.
10. **Provenance table:** add the `method` column (fgsea/gsva/both) and ensure the
    full table collection name is carried in the description field.

---

## Deferred / out of scope (recorded so it is not a later scramble)

- **Pass 2 TFs:** PGR, STAT5, +/- GATA3/FOXA1/ELF5 - only if the analysis extends
  past the 6-12wk window.
- **Pal 2021 (BCR) TEB re-analysis** via the Cheng 2023 workflow - bigger lift,
  out of Phase C.
- **GRN confirmation** (SCENIC/pySCENIC, decoupleR) - out of scope for this paper.
- **Phase C-human sub-section:** Reed 2024, Kumar 2023, Pal 2021 EMBO, Chen 2022,
  + the Hannon dataset (awaiting raw count matrices).
- **Deconvolution** (CIBERSORTx / MuSiC / Bisque) - the only way to resolve the
  composition-vs-regulation ambiguity in the GSVA scores; flagged, not scheduled.

## Adjacent literature (manuscript ref candidate, not Phase C)

- Mycn 2025 (Cells, GSE251933, Westlake group), "Mycn Is Essential for Pubertal
  Mammary Gland Development". Mycn != c-Myc, but adjacent to a MYC-driven mammary
  tumorigenesis manuscript.

## Citations / DOIs (carried from the review note)

- Pal 2017: 10.1038/s41467-017-01560-x
- Bach 2017: 10.1038/s41467-017-02001-5
- Giraddi 2018: 10.1016/j.celrep.2018.07.025
- Sun 2018: 10.1074/jbc.RA118.002297
- Pal 2021 (BCR): 10.1186/s13058-021-01445-4
- Garcia Sola 2021 (= Gutierrez 2021): 10.1007/s10911-021-09488-1
- Saeki 2021: 10.1038/s42003-021-02201-2
- Gray 2023: 10.1016/j.celrep.2023.113293
- Scheele 2017: 10.1038/nature21046
- Henry 2021: 10.1007/s10911-021-09486-3
- Chung 2019: 10.1016/j.celrep.2019.08.089
- Cheng 2023: 10.12688/f1000research.134078.2
- Chen 2022: 10.1038/s41597-022-01236-2
- Reed 2024: 10.1038/s41588-024-01688-9
- Kumar 2023: Nature 620:181-191 (2023)
- Pal 2021 (EMBO J): EMBO J 40:e107333
