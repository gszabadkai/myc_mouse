---
date: 2026-06-15
tags: [project/mammary_geneset_library, planning, phase-c, dataset-review]
supersedes_section: "Phase C of PLAN.md (Pal/Bach/Giraddi/Sun shortlist)"
status: input-for-replanning
---

# 2026-06-15 - Phase C dataset review (pre-replanning summary)

Working notes from a review session on the Phase C manual dataset downloads
for `mammary_geneset_library`. Purpose: capture the rationale, the candidate
datasets found, and the open decisions, so the final Phase C re-planning can
start from a clean slate. This note is INPUT to that re-plan, not the plan
itself.

## What Phase C is for (reminder)

Script 06 needs differential-expression gene lists to emit stage-specific
developmental GMT sets: `MG_EARLY_POSTNATAL_UP`, `MG_PUBERTY_UP`,
`MG_PUBERTY_TO_ADULT_UP`, `MG_TERMINAL_END_BUD_UP`, `MG_MATURE_DUCTAL_UP`
(plus `_DOWN` counterparts), filtered to top 100-300 genes by logFC at
FDR < 0.05, and consensus sets (e.g. `MG_LUMINAL_PROGENITOR_CONSENSUS`,
intersection across sources, size >= 15). The whole point of using multiple
independent atlases is to make `_CONSENSUS` filtering meaningful - three-plus
sources voting on the same cell-type/stage signature.

## Original shortlist and why it was chosen

The v1/v2 plan named four mouse scRNA-seq atlases, each covering a different
slice of the developmental axis:

- **Pal 2017** (Nat Commun 8:1627) - temporal: pre-puberty -> puberty -> adult
- **Bach 2017** (Nat Commun 8:2128) - reproductive cycle: nulliparous ->
  gestation -> lactation -> post-involution; 15-cluster luminal/basal structure
- **Giraddi 2018** (Cell Reports) - early developmental: embryonic ->
  postnatal -> adult; fMaSC identity
- **Sun 2018** (J Biol Chem 293:8315) - adult virgin (3mo) + pregnant (P12)

Pommier (already loaded via script 04) supplies the human MaSC/LP/mL axis as
an ortholog-mapped fifth anchor.

## Problems found with the original shortlist

1. **Sun 2018 is mislabelled.** PROJECT_CONTEXT.md calls it "pubertal mammary
   gland transcriptomics", but the Heng Sun / Chu-Xia Deng paper actually
   profiles 3-month-old virgin and P12 pregnant mice - both post-pubertal.
   It does NOT cover puberty. So the original four leaned on Pal 2017 and
   Giraddi 2018 alone for the puberty signal.

2. **The four are old (2017-2018).** Newer atlases exist that cover the same
   biology with more cells, more stages, and - critically - direct
   micro-dissected TEB-vs-duct contrasts.

3. **Phase C was entirely mouse-side.** No human-native references staged for
   the downstream Hannon DCIS-IDC analysis.

## Newer candidate datasets found (mouse)

- **Pal 2021** (Breast Cancer Research, doi:10.1186/s13058-021-01445-4).
  Same Visvader/Smyth group as Pal 2017. 132,599 cells across 9 developmental
  stages (late embryogenesis, early postnatal, prepuberty, adult,
  mid/late-pregnancy, post-involution) PLUS micro-dissected pubertal TEBs and
  subtending ducts. The TEB-vs-duct micro-dissection gives a DIRECT contrast
  for `MG_TERMINAL_END_BUD_UP` that Pal 2017 lacks. Same edgeR/Smyth
  conventions. Effectively a strict superset of Pal 2017 for our purposes.
  => Strong candidate to REPLACE Pal 2017 as the transcriptome backbone.

- **Gutierrez 2021** (J Mammary Gland Biol Neoplasia 26(1),
  doi:10.1007/s10911-021-09488-1). Pre-built integration of exactly
  Pal/Bach/Giraddi (the three originals), 53,686 cells, post-natal focus
  (puberty -> post-involution). Either a drop-in replacement for the manual
  consensus step or a cross-check. Worth inspecting BEFORE downloading the
  three originals separately - may make some of them redundant.

- **Saeki 2021** (Communications Biology, doi:10.1038/s42003-021-02201-2).
  Integrated ~50K mouse + ~24K human mammary epithelial atlas; explicitly maps
  mouse cell types to human breast cancer subtypes via TCGA. This is the
  mouse-human BRIDGE the broader project needs (links myc_mouse to the Hannon
  human data). => Add as the cross-species bridge.

- **Gray 2023** (Cell Reports 42:113293, doi:10.1016/j.celrep.2023.113293).
  Brugge lab. The "most relevant" candidate flagged in review. CAVEAT: this is
  PRIMARILY a mass cytometry (CyTOF, ~40-marker panel) + cyclic-IF study, with
  scRNA-seq as a complement on smaller cell numbers - NOT transcriptome-wide
  scRNA by design. Covers puberty, estrous cycle, sex directly and spatially.
  Headline finding: proliferative expansion systematically suppresses lineage
  programs (protein + RNA), with ribosomal-protein enrichment during
  progesterone-driven expansion.
  => ROLE: strong puberty-stage CONSENSUS VOTER and a MANUSCRIPT DISCUSSION
     reference (the proliferation-vs-lineage-suppression concept mirrors the
     MYC story; ribosomal/biosynthetic signature rhymes with the "metabolic
     not apoptotic" reframing). Probably NOT a GMT backbone - the scRNA DE
     tables are likely thin. MUST VERIFY supplement before promoting (see below).

- **Henry 2021** (J Mammary Gland Biol Neoplasia 26(1):43-66, Dos Santos lab).
  Gene-expression signatures for cellular heterogeneity in developing gland.
  Smaller/targeted; lower priority than Pal 2021.

## Human-side candidates (for a separate Phase C-human sub-section)

For the downstream Hannon DCIS-IDC analysis - not the immediate download round,
but worth marking now so it doesn't become a later scramble:

- **Reed 2024** (Nat Genet 56:652-662, doi:10.1038/s41588-024-01688-9) - HBCA,
  ~800K cells / 55 donors; iHBCA integration has CellTypist models for label
  transfer.
- **Kumar 2023** (Nature 620:181-191) - 535,941 cells / 62 women + 120,024
  nuclei / 20 women; 11 cell types, 53 cell states; spatially resolved.
- **Pal 2021 EMBO J** (e107333) - normal, preneoplastic AND tumorigenic human
  states; the preneoplastic component is the closest direct parallel to the
  Hannon normal-DCIS-IDC continuum.

## Proposed revision (to be finalised in the re-plan, NOT yet decided)

1. Replace Pal 2017 -> **Pal 2021 (BCR)** as transcriptome backbone.
2. Keep **Bach 2017** and **Giraddi 2018** for distinct biology (pregnancy
   cycle; embryonic/fMaSC) - but check Gutierrez 2021 first; it may subsume them.
3. Drop **Sun 2018** (mislabelled; now dominated by Pal 2021 + Bach 2017).
4. Add **Saeki 2021** as the mouse-human bridge.
5. Add **Gray 2023** as a puberty consensus voter + manuscript reference
   (pending supplement check).
6. Open a **Phase C-human** sub-section: Reed 2024, Kumar 2023, Pal 2021 EMBO.

## Open questions to resolve in the re-plan

- **Gray 2023 supplement check (blocking its role):** does the scRNA-seq
  component ship transcriptome-wide DE tables (=> could be backbone) or only
  the CyTOF marker panel + cluster annotations (=> voter only)? Check STAR
  Methods / Supplemental info; Cell Reports uses `mmc2.xlsx`, `mmc3.xlsx`, etc.
  Also note GEO accession in case re-derivation from raw counts is ever needed
  (bigger lift, likely out of Phase C scope).
- **Gutierrez 2021 vs the three originals:** does the integrated atlas provide
  adequate per-stage DE to replace downloading Pal/Bach/Giraddi individually?
- **Species handling in script 06:** do human atlases feed the same
  `mammary_geneset_library` GMTs (script 06 handles both species) or a sibling
  human-only collection (parallel script)? Decide before writing 06.
- **Supplement layouts unverified** for Pal 2021, Gutierrez 2021, Saeki 2021,
  Gray 2023 - each needs the same open-and-check pass (which sheet = which
  contrast; column schema) before committing to swap in.

## Knock-on edits if the swap is adopted

- PLAN.md pre-flight checklist (lines ~122-139) and file list in section 1
  (lines ~80-82) both reference `pal2017/bach2017/giraddi2018_supplementary.xlsx`
  and omit `sun2018` - update filenames to match the revised set.
- `data/raw/README.md` - one-line provenance entry per new file (source URL,
  download date, purpose, sheet-of-interest), mirroring
  `data/raw/from_myc_mouse/README.md`.
- Capture the decision trail in a new dated note that supersedes the Phase C
  section of PLAN.md (same pattern as v2 superseding v1).

## Adjacent literature note (not Phase C)

- **Mycn 2025** (Cells, GSE251933, Westlake group) - "Mycn Is Essential for
  Pubertal Mammary Gland Development...". Mycn != c-Myc, but directly adjacent
  to a MYC-driven mammary tumorigenesis manuscript. Candidate for the
  MK_myc_paper reference list regardless of Phase C.

## Citations / DOIs

- Pal 2017: 10.1038/s41467-017-01560-x
- Bach 2017: 10.1038/s41467-017-02001-5
- Giraddi 2018: 10.1016/j.celrep.2018.08.069
- Sun 2018: 10.1074/jbc.RA118.002297
- Pal 2021 (BCR): 10.1186/s13058-021-01445-4
- Gutierrez 2021: 10.1007/s10911-021-09488-1
- Saeki 2021: 10.1038/s42003-021-02201-2
- Gray 2023: 10.1016/j.celrep.2023.113293
- Henry 2021: 10.1007/s10911-021-09486-3
- Reed 2024: 10.1038/s41588-024-01688-9
- Kumar 2023: (Nature 620:181-191, 2023)
- Pal 2021 (EMBO J): EMBO J 40:e107333
