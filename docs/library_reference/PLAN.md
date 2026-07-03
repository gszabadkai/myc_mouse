# Custom Gene Set Library for fGSEA — Project Plan (v3)

**Project name:** `mammary_geneset_library`
**Output:** versioned master `.gmt` (`mammary_mito_myc_metab_v1_mouse.gmt`) plus
nine **category-stratified `.gmt` files** so fGSEA can be run on any single
functional group in isolation. Plus a parallel human-symbol GMT family.
**Identifier convention:** mouse MGI symbols primary; HGNC symbols for export.

This is v3 of the plan. v2 amended the from-scratch build into a
**consolidate-and-extend** effort after the gene_set_history note. v3 folds in the
decisions of the 2026-06-15 Phase C re-plan note
(`2026-06-15_phase_c_replan_and_library_decisions.md`): the Phase C dataset roster is
finalised, bespoke per-file atlas readers replace the old DE-table scheme, GSVA is
added alongside fGSEA as a second scoring layer, Categories 6 and 7 are rescoped
(two TF lanes; "oncogenic vs developmental biogenesis"), and per-set method tagging
is introduced. Scripts 00-05 are unchanged from v2.

---

## 0. Scope & non-goals

**In scope**
- Build a modular GMT covering nine functional categories (see §3 below).
- Reuse curated assets from `myc_mouse`: ortholog table, Felsher MYC sets,
  Tang et al. cell death sets, consolidated cell death gene catalog.
- Preserve the MitoCarta modifications already in use in
  `08_mitoPPS_analysis.R`: nu/mt OXPHOS split (overall + per complex),
  Pro/Anti apoptotic split.
- Extend with: MSigDB Hallmark+C2+C5, mammary developmental signatures,
  TF target sets, MYC and PGC1 intersections, METABRIC and Pommier
  human-derived sets, GS_metabolic compendium.
- Output: one master GMT + nine category GMTs, plus human-ortholog twins.
- **Two complementary scoring layers, both supported.** fGSEA stays the
  contrast-level competitive layer (signed NES on a ranked DESeq2 Wald `stat`;
  size 15-500). **GSVA/ssGSEA is added** as the per-sample scorer for the large
  developmental/lineage signatures that exceed the fGSEA size ceiling; its scores
  feed the same genotype x timepoint interaction model as the gene-level analysis.
  Every set therefore carries a **method tag** (`fgsea` / `gsva` / `both`) in the
  provenance table, and the QC size filter gains a **per-set exemption flag** so
  `gsva`-tagged large sets are not dropped (see §4 step 16).
- **Framing.** Category 7 is the discrimination of **oncogenic vs developmental**
  mitochondrial biogenesis (MYC-driven vs ER/PGC1a-driven), at the
  hypothesis-generation level.
- QC: deduplication, size filtering (15 <= |set| <= 500, with method exemptions),
  Jaccard overlap, provenance/tier annotation.

**Out of scope (for v1)**
- scRNA-seq de novo signature derivation.
- Modifying anything in `myc_mouse` itself (Phase E touches it separately).
- **GRN inference / confirmation (SCENIC, pySCENIC, decoupleR).** The TF
  work here is hypothesis-level (generic regulons + Gray's mammary-filtered CHEA);
  GRN confirmation is explicitly deferred and is not part of this paper.

---

## 1. Project layout

```
mammary_geneset_library/
├── CLAUDE.md
├── PLAN.md
├── PROMPT_DAY1.md
├── gene_set_history.md                    # user's note on prior work
├── README.md
├── R/
│   ├── 00_setup_packages.R
│   ├── 01_pull_msigdb.R                   # MSigDB Hallmark/C2/C5/C8 (mouse)
│   ├── 02_load_mitocarta.R                # MitoCarta3.0 + nu/mt + Pro/Anti splits
│   ├── 03_ortholog_utilities.R            # inner_join on cached biomaRt table
│   ├── 04_load_user_curated_sets.R        # GS_metabolic + METABRIC + Pommier
│   ├── 05_load_myc_mouse_assets.R         # Felsher MYC + Tang cell-death +
│   │                                       #   cell_death_genes_consolidated.rds
│   ├── 06_build_mammary_dev_sets.R        # atlas readers -> MG_* cell-state sets
│   ├── 07_build_tf_target_sets.R          # two TF lanes (generic + Gray-CHEA) incl ESR1
│   ├── 08_build_intersection_sets.R       # oncogenic-vs-developmental mito algebra
│   ├── 09_build_metabolism_sets.R         # MSigDB + GS_metabolic
│   ├── 10_build_apoptosis_sets.R          # MitoCarta Pro/Anti + Tang Apoptosis +
│   │                                       #   consolidated; intersections w/ biogen.
│   ├── 11_build_proliferation_controls.R  # Hallmark E2F/G2M/etc.
│   ├── 12_qc_filter_and_dedup.R
│   ├── 13_merge_with_existing.R           # join with gene_sets_final_annotated.csv
│   └── 14_export_gmts.R                   # master + 9 category GMTs (mouse+human)
├── data/
│   ├── raw/
│   │   ├── user_curated/                  # the three xlsx
│   │   │   ├── GS_metabolic_genes_list.xlsx
│   │   │   ├── METABRIC_…biclustergenelists….xlsx
│   │   │   └── Pommier_dev_genesets.xlsx
│   │   ├── from_myc_mouse/                # snapshot — see README.md inside
│   │   │   ├── README.md
│   │   │   ├── ortholog_table.rds
│   │   │   ├── gene_sets_list.rds
│   │   │   ├── cell_death_genes_consolidated.rds
│   │   │   ├── myc_signature_genesets.gmx
│   │   │   ├── felsher_integrative_signature.csv
│   │   │   ├── mitocarta_pathways.csv     # legacy custom subset (reference only)
│   │   │   └── cell-death/                # Tang et al 15 modalities
│   │   ├── mammary_dev_atlases/           # Phase C atlas supplements (manual download 2026-06-15)
│   │   │   ├── README.md                  # resources table + prepended provenance header
│   │   │   ├── Gray-et-al_2023_Suppl1.xlsx     # BA/HS/AP major-type (gsva)
│   │   │   ├── Gray-et-al_2023_Suppl2.xlsx     # HE vs LE subtype DE (fgsea)
│   │   │   ├── Gray-et-al_2023_Suppl3.xlsx     # CHEA subtype TF lists (gsva)
│   │   │   ├── Gray-et-al_2023_Suppl4.xlsx     # TEB vs ductal DE (fgsea)
│   │   │   ├── Gray-et-al_2023_Suppl5.xlsx     # CHEA TEB/ductal TF lists (gsva)
│   │   │   ├── Garcia-Sola-et-al_2021_Suppl1.xlsx  # 20 global clusters
│   │   │   ├── Garcia-Sola-et-al_2021_Suppl2.xlsx  # HS trajectory clusters
│   │   │   ├── Garcia-Sola-et-al_2021_Suppl3.xlsx  # Alv trajectory clusters
│   │   │   ├── Scheele-et-al_2017_Suppl1.xlsx      # 9 pubertal SC clusters
│   │   │   ├── Saeki-et-al_2021_Suppl4.xlsx        # C1-C6 cluster markers (fgsea)
│   │   │   ├── Saeki-et-al_2021_Suppl5.xlsx        # pseudotime MSigDB pathways (ref)
│   │   │   ├── Saeki-et-al_2021_Suppl6.xlsx        # scGSVA Stem/Basal/Alv/Hor (gsva)
│   │   │   ├── Pal-et-al_2017_Suppl2.xlsx          # timeline clusters I-VI
│   │   │   ├── Pal-et-al_2017_Suppl3.xlsx          # basal/luminal signatures (fgsea)
│   │   │   ├── Pal-et-al_2017_Suppl4.xlsx          # adult vs pubertal clusters
│   │   │   ├── Henry-et-al_2021_Suppl1.xlsx        # consensus lineage signatures
│   │   │   ├── Henry-et-al_2021_Suppl2.xlsx        # mEC cluster signatures
│   │   │   ├── Henry-et-al_2021_Suppl5.xlsx        # mNEC cluster signatures
│   │   │   ├── Giraddi-et-al_2018_Suppl3.xlsx      # NMF fMaSC/balancer sigs (gsva)
│   │   │   ├── Chung-et-al_2019_Suppl2.xlsx        # ATAC TF signatures
│   │   │   └── Chung-et-al_2019_Suppl3.xlsx        # ATAC accessible gene sigs
│   │   └── mitocarta3_mouse.xlsx          # full MitoCarta3.0 (Broad)
│   └── processed/
├── docs/
│   └── reference/myc_mouse/               # reference R scripts (read-only)
│       ├── README.md
│       ├── 01_load_data.R
│       ├── 08_mitoPPS_analysis.R
│       ├── 09_mitoPPS_vs_fgsea_comparison.R
│       ├── 09_cell_death_pathway_summary.R
│       └── Myc_timecourse_analysis_GS_sandbox.R
├── results/                               # .rds intermediates (gitignored)
├── outputs/
│   ├── gmt/
│   │   ├── mammary_mito_myc_metab_v1_mouse.gmt           # master
│   │   ├── by_category/
│   │   │   ├── 01_mitocarta_mouse.gmt
│   │   │   ├── 02_myc_signatures_mouse.gmt
│   │   │   ├── 03_mammary_development_mouse.gmt
│   │   │   ├── 04_metabolism_mouse.gmt
│   │   │   ├── 05_proliferation_mouse.gmt
│   │   │   ├── 06_tf_targets_mouse.gmt
│   │   │   ├── 07_biogenesis_discrimination_mouse.gmt
│   │   │   ├── 08_apoptosis_mouse.gmt
│   │   │   └── 09_biogenesis_apoptosis_intersections_mouse.gmt
│   │   └── human/                                         # parallel HGNC versions
│   ├── qc/
│   └── README.md
└── .gitignore
```

Removed from `data/raw/` in v3: `pal2017_supplementary.xlsx`,
`bach2017_supplementary.xlsx`, `giraddi2018_supplementary.xlsx` (superseded by the
`mammary_dev_atlases/` roster). Sun 2018 dropped entirely (mislabelled as pubertal).

---

## 2. Pre-flight checkpoints (before starting Claude Code)

Each item must be checked off before pasting `PROMPT_DAY1.md`. The
checkpoint list is mirrored at the end of `PROMPT_DAY1.md` for tomorrow
morning.

- [ ] Project skeleton created (`mkdir` commands above).
- [ ] Git initialised, on `dev` branch, NOT `main`.
- [ ] **`data/raw/user_curated/` contains all three xlsx files.**
- [ ] **`data/raw/from_myc_mouse/` contains all seven items** (ortholog_table.rds,
      gene_sets_list.rds, cell_death_genes_consolidated.rds,
      myc_signature_genesets.gmx, felsher_integrative_signature.csv,
      mitocarta_pathways.csv, cell-death/*.csv).
- [ ] **`data/raw/from_myc_mouse/README.md` filled in** with copy date and
      source git commit hash from myc_mouse.
- [ ] `docs/reference/myc_mouse/` contains the five reference scripts with
      its own `README.md` flagging them as read-only.
- [ ] `data/raw/mitocarta3_mouse.xlsx` present (the full MitoCarta3.0 Sheet 4
      with MitoPathways — separate from the legacy custom subset).
- [ ] **`data/raw/mammary_dev_atlases/` contains the Phase C roster**
      (download_date 2026-06-15): Gray Suppl1-5; Garcia-Sola Suppl1-3;
      Scheele Suppl1; Saeki Suppl4/5/6; Pal-2017 Suppl2/3/4; Henry Suppl1/2/5;
      Giraddi Suppl3; Chung Suppl2/3 — plus its `README.md` with the prepended
      provenance header.
- [ ] `PLAN.md`, `CLAUDE.md`, `PROMPT_DAY1.md`, `gene_set_history.md` all in
      project root.
- [ ] Initial `.gitignore` written and `git add` + first commit done so
      Claude Code starts from a clean tree.
- [ ] Positron open in project folder; one R session running, no packages
      sourced yet.
- [ ] Claude Code launched with `claude --model opus`.

---

## 3. The nine functional categories (output structure)

This is the canonical category list. Every gene set in the final library
maps to exactly one `Category_Primary`. fGSEA can then be run by loading
only that category's GMT.

| # | Category_Primary                            | Contents (examples)                                                                                                                                              |
|---|---------------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| 1 | `MitoCarta`                                 | `MITOCARTA_OXPHOS`, `MITOCARTA_OXPHOS_NU`, `MITOCARTA_OXPHOS_MT`, `MITOCARTA_COMPLEX_1_NU`, `MITOCARTA_COMPLEX_1_MT`, …, TCA, FAO, ribosome, dynamics, mitophagy. |
| 2 | `MYC_signatures`                            | `HALLMARK_MYC_TARGETS_V1/V2`, Dang/Schuhmacher/Coller/Kim/Zeller, Felsher GMX-derived sets, **`MYC_FELSHER_INTEGRATIVE_SIGNATURE`**.                              |
| 3 | `Mammary_development`                       | GO mammary terms; cell-state/lineage `MG_*` sets from the Phase C atlas roster (Gray, Garcia-Sola, Scheele, Saeki, Pal-2017, Henry, Giraddi, Chung); Pommier; METABRIC; hormone/EGFR/IGF signalling. |
| 4 | `Metabolism`                                | MSigDB Hallmark + Reactome + KEGG metabolism, GS_metabolic compendium, granular AA/lipid/nucleotide/one-carbon.                                                  |
| 5 | `Proliferation`                             | Hallmark `E2F_TARGETS`, `G2M_CHECKPOINT`, `MITOTIC_SPINDLE`, `MTORC1_SIGNALING` (confounder controls).                                                            |
| 6 | `TF_targets`                                | **Two lanes, side by side, not merged.** Generic regulon lane (DoRothEA A/B + ChIP-Atlas) for ESR1, NRF1, GABPA, ESRRA, PPARGC1A-axis, MYC, E2F -> used for Category-7 set algebra (`fgsea`): `TFT_<TF>_DOROTHEA_AB`, `TFT_<TF>_CHIPATLAS`. Mammary-context lane (Gray Suppl3/Suppl5 CHEA) for per-sample TF-program scoring (`gsva`): `TFT_<TF>_GRAY_<CONTEXT>`, CONTEXT in {HS_HE, HS_LE, AP_HE, AP_LE, BA_HE, BA_LE, TEB, DUCTAL}. |
| 7 | `Biogenesis_discrimination` (intersections) | **Rescoped: oncogenic vs developmental biogenesis.** Pass-1 mito algebra (Decision 4): `DEVELOPMENTAL_MITO` (ESR1 ∩ MitoCarta), `CORE_MITO` (PGC1a axis = union(ESRRA,NRF1,GABPA) ∩ MitoCarta), `MYC_MITO`, and complements `MYC_SPECIFIC_MITO`, `ESR1_NOT_CORE_MITO` (discovery), `ESR1_AND_CORE_MITO`, `MYC_AND_DEVELOPMENTAL_MITO`, `MYC_AND_CORE_MITO` (hijack), `MYC_DEV_CORE_TRIPLE`. |
| 8 | `Apoptosis`                                 | Tang Apoptosis, MitoCarta `Apoptosis-PRO`/`Apoptosis-ANTI`, cell_death_consolidated Pro/Anti/Mito/Core, Hallmark APOPTOSIS.                                       |
| 9 | `Biogenesis_apoptosis_intersections`        | `MYC_MITO_∩_APOPTOSIS_PRO`, `MYC_MITO_∩_APOPTOSIS_ANTI`, `PPARGC1A_MITO_∩_APOPTOSIS_PRO/ANTI`, MITO ∩ Tang Apoptosis subsets.                                     |

Category count stays at **9** (Category 7 is rescoped, not added to). A `Category`
(finer) field gives sub-grouping within each Primary (e.g. for `MitoCarta`:
`OXPHOS`, `OXPHOS_mt`, `OXPHOS_nu`, `Complex_I_nu`, `TCA`, `FAO`, `Apoptosis_PRO`,
etc.). Sets that legitimately belong to multiple sub-categories use pipe-delimited
values (e.g. `Category = "Apoptosis_PRO|Complex_I_nu"`).

---

## 4. Action points

### Phase A — Setup (Claude Code; trivial, completed in one short session)

1. Verify `data/raw/from_myc_mouse/` snapshot is present and READMEs are
   filled in. Refuse to proceed if missing — these are non-negotiable inputs.
2. Write `00_setup_packages.R` declaring all packages once:
   `msigdbr`, `dplyr`, `readr`, `tidyr`, `readxl`, `stringr`, `fgsea`,
   **`GSVA`**, `gprofiler2`, `org.Mm.eg.db`, `org.Hs.eg.db`, `biomaRt` (only if
   rebuilding ortholog table is needed), `here`, `purrr`; optionally `babelgene` /
   `orthogene` for sandbox checks.

### Phase B — Build sets from automatable + snapshot sources (Claude Code)

3. **`01_pull_msigdb.R`** — pull mouse MSigDB:
   - Hallmark (all).
   - C2:CP:REACTOME, C2:CP:KEGG_MEDICUS, C2:CP:WIKIPATHWAYS.
   - C2:CGP mammary/breast (Charafe, Landis, Stein ERBB2…) and MYC
     (Dang, Schuhmacher, Coller, Kim, Zeller).
   - C3:TFT for MYC/NRF1/GABPA/ESRRA.
   - C5:GOBP mammary/development/branching/hormone response.
   - Save `results/msigdb_sets.rds`.

4. **`02_load_mitocarta.R`** — load full MitoCarta3.0 mouse, build sets
   with the EXACT modifications used in `08_mitoPPS_analysis.R`:
   - `MITOCARTA_ALL`, `MITOCARTA_NUCLEAR_ENCODED`, `MITOCARTA_MTDNA_ENCODED`.
   - Pathways from MitoPathways column: TCA, pyruvate, FAO, AA metab,
     mito ribosome, mito translation, mtDNA maintenance, import, dynamics,
     mitophagy, ROS detox, Fe-S, heme, CoQ.
   - **OXPHOS splits (mandatory):**
       - `MITOCARTA_OXPHOS_NU` — all OXPHOS subunits MINUS mt-* genes
       - `MITOCARTA_OXPHOS_MT` — only the 13 mt-* protein-coding genes
       - `MITOCARTA_COMPLEX_I_NU`, `MITOCARTA_COMPLEX_I_MT`,
         `MITOCARTA_COMPLEX_III_NU`, `MITOCARTA_COMPLEX_III_MT`,
         `MITOCARTA_COMPLEX_IV_NU`, `MITOCARTA_COMPLEX_IV_MT`,
         `MITOCARTA_COMPLEX_V_NU`, `MITOCARTA_COMPLEX_V_MT`
         (Complex II has no mt-* genes — only nu version)
   - **Apoptosis splits (mandatory):**
       - `MITOCARTA_APOPTOSIS_PRO`, `MITOCARTA_APOPTOSIS_ANTI`
   - Read the existing convention from `docs/reference/myc_mouse/08_mitoPPS_analysis.R`
     to confirm exact membership. Do not regenerate from scratch if it
     conflicts with that script — the existing convention is authoritative.

5. **`03_ortholog_utilities.R`** — central bidirectional mapping utility.
   - Load `data/raw/from_myc_mouse/ortholog_table.rds` as the cached table.
     Do NOT rebuild via biomaRt unless the cached file is missing.
   - Function `human_to_mouse(genes, expand_one_to_many = TRUE)`: returns a
     tibble `{human_symbol, mouse_symbol, orthology_type}`. **Default expands
     one-to-many.** Use `dplyr::inner_join` on the cached table, NOT
     `deframe()` + named-vector lookup.
   - Function `mouse_to_human(genes, collapse_many_to_one = "union")`.
   - Append per-gene results to `outputs/qc/ortholog_mapping_report.csv`.
   - **Unit test inside the script** (in an `if (FALSE)` block I run
     manually): assert that a known one-to-many human gene returns ≥ 2 mouse
     symbols. Suggest a test gene from the cached table during planning.
   - Sandbox-only: `babelgene::orthologs()` and/or `orthogene` for spot
     checks; not in production pipeline.

6. **`04_load_user_curated_sets.R`** — three sub-routines:
   - **GS_metabolic** — long-format human × Classification; join via the
     file's own built-in mouse orthologs sheet. Tag provenance
     `source = "GS_metabolic (pre-mapped)"`. Output `GS_METAB_*` sets,
     Category_Primary = `Metabolism`.
   - **METABRIC biclusters** — process 18 sheets matching
     `^MB[1-3]\.(bicluster|hi\.cv|metab)\.group[12]$`. Map human → mouse via
     `03_ortholog_utilities.R::human_to_mouse`. Output `METABRIC_*` sets,
     Category_Primary = `Mammary_development` (Category = `Human BC subtype`).
   - **Pommier** — six MaSC/LP/mL up/down sets, map human → mouse. Output
     `POMMIER_*` sets, Category_Primary = `Mammary_development`.

7. **`05_load_myc_mouse_assets.R`** — the critical reuse layer.
   - Load `data/raw/from_myc_mouse/gene_sets_list.rds` and extract ONLY the
     MYC subset (names starting `MYC_`). Drop MitoCarta and Hallmark portions
     — those are regenerated in steps 4 and 3 respectively.
   - Verify `MYC_felsher_integrative_signature` is present and non-empty.
   - Load the Felsher MYC GMX directly from
     `data/raw/from_myc_mouse/myc_signature_genesets.gmx` and ortholog-map
     human → mouse (via `03_ortholog_utilities.R`). Cross-check vs
     gene_sets_list.rds; flag any drift.
   - Load `data/raw/from_myc_mouse/cell-death/*.csv` (Tang et al. 15 cell
     death modalities; `gene` column already mouse). One set per file:
     `TANG_<DEATHTYPE>`. Category_Primary = `Apoptosis` for all 15 (Category =
     death type name).
   - Load `data/raw/from_myc_mouse/cell_death_genes_consolidated.rds`
     (own MyGene+HomoloGene orthology — **do NOT remap**). Build `CDC_*` sets
     stratified by `pathway × pro_anti × mitochondrial × core`.
     Category_Primary = `Apoptosis`.

8. **`06_build_mammary_dev_sets.R`** — atlas-roster cell-state /
   lineage signatures. **The old `MG_*_STAGE_UP/DOWN` DE-contrast scheme is
   REPLACED** (the supplements give cell-state/lineage signatures, not clean
   per-stage DE contrasts; the only genuine contrast is Gray TEB-vs-ductal).

   - **Workflow (Decision 2), one file at a time with a stop after each reader:**
     Claude Code writes a per-file inspection block (sheet names, dims, head) ->
     user runs it in Positron -> Claude Code writes **one bespoke reader for that
     file** -> **STOP for confirmation**: the user runs the reader in Positron and
     reviews the result before Claude Code moves on to the next file. Structures
     genuinely differ, so readers are never batched. All readers emit one tidy
     schema: `{set_name, gene, logFC, padj, source, stage, method}`.
   - **Naming grammar (Decision 6):**
     `<PREFIX>_<CELLSTATE-or-CONTRAST>_<SOURCE>[_<UP|DN>][_<SUBTOKEN>]`.
     PREFIX `MG_`; SOURCE token in {GRAY, GARCIASOLA, SCHEELE, SAEKI, PAL2017,
     HENRY, GIRADDI, CHUNG}; SUBTOKEN only for disambiguation
     (`_GSVA`/`_C` for same state via two analyses in one source; `_ATAC` for
     Chung). Stage lives in the **description**, not the name (revisable);
     only genuine contrasts carry stage in the name (Gray `TEB_VS_DUCTAL`).
     Description field carries the verbatim full table collection name + Suppl#
     + stage annotation + method.
   - **Method assignment (Decision 3 table) - provisional until confirmed per set.**
     The assignment below is the starting guide; the actual method for each set is
     decided only after that set's file/sheet is inspected at the reader-confirmation
     stop above (size, presence of source DE stats, signed vs unsigned):
     - `fgsea` (in-range, signed/truncatable): Gray Suppl2 (HE/LE),
       Gray Suppl4 (`MG_TEB_VS_DUCTAL_GRAY_UP/_DN`), Garcia-Sola clusters
       (truncated), Pal-2017 Suppl3, Scheele clusters, Henry mEC/mNEC.
     - `gsva` (large, used whole): Giraddi fMaSC/balancers, Saeki scGSVA
       Stem/Basal/Alv/Hor, Gray major-type (BA/HS/AP).
     - `both` where a program exists in both a tight and a large form
       (concordance test).
   - **No graded inspection tiers.** Every reader is confirmed with a full
     'real look' (run in Positron + review the emitted sets) before its structure
     and method are finalized - including the files the earlier draft would have
     generated straight. The cluster-naming and stage annotations for the harder
     files (Garcia-Sola ~39 numbered clusters across 3 sheets, Scheele 9 mixed
     clusters) are decided at that same per-file stop.
   - **Consensus sets** vote across the new roster:
     `MG_<STATE>_CONSENSUS` = intersection across sources, size >= 15
     (the source token is what makes voting possible, e.g.
     `MG_BASAL_GRAY` vs `MG_BASAL_PAL2017` vs `MG_BASAL_SAEKI`).
   - Add signalling pathway sets relevant to mammary epithelium
     (EGFR, IGF1/IGF1R, ER/PR/PRL hormone response) pulled from
     msigdb_sets.rds, retagged Category_Primary = `Mammary_development`.

9. **`07_build_tf_target_sets.R`** — two TF lanes, side by side,
   **do NOT merge** (same object type, different method):
   - **Generic regulon lane** (`fgsea`; for Category-7 set algebra): DoRothEA
     A/B + ChIP-Atlas. TF roster **now includes ESR1** (ESR1 != ESRRA; the
     original roster lacked a ligand-activated estrogen receptor). PGC1a is a
     coactivator, not DNA-binding -> "PGC1a axis" = `union(ESRRA, NRF1, GABPA)`.
     Output `TFT_<TF>_DOROTHEA_AB`, `TFT_<TF>_CHIPATLAS`.
   - **Mammary-context lane** (`gsva`; per-sample TF-program scoring): Gray
     Suppl3 (adult HE/LE CHEA) + Suppl5 (TEB/ductal CHEA). Output
     `TFT_<TF>_GRAY_<CONTEXT>`.
   - **Developmental arm pass 1 = ESR1 only** (mouse model is 6-12 wk, active
     estrous cycle, no pregnancy). PGR/STAT5 (+/- GATA3/FOXA1/ELF5) deferred to
     pass 2.
   - **Gray-CHEA mito-overlap TF search (Decision 5).** For each Gray
     CHEA TF, intersect its overlapping signature genes with the **full
     MitoCarta (~1140 genes, not the biogenesis subset)** — captures TFs acting
     on apoptotic/metabolic mito aspects, relevant to the apoptosis-evasion
     thread. Annotate each overlap gene by MitoCarta sub-pathway (from script
     02's splits); **no per-aspect cap** - overlap counts are expected to be
     small, so the full annotated list is kept; revisit only if the use review
     shows it is unwieldy. Artifact:
     `outputs/qc/gray_chea_mito_tf_shortlist.csv` with columns
     `TF | n_mito_total | n_oxphos | n_mitoribo | n_apoptosis_pro |
     n_apoptosis_anti | n_metab_other | frac | hyperg_p`.
     Validation: search should recover ESRRA/NRF1/GABPA in-context; novel
     mito-overlap TFs -> GSVA panel; optionally emit small
     `TFT_<TF>_GRAY_<CONTEXT>_MITO` (gsva-only/descriptive). This is a
     stop-and-check the user reviews on real Suppl3/Suppl5 before the panel locks.

10. **`08_build_intersection_sets.R`** — Category 7, rescoped to
    **oncogenic vs developmental biogenesis** (Decision 4 pass-1 algebra; all
    gated size >= 15). `MITO` = full MitoCarta (~1140); generic `*_targets`
    from the generic regulon lane:
    ```
    DEVELOPMENTAL_MITO         = intersect(ESR1_targets, MITO)
    CORE_MITO                  = intersect(union(ESRRA,NRF1,GABPA)_targets, MITO)  # PGC1a axis
    MYC_MITO                   = intersect(MYC_targets, MITO)
    MYC_SPECIFIC_MITO          = setdiff(MYC_MITO, union(DEVELOPMENTAL_MITO, CORE_MITO))  # key set
    ESR1_NOT_CORE_MITO         = setdiff(DEVELOPMENTAL_MITO, CORE_MITO)                   # discovery set
    ESR1_AND_CORE_MITO         = intersect(DEVELOPMENTAL_MITO, CORE_MITO)                 # canonical route
    MYC_AND_DEVELOPMENTAL_MITO = intersect(MYC_MITO, DEVELOPMENTAL_MITO)
    MYC_AND_CORE_MITO          = intersect(MYC_MITO, CORE_MITO)                           # hijack set
    MYC_DEV_CORE_TRIPLE        = intersect(intersect(MYC_MITO, DEVELOPMENTAL_MITO), CORE_MITO)
    ```
    `ESR1_NOT_CORE_MITO` is the discovery set (ER-reachable mito genes the PGC1a
    axis does not explain); a near-empty result is itself evidence for the
    exclusivity model. Caveat to record: ER's broad cistrome makes
    `DEVELOPMENTAL_MITO` large/noisier — read fGSEA enrichment as
    "ER-associated", not "ER-driven".
    - **METABRIC integration:** `MYC_TARGETS_∩_METABRIC_MB<N>_BICLUSTER_GROUP<G>`
      where size >= 15.

11. **`09_build_metabolism_sets.R`** — granular metabolism from MSigDB + GS_metabolic:
    - `METAB_GLYCOLYSIS`, `METAB_TCA`, `METAB_PPP`, `METAB_OXPHOS_MSIGDB`.
    - AA: serine-glycine-1C, glutamine, BCAA, arginine-proline, methionine-cysteine.
    - Lipid: FA synthesis, β-oxidation, cholesterol, phospholipid.
    - Nucleotide: purine de novo, pyrimidine de novo, salvage, dNTP.
    - GS_METAB_* alternatives kept alongside as a separate curation.

12. **`10_build_apoptosis_sets.R`** — Categories 8 and 9:
    - Collect from step 4 (MitoCarta Pro/Anti), step 7 (Tang Apoptosis +
      cell_death_consolidated), MSigDB Hallmark `APOPTOSIS`.
    - Category 9 intersections:
      `MYC_MITO_∩_APOPTOSIS_PRO`, `MYC_MITO_∩_APOPTOSIS_ANTI`,
      `CORE_MITO_∩_APOPTOSIS_PRO`, `CORE_MITO_∩_APOPTOSIS_ANTI`,
      `MITOCARTA_NU_∩_TANG_APOPTOSIS`, `MITOCARTA_NU_∩_CDC_APOPTOSIS_PRO`.
    - Each only emitted if size >= 15.

13. **`11_build_proliferation_controls.R`** — Hallmark E2F/G2M etc., already
    in msigdb_sets.rds; relabel and re-tag Category_Primary = `Proliferation`.

### Phase C — Manual downloads (you) — COMPLETED 2026-06-15

14. **Atlas supplement roster** staged into
    `data/raw/mammary_dev_atlases/` (download_date 2026-06-15): Gray Suppl1-5,
    Garcia-Sola Suppl1-3, Scheele Suppl1, Saeki Suppl4/5/6, Pal-2017 Suppl2/3/4,
    Henry Suppl1/2/5, Giraddi Suppl3, Chung Suppl2/3. README provenance header
    prepended. Roster rationale (Decision 1): Gray/Garcia-Sola/Scheele backbone;
    Saeki cross-species bridge; Pal-2017/Henry consensus voters; Giraddi/Chung
    embryonic-stem slice. **Pal-2021 (BCR) PARKED** (no usable supplement);
    **Sun-2018 DROPPED** (mislabelled as pubertal); Bach-2017 reachable via
    Garcia-Sola/Saeki integrations. Human datasets (Reed/Kumar/Pal-2021-EMBO/
    Chen) deferred to a Phase C-human sub-section.
15. **MYC ChIP target lists** — *only if* generic-regulon `MYC ∩ MitoCarta` in
    step 10 yields < 15 genes. Sources: Sabò 2014, Walz 2014, ChIP-Atlas.

### Phase D — QC, merge, export (Claude Code)

16. **`12_qc_filter_and_dedup.R`** — size filter (15-500), within-set
    dedup, pairwise Jaccard, flag pairs J > 0.7, write heatmap and size
    distribution figures. **Add a per-set size-exemption flag keyed on the
    method tag** so `gsva`-tagged large sets (Giraddi ~477-937, Saeki scGSVA
    160/240/500/200, Gray major-type) survive the 15-500 filter; `<15` is still
    dropped for either method.

17. **`13_merge_with_existing.R`** — merge with existing
    `gene_sets_final_annotated.csv` from myc_mouse. Resolve name collisions
    (existing names win; append `_v2` if intentional). Apply
    Category_Primary / Category labels per §3.

18. **`14_export_gmts.R`** — write:
    - Master GMT: `outputs/gmt/mammary_mito_myc_metab_v1_mouse.gmt`.
    - Nine category GMTs: `outputs/gmt/by_category/0N_<category>_mouse.gmt`.
    - Mouse → human ortholog mapping → parallel human GMTs in `outputs/gmt/human/`.
    - **`outputs/qc/provenance_table.csv`** (one row per set) — add a
      **`method` column** (`fgsea` / `gsva` / `both`) and ensure the **full table
      collection name is carried verbatim in the description field** (alongside
      set_name, source, tier, category_primary, category, size_mouse, size_human,
      n_unmapped_to_human, n_unmapped_from_human, notes).

### Phase E — Integration with `myc_mouse` (Claude Code, in myc_mouse repo)

19. Update `01_load_data.R` to load the new GMT (master + relevant category
    GMTs). Keep old path as fallback.
20. Re-run scripts 04, 06, 10, 12 with the new GMT; diff report.
21. **New cross-species fGSEA:** mouse 6W/12W contrasts against METABRIC
    bicluster signatures.
22. Quarto report: old vs new GMT diff for paper.

---

## 5. fGSEA + GSVA execution conventions (downstream use)

- **fGSEA:** rank by signed DESeq2 Wald `stat`; `fgseaMultilevel`,
  `minSize = 15`, `maxSize = 500`; `collapsePathways()` for redundancy in
  reports, keep full results. Category-stratified: load one category GMT.
- **GSVA/ssGSEA:** per-sample scores for `gsva`-tagged sets feed the
  genotype x timepoint interaction model. **Input scale is OPPOSITE to mitoPPS**
  — use log-CPM or VST with `kcdf="Gaussian"`, OR raw counts with
  `kcdf="Poisson"`; do NOT carry the mitoPPS linear-scale choice over. Scores
  are **cohort-relative**: score ALL samples together in one run, never in two
  batches. Interpretive ceiling: a bulk signature score conflates composition vs
  per-cell regulation; report accordingly (deconvolution out of scope).

---

## 6. Open decisions to confirm in Day 1 planning

- **D1 (RESOLVED):** Ortholog source — cached biomaRt table; `inner_join`,
  expand one-to-many.
- **D2 (RESOLVED):** Keep all 18 METABRIC sets.
- **D3 (RESOLVED):** GS_metabolic uses its own pre-mapped sheet; tag in provenance.
- **D4 (RESOLVED):** cell_death_genes_consolidated uses MyGene+HomoloGene, do NOT remap.
- **D5 (RESOLVED):** METABRIC sits in `Mammary_development` with
  `Category = "Human BC subtype"` (preserve nine canonical groups).
- **D6 (RESOLVED):** All 15 Tang cell death modalities under `Apoptosis` Primary
  with Category = death type (avoids 10-category drift).
- **D7 (RESOLVED):** Version master GMT, with CHANGELOG.
- **D8 (RESOLVED by re-plan note):** Category count stays 9; Category 7
  rescoped to oncogenic-vs-developmental biogenesis; ESR1 added to TF roster and
  kept inside the mito algebra; GRN confirmation out of scope; pass-2 TFs
  (PGR/STAT5/GATA3/FOXA1/ELF5) deferred.
