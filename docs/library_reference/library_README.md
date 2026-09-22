# mammary_geneset_library

A versioned, modular, **mouse-native** `.gmt` gene-set collection (plus nine
category-stratified GMTs and parallel human-symbol versions) for fGSEA / GSVA
across the `myc_mouse` and human mitoPPS projects. It consolidates and extends
the gene-set assets built in `myc_mouse`, adding mammary developmental
signatures, TF target sets, MYC vs PGC1a biogenesis-discrimination
intersections, biogenesis-apoptosis intersections, and curated metabolism /
proliferation / apoptosis panels.

See `PLAN.md` for the full action list and `CLAUDE.md` for project rules.

## Output (v1)

Produced by `R/14_export_gmts.R` into `outputs/gmt/`:

- `mammary_mito_myc_metab_v1_mouse.gmt` — master (986 sets, MGI symbols)
- `by_category/0N_<category>_mouse.gmt` — nine category GMTs (sum to the master)
- `human/` — parallel human-symbol GMTs (mouse to human via biomaRt orthologs)
- `outputs/qc/provenance_table.csv` — one row per set (source, tier, category,
  method, fGSEA/GSVA eligibility, mouse/human sizes, unmapped counts, description)
- `CHANGELOG.md` — versioned release notes

### Nine categories (each set has exactly one `Category_Primary`)

| # | Category | Sets | Prefixes |
|---|----------|------|----------|
| 01 | MitoCarta | 127 | `MITOCARTA_*` |
| 02 | MYC_signatures | 17 | `MYC_*`, `MYC_FELSHER_*` |
| 03 | Mammary_development | 203 | `MG_*`, `MG_SIG_*`, `POMMIER_*`, `METABRIC_*` |
| 04 | Metabolism | 103 | `METAB_*`, `GS_METAB_*` |
| 05 | Proliferation | 14 | `PROLIF_*` |
| 06 | TF_targets | 444 | `TFT_*` |
| 07 | Biogenesis_discrimination | 16 | MYC/ER/PGC1a-axis mito set algebra |
| 08 | Apoptosis | 36 | `APOP_*`, `TANG_*`, `CDC_*` |
| 09 | Biogenesis_apoptosis_intersections | 26 | biogenesis x cell-death |

## Pipeline

Numbered R scripts in `R/`, run in order; each writes `.rds` intermediates to
`results/`. Convenience: `source("source_all.R")` runs 00 to 14.

| Script | Builds |
|--------|--------|
| `00_setup_packages.R` | package loader (sourced by every script) |
| `01_pull_msigdb.R` | mouse MSigDB (Hallmark, C2:CP, GO:BP, C3:TFT) |
| `02_load_mitocarta.R` | MitoCarta3.0 (full 1140 inventory; nu/mt + apoptosis splits) |
| `03_ortholog_utilities.R` | bidirectional human/mouse ortholog mapping |
| `04_load_user_curated_sets.R` | GS_metabolic, METABRIC, Pommier |
| `05_load_myc_mouse_assets.R` | MYC signatures, Tang cell-death, CDC |
| `06_build_mammary_dev_sets.R` | atlas cell-state + signalling sets |
| `07_build_tf_target_sets.R` | TF regulons (DoRothEA / ChIP-Atlas / MSigDB) |
| `08_build_intersection_sets.R` | Category 7 biogenesis discrimination |
| `09_build_metabolism_sets.R` | curated MSigDB metabolism panel |
| `10_build_apoptosis_sets.R` | Categories 8 and 9 |
| `11_build_proliferation_controls.R` | Category 5 cell-cycle controls |
| `12_qc_filter_and_dedup.R` | assemble + QC (size filter, Jaccard) -> `library_*` |
| `13_finalize_library.R` | finalize (no external merge) -> `final_library_*` |
| `14_export_gmts.R` | write the GMTs + provenance + CHANGELOG |

## Conventions

- **Identifiers:** mouse = MGI symbols; human export = HGNC via one-to-one
  ortholog mapping (unmapped recorded per set).
- **Size policy:** hard-drop sets < 5 genes; `method` tag is size-aware
  (`both` 15-500, `gsva` 5-14); GSVA is upper-bound-exempt, fGSEA window 15-500
  captured as `fgsea_eligible`.
- **Provenance:** every set carries source, tier, category, method, and a
  verbatim collection-name description.
- **Versioning:** never overwrite; bump the version and update `CHANGELOG.md`.

## Development

Developed with Claude Code in the Positron terminal; scripts are sourced
interactively in Positron (not auto-run end-to-end). Reference scripts from
`myc_mouse` live in `docs/reference/myc_mouse/` (read-only context).
