# Dropped gene sets (< 5 genes)

Generated 2026-06-19 from the builder outputs in `results/`.
These 40 sets fall below the hard floor of 5 genes and are removed by
`R/12_qc_filter_and_dedup.R` (user decision 2026-06-18: keep >= 5, hard-drop < 5).
They are unusable for fGSEA (needs 15-500) and unstable for GSVA (floor 5).
All come from two legacy builders; no `MG_*`, `METAB_*`, `APOP_*`, `PROLIF_*`,
`TFT_*`, or intersection sets are affected.

Script 12 also writes this list to `outputs/qc/dropped_small_sets.csv`.

## gs_metabolic (10 sets)

Source: GS_metabolic curated compendium (pre-mapped mouse orthologs).

| set_name | n_genes | category |
|---|---|---|
| `GS_METAB_B12` | 4 | B12 |
| `GS_METAB_BH4` | 4 | BH4 |
| `GS_METAB_COMPLEX_II` | 4 | COMPLEX_II |
| `GS_METAB_LYSINE` | 4 | LYSINE |
| `GS_METAB_PIGMENT` | 4 | PIGMENT |
| `GS_METAB_PROPANOATE` | 4 | PROPANOATE |
| `GS_METAB_NUCLEOTIDE_SUGAR` | 3 | NUCLEOTIDE_SUGAR |
| `GS_METAB_VITAMIN_B6` | 3 | VITAMIN_B6 |
| `GS_METAB_SERINE` | 2 | SERINE |
| `GS_METAB_XYULOSE` | 2 | XYULOSE |

## mitocarta (30 sets)

Source: MitoCarta3.0 per-pathway sets (Sheet 4), after mt-* removal.

| set_name | n_genes | category |
|---|---|---|
| `MITOCARTA_AMIDOXIME_REDUCING_COMPLEX` | 4 | MitoCarta pathway |
| `MITOCARTA_CII_ASSEMBLY_FACTORS` | 4 | MitoCarta pathway |
| `MITOCARTA_CII_SUBUNITS` | 4 | MitoCarta pathway |
| `MITOCARTA_CYTOCHROMES` | 4 | MitoCarta pathway |
| `MITOCARTA_GLYCINE_CLEAVAGE_SYSTEM` | 4 | MitoCarta pathway |
| `MITOCARTA_LIPOATE_INSERTION` | 4 | MitoCarta pathway |
| `MITOCARTA_MT_MRNA_MODIFICATIONS` | 4 | MitoCarta pathway |
| `MITOCARTA_PROPANOATE_METABOLISM` | 4 | MitoCarta pathway |
| `MITOCARTA_RESPIRASOME_ASSEMBLY` | 4 | MitoCarta pathway |
| `MITOCARTA_VITAMIN_B12_METABOLISM` | 4 | MitoCarta pathway |
| `MITOCARTA_VITAMIN_D_METABOLISM` | 4 | MitoCarta pathway |
| `MITOCARTA_CATECHOL_METABOLISM` | 3 | MitoCarta pathway |
| `MITOCARTA_CREATINE_METABOLISM` | 3 | MitoCarta pathway |
| `MITOCARTA_CYTOCHROME_C` | 3 | MitoCarta pathway |
| `MITOCARTA_EICOSANOID_METABOLISM` | 3 | MitoCarta pathway |
| `MITOCARTA_FMET_PROCESSING` | 3 | MitoCarta pathway |
| `MITOCARTA_KYNURENINE_METABOLISM` | 3 | MitoCarta pathway |
| `MITOCARTA_MIA40` | 3 | MitoCarta pathway |
| `MITOCARTA_MTDNA_STABILITY_AND_DECAY` | 3 | MitoCarta pathway |
| `MITOCARTA_SAM` | 3 | MitoCarta pathway |
| `MITOCARTA_TETRAHYDROBIOPTERIN_SYNTHESIS` | 3 | MitoCarta pathway |
| `MITOCARTA_TIM22_CARRIER_PATHWAY` | 3 | MitoCarta pathway |
| `MITOCARTA_VITAMIN_B2_METABOLISM` | 3 | MitoCarta pathway |
| `MITOCARTA_CHOLESTEROL_ASSOCIATED` | 2 | MitoCarta pathway |
| `MITOCARTA_VITAMIN_B6_METABOLISM` | 2 | MitoCarta pathway |
| `MITOCARTA_VITAMIN_C_METABOLISM` | 2 | MitoCarta pathway |
| `MITOCARTA_GLYCEROL_PHOSPHATE_SHUTTLE` | 1 | MitoCarta pathway |
| `MITOCARTA_MTDNA_MODIFICATIONS` | 1 | MitoCarta pathway |
| `MITOCARTA_OXA` | 1 | MitoCarta pathway |
| `MITOCARTA_VITAMIN_B1_METABOLISM` | 1 | MitoCarta pathway |

