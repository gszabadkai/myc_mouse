# scripts/06_fgsea_cross_sectional.R
# =============================================================================
# fGSEA Cross-Sectional Analysis: Myc+ vs Myc- at each timepoint
# =============================================================================
#
# Purpose: Understand ABSOLUTE pathway differences between Myc+ and Myc- at
#   each timepoint, complementing the temporal analysis in script 04.
#
# Questions addressed:
#   Q1. What pathways are enriched in Myc+ vs Myc- at 6W (early)?
#   Q2. What pathways are enriched in Myc+ vs Myc- at 12W (late)?
#   Q3. Which pathway differences are stable vs changing over time?
#
# Strategy:
#   - Use existing DESeq2 contrasts: myc_6W_raw (Myc effect at 6W)
#     and myc_12W_raw (Myc effect at 12W) from the interaction model
#   - Rank by Wald statistic (same as script 04)
#   - Classify pathways by significance pattern across timepoints
#
# Contrast interpretation:
#   Positive NES = enriched in Myc+ relative to Myc-
#   Negative NES = depleted in Myc+ relative to Myc-
#
# =============================================================================

source(here::here("scripts", "00_setup_packages.R"))

# === Create output directory ===
fgsea_xs_dir <- here("outputs", "fgsea_cross_sectional")
dir.create(fgsea_xs_dir, showWarnings = FALSE, recursive = TRUE)

# === Load data ===
interaction_results <- readRDS(here("results", "interaction_results.rds"))
gene_sets <- readRDS(here("results", "gene_sets_list.rds"))
ortholog_table <- readRDS(here("results", "ortholog_table.rds"))

# Cross-sectional contrasts: Myc+ vs Myc- at each timepoint
res_myc_at_6W <- interaction_results$myc_6W_raw
res_myc_at_12W <- interaction_results$myc_12W_raw

# Ensembl to symbol mapping
ensembl_to_symbol <- setNames(
  ortholog_table$external_gene_name,
  ortholog_table$ensembl_gene_id
)

# =============================================================================
# PART 1: PREPARE RANKED GENE LISTS
# =============================================================================

#' Create a ranking metric for GSEA using Wald statistic
#' (same approach as script 04)
create_ranks <- function(res, id_to_symbol) {
  res_df <- as.data.frame(res) |>
    rownames_to_column("ensembl_id") |>
    filter(!is.na(stat)) |>
    mutate(
      gene_symbol = id_to_symbol[ensembl_id],
      rank = stat
    ) |>
    filter(!is.na(gene_symbol) & gene_symbol != "") |>
    group_by(gene_symbol) |>
    slice_max(abs(rank), n = 1, with_ties = FALSE) |>
    ungroup() |>
    arrange(desc(rank))

  setNames(res_df$rank, res_df$gene_symbol)
}

ranks_6W <- create_ranks(res_myc_at_6W, ensembl_to_symbol)
ranks_12W <- create_ranks(res_myc_at_12W, ensembl_to_symbol)

message("Ranked gene lists created (Myc+ vs Myc-, Wald statistic):")
message("  6W:  ", length(ranks_6W), " genes")
message("  12W: ", length(ranks_12W), " genes")

# =============================================================================
# PART 2: RUN fGSEA
# =============================================================================

set.seed(42)

message("\nRunning fGSEA (cross-sectional)...")

fgsea_6W <- fgsea(
  pathways = gene_sets,
  stats = ranks_6W,
  minSize = 10,
  maxSize = 500,
  nPermSimple = 10000
)

fgsea_12W <- fgsea(
  pathways = gene_sets,
  stats = ranks_12W,
  minSize = 10,
  maxSize = 500,
  nPermSimple = 10000
)

message("fGSEA complete:")
message("  6W  significant (padj < 0.05): ",
        sum(fgsea_6W$padj < 0.05, na.rm = TRUE), " pathways")
message("  12W significant (padj < 0.05): ",
        sum(fgsea_12W$padj < 0.05, na.rm = TRUE), " pathways")

# =============================================================================
# PART 3: COMBINE AND CLASSIFY
# =============================================================================

fgsea_xs_combined <- fgsea_6W |>
  dplyr::select(pathway, NES_6W = NES, padj_6W = padj, size) |>
  left_join(
    fgsea_12W |> dplyr::select(pathway, NES_12W = NES, padj_12W = padj),
    by = "pathway"
  ) |>
  mutate(
    NES_diff = NES_12W - NES_6W,
    sig_6W = padj_6W < 0.05,
    sig_12W = padj_12W < 0.05,
    category = case_when(
      sig_6W & sig_12W & sign(NES_6W) == sign(NES_12W) ~ "Stable (both timepoints)",
      sig_6W & sig_12W & sign(NES_6W) != sign(NES_12W) ~ "Reversed",
      sig_6W & !sig_12W ~ "6W only (lost by 12W)",
      !sig_6W & sig_12W ~ "12W only (gained by 12W)",
      TRUE ~ "Not significant"
    )
  ) |>
  arrange(desc(abs(NES_diff)))

message("\nCross-sectional pathway classification:")
print(table(fgsea_xs_combined$category))

# =============================================================================
# PART 4: EXTRACT PATHWAY CATEGORIES
# =============================================================================

filter_pathways <- function(df, pattern) {
  df |> filter(grepl(pattern, pathway, ignore.case = TRUE))
}

xs_myc_pathways <- filter_pathways(fgsea_xs_combined, "MYC")
xs_mito_pathways <- filter_pathways(fgsea_xs_combined, "^MC_")
xs_apoptosis_pathways <- filter_pathways(fgsea_xs_combined, "apoptosis")
xs_hallmark_pathways <- filter_pathways(fgsea_xs_combined, "^MSigDB_HALLMARK")

# =============================================================================
# PART 5: SAVE RESULTS
# =============================================================================

message("\nSaving results...")

fgsea_xs_results <- list(
  fgsea_6W = fgsea_6W,
  fgsea_12W = fgsea_12W,
  combined = fgsea_xs_combined,
  pathways_myc = xs_myc_pathways,
  pathways_mito = xs_mito_pathways,
  pathways_apoptosis = xs_apoptosis_pathways,
  pathways_hallmark = xs_hallmark_pathways,
  ranks_6W = ranks_6W,
  ranks_12W = ranks_12W,
  analysis_date = Sys.Date(),
  description = "Cross-sectional fGSEA: Myc+ vs Myc- at 6W and 12W"
)

saveRDS(fgsea_xs_results, here("results", "fgsea_xs_results.rds"))
message("Results saved to: results/fgsea_xs_results.rds")

# =============================================================================
# PART 6: SUMMARY TABLES
# =============================================================================

message("\n", strrep("=", 70))
message("CROSS-SECTIONAL fGSEA RESULTS SUMMARY")
message(strrep("=", 70))

message("\n--- MYC Signatures (Myc+ vs Myc-) ---")
xs_myc_pathways |>
  dplyr::select(pathway, NES_6W, padj_6W, NES_12W, padj_12W, category) |>
  arrange(padj_6W) |>
  print()

message("\n--- MitoCarta Pathways ---")
xs_mito_pathways |>
  dplyr::select(pathway, NES_6W, padj_6W, NES_12W, padj_12W, category) |>
  arrange(padj_6W) |>
  head(15) |>
  print()

message("\n--- Apoptosis Pathways ---")
xs_apoptosis_pathways |>
  dplyr::select(pathway, NES_6W, padj_6W, NES_12W, padj_12W, category) |>
  print()

message("\n--- Top Hallmark Pathways (by |NES_diff|) ---")
xs_hallmark_pathways |>
  dplyr::select(pathway, NES_6W, padj_6W, NES_12W, padj_12W, NES_diff, category) |>
  arrange(desc(abs(NES_diff))) |>
  head(15) |>
  print()

message("\n", strrep("=", 70))
message("CROSS-SECTIONAL fGSEA COMPLETE")
message(strrep("=", 70))
