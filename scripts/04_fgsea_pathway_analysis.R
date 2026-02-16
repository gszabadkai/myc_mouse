# scripts/04_fgsea_pathway_analysis.R
# =============================================================================
# fGSEA Pathway Analysis
# =============================================================================
#
# Purpose: Perform gene set enrichment analysis to address biological questions:
#   Q1. How do Myc+ tumours change from 6W to 12W?
#   Q2. Is a change Myc-specific vs developmental?
#   Q3. What's the baseline developmental effect?
#
# Strategy:
#   - Run fGSEA on 12W_pos vs 6W_pos (Q1: Myc+ progression)
#   - Run fGSEA on 12W_neg vs 6W_neg (Q3: developmental baseline)
#   - Compare enrichment between comparisons to identify Myc-specific effects (Q2)
#
# Ranking metric: Wald statistic
#   We use the Wald statistic (log2FoldChange / lfcSE) from DESeq2 for ranking:
#   1. It's a single quantity directly from the model, not a derived metric
#   2. Already incorporates both effect size and precision
#   3. Signed, so preserves direction of change
#   4. Compared to sign(LFC) × -log10(p), Wald is slightly more sensitive
#      for detecting pathways with moderate but consistent effects (r = 0.95
#      correlation between metrics, but Wald detected 11 additional pathways)
#
# Note: We also keep shrunken LFC results available for later visualisation.
#
# Gene sets analysed (89 total):
#   - MitoCarta pathways (22 sets)
#   - MYC signatures (17 sets)
#   - Apoptosis pathways (2 sets)
#   - MSigDB Hallmark pathways (50 sets)
#
# =============================================================================

source(here::here("scripts", "00_setup_packages.R"))

# === Create output directory ===
fgsea_dir <- here("outputs", "fgsea")
dir.create(fgsea_dir, showWarnings = FALSE, recursive = TRUE)

# === Load data ===
group_results <- readRDS(here("results", "group_results.rds"))
gene_sets <- readRDS(here("results", "gene_sets_list.rds"))
ortholog_table <- readRDS(here("results", "ortholog_table.rds"))

# Use raw results for Wald statistic ranking
res_myc_pos <- group_results$pos_12W_vs_6W_raw
res_myc_neg <- group_results$neg_12W_vs_6W_raw

# Create Ensembl to gene symbol mapping
ensembl_to_symbol <- setNames(
  ortholog_table$external_gene_name,
  ortholog_table$ensembl_gene_id
)

# =============================================================================
# PART 1: PREPARE RANKED GENE LISTS
# =============================================================================

#' Create a ranking metric for GSEA using Wald statistic
#' Wald = log2FoldChange / lfcSE
#' Already signed, incorporates both effect size and precision
#'
#' @param res DESeq2 results object (raw, not shrunk)
#' @param id_to_symbol Named vector mapping Ensembl IDs to gene symbols
#' @return Named numeric vector of ranks (named by gene symbol)
create_ranks <- function(res, id_to_symbol) {
  res_df <- as.data.frame(res) |>
    rownames_to_column("ensembl_id") |>
    filter(!is.na(stat)) |>
    mutate(
      gene_symbol = id_to_symbol[ensembl_id],
      rank = stat  # Wald statistic
    ) |>
    filter(!is.na(gene_symbol) & gene_symbol != "") |>
    # Handle duplicate symbols by keeping highest absolute rank
    group_by(gene_symbol) |>
    slice_max(abs(rank), n = 1, with_ties = FALSE) |>
    ungroup() |>
    arrange(desc(rank))
  
  setNames(res_df$rank, res_df$gene_symbol)
}

ranks_myc_pos <- create_ranks(res_myc_pos, ensembl_to_symbol)
ranks_myc_neg <- create_ranks(res_myc_neg, ensembl_to_symbol)

message("Ranked gene lists created (using Wald statistic):")
message("  Myc+ (12W vs 6W): ", length(ranks_myc_pos), " genes")
message("  Myc- (12W vs 6W): ", length(ranks_myc_neg), " genes")

# Verify gene name overlap with gene sets
example_overlap <- sum(names(ranks_myc_pos) %in% gene_sets[["MSigDB_HALLMARK_MYC_TARGETS_V1"]])
message("  Example overlap (MYC_TARGETS_V1): ", example_overlap, " / ", 
        length(gene_sets[["MSigDB_HALLMARK_MYC_TARGETS_V1"]]), " genes")

# =============================================================================
# PART 2: RUN fGSEA
# =============================================================================

set.seed(42)

message("\nRunning fGSEA...")

fgsea_myc_pos <- fgsea(
  pathways = gene_sets,
  stats = ranks_myc_pos,
  minSize = 10,
  maxSize = 500,
  nPermSimple = 10000
)

fgsea_myc_neg <- fgsea(
  pathways = gene_sets,
  stats = ranks_myc_neg,
  minSize = 10,
  maxSize = 500,
  nPermSimple = 10000
)

message("fGSEA complete:")
message("  Myc+ significant (padj < 0.05): ", sum(fgsea_myc_pos$padj < 0.05, na.rm = TRUE), " pathways")
message("  Myc- significant (padj < 0.05): ", sum(fgsea_myc_neg$padj < 0.05, na.rm = TRUE), " pathways")

# =============================================================================
# PART 3: COMBINE RESULTS FOR COMPARISON (Q2: Myc-specific effects)
# =============================================================================

fgsea_combined <- fgsea_myc_pos |>
  dplyr::select(pathway, NES_pos = NES, padj_pos = padj, size) |>
  left_join(
    fgsea_myc_neg |> dplyr::select(pathway, NES_neg = NES, padj_neg = padj),
    by = "pathway"
  ) |>
  mutate(
    NES_diff = NES_pos - NES_neg,
    sig_pos = padj_pos < 0.05,
    sig_neg = padj_neg < 0.05,
    category = case_when(
      sig_pos & !sig_neg ~ "Myc+ specific",
      !sig_pos & sig_neg ~ "Developmental only",
      sig_pos & sig_neg & sign(NES_pos) == sign(NES_neg) ~ "Shared (same direction)",
      sig_pos & sig_neg & sign(NES_pos) != sign(NES_neg) ~ "Opposite effects",
      TRUE ~ "Not significant"
    )
  ) |>
  arrange(desc(abs(NES_diff)))

message("\nPathway classification summary:")
print(table(fgsea_combined$category))

# =============================================================================
# PART 4: EXTRACT PATHWAY CATEGORIES FOR FOCUSED REPORTING
# =============================================================================

filter_pathways <- function(df, pattern) {
  df |> filter(grepl(pattern, pathway, ignore.case = TRUE))
}

myc_pathways <- filter_pathways(fgsea_combined, "MYC")
mito_pathways <- filter_pathways(fgsea_combined, "^MC_")
apoptosis_pathways <- filter_pathways(fgsea_combined, "apoptosis")
hallmark_pathways <- filter_pathways(fgsea_combined, "^MSigDB_HALLMARK")

# =============================================================================
# PART 5: SAVE RESULTS
# =============================================================================

message("\nSaving results...")

fgsea_results <- list(
  myc_pos = fgsea_myc_pos,
  myc_neg = fgsea_myc_neg,
  combined = fgsea_combined,
  pathways_myc = myc_pathways,
  pathways_mito = mito_pathways,
  pathways_apoptosis = apoptosis_pathways,
  pathways_hallmark = hallmark_pathways,
  ranks_myc_pos = ranks_myc_pos,
  ranks_myc_neg = ranks_myc_neg,
  analysis_date = Sys.Date(),
  gene_set_source = "gene_sets_list.rds (MitoCarta + MYC + Apoptosis + Hallmark)"
)

saveRDS(fgsea_results, here("results", "fgsea_results.rds"))
message("Results saved to: results/fgsea_results.rds")

# =============================================================================
# PART 6: SUMMARY TABLES FOR KEY PATHWAYS
# =============================================================================

message("\n", strrep("=", 70))
message("KEY RESULTS SUMMARY")
message(strrep("=", 70))

message("\n--- MYC Signatures ---")
myc_pathways |>
  dplyr::select(pathway, NES_pos, padj_pos, NES_neg, padj_neg, category) |>
  arrange(padj_pos) |
  head(20) |
  print()

message("\n--- MitoCarta Pathways ---")
mito_pathways |>
  dplyr::select(pathway, NES_pos, padj_pos, NES_neg, padj_neg, category) |>
  arrange(padj_pos) |
  head(25) |
  print()

message("\n--- Apoptosis Pathways ---")
apoptosis_pathways |>
  dplyr::select(pathway, NES_pos, padj_pos, NES_neg, padj_neg, category) |>
  print()

message("\n--- Top Hallmark Pathways (by |NES_diff|) ---")
hallmark_pathways |>
  dplyr::select(pathway, NES_pos, padj_pos, NES_neg, padj_neg, NES_diff, category) |>
  arrange(desc(abs(NES_diff))) |
  head(15) |
  print()

message("\n", strrep("=", 70))
message("fGSEA PATHWAY ANALYSIS COMPLETE")
message(strrep("=", 70))
