# scripts/11_interaction_gene_characterisation.R
# =============================================================================
# Characterisation of Interaction Gene Sets: Gene-Level Dissection of the
# Fading Myc Effect
# =============================================================================
#
# Background:
#   Script 03 identified three classes of temporally regulated genes based on
#   whether their timepoint effect (padj < 0.1) is significant in Myc+, Myc−,
#   or both:
#     - Unique to Myc+:  ~1,662 genes (temporal change only in Myc+ samples)
#     - Unique to Myc−:  ~906 genes  (temporal change only in Myc− samples)
#     - Overlapping:     ~961 genes  (temporal change in both genotypes)
#
#   Among the unique-to-Myc+ set, ~357 genes have nominal interaction p < 0.05,
#   well above the ~83 expected by chance — suggesting a weak but real signal
#   that fails FDR correction. Script 10 confirmed pathway-level interaction
#   signal via fGSEA (52/88 curated sets, 40/142 MitoCarta pathways).
#
# This script asks:
#   1. DIRECTION (Part 2): Are the interaction genes predominantly negative
#      (Myc effect weakening at 12W) or mixed?
#
#   2. FUNCTIONAL ENRICHMENT (Part 3): What biological processes are enriched?
#      Two complementary g:Profiler approaches:
#        a) Ordered query on all unique-to-Myc+ genes ranked by |Wald stat|
#           (more powerful, uses the continuous ranking)
#        b) Unordered ORA on the 357 p < 0.05 genes, split by direction
#      Also run on unique-to-Myc− and overlapping sets for comparison.
#
#   3. LEADING EDGE OVERLAP (Part 4): Do these 357 genes drive the significant
#      fGSEA interaction pathways from script 10? Hypergeometric test of
#      enrichment within leading edges.
#
#   4. SELECTION vs REMODELLING (Part 5B): If Bbc3/PUMA (a pro-apoptotic gene)
#      shows a fading Myc effect, is this clonal selection against apoptosis-
#      prone cells, or transcriptional remodelling of Myc's promoter targets?
#      Selection predicts a coordinated shift across the entire pro-apoptotic
#      module; remodelling predicts isolated gene-level changes. Tested with:
#        - Forest plots: Apoptosis-PRO (25 genes), Apoptosis-ANTI (9 genes)
#        - Volcano plots: all MitoCarta genes, Felsher integrative MYC signature
#        - One-sample t-tests for group-level shifts
#
#   5. THREE-SET COMPARISON (Part 5C): Parallel analysis of the unique-to-Myc−
#      and overlapping gene sets. The overlapping set serves as a negative
#      control (shared developmental genes should show minimal interaction).
#      The unique-to-Myc− set tests whether Myc buffers developmental changes
#      that occur in controls. Comparative summary table quantifies interaction
#      signal enrichment, directional bias, and functional themes across all
#      three sets.
#
# Input:
#   - results/interaction_results.rds       (DESeq2 results, all contrasts)
#   - results/interaction_fgsea_mitopps.rds (script 10 fGSEA + mitoPPS interaction)
#   - results/mitopps_scores.rds            (gene-to-pathway mapping)
#   - results/ortholog_table.rds
#   - results/gene_sets_list.rds
#
# Output:
#   - results/interaction_gene_characterisation.rds
#   - outputs/interaction_analysis/:
#       Direction:    histogram_interaction_gene_direction.pdf
#                     histogram_interaction_gene_basemean.pdf
#                     volcano_interaction_unique_to_pos.pdf
#       g:Profiler:   gprofiler_manhattan_ordered.pdf
#                     gprofiler_manhattan_unordered.pdf
#                     dotplot_gprofiler_ordered.pdf
#                     gprofiler_ordered_sig.csv
#                     gprofiler_unordered_sig.csv
#       Leading edge: dotplot_le_overlap_curated.pdf
#                     dotplot_le_overlap_mitocarta.pdf
#                     le_overlap_curated.csv, le_overlap_mitocarta.csv
#       Selection:    focused_interaction_selection_vs_remodelling.pdf (4-panel)
#                     forest_apoptosis_pro_interaction.pdf
#                     forest_apoptosis_anti_interaction.pdf
#                     volcano_mitocarta_interaction.pdf
#                     volcano_felsher_myc_interaction.pdf
#       Three-set:    histogram_interaction_lfc_three_sets.pdf
#                     volcano_interaction_three_sets.pdf (3-panel)
#                     volcano_interaction_unique_to_neg.pdf
#                     volcano_interaction_overlapping.pdf
#                     interaction_three_sets_summary.csv
#                     interaction_sig_genes_p05.csv
#
# =============================================================================

source(here::here("scripts", "00_setup_packages.R"))

# === Output directory (shared with script 10) ===
int_dir <- here("outputs", "interaction_analysis")
dir.create(int_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# PART 1: LOAD DATA AND RECONSTRUCT THE GENE SET
# =============================================================================

message("Loading data...")

interaction_results <- readRDS(here("results", "interaction_results.rds"))
interaction_analysis <- readRDS(here("results", "interaction_fgsea_mitopps.rds"))
ortholog_table <- readRDS(here("results", "ortholog_table.rds"))
gene_sets <- readRDS(here("results", "gene_sets_list.rds"))

# Ensembl → symbol mapping
ensembl_to_symbol <- setNames(
  ortholog_table$external_gene_name,
  ortholog_table$ensembl_gene_id
)

# --- Reconstruct "unique to Myc+" genes (from script 03 logic) ---
sig_timepoint_pos <- interaction_results$timepoint_pos_raw |>
  as.data.frame() |>
  filter(padj < 0.1) |>
  rownames()

sig_timepoint_neg <- interaction_results$timepoint_neg_raw |>
  as.data.frame() |>
  filter(padj < 0.1) |>
  rownames()

unique_to_pos <- setdiff(sig_timepoint_pos, sig_timepoint_neg)

message(sprintf("Genes with significant timepoint effect only in Myc+: %d", length(unique_to_pos)))

# --- Extract interaction statistics for these genes ---
interaction_raw_df <- interaction_results$interaction_raw |>
  as.data.frame() |>
  rownames_to_column("ensembl_id") |>
  mutate(gene_symbol = ensembl_to_symbol[ensembl_id])

unique_pos_interaction <- interaction_raw_df |>
  filter(ensembl_id %in% unique_to_pos)

# --- The key gene set: nominal interaction p < 0.05 ---
int_sig_genes <- unique_pos_interaction |>
  filter(pvalue < 0.05)

message(sprintf("Of these, %d have interaction p < 0.05 (expected by chance: %d)",
                nrow(int_sig_genes),
                round(length(unique_to_pos) * 0.05)))

# =============================================================================
# PART 2: DIRECTION ANALYSIS
# =============================================================================
# The interaction LFC sign tells us:
#   Negative = Myc effect weakens at 12W vs 6W (Myc drives gene UP at 6W,
#              but the UP-regulation fades or reverses at 12W)
#   Positive = Myc effect strengthens at 12W

message("\n", strrep("=", 70))
message("PART 2: DIRECTION ANALYSIS")
message(strrep("=", 70))

n_neg <- sum(int_sig_genes$log2FoldChange < 0)
n_pos <- sum(int_sig_genes$log2FoldChange > 0)
pct_neg <- round(100 * n_neg / nrow(int_sig_genes), 1)

message(sprintf("\nInteraction direction for %d genes (p < 0.05):", nrow(int_sig_genes)))
message(sprintf("  Negative (Myc effect weakens at 12W): %d (%.1f%%)", n_neg, pct_neg))
message(sprintf("  Positive (Myc effect strengthens at 12W): %d (%.1f%%)", n_pos, 100 - pct_neg))

# Binomial test: is the direction bias significant?
binom_result <- binom.test(n_neg, nrow(int_sig_genes), p = 0.5)
message(sprintf("  Binomial test p-value: %.2e", binom_result$p.value))

# Summary statistics on LFC magnitude
message(sprintf("\n  Interaction LFC — median: %.3f, mean: %.3f, range: [%.3f, %.3f]",
                median(int_sig_genes$log2FoldChange),
                mean(int_sig_genes$log2FoldChange),
                min(int_sig_genes$log2FoldChange),
                max(int_sig_genes$log2FoldChange)))

# Also characterise the broader set (all 1,662 unique_to_pos)
n_neg_all <- sum(unique_pos_interaction$log2FoldChange < 0, na.rm = TRUE)
n_pos_all <- sum(unique_pos_interaction$log2FoldChange > 0, na.rm = TRUE)
message(sprintf("\n  For context — all %d unique-to-Myc+ genes:",
                nrow(unique_pos_interaction)))
message(sprintf("  Negative interaction LFC: %d (%.1f%%)", n_neg_all,
                100 * n_neg_all / nrow(unique_pos_interaction)))

# --- Direction analysis plot ---
p_direction <- ggplot(int_sig_genes, aes(x = log2FoldChange)) +
  geom_histogram(bins = 40, fill = "steelblue", colour = "white", boundary = 0) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "red", linewidth = 0.6) +
  geom_vline(xintercept = median(int_sig_genes$log2FoldChange),
             linetype = "dotted", colour = "darkblue", linewidth = 0.6) +
  annotate("text", x = Inf, y = Inf,
           label = sprintf("n = %d\n%d negative (%.0f%%)\nMedian LFC = %.3f\nBinomial p = %.2e",
                           nrow(int_sig_genes), n_neg, pct_neg,
                           median(int_sig_genes$log2FoldChange),
                           binom_result$p.value),
           hjust = 1.1, vjust = 1.3, size = 3.5) +
  labs(
    title = "Interaction LFC direction: genes unique to Myc+ timepoint (p < 0.05)",
    subtitle = paste0(
      "Negative LFC = Myc effect weakens at 12W  |  ",
      "Positive LFC = Myc effect strengthens at 12W\n",
      "Red dashed = zero  |  Blue dotted = median"
    ),
    x = "Interaction log2FC (timepoint12W:myc_statuspos)",
    y = "Count"
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(int_dir, "histogram_interaction_gene_direction.pdf"),
       p_direction, width = 8, height = 5)

# --- Expression level of these genes ---
# Are they predominantly low- or high-expression?
p_basemean <- ggplot(unique_pos_interaction,
                     aes(x = log10(baseMean + 1), fill = pvalue < 0.05)) +
  geom_histogram(bins = 50, colour = "white", alpha = 0.7, position = "identity") +
  scale_fill_manual(
    values = c("TRUE" = "#B2182B", "FALSE" = "grey70"),
    labels = c("TRUE" = "p < 0.05 (interaction signal)", "FALSE" = "p >= 0.05"),
    name = "Interaction"
  ) +
  labs(
    title = "Expression level of unique-to-Myc+ genes",
    subtitle = "Red: genes with nominal interaction signal  |  Grey: non-significant",
    x = "log10(baseMean + 1)",
    y = "Count"
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "bottom")

ggsave(file.path(int_dir, "histogram_interaction_gene_basemean.pdf"),
       p_basemean, width = 8, height = 5)

# =============================================================================
# PART 3: FUNCTIONAL ENRICHMENT WITH g:Profiler
# =============================================================================
# g:Profiler performs over-representation analysis against multiple databases
# simultaneously: GO (BP, MF, CC), KEGG, Reactome, WikiPathways, TRANSFAC,
# miRTarBase, HP, CORUM.
#
# Two complementary approaches:
#
# A) ORDERED QUERY on all ~1,662 unique-to-Myc+ genes ranked by interaction
#    Wald statistic (absolute value, descending). g:Profiler's ordered mode
#    applies a minimum hypergeometric (mHG) test: for each functional term,
#    it walks down the ranked list and finds the rank threshold that gives the
#    strongest enrichment. This is strictly more powerful than a fixed-cutoff
#    ORA because it uses the continuous ranking rather than an arbitrary binary
#    split, and can detect enrichment driven by genes in the p = 0.05–0.2
#    range that a hard cutoff would miss.
#
# B) UNORDERED ORA on the 357 genes (interaction p < 0.05), split by direction.
#    This directly characterises the specific gene set identified in script 03
#    and tests whether positive- and negative-direction genes map to different
#    biology.
#
# If the same terms emerge from both approaches, that's strong convergent
# evidence. Terms unique to the ordered query are driven by sub-threshold genes.

message("\n", strrep("=", 70))
message("PART 3: FUNCTIONAL ENRICHMENT (g:Profiler)")
message(strrep("=", 70))

# --- Prepare gene lists ---

# Map all unique-to-Myc+ genes to symbols with interaction Wald stats
unique_pos_with_stats <- interaction_raw_df |>
  filter(ensembl_id %in% unique_to_pos,
         !is.na(gene_symbol), gene_symbol != "",
         !is.na(stat)) |>
  # Deduplicate symbols: keep strongest interaction signal
  group_by(gene_symbol) |>
  slice_max(abs(stat), n = 1, with_ties = FALSE) |>
  ungroup()

# Ordered list: rank by absolute Wald statistic (strongest interaction first)
# Using absolute value because g:Profiler ordered mode tests enrichment at
# the TOP of the list — we want genes with the strongest interaction
# (regardless of direction) ranked first
ordered_symbols <- unique_pos_with_stats |>
  arrange(desc(abs(stat))) |>
  pull(gene_symbol)

# The 357 gene set (p < 0.05)
int_sig_symbols <- int_sig_genes |>
  filter(!is.na(gene_symbol) & gene_symbol != "") |>
  pull(gene_symbol) |>
  unique()

# Background: all genes in the DESeq2 model, mapped to symbols
all_tested_ensembl <- rownames(interaction_results$interaction_raw)
all_tested_symbols <- ensembl_to_symbol[all_tested_ensembl]
all_tested_symbols <- unique(na.omit(all_tested_symbols[all_tested_symbols != ""]))

# Direction subsets
int_sig_neg_symbols <- int_sig_genes |>
  filter(log2FoldChange < 0, !is.na(gene_symbol), gene_symbol != "") |>
  pull(gene_symbol) |> unique()

int_sig_pos_symbols <- int_sig_genes |>
  filter(log2FoldChange > 0, !is.na(gene_symbol), gene_symbol != "") |>
  pull(gene_symbol) |> unique()

message(sprintf("  Ordered query: %d unique symbols (all unique-to-Myc+, ranked by |Wald|)",
                length(ordered_symbols)))
message(sprintf("  Unordered query: %d symbols (interaction p < 0.05)",
                length(int_sig_symbols)))
message(sprintf("    Negative direction: %d", length(int_sig_neg_symbols)))
message(sprintf("    Positive direction: %d", length(int_sig_pos_symbols)))
message(sprintf("  Background: %d unique gene symbols", length(all_tested_symbols)))

# -------------------------------------------------------------------------
# 3A. ORDERED QUERY: full ranked list of unique-to-Myc+ genes
# -------------------------------------------------------------------------

message("\n--- 3A: Ordered query (all unique-to-Myc+ genes, ranked by |Wald|) ---")

gost_ordered <- gost(
  query            = ordered_symbols,
  organism         = "mmusculus",
  ordered_query    = TRUE,
  significant      = FALSE,
  custom_bg        = all_tested_symbols,
  sources          = c("GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC", "WP"),
  evcodes          = TRUE
)

if (!is.null(gost_ordered) && nrow(gost_ordered$result) > 0) {
  n_sig_ordered <- sum(gost_ordered$result$p_value < 0.05)
  message(sprintf("  Significant terms (g:SCS p < 0.05): %d out of %d tested",
                  n_sig_ordered, nrow(gost_ordered$result)))

  message("\n  Top 20 enriched terms (ordered query):")
  gost_ordered$result |>
    filter(p_value < 0.05) |>
    arrange(p_value) |>
    dplyr::select(source, term_name, term_size, intersection_size, p_value) |>
    head(20) |>
    print()
} else {
  message("  g:Profiler ordered query: no results returned")
}

# -------------------------------------------------------------------------
# 3B. UNORDERED ORA: 357 genes (all, then split by direction)
# -------------------------------------------------------------------------

message("\n--- 3B: Unordered ORA (357 genes with interaction p < 0.05) ---")

# All 357 combined
gost_all <- gost(
  query            = int_sig_symbols,
  organism         = "mmusculus",
  ordered_query    = FALSE,
  significant      = FALSE,
  custom_bg        = all_tested_symbols,
  sources          = c("GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC", "WP"),
  evcodes          = TRUE
)

if (!is.null(gost_all) && nrow(gost_all$result) > 0) {
  n_sig_gost <- sum(gost_all$result$p_value < 0.05)
  message(sprintf("  Unordered (all): %d significant terms", n_sig_gost))

  message("\n  Top 20 enriched terms (all 357 genes, unordered):")
  gost_all$result |>
    filter(p_value < 0.05) |>
    arrange(p_value) |>
    dplyr::select(source, term_name, term_size, intersection_size, p_value) |>
    head(20) |>
    print()
} else {
  message("  g:Profiler unordered (all): no results returned")
}

# Negative-direction genes
gost_neg <- NULL
if (length(int_sig_neg_symbols) >= 10) {
  gost_neg <- gost(
    query            = int_sig_neg_symbols,
    organism         = "mmusculus",
    ordered_query    = FALSE,
    significant      = FALSE,
    custom_bg        = all_tested_symbols,
    sources          = c("GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC", "WP"),
    evcodes          = TRUE
  )

  if (!is.null(gost_neg) && nrow(gost_neg$result) > 0) {
    n_sig_neg <- sum(gost_neg$result$p_value < 0.05)
    message(sprintf("\n  Unordered (negative direction): %d significant terms", n_sig_neg))

    message("  Top 20 enriched terms (negative interaction — Myc effect weakens):")
    gost_neg$result |>
      filter(p_value < 0.05) |>
      arrange(p_value) |>
      dplyr::select(source, term_name, term_size, intersection_size, p_value) |>
      head(20) |>
      print()
  }
}

# Positive-direction genes
gost_pos <- NULL
if (length(int_sig_pos_symbols) >= 10) {
  gost_pos <- gost(
    query            = int_sig_pos_symbols,
    organism         = "mmusculus",
    ordered_query    = FALSE,
    significant      = FALSE,
    custom_bg        = all_tested_symbols,
    sources          = c("GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC", "WP"),
    evcodes          = TRUE
  )

  if (!is.null(gost_pos) && nrow(gost_pos$result) > 0) {
    n_sig_pos_gost <- sum(gost_pos$result$p_value < 0.05)
    message(sprintf("\n  Unordered (positive direction): %d significant terms", n_sig_pos_gost))

    message("  Top 20 enriched terms (positive interaction — Myc effect strengthens):")
    gost_pos$result |>
      filter(p_value < 0.05) |>
      arrange(p_value) |>
      dplyr::select(source, term_name, term_size, intersection_size, p_value) |>
      head(20) |>
      print()
  }
} else {
  message(sprintf("\n  Positive direction: only %d genes — skipping g:Profiler",
                  length(int_sig_pos_symbols)))
}

# -------------------------------------------------------------------------
# 3C. COMPARE ORDERED vs UNORDERED results
# -------------------------------------------------------------------------

message("\n--- 3C: Comparing ordered vs unordered results ---")

if (!is.null(gost_ordered) && !is.null(gost_all) &&
    nrow(gost_ordered$result) > 0 && nrow(gost_all$result) > 0) {

  sig_ordered_terms <- gost_ordered$result |>
    filter(p_value < 0.05) |>
    pull(term_id)
  sig_unordered_terms <- gost_all$result |>
    filter(p_value < 0.05) |>
    pull(term_id)

  shared   <- intersect(sig_ordered_terms, sig_unordered_terms)
  only_ord <- setdiff(sig_ordered_terms, sig_unordered_terms)
  only_uno <- setdiff(sig_unordered_terms, sig_ordered_terms)

  message(sprintf("  Significant terms (p < 0.05):"))
  message(sprintf("    Ordered:   %d", length(sig_ordered_terms)))
  message(sprintf("    Unordered: %d", length(sig_unordered_terms)))
  message(sprintf("    Shared:    %d", length(shared)))
  message(sprintf("    Ordered only:   %d (driven by sub-threshold genes)", length(only_ord)))
  message(sprintf("    Unordered only: %d", length(only_uno)))

  if (length(only_ord) > 0) {
    message("\n  Terms found ONLY in ordered query (top 10):")
    gost_ordered$result |>
      filter(term_id %in% only_ord, p_value < 0.05) |>
      arrange(p_value) |>
      dplyr::select(source, term_name, term_size, intersection_size, p_value) |>
      head(10) |>
      print()
  }
}

# --- g:Profiler visualisations ---

# Manhattan plot — ordered query (primary result)
if (!is.null(gost_ordered) && nrow(filter(gost_ordered$result, p_value < 0.05)) > 0) {
  p_gost_manhattan_ord <- gostplot(gost_ordered, capped = TRUE, interactive = FALSE)
  ggsave(file.path(int_dir, "gprofiler_manhattan_ordered.pdf"),
         p_gost_manhattan_ord, width = 14, height = 7)
}

# Manhattan plot — unordered (for comparison)
if (!is.null(gost_all) && nrow(filter(gost_all$result, p_value < 0.05)) > 0) {
  p_gost_manhattan_uno <- gostplot(gost_all, capped = TRUE, interactive = FALSE)
  ggsave(file.path(int_dir, "gprofiler_manhattan_unordered.pdf"),
         p_gost_manhattan_uno, width = 14, height = 7)
}

# Custom dotplot of top terms from ordered query
if (!is.null(gost_ordered) && nrow(filter(gost_ordered$result, p_value < 0.05)) > 0) {
  top_terms <- gost_ordered$result |>
    filter(p_value < 0.05) |>
    arrange(p_value) |>
    group_by(source) |>
    slice_head(n = 8) |>
    ungroup() |>
    head(30)

  # Flag terms also found in unordered ORA
  if (!is.null(gost_all) && nrow(gost_all$result) > 0) {
    sig_unordered_ids <- gost_all$result |>
      filter(p_value < 0.05) |>
      pull(term_id)
    top_terms <- top_terms |>
      mutate(also_in_unordered = term_id %in% sig_unordered_ids)
  } else {
    top_terms$also_in_unordered <- FALSE
  }

  p_gost_dot <- ggplot(top_terms,
                       aes(x = -log10(p_value),
                           y = reorder(term_name, -log10(p_value)))) +
    geom_point(aes(size = intersection_size, colour = source,
                   shape = also_in_unordered),
               alpha = 0.8) +
    scale_size_continuous(range = c(2, 7), name = "Overlap") +
    scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1),
                       labels = c("TRUE" = "Also in unordered ORA",
                                  "FALSE" = "Ordered only"),
                       name = "Convergence") +
    labs(
      title = "g:Profiler: ordered query on unique-to-Myc+ genes",
      subtitle = paste0(
        length(ordered_symbols),
        " genes ranked by |interaction Wald stat| (strongest interaction first)\n",
        "Top terms per source (g:SCS p < 0.05)  |  ",
        "Filled = also significant in unordered ORA"
      ),
      x = expression(-log[10](p)),
      y = NULL,
      colour = "Source"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.y  = element_text(size = 7),
      plot.title   = element_text(face = "bold"),
      legend.position = "right"
    )

  ggsave(file.path(int_dir, "dotplot_gprofiler_ordered.pdf"),
         p_gost_dot, width = 12, height = max(6, nrow(top_terms) * 0.3))
}

# =============================================================================
# PART 4: LEADING EDGE OVERLAP WITH fGSEA INTERACTION PATHWAYS
# =============================================================================
# fGSEA stores the "leading edge" genes for each pathway — the subset of
# pathway genes that drive the enrichment signal. If our 357 interaction genes
# are concentrated in the leading edges of the significant interaction
# pathways, that confirms they ARE the molecular drivers of the pathway-level
# signal detected in script 10.

message("\n", strrep("=", 70))
message("PART 4: LEADING EDGE OVERLAP")
message(strrep("=", 70))

# --- Map 357 genes to symbols for overlap ---
int_sig_ensembl <- int_sig_genes$ensembl_id
int_sig_symbols_all <- ensembl_to_symbol[int_sig_ensembl]
int_sig_symbols_set <- unique(na.omit(int_sig_symbols_all[int_sig_symbols_all != ""]))

# --- Extract leading edges from curated fGSEA interaction ---
fgsea_curated <- interaction_analysis$fgsea_interaction_curated
fgsea_mc      <- interaction_analysis$fgsea_interaction_mitocarta

# Function to compute overlap statistics
compute_leading_edge_overlap <- function(fgsea_results, query_genes,
                                         gene_sets, padj_threshold = 0.05) {
  sig_pathways <- fgsea_results |>
    filter(padj < padj_threshold) |>
    arrange(pval)

  if (nrow(sig_pathways) == 0) {
    message("  No significant pathways at padj < ", padj_threshold)
    return(NULL)
  }

  overlap_results <- lapply(seq_len(nrow(sig_pathways)), function(i) {
    pw <- sig_pathways$pathway[i]
    le_genes <- unlist(sig_pathways$leadingEdge[i])
    pw_genes <- gene_sets[[pw]]

    n_le        <- length(le_genes)
    n_pw        <- length(pw_genes)
    le_in_query <- intersect(le_genes, query_genes)
    pw_in_query <- intersect(pw_genes, query_genes)

    # Hypergeometric test: are query genes over-represented in leading edge?
    # Drawing from pathway genes, asking if query genes are enriched in LE
    # This is a within-pathway enrichment test
    # Universe = pathway genes, success in universe = query genes in pathway,
    # sample = leading edge, success in sample = query genes in LE
    n_success_pop <- length(pw_in_query)  # query genes in pathway
    n_fail_pop    <- n_pw - n_success_pop # non-query genes in pathway
    n_drawn       <- n_le                 # leading edge size
    n_success_sample <- length(le_in_query)

    hyper_p <- NA
    if (n_pw > 0 && n_le > 0 && n_success_pop > 0) {
      hyper_p <- phyper(n_success_sample - 1, n_success_pop, n_fail_pop,
                        n_drawn, lower.tail = FALSE)
    }

    data.frame(
      pathway          = pw,
      NES              = sig_pathways$NES[i],
      padj_fgsea       = sig_pathways$padj[i],
      pathway_size     = n_pw,
      leading_edge_size = n_le,
      query_in_pathway  = n_success_pop,
      query_in_le       = n_success_sample,
      pct_le_from_query = ifelse(n_le > 0, round(100 * n_success_sample / n_le, 1), 0),
      pct_query_in_le   = ifelse(n_success_pop > 0,
                                  round(100 * n_success_sample / n_success_pop, 1), 0),
      hyper_p           = hyper_p,
      le_query_genes    = paste(le_in_query, collapse = ", "),
      stringsAsFactors  = FALSE
    )
  }) |> bind_rows()

  overlap_results |>
    mutate(hyper_padj = p.adjust(hyper_p, method = "BH"))
}

# --- Curated gene sets ---
message("\n--- Leading edge overlap: curated gene sets (padj < 0.05) ---")
le_overlap_curated <- compute_leading_edge_overlap(
  fgsea_curated, int_sig_symbols_set, gene_sets, padj_threshold = 0.05
)

if (!is.null(le_overlap_curated)) {
  message(sprintf("  %d significant pathways examined", nrow(le_overlap_curated)))
  message(sprintf("  %d show enrichment of interaction genes in LE (hyper padj < 0.05)",
                  sum(le_overlap_curated$hyper_padj < 0.05, na.rm = TRUE)))

  message("\n  Top pathways by LE overlap:")
  le_overlap_curated |>
    arrange(hyper_p) |>
    dplyr::select(pathway, NES, leading_edge_size, query_in_le,
                  pct_le_from_query, hyper_p, hyper_padj) |>
    head(20) |>
    print()
}

# --- MitoCarta gene sets ---
message("\n--- Leading edge overlap: MitoCarta pathways (padj < 0.05) ---")

# Need MitoCarta gene sets from script 10 — reconstruct from gene_to_pathway
mitopps_results <- readRDS(here("results", "mitopps_scores.rds"))
gene_to_pathway <- mitopps_results$gene_to_pathway

mitocarta_gene_sets <- gene_to_pathway |>
  group_by(Pathway) |>
  summarise(genes = list(unique(Gene)), .groups = "drop") |>
  deframe()

le_overlap_mc <- compute_leading_edge_overlap(
  fgsea_mc, int_sig_symbols_set, mitocarta_gene_sets, padj_threshold = 0.05
)

if (!is.null(le_overlap_mc)) {
  message(sprintf("  %d significant MitoCarta pathways examined", nrow(le_overlap_mc)))
  message(sprintf("  %d show enrichment of interaction genes in LE (hyper padj < 0.05)",
                  sum(le_overlap_mc$hyper_padj < 0.05, na.rm = TRUE)))

  message("\n  Top MitoCarta pathways by LE overlap:")
  le_overlap_mc |>
    arrange(hyper_p) |>
    dplyr::select(pathway, NES, leading_edge_size, query_in_le,
                  pct_le_from_query, hyper_p, hyper_padj) |>
    head(20) |>
    print()
}

# --- Global leading edge membership ---
# What fraction of the 357 genes appear in ANY significant leading edge?
all_sig_le_curated <- fgsea_curated |>
  filter(padj < 0.05) |>
  pull(leadingEdge) |>
  unlist() |>
  unique()

all_sig_le_mc <- fgsea_mc |>
  filter(padj < 0.05) |>
  pull(leadingEdge) |>
  unlist() |>
  unique()

all_sig_le <- union(all_sig_le_curated, all_sig_le_mc)

n_in_any_le <- length(intersect(int_sig_symbols_set, all_sig_le))
pct_in_le <- round(100 * n_in_any_le / length(int_sig_symbols_set), 1)

message(sprintf("\n--- Global leading edge membership ---"))
message(sprintf("  %d / %d interaction genes (%.1f%%) appear in at least one",
                n_in_any_le, length(int_sig_symbols_set), pct_in_le))
message(sprintf("  significant leading edge (curated or MitoCarta)"))

# Genes NOT in any leading edge — these may represent novel biology
not_in_le <- setdiff(int_sig_symbols_set, all_sig_le)
message(sprintf("  %d interaction genes are NOT in any significant leading edge",
                length(not_in_le)))

# =============================================================================
# PART 5: VISUALISATIONS
# =============================================================================

message("\nGenerating visualisations...")

# --- 5a. Leading edge enrichment dotplot (curated) ---
if (!is.null(le_overlap_curated) && nrow(le_overlap_curated) > 0) {
  le_plot_data <- le_overlap_curated |>
    filter(query_in_le > 0) |>
    mutate(
      pathway_label = str_replace_all(pathway, "^MSigDB_HALLMARK_|^MC_|^MYC_", "") |>
        str_replace_all("_", " "),
      sig_le = hyper_padj < 0.05
    )

  if (nrow(le_plot_data) > 0) {
    p_le_curated <- ggplot(le_plot_data,
                           aes(x = pct_le_from_query,
                               y = reorder(pathway_label, pct_le_from_query))) +
      geom_point(aes(size = query_in_le,
                     colour = -log10(pmax(hyper_padj, 1e-10))),
                 alpha = 0.85) +
      scale_colour_gradient(low = "grey70", high = "#B2182B",
                            name = expression(-log[10](padj[hyper]))) +
      scale_size_continuous(range = c(2, 7), name = "Interaction\ngenes in LE") +
      geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
      labs(
        title = "Interaction genes in fGSEA leading edges (curated sets)",
        subtitle = paste0(
          "% of each pathway's leading edge composed of interaction genes\n",
          "Colour: hypergeometric enrichment significance"
        ),
        x = "% of leading edge from interaction gene set",
        y = NULL
      ) +
      theme_minimal(base_size = 11) +
      theme(
        axis.text.y  = element_text(size = 7),
        plot.title   = element_text(face = "bold"),
        legend.position = "right"
      )

    ggsave(file.path(int_dir, "dotplot_le_overlap_curated.pdf"),
           p_le_curated, width = 12, height = max(6, nrow(le_plot_data) * 0.3))
  }
}

# --- 5b. Leading edge enrichment dotplot (MitoCarta) ---
if (!is.null(le_overlap_mc) && nrow(le_overlap_mc) > 0) {
  le_mc_plot <- le_overlap_mc |>
    filter(query_in_le > 0) |>
    mutate(
      pathway_label = str_replace_all(pathway, "_", " "),
      tier1 = mitopps_results$pathway_tier1_map[pathway]
    )

  if (nrow(le_mc_plot) > 0) {
    tier1_colours <- c(
      "Metabolism"                                = "#2166AC",
      "OXPHOS"                                    = "#B2182B",
      "Protein import, sorting, and homeostasis"  = "#1B9E77",
      "Mitochondrial central dogma"               = "#D95F02",
      "Mitochondrial dynamics and surveillance"   = "#7570B3",
      "Small molecule transport"                  = "#E7298A",
      "Signaling"                                 = "#66A61E"
    )

    p_le_mc <- ggplot(le_mc_plot,
                      aes(x = pct_le_from_query,
                          y = reorder(pathway_label, pct_le_from_query))) +
      geom_point(aes(size = query_in_le, colour = tier1), alpha = 0.85) +
      scale_colour_manual(values = tier1_colours, na.value = "grey50",
                          name = "MitoCarta Tier 1") +
      scale_size_continuous(range = c(2, 7), name = "Interaction\ngenes in LE") +
      labs(
        title = "Interaction genes in fGSEA leading edges (MitoCarta)",
        subtitle = paste0(
          "% of each pathway's leading edge composed of interaction genes\n",
          "Pathways significant at padj < 0.05 in interaction fGSEA"
        ),
        x = "% of leading edge from interaction gene set",
        y = NULL
      ) +
      theme_minimal(base_size = 11) +
      theme(
        axis.text.y  = element_text(size = 7),
        plot.title   = element_text(face = "bold"),
        legend.position = "right"
      )

    ggsave(file.path(int_dir, "dotplot_le_overlap_mitocarta.pdf"),
           p_le_mc, width = 12, height = max(6, nrow(le_mc_plot) * 0.3))
  }
}

# --- 5c. Volcano-style plot: interaction genes ---
# Show all unique-to-Myc+ genes with interaction stats, highlight the 357
p_volcano <- ggplot(unique_pos_interaction,
                    aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(data = filter(unique_pos_interaction, pvalue >= 0.05),
             colour = "grey80", size = 1, alpha = 0.5) +
  geom_point(data = filter(unique_pos_interaction, pvalue < 0.05),
             aes(colour = log2FoldChange < 0),
             size = 1.5, alpha = 0.7) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey50") +
  scale_colour_manual(
    values = c("TRUE" = "#2166AC", "FALSE" = "#B2182B"),
    labels = c("TRUE" = "Negative (Myc weakens)", "FALSE" = "Positive (Myc strengthens)"),
    name = "Direction"
  ) +
  labs(
    title = "Interaction statistics: genes unique to Myc+ timepoint effect",
    subtitle = sprintf(
      "%d genes total  |  %d with p < 0.05 (%d negative, %d positive)",
      nrow(unique_pos_interaction), nrow(int_sig_genes), n_neg, n_pos
    ),
    x = "Interaction log2FC (timepoint12W:myc_statuspos)",
    y = expression(-log[10](p))
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "bottom")

ggsave(file.path(int_dir, "volcano_interaction_unique_to_pos.pdf"),
       p_volcano, width = 8, height = 6)

# --- 5d. Overlap with curated gene sets (bar plot) ---
# How many of the 357 genes fall in each of our curated gene set categories?
gene_set_membership <- lapply(names(gene_sets), function(gs_name) {
  gs_genes <- gene_sets[[gs_name]]
  overlap  <- intersect(int_sig_symbols_set, gs_genes)
  if (length(overlap) > 0) {
    data.frame(
      gene_set = gs_name,
      n_overlap = length(overlap),
      pct_of_query = round(100 * length(overlap) / length(int_sig_symbols_set), 1),
      pct_of_set   = round(100 * length(overlap) / length(gs_genes), 1)
    )
  } else {
    NULL
  }
}) |> bind_rows()

# Categorise and summarise
gene_set_membership <- gene_set_membership |>
  mutate(
    category = case_when(
      grepl("^MC_", gene_set)     ~ "MitoCarta",
      grepl("^MYC_", gene_set)    ~ "MYC signature",
      grepl("^MSigDB_", gene_set) ~ "MSigDB Hallmark"
    )
  ) |>
  arrange(desc(n_overlap))

message("\n--- Overlap with curated gene sets (top 20) ---")
gene_set_membership |>
  head(20) |>
  print()

# =============================================================================
# PART 5B: FOCUSED INTERACTION PLOTS — SELECTION vs REMODELLING
# =============================================================================
# If the fading Myc effect on Bbc3 (and other apoptotic genes) reflects clonal
# selection against apoptosis-prone cells, we'd expect a COORDINATED shift
# across pro-apoptotic genes — because selection removes entire cells, not
# individual transcripts. If instead Myc's transcriptional engagement with
# specific promoters changes, we'd expect ISOLATED gene-level changes.
#
# Four panels test this:
#   1. Apoptosis-PRO (25 genes)  — forest plot: is Bbc3 an outlier?
#   2. Apoptosis-ANTI (9 genes)  — forest plot: do anti-apoptotic genes change?
#   3. MitoCarta (all genes)     — volcano: broader mitochondrial landscape
#   4. Felsher integrative MYC   — volcano: is the fading Myc-specific?
#
# Forest plots for small sets: show every gene with LFC ± SE, highlight p<0.05
# Volcano plots for large sets: LFC vs -log10(p), label significant genes

message("\n", strrep("=", 70))
message("PART 5B: FOCUSED INTERACTION PLOTS (SELECTION vs REMODELLING)")
message(strrep("=", 70))

# --- Prepare interaction statistics for all genes (mapped to symbols) ---
# interaction_raw_df already exists from Part 1 (ensembl_id, gene_symbol, LFC, etc.)

# Helper: extract interaction stats for a gene set (by symbol)
get_interaction_stats <- function(symbols, interaction_df) {
  interaction_df |>
    filter(gene_symbol %in% symbols) |>
    # Deduplicate: keep strongest signal per symbol
    group_by(gene_symbol) |>
    slice_min(pvalue, n = 1, with_ties = FALSE) |>
    ungroup() |>
    mutate(
      sig = !is.na(pvalue) & pvalue < 0.05,
      direction = ifelse(log2FoldChange < 0, "Negative", "Positive")
    )
}

# --- Gene sets ---
apoptosis_pro_symbols <- gene_sets[["MC_Apoptosis_Pro"]]
apoptosis_anti_symbols <- gene_sets[["MC_Apoptosis_Anti"]]
felsher_symbols <- gene_sets[["MYC_felsher_integrative_signature"]]

# All MitoCarta genes (union across all pathways from gene_to_pathway)
mitocarta_all_symbols <- unique(gene_to_pathway$Gene)

# Get interaction stats for each set
stats_apop_pro  <- get_interaction_stats(apoptosis_pro_symbols, interaction_raw_df)
stats_apop_anti <- get_interaction_stats(apoptosis_anti_symbols, interaction_raw_df)
stats_mitocarta <- get_interaction_stats(mitocarta_all_symbols, interaction_raw_df)
stats_felsher   <- get_interaction_stats(felsher_symbols, interaction_raw_df)

message(sprintf("  Apoptosis-PRO:  %d / %d genes found, %d with interaction p < 0.05",
                nrow(stats_apop_pro), length(apoptosis_pro_symbols),
                sum(stats_apop_pro$sig, na.rm = TRUE)))
message(sprintf("  Apoptosis-ANTI: %d / %d genes found, %d with interaction p < 0.05",
                nrow(stats_apop_anti), length(apoptosis_anti_symbols),
                sum(stats_apop_anti$sig, na.rm = TRUE)))
message(sprintf("  MitoCarta:      %d / %d genes found, %d with interaction p < 0.05",
                nrow(stats_mitocarta), length(mitocarta_all_symbols),
                sum(stats_mitocarta$sig, na.rm = TRUE)))
message(sprintf("  Felsher MYC:    %d / %d genes found, %d with interaction p < 0.05",
                nrow(stats_felsher), length(felsher_symbols),
                sum(stats_felsher$sig, na.rm = TRUE)))

# --- 5B-1. Forest plot: Apoptosis-PRO ---
p_forest_pro <- ggplot(stats_apop_pro,
                       aes(x = log2FoldChange,
                           y = reorder(gene_symbol, log2FoldChange))) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbar(aes(xmin = log2FoldChange - 1.96 * lfcSE,
                    xmax = log2FoldChange + 1.96 * lfcSE,
                     colour = sig),
                 height = 0.3, linewidth = 0.5, orientation = "y") +
  geom_point(aes(colour = sig, size = sig)) +
  scale_colour_manual(values = c("TRUE" = "#B2182B", "FALSE" = "grey60"),
                      labels = c("TRUE" = "p < 0.05", "FALSE" = "n.s."),
                      name = "Interaction") +
  scale_size_manual(values = c("TRUE" = 3, "FALSE" = 2), guide = "none") +
  labs(
    title = "Apoptosis-PRO: interaction LFC (Myc effect change 6W → 12W)",
    subtitle = paste0(
      "Negative = Myc effect weakens at 12W  |  Error bars = 95% CI\n",
      nrow(stats_apop_pro), " genes  |  ",
      sum(stats_apop_pro$sig, na.rm = TRUE), " significant (p < 0.05)"
    ),
    x = "Interaction log2FC (timepoint12W:myc_statuspos)",
    y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "bottom")

# --- 5B-2. Forest plot: Apoptosis-ANTI ---
p_forest_anti <- ggplot(stats_apop_anti,
                        aes(x = log2FoldChange,
                            y = reorder(gene_symbol, log2FoldChange))) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbar(aes(xmin = log2FoldChange - 1.96 * lfcSE,
                    xmax = log2FoldChange + 1.96 * lfcSE,
                     colour = sig),
                 height = 0.3, linewidth = 0.5, orientation = "y") +
  geom_point(aes(colour = sig, size = sig)) +
  scale_colour_manual(values = c("TRUE" = "#B2182B", "FALSE" = "grey60"),
                      labels = c("TRUE" = "p < 0.05", "FALSE" = "n.s."),
                      name = "Interaction") +
  scale_size_manual(values = c("TRUE" = 3, "FALSE" = 2), guide = "none") +
  labs(
    title = "Apoptosis-ANTI: interaction LFC (Myc effect change 6W → 12W)",
    subtitle = paste0(
      "Negative = Myc effect weakens at 12W  |  Error bars = 95% CI\n",
      nrow(stats_apop_anti), " genes  |  ",
      sum(stats_apop_anti$sig, na.rm = TRUE), " significant (p < 0.05)"
    ),
    x = "Interaction log2FC (timepoint12W:myc_statuspos)",
    y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "bottom")

# --- 5B-3. Volcano plot: MitoCarta ---
# Label genes with p < 0.05 (or top N if too many)
n_sig_mc_vol <- sum(stats_mitocarta$sig, na.rm = TRUE)
label_mc <- if (n_sig_mc_vol <= 40) {
  filter(stats_mitocarta, sig)
} else {
  stats_mitocarta |> arrange(pvalue) |> head(40)
}

p_volcano_mc <- ggplot(stats_mitocarta,
                       aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(data = filter(stats_mitocarta, !sig),
             colour = "grey80", size = 1, alpha = 0.5) +
  geom_point(data = filter(stats_mitocarta, sig),
             aes(colour = direction), size = 2, alpha = 0.7) +
  ggrepel::geom_text_repel(
    data = label_mc,
    aes(label = gene_symbol, colour = direction),
    size = 2.5, max.overlaps = 30,
    segment.colour = "grey60", segment.size = 0.3,
    show.legend = FALSE
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey50") +
  scale_colour_manual(values = c("Negative" = "#2166AC", "Positive" = "#B2182B"),
                      name = "Direction") +
  labs(
    title = "MitoCarta genes: interaction LFC",
    subtitle = sprintf(
      "%d mitochondrial genes  |  %d with p < 0.05 (%d neg, %d pos)",
      nrow(stats_mitocarta), n_sig_mc_vol,
      sum(stats_mitocarta$sig & stats_mitocarta$direction == "Negative", na.rm = TRUE),
      sum(stats_mitocarta$sig & stats_mitocarta$direction == "Positive", na.rm = TRUE)
    ),
    x = "Interaction log2FC",
    y = expression(-log[10](p))
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "bottom")

# --- 5B-4. Volcano plot: Felsher integrative MYC signature ---
n_sig_felsher <- sum(stats_felsher$sig, na.rm = TRUE)
label_felsher <- if (n_sig_felsher <= 40) {
  filter(stats_felsher, sig)
} else {
  stats_felsher |> arrange(pvalue) |> head(40)
}

p_volcano_felsher <- ggplot(stats_felsher,
                            aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(data = filter(stats_felsher, !sig),
             colour = "grey80", size = 1, alpha = 0.5) +
  geom_point(data = filter(stats_felsher, sig),
             aes(colour = direction), size = 2, alpha = 0.7) +
  ggrepel::geom_text_repel(
    data = label_felsher,
    aes(label = gene_symbol, colour = direction),
    size = 2.5, max.overlaps = 30,
    segment.colour = "grey60", segment.size = 0.3,
    show.legend = FALSE
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey50") +
  scale_colour_manual(values = c("Negative" = "#2166AC", "Positive" = "#B2182B"),
                      name = "Direction") +
  labs(
    title = "Felsher integrative MYC signature: interaction LFC",
    subtitle = sprintf(
      "%d MYC signature genes  |  %d with p < 0.05 (%d neg, %d pos)",
      nrow(stats_felsher), n_sig_felsher,
      sum(stats_felsher$sig & stats_felsher$direction == "Negative", na.rm = TRUE),
      sum(stats_felsher$sig & stats_felsher$direction == "Positive", na.rm = TRUE)
    ),
    x = "Interaction log2FC",
    y = expression(-log[10](p))
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "bottom")

# --- Combined PDF: 4 panels ---
pdf(file.path(int_dir, "focused_interaction_selection_vs_remodelling.pdf"),
    width = 14, height = 16)

gridExtra::grid.arrange(
  p_forest_pro, p_forest_anti,
  p_volcano_mc, p_volcano_felsher,
  ncol = 2,
  top = grid::textGrob(
    "Selection vs remodelling: interaction LFC across gene sets",
    gp = grid::gpar(fontsize = 14, fontface = "bold")
  )
)

dev.off()

# Also save individual panels
ggsave(file.path(int_dir, "forest_apoptosis_pro_interaction.pdf"),
       p_forest_pro, width = 8, height = 7)
ggsave(file.path(int_dir, "forest_apoptosis_anti_interaction.pdf"),
       p_forest_anti, width = 8, height = 5)
ggsave(file.path(int_dir, "volcano_mitocarta_interaction.pdf"),
       p_volcano_mc, width = 10, height = 8)
ggsave(file.path(int_dir, "volcano_felsher_myc_interaction.pdf"),
       p_volcano_felsher, width = 10, height = 8)

# --- Summary statistics for the selection vs remodelling argument ---
message("\n--- Selection vs remodelling evidence ---")

# Apoptosis-PRO: is Bbc3 an outlier?
apop_pro_median_lfc <- median(stats_apop_pro$log2FoldChange, na.rm = TRUE)
bbc3_lfc <- stats_apop_pro |> filter(gene_symbol == "Bbc3") |> pull(log2FoldChange)
bbc3_rank <- which(sort(stats_apop_pro$log2FoldChange) == bbc3_lfc)

message(sprintf("  Apoptosis-PRO median interaction LFC: %.4f", apop_pro_median_lfc))
if (length(bbc3_lfc) > 0) {
  message(sprintf("  Bbc3 interaction LFC: %.4f (rank %d / %d in pathway)",
                  bbc3_lfc, bbc3_rank[1], nrow(stats_apop_pro)))
}

# Test: is the Apoptosis-PRO set shifted as a group?
apop_pro_ttest <- t.test(stats_apop_pro$log2FoldChange, mu = 0)
message(sprintf("  Apoptosis-PRO one-sample t-test (mean ≠ 0): p = %.4f, mean = %.4f",
                apop_pro_ttest$p.value, apop_pro_ttest$estimate))

apop_anti_ttest <- t.test(stats_apop_anti$log2FoldChange, mu = 0)
message(sprintf("  Apoptosis-ANTI one-sample t-test (mean ≠ 0): p = %.4f, mean = %.4f",
                apop_anti_ttest$p.value, apop_anti_ttest$estimate))

message("\n  Interpretation:")
if (apop_pro_ttest$p.value > 0.05 && length(bbc3_lfc) > 0 && abs(bbc3_lfc) > abs(apop_pro_median_lfc) * 2) {
  message("  → Apoptosis-PRO is NOT shifted as a group; Bbc3 is an outlier")
  message("  → FAVOURS transcriptional remodelling over clonal selection")
} else if (apop_pro_ttest$p.value < 0.05) {
  message("  → Apoptosis-PRO IS shifted as a group")
  message("  → CONSISTENT WITH clonal selection (or coordinated transcriptional change)")
} else {
  message("  → No clear pattern — examine forest plot for details")
}

# =============================================================================
# PART 5C: PARALLEL ANALYSIS — UNIQUE-TO-MYC− AND OVERLAPPING GENE SETS
# =============================================================================
# The unique-to-Myc+ analysis above characterises genes where the temporal
# change is only significant in Myc+ samples. Two complementary gene sets
# complete the picture:
#
# UNIQUE-TO-MYC−: ~906 genes with significant timepoint effect only in Myc−.
#   These change developmentally but NOT in Myc+ samples. If they show
#   positive interaction (Myc effect strengthens at 12W), it means Myc
#   BUFFERS against a developmental change that occurs in controls.
#   If no interaction signal: the failure to reach significance in Myc+
#   is just a power issue, not a real difference.
#
# OVERLAPPING: ~961 genes significant in both genotypes.
#   These are the core shared developmental programme. Their interaction
#   LFCs should centre on zero (negative control). Any genes with
#   interaction signal here are particularly credible — robustly
#   temporally regulated genes where Myc modulates the MAGNITUDE of a
#   change that occurs in both genotypes.

message("\n", strrep("=", 70))
message("PART 5C: UNIQUE-TO-MYC− AND OVERLAPPING GENE SETS")
message(strrep("=", 70))

# --- Define the three gene sets ---
unique_to_neg <- setdiff(sig_timepoint_neg, sig_timepoint_pos)
overlap_genes <- intersect(sig_timepoint_pos, sig_timepoint_neg)

message(sprintf("  Unique to Myc+:  %d genes", length(unique_to_pos)))
message(sprintf("  Unique to Myc−:  %d genes", length(unique_to_neg)))
message(sprintf("  Overlapping:     %d genes", length(overlap_genes)))

# --- Extract interaction stats for each set ---
unique_neg_interaction <- interaction_raw_df |>
  filter(ensembl_id %in% unique_to_neg)

overlap_interaction <- interaction_raw_df |>
  filter(ensembl_id %in% overlap_genes)

# --- Helper function for characterising each set ---
characterise_set <- function(int_df, set_name) {
  message(sprintf("\n--- %s ---", set_name))

  # Genes with nominal interaction p < 0.05
  sig_int <- int_df |> filter(!is.na(pvalue) & pvalue < 0.05)
  expected <- round(sum(!is.na(int_df$pvalue)) * 0.05)

  message(sprintf("  Total genes: %d", nrow(int_df)))
  message(sprintf("  Interaction p < 0.05: %d (expected by chance: %d, ratio: %.1fx)",
                  nrow(sig_int), expected,
                  nrow(sig_int) / max(expected, 1)))

  # Direction
  n_neg_dir <- sum(sig_int$log2FoldChange < 0, na.rm = TRUE)
  n_pos_dir <- sum(sig_int$log2FoldChange > 0, na.rm = TRUE)
  pct_neg_dir <- if (nrow(sig_int) > 0) round(100 * n_neg_dir / nrow(sig_int), 1) else NA

  message(sprintf("  Direction (p < 0.05): %d negative (%.1f%%), %d positive (%.1f%%)",
                  n_neg_dir, pct_neg_dir, n_pos_dir, 100 - pct_neg_dir))

  if (nrow(sig_int) > 1) {
    binom <- binom.test(n_neg_dir, nrow(sig_int), p = 0.5)
    message(sprintf("  Binomial test: p = %.2e", binom$p.value))
  }

  # Global distribution (all genes in set, not just p < 0.05)
  all_median <- median(int_df$log2FoldChange, na.rm = TRUE)
  all_mean   <- mean(int_df$log2FoldChange, na.rm = TRUE)
  ttest_all  <- tryCatch(
    t.test(int_df$log2FoldChange, mu = 0),
    error = function(e) NULL
  )

  message(sprintf("  All genes — median LFC: %.4f, mean LFC: %.4f", all_median, all_mean))
  if (!is.null(ttest_all)) {
    message(sprintf("  One-sample t-test (mean != 0): p = %.2e", ttest_all$p.value))
  }

  list(
    int_df      = int_df,
    sig_int     = sig_int,
    n_sig       = nrow(sig_int),
    expected    = expected,
    n_neg       = n_neg_dir,
    n_pos       = n_pos_dir,
    pct_neg     = pct_neg_dir,
    median_lfc  = all_median,
    mean_lfc    = all_mean,
    ttest_p     = if (!is.null(ttest_all)) ttest_all$p.value else NA
  )
}

char_unique_neg <- characterise_set(unique_neg_interaction, "Unique to Myc−")
char_overlap    <- characterise_set(overlap_interaction, "Overlapping (both genotypes)")
# For reference, recompute for unique-to-Myc+ with same framework
char_unique_pos <- characterise_set(unique_pos_interaction, "Unique to Myc+ (recap)")

# --- Combined direction histogram ---
all_three_sets <- bind_rows(
  unique_pos_interaction |>
    mutate(gene_set = sprintf("Unique to Myc+ (n=%d)", length(unique_to_pos))),
  unique_neg_interaction |>
    mutate(gene_set = sprintf("Unique to Myc− (n=%d)", length(unique_to_neg))),
  overlap_interaction |>
    mutate(gene_set = sprintf("Overlapping (n=%d)", length(overlap_genes)))
) |>
  filter(!is.na(log2FoldChange))

# Order factor so Myc+ is first
all_three_sets$gene_set <- factor(all_three_sets$gene_set,
                                   levels = unique(all_three_sets$gene_set))

p_three_hist <- ggplot(all_three_sets, aes(x = log2FoldChange)) +
  geom_histogram(bins = 50, fill = "steelblue", colour = "white", boundary = 0) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "red", linewidth = 0.5) +
  facet_wrap(~ gene_set, ncol = 1, scales = "free_y") +
  labs(
    title = "Interaction LFC distribution by gene set",
    subtitle = paste0(
      "Genes with significant timepoint effect (padj < 0.1) in one or both genotypes\n",
      "Negative LFC = Myc effect weakens at 12W  |  Red dashed = zero"
    ),
    x = "Interaction log2FC (timepoint12W:myc_statuspos)",
    y = "Count"
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"))

ggsave(file.path(int_dir, "histogram_interaction_lfc_three_sets.pdf"),
       p_three_hist, width = 8, height = 10)

# --- Volcano plots: unique-to-Myc− and overlapping ---

# Helper to make volcano for any set
make_interaction_volcano <- function(int_df, title, subtitle_extra = "") {
  int_df <- int_df |>
    mutate(
      sig = !is.na(pvalue) & pvalue < 0.05,
      direction = ifelse(log2FoldChange < 0, "Negative", "Positive")
    )

  n_sig <- sum(int_df$sig, na.rm = TRUE)

  label_genes <- if (n_sig <= 40) {
    filter(int_df, sig)
  } else {
    int_df |> filter(sig) |> arrange(pvalue) |> head(40)
  }

  # Map Ensembl to symbols for labels
  label_genes <- label_genes |>
    mutate(label = ifelse(!is.na(gene_symbol) & gene_symbol != "",
                          gene_symbol, ensembl_id))

  p <- ggplot(int_df, aes(x = log2FoldChange, y = -log10(pvalue))) +
    geom_point(data = filter(int_df, !sig),
               colour = "grey80", size = 1, alpha = 0.5) +
    geom_point(data = filter(int_df, sig),
               aes(colour = direction), size = 1.5, alpha = 0.7) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey50") +
    scale_colour_manual(values = c("Negative" = "#2166AC", "Positive" = "#B2182B"),
                        name = "Direction")

  if (nrow(label_genes) > 0) {
    p <- p + ggrepel::geom_text_repel(
      data = label_genes,
      aes(label = label, colour = direction),
      size = 2.5, max.overlaps = 30,
      segment.colour = "grey60", segment.size = 0.3,
      show.legend = FALSE
    )
  }

  n_neg_sig <- sum(int_df$sig & int_df$direction == "Negative", na.rm = TRUE)
  n_pos_sig <- sum(int_df$sig & int_df$direction == "Positive", na.rm = TRUE)

  p + labs(
    title = title,
    subtitle = sprintf(
      "%d genes  |  %d with p < 0.05 (%d neg, %d pos)%s",
      nrow(int_df), n_sig, n_neg_sig, n_pos_sig, subtitle_extra
    ),
    x = "Interaction log2FC (timepoint12W:myc_statuspos)",
    y = expression(-log[10](p))
  ) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "bottom")
}

p_volcano_neg <- make_interaction_volcano(
  unique_neg_interaction,
  "Interaction: genes unique to Myc- timepoint effect",
  sprintf("\nExpected by chance: %d", char_unique_neg$expected)
)

p_volcano_overlap <- make_interaction_volcano(
  overlap_interaction,
  "Interaction: genes with timepoint effect in both genotypes",
  sprintf("\nExpected by chance: %d", char_overlap$expected)
)

ggsave(file.path(int_dir, "volcano_interaction_unique_to_neg.pdf"),
       p_volcano_neg, width = 10, height = 8)
ggsave(file.path(int_dir, "volcano_interaction_overlapping.pdf"),
       p_volcano_overlap, width = 10, height = 8)

# --- Combined 3-panel volcano ---
pdf(file.path(int_dir, "volcano_interaction_three_sets.pdf"),
    width = 10, height = 20)
gridExtra::grid.arrange(
  # Reuse the existing volcano from Part 5c for unique-to-Myc+ (p_volcano)
  p_volcano +
    labs(title = "Interaction: genes unique to Myc+ timepoint effect"),
  p_volcano_neg,
  p_volcano_overlap,
  ncol = 1,
  top = grid::textGrob(
    "Interaction signal across three temporal gene sets",
    gp = grid::gpar(fontsize = 14, fontface = "bold")
  )
)
dev.off()

# --- g:Profiler on unique-to-Myc− interaction genes (if enough) ---
neg_sig_int <- unique_neg_interaction |>
  filter(!is.na(pvalue) & pvalue < 0.05)

neg_sig_symbols <- neg_sig_int |>
  filter(!is.na(gene_symbol) & gene_symbol != "") |>
  pull(gene_symbol) |>
  unique()

gost_neg_set <- NULL
if (length(neg_sig_symbols) >= 10) {
  message(sprintf("\n  Running g:Profiler on %d unique-to-Myc− interaction genes...",
                  length(neg_sig_symbols)))

  gost_neg_set <- gost(
    query         = neg_sig_symbols,
    organism      = "mmusculus",
    ordered_query = FALSE,
    significant   = FALSE,
    custom_bg     = all_tested_symbols,
    sources       = c("GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC", "WP"),
    evcodes       = TRUE
  )

  if (!is.null(gost_neg_set) && nrow(gost_neg_set$result) > 0) {
    n_sig_gost_neg <- sum(gost_neg_set$result$p_value < 0.05)
    message(sprintf("  g:Profiler (unique-to-Myc−): %d significant terms", n_sig_gost_neg))

    if (n_sig_gost_neg > 0) {
      message("  Top 15:")
      gost_neg_set$result |>
        filter(p_value < 0.05) |>
        arrange(p_value) |>
        dplyr::select(source, term_name, term_size, intersection_size, p_value) |>
        head(15) |>
        print()
    }
  } else {
    message("  No significant terms")
  }
} else {
  message(sprintf("\n  Unique-to-Myc− interaction genes: only %d — skipping g:Profiler",
                  length(neg_sig_symbols)))
}

# --- g:Profiler on overlapping set interaction genes (if enough) ---
overlap_sig_int <- overlap_interaction |>
  filter(!is.na(pvalue) & pvalue < 0.05)

overlap_sig_symbols <- overlap_sig_int |>
  filter(!is.na(gene_symbol) & gene_symbol != "") |>
  pull(gene_symbol) |>
  unique()

gost_overlap_set <- NULL
if (length(overlap_sig_symbols) >= 10) {
  message(sprintf("\n  Running g:Profiler on %d overlapping-set interaction genes...",
                  length(overlap_sig_symbols)))

  gost_overlap_set <- gost(
    query         = overlap_sig_symbols,
    organism      = "mmusculus",
    ordered_query = FALSE,
    significant   = FALSE,
    custom_bg     = all_tested_symbols,
    sources       = c("GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC", "WP"),
    evcodes       = TRUE
  )

  if (!is.null(gost_overlap_set) && nrow(gost_overlap_set$result) > 0) {
    n_sig_gost_ov <- sum(gost_overlap_set$result$p_value < 0.05)
    message(sprintf("  g:Profiler (overlapping): %d significant terms", n_sig_gost_ov))

    if (n_sig_gost_ov > 0) {
      message("  Top 15:")
      gost_overlap_set$result |>
        filter(p_value < 0.05) |>
        arrange(p_value) |>
        dplyr::select(source, term_name, term_size, intersection_size, p_value) |>
        head(15) |>
        print()
    }
  } else {
    message("  No significant terms")
  }
} else {
  message(sprintf("\n  Overlapping interaction genes: only %d — skipping g:Profiler",
                  length(overlap_sig_symbols)))
}

# --- Comparative summary table ---
comparison_table <- data.frame(
  gene_set = c("Unique to Myc+", "Unique to Myc-", "Overlapping"),
  total_genes = c(length(unique_to_pos), length(unique_to_neg), length(overlap_genes)),
  int_p05 = c(nrow(int_sig_genes), char_unique_neg$n_sig, char_overlap$n_sig),
  expected = c(round(length(unique_to_pos) * 0.05),
               char_unique_neg$expected,
               char_overlap$expected),
  ratio = c(nrow(int_sig_genes) / max(round(length(unique_to_pos) * 0.05), 1),
            char_unique_neg$n_sig / max(char_unique_neg$expected, 1),
            char_overlap$n_sig / max(char_overlap$expected, 1)),
  pct_negative = c(pct_neg, char_unique_neg$pct_neg, char_overlap$pct_neg),
  median_lfc = c(median(int_sig_genes$log2FoldChange),
                 char_unique_neg$median_lfc,
                 char_overlap$median_lfc),
  ttest_p = c(t.test(unique_pos_interaction$log2FoldChange, mu = 0)$p.value,
              char_unique_neg$ttest_p,
              char_overlap$ttest_p)
)

message("\n--- Comparative summary: interaction signal across gene sets ---")
print(comparison_table)

write.csv(comparison_table,
          file.path(int_dir, "interaction_three_sets_summary.csv"),
          row.names = FALSE)

# =============================================================================
# PART 6: SUMMARY
# =============================================================================

message("\n", strrep("=", 70))
message("SUMMARY: INTERACTION GENE SET CHARACTERISATION")
message(strrep("=", 70))

message(sprintf("\n1. GENE SET SIZE: %d genes (from %d unique-to-Myc+ with interaction p < 0.05)",
                nrow(int_sig_genes), length(unique_to_pos)))

message(sprintf("\n2. DIRECTION: %d negative (%.1f%%), %d positive (%.1f%%)",
                n_neg, pct_neg, n_pos, 100 - pct_neg))
message(sprintf("   Binomial test: p = %.2e", binom_result$p.value))
if (pct_neg > 60) {
  message("   → Predominantly negative: Myc effect WEAKENS at 12W for these genes")
} else if (pct_neg < 40) {
  message("   → Predominantly positive: Myc effect STRENGTHENS at 12W for these genes")
} else {
  message("   → Mixed directions: interaction signal is heterogeneous")
}

if (!is.null(gost_ordered) && nrow(filter(gost_ordered$result, p_value < 0.05)) > 0) {
  n_ordered <- sum(gost_ordered$result$p_value < 0.05)
  n_unordered <- if (!is.null(gost_all)) sum(gost_all$result$p_value < 0.05) else 0
  top_source <- gost_ordered$result |>
    filter(p_value < 0.05) |>
    dplyr::count(source) |>
    arrange(desc(n)) |>
    head(3)
  message(sprintf("\n3. FUNCTIONAL ENRICHMENT:"))
  message(sprintf("   Ordered query (ranked by |Wald|): %d significant terms", n_ordered))
  message(sprintf("   Unordered ORA (p < 0.05 cutoff):  %d significant terms", n_unordered))
  if (exists("shared")) {
    message(sprintf("   Convergent (both): %d  |  Ordered only: %d  |  Unordered only: %d",
                    length(shared), length(only_ord), length(only_uno)))
  }
  message("   Top sources: ", paste(sprintf("%s (%d)", top_source$source, top_source$n),
                                     collapse = ", "))
} else if (!is.null(gost_all) && nrow(filter(gost_all$result, p_value < 0.05)) > 0) {
  top_source <- gost_all$result |>
    filter(p_value < 0.05) |>
    dplyr::count(source) |>
    arrange(desc(n)) |>
    head(3)
  message(sprintf("\n3. FUNCTIONAL ENRICHMENT: %d significant g:Profiler terms (unordered)",
                  sum(gost_all$result$p_value < 0.05)))
  message("   Top sources: ", paste(sprintf("%s (%d)", top_source$source, top_source$n),
                                     collapse = ", "))
} else {
  message("\n3. FUNCTIONAL ENRICHMENT: No significant terms from g:Profiler")
}

message(sprintf("\n4. LEADING EDGE OVERLAP: %d / %d (%.1f%%) interaction genes appear in",
                n_in_any_le, length(int_sig_symbols_set), pct_in_le))
message("   at least one significant fGSEA leading edge")
if (!is.null(le_overlap_curated)) {
  n_enriched <- sum(le_overlap_curated$hyper_padj < 0.05, na.rm = TRUE)
  message(sprintf("   %d curated pathways show significant LE enrichment for these genes",
                  n_enriched))
}

message(sprintf("\n5. SELECTION vs REMODELLING:"))
message(sprintf("   Apoptosis-PRO: %d / %d genes with interaction p < 0.05",
                sum(stats_apop_pro$sig, na.rm = TRUE), nrow(stats_apop_pro)))
message(sprintf("   Apoptosis-PRO group shift t-test: p = %.4f (mean LFC = %.4f)",
                apop_pro_ttest$p.value, apop_pro_ttest$estimate))
if (length(bbc3_lfc) > 0) {
  message(sprintf("   Bbc3 LFC = %.4f vs pathway median = %.4f",
                  bbc3_lfc, apop_pro_median_lfc))
}
message(sprintf("   MitoCarta: %d / %d with p < 0.05  |  Felsher MYC: %d / %d with p < 0.05",
                sum(stats_mitocarta$sig, na.rm = TRUE), nrow(stats_mitocarta),
                sum(stats_felsher$sig, na.rm = TRUE), nrow(stats_felsher)))

message("\n6. THREE-SET COMPARISON:")
message("   Gene set          | Total | Int p<0.05 | Expected | Ratio | %Neg  | Median LFC")
for (i in seq_len(nrow(comparison_table))) {
  message(sprintf("   %-18s | %5d | %10d | %8d | %5.1fx | %5.1f%% | %.4f",
                  comparison_table$gene_set[i],
                  comparison_table$total_genes[i],
                  comparison_table$int_p05[i],
                  comparison_table$expected[i],
                  comparison_table$ratio[i],
                  comparison_table$pct_negative[i],
                  comparison_table$median_lfc[i]))
}

# =============================================================================
# PART 7: SAVE RESULTS
# =============================================================================

message("\nSaving results...")

interaction_gene_results <- list(
  # The gene set
  int_sig_genes      = int_sig_genes,
  int_sig_symbols    = int_sig_symbols_set,
  unique_to_pos      = unique_to_pos,

  # Direction analysis
  direction = list(
    n_negative  = n_neg,
    n_positive  = n_pos,
    pct_negative = pct_neg,
    binom_p     = binom_result$p.value,
    median_lfc  = median(int_sig_genes$log2FoldChange),
    mean_lfc    = mean(int_sig_genes$log2FoldChange)
  ),

  # g:Profiler results
  gprofiler_ordered = gost_ordered,
  gprofiler_all = gost_all,
  gprofiler_neg = gost_neg,
  gprofiler_pos = gost_pos,
  gprofiler_comparison = if (exists("shared")) {
    list(shared = shared, ordered_only = only_ord, unordered_only = only_uno)
  } else { NULL },

  # Leading edge overlap
  le_overlap_curated = le_overlap_curated,
  le_overlap_mc      = le_overlap_mc,
  le_membership = list(
    n_in_any_le   = n_in_any_le,
    pct_in_le     = pct_in_le,
    not_in_le     = not_in_le
  ),

  # Focused interaction stats (selection vs remodelling)
  interaction_by_geneset = list(
    apoptosis_pro  = stats_apop_pro,
    apoptosis_anti = stats_apop_anti,
    mitocarta      = stats_mitocarta,
    felsher_myc    = stats_felsher
  ),

  # Three-set comparison (unique-to-Myc+, unique-to-Myc−, overlapping)
  three_set_analysis = list(
    unique_to_pos  = char_unique_pos,
    unique_to_neg  = char_unique_neg,
    overlapping    = char_overlap,
    comparison_table = comparison_table
  ),
  gprofiler_unique_neg = gost_neg_set,
  gprofiler_overlap    = gost_overlap_set,

  # Gene set membership
  gene_set_membership = gene_set_membership,

  # Metadata
  analysis_date = Sys.Date(),
  description = paste(
    "Characterisation of ~357 genes unique to Myc+ timepoint effect with",
    "nominal interaction p < 0.05. Direction analysis, g:Profiler ORA,",
    "and leading edge overlap with fGSEA interaction pathways from script 10."
  )
)

saveRDS(interaction_gene_results,
        here("results", "interaction_gene_characterisation.rds"))

# CSV exports
int_sig_genes |>
  arrange(pvalue) |>
  write.csv(file.path(int_dir, "interaction_sig_genes_p05.csv"),
            row.names = FALSE)

if (!is.null(gost_ordered) && nrow(gost_ordered$result) > 0) {
  gost_ordered$result |>
    filter(p_value < 0.05) |>
    arrange(p_value) |>
    dplyr::select(source, term_id, term_name, term_size,
                  query_size, intersection_size, p_value) |>
    write.csv(file.path(int_dir, "gprofiler_ordered_sig.csv"),
              row.names = FALSE)
}

if (!is.null(gost_all) && nrow(gost_all$result) > 0) {
  gost_all$result |>
    filter(p_value < 0.05) |>
    arrange(p_value) |>
    dplyr::select(source, term_id, term_name, term_size,
                  query_size, intersection_size, p_value) |>
    write.csv(file.path(int_dir, "gprofiler_unordered_sig.csv"),
              row.names = FALSE)
}

if (!is.null(le_overlap_curated)) {
  le_overlap_curated |>
    arrange(hyper_p) |>
    write.csv(file.path(int_dir, "le_overlap_curated.csv"), row.names = FALSE)
}

if (!is.null(le_overlap_mc)) {
  le_overlap_mc |>
    arrange(hyper_p) |>
    write.csv(file.path(int_dir, "le_overlap_mitocarta.csv"), row.names = FALSE)
}

message("\n", strrep("=", 70))
message("INTERACTION GENE CHARACTERISATION COMPLETE")
message(strrep("=", 70))
message(sprintf("Results saved to: results/interaction_gene_characterisation.rds"))
message(sprintf("Figures and tables saved to: %s", int_dir))
message(strrep("=", 70))
