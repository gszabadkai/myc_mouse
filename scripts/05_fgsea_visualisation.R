# scripts/05_fgsea_visualisation.R
# =============================================================================
# fGSEA Visualisation
# =============================================================================
#
# Visualisations:
#   1. Dot plots - NES vs pathway, sized by gene set size, coloured by padj
#   2. Comparative bar plots - NES_pos vs NES_neg side by side
#   3. Enrichment plots - running score for key pathways
#
# =============================================================================

source(here::here("scripts", "00_setup_packages.R"))

# === Create output directory ===
fig_dir <- here("outputs", "fgsea")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# === Load results ===
fgsea_results <- readRDS(here("results", "fgsea_results.rds"))

# Extract components
fgsea_combined <- fgsea_results$combined
myc_pathways <- fgsea_results$pathways_myc
mito_pathways <- fgsea_results$pathways_mito
apoptosis_pathways <- fgsea_results$pathways_apoptosis
hallmark_pathways <- fgsea_results$pathways_hallmark
ranks_myc_pos <- fgsea_results$ranks_myc_pos
ranks_myc_neg <- fgsea_results$ranks_myc_neg
gene_sets <- readRDS(here("results", "gene_sets_list.rds"))

# =============================================================================
# PART 1: DOT PLOTS
# =============================================================================

#' Create a dot plot for fGSEA results
create_dotplot <- function(data, title, nes_col = "NES_pos", padj_col = "padj_pos") {
  data |>
    as_tibble() |>
    mutate(
      pathway = str_remove(pathway, "^(MSigDB_HALLMARK_|MYC_|MC_)"),
      pathway = str_replace_all(pathway, "_", " "),
      pathway = forcats::fct_reorder(pathway, .data[[nes_col]])
    ) |>
    ggplot(aes(x = .data[[nes_col]], y = pathway)) +
    geom_point(aes(size = size, colour = .data[[padj_col]])) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    scale_colour_gradient(low = "red", high = "grey80", limits = c(0, 0.1)) +
    scale_size_continuous(range = c(2, 6)) +
    labs(
      title = title,
      x = "Normalized Enrichment Score (NES)",
      y = NULL,
      colour = "padj",
      size = "Gene set size"
    ) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 9),
      plot.title = element_text(size = 11, face = "bold")
    )
}

# MYC pathways dot plot
p_myc_dot <- create_dotplot(
  filter(myc_pathways, padj_pos < 0.1),
  "MYC Signatures: Myc+ tumours (12W vs 6W)"
)

# MitoCarta pathways dot plot
p_mito_dot <- create_dotplot(
  filter(mito_pathways, padj_pos < 0.1),
  "MitoCarta Pathways: Myc+ tumours (12W vs 6W)"
)

# Hallmark pathways dot plot (top 20 by padj)
p_hallmark_dot <- create_dotplot(
  hallmark_pathways |> slice_min(padj_pos, n = 20),
  "Top Hallmark Pathways: Myc+ tumours (12W vs 6W)"
)

# Save dot plots
ggsave(
  file.path(fig_dir, "dotplot_myc_signatures.pdf"),
  p_myc_dot, width = 10, height = 6
)
ggsave(
  file.path(fig_dir, "dotplot_mitocarta.pdf"),
  p_mito_dot, width = 10, height = 7
)
ggsave(
  file.path(fig_dir, "dotplot_hallmark.pdf"),
  p_hallmark_dot, width = 10, height = 8
)

message("Dot plots saved")

# =============================================================================
# PART 2: COMPARATIVE BAR PLOTS (Myc+ vs Myc-)
# =============================================================================

#' Create comparative bar plot showing NES for both comparisons
create_comparative_barplot <- function(data, title) {
  data_long <- data |>
    as_tibble() |>
    mutate(
      pathway = str_remove(pathway, "^(MSigDB_HALLMARK_|MYC_|MC_)"),
      pathway = str_replace_all(pathway, "_", " ")
    ) |>
    dplyr::select(pathway, NES_pos, NES_neg, category) |>
    pivot_longer(
      cols = c(NES_pos, NES_neg),
      names_to = "comparison",
      values_to = "NES"
    ) |>
    mutate(
      comparison = ifelse(comparison == "NES_pos", "Myc+ (12W vs 6W)", "Myc- (12W vs 6W)"),
      pathway = forcats::fct_reorder(pathway, NES, .fun = function(x) x[1])
    )
  
  ggplot(data_long, aes(x = NES, y = pathway, fill = comparison)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    scale_fill_manual(values = c("Myc+ (12W vs 6W)" = "#E41A1C", "Myc- (12W vs 6W)" = "#377EB8")) +
    labs(
      title = title,
      x = "Normalized Enrichment Score (NES)",
      y = NULL,
      fill = "Comparison"
    ) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 9),
      plot.title = element_text(size = 11, face = "bold"),
      legend.position = "bottom"
    )
}

# MYC signatures comparison
p_myc_bar <- create_comparative_barplot(
  filter(myc_pathways, padj_pos < 0.1 | padj_neg < 0.1),
  "MYC Signatures: Myc+ vs Myc- comparison"
)

# MitoCarta comparison
p_mito_bar <- create_comparative_barplot(
  filter(mito_pathways, padj_pos < 0.1 | padj_neg < 0.1),
  "MitoCarta Pathways: Myc+ vs Myc- comparison"
)

# Hallmark comparison (significant in either)
p_hallmark_bar <- create_comparative_barplot(
  hallmark_pathways |> filter(padj_pos < 0.05 | padj_neg < 0.05),
  "Hallmark Pathways: Myc+ vs Myc- comparison"
)

# Save comparative bar plots
ggsave(
  file.path(fig_dir, "barplot_myc_comparison.pdf"),
  p_myc_bar, width = 10, height = 6
)
ggsave(
  file.path(fig_dir, "barplot_mito_comparison.pdf"),
  p_mito_bar, width = 10, height = 8
)
ggsave(
  file.path(fig_dir, "barplot_hallmark_comparison.pdf"),
  p_hallmark_bar, width = 10, height = 10
)

message("Comparative bar plots saved")

# =============================================================================
# PART 3: CATEGORY SUMMARY
# =============================================================================

# Summary by category
category_summary <- fgsea_combined |>
  as_tibble() |>
  dplyr::count(category) |>
  mutate(category = forcats::fct_reorder(category, n))

p_category <- ggplot(category_summary, aes(x = n, y = category, fill = category)) +
  geom_col() +
  geom_text(aes(label = n), hjust = -0.2) +
  scale_fill_manual(values = c(
    "Myc+ specific" = "#E41A1C",
    "Developmental only" = "#377EB8",
    "Shared (same direction)" = "#984EA3",
    "Opposite effects" = "#FF7F00",
    "Not significant" = "grey70"
  )) +
  labs(
    title = "Pathway classification summary",
    subtitle = "Based on significance (padj < 0.05) in each comparison",
    x = "Number of pathways",
    y = NULL
  ) +
  theme_minimal() +
  guides(fill = "none") +
  xlim(0, max(category_summary$n) * 1.15)

ggsave(
  file.path(fig_dir, "category_summary.pdf"),
  p_category, width = 8, height = 5
)

message("Category summary saved")

# =============================================================================
# PART 4: ENRICHMENT PLOTS (Running score)
# =============================================================================

# Key pathways to plot
key_pathways <- c(
  "MSigDB_HALLMARK_MYC_TARGETS_V1",
  "MSigDB_HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "MSigDB_HALLMARK_APOPTOSIS",
  "MC_OXPHOS_subunits",
  "MYC_felsher_integrative_signature"
)

# Create enrichment plots for key pathways
for (pw in key_pathways) {
  if (pw %in% names(gene_sets)) {
    # Myc+ enrichment
    p_pos <- plotEnrichment(gene_sets[[pw]], ranks_myc_pos) +
      labs(title = paste0(pw, "\nMyc+ (12W vs 6W)")) +
      theme_minimal()
    
    # Myc- enrichment
    p_neg <- plotEnrichment(gene_sets[[pw]], ranks_myc_neg) +
      labs(title = paste0(pw, "\nMyc- (12W vs 6W)")) +
      theme_minimal()
    
    # Combine
    p_combined <- gridExtra::grid.arrange(p_pos, p_neg, ncol = 2)
    
    ggsave(
      file.path(fig_dir, paste0("enrichment_", gsub("MSigDB_HALLMARK_|MYC_|MC_", "", pw), ".pdf")),
      p_combined, width = 12, height = 5
    )
  }
}

message("Enrichment plots saved")

# =============================================================================
# PART 5: HEATMAP OF TOP PATHWAYS BY CATEGORY
# =============================================================================

# Get top pathways from each category
top_pathways <- fgsea_combined |>
  as_tibble() |>
  filter(category != "Not significant") |>
  group_by(category) |>
  slice_min(padj_pos, n = 5, with_ties = FALSE) |>
  ungroup() |>
  mutate(
    pathway = str_remove(pathway, "^(MSigDB_HALLMARK_|MYC_|MC_)"),
    pathway = str_replace_all(pathway, "_", " "),
    pathway = forcats::fct_reorder(pathway, NES_diff)
  )

p_heatmap <- top_pathways |>
  dplyr::select(pathway, category, NES_pos, NES_neg) |>
  pivot_longer(cols = c(NES_pos, NES_neg), names_to = "comparison", values_to = "NES") |>
  mutate(comparison = ifelse(comparison == "NES_pos", "Myc+", "Myc-")) |>
  ggplot(aes(x = comparison, y = pathway, fill = NES)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
  facet_grid(category ~ ., scales = "free_y", space = "free_y") +
  labs(
    title = "NES comparison: Myc+ vs Myc- (12W vs 6W)",
    subtitle = "Top 5 pathways per category",
    x = NULL,
    y = NULL,
    fill = "NES"
  ) +
  theme_minimal() +
  theme(
    strip.text.y = element_text(angle = 0, hjust = 0),
    axis.text.y = element_text(size = 8),
    panel.grid = element_blank()
  )

ggsave(
  file.path(fig_dir, "heatmap_top_pathways.pdf"),
  p_heatmap, width = 8, height = 12
)

message("Heatmap saved")

message("\n", strrep("=", 70))
message("fGSEA VISUALISATION COMPLETE")
message("Figures saved to: ", fig_dir)
message(strrep("=", 70))
