# scripts/07_fgsea_xs_visualisation.R
# =============================================================================
# Cross-Sectional fGSEA Visualisation
# =============================================================================
#
# Simple category summary for Myc+ vs Myc- at 6W and 12W timepoints
#
# =============================================================================

source(here::here("scripts", "00_setup_packages.R"))

# === Create output directory ===
fig_dir <- here("outputs", "fgsea_cross_sectional")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# === Load results ===
fgsea_xs_results <- readRDS(here("results", "fgsea_xs_results.rds"))
fgsea_xs_combined <- fgsea_xs_results$combined

# Extract pathway subsets
xs_myc_pathways <- fgsea_xs_results$pathways_myc
xs_mito_pathways <- fgsea_xs_results$pathways_mito
xs_hallmark_pathways <- fgsea_xs_results$pathways_hallmark

# =============================================================================
# CATEGORY SUMMARY
# =============================================================================

category_summary <- fgsea_xs_combined |>
  as_tibble() |>
  dplyr::count(category) |>
  mutate(category = forcats::fct_reorder(category, n))

p_category <- ggplot(category_summary, aes(x = n, y = category, fill = category)) +
  geom_col() +
  geom_text(aes(label = n), hjust = -0.2) +
  scale_fill_manual(values = c(
    "Stable (both timepoints)" = "#4DAF4A",
    "6W only (lost by 12W)" = "#377EB8",
    "12W only (gained by 12W)" = "#E41A1C",
    "Reversed" = "#FF7F00",
    "Not significant" = "grey70"
  )) +
  labs(
    title = "Cross-sectional pathway classification",
    subtitle = "Myc+ vs Myc- enrichment at 6W and 12W (padj < 0.05)",
    x = "Number of pathways",
    y = NULL
  ) +
  theme_minimal() +
  guides(fill = "none") +
  xlim(0, max(category_summary$n) * 1.15)

ggsave(
  file.path(fig_dir, "xs_category_summary.pdf"),
  p_category, width = 8, height = 5
)

message("Category summary saved")

# =============================================================================
# DOT PLOTS
# =============================================================================

#' Create a dot plot for cross-sectional fGSEA results
create_dotplot <- function(data, title, nes_col = "NES_6W", padj_col = "padj_6W") {
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

# MYC pathways dot plots (6W and 12W)
p_myc_6W <- create_dotplot(
 dplyr::filter(xs_myc_pathways, padj_6W < 0.1),
  "MYC Signatures: Myc+ vs Myc- at 6W",
  nes_col = "NES_6W", padj_col = "padj_6W"
)

p_myc_12W <- create_dotplot(
 dplyr::filter(xs_myc_pathways, padj_12W < 0.1),
  "MYC Signatures: Myc+ vs Myc- at 12W",
  nes_col = "NES_12W", padj_col = "padj_12W"
)

# MitoCarta pathways dot plots
p_mito_6W <- create_dotplot(
 dplyr::filter(xs_mito_pathways, padj_6W < 0.1),
  "MitoCarta Pathways: Myc+ vs Myc- at 6W",
  nes_col = "NES_6W", padj_col = "padj_6W"
)

p_mito_12W <- create_dotplot(
 dplyr::filter(xs_mito_pathways, padj_12W < 0.1),
  "MitoCarta Pathways: Myc+ vs Myc- at 12W",
  nes_col = "NES_12W", padj_col = "padj_12W"
)

# Hallmark pathways dot plots (top 20 by padj)
p_hallmark_6W <- create_dotplot(
  xs_hallmark_pathways |> dplyr::slice_min(padj_6W, n = 20),
  "Top Hallmark Pathways: Myc+ vs Myc- at 6W",
  nes_col = "NES_6W", padj_col = "padj_6W"
)

p_hallmark_12W <- create_dotplot(
  xs_hallmark_pathways |> dplyr::slice_min(padj_12W, n = 20),
  "Top Hallmark Pathways: Myc+ vs Myc- at 12W",
  nes_col = "NES_12W", padj_col = "padj_12W"
)

# Save dot plots
ggsave(
  file.path(fig_dir, "dotplot_myc_6W.pdf"),
  p_myc_6W, width = 10, height = 6
)
ggsave(
  file.path(fig_dir, "dotplot_myc_12W.pdf"),
  p_myc_12W, width = 10, height = 6
)
ggsave(
  file.path(fig_dir, "dotplot_mito_6W.pdf"),
  p_mito_6W, width = 10, height = 7
)
ggsave(
  file.path(fig_dir, "dotplot_mito_12W.pdf"),
  p_mito_12W, width = 10, height = 7
)
ggsave(
  file.path(fig_dir, "dotplot_hallmark_6W.pdf"),
  p_hallmark_6W, width = 10, height = 8
)
ggsave(
  file.path(fig_dir, "dotplot_hallmark_12W.pdf"),
  p_hallmark_12W, width = 10, height = 8
)

message("Dot plots saved")

# =============================================================================
# NES CORRELATION PLOT (Temporal vs Cross-Sectional)
# =============================================================================

# Load temporal results for comparison
fgsea_results <- readRDS(here("results", "fgsea_results.rds"))
fgsea_combined_temporal <- fgsea_results$combined

# Combine temporal and cross-sectional results
fgsea_correlation_data <- fgsea_combined_temporal |>
  as_tibble() |>
  dplyr::select(pathway, NES_pos, NES_neg, padj_pos, padj_neg, size) |>
  left_join(
    fgsea_xs_combined |>
      as_tibble() |>
      dplyr::select(pathway, NES_6W, NES_12W, padj_6W, padj_12W),
    by = "pathway"
  ) |>
  mutate(
    # X: developmental baseline effect (Myc- 12W vs 6W)
    x_dev_baseline = NES_neg,

    # Y: maintenance of the Myc effect (change in Myc+ vs Myc- enrichment over time)
    y_myc_maintenance = NES_12W - NES_6W,

    # Size: highest significance across all comparisons
    min_padj = pmin(padj_pos, padj_neg, padj_6W, padj_12W, na.rm = TRUE),
    sig_level = -log10(min_padj),

    # Color: 1 / (NES_pos - NES_neg) = inverse of Myc-specific temporal effect
    nes_diff_temporal = NES_pos - NES_neg,
    color_dev_contribution = 1 / nes_diff_temporal,

    # Clean pathway names for labeling
    pathway_label = str_remove(pathway, "^(MSigDB_HALLMARK_|MYC_|MC_)") |>
      str_replace_all("_", " ")
  ) |>
  mutate(
    color_dev_contribution = case_when(
      is.infinite(color_dev_contribution) ~ NA_real_,
      abs(color_dev_contribution) > 10 ~ sign(color_dev_contribution) * 10,
      TRUE ~ color_dev_contribution
    )
  )

p_correlation <- ggplot(
  fgsea_correlation_data,
  aes(x = x_dev_baseline, y = y_myc_maintenance)
) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_point(
    aes(size = sig_level, colour = color_dev_contribution),
    alpha = 0.7
  ) +
  ggrepel::geom_text_repel(
    data = fgsea_correlation_data |>
      filter(sig_level > 5 & (abs(y_myc_maintenance) > 1 | abs(x_dev_baseline) > 1.5)),
    aes(label = pathway_label),
    size = 2.5,
    max.overlaps = 15,
    segment.alpha = 0.5
  ) +
  scale_colour_gradient2(
    low = "#2166AC", mid = "grey90", high = "#B2182B",
    midpoint = 0,
    limits = c(-5, 5),
    oob = scales::squish,
    na.value = "grey50"
  ) +
  scale_size_continuous(range = c(1, 8), breaks = c(2, 5, 10, 20)) +
  labs(
    title = "fGSEA NES correlations",
    x = "Developmental baseline effect\n(NES: Myc- 12W vs 6W)",
    y = "Maintenance of the Myc effect\n(NES Myc+vsMyc- 12W - NES Myc+vsMyc- 6W)",
    colour = "Developmental\ncontribution\n(1/dNES)",
    size = "-log10(padj)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    legend.position = "right"
  )

ggsave(
  file.path(fig_dir, "xs_nes_correlation.pdf"),
  p_correlation, width = 10, height = 8
)

message("NES correlation plot saved")

# --- Plot 2: Developmental contribution (clearer interpretation) ---
# X: Developmental baseline (NES Myc- 12W vs 6W)
# Y: Developmental contribution (1/dNES) - high absolute values = development explains Myc effect
# Colour: Myc+ trajectory (NES Myc+ 12W vs 6W)
# Size: Significance

p_dev_contribution <- ggplot(
  fgsea_correlation_data,
  aes(x = x_dev_baseline, y = color_dev_contribution)
) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_point(
    aes(size = sig_level, colour = NES_pos),
    alpha = 0.7
  ) +
  ggrepel::geom_text_repel(
    data = fgsea_correlation_data |>
      filter(sig_level > 5 & (abs(color_dev_contribution) > 1 | abs(x_dev_baseline) > 1.5)),
    aes(label = pathway_label),
    size = 2.5,
    max.overlaps = 20,
    segment.alpha = 0.5
  ) +
  scale_colour_gradient2(
    low = "#2166AC", mid = "grey90", high = "#B2182B",
    midpoint = 0,
    limits = c(-3, 3),
    oob = scales::squish
  ) +
  scale_size_continuous(range = c(1, 8), breaks = c(2, 5, 10, 20)) +
  labs(
    title = "fGSEA: Developmental contribution to Myc enrichment",
    x = "Developmental baseline effect\n(NES: Myc- 12W vs 6W)",
    y = "Developmental contribution\n(1 / [NES Myc+ - NES Myc-])",
    colour = "Myc+ trajectory\n(NES 12W vs 6W)",
    size = "-log10(padj)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    legend.position = "right"
  )

ggsave(
  file.path(fig_dir, "xs_dev_contribution.pdf"),
  p_dev_contribution, width = 12, height = 9
)

message("Developmental contribution plot saved")

# =============================================================================
# SIMPLIFIED PLOTS
# =============================================================================

# --- Plot 3: NES timepoint_neg vs NES timepoint_pos ---
# X: NES (Myc- 12W vs 6W) - developmental baseline
# Y: NES (Myc+ 12W vs 6W) - Myc+ trajectory
# Size: significance, Colour: pathway category

fgsea_simple_data <- fgsea_combined_temporal |>
  as_tibble() |>
  mutate(
    min_padj = pmin(padj_pos, padj_neg, na.rm = TRUE),
    sig_level = -log10(min_padj),
    pathway_category = case_when(
      str_detect(pathway, "^MC_|^Mitochondrial|^Metabolism|^Protein import|^Small molecule|^Signaling|^OXPHOS") ~ "MitoCarta",
      str_detect(pathway, "^MYC_") ~ "MYC Signature",
      str_detect(pathway, "^MSigDB_HALLMARK_") ~ "Hallmark",
      TRUE ~ "Other"
    ),
    pathway_label = str_remove(pathway, "^(MSigDB_HALLMARK_|MYC_|MC_)") |>
      str_replace_all("_", " ")
  )

p_nes_comparison <- ggplot(
  fgsea_simple_data,
  aes(x = NES_neg, y = NES_pos)
) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", colour = "grey40") +
  geom_point(
    aes(size = sig_level, colour = pathway_category),
    alpha = 0.7
  ) +
  ggrepel::geom_text_repel(
    aes(label = pathway_label),
    size = 2.5,
    max.overlaps = Inf,
    segment.alpha = 0.3,
    segment.size = 0.2,
    force = 2,
    box.padding = 0.2
  ) +
  scale_colour_manual(values = c(
    "MitoCarta" = "#E41A1C",
    "MYC Signature" = "#377EB8",
    "Hallmark" = "#4DAF4A"
  )) +
  scale_size_continuous(range = c(2, 10), breaks = c(2, 5, 10, 20)) +
  labs(
    title = "Temporal pathway enrichment: Myc+ vs Myc-",
    x = "NES (Myc- 12W vs 6W)",
    y = "NES (Myc+ 12W vs 6W)",
    colour = "Category",
    size = "-log10(padj)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    legend.position = "right"
  )

ggsave(
  file.path(fig_dir, "xs_nes_comparison_simple.pdf"),
  p_nes_comparison, width = 14, height = 10
)

message("Simplified NES comparison plot saved")

# --- Plot 4: NES 6W vs NES 12W (cross-sectional Myc effect) ---
# X: NES (Myc+ vs Myc- at 6W)
# Y: NES (Myc+ vs Myc- at 12W)
# Size: significance, Colour: pathway category

fgsea_simple_xs <- fgsea_xs_combined |>
  as_tibble() |>
  mutate(
    min_padj = pmin(padj_6W, padj_12W, na.rm = TRUE),
    sig_level = -log10(min_padj),
    pathway_category = case_when(
      str_detect(pathway, "^MC_") ~ "MitoCarta",
      str_detect(pathway, "^MYC_") ~ "MYC Signature",
      str_detect(pathway, "^MSigDB_HALLMARK_") ~ "Hallmark",
      TRUE ~ "Other"
    ),
    pathway_label = str_remove(pathway, "^(MSigDB_HALLMARK_|MYC_|MC_)") |>
      str_replace_all("_", " ")
  )

p_nes_xs <- ggplot(
  fgsea_simple_xs,
  aes(x = NES_6W, y = NES_12W)
) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", colour = "grey40") +
  geom_point(
    aes(size = sig_level, colour = pathway_category),
    alpha = 0.7
  ) +
  ggrepel::geom_text_repel(
    aes(label = pathway_label),
    size = 2.5,
    max.overlaps = Inf,
    segment.alpha = 0.3,
    segment.size = 0.2,
    force = 2,
    box.padding = 0.2
  ) +
  scale_colour_manual(values = c(
    "MitoCarta" = "#E41A1C",
    "MYC Signature" = "#377EB8",
    "Hallmark" = "#4DAF4A"
  )) +
  scale_size_continuous(range = c(2, 10), breaks = c(2, 5, 10, 20)) +
  labs(
    title = "Cross-sectional Myc effect: 6W vs 12W",
    x = "NES (Myc+ vs Myc- at 6W)",
    y = "NES (Myc+ vs Myc- at 12W)",
    colour = "Category",
    size = "-log10(padj)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    legend.position = "right"
  )

ggsave(
  file.path(fig_dir, "xs_nes_crosssectional_comparison.pdf"),
  p_nes_xs, width = 14, height = 10
)

message("Cross-sectional NES comparison plot saved")

message("\n", strrep("=", 70))
message("CROSS-SECTIONAL fGSEA VISUALISATION COMPLETE")
message("Figures saved to: ", fig_dir)
message(strrep("=", 70))
