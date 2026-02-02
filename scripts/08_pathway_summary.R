# scripts/08_pathway_summary.R
# MitoCarta pathway-level summary analysis
# Addresses Q1: How does the mitochondrial transcriptome change across conditions?

source("scripts/00_setup_packages.R")
library(tidyr)

# === Load data ===
combined_df <- readRDS("results/combined_df_annotated.rds")

# === Parse MitoCarta pathways ===
# Format: pathway name in column 1, comma-separated genes in column 2
mitocarta_raw <- read_csv("data/mitocarta_pathways.csv", col_names = FALSE, show_col_types = FALSE)

mitocarta_long <- mitocarta_raw |>
  dplyr::rename(pathway = X1, genes_str = X2) |>
  separate_longer_delim(genes_str, delim = ", ") |>
  dplyr::rename(gene = genes_str) |>
  mutate(gene = str_trim(gene)) |>
  filter(gene != "")

# Save for reuse
saveRDS(mitocarta_long, "results/mitocarta_long.rds")

message(sprintf("Parsed %d gene-pathway pairs across %d pathways", 
                nrow(mitocarta_long), n_distinct(mitocarta_long$pathway)))

# === Join pathways with DE results ===
mito_de <- combined_df |>
  inner_join(mitocarta_long, by = c("mgi_symbol" = "gene"), relationship = "many-to-many")

message(sprintf("Matched %d gene-pathway pairs from %d unique genes", 
                nrow(mito_de), n_distinct(mito_de$mgi_symbol)))

# === Pathway-level summary function ===
summarize_pathway <- function(df, lfc_col, padj_col = NULL, stat_col = NULL, 
                               lfc_thresh = 0.3, padj_thresh = 0.1, stat_thresh = 2) {
  lfc_sym <- sym(lfc_col)
  
  # Determine significance: use padj if available, otherwise use stat
  if (!is.null(padj_col) && padj_col %in% names(df)) {
    padj_sym <- sym(padj_col)
    df <- df |>
      mutate(.sig = !is.na(!!padj_sym) & !!padj_sym < padj_thresh)
  } else if (!is.null(stat_col) && stat_col %in% names(df)) {
    stat_sym <- sym(stat_col)
    df <- df |>
      mutate(.sig = !is.na(!!stat_sym) & abs(!!stat_sym) > stat_thresh)
  } else {
    df <- df |>
      mutate(.sig = FALSE)
  }
  
  df |>
    group_by(pathway) |>
    dplyr::summarize(
      n_genes = n(),
      mean_lfc = mean(!!lfc_sym, na.rm = TRUE),
      median_lfc = median(!!lfc_sym, na.rm = TRUE),
      sd_lfc = sd(!!lfc_sym, na.rm = TRUE),
      n_up = sum(!!lfc_sym > lfc_thresh & .sig, na.rm = TRUE),
      n_down = sum(!!lfc_sym < -lfc_thresh & .sig, na.rm = TRUE),
      n_sig = n_up + n_down,
      pct_sig = 100 * n_sig / n(),
      # One-sample t-test against 0
      t_pval = if (n() > 2 & sd(!!lfc_sym, na.rm = TRUE) > 0) {
        t.test(!!lfc_sym, mu = 0)$p.value
      } else NA_real_,
      .groups = "drop"
    ) |>
    mutate(
      direction = case_when(
        mean_lfc > 0.1 ~ "up",
        mean_lfc < -0.1 ~ "down",
        TRUE ~ "unchanged"
      )
    )
}

# === Generate summaries for each contrast ===

# Baseline aging: 6W -> 12W in MYC- condition
pathway_baseline <- summarize_pathway(
  mito_de, 
  lfc_col = "timepoint_neg_log2FC", 
  padj_col = "tau_padj",
  stat_col = "tau_stat"
) |>
  dplyr::rename_with(~ paste0("baseline_", .), -pathway)

# MYC+ aging: 6W -> 12W in MYC+ condition  
pathway_myc_aging <- summarize_pathway(
  mito_de,
  lfc_col = "timepoint_pos_log2FC",
  stat_col = "delta_stat"  # No direct padj available
) |>
  dplyr::rename_with(~ paste0("myc_aging_", .), -pathway)

# MYC effect at 6W
pathway_myc_6W <- summarize_pathway(
  mito_de,
  lfc_col = "myc_6W_log2FC",
  padj_col = "mu6_padj",
  stat_col = "mu6_stat"
) |>
  dplyr::rename_with(~ paste0("myc_6W_", .), -pathway)

# MYC effect at 12W
pathway_myc_12W <- summarize_pathway(
  mito_de,
  lfc_col = "myc_12W_log2FC",
  padj_col = "mu12_padj",
  stat_col = "mu12_stat"
) |>
  dplyr::rename_with(~ paste0("myc_12W_", .), -pathway)

# === Combine all summaries ===
pathway_summary <- pathway_baseline |>
  left_join(pathway_myc_aging, by = "pathway") |>
  left_join(pathway_myc_6W, by = "pathway") |>
  left_join(pathway_myc_12W, by = "pathway")

# Save results
write_csv(pathway_summary, "results/mitocarta_pathway_summary.csv")
saveRDS(pathway_summary, "results/mitocarta_pathway_summary.rds")
saveRDS(mito_de, "results/mito_de.rds")

message("Pathway summaries saved to results/")

# === Visualization: Pathway-level heatmap ===
dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)
# Prepare matrix of mean LFCs
pathway_lfc_matrix <- pathway_summary |>
  dplyr::select(pathway, 
         `Baseline (6W-12W)` = baseline_mean_lfc,
         `MYC+ (6W-12W)` = myc_aging_mean_lfc,
         `MYC effect 6W` = myc_6W_mean_lfc,
         `MYC effect 12W` = myc_12W_mean_lfc) |>
  column_to_rownames("pathway") |>
  as.matrix()

# Shorten pathway names for display
rownames(pathway_lfc_matrix) <- rownames(pathway_lfc_matrix) |>
  str_replace("Mitochondrial central dogma > ", "Central dogma: ") |>
  str_replace("Protein import, sorting and homeostasis > ", "") |>
  str_replace("Protein homeostasis > ", "") |>
  str_replace("Metabolism > ", "") |>
  str_replace("OXPHOS > ", "OXPHOS: ")

# Color palette
lfc_colors <- colorRamp2(
  c(-0.3, -0.15, 0, 0.15, 0.3),
  c("#2166AC", "#92C5DE", "white", "#F4A582", "#B2182B")
)

# Generate heatmap
ht <- Heatmap(
  pathway_lfc_matrix,
  name = "Mean LFC",
  col = lfc_colors,
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  row_names_gp = gpar(fontsize = 9),
  column_names_gp = gpar(fontsize = 10),
  column_names_rot = 45,
  row_title = "MitoCarta Pathways",
  column_title = "Mitochondrial Pathway Changes",
  heatmap_legend_param = list(
    title = "Mean Log2FC",
    at = c(-0.3, -0.15, 0, 0.15, 0.3)
  ),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.2f", pathway_lfc_matrix[i, j]), x, y, gp = gpar(fontsize = 8))
  }
)

pdf("results/figures/mitocarta_pathway_heatmap.pdf", width = 9, height = 10)
draw(ht)
dev.off()

# === Visualization: Dot plot (alternative) ===
pathway_long <- pathway_summary |>
  dplyr::select(pathway, 
         baseline_mean_lfc, baseline_pct_sig,
         myc_aging_mean_lfc, myc_aging_pct_sig,
         myc_6W_mean_lfc, myc_6W_pct_sig,
         myc_12W_mean_lfc, myc_12W_pct_sig) |>
  pivot_longer(
    cols = -pathway,
    names_to = c("contrast", ".value"),
    names_pattern = "(.+)_(mean_lfc|pct_sig)"
  ) |>
  mutate(
    contrast = case_when(
      contrast == "baseline" ~ "Baseline (6W→12W)",
      contrast == "myc_aging" ~ "MYC+ (6W→12W)",
      contrast == "myc_6W" ~ "MYC effect @ 6W",
      contrast == "myc_12W" ~ "MYC effect @ 12W"
    ),
    contrast = factor(contrast, levels = c(
      "Baseline (6W→12W)", "MYC+ (6W→12W)", 
      "MYC effect @ 6W", "MYC effect @ 12W"
    ))
  )

p_dotplot <- ggplot(pathway_long, aes(x = contrast, y = pathway)) +

  geom_point(aes(size = pct_sig, color = mean_lfc)) +
  scale_color_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B",
    midpoint = 0, name = "Mean LFC",
    limits = c(-1, 1), oob = scales::squish
  ) +
  scale_size_continuous(name = "% Significant", range = c(1, 8)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 7),
    panel.grid.major = element_line(color = "grey90")
  ) +
  labs(
    x = NULL, y = NULL,
    title = "MitoCarta Pathway Changes",
    subtitle = "Dot size = % significant genes, Color = mean log2FC"
  )

ggsave("results/figures/mitocarta_pathway_dotplot.pdf", p_dotplot, 
       width = 10, height = 14)

message("Pathway summary complete. Results saved to results/")
