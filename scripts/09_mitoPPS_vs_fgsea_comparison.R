# scripts/09_mitoPPS_vs_fgsea_comparison.R
# =============================================================================
# MitoPPS vs fGSEA Comparison: Full MitoCarta Pathway Analysis
# =============================================================================
#
# Purpose:
#   1. Run fGSEA on ALL MitoCarta3.0 pathways for all four contrasts
#      (matching exactly the four contrasts in the mitoPPS analysis)
#   2. Compare fGSEA NES with mitoPPS differences per pathway per contrast
#   3. Classify pathways by agreement/disagreement between methods
#   4. Draw biological conclusions about Myc-driven mitochondrial remodelling
#
# Contrasts (matching 08_mitoPPS_analysis.R):
#   - Myc+ vs Myc- at 6W    (cross-sectional)
#   - Myc+ vs Myc- at 12W   (cross-sectional)
#   - 12W vs 6W in Myc-     (temporal, baseline ageing)
#   - 12W vs 6W in Myc+     (temporal, Myc-driven)
#
# Key question:
#   Do fGSEA and mitoPPS agree? Disagreement reveals whether Myc drives
#   genuine mitochondrial REPRIORITISATION vs global transcriptional amplification.
#
# Interpretation framework:
#   fGSEA sig + mitoPPS sig (same direction) → absolute AND relative change
#   fGSEA sig + mitoPPS NOT sig              → absolute change, proportional
#                                              (global mito amplification)
#   fGSEA NOT sig + mitoPPS sig              → relative reallocation without
#                                              strong gene-level DE signal
#   Neither significant                      → no detected change
#
# Input:
#   - results/mitopps_scores.rds       (from 08_mitoPPS_analysis.R)
#   - results/interaction_results.rds  (DESeq2 Myc effect contrasts)
#   - results/group_results.rds        (DESeq2 temporal contrasts)
#   - results/ortholog_table.rds
#
# Output:
#   - results/mitopps_fgsea_comparison.rds
#   - outputs/mitopps_fgsea/
#
# =============================================================================

source(here::here("scripts", "00_setup_packages.R"))

# === Output directory ===
comp_dir <- here("outputs", "mitopps_fgsea")
dir.create(comp_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# PART 1: LOAD DATA
# =============================================================================

message("Loading data...")

# mitoPPS results — pairwise comparison diffs and pathway annotations
mitopps_results      <- readRDS(here("results", "mitopps_scores.rds"))
mitopps_pairwise     <- mitopps_results$mitopps_pairwise
gene_to_pathway      <- mitopps_results$gene_to_pathway
pathway_tier1_map    <- mitopps_results$pathway_tier1_map

# DESeq2 contrasts
interaction_results  <- readRDS(here("results", "interaction_results.rds"))
group_results        <- readRDS(here("results", "group_results.rds"))
ortholog_table       <- readRDS(here("results", "ortholog_table.rds"))

# Ensembl → symbol mapping
ensembl_to_symbol <- setNames(
  ortholog_table$external_gene_name,
  ortholog_table$ensembl_gene_id
)

# =============================================================================
# PART 2: BUILD FULL MitoCarta GENE SETS FOR fGSEA
# =============================================================================
# IMPORTANT: Use the MODIFIED gene_to_pathway from mitopps_scores.rds, which:
#   1. Has all mt-* genes REMOVED from their original pathways and placed in
#      the synthetic "mtDNA-encoded OXPHOS subunits" pathway
#   2. Has "Apoptosis-PRO" and "Apoptosis-ANTI" added as additional pathways
# This ensures fGSEA and mitoPPS use identical pathway gene set definitions.

message("Building MitoCarta gene sets for fGSEA (using modified gene_to_pathway)...")

# gene_to_pathway is already the modified version from mitopps_scores.rds
# (mtDNA separated, apoptosis subgroups added)
mitocarta_gene_sets <- gene_to_pathway %>%
  group_by(Pathway) %>%
  summarise(genes = list(unique(Gene)), .groups = "drop") %>%
  deframe()  # named list: pathway name → character vector of gene symbols

# Also update pathway_tier1_map with the custom pathways saved in results
MTDNA_PATHWAY_NAME  <- mitopps_results$mtdna_pathway_name
APOPTOSIS_PRO_NAME  <- "Apoptosis-PRO"
APOPTOSIS_ANTI_NAME <- "Apoptosis-ANTI"

# These should already be in pathway_tier1_map loaded from mitopps_results,
# but make explicit in case of any mismatch
pathway_tier1_map[MTDNA_PATHWAY_NAME]  <- "OXPHOS"
pathway_tier1_map[APOPTOSIS_PRO_NAME]  <- "Mitochondrial dynamics and surveillance"
pathway_tier1_map[APOPTOSIS_ANTI_NAME] <- "Mitochondrial dynamics and surveillance"

message(sprintf("  %d gene sets built (incl. '%s', '%s', '%s')",
                length(mitocarta_gene_sets),
                MTDNA_PATHWAY_NAME, APOPTOSIS_PRO_NAME, APOPTOSIS_ANTI_NAME))

# Verify the custom pathways are present
for (pw in c(MTDNA_PATHWAY_NAME, APOPTOSIS_PRO_NAME, APOPTOSIS_ANTI_NAME)) {
  if (pw %in% names(mitocarta_gene_sets)) {
    message(sprintf("  ✓ '%s': %d genes", pw, length(mitocarta_gene_sets[[pw]])))
  } else {
    warning(sprintf("  ✗ '%s' NOT FOUND in gene sets — check mitopps_scores.rds", pw))
  }
}

# =============================================================================
# PART 3: PREPARE RANKED GENE LISTS FOR ALL FOUR CONTRASTS
# =============================================================================

#' Create Wald-statistic ranked gene list from DESeq2 results
#' Identical to create_ranks() in scripts 04 and 06
create_ranks <- function(res, id_to_symbol) {
  as.data.frame(res) %>%
    rownames_to_column("ensembl_id") %>%
    filter(!is.na(stat)) %>%
    mutate(
      gene_symbol = id_to_symbol[ensembl_id],
      rank        = stat
    ) %>%
    filter(!is.na(gene_symbol) & gene_symbol != "") %>%
    group_by(gene_symbol) %>%
    slice_max(abs(rank), n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(desc(rank)) %>%
    { setNames(.$rank, .$gene_symbol) }
}

message("Creating ranked gene lists (Wald statistic)...")

ranks <- list(
  Myc_effect_6W   = create_ranks(interaction_results$myc_6W_raw,    ensembl_to_symbol),
  Myc_effect_12W  = create_ranks(interaction_results$myc_12W_raw,   ensembl_to_symbol),
  Temporal_Myc_neg = create_ranks(group_results$neg_12W_vs_6W_raw,  ensembl_to_symbol),
  Temporal_Myc_pos = create_ranks(group_results$pos_12W_vs_6W_raw,  ensembl_to_symbol)
)

for (nm in names(ranks)) {
  message(sprintf("  %-20s: %d genes ranked", nm, length(ranks[[nm]])))
}

# =============================================================================
# PART 4: RUN fGSEA ON ALL MITOCARTA PATHWAYS, ALL FOUR CONTRASTS
# =============================================================================

message("\nRunning fGSEA on all MitoCarta pathways...")
set.seed(42)

fgsea_results <- lapply(names(ranks), function(contrast_name) {
  message(sprintf("  Running: %s", contrast_name))
  fgsea(
    pathways    = mitocarta_gene_sets,
    stats       = ranks[[contrast_name]],
    minSize     = 3,    # match mitoPPS min_genes = 3
    maxSize     = 500,
    nPermSimple = 10000
  ) %>%
    mutate(contrast = contrast_name)
})

names(fgsea_results) <- names(ranks)

# Summary
for (nm in names(fgsea_results)) {
  n_sig <- sum(fgsea_results[[nm]]$padj < 0.05, na.rm = TRUE)
  message(sprintf("  %-20s: %d pathways significant (padj < 0.05)", nm, n_sig))
}

# =============================================================================
# PART 5: ALIGN fGSEA AND mitoPPS RESULTS
# =============================================================================
# Map contrast names between the two methods, then join on pathway name.
# fGSEA gives NES (normalised enrichment score, signed)
# mitoPPS gives diff (mean mitoPPS difference between groups, signed)
# Both are signed the same way (positive = up in condition B vs A)

message("\nAligning fGSEA and mitoPPS results...")

# mitoPPS contrast name mapping to fGSEA names
contrast_map <- c(
  "Myc_effect_6W"    = "Myc_effect_6W",
  "Myc_effect_12W"   = "Myc_effect_12W",
  "Temporal_Myc-"    = "Temporal_Myc_neg",
  "Temporal_Myc+"    = "Temporal_Myc_pos"
)

# Combine all fGSEA results into one long table
fgsea_long <- bind_rows(fgsea_results) %>%
  dplyr::select(pathway, contrast, NES, padj_fgsea = padj, size) %>%
  mutate(sig_fgsea = padj_fgsea < 0.05)

# Prepare mitoPPS pairwise results — rename contrast to match fGSEA
mitopps_long <- mitopps_pairwise %>%
  mutate(contrast_fgsea = contrast_map[contrast]) %>%
  filter(!is.na(contrast_fgsea)) %>%
  dplyr::select(pathway, contrast = contrast_fgsea,
                mitopps_diff = diff, padj_mitopps = padj) %>%
  mutate(sig_mitopps = padj_mitopps < 0.1)  # mitoPPS uses padj < 0.1 threshold

# Join on pathway + contrast
comparison <- fgsea_long %>%
  inner_join(mitopps_long, by = c("pathway", "contrast")) %>%
  mutate(
    tier1 = pathway_tier1_map[pathway],
    pathway_label = str_replace_all(pathway, "_", " "),
    # Classification by agreement
    category = case_when(
      sig_fgsea & sig_mitopps & sign(NES) == sign(mitopps_diff) ~ "Both: concordant",
      sig_fgsea & sig_mitopps & sign(NES) != sign(mitopps_diff) ~ "Both: discordant",
      sig_fgsea & !sig_mitopps ~ "fGSEA only\n(absolute change, proportional)",
      !sig_fgsea & sig_mitopps ~ "mitoPPS only\n(reallocation without bulk DE)",
      TRUE ~ "Neither significant"
    )
  )

message("\nPathway classification per contrast:")
comparison %>%
  group_by(contrast, category) %>%
  summarise(n = n(), .groups = "drop") %>%
  print(n = 40)

# =============================================================================
# PART 6: VISUALISATIONS
# =============================================================================

message("\nGenerating visualisations...")

tier1_colours <- c(
  "Metabolism"                                      = "#2166AC",
  "OXPHOS"                                          = "#B2182B",
  "Protein import, sorting, and homeostasis"        = "#1B9E77",
  "Mitochondrial central dogma"                     = "#D95F02",
  "Mitochondrial dynamics and surveillance"         = "#7570B3",
  "Small molecule transport"                        = "#E7298A",
  "Signaling"                                       = "#66A61E"
)

category_colours <- c(
  "Both: concordant"                        = "#B2182B",
  "Both: discordant"                        = "#FF7F00",
  "fGSEA only\n(absolute change, proportional)" = "#2166AC",
  "mitoPPS only\n(reallocation without bulk DE)" = "#4DAF4A",
  "Neither significant"                     = "grey80"
)

# Friendly contrast labels for plot facets
contrast_labels <- c(
  "Myc_effect_6W"    = "Myc+ vs Myc\u2212  (6W)",
  "Myc_effect_12W"   = "Myc+ vs Myc\u2212  (12W)",
  "Temporal_Myc_neg" = "12W vs 6W  (Myc\u2212)",
  "Temporal_Myc_pos" = "12W vs 6W  (Myc+)"
)

comparison <- comparison %>%
  mutate(contrast_label = factor(contrast_labels[contrast],
                                 levels = contrast_labels))

# --- 6a. Scatterplot: NES vs mitoPPS diff, faceted by contrast ---
# This is the core comparison plot.
# x = fGSEA NES (absolute transcriptional enrichment)
# y = mitoPPS diff (relative mitochondrial resource reallocation)
# Colour = classification category
# Size = mean significance (-log10 of geometric mean of both padj)

comparison_plot <- comparison %>%
  mutate(
    combined_sig = sqrt(-log10(padj_fgsea + 1e-10) * -log10(padj_mitopps + 1e-10)),
    label_this   = (sig_fgsea | sig_mitopps) & category != "Neither significant"
  )

# Compute per-contrast correlations BEFORE the ggplot call
cor_labels <- comparison %>%
  group_by(contrast_label) %>%
  summarise(
    r     = cor(NES, mitopps_diff, use = "complete.obs"),
    p_val = cor.test(NES, mitopps_diff)$p.value,
    .groups = "drop"
  ) %>%
  mutate(label = sprintf("r = %.2f\np = %.2e", r, p_val))

p_scatter_main <- ggplot(comparison_plot,
                         aes(x = NES, y = mitopps_diff)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey70", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey70", linewidth = 0.3) +
  # Regression line through ALL points — goes before geom_point so points sit on top
  geom_smooth(method = "lm", formula = y ~ x,
              colour = "grey30", fill = "grey80",
              linewidth = 0.5, alpha = 0.15, se = TRUE) +
  # Background: non-significant
  geom_point(data = filter(comparison_plot, category == "Neither significant"),
             colour = "grey85", size = 1.5, alpha = 0.5) +
  # Foreground: classified pathways
  geom_point(data = filter(comparison_plot, category != "Neither significant"),
             aes(colour = category, size = combined_sig), alpha = 0.85) +
  geom_text_repel(
    data          = filter(comparison_plot, label_this),
    aes(label = pathway_label, colour = category),
    size          = 2.2,
    max.overlaps  = 20,
    segment.colour = "grey60",
    segment.size  = 0.3,
    show.legend   = FALSE
  ) +
  # Correlation annotation — goes after geom_text_repel, before scales
  geom_text(
    data        = cor_labels,
    aes(x = -Inf, y = Inf, label = label),
    hjust       = -0.1,
    vjust       = 1.3,
    size        = 3,
    colour      = "grey20",
    inherit.aes = FALSE
  ) +
  facet_wrap(~ contrast_label, ncol = 2) +
  scale_colour_manual(values = category_colours, name = "Classification") +
  scale_size_continuous(range = c(1.5, 6),
                        name  = expression(sqrt(-log[10](padj)))) +
  labs(
    title    = "fGSEA vs mitoPPS: MitoCarta pathway comparison",
    subtitle = paste0(
      "x-axis: fGSEA NES (absolute transcriptional enrichment)\n",
      "y-axis: \u0394 mitoPPS (relative mitochondrial resource reallocation)\n",
      "Grey band: 95% CI of regression across all pathways"
    ),
    x = "fGSEA NES",
    y = "\u0394 mitoPPS (Myc+ \u2212 Myc\u2212  or  12W \u2212 6W)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title      = element_text(face = "bold"),
    strip.text      = element_text(face = "bold", size = 10),
    legend.position = "bottom",
    legend.text     = element_text(size = 8)
  )

ggsave(file.path(comp_dir, "scatter_fgsea_vs_mitopps_all_contrasts.pdf"),
       p_scatter_main, width = 14, height = 12)

# --- 6b. Classification summary barplot ---
summary_counts <- comparison %>%
  group_by(contrast_label, category) %>%
  summarise(n = n(), .groups = "drop") %>%
  filter(category != "Neither significant")

p_barplot <- ggplot(summary_counts,
                    aes(x = contrast_label, y = n, fill = category)) +
  geom_col(position = "dodge", alpha = 0.85) +
  scale_fill_manual(values = category_colours, name = NULL) +
  labs(
    title    = "Pathway classification: fGSEA vs mitoPPS",
    subtitle = "Number of pathways in each agreement category per contrast",
    x        = NULL,
    y        = "Number of pathways"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title      = element_text(face = "bold"),
    axis.text.x     = element_text(angle = 30, hjust = 1),
    legend.position = "bottom",
    legend.text     = element_text(size = 8)
  )

ggsave(file.path(comp_dir, "barplot_classification_summary.pdf"),
       p_barplot, width = 10, height = 6)

# --- 6c. Correlation between NES and mitoPPS diff per contrast ---
cor_by_contrast <- comparison %>%
  group_by(contrast, contrast_label) %>%
  summarise(
    r       = cor(NES, mitopps_diff, use = "complete.obs"),
    p_value = cor.test(NES, mitopps_diff)$p.value,
    n       = n(),
    .groups = "drop"
  )

message("\nCorrelation between fGSEA NES and mitoPPS diff per contrast:")
print(cor_by_contrast)

p_cor_scatter <- ggplot(comparison,
                        aes(x = NES, y = mitopps_diff, colour = tier1)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey70", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey70", linewidth = 0.3) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", formula = y ~ x, colour = "grey30",
              fill = "grey80", linewidth = 0.7, alpha = 0.25, se = TRUE) +
  geom_text(
    data = cor_by_contrast,
    aes(x = -Inf, y = Inf,
        label = sprintf("r = %.2f\np = %.2e", r, p_value)),
    hjust = -0.1, vjust = 1.3,
    colour = "grey20", size = 3, inherit.aes = FALSE
  ) +
  facet_wrap(~ contrast_label, ncol = 2) +
  scale_colour_manual(values = tier1_colours, na.value = "grey60",
                      name = "MitoCarta Tier 1") +
  labs(
    title    = "fGSEA NES vs \u0394 mitoPPS correlation per contrast",
    subtitle = "Grey band: 95% CI of linear regression",
    x        = "fGSEA NES",
    y        = "\u0394 mitoPPS"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title      = element_text(face = "bold"),
    strip.text      = element_text(face = "bold", size = 10),
    legend.position = "right"
  )

ggsave(file.path(comp_dir, "scatter_nes_vs_mitopps_correlation.pdf"),
       p_cor_scatter, width = 14, height = 12)

# --- 6d. Dotplot: pathways significant in mitoPPS only (reallocation signal) ---
# These are the most interesting: mitochondrial reprioritisation NOT captured
# by conventional enrichment analysis
mitopps_only <- comparison %>%
  filter(category == "mitoPPS only\n(reallocation without bulk DE)") %>%
  mutate(pathway_label = str_replace_all(pathway, "_", " "))

if (nrow(mitopps_only) > 0) {
  p_mitopps_only <- ggplot(mitopps_only,
                           aes(x = mitopps_diff, y = reorder(pathway_label, mitopps_diff))) +
    geom_point(aes(size = -log10(padj_mitopps), colour = tier1), alpha = 0.85) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    facet_wrap(~ contrast_label, ncol = 2) +
    scale_colour_manual(values = tier1_colours, na.value = "grey50",
                        name = "MitoCarta Tier 1") +
    scale_size_continuous(range = c(2, 8), name = expression(-log[10](padj))) +
    labs(
      title    = "mitoPPS-only pathways: reallocation without bulk DE",
      subtitle = "Significant in mitoPPS (padj < 0.1) but NOT in fGSEA (padj \u2265 0.05)\nThese pathways are reprioritised independently of global transcriptional changes",
      x        = "\u0394 mitoPPS",
      y        = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.y     = element_text(size = 7),
      plot.title      = element_text(face = "bold"),
      strip.text      = element_text(face = "bold", size = 10),
      legend.position = "right"
    )

  ggsave(file.path(comp_dir, "dotplot_mitopps_only_pathways.pdf"),
         p_mitopps_only,
         width  = 14,
         height = max(5, length(unique(mitopps_only$pathway_label)) * 0.3))
}

# --- 6e. Dotplot: pathways significant in fGSEA only (proportional amplification) ---
fgsea_only <- comparison %>%
  filter(category == "fGSEA only\n(absolute change, proportional)") %>%
  mutate(pathway_label = str_replace_all(pathway, "_", " "))

if (nrow(fgsea_only) > 0) {
  p_fgsea_only <- ggplot(fgsea_only,
                         aes(x = NES, y = reorder(pathway_label, NES))) +
    geom_point(aes(size = -log10(padj_fgsea), colour = tier1), alpha = 0.85) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    facet_wrap(~ contrast_label, ncol = 2) +
    scale_colour_manual(values = tier1_colours, na.value = "grey50",
                        name = "MitoCarta Tier 1") +
    scale_size_continuous(range = c(2, 8), name = expression(-log[10](padj))) +
    labs(
      title    = "fGSEA-only pathways: absolute change, proportional across mito compartment",
      subtitle = "Significant in fGSEA (padj < 0.05) but NOT in mitoPPS (padj \u2265 0.1)\nLikely driven by global mitochondrial content change rather than selective reallocation",
      x        = "fGSEA NES",
      y        = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.y     = element_text(size = 7),
      plot.title      = element_text(face = "bold"),
      strip.text      = element_text(face = "bold", size = 10),
      legend.position = "right"
    )

  ggsave(file.path(comp_dir, "dotplot_fgsea_only_pathways.pdf"),
         p_fgsea_only,
         width  = 14,
         height = max(5, length(unique(fgsea_only$pathway_label)) * 0.3))
}

# =============================================================================
# PART 7: SAVE RESULTS
# =============================================================================

message("\nSaving results...")

comparison_results <- list(
  comparison        = comparison,
  fgsea_results     = fgsea_results,
  cor_by_contrast   = cor_by_contrast,
  mitocarta_gene_sets = mitocarta_gene_sets,
  analysis_date     = Sys.Date(),
  description       = paste(
    "Full MitoCarta fGSEA vs mitoPPS comparison.",
    "fGSEA ranked by DESeq2 Wald statistic, minSize = 3.",
    "mitoPPS threshold padj < 0.1; fGSEA threshold padj < 0.05."
  )
)

saveRDS(comparison_results, here("results", "mitopps_fgsea_comparison.rds"))

comparison %>%
  dplyr::select(pathway, pathway_label, tier1, contrast, contrast_label,
                NES, padj_fgsea, sig_fgsea,
                mitopps_diff, padj_mitopps, sig_mitopps,
                category) %>%
  arrange(contrast, category, padj_fgsea) %>%
  write.csv(file.path(comp_dir, "fgsea_vs_mitopps_full_comparison.csv"),
            row.names = FALSE)

message("\n", strrep("=", 70))
message("fGSEA vs mitoPPS COMPARISON COMPLETE")
message(strrep("=", 70))
message(sprintf("Results saved to: results/mitopps_fgsea_comparison.rds"))
message(sprintf("Figures saved to: %s", comp_dir))
message(strrep("=", 70))