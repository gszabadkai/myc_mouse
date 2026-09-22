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
#   fGSEA sig + mitoPPS sig (same direction) → absolute AND relative change:
#                                              pathway specifically prioritised
#   fGSEA sig + mitoPPS NOT sig              → absolute change, proportional
#                                              (global mito amplification, e.g. biogenesis)
#   fGSEA sig + mitoPPS sig (opposite dir.)  → pathway increased absolutely but
#                                              deprioritised relative to others
#                                              (left behind by biogenesis programme)
#   fGSEA NOT sig + mitoPPS sig              → relative reallocation with small
#                                              but consistent gene-level shifts,
#                                              invisible against genome-wide background
#   Neither significant                      → no detected change
#
# Dotplot layout:
#   Each classification category produces ONE combined PDF with four panels:
#     Top row:    Myc+ vs Myc- at 6W  |  Myc+ vs Myc- at 12W
#     Bottom row: 12W vs 6W in Myc-   |  12W vs 6W in Myc+
#   Y-axis: shared within top row (mean effect at 6W+12W) and within bottom
#   row (mean effect in Myc-+Myc+) independently, so the Myc-effect and
#   temporal comparisons are directly visible side-by-side in one figure.
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

mitocarta_gene_sets <- gene_to_pathway %>%
  group_by(Pathway) %>%
  summarise(genes = list(unique(Gene)), .groups = "drop") %>%
  deframe()  # named list: pathway name → character vector of gene symbols

# Custom pathway names from script 08
MTDNA_PATHWAY_NAME  <- mitopps_results$mtdna_pathway_name
APOPTOSIS_PRO_NAME  <- "Apoptosis-PRO"
APOPTOSIS_ANTI_NAME <- "Apoptosis-ANTI"

# Ensure custom pathways are in tier1 map
pathway_tier1_map[MTDNA_PATHWAY_NAME]  <- "OXPHOS"
pathway_tier1_map[APOPTOSIS_PRO_NAME]  <- "Mitochondrial dynamics and surveillance"
pathway_tier1_map[APOPTOSIS_ANTI_NAME] <- "Mitochondrial dynamics and surveillance"

message(sprintf("  %d gene sets built (incl. '%s', '%s', '%s')",
                length(mitocarta_gene_sets),
                MTDNA_PATHWAY_NAME, APOPTOSIS_PRO_NAME, APOPTOSIS_ANTI_NAME))

for (pw in c(MTDNA_PATHWAY_NAME, APOPTOSIS_PRO_NAME, APOPTOSIS_ANTI_NAME)) {
  if (pw %in% names(mitocarta_gene_sets)) {
    message(sprintf("  \u2713 '%s': %d genes", pw, length(mitocarta_gene_sets[[pw]])))
  } else {
    warning(sprintf("  \u2717 '%s' NOT FOUND in gene sets — check mitopps_scores.rds", pw))
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
  Myc_effect_6W    = create_ranks(interaction_results$myc_6W_raw,   ensembl_to_symbol),
  Myc_effect_12W   = create_ranks(interaction_results$myc_12W_raw,  ensembl_to_symbol),
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

for (nm in names(fgsea_results)) {
  n_sig <- sum(fgsea_results[[nm]]$padj < 0.05, na.rm = TRUE)
  message(sprintf("  %-20s: %d pathways significant (padj < 0.05)", nm, n_sig))
}

# =============================================================================
# PART 5: ALIGN fGSEA AND mitoPPS RESULTS
# =============================================================================

message("\nAligning fGSEA and mitoPPS results...")

# Contrast name mapping: mitoPPS → fGSEA
contrast_map <- c(
  "Myc_effect_6W"   = "Myc_effect_6W",
  "Myc_effect_12W"  = "Myc_effect_12W",
  "Temporal_Myc-"   = "Temporal_Myc_neg",
  "Temporal_Myc+"   = "Temporal_Myc_pos"
)

# Friendly labels for facets — used consistently across all plots
contrast_labels <- c(
  "Myc_effect_6W"    = "Myc+ vs Myc\u2212  (6W)",
  "Myc_effect_12W"   = "Myc+ vs Myc\u2212  (12W)",
  "Temporal_Myc_neg" = "12W vs 6W  (Myc\u2212)",
  "Temporal_Myc_pos" = "12W vs 6W  (Myc+)"
)

# Contrast group membership — shared y-axis ordering within each group
myc_effect_contrasts <- c("Myc+ vs Myc\u2212  (6W)", "Myc+ vs Myc\u2212  (12W)")
temporal_contrasts   <- c("12W vs 6W  (Myc\u2212)", "12W vs 6W  (Myc+)")

# Combine all fGSEA results
fgsea_long <- bind_rows(fgsea_results) %>%
  dplyr::select(pathway, contrast, NES, padj_fgsea = padj, size) %>%
  mutate(sig_fgsea = padj_fgsea < 0.05)

# Prepare mitoPPS pairwise results
mitopps_long <- mitopps_pairwise %>%
  mutate(contrast_fgsea = contrast_map[contrast]) %>%
  filter(!is.na(contrast_fgsea)) %>%
  dplyr::select(pathway, contrast = contrast_fgsea,
                mitopps_diff = diff, padj_mitopps = padj) %>%
  mutate(sig_mitopps = padj_mitopps < 0.1)

# Join and classify
comparison <- fgsea_long %>%
  inner_join(mitopps_long, by = c("pathway", "contrast")) %>%
  mutate(
    tier1         = pathway_tier1_map[pathway],
    pathway_label = str_replace_all(pathway, "_", " "),
    category      = case_when(
      sig_fgsea & sig_mitopps & sign(NES) == sign(mitopps_diff) ~ "Both: concordant",
      sig_fgsea & sig_mitopps & sign(NES) != sign(mitopps_diff) ~ "Both: discordant",
      sig_fgsea  & !sig_mitopps ~ "fGSEA only\n(absolute change, proportional)",
      !sig_fgsea & sig_mitopps  ~ "mitoPPS only\n(reallocation without bulk DE)",
      TRUE ~ "Neither significant"
    ),
    contrast_label = factor(contrast_labels[contrast], levels = contrast_labels)
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
  "Metabolism"                                = "#2166AC",
  "OXPHOS"                                    = "#B2182B",
  "Protein import, sorting, and homeostasis"  = "#1B9E77",
  "Mitochondrial central dogma"               = "#D95F02",
  "Mitochondrial dynamics and surveillance"   = "#7570B3",
  "Small molecule transport"                  = "#E7298A",
  "Signaling"                                 = "#66A61E"
)

category_colours <- c(
  "Both: concordant"                             = "#B2182B",
  "Both: discordant"                             = "#FF7F00",
  "fGSEA only\n(absolute change, proportional)"  = "#2166AC",
  "mitoPPS only\n(reallocation without bulk DE)" = "#4DAF4A",
  "Neither significant"                          = "grey80"
)

# =============================================================================
# Helper: four-panel combined dotplot for one classification category.
#
# Layout (facet_wrap ncol = 2, scales = "free_y"):
#   Panel 1 (top-left):  Myc+ vs Myc- (6W)
#   Panel 2 (top-right): Myc+ vs Myc- (12W)
#   Panel 3 (bot-left):  12W vs 6W (Myc-)
#   Panel 4 (bot-right): 12W vs 6W (Myc+)
#
# Y-axis ordering:
#   Top row (Myc-effect): shared order by mean effect at 6W + 12W
#   Bottom row (temporal): shared order by mean effect in Myc- + Myc+
#   Implemented via [M]/[T] tagged factor levels that are stripped from
#   display labels by scale_y_discrete(labels = ...).
# =============================================================================

make_combined_dotplot <- function(dat, x_var, size_var, title, subtitle, x_label) {

  if (nrow(dat) == 0) return(NULL)

  # Shared ordering within each row
  order_myc <- dat %>%
    filter(contrast_label %in% myc_effect_contrasts) %>%
    group_by(pathway_label) %>%
    summarise(mean_eff = mean(.data[[x_var]], na.rm = TRUE), .groups = "drop") %>%
    arrange(mean_eff) %>%
    pull(pathway_label)

  order_temporal <- dat %>%
    filter(contrast_label %in% temporal_contrasts) %>%
    group_by(pathway_label) %>%
    summarise(mean_eff = mean(.data[[x_var]], na.rm = TRUE), .groups = "drop") %>%
    arrange(mean_eff) %>%
    pull(pathway_label)

  # Tag labels so the two rows have independent factor levels
  dat <- dat %>%
    mutate(
      pathway_tagged = case_when(
        contrast_label %in% myc_effect_contrasts ~ paste0(pathway_label, " [M]"),
        TRUE                                     ~ paste0(pathway_label, " [T]")
      ),
      pathway_tagged = factor(
        pathway_tagged,
        levels = c(paste0(order_myc,      " [M]"),
                   paste0(order_temporal, " [T]"))
      )
    )

  strip_tag <- function(x) sub(" \\[[MT]\\]$", "", x)

  ggplot(dat, aes(x = .data[[x_var]], y = pathway_tagged)) +
    geom_point(aes(size = -log10(.data[[size_var]]), colour = tier1), alpha = 0.85) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    facet_wrap(~ contrast_label, ncol = 2, scales = "free_y") +
    scale_y_discrete(labels = strip_tag) +
    scale_colour_manual(values = tier1_colours, na.value = "grey50",
                        name = "MitoCarta Tier 1") +
    scale_size_continuous(range = c(2, 8), name = expression(-log[10](padj))) +
    labs(title = title, subtitle = subtitle, x = x_label, y = NULL) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.y     = element_text(size = 7),
      plot.title      = element_text(face = "bold"),
      strip.text      = element_text(face = "bold", size = 10),
      legend.position = "right"
    )
}

# --- 6a. Scatterplot: NES vs mitoPPS diff, faceted by contrast ---
# Core comparison: all pathways plotted, coloured by classification category,
# with regression line and per-contrast correlation annotation.

comparison_plot <- comparison %>%
  mutate(
    combined_sig = sqrt(-log10(padj_fgsea + 1e-10) * -log10(padj_mitopps + 1e-10)),
    label_this   = (sig_fgsea | sig_mitopps) & category != "Neither significant"
  )

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
  geom_smooth(method = "lm", formula = y ~ x,
              colour = "grey30", fill = "grey80",
              linewidth = 0.5, alpha = 0.15, se = TRUE) +
  geom_point(data = filter(comparison_plot, category == "Neither significant"),
             colour = "grey85", size = 1.5, alpha = 0.5) +
  geom_point(data = filter(comparison_plot, category != "Neither significant"),
             aes(colour = category, size = combined_sig), alpha = 0.85) +
  geom_text_repel(
    data           = filter(comparison_plot, label_this),
    aes(label = pathway_label, colour = category),
    size           = 2.2,
    max.overlaps   = 20,
    segment.colour = "grey60",
    segment.size   = 0.3,
    show.legend    = FALSE
  ) +
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

# --- 6c. Correlation scatter: NES vs mitoPPS diff, coloured by tier1 ---
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
    subtitle = "Grey band: 95% CI of linear regression  |  Coloured by MitoCarta Tier 1",
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

# --- 6d-6g. Combined four-panel dotplots, one PDF per classification category ---
#
# Each PDF layout (facet_wrap ncol = 2, scales = "free_y"):
#   Top row:    Myc+ vs Myc- at 6W  |  Myc+ vs Myc- at 12W
#   Bottom row: 12W vs 6W in Myc-   |  12W vs 6W in Myc+
#
# Y-axis: shared within top row and shared within bottom row independently.

dotplot_specs <- list(

  list(
    category = "mitoPPS only\n(reallocation without bulk DE)",
    x_var    = "mitopps_diff",
    size_var = "padj_mitopps",
    title    = "mitoPPS-only pathways: relative reallocation without bulk DE signal",
    subtitle = paste0(
      "Significant in mitoPPS (padj < 0.1) but NOT fGSEA (padj \u2265 0.05)\n",
      "Small consistent shifts visible within mito compartment, sub-threshold genome-wide\n",
      "Top row: Myc effect (6W | 12W)  |  Bottom row: Temporal (Myc\u2212 | Myc+)  |  ",
      "Y-axis: shared order within each row by mean \u0394 mitoPPS"
    ),
    x_label  = "\u0394 mitoPPS",
    filename = "dotplot_mitopps_only.pdf"
  ),

  list(
    category = "fGSEA only\n(absolute change, proportional)",
    x_var    = "NES",
    size_var = "padj_fgsea",
    title    = "fGSEA-only pathways: absolute change proportional across mito compartment",
    subtitle = paste0(
      "Significant in fGSEA (padj < 0.05) but NOT mitoPPS (padj \u2265 0.1)\n",
      "Consistent with global mitochondrial biogenesis — all pathways amplified equally\n",
      "Top row: Myc effect (6W | 12W)  |  Bottom row: Temporal (Myc\u2212 | Myc+)  |  ",
      "Y-axis: shared order within each row by mean NES"
    ),
    x_label  = "fGSEA NES",
    filename = "dotplot_fgsea_only.pdf"
  ),

  list(
    category = "Both: concordant",
    x_var    = "mitopps_diff",
    size_var = "padj_mitopps",
    title    = "Concordant pathways: significant in both fGSEA and mitoPPS (same direction)",
    subtitle = paste0(
      "Absolute AND relative change — strongest evidence for selective mitochondrial targeting\n",
      "Top row: Myc effect (6W | 12W)  |  Bottom row: Temporal (Myc\u2212 | Myc+)  |  ",
      "Y-axis: shared order within each row by mean \u0394 mitoPPS"
    ),
    x_label  = "\u0394 mitoPPS",
    filename = "dotplot_concordant.pdf"
  ),

  list(
    category = "Both: discordant",
    x_var    = "mitopps_diff",
    size_var = "padj_mitopps",
    title    = "Discordant pathways: significant in both fGSEA and mitoPPS (opposite direction)",
    subtitle = paste0(
      "Pathway gains absolute transcript abundance but is relatively deprioritised\n",
      "Consistent with being left behind by global mitochondrial biogenesis\n",
      "Top row: Myc effect (6W | 12W)  |  Bottom row: Temporal (Myc\u2212 | Myc+)  |  ",
      "Y-axis: shared order within each row by mean \u0394 mitoPPS"
    ),
    x_label  = "\u0394 mitoPPS",
    filename = "dotplot_discordant.pdf"
  )
)

for (spec in dotplot_specs) {

  dat <- comparison %>% filter(category == spec$category)

  if (nrow(dat) == 0) {
    message(sprintf("  No pathways in '%s' — skipping", gsub("\n", " ", spec$category)))
    next
  }

  p <- make_combined_dotplot(
    dat      = dat,
    x_var    = spec$x_var,
    size_var = spec$size_var,
    title    = spec$title,
    subtitle = spec$subtitle,
    x_label  = spec$x_label
  )

  if (!is.null(p)) {
    n_myc <- dat %>%
      filter(contrast_label %in% myc_effect_contrasts) %>%
      pull(pathway_label) %>% unique() %>% length()
    n_temporal <- dat %>%
      filter(contrast_label %in% temporal_contrasts) %>%
      pull(pathway_label) %>% unique() %>% length()

    ggsave(file.path(comp_dir, spec$filename), p,
           width  = 14,
           height = max(6, (n_myc + n_temporal) * 0.3))
    message(sprintf("  Saved: %s", spec$filename))
  }
}

# =============================================================================
# PART 7: SAVE RESULTS
# =============================================================================

message("\nSaving results...")

comparison_results <- list(
  comparison          = comparison,
  fgsea_results       = fgsea_results,
  cor_by_contrast     = cor_by_contrast,
  mitocarta_gene_sets = mitocarta_gene_sets,
  analysis_date       = Sys.Date(),
  description         = paste(
    "Full MitoCarta fGSEA vs mitoPPS comparison.",
    "fGSEA ranked by DESeq2 Wald statistic, minSize = 3.",
    "mitoPPS threshold padj < 0.1; fGSEA threshold padj < 0.05.",
    "Dotplots: one combined four-panel PDF per classification category.",
    "Top row: Myc effect (6W | 12W); bottom row: temporal (Myc- | Myc+).",
    "Y-axis ordering shared within each row independently."
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
