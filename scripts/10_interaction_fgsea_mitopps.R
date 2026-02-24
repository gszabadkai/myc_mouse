# scripts/10_interaction_fgsea_mitopps.R
# =============================================================================
# Interaction Analysis: Does Myc's Effect on Mitochondrial Pathways
# Evolve Between 6W and 12W?
# =============================================================================
#
# Motivation:
#   The DESeq2 interaction term (timepoint12W:myc_statuspos) found NO
#   individually significant genes (minimum padj ~ 0.16). However, the
#   interaction term still carries Wald statistics for every gene — and fGSEA
#   can detect coordinated pathway-level shifts even when no single gene
#   passes FDR. This is exactly the scenario where pathway analysis adds
#   power beyond gene-level testing.
#
# This script asks:
#   1. Does fGSEA on the interaction Wald statistics detect any pathways
#      whose Myc effect changes from 6W → 12W?
#   2. Does the mitoPPS interaction (Δ mitoPPS_Myc_12W − Δ mitoPPS_Myc_6W)
#      agree with fGSEA?
#   3. Are the same or different pathways flagged?
#
# Interpretation of the interaction term:
#   The DESeq2 interaction coefficient tests:
#     (Myc_effect_at_12W) − (Myc_effect_at_6W)
#   Positive Wald stat → Myc effect is STRONGER (or more positive) at 12W
#   Negative Wald stat → Myc effect is WEAKER (or more negative) at 12W
#
#   Similarly, the mitoPPS interaction is:
#     Δ mitoPPS_Myc_12W − Δ mitoPPS_Myc_6W
#     where Δ mitoPPS_Myc = mean(Myc+) − mean(Myc−) at that timepoint
#   Positive → pathway reprioritisation by Myc increases over time
#   Negative → pathway reprioritisation by Myc decreases over time
#
# Biological question:
#   Does Myc progressively remodel mitochondrial resource allocation, or is
#   the reprioritisation established early (6W) and stable? Pathways with
#   significant interaction in BOTH fGSEA and mitoPPS represent the strongest
#   evidence for time-dependent Myc-driven mitochondrial remodelling.
#
# Input:
#   - results/interaction_results.rds   (DESeq2, contains interaction_raw)
#   - results/mitopps_scores.rds        (mitoPPS scores + pairwise results)
#   - results/gene_sets_list.rds        (Hallmark/MYC/MitoCarta sets)
#   - results/ortholog_table.rds
#
# Output:
#   - results/interaction_fgsea_mitopps.rds
#   - outputs/interaction_analysis/
#
# =============================================================================

source(here::here("scripts", "00_setup_packages.R"))

# === Output directory ===
int_dir <- here("outputs", "interaction_analysis")
dir.create(int_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# PART 1: LOAD DATA
# =============================================================================

message("Loading data...")

# DESeq2 interaction results
interaction_results <- readRDS(here("results", "interaction_results.rds"))
res_interaction     <- interaction_results$interaction_raw

# mitoPPS results
mitopps_results   <- readRDS(here("results", "mitopps_scores.rds"))
mitopps_pairwise  <- mitopps_results$mitopps_pairwise
mitopps_anova     <- mitopps_results$mitopps_anova
gene_to_pathway   <- mitopps_results$gene_to_pathway
pathway_tier1_map <- mitopps_results$pathway_tier1_map

# Gene sets — both the curated set (Hallmark/MYC/MitoCarta) and full MitoCarta
gene_sets     <- readRDS(here("results", "gene_sets_list.rds"))
ortholog_table <- readRDS(here("results", "ortholog_table.rds"))

# Ensembl → symbol mapping
ensembl_to_symbol <- setNames(
  ortholog_table$external_gene_name,
  ortholog_table$ensembl_gene_id
)

# =============================================================================
# PART 2: RANK GENES BY INTERACTION WALD STATISTIC
# =============================================================================

message("Creating ranked gene list from interaction Wald statistics...")

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

ranks_interaction <- create_ranks(res_interaction, ensembl_to_symbol)

message(sprintf("  %d genes ranked by interaction Wald statistic", length(ranks_interaction)))
message(sprintf("  Range: [%.2f, %.2f]", min(ranks_interaction), max(ranks_interaction)))

# Quick diagnostic: is the ranking informative?
# Under the null (no interaction), the Wald stats should be ~N(0,1)
# A fat tail or skew suggests real signal even if no gene passes FDR
ks_test <- ks.test(ranks_interaction, "pnorm", mean = 0, sd = 1)
message(sprintf("  KS test vs N(0,1): D = %.4f, p = %.2e (deviation from null)",
                ks_test$statistic, ks_test$p.value))

# =============================================================================
# PART 3: fGSEA ON INTERACTION — CURATED GENE SETS
# =============================================================================
# Same gene sets as scripts 04/06: Hallmark, MYC, MitoCarta, Apoptosis

message("\nRunning fGSEA on interaction term — curated gene sets...")
set.seed(42)

fgsea_interaction_curated <- fgsea(
  pathways    = gene_sets,
  stats       = ranks_interaction,
  minSize     = 10,
  maxSize     = 500,
  nPermSimple = 10000
)

n_sig_curated <- sum(fgsea_interaction_curated$padj < 0.05, na.rm = TRUE)
n_sig_curated_10 <- sum(fgsea_interaction_curated$padj < 0.10, na.rm = TRUE)
message(sprintf("  Significant: %d (padj < 0.05), %d (padj < 0.10) out of %d sets",
                n_sig_curated, n_sig_curated_10, nrow(fgsea_interaction_curated)))

# =============================================================================
# PART 4: fGSEA ON INTERACTION — FULL MitoCarta PATHWAYS
# =============================================================================
# Same modified gene sets as script 09 (with mtDNA separation and
# Apoptosis-PRO/ANTI), ensuring identical pathway definitions to mitoPPS

message("Running fGSEA on interaction term — full MitoCarta pathways...")

mitocarta_gene_sets <- gene_to_pathway %>%
  group_by(Pathway) %>%
  summarise(genes = list(unique(Gene)), .groups = "drop") %>%
  deframe()

# Custom pathway names from script 08
MTDNA_PATHWAY_NAME  <- mitopps_results$mtdna_pathway_name
APOPTOSIS_PRO_NAME  <- "Apoptosis-PRO"
APOPTOSIS_ANTI_NAME <- "Apoptosis-ANTI"

# Ensure custom pathways are in tier1 map
pathway_tier1_map[MTDNA_PATHWAY_NAME]  <- "OXPHOS"
pathway_tier1_map[APOPTOSIS_PRO_NAME]  <- "Mitochondrial dynamics and surveillance"
pathway_tier1_map[APOPTOSIS_ANTI_NAME] <- "Mitochondrial dynamics and surveillance"

message(sprintf("  %d MitoCarta gene sets", length(mitocarta_gene_sets)))

fgsea_interaction_mitocarta <- fgsea(
  pathways    = mitocarta_gene_sets,
  stats       = ranks_interaction,
  minSize     = 3,      # match mitoPPS min_genes = 3
  maxSize     = 500,
  nPermSimple = 10000
)

n_sig_mc <- sum(fgsea_interaction_mitocarta$padj < 0.05, na.rm = TRUE)
n_sig_mc_10 <- sum(fgsea_interaction_mitocarta$padj < 0.10, na.rm = TRUE)
message(sprintf("  Significant: %d (padj < 0.05), %d (padj < 0.10) out of %d sets",
                n_sig_mc, n_sig_mc_10, nrow(fgsea_interaction_mitocarta)))

# =============================================================================
# PART 5: COMPUTE mitoPPS INTERACTION TERM
# =============================================================================
# The mitoPPS interaction is the difference-of-differences:
#   Δ_interaction = (mean_12W_pos − mean_12W_neg) − (mean_6W_pos − mean_6W_neg)
#
# This is computed from the existing pairwise results:
#   Δ_Myc_12W = diff from Myc_effect_12W contrast
#   Δ_Myc_6W  = diff from Myc_effect_6W contrast
#   Interaction = Δ_Myc_12W − Δ_Myc_6W
#
# For the p-value, we use the two-way ANOVA interaction term already computed
# in script 08 (mitopps_anova, effect == "timepoint:myc_status").

message("\nComputing mitoPPS interaction (difference of Myc effects: 12W − 6W)...")

pps_myc_6W <- mitopps_pairwise %>%
  filter(contrast == "Myc_effect_6W") %>%
  dplyr::select(pathway, diff_6W = diff, padj_6W = padj)

pps_myc_12W <- mitopps_pairwise %>%
  filter(contrast == "Myc_effect_12W") %>%
  dplyr::select(pathway, diff_12W = diff, padj_12W = padj)

pps_interaction <- pps_myc_6W %>%
  inner_join(pps_myc_12W, by = "pathway") %>%
  mutate(
    interaction_diff = diff_12W - diff_6W  # positive = Myc effect grows over time
  )

# Add the ANOVA interaction p-value from script 08
anova_interaction <- mitopps_anova %>%
  filter(effect == "timepoint:myc_status") %>%
  dplyr::select(pathway, anova_F = F_value, anova_p = p_value, anova_padj = padj)

pps_interaction <- pps_interaction %>%
  left_join(anova_interaction, by = "pathway") %>%
  mutate(sig_anova = anova_padj < 0.10)

n_sig_anova <- sum(pps_interaction$sig_anova, na.rm = TRUE)
message(sprintf("  Pathways with significant ANOVA interaction (padj < 0.10): %d / %d",
                n_sig_anova, nrow(pps_interaction)))

# =============================================================================
# PART 6: ALIGN fGSEA AND mitoPPS INTERACTION RESULTS
# =============================================================================

message("Aligning fGSEA interaction (MitoCarta) with mitoPPS interaction...")

fgsea_int_mc <- fgsea_interaction_mitocarta %>%
  dplyr::select(pathway, NES, padj_fgsea = padj, size) %>%
  mutate(sig_fgsea = padj_fgsea < 0.05)

comparison <- fgsea_int_mc %>%
  inner_join(pps_interaction, by = "pathway") %>%
  mutate(
    tier1         = pathway_tier1_map[pathway],
    pathway_label = str_replace_all(pathway, "_", " "),
    # Classification: same framework as script 09 but for interaction
    category = case_when(
      sig_fgsea & sig_anova & sign(NES) == sign(interaction_diff)
        ~ "Both: concordant\n(evolving Myc reprioritisation)",
      sig_fgsea & sig_anova & sign(NES) != sign(interaction_diff)
        ~ "Both: discordant",
      sig_fgsea  & !sig_anova
        ~ "fGSEA only\n(transcriptional interaction, stable mitoPPS)",
      !sig_fgsea & sig_anova
        ~ "mitoPPS only\n(reallocation shift, sub-threshold genome-wide)",
      TRUE ~ "Neither significant"
    )
  )

message("\nInteraction pathway classification:")
comparison %>%
  group_by(category) %>%
  summarise(n = n(), .groups = "drop") %>%
  print()

# =============================================================================
# PART 7: COMPARE WITH ORIGINAL FOUR CONTRASTS
# =============================================================================
# For context: show the original Myc effect at each timepoint alongside
# the interaction, so we can see whether pathways that interact are
# "emerging" (small at 6W, large at 12W) vs "fading" vs "reversing"

message("Adding original pairwise context...")

comparison <- comparison %>%
  mutate(
    myc_effect_pattern = case_when(
      abs(diff_6W) < 0.005 & abs(diff_12W) >= 0.005
        ~ "Emerging (6W ≈ 0, 12W active)",
      abs(diff_6W) >= 0.005 & abs(diff_12W) < 0.005
        ~ "Fading (6W active, 12W ≈ 0)",
      sign(diff_6W) != sign(diff_12W) & abs(diff_6W) >= 0.005 & abs(diff_12W) >= 0.005
        ~ "Reversing (opposite at 6W vs 12W)",
      sign(diff_6W) == sign(diff_12W) & abs(interaction_diff) >= 0.005
        ~ "Amplifying (same direction, magnitude changes)",
      TRUE ~ "Stable (minimal interaction)"
    )
  )

message("\nMyc effect temporal patterns:")
comparison %>%
  group_by(myc_effect_pattern) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(desc(n)) %>%
  print()

# =============================================================================
# PART 8: VISUALISATIONS
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
  "Both: concordant\n(evolving Myc reprioritisation)" = "#B2182B",
  "Both: discordant"                                   = "#FF7F00",
  "fGSEA only\n(transcriptional interaction, stable mitoPPS)" = "#2166AC",
  "mitoPPS only\n(reallocation shift, sub-threshold genome-wide)" = "#4DAF4A",
  "Neither significant"                                = "grey80"
)

# --- 8a. Scatter: fGSEA NES vs mitoPPS interaction diff ---

cor_result <- cor.test(comparison$NES, comparison$interaction_diff,
                       use = "complete.obs")

p_scatter <- ggplot(comparison, aes(x = NES, y = interaction_diff)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey70", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey70", linewidth = 0.3) +
  geom_smooth(method = "lm", formula = y ~ x,
              colour = "grey30", fill = "grey80",
              linewidth = 0.5, alpha = 0.15, se = TRUE) +
  geom_point(data = filter(comparison, category == "Neither significant"),
             colour = "grey85", size = 2, alpha = 0.5) +
  geom_point(data = filter(comparison, category != "Neither significant"),
             aes(colour = category, size = -log10(pmin(padj_fgsea, anova_padj))),
             alpha = 0.85) +
  ggrepel::geom_text_repel(
    data           = filter(comparison, sig_fgsea | sig_anova),
    aes(label = pathway_label, colour = category),
    size           = 2.5,
    max.overlaps   = 25,
    segment.colour = "grey60",
    segment.size   = 0.3,
    show.legend    = FALSE
  ) +
  annotate("text",
           x = -Inf, y = Inf,
           label = sprintf("r = %.2f\np = %.2e", cor_result$estimate, cor_result$p.value),
           hjust = -0.1, vjust = 1.3,
           size = 3.5, colour = "grey20") +
  scale_colour_manual(values = category_colours, name = "Classification") +
  scale_size_continuous(range = c(2, 7), name = expression(-log[10](padj))) +
  labs(
    title = "Interaction: does Myc\u2019s mitochondrial effect evolve from 6W to 12W?",
    subtitle = paste0(
      "x-axis: fGSEA NES of DESeq2 interaction term (timepoint \u00D7 Myc)\n",
      "y-axis: \u0394 mitoPPS interaction (\u0394Myc_12W \u2212 \u0394Myc_6W)\n",
      "Positive = Myc effect strengthens over time  |  Negative = Myc effect weakens"
    ),
    x = "fGSEA NES (interaction Wald statistic)",
    y = "\u0394 mitoPPS interaction (\u0394Myc_12W \u2212 \u0394Myc_6W)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title      = element_text(face = "bold"),
    legend.position = "bottom",
    legend.text     = element_text(size = 8)
  )

ggsave(file.path(int_dir, "scatter_interaction_fgsea_vs_mitopps.pdf"),
       p_scatter, width = 12, height = 10)

# --- 8b. Dotplot: significant interaction pathways from fGSEA (curated sets) ---
# This shows ALL gene sets (Hallmark + MYC + MitoCarta + Apoptosis), not
# just MitoCarta, to contextualise mitochondrial changes within the broader
# transcriptional landscape.

fgsea_int_sig <- fgsea_interaction_curated %>%
  filter(padj < 0.10) %>%
  mutate(pathway_label = str_replace_all(pathway, "^MSigDB_HALLMARK_|^MC_", "") %>%
           str_replace_all("_", " "))

if (nrow(fgsea_int_sig) > 0) {
  p_curated <- ggplot(fgsea_int_sig,
                      aes(x = NES, y = reorder(pathway_label, NES))) +
    geom_point(aes(size = -log10(padj)), colour = "#B2182B", alpha = 0.8) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    scale_size_continuous(range = c(2, 7), name = expression(-log[10](padj))) +
    labs(
      title = "fGSEA interaction: pathways with evolving Myc effect (6W \u2192 12W)",
      subtitle = paste0(
        "Curated gene sets (Hallmark + MYC + MitoCarta + Apoptosis)  |  padj < 0.10\n",
        "Positive NES: Myc effect strengthens at 12W  |  Negative NES: Myc effect weakens"
      ),
      x = "fGSEA NES (interaction)",
      y = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.y  = element_text(size = 8),
      plot.title   = element_text(face = "bold")
    )

  ggsave(file.path(int_dir, "dotplot_fgsea_interaction_curated.pdf"),
         p_curated, width = 10, height = max(5, nrow(fgsea_int_sig) * 0.3))
} else {
  message("  No curated gene sets significant at padj < 0.10 — skipping dotplot")
}

# --- 8c. Dotplot: significant interaction pathways from fGSEA (MitoCarta) ---

fgsea_int_mc_sig <- fgsea_interaction_mitocarta %>%
  filter(padj < 0.10) %>%
  mutate(
    pathway_label = str_replace_all(pathway, "_", " "),
    tier1 = pathway_tier1_map[pathway]
  )

if (nrow(fgsea_int_mc_sig) > 0) {
  p_mc <- ggplot(fgsea_int_mc_sig,
                 aes(x = NES, y = reorder(pathway_label, NES))) +
    geom_point(aes(size = -log10(padj), colour = tier1), alpha = 0.8) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    scale_colour_manual(values = tier1_colours, na.value = "grey50",
                        name = "MitoCarta Tier 1") +
    scale_size_continuous(range = c(2, 7), name = expression(-log[10](padj))) +
    labs(
      title = "fGSEA interaction: MitoCarta pathways with evolving Myc effect",
      subtitle = paste0(
        "Full MitoCarta gene sets (modified: mtDNA separated, Apoptosis subgroups)\n",
        "padj < 0.10  |  Positive NES: Myc effect strengthens at 12W"
      ),
      x = "fGSEA NES (interaction)",
      y = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.y  = element_text(size = 7),
      plot.title   = element_text(face = "bold"),
      legend.position = "right"
    )

  ggsave(file.path(int_dir, "dotplot_fgsea_interaction_mitocarta.pdf"),
         p_mc, width = 12, height = max(5, nrow(fgsea_int_mc_sig) * 0.3))
} else {
  message("  No MitoCarta pathways significant at padj < 0.10 — skipping dotplot")
}

# --- 8d. Paired dotplot: Myc effect at 6W vs 12W for interaction pathways ---
# For pathways with significant interaction (either method), show the actual
# Δ mitoPPS at each timepoint side by side — makes the interaction visible.

int_pathways <- comparison %>%
  filter(sig_fgsea | sig_anova) %>%
  pull(pathway)

if (length(int_pathways) > 0) {
  paired_data <- pps_interaction %>%
    filter(pathway %in% int_pathways) %>%
    dplyr::select(pathway, diff_6W, diff_12W) %>%
    pivot_longer(cols = c(diff_6W, diff_12W),
                 names_to = "timepoint",
                 values_to = "delta_mitopps") %>%
    mutate(
      timepoint     = ifelse(timepoint == "diff_6W", "6W", "12W"),
      pathway_label = str_replace_all(pathway, "_", " "),
      tier1         = pathway_tier1_map[pathway]
    )

  # Order by interaction magnitude
  pw_order <- pps_interaction %>%
    filter(pathway %in% int_pathways) %>%
    arrange(interaction_diff) %>%
    pull(pathway) %>%
    str_replace_all("_", " ")

  paired_data$pathway_label <- factor(paired_data$pathway_label, levels = pw_order)

  p_paired <- ggplot(paired_data,
                     aes(x = delta_mitopps, y = pathway_label, colour = timepoint)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_line(aes(group = pathway_label), colour = "grey60", linewidth = 0.4) +
    geom_point(size = 3.5, alpha = 0.85) +
    scale_colour_manual(values = c("6W" = "#377EB8", "12W" = "#E41A1C"),
                        name = "Timepoint") +
    labs(
      title = "Evolving Myc effect on mitoPPS: 6W vs 12W",
      subtitle = paste0(
        "Pathways with significant interaction (fGSEA padj < 0.05 or ANOVA padj < 0.10)\n",
        "\u0394 mitoPPS = mean(Myc+) \u2212 mean(Myc\u2212)  |  Connecting lines show interaction"
      ),
      x = "\u0394 mitoPPS (Myc+ \u2212 Myc\u2212)",
      y = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.y  = element_text(size = 8),
      plot.title   = element_text(face = "bold"),
      legend.position = "right"
    )

  ggsave(file.path(int_dir, "dotplot_paired_myc_effect_6W_vs_12W.pdf"),
         p_paired, width = 10, height = max(5, length(int_pathways) * 0.35))
} else {
  message("  No pathways with significant interaction — skipping paired dotplot")
}

# --- 8e. P-value histogram: fGSEA interaction (diagnostic) ---
# If there is pathway-level signal in the interaction, we expect
# enrichment of small p-values (left-skew) beyond the uniform null.

p_pval_hist <- ggplot(fgsea_interaction_mitocarta, aes(x = pval)) +
  geom_histogram(bins = 30, fill = "steelblue", colour = "white", boundary = 0) +
  geom_hline(yintercept = nrow(fgsea_interaction_mitocarta) / 30,
             colour = "red", linetype = "dashed") +
  labs(
    title = "fGSEA interaction p-value distribution (MitoCarta pathways)",
    subtitle = "Red line = expected under null (uniform). Left-skew = real signal.",
    x = "fGSEA p-value",
    y = "Count"
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(int_dir, "histogram_fgsea_interaction_pvalues.pdf"),
       p_pval_hist, width = 8, height = 5)

# --- 8f. Classification summary barplot ---
summary_counts <- comparison %>%
  group_by(category) %>%
  summarise(n = n(), .groups = "drop") %>%
  filter(category != "Neither significant")

if (nrow(summary_counts) > 0) {
  p_bar <- ggplot(summary_counts, aes(x = reorder(category, n), y = n, fill = category)) +
    geom_col(alpha = 0.85) +
    coord_flip() +
    scale_fill_manual(values = category_colours) +
    labs(
      title    = "Interaction classification: fGSEA vs mitoPPS",
      subtitle = "Number of MitoCarta pathways per category",
      x        = NULL,
      y        = "Number of pathways"
    ) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "none")

  ggsave(file.path(int_dir, "barplot_interaction_classification.pdf"),
         p_bar, width = 10, height = 5)
}

# =============================================================================
# PART 9: SUMMARY TABLE
# =============================================================================

message("\n", strrep("=", 70))
message("KEY RESULTS: INTERACTION ANALYSIS")
message(strrep("=", 70))

message("\n--- fGSEA interaction: top MitoCarta pathways (by p-value) ---")
fgsea_interaction_mitocarta %>%
  mutate(pathway_label = str_replace_all(pathway, "_", " "),
         tier1 = pathway_tier1_map[pathway]) %>%
  arrange(pval) %>%
  dplyr::select(pathway_label, tier1, NES, pval, padj, size) %>%
  head(20) %>%
  print()

message("\n--- fGSEA interaction: top curated pathways (by p-value) ---")
fgsea_interaction_curated %>%
  arrange(pval) %>%
  dplyr::select(pathway, NES, pval, padj, size) %>%
  head(20) %>%
  print()

if (any(comparison$sig_fgsea | comparison$sig_anova, na.rm = TRUE)) {
  message("\n--- Pathways with significant interaction (either method) ---")
  comparison %>%
    filter(sig_fgsea | sig_anova) %>%
    arrange(pmin(padj_fgsea, anova_padj)) %>%
    dplyr::select(pathway_label, tier1, NES, padj_fgsea,
                  diff_6W, diff_12W, interaction_diff, anova_padj,
                  category, myc_effect_pattern) %>%
    as_tibble() |>
    print(n = 40)
} else {
  message("\n  No pathways reached significance in either method.")
  message("  This is consistent with the Myc reprioritisation being")
  message("  established early (by 6W) and largely stable through 12W.")
}

# =============================================================================
# PART 10: SAVE RESULTS
# =============================================================================

message("\nSaving results...")

interaction_analysis <- list(
  # fGSEA results
  fgsea_interaction_curated   = fgsea_interaction_curated,
  fgsea_interaction_mitocarta = fgsea_interaction_mitocarta,

  # mitoPPS interaction
  pps_interaction = pps_interaction,

  # Combined comparison
  comparison = comparison,

  # Correlation
  correlation = list(
    r       = cor_result$estimate,
    p_value = cor_result$p.value,
    method  = "Pearson"
  ),

  # Ranking
  ranks_interaction = ranks_interaction,

  # Metadata
  analysis_date = Sys.Date(),
  description = paste(
    "Interaction analysis: does the Myc effect on mitochondrial pathways",
    "evolve between 6W and 12W?",
    "fGSEA on DESeq2 interaction Wald statistics (timepoint12W:myc_statuspos).",
    "mitoPPS interaction = delta_Myc_12W - delta_Myc_6W.",
    "Curated gene sets: Hallmark + MYC + MitoCarta + Apoptosis.",
    "MitoCarta sets: modified with mtDNA separation and Apoptosis subgroups."
  )
)

saveRDS(interaction_analysis, here("results", "interaction_fgsea_mitopps.rds"))

# CSV export
comparison %>%
  dplyr::select(pathway, pathway_label, tier1,
                NES, padj_fgsea, sig_fgsea,
                diff_6W, diff_12W, interaction_diff,
                anova_F, anova_p, anova_padj, sig_anova,
                category, myc_effect_pattern) %>%
  arrange(pmin(padj_fgsea, anova_padj, na.rm = TRUE)) %>%
  write.csv(file.path(int_dir, "interaction_fgsea_vs_mitopps.csv"),
            row.names = FALSE)

fgsea_interaction_curated %>%
  arrange(pval) %>%
  mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ";")) %>%
  write.csv(file.path(int_dir, "fgsea_interaction_curated_all.csv"),
            row.names = FALSE)

message("\n", strrep("=", 70))
message("INTERACTION ANALYSIS COMPLETE")
message(strrep("=", 70))
message(sprintf("Results saved to: results/interaction_fgsea_mitopps.rds"))
message(sprintf("Figures saved to: %s", int_dir))
message(strrep("=", 70))
