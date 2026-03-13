# scripts/09_cell_death_analysis.R
# =============================================================================
# Cell Death Gene Analysis for MYC Timecourse Study
# =============================================================================
# 
# PURPOSE: Test whether cell death gene expression supports the experimental
#          observation that MYC causes more apoptosis at 6W than 12W
#
# HYPOTHESIS: 
#   - Pro-death genes are induced MORE by MYC at 6W than 12W
#   - Pro-survival genes are induced MORE by MYC at 12W than 6W
#
# OUTPUT:
#   - results/cell_death_hypothesis_results.csv (summary by category)
#   - results/cell_death_genes_full.csv (all genes with scores)
#   - results/cell_death_supporting_genes.csv
#   - results/cell_death_opposing_genes.csv
#   - results/figures/hypothesis_all_categories.pdf
#   - results/figures/cell_death_hypothesis_scatter.pdf
#
# =============================================================================

source("scripts/00_setup_packages.R")
library(tidyr)

# =============================================================================
# PART 1: LOAD DATA
# =============================================================================

combined_df <- readRDS("results/combined_df_annotated.rds")
cell_death_genes <- read_csv("data/cell_death_genes_consolidated.csv", show_col_types = FALSE)

message(sprintf("Loaded %d genes from DE results", nrow(combined_df)))
message(sprintf("Loaded %d cell death genes", nrow(cell_death_genes)))

# =============================================================================
# PART 2: JOIN CELL DEATH GENES WITH DE RESULTS
# =============================================================================

cell_death_de <- combined_df |>
  inner_join(
    cell_death_genes |> 
      dplyr::select(mouse_symbol, effect, pathway, is_core, is_mitochondrial,
                    in_GO, in_KEGG, in_Reactome, in_Hallmark),
    by = c("mgi_symbol" = "mouse_symbol")
  )

# Remove duplicates (keep unique gene-effect combinations)
cell_death_de <- cell_death_de |>
  distinct(mgi_symbol, .keep_all = TRUE)

message(sprintf("Matched %d unique cell death genes in DE results", nrow(cell_death_de)))
message(sprintf("  Pro-death: %d", sum(cell_death_de$effect == "pro-death", na.rm = TRUE)))
message(sprintf("  Pro-survival: %d", sum(cell_death_de$effect == "pro-survival", na.rm = TRUE)))

# Save joined data
saveRDS(cell_death_de, "results/cell_death_de.rds")

# =============================================================================
# PART 3: CALCULATE HYPOTHESIS SCORES
# =============================================================================

# Hypothesis score calculation:
# For pro-death genes: positive score if MYC effect stronger at 6W
# For pro-survival genes: positive score if MYC effect stronger at 12W (protective)

cell_death_full <- cell_death_de |>
  filter(effect %in% c("pro-death", "pro-survival")) |>
  mutate(
    # Raw difference in MYC effect
    delta_myc_effect = myc_6W_log2FC - myc_12W_log2FC,
    
    # Hypothesis score (aligned so positive = supports hypothesis)
    hypothesis_score = case_when(
      effect == "pro-death" ~ delta_myc_effect,      # Pro-death: higher at 6W supports
      effect == "pro-survival" ~ -delta_myc_effect,  # Pro-survival: higher at 12W supports
      TRUE ~ NA_real_
    ),
    
    # Direction classification (threshold = 0.2 LFC difference)
    hypothesis_direction = case_when(
      hypothesis_score > 0.2 ~ "SUPPORTING",
      hypothesis_score < -0.2 ~ "OPPOSING",
      TRUE ~ "NEUTRAL"
    )
  )

message(sprintf("\nGenes with hypothesis scores: %d", nrow(cell_death_full)))
message(sprintf("  Supporting: %d", sum(cell_death_full$hypothesis_direction == "SUPPORTING")))
message(sprintf("  Opposing: %d", sum(cell_death_full$hypothesis_direction == "OPPOSING")))
message(sprintf("  Neutral: %d", sum(cell_death_full$hypothesis_direction == "NEUTRAL")))

# =============================================================================
# PART 4: TEST HYPOTHESIS ACROSS CATEGORIES
# =============================================================================

# Function to test hypothesis for a subset of genes
test_hypothesis_subset <- function(df, subset_name) {
  supporting <- sum(df$hypothesis_direction == "SUPPORTING", na.rm = TRUE)
  opposing <- sum(df$hypothesis_direction == "OPPOSING", na.rm = TRUE)
  total <- supporting + opposing
  
  if (total == 0) {
    return(tibble(
      category = subset_name,
      n_genes = nrow(df),
      n_supporting = 0,
      n_opposing = 0,
      pct_supporting = NA_real_,
      binom_p = NA_real_
    ))
  }
  
  # Binomial test: is there a bias toward supporting or opposing?
  binom_result <- binom.test(supporting, total, p = 0.5)
  
  tibble(
    category = subset_name,
    n_genes = total,
    n_supporting = supporting,
    n_opposing = opposing,
    pct_supporting = round(100 * supporting / total, 1),
    binom_p = round(binom_result$p.value, 4)
  )
}

# Test across all categories
results <- bind_rows(
  # Overall
  test_hypothesis_subset(cell_death_full, "All genes"),
  
  # By pathway
  test_hypothesis_subset(cell_death_full |> filter(pathway == "apoptosis"), "Apoptosis pathway"),
  test_hypothesis_subset(cell_death_full |> filter(pathway %in% c("CICD", "both")), "CICD pathway"),
  
  # By core status
  test_hypothesis_subset(cell_death_full |> filter(is_core), "Core genes"),
  test_hypothesis_subset(cell_death_full |> filter(!is_core), "Non-core genes"),
  
  # By mitochondrial status
  test_hypothesis_subset(cell_death_full |> filter(is_mitochondrial), "Mitochondrial"),
  test_hypothesis_subset(cell_death_full |> filter(!is_mitochondrial), "Non-mitochondrial"),
  
  # By database source
  test_hypothesis_subset(cell_death_full |> filter(in_GO), "GO database"),
  test_hypothesis_subset(cell_death_full |> filter(in_KEGG), "KEGG database"),
  test_hypothesis_subset(cell_death_full |> filter(in_Reactome), "Reactome database"),
  test_hypothesis_subset(cell_death_full |> filter(in_Hallmark), "Hallmark (MSigDB)")
)

# Add direction interpretation
results <- results |>
  mutate(
    direction = case_when(
      binom_p < 0.01 & n_supporting > n_opposing ~ "SUPPORTS ***",
      binom_p < 0.1 & n_supporting > n_opposing ~ "SUPPORTS *",
      binom_p < 0.01 & n_opposing > n_supporting ~ "OPPOSES ***",
      binom_p < 0.1 & n_opposing > n_supporting ~ "OPPOSES *",
      TRUE ~ "No bias"
    )
  )

# Print results
cat("\n")
cat("================================================================================\n")
cat("         HYPOTHESIS TEST: MYC causes more apoptosis at 6W than 12W\n")
cat("================================================================================\n\n")
print(results)

# =============================================================================
# PART 5: IDENTIFY KEY GENES
# =============================================================================

# Top genes supporting hypothesis
supporting_genes <- cell_death_full |>
  filter(hypothesis_direction == "SUPPORTING") |>
  arrange(desc(hypothesis_score)) |>
  dplyr::select(mgi_symbol, effect, pathway, is_core, is_mitochondrial,
                myc_6W_log2FC, myc_12W_log2FC, hypothesis_score)

# Top genes opposing hypothesis  
opposing_genes <- cell_death_full |>
  filter(hypothesis_direction == "OPPOSING") |>
  arrange(hypothesis_score) |>
  dplyr::select(mgi_symbol, effect, pathway, is_core, is_mitochondrial,
                myc_6W_log2FC, myc_12W_log2FC, hypothesis_score)

cat("\n=== TOP 15 GENES SUPPORTING HYPOTHESIS ===\n\n")
print(head(supporting_genes, 15))

cat("\n=== TOP 15 GENES OPPOSING HYPOTHESIS ===\n\n")
print(head(opposing_genes, 15))

# =============================================================================
# PART 6: VISUALIZATIONS
# =============================================================================

dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)

# --- Plot 1: Category comparison bar plot ---
plot_results <- results |>
  mutate(
    net_support = n_supporting - n_opposing,
    category = factor(category, levels = category[order(binom_p)])
  )

p_categories <- ggplot(plot_results, aes(x = reorder(category, -binom_p), y = net_support)) +
  geom_col(aes(fill = direction), width = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey30") +
  geom_text(aes(label = paste0(n_supporting, "/", n_supporting + n_opposing)), 
            vjust = ifelse(plot_results$net_support >= 0, -0.5, 1.5), size = 3) +
  scale_fill_manual(
    values = c("SUPPORTS *" = "#2E7D32", "SUPPORTS ***" = "#1B5E20", 
               "OPPOSES *" = "#C62828", "OPPOSES ***" = "#B71C1C",
               "No bias" = "grey60"),
    name = "Direction"
  ) +
  coord_flip() +
  theme_minimal(base_size = 11) +
  labs(
    x = NULL,
    y = "Net support (Supporting - Opposing genes)",
    title = "Hypothesis Test Across Gene Categories",
    subtitle = "Hypothesis: MYC causes more apoptosis at 6W than 12W\nNumbers show supporting/total genes with |delta LFC| > 0.2"
  )

ggsave("results/figures/hypothesis_all_categories.pdf", p_categories, 
       width = 10, height = 6)

# --- Plot 2: Scatter plot of MYC effects ---
p_scatter <- ggplot(cell_death_full, aes(x = myc_6W_log2FC, y = myc_12W_log2FC)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(aes(color = effect, shape = hypothesis_direction), alpha = 0.7, size = 2) +
  scale_color_manual(
    values = c("pro-death" = "#D73027", "pro-survival" = "#4575B4"),
    name = "Effect"
  ) +
  scale_shape_manual(
    values = c("SUPPORTING" = 17, "OPPOSING" = 25, "NEUTRAL" = 16),
    name = "Hypothesis"
  ) +
  theme_minimal(base_size = 11) +
  labs(
    x = "MYC effect at 6W (log2FC)",
    y = "MYC effect at 12W (log2FC)",
    title = "Cell Death Genes: MYC Effect Comparison",
    subtitle = "Points above diagonal: MYC effect stronger at 12W"
  )

# Add labels for key genes
label_genes <- cell_death_full |>
  filter(is_core | abs(hypothesis_score) > 0.4) |>
  filter(hypothesis_direction != "NEUTRAL")

p_scatter_labeled <- p_scatter +
  ggrepel::geom_text_repel(
    data = label_genes,
    aes(label = mgi_symbol),
    size = 2.5,
    max.overlaps = 20,
    segment.color = "grey60",
    segment.size = 0.3
  )

ggsave("results/figures/cell_death_hypothesis_scatter.pdf", p_scatter_labeled, 
       width = 10, height = 8)

# =============================================================================
# PART 7: SAVE ALL RESULTS
# =============================================================================

# Save summary results
write_csv(results, "results/cell_death_hypothesis_results.csv")

# Save full gene table
write_csv(
  cell_death_full |>
    dplyr::select(mgi_symbol, gene, effect, pathway, is_core, is_mitochondrial,
                  in_GO, in_KEGG, in_Reactome, in_Hallmark,
                  myc_6W_log2FC, myc_12W_log2FC, delta_myc_effect,
                  hypothesis_score, hypothesis_direction),
  "results/cell_death_genes_full.csv"
)

# Save supporting and opposing gene lists
write_csv(supporting_genes, "results/cell_death_supporting_genes.csv")
write_csv(opposing_genes, "results/cell_death_opposing_genes.csv")

message("\n=============================================================================")
message("Cell death analysis complete.")
message("Results saved to results/")
message("Figures saved to results/figures/")
message("=============================================================================")