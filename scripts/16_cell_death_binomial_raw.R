# scripts/16_cell_death_binomial_raw.R
# =============================================================================
# Cell Death directional binomial test (branch 1) -- RAW / unshrunken LFCs
# Ported from archive_main_pipeline/09_cell_death_pathway_summary.R
# =============================================================================
#
# PURPOSE: Test whether cell death gene expression supports the observation
#   that MYC causes more apoptosis at 6W than 12W. Directional-hypothesis lens
#   (the complement of the broad-enrichment branch 2, script 12). Both branches
#   go into the final analysis; report convergence and divergence (CLAUDE.md
#   "Cell death -- both branches, neither canonical").
#
# HYPOTHESIS:
#   - Pro-death genes are induced MORE by MYC at 6W than 12W
#   - Pro-survival genes are induced MORE by MYC at 12W than 6W
#
# WHY THIS PORT EXISTS (the single analytical change vs archive 09):
#   The metric is a DELTA-LFC difference (myc_6W_log2FC - myc_12W_log2FC).
#   Per CLAUDE.md "Shrunken vs raw LFCs", delta-LFC methods MUST run on RAW
#   (unshrunken MLE) LFCs -- shrinkage distorts the difference. The archived
#   script read the SHRUNKEN combined_df_annotated.rds; this port reads
#   combined_df_annotated_raw.rds instead. Everything else (binomial test,
#   category breakdown, plots) is faithful to the original.
#
# COLUMN NOTE (verified 2026-07-06): the raw df carries suffixed columns
#   myc_6W_log2FC_raw / myc_12W_log2FC_raw and has NO `gene` column. We alias
#   the two LFC columns to the archived names right after the read so the ported
#   logic below is unchanged, and drop `gene` from the final select.
#
# Outputs are all suffixed _raw so they never clobber the shrunken run.
#
# Input:
#   - results/combined_df_annotated_raw.rds     (RAW LFCs; keys mgi_symbol)
#   - data/cell_death_genes_consolidated.csv    (effect / pathway annotations)
#
# Output:
#   - results/cell_death_de_raw.rds
#   - results/cell_death_hypothesis_results_raw.csv
#   - results/cell_death_genes_full_raw.csv
#   - results/cell_death_supporting_genes_raw.csv
#   - results/cell_death_opposing_genes_raw.csv
#   - outputs/cell_death_raw/hypothesis_all_categories.pdf
#   - outputs/cell_death_raw/cell_death_hypothesis_scatter.pdf
# =============================================================================

source(here::here("scripts", "00_setup_packages.R"))

# =============================================================================
# PART 1: LOAD DATA  (RAW file; alias the _raw LFC columns)
# =============================================================================

combined_df <- readRDS(here::here("results", "combined_df_annotated_raw.rds"))
cell_death_genes <- readr::read_csv(
  here::here("data", "cell_death_genes_consolidated.csv"), show_col_types = FALSE
)

# Alias raw-suffixed columns to the names the ported logic expects.
combined_df <- combined_df |>
  dplyr::mutate(
    myc_6W_log2FC  = myc_6W_log2FC_raw,
    myc_12W_log2FC = myc_12W_log2FC_raw
  )

message(sprintf("Loaded %d genes from RAW DE results", nrow(combined_df)))
message(sprintf("Loaded %d cell death genes", nrow(cell_death_genes)))

# =============================================================================
# PART 2: JOIN CELL DEATH GENES WITH DE RESULTS
# =============================================================================

cell_death_de <- combined_df |>
  dplyr::inner_join(
    cell_death_genes |>
      dplyr::select(mouse_symbol, effect, pathway, is_core, is_mitochondrial,
                    in_GO, in_KEGG, in_Reactome, in_Hallmark),
    by = c("mgi_symbol" = "mouse_symbol")
  )

# Remove duplicates (keep unique gene-effect combinations)
cell_death_de <- cell_death_de |>
  dplyr::distinct(mgi_symbol, .keep_all = TRUE)

message(sprintf("Matched %d unique cell death genes in DE results", nrow(cell_death_de)))
message(sprintf("  Pro-death: %d", sum(cell_death_de$effect == "pro-death", na.rm = TRUE)))
message(sprintf("  Pro-survival: %d", sum(cell_death_de$effect == "pro-survival", na.rm = TRUE)))

saveRDS(cell_death_de, here::here("results", "cell_death_de_raw.rds"))

# =============================================================================
# PART 3: CALCULATE HYPOTHESIS SCORES
# =============================================================================
# For pro-death genes:    positive score if MYC effect stronger at 6W
# For pro-survival genes:  positive score if MYC effect stronger at 12W

cell_death_full <- cell_death_de |>
  dplyr::filter(effect %in% c("pro-death", "pro-survival")) |>
  dplyr::mutate(
    delta_myc_effect = myc_6W_log2FC - myc_12W_log2FC,
    hypothesis_score = dplyr::case_when(
      effect == "pro-death"    ~ delta_myc_effect,   # higher at 6W supports
      effect == "pro-survival" ~ -delta_myc_effect,  # higher at 12W supports
      TRUE ~ NA_real_
    ),
    hypothesis_direction = dplyr::case_when(
      hypothesis_score > 0.2  ~ "SUPPORTING",
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

test_hypothesis_subset <- function(df, subset_name) {
  supporting <- sum(df$hypothesis_direction == "SUPPORTING", na.rm = TRUE)
  opposing   <- sum(df$hypothesis_direction == "OPPOSING", na.rm = TRUE)
  total      <- supporting + opposing

  if (total == 0) {
    return(tibble::tibble(
      category = subset_name, n_genes = nrow(df),
      n_supporting = 0, n_opposing = 0,
      pct_supporting = NA_real_, binom_p = NA_real_
    ))
  }

  binom_result <- binom.test(supporting, total, p = 0.5)

  tibble::tibble(
    category = subset_name,
    n_genes = total,
    n_supporting = supporting,
    n_opposing = opposing,
    pct_supporting = round(100 * supporting / total, 1),
    binom_p = round(binom_result$p.value, 4)
  )
}

results <- dplyr::bind_rows(
  test_hypothesis_subset(cell_death_full, "All genes"),
  test_hypothesis_subset(cell_death_full |> dplyr::filter(pathway == "apoptosis"), "Apoptosis pathway"),
  test_hypothesis_subset(cell_death_full |> dplyr::filter(pathway %in% c("CICD", "both")), "CICD pathway"),
  test_hypothesis_subset(cell_death_full |> dplyr::filter(is_core), "Core genes"),
  test_hypothesis_subset(cell_death_full |> dplyr::filter(!is_core), "Non-core genes"),
  test_hypothesis_subset(cell_death_full |> dplyr::filter(is_mitochondrial), "Mitochondrial"),
  test_hypothesis_subset(cell_death_full |> dplyr::filter(!is_mitochondrial), "Non-mitochondrial"),
  test_hypothesis_subset(cell_death_full |> dplyr::filter(in_GO), "GO database"),
  test_hypothesis_subset(cell_death_full |> dplyr::filter(in_KEGG), "KEGG database"),
  test_hypothesis_subset(cell_death_full |> dplyr::filter(in_Reactome), "Reactome database"),
  test_hypothesis_subset(cell_death_full |> dplyr::filter(in_Hallmark), "Hallmark (MSigDB)")
)

results <- results |>
  dplyr::mutate(
    direction = dplyr::case_when(
      binom_p < 0.01 & n_supporting > n_opposing ~ "SUPPORTS ***",
      binom_p < 0.1  & n_supporting > n_opposing ~ "SUPPORTS *",
      binom_p < 0.01 & n_opposing > n_supporting ~ "OPPOSES ***",
      binom_p < 0.1  & n_opposing > n_supporting ~ "OPPOSES *",
      TRUE ~ "No bias"
    )
  )

cat("\n")
cat("================================================================================\n")
cat("     HYPOTHESIS TEST (RAW LFCs): MYC causes more apoptosis at 6W than 12W\n")
cat("================================================================================\n\n")
print(results)

# =============================================================================
# PART 5: IDENTIFY KEY GENES
# =============================================================================

supporting_genes <- cell_death_full |>
  dplyr::filter(hypothesis_direction == "SUPPORTING") |>
  dplyr::arrange(dplyr::desc(hypothesis_score)) |>
  dplyr::select(mgi_symbol, effect, pathway, is_core, is_mitochondrial,
                myc_6W_log2FC, myc_12W_log2FC, hypothesis_score)

opposing_genes <- cell_death_full |>
  dplyr::filter(hypothesis_direction == "OPPOSING") |>
  dplyr::arrange(hypothesis_score) |>
  dplyr::select(mgi_symbol, effect, pathway, is_core, is_mitochondrial,
                myc_6W_log2FC, myc_12W_log2FC, hypothesis_score)

cat("\n=== TOP 15 GENES SUPPORTING HYPOTHESIS ===\n\n")
supporting_genes |> head(15) |> print()

cat("\n=== TOP 15 GENES OPPOSING HYPOTHESIS ===\n\n")
opposing_genes |> head(15) |> print()

# =============================================================================
# PART 6: VISUALIZATIONS
# =============================================================================

fig_dir <- here::here("outputs", "cell_death_raw")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# --- Plot 1: Category comparison bar plot ---
plot_results <- results |>
  dplyr::mutate(
    net_support = n_supporting - n_opposing,
    category = factor(category, levels = category[order(binom_p)])
  )

p_categories <- ggplot2::ggplot(plot_results,
    ggplot2::aes(x = reorder(category, -binom_p), y = net_support)) +
  ggplot2::geom_col(ggplot2::aes(fill = direction), width = 0.7) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey30") +
  ggplot2::geom_text(ggplot2::aes(label = paste0(n_supporting, "/", n_supporting + n_opposing)),
            vjust = ifelse(plot_results$net_support >= 0, -0.5, 1.5), size = 3) +
  ggplot2::scale_fill_manual(
    values = c("SUPPORTS *" = "#2E7D32", "SUPPORTS ***" = "#1B5E20",
               "OPPOSES *" = "#C62828", "OPPOSES ***" = "#B71C1C",
               "No bias" = "grey60"),
    name = "Direction"
  ) +
  ggplot2::coord_flip() +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::labs(
    x = NULL,
    y = "Net support (Supporting - Opposing genes)",
    title = "Hypothesis Test Across Gene Categories (RAW LFCs)",
    subtitle = "Hypothesis: MYC causes more apoptosis at 6W than 12W\nNumbers show supporting/total genes with |delta LFC| > 0.2"
  )

ggplot2::ggsave(file.path(fig_dir, "hypothesis_all_categories.pdf"), p_categories,
       width = 10, height = 6)

# --- Plot 2: Scatter plot of MYC effects ---
p_scatter <- ggplot2::ggplot(cell_death_full,
    ggplot2::aes(x = myc_6W_log2FC, y = myc_12W_log2FC)) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  ggplot2::geom_point(ggplot2::aes(color = effect, shape = hypothesis_direction),
             alpha = 0.7, size = 2) +
  ggplot2::scale_color_manual(
    values = c("pro-death" = "#D73027", "pro-survival" = "#4575B4"),
    name = "Effect"
  ) +
  ggplot2::scale_shape_manual(
    values = c("SUPPORTING" = 17, "OPPOSING" = 25, "NEUTRAL" = 16),
    name = "Hypothesis"
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::labs(
    x = "MYC effect at 6W (log2FC, raw)",
    y = "MYC effect at 12W (log2FC, raw)",
    title = "Cell Death Genes: MYC Effect Comparison (RAW LFCs)",
    subtitle = "Points above diagonal: MYC effect stronger at 12W"
  )

label_genes <- cell_death_full |>
  dplyr::filter(is_core | abs(hypothesis_score) > 0.4) |>
  dplyr::filter(hypothesis_direction != "NEUTRAL")

p_scatter_labeled <- p_scatter +
  ggrepel::geom_text_repel(
    data = label_genes,
    ggplot2::aes(label = mgi_symbol),
    size = 2.5, max.overlaps = 20,
    segment.color = "grey60", segment.size = 0.3
  )

ggplot2::ggsave(file.path(fig_dir, "cell_death_hypothesis_scatter.pdf"), p_scatter_labeled,
       width = 10, height = 8)

# =============================================================================
# PART 7: SAVE ALL RESULTS
# =============================================================================
# NOTE: dropped `gene` from the full-table select -- the raw df has no `gene`
# column (mgi_symbol is the gene identifier).

readr::write_csv(results, here::here("results", "cell_death_hypothesis_results_raw.csv"))

readr::write_csv(
  cell_death_full |>
    dplyr::select(mgi_symbol, effect, pathway, is_core, is_mitochondrial,
                  in_GO, in_KEGG, in_Reactome, in_Hallmark,
                  myc_6W_log2FC, myc_12W_log2FC, delta_myc_effect,
                  hypothesis_score, hypothesis_direction),
  here::here("results", "cell_death_genes_full_raw.csv")
)

readr::write_csv(supporting_genes, here::here("results", "cell_death_supporting_genes_raw.csv"))
readr::write_csv(opposing_genes, here::here("results", "cell_death_opposing_genes_raw.csv"))

message("\n=============================================================================")
message("Cell death analysis (RAW) complete.")
message("Results saved to results/ (suffix _raw); figures to outputs/cell_death_raw/")
message("=============================================================================")

# =============================================================================
# SANDBOX (skipped by source(); run line-by-line in Positron)
# =============================================================================
if (FALSE) {

  # FIRST LINE CHECK -- confirm the raw df columns before trusting the join.
  cdr <- readRDS(here::here("results", "combined_df_annotated_raw.rds"))
  head(cdr) |> print()
  colnames(cdr)
  stopifnot(all(c("mgi_symbol", "myc_6W_log2FC_raw", "myc_12W_log2FC_raw") %in% colnames(cdr)))

  # Compare the RAW binomial output against the SHRUNKEN run (archive 09).
  # Expect: same direction of bias, magnitude typically sharper on raw.
  raw_res <- readr::read_csv(here::here("results", "cell_death_hypothesis_results_raw.csv"),
                             show_col_types = FALSE)
  shrunk_path <- here::here("results", "cell_death_hypothesis_results.csv")
  if (file.exists(shrunk_path)) {
    shrunk_res <- readr::read_csv(shrunk_path, show_col_types = FALSE)
    dplyr::full_join(
      raw_res |> dplyr::select(category, raw_pct = pct_supporting, raw_p = binom_p),
      shrunk_res |> dplyr::select(category, shrunk_pct = pct_supporting, shrunk_p = binom_p),
      by = "category"
    ) |> print(n = Inf)
  }

  # Cross-check with Gate 2 (script 14): PRO module coordinately shifted negative.
  # Here the pro-death SUPPORTING bias is the delta-LFC analogue of that shift.
  results |> dplyr::filter(category %in% c("All genes", "Apoptosis pathway")) |> print()

  print(p_categories)
  print(p_scatter_labeled)
}
