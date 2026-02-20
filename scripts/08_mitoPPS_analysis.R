# scripts/08_mitoPPS_analysis.R
# =============================================================================
# MitoPPS Analysis: Mitochondrial Pathway Prioritization Scores
# =============================================================================
#
# Based on: Monzel et al. (2025) "A Quantitative Approach to Mapping
#   Mitochondrial Specialization and Plasticity"
#   doi: 10.1101/2025.02.03.635951
#   Code: https://github.com/annamonzel/mitotyping (branch: new_code)
#
# Purpose:
#   1. Compute per-sample MitoPathway scores (raw) for all MitoCarta3.0
#      pathways — reflects absolute mitochondrial pathway expression
#   2. Compute per-sample mitoPPS via pairwise ratio normalisation —
#      reflects which pathways are *prioritised* independent of total
#      mito content and intrinsic pathway-scale differences
#   3. Statistical comparison of both metrics across four experimental groups
#
# mitoPPS algorithm (from Monzel et al. source code):
#   Step 1: For each sample, compute ratio_ij = score_i / score_j for all
#           pathway pairs (i ≠ j)
#   Step 2: For each pair (i,j), compute average_ratio across all samples,
#           then corrected_ij = ratio_ij / average_ratio_ij
#   Step 3: For each sample and pathway i,
#           mitoPPS_i = mean(corrected_ij) across all j ≠ i
#
#   Interpretation:
#     mitoPPS ≈ 1.0 → pathway prioritised at dataset average
#     mitoPPS > 1.0 → pathway selectively up-prioritised
#     mitoPPS < 1.0 → pathway selectively down-prioritised
#
# Input:
#   - DESeq2 object (dds_int_run.rds) — normalised counts extracted directly
#   - MitoCarta3.0 mouse annotation file (Mouse.MitoCarta3.0.xls), Sheet 4
#     Download from: https://www.broadinstitute.org/mitocarta
#   - Sample metadata (coldata)
#
# Output:
#   - results/mitopps_scores.rds       (all scores and metadata)
#   - outputs/mitopps/                  (visualisations)
#
# =============================================================================

source(here::here("scripts", "00_setup_packages.R"))

# Additional packages for this script
if (!requireNamespace("readxl", quietly = TRUE)) install.packages("readxl")
if (!requireNamespace("splitstackshape", quietly = TRUE)) install.packages("splitstackshape")
library(readxl)
library(splitstackshape)

# === Create output directories ===
mitopps_fig_dir <- here("outputs", "mitopps")
dir.create(mitopps_fig_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# PART 1: LOAD DATA
# =============================================================================

message("Loading data...")

# Load DESeq2 object — extract SIZE-FACTOR NORMALISED COUNTS (linear scale)
# This is the correct input for mitoPPS: linear-scale expression values,
# analogous to CPM. The Monzel et al. ROSMAP code back-transforms log2CPM
# with 2^value before computing pathway means (mitoPPS_RM.R, line 49).
# Using normalised counts directly avoids any log-space approximation.
dds_int <- readRDS(here("results", "dds_int_run.rds"))
expr_matrix <- counts(dds_int, normalized = TRUE)  # genes (Ensembl IDs) × samples

# Load metadata
coldata <- readRDS(here("results", "coldata.rds"))

# Load ortholog table for Ensembl → gene symbol mapping
ortholog_table <- readRDS(here("results", "ortholog_table.rds"))

# Create Ensembl → mouse gene symbol map (one-to-one where possible)
ensembl_to_symbol <- ortholog_table %>%
  dplyr::select(ensembl_gene_id, external_gene_name) %>%
  distinct() %>%
  filter(external_gene_name != "") %>%
  # If multiple symbols per Ensembl, keep first
  group_by(ensembl_gene_id) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  deframe()

# =============================================================================
# PART 2: LOAD MitoCarta 3.0 PATHWAY ANNOTATIONS (SHEET 4)
# =============================================================================
# Following Monzel et al. source code exactly:
#   Sheet 4 contains columns: MitoPathway (pathway name), Genes (comma-separated)
#   The gene-to-pathway mapping uses splitstackshape::cSplit to expand the
#   comma-separated gene lists.
#   Pathway hierarchy is in a separate column "MitoPathway Hierarchy" (mouse)
#   parsed with separate(..., sep = " > ") for annotation/colouring only.

message("Loading MitoCarta 3.0 annotations (Sheet 4)...")

mitocarta_path <- here("data", "Mouse.MitoCarta3.0.xls")

if (!file.exists(mitocarta_path)) {
  stop(
    "MitoCarta3.0 mouse file not found at:\n  ", mitocarta_path, "\n",
    "Download from: https://www.broadinstitute.org/mitocarta\n",
    "Place Mouse.MitoCarta3.0.xls in the data/ directory."
  )
}

# --- Gene-to-pathway mapping from Sheet 4 ---
# Exactly as in mitoPPS_RM.R lines 30–39, mitoPPS_MBM.R lines 5–14,
# MitoPPS_calculation.R lines 50–59
mitocarta_sheet4 <- read_xls(mitocarta_path, sheet = 4) %>%
  dplyr::select(MitoPathway, Genes) %>%
  na.omit()

gene_to_pathway <- cSplit(mitocarta_sheet4, "Genes", ",") %>%
  column_to_rownames("MitoPathway") %>%
  t() %>%
  as.data.frame()
gene_to_pathway <- gene_to_pathway %>%
  pivot_longer(cols = colnames(gene_to_pathway),
               names_to = "Pathway", values_to = "Gene") %>%
  na.omit() %>%
  mutate(Gene = as.character(Gene))  # ensure character type

n_pathways_total <- length(unique(gene_to_pathway$Pathway))
n_genes_mc <- length(unique(gene_to_pathway$Gene))
message(sprintf("MitoCarta3.0 Sheet 4: %d genes mapped to %d pathways",
                n_genes_mc, n_pathways_total))

# --- Pathway hierarchy from Sheet 4 (for annotation/colouring) ---
# As in MitoPPS_calculation.R lines 64–73
# Mouse file uses "MitoPathway Hierarchy" (not "MitoPathways Hierarchy")
sheet4_cols <- colnames(read_xls(mitocarta_path, sheet = 4))
hierarchy_col <- intersect(
  c("MitoPathway Hierarchy", "MitoPathways Hierarchy"),
  sheet4_cols
)
if (length(hierarchy_col) == 0) {
  warning("No hierarchy column found in Sheet 4. Tier1 annotations unavailable.")
  pathway_levels <- data.frame(
    Pathway = unique(gene_to_pathway$Pathway),
    Pathway_Level1 = NA_character_,
    Pathway_Level2 = NA_character_,
    Pathway_Level3 = NA_character_,
    Level = NA_character_
  )
} else {
  pathway_levels <- read_xls(mitocarta_path, sheet = 4) %>%
    dplyr::select(MitoPathway, all_of(hierarchy_col[1])) %>%
    dplyr::rename(Hierarchy = 2) %>%
    separate(Hierarchy,
             into = c("Pathway_Level1", "Pathway_Level2", "Pathway_Level3"),
             sep = " > ", fill = "right") %>%
    mutate(Level = case_when(
      MitoPathway == Pathway_Level1 ~ "Pathway_Level1",
      MitoPathway == Pathway_Level2 ~ "Pathway_Level2",
      MitoPathway == Pathway_Level3 ~ "Pathway_Level3"
    )) %>%
    unique() %>%
    dplyr::rename(Pathway = MitoPathway)
}

# Create tier1 lookup: pathway name → top-level category
pathway_tier1_map <- pathway_levels %>%
  dplyr::select(Pathway, Pathway_Level1) %>%
  distinct() %>%
  group_by(Pathway) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  deframe()

# =============================================================================
# PART 3: MAP EXPRESSION MATRIX TO GENE SYMBOLS
# =============================================================================

message("Mapping expression matrix to gene symbols...")

# Convert Ensembl-indexed expression matrix to gene symbol-indexed
mapped_genes <- intersect(rownames(expr_matrix), names(ensembl_to_symbol))
expr_symbols <- expr_matrix[mapped_genes, ]
rownames(expr_symbols) <- ensembl_to_symbol[mapped_genes]

# Handle duplicate symbols: keep the one with highest mean expression
mean_expr <- rowMeans(expr_symbols)
expr_symbols <- expr_symbols[order(-mean_expr), ]
expr_symbols <- expr_symbols[!duplicated(rownames(expr_symbols)), ]

# Check overlap with MitoCarta
mc_symbols <- unique(gene_to_pathway$Gene)
overlap <- intersect(rownames(expr_symbols), mc_symbols)
message(sprintf("Gene overlap: %d / %d MitoCarta genes found in expression data (%.1f%%)",
                length(overlap), length(mc_symbols),
                100 * length(overlap) / length(mc_symbols)))

# =============================================================================
# PART 4: COMPUTE RAW MitoPathway SCORES
# =============================================================================
# For each pathway and each sample: mean NORMALISED COUNT expression of
# member genes. This reflects ABSOLUTE pathway activity (confounded by mito
# content). Uses linear-scale values as per Monzel et al.

message("Computing raw MitoPathway scores...")

compute_pathway_scores <- function(expr_mat, gene_to_pathway_df, min_genes = 3) {
  pathways <- unique(gene_to_pathway_df$Pathway)

  scores <- sapply(pathways, function(pw) {
    genes <- gene_to_pathway_df %>%
      filter(Pathway == pw) %>%
      pull(Gene) %>%
      unique()
    genes_found <- intersect(genes, rownames(expr_mat))

    if (length(genes_found) < min_genes) {
      return(rep(NA, ncol(expr_mat)))
    }
    colMeans(expr_mat[genes_found, , drop = FALSE])
  })

  # samples × pathways matrix
  scores <- as.data.frame(scores)
  rownames(scores) <- colnames(expr_mat)
  scores
}

raw_pathway_scores <- compute_pathway_scores(expr_symbols, gene_to_pathway, min_genes = 3)

# Remove pathways with insufficient genes (NA columns)
n_before <- ncol(raw_pathway_scores)
raw_pathway_scores <- raw_pathway_scores[, colSums(is.na(raw_pathway_scores)) == 0]
n_after <- ncol(raw_pathway_scores)
message(sprintf("Pathway scores computed: %d pathways (dropped %d with < 3 genes in data)",
                n_after, n_before - n_after))

pathway_names <- colnames(raw_pathway_scores)

# =============================================================================
# PART 5: COMPUTE mitoPPS (PAIRWISE RATIO NORMALISATION)
# =============================================================================
# This is the actual mitoPPS algorithm from Monzel et al., verified against
# the source code in mitoPPS_RM.R, mitoPPS_MBM.R, and mitoPPS_calculations.R.
#
# NOT a simple score/sum normalisation. Instead:
#   1. For each sample, compute all pairwise ratios between pathways
#   2. Normalise each ratio by its global average across samples
#   3. Aggregate corrected ratios per pathway per sample
#
# This cancels out both total mitochondrial content AND intrinsic scale
# differences between pathways.

message("Computing mitoPPS (pairwise ratio normalisation)...")

compute_mitopps <- function(raw_scores) {
  # Convert to long format: one row per (sample, pathway)
  data_long <- raw_scores %>%
    as.data.frame() %>%
    rownames_to_column("Sample") %>%
    pivot_longer(cols = -Sample, names_to = "Pathway", values_to = "score")

  # Step 1: For each sample, compute all pairwise ratios
  message("  Step 1: Computing pairwise ratios...")
  data_ratios <- data_long %>%
    group_by(Sample) %>%
    nest() %>%
    mutate(pairs = map(data, ~ {
      df <- .x
      expand.grid(
        Pathway1 = df$Pathway,
        Pathway2 = df$Pathway,
        stringsAsFactors = FALSE
      ) %>%
        filter(Pathway1 != Pathway2) %>%
        left_join(df, by = c("Pathway1" = "Pathway")) %>%
        dplyr::rename(value1 = score) %>%
        left_join(df, by = c("Pathway2" = "Pathway")) %>%
        dplyr::rename(value2 = score) %>%
        mutate(ratio = value1 / value2)
    })) %>%
    dplyr::select(-data) %>%
    unnest(pairs)

  # Step 2: Normalise each ratio by its global average
  message("  Step 2: Normalising by global average ratios...")
  data_ratios <- data_ratios %>%
    group_by(Pathway1, Pathway2) %>%
    mutate(average_ratio = mean(ratio, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(corrected = ratio / average_ratio)

  # Step 3: Aggregate corrected ratios per sample per pathway
  message("  Step 3: Aggregating corrected ratios...")
  mitopps <- data_ratios %>%
    filter(!is.na(corrected)) %>%
    group_by(Sample, Pathway1) %>%
    summarize(mitoPPS = mean(corrected, na.rm = TRUE), .groups = "drop") %>%
    dplyr::rename(Pathway = Pathway1)

  # Convert back to wide format: samples × pathways
  mitopps_wide <- mitopps %>%
    pivot_wider(names_from = Pathway, values_from = mitoPPS) %>%
    column_to_rownames("Sample")

  # Ensure same column order as input
  mitopps_wide <- mitopps_wide[rownames(raw_scores), colnames(raw_scores)]

  mitopps_wide
}

mitopps_scores <- compute_mitopps(raw_pathway_scores)

# Verify: mitoPPS values should centre around 1.0
mitopps_mean <- mean(as.matrix(mitopps_scores), na.rm = TRUE)
mitopps_range <- range(as.matrix(mitopps_scores), na.rm = TRUE)
message(sprintf("mitoPPS: global mean = %.4f (should be ~1.0), range = [%.4f, %.4f]",
                mitopps_mean, mitopps_range[1], mitopps_range[2]))

# =============================================================================
# PART 6: ANNOTATE WITH METADATA
# =============================================================================

# Add group information to both score matrices
annotate_scores <- function(scores_df, coldata) {
  scores_df %>%
    as.data.frame() %>%
    rownames_to_column("sample") %>%
    left_join(coldata %>% rownames_to_column("sample") %>%
                dplyr::select(sample, group, timepoint, myc_status),
              by = "sample")
}

raw_pathway_scores_ann <- annotate_scores(raw_pathway_scores, coldata)
mitopps_scores_ann <- annotate_scores(mitopps_scores, coldata)

# =============================================================================
# PART 7: STATISTICAL ANALYSIS — RAW PATHWAY SCORES
# =============================================================================

message("Running statistical tests on raw pathway scores...")

run_anova <- function(scores_ann, pathway_names) {
  lapply(pathway_names, function(pw) {
    df <- scores_ann %>%
      dplyr::select(value = all_of(pw), timepoint, myc_status, group)

    fit <- tryCatch(
      aov(value ~ timepoint * myc_status, data = df),
      error = function(e) NULL
    )

    if (is.null(fit)) {
      return(data.frame(pathway = pw, effect = NA, F_value = NA, p_value = NA))
    }

    anova_table <- summary(fit)[[1]]

    data.frame(
      pathway = pw,
      effect = c("timepoint", "myc_status", "timepoint:myc_status"),
      F_value = anova_table[1:3, "F value"],
      p_value = anova_table[1:3, "Pr(>F)"]
    )
  }) %>%
    bind_rows() %>%
    filter(!is.na(p_value)) %>%
    group_by(effect) %>%
    mutate(padj = p.adjust(p_value, method = "BH")) %>%
    ungroup()
}

raw_stats <- run_anova(raw_pathway_scores_ann, pathway_names)

raw_stats_summary <- raw_stats %>%
  group_by(effect) %>%
  summarise(sig_005 = sum(padj < 0.05), sig_01 = sum(padj < 0.1), .groups = "drop")

message("\nRaw pathway scores — significant pathways (two-way ANOVA):")
print(raw_stats_summary)

# Group means for heatmap
raw_group_means <- raw_pathway_scores_ann %>%
  pivot_longer(cols = all_of(pathway_names), names_to = "pathway", values_to = "score") %>%
  group_by(pathway, group) %>%
  summarise(mean_score = mean(score), sd_score = sd(score), .groups = "drop")

# =============================================================================
# PART 8: STATISTICAL ANALYSIS — mitoPPS
# =============================================================================

message("Running statistical tests on mitoPPS scores...")

mitopps_stats <- run_anova(mitopps_scores_ann, pathway_names)

mitopps_stats_summary <- mitopps_stats %>%
  group_by(effect) %>%
  summarise(sig_005 = sum(padj < 0.05), sig_01 = sum(padj < 0.1), .groups = "drop")

message("\nmitoPPS — significant pathways (two-way ANOVA):")
print(mitopps_stats_summary)

# Group means for mitoPPS
mitopps_group_means <- mitopps_scores_ann %>%
  pivot_longer(cols = all_of(pathway_names), names_to = "pathway", values_to = "mitopps") %>%
  group_by(pathway, group) %>%
  summarise(mean_mitopps = mean(mitopps), sd_mitopps = sd(mitopps), .groups = "drop")

# =============================================================================
# PART 9: PAIRWISE COMPARISONS
# =============================================================================

message("Running pairwise comparisons...")

run_pairwise <- function(scores_ann, pathway_names, group_a, group_b, contrast_name) {
  results <- lapply(pathway_names, function(pw) {
    vals_a <- scores_ann %>% filter(group == group_a) %>% pull(all_of(pw))
    vals_b <- scores_ann %>% filter(group == group_b) %>% pull(all_of(pw))

    tt <- tryCatch(t.test(vals_b, vals_a), error = function(e) NULL)

    if (is.null(tt)) {
      return(data.frame(pathway = pw, contrast = contrast_name,
                        mean_a = NA, mean_b = NA, diff = NA, p_value = NA))
    }

    data.frame(
      pathway = pw,
      contrast = contrast_name,
      mean_a = mean(vals_a),
      mean_b = mean(vals_b),
      diff = mean(vals_b) - mean(vals_a),
      p_value = tt$p.value
    )
  }) %>% bind_rows()

  results %>% mutate(padj = p.adjust(p_value, method = "BH"))
}

# Raw pathway scores — pairwise
raw_pw_myc6W    <- run_pairwise(raw_pathway_scores_ann, pathway_names, "6W_neg", "6W_pos", "Myc_effect_6W")
raw_pw_myc12W   <- run_pairwise(raw_pathway_scores_ann, pathway_names, "12W_neg", "12W_pos", "Myc_effect_12W")
raw_pw_time_pos <- run_pairwise(raw_pathway_scores_ann, pathway_names, "6W_pos", "12W_pos", "Temporal_Myc+")
raw_pw_time_neg <- run_pairwise(raw_pathway_scores_ann, pathway_names, "6W_neg", "12W_neg", "Temporal_Myc-")
raw_pairwise <- bind_rows(raw_pw_myc6W, raw_pw_myc12W, raw_pw_time_pos, raw_pw_time_neg)

# mitoPPS — pairwise
pps_pw_myc6W    <- run_pairwise(mitopps_scores_ann, pathway_names, "6W_neg", "6W_pos", "Myc_effect_6W")
pps_pw_myc12W   <- run_pairwise(mitopps_scores_ann, pathway_names, "12W_neg", "12W_pos", "Myc_effect_12W")
pps_pw_time_pos <- run_pairwise(mitopps_scores_ann, pathway_names, "6W_pos", "12W_pos", "Temporal_Myc+")
pps_pw_time_neg <- run_pairwise(mitopps_scores_ann, pathway_names, "6W_neg", "12W_neg", "Temporal_Myc-")
pps_pairwise <- bind_rows(pps_pw_myc6W, pps_pw_myc12W, pps_pw_time_pos, pps_pw_time_neg)

pairwise_summary <- pps_pairwise %>%
  group_by(contrast) %>%
  summarise(sig_005 = sum(padj < 0.05, na.rm = TRUE),
            sig_01 = sum(padj < 0.1, na.rm = TRUE),
            .groups = "drop")

message("\nmitoPPS pairwise comparisons — significant pathways:")
print(pairwise_summary)

# =============================================================================
# PART 10: VISUALISATIONS
# =============================================================================

message("Generating visualisations...")

group_colours <- c("6W_neg" = "#377EB8", "6W_pos" = "#E41A1C",
                   "12W_neg" = "#4DAF4A", "12W_pos" = "#FF7F00")
group_shapes  <- c("6W_neg" = 16, "6W_pos" = 17,
                   "12W_neg" = 15, "12W_pos" = 18)
tier1_colours <- c(
  "Metabolism" = "#2166AC",
  "OXPHOS" = "#B2182B",
  "Protein import, sorting, and homeostasis" = "#1B9E77",
  "Mitochondrial central dogma" = "#D95F02",
  "Mitochondrial dynamics and surveillance" = "#7570B3",
  "Small molecule transport" = "#E7298A",
  "Signaling" = "#66A61E"
)

# --- 10a. PCA of raw MitoPathway scores ---
pca_raw <- prcomp(raw_pathway_scores_ann[, pathway_names], center = TRUE, scale. = TRUE)
pca_raw_df <- as.data.frame(pca_raw$x[, 1:3]) %>%
  mutate(sample = raw_pathway_scores_ann$sample,
         group = raw_pathway_scores_ann$group,
         timepoint = raw_pathway_scores_ann$timepoint,
         myc_status = raw_pathway_scores_ann$myc_status)

var_explained_raw <- round(100 * summary(pca_raw)$importance[2, 1:3], 1)

p_pca_raw <- ggplot(pca_raw_df, aes(x = PC1, y = PC2, colour = group, shape = group)) +
  geom_point(size = 4, alpha = 0.8) +
  scale_colour_manual(values = group_colours) +
  scale_shape_manual(values = group_shapes) +
  labs(title = "PCA of raw MitoPathway scores",
       subtitle = "Based on mean normalised counts of MitoCarta3.0 pathway genes",
       x = sprintf("PC1 (%.1f%%)", var_explained_raw[1]),
       y = sprintf("PC2 (%.1f%%)", var_explained_raw[2])) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(mitopps_fig_dir, "pca_raw_pathway_scores.pdf"),
       p_pca_raw, width = 8, height = 6)

# --- 10b. PCA of mitoPPS ---
# Following Monzel et al. Fig3D: PCA with scale. = TRUE
pca_pps <- prcomp(mitopps_scores_ann[, pathway_names], center = TRUE, scale. = TRUE)
pca_pps_df <- as.data.frame(pca_pps$x[, 1:3]) %>%
  mutate(sample = mitopps_scores_ann$sample,
         group = mitopps_scores_ann$group,
         timepoint = mitopps_scores_ann$timepoint,
         myc_status = mitopps_scores_ann$myc_status)

var_explained_pps <- round(100 * summary(pca_pps)$importance[2, 1:3], 1)

p_pca_pps <- ggplot(pca_pps_df, aes(x = PC1, y = PC2, colour = group, shape = group)) +
  geom_point(size = 4, alpha = 0.8) +
  scale_colour_manual(values = group_colours) +
  scale_shape_manual(values = group_shapes) +
  labs(title = "PCA of mitoPPS (pathway prioritisation scores)",
       subtitle = "Pairwise-ratio normalised: independent of total mito content",
       x = sprintf("PC1 (%.1f%%)", var_explained_pps[1]),
       y = sprintf("PC2 (%.1f%%)", var_explained_pps[2])) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(mitopps_fig_dir, "pca_mitopps.pdf"),
       p_pca_pps, width = 8, height = 6)

# --- 10c. Heatmap: mitoPPS (group means, log10-transformed) ---
# Following Monzel et al. Fig3C: log10 of mitoPPS values for heatmap
pps_means_wide <- mitopps_group_means %>%
  dplyr::select(pathway, group, mean_mitopps) %>%
  pivot_wider(names_from = group, values_from = mean_mitopps) %>%
  column_to_rownames("pathway") %>%
  as.matrix()

# Log10 transform as in Fig3C_Heatmap.R line 74
pps_means_log10 <- log10(pps_means_wide)

# Z-score across groups for visualisation
pps_means_z <- t(scale(t(pps_means_log10)))

pps_row_annotation <- data.frame(
  Category = pathway_tier1_map[rownames(pps_means_z)],
  row.names = rownames(pps_means_z)
)
valid_rows_pps <- !is.na(pps_row_annotation$Category)

pdf(file.path(mitopps_fig_dir, "heatmap_mitopps.pdf"), width = 8, height = 18)
pheatmap(pps_means_z[valid_rows_pps, ],
         color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
         cluster_cols = FALSE,
         annotation_row = pps_row_annotation[valid_rows_pps, , drop = FALSE],
         annotation_colors = list(Category = tier1_colours),
         fontsize_row = 6,
         fontsize_col = 10,
         main = "mitoPPS (log10, z-scored group means)")
dev.off()

# --- 10d. Heatmap: raw pathway scores (group means, z-scored) ---
raw_means_wide <- raw_group_means %>%
  dplyr::select(pathway, group, mean_score) %>%
  pivot_wider(names_from = group, values_from = mean_score) %>%
  column_to_rownames("pathway") %>%
  as.matrix()

raw_means_z <- t(scale(t(raw_means_wide)))

raw_row_annotation <- data.frame(
  Category = pathway_tier1_map[rownames(raw_means_z)],
  row.names = rownames(raw_means_z)
)
valid_rows_raw <- !is.na(raw_row_annotation$Category)

pdf(file.path(mitopps_fig_dir, "heatmap_raw_pathway_scores.pdf"), width = 8, height = 18)
pheatmap(raw_means_z[valid_rows_raw, ],
         color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
         cluster_cols = FALSE,
         annotation_row = raw_row_annotation[valid_rows_raw, , drop = FALSE],
         annotation_colors = list(Category = tier1_colours),
         fontsize_row = 6,
         fontsize_col = 10,
         main = "Raw MitoPathway scores (z-scored group means)")
dev.off()

# --- 10e. Dot plot: significantly reprioritised pathways (mitoPPS, Myc effect) ---
myc_effect_pps <- bind_rows(pps_pw_myc6W, pps_pw_myc12W) %>%
  filter(padj < 0.1) %>%
  mutate(
    pathway_label = str_replace_all(pathway, "_", " "),
    tier1 = pathway_tier1_map[pathway],
    direction = ifelse(diff > 0, "Prioritised in Myc+", "Deprioritised in Myc+")
  )

if (nrow(myc_effect_pps) > 0) {
  p_myc_pps <- ggplot(myc_effect_pps,
                       aes(x = diff, y = reorder(pathway_label, diff))) +
    geom_point(aes(size = -log10(padj), colour = tier1), alpha = 0.8) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    facet_wrap(~ contrast, ncol = 2) +
    scale_colour_manual(values = tier1_colours, na.value = "grey50") +
    labs(title = "mitoPPS: Myc-driven pathway reprioritisation",
         subtitle = "Pathways significantly reprioritised (padj < 0.1) in Myc+ vs Myc-",
         x = "Difference in mitoPPS (Myc+ - Myc-)",
         y = NULL,
         colour = "Category",
         size = "-log10(padj)") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 7),
          plot.title = element_text(face = "bold"))

  ggsave(file.path(mitopps_fig_dir, "dotplot_mitopps_myc_effect.pdf"),
         p_myc_pps, width = 14, height = max(6, nrow(myc_effect_pps) * 0.2))
}

# --- 10f. Box plots for top pathways ---
top_myc_pathways <- mitopps_stats %>%
  filter(effect == "myc_status") %>%
  arrange(p_value) %>%
  head(12) %>%
  pull(pathway)

if (length(top_myc_pathways) > 0) {
  boxplot_data <- mitopps_scores_ann %>%
    pivot_longer(cols = all_of(top_myc_pathways),
                 names_to = "pathway", values_to = "mitopps") %>%
    mutate(pathway_label = str_replace_all(pathway, "_", " "))

  p_box <- ggplot(boxplot_data, aes(x = group, y = mitopps, fill = group)) +
    geom_boxplot(alpha = 0.7, outlier.size = 1) +
    geom_point(size = 1.2, alpha = 0.6, position = position_jitter(width = 0.15)) +
    facet_wrap(~ pathway_label, scales = "free_y", ncol = 4) +
    scale_fill_manual(values = group_colours) +
    labs(title = "Top Myc-affected mitoPPS pathways",
         subtitle = "Top 12 pathways by Myc effect (two-way ANOVA p-value)",
         x = NULL, y = "mitoPPS") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
          strip.text = element_text(size = 7),
          plot.title = element_text(face = "bold"),
          legend.position = "none")

  ggsave(file.path(mitopps_fig_dir, "boxplots_top_mitopps_myc.pdf"),
         p_box, width = 12, height = 10)
}

# --- 10g. Total mito expression by group ---
total_mito <- data.frame(
  sample = raw_pathway_scores_ann$sample,
  group = raw_pathway_scores_ann$group,
  total_mito_score = rowSums(raw_pathway_scores_ann[, pathway_names])
)

p_total_mito <- ggplot(total_mito, aes(x = group, y = total_mito_score, fill = group)) +
  geom_boxplot(alpha = 0.7) +
  geom_point(size = 2, alpha = 0.6, position = position_jitter(width = 0.15)) +
  scale_fill_manual(values = group_colours) +
  labs(title = "Total mitochondrial pathway expression",
       subtitle = "Sum of all raw MitoPathway scores per sample (reflects mito content)",
       x = NULL, y = "Total MitoPathway score") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"), legend.position = "none")

ggsave(file.path(mitopps_fig_dir, "boxplot_total_mito_expression.pdf"),
       p_total_mito, width = 6, height = 5)

# =============================================================================
# PART 11: SAVE RESULTS
# =============================================================================

message("Saving results...")

mitopps_results <- list(
  # Per-sample scores (samples × pathways + metadata)
  raw_pathway_scores = raw_pathway_scores_ann,
  mitopps_scores = mitopps_scores_ann,

  # Group-level summaries
  raw_group_means = raw_group_means,
  mitopps_group_means = mitopps_group_means,

  # Statistical results
  raw_anova = raw_stats,
  mitopps_anova = mitopps_stats,
  raw_pairwise = raw_pairwise,
  mitopps_pairwise = pps_pairwise,

  # Pathway annotations
  gene_to_pathway = gene_to_pathway,
  pathway_levels = pathway_levels,
  pathway_tier1_map = pathway_tier1_map,

  # PCA objects
  pca_raw = pca_raw,
  pca_mitopps = pca_pps,

  # Metadata
  n_pathways = n_after,
  n_mitocarta_genes_found = length(overlap),
  analysis_date = Sys.Date(),
  description = paste(
    "MitoPPS analysis using pairwise ratio normalisation.",
    "Raw pathway scores from DESeq2 normalised counts (linear scale).",
    "mitoPPS computed per Monzel et al. (2025) algorithm."
  )
)

saveRDS(mitopps_results, here("results", "mitopps_scores.rds"))

# Save summary tables as CSV
raw_stats %>%
  filter(padj < 0.1) %>%
  arrange(effect, padj) %>%
  write.csv(file.path(mitopps_fig_dir, "raw_pathway_anova_sig.csv"), row.names = FALSE)

mitopps_stats %>%
  filter(padj < 0.1) %>%
  arrange(effect, padj) %>%
  write.csv(file.path(mitopps_fig_dir, "mitopps_anova_sig.csv"), row.names = FALSE)

pps_pairwise %>%
  filter(padj < 0.1) %>%
  arrange(contrast, padj) %>%
  write.csv(file.path(mitopps_fig_dir, "mitopps_pairwise_sig.csv"), row.names = FALSE)

message("\n", strrep("=", 70))
message("mitoPPS ANALYSIS COMPLETE")
message(strrep("=", 70))
message(sprintf("Pathways analysed: %d", n_after))
message(sprintf("Results saved to: results/mitopps_scores.rds"))
message(sprintf("Figures saved to: %s", mitopps_fig_dir))
message(strrep("=", 70))
