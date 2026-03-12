# scripts/12_cell_death_pathway_analysis.R
# =============================================================================
# Cell Death Pathway Analysis: How Do 15 Regulated Cell Death Modalities
# Change Between 6W and 12W in Myc+ vs Myc- Samples?
# =============================================================================
#
# Motivation:
#   We want to understand whether cell death pathway activity changes over
#   time (6W → 12W) and whether these changes are Myc-specific or
#   developmental. Using gene sets from Tang et al. (2024) covering 15
#   regulated cell death modalities, we run fGSEA on five DESeq2 contrasts
#   and extract leading edge genes to determine whether samples become more
#   resistant or more susceptible to each death modality over time.
#
# Reference:
#   Tang et al. (2024) Comput Struct Biotechnol J
#   https://doi.org/10.1016/j.csbj.2024.08.012
#
# Approach:
#   1. Load 15 cell death gene sets (human) and convert to mouse orthologs
#   2. Create ranked gene lists (Wald statistic) for 5 contrasts:
#      - Myc+ vs Myc- at 6W (cross-sectional early Myc effect)
#      - Myc+ vs Myc- at 12W (cross-sectional late Myc effect)
#      - 12W vs 6W in Myc- (developmental baseline)
#      - 12W vs 6W in Myc+ (Myc+ temporal progression)
#      - Interaction (timepoint × myc_status)
#   3. Run fGSEA on all five contrasts
#   4. Extract leading edge genes, annotate with pro-/anti-death roles
#   5. Interpret whether 12W samples are more resistant or susceptible
#
# Input:
#   - data/cell-death/*.csv           (15 cell death gene set CSVs)
#   - results/interaction_results.rds (DESeq2 results, all contrasts)
#   - results/ortholog_table.rds      (human-mouse ortholog mapping)
#
# Output:
#   - results/cell_death_fgsea.rds
#   - outputs/cell_death/
#
# =============================================================================

source(here::here("scripts", "00_setup_packages.R"))

# === Output directory ===
cd_dir <- here("outputs", "cell_death")
dir.create(cd_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# PART 1: LOAD AND CONVERT CELL DEATH GENE SETS
# =============================================================================

message("Loading cell death gene sets from CSV files...")

# --- 1a. Read all CSV files ---
cd_path <- here("data", "cell-death")
cd_files <- list.files(cd_path, pattern = "\\.csv$", full.names = TRUE)

message(sprintf("  Found %d cell death CSV files", length(cd_files)))

# Read all files into a single data frame
cd_all <- lapply(cd_files, function(f) {
  # Some files are cp1252, others UTF-8 — read raw and convert
  df <- read.csv(f, stringsAsFactors = FALSE, fileEncoding = "latin1")
  # Ensure all character columns are valid UTF-8
  df[] <- lapply(df, function(col) {
    if (is.character(col)) iconv(col, from = "latin1", to = "UTF-8", sub = "")
    else col
  })
  df
}) %>% bind_rows()

message(sprintf("  Total rows across all files: %d", nrow(cd_all)))
message(sprintf("  Unique death types: %s", paste(unique(cd_all$deathtype), collapse = ", ")))

# --- 1b. Create gene sets: death type → unique human gene symbols ---
cd_human_sets <- cd_all %>%
  group_by(deathtype) %>%
  summarise(genes = list(unique(gene)), .groups = "drop") %>%
  deframe()

message("\nHuman gene set sizes:")
for (nm in sort(names(cd_human_sets))) {
  message(sprintf("  %-35s: %d genes", nm, length(cd_human_sets[[nm]])))
}

# --- 1c. Convert human → mouse orthologs ---
ortholog_table <- readRDS(here("results", "ortholog_table.rds"))

human_to_mouse_map <- deframe(
  ortholog_table[, c("hsapiens_homolog_associated_gene_name", "external_gene_name")]
)

cd_mouse_sets <- lapply(cd_human_sets, function(genes) {
  mapped <- human_to_mouse_map[genes]
  unique(na.omit(mapped))
})

# Prefix for clarity
names(cd_mouse_sets) <- paste0("CD_", names(cd_mouse_sets))

# Report coverage
message("\nMouse ortholog mapping coverage:")
for (i in seq_along(cd_human_sets)) {
  nm_h <- names(cd_human_sets)[i]
  nm_m <- paste0("CD_", nm_h)
  n_h <- length(cd_human_sets[[nm_h]])
  n_m <- length(cd_mouse_sets[[nm_m]])
  message(sprintf("  %-35s: %d / %d mapped (%.0f%%)", nm_h, n_m, n_h,
                  100 * n_m / n_h))
}

# Also keep a lookup: mouse symbol → human descriptions (for leading edge)
# For genes with multiple entries per death type, collapse descriptions
cd_descriptions <- cd_all %>%
  mutate(
    mouse_symbol = human_to_mouse_map[gene]
  ) %>%
  filter(!is.na(mouse_symbol)) %>%
  dplyr::select(deathtype, human_gene = gene, mouse_symbol, description, comment)

# =============================================================================
# PART 2: CREATE RANKED GENE LISTS FOR 5 CONTRASTS
# =============================================================================

message("\nCreating ranked gene lists from DESeq2 Wald statistics...")

interaction_results <- readRDS(here("results", "interaction_results.rds"))

# Ensembl → symbol mapping
ensembl_to_symbol <- setNames(
  ortholog_table$external_gene_name,
  ortholog_table$ensembl_gene_id
)

# Rank creation function (from script 10)
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

# Five contrasts
ranks_list <- list(
  myc_effect_6W   = create_ranks(interaction_results$myc_6W_raw, ensembl_to_symbol),
  myc_effect_12W  = create_ranks(interaction_results$myc_12W_raw, ensembl_to_symbol),
  temporal_neg    = create_ranks(interaction_results$timepoint_neg_raw, ensembl_to_symbol),
  temporal_pos    = create_ranks(interaction_results$timepoint_pos_raw, ensembl_to_symbol),
  interaction     = create_ranks(interaction_results$interaction_raw, ensembl_to_symbol)
)

for (nm in names(ranks_list)) {
  message(sprintf("  %-18s: %d genes ranked, range [%.2f, %.2f]",
                  nm, length(ranks_list[[nm]]),
                  min(ranks_list[[nm]]), max(ranks_list[[nm]])))
}

# =============================================================================
# PART 3: RUN fGSEA ON ALL 5 CONTRASTS
# =============================================================================

message("\nRunning fGSEA on cell death gene sets across 5 contrasts...")

set.seed(42)

fgsea_results <- lapply(names(ranks_list), function(contrast_name) {
  message(sprintf("  Running fGSEA: %s", contrast_name))
  res <- fgsea(
    pathways    = cd_mouse_sets,
    stats       = ranks_list[[contrast_name]],
    minSize     = 5,
    maxSize     = 2000,  # some sets are large (Autophagy ~1195 human genes)
    nPermSimple = 10000
  )
  res$contrast <- contrast_name
  res
})

names(fgsea_results) <- names(ranks_list)

# Summary
message("\nfGSEA results summary (padj < 0.05 / padj < 0.10):")
for (nm in names(fgsea_results)) {
  n05 <- sum(fgsea_results[[nm]]$padj < 0.05, na.rm = TRUE)
  n10 <- sum(fgsea_results[[nm]]$padj < 0.10, na.rm = TRUE)
  message(sprintf("  %-18s: %d (padj < 0.05), %d (padj < 0.10) out of %d",
                  nm, n05, n10, nrow(fgsea_results[[nm]])))
}

# =============================================================================
# PART 4: COMBINED RESULTS TABLE AND NES HEATMAP
# =============================================================================

message("\nBuilding combined results table...")

# --- 4a. Wide NES table for heatmap ---
fgsea_combined <- bind_rows(fgsea_results) %>%
  mutate(pathway_label = str_replace_all(pathway, "^CD_", "") %>%
           str_replace_all("_", " "))

nes_wide <- fgsea_combined %>%
  dplyr::select(pathway_label, contrast, NES) %>%
  pivot_wider(names_from = contrast, values_from = NES) %>%
  column_to_rownames("pathway_label")

padj_wide <- fgsea_combined %>%
  dplyr::select(pathway_label, contrast, padj) %>%
  pivot_wider(names_from = contrast, values_from = padj) %>%
  column_to_rownames("pathway_label")

# Significance stars for annotation
sig_mat <- matrix("", nrow = nrow(padj_wide), ncol = ncol(padj_wide),
                  dimnames = dimnames(padj_wide))
sig_mat[padj_wide < 0.10] <- "."
sig_mat[padj_wide < 0.05] <- "*"
sig_mat[padj_wide < 0.01] <- "**"
sig_mat[padj_wide < 0.001] <- "***"

# --- 4b. NES Heatmap ---
message("  Generating NES heatmap...")

# Column labels
contrast_labels <- c(
  myc_effect_6W  = "Myc effect\n6W",
  myc_effect_12W = "Myc effect\n12W",
  temporal_neg   = "Temporal\nMyc-",
  temporal_pos   = "Temporal\nMyc+",
  interaction    = "Interaction\n(Time x Myc)"
)

# Reorder columns
col_order <- c("myc_effect_6W", "myc_effect_12W", "temporal_neg",
               "temporal_pos", "interaction")
nes_mat <- as.matrix(nes_wide[, col_order])
sig_mat_ordered <- sig_mat[, col_order]

# Column annotation: contrast type
col_anno <- HeatmapAnnotation(
  Type = c("Cross-sectional", "Cross-sectional", "Temporal", "Temporal", "Interaction"),
  col = list(Type = c("Cross-sectional" = "#4DAF4A",
                       "Temporal" = "#377EB8",
                       "Interaction" = "#E41A1C")),
  annotation_name_side = "left"
)

# Colour scale
max_abs <- max(abs(nes_mat), na.rm = TRUE)
col_fun <- circlize::colorRamp2(
  c(-max_abs, 0, max_abs),
  c("#2166AC", "white", "#B2182B")
)

pdf(file.path(cd_dir, "heatmap_nes_cell_death_pathways.pdf"),
    width = 10, height = 8)
ht <- ComplexHeatmap::Heatmap(
  nes_mat,
  name = "NES",
  col  = col_fun,
  cluster_columns = FALSE,
  cluster_rows    = TRUE,
  column_labels   = contrast_labels[col_order],
  row_names_side  = "left",
  row_names_gp    = gpar(fontsize = 9),
  column_names_gp = gpar(fontsize = 10),
  top_annotation  = col_anno,
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sig_mat_ordered[i, j], x, y, gp = gpar(fontsize = 8))
  },
  column_title = "Cell Death Pathway Enrichment (fGSEA NES) Across Contrasts",
  column_title_gp = gpar(fontsize = 12, fontface = "bold"),
  heatmap_legend_param = list(
    title = "NES",
    title_gp = gpar(fontsize = 10),
    labels_gp = gpar(fontsize = 9)
  )
)
draw(ht)
dev.off()

# =============================================================================
# PART 5: DOTPLOTS FOR TEMPORAL AND INTERACTION CONTRASTS
# =============================================================================

message("Generating dotplots...")

# --- 5a. Temporal comparison dotplot (Myc+ vs Myc-) ---
temporal_combined <- fgsea_combined %>%
  filter(contrast %in% c("temporal_neg", "temporal_pos")) %>%
  mutate(
    contrast_label = ifelse(contrast == "temporal_neg",
                            "Myc- (12W vs 6W)",
                            "Myc+ (12W vs 6W)"),
    sig = padj < 0.05
  )

p_temporal <- ggplot(temporal_combined,
                     aes(x = NES, y = reorder(pathway_label, NES),
                         colour = contrast_label)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(aes(size = -log10(padj),
                 shape = sig),
             alpha = 0.8) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1),
                     name = "padj < 0.05") +
  scale_colour_manual(values = c("Myc- (12W vs 6W)" = "#377EB8",
                                  "Myc+ (12W vs 6W)" = "#E41A1C"),
                      name = "Contrast") +
  scale_size_continuous(range = c(1.5, 5), name = expression(-log[10](padj))) +
  labs(
    title = "Cell Death Pathways: Temporal Change (12W vs 6W)",
    subtitle = paste0(
      "Positive NES = pathway activity increases at 12W\n",
      "Negative NES = pathway activity decreases at 12W\n",
      "Filled = padj < 0.05  |  Open = not significant"
    ),
    x = "fGSEA NES",
    y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.y  = element_text(size = 9),
    plot.title   = element_text(face = "bold"),
    legend.position = "right"
  )

ggsave(file.path(cd_dir, "dotplot_temporal_cell_death.pdf"),
       p_temporal, width = 11, height = 7)

# --- 5b. Cross-sectional comparison dotplot (Myc effect at 6W vs 12W) ---
xs_combined <- fgsea_combined %>%
  filter(contrast %in% c("myc_effect_6W", "myc_effect_12W")) %>%
  mutate(
    contrast_label = ifelse(contrast == "myc_effect_6W",
                            "Myc+ vs Myc- (6W)",
                            "Myc+ vs Myc- (12W)"),
    sig = padj < 0.05
  )

p_xs <- ggplot(xs_combined,
               aes(x = NES, y = reorder(pathway_label, NES),
                   colour = contrast_label)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(aes(size = -log10(padj),
                 shape = sig),
             alpha = 0.8) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1),
                     name = "padj < 0.05") +
  scale_colour_manual(values = c("Myc+ vs Myc- (6W)" = "#377EB8",
                                  "Myc+ vs Myc- (12W)" = "#E41A1C"),
                      name = "Contrast") +
  scale_size_continuous(range = c(1.5, 5), name = expression(-log[10](padj))) +
  labs(
    title = "Cell Death Pathways: Myc Effect at Each Timepoint",
    subtitle = paste0(
      "Positive NES = pathway upregulated in Myc+\n",
      "Negative NES = pathway downregulated in Myc+"
    ),
    x = "fGSEA NES",
    y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.y  = element_text(size = 9),
    plot.title   = element_text(face = "bold"),
    legend.position = "right"
  )

ggsave(file.path(cd_dir, "dotplot_crosssectional_cell_death.pdf"),
       p_xs, width = 11, height = 7)

# --- 5c. Interaction dotplot ---
int_data <- fgsea_combined %>%
  filter(contrast == "interaction") %>%
  mutate(sig = padj < 0.05)

p_int <- ggplot(int_data,
                aes(x = NES, y = reorder(pathway_label, NES))) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(aes(size = -log10(padj),
                 colour = NES,
                 shape = sig),
             alpha = 0.85) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1),
                     name = "padj < 0.05") +
  scale_colour_gradient2(low = "#2166AC", mid = "grey80", high = "#B2182B",
                         midpoint = 0, name = "NES") +
  scale_size_continuous(range = c(2, 6), name = expression(-log[10](padj))) +
  labs(
    title = "Cell Death Pathways: Interaction (Timepoint x Myc Status)",
    subtitle = paste0(
      "Positive NES = Myc effect on pathway STRENGTHENS at 12W\n",
      "Negative NES = Myc effect on pathway WEAKENS at 12W"
    ),
    x = "fGSEA NES (interaction)",
    y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.y  = element_text(size = 9),
    plot.title   = element_text(face = "bold"),
    legend.position = "right"
  )

ggsave(file.path(cd_dir, "dotplot_interaction_cell_death.pdf"),
       p_int, width = 10, height = 7)

# =============================================================================
# PART 6: PAIRED NES COMPARISON (TEMPORAL: Myc+ vs Myc-)
# =============================================================================
# For each death pathway, show NES for temporal_pos and temporal_neg
# side by side — highlights Myc-specific vs shared temporal changes

message("Generating paired NES comparison...")

paired_nes_data <- fgsea_combined %>%
  filter(contrast %in% c("temporal_neg", "temporal_pos")) %>%
  dplyr::select(pathway_label, contrast, NES, padj)

paired_nes_long <- paired_nes_data %>%
  mutate(
    genotype = ifelse(contrast == "temporal_neg", "Myc-", "Myc+"),
    sig = padj < 0.05
  )

# Order pathways by mean NES
pw_order <- paired_nes_data %>%
  group_by(pathway_label) %>%
  summarise(mean_nes = mean(NES), .groups = "drop") %>%
  arrange(mean_nes) %>%
  pull(pathway_label)

paired_nes_long$pathway_label <- factor(paired_nes_long$pathway_label,
                                         levels = pw_order)

p_paired <- ggplot(paired_nes_long,
                   aes(x = NES, y = pathway_label, colour = genotype)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_line(aes(group = pathway_label), colour = "grey60", linewidth = 0.4) +
  geom_point(aes(shape = sig), size = 3.5, alpha = 0.85) +
  scale_colour_manual(values = c("Myc-" = "#377EB8", "Myc+" = "#E41A1C"),
                      name = "Genotype") +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1),
                     name = "padj < 0.05") +
  labs(
    title = "Temporal Change in Cell Death Pathways: Myc+ vs Myc-",
    subtitle = paste0(
      "NES from fGSEA (12W vs 6W)  |  Connecting lines show divergence\n",
      "Filled = padj < 0.05  |  Open = not significant"
    ),
    x = "fGSEA NES (12W vs 6W)",
    y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.y  = element_text(size = 9),
    plot.title   = element_text(face = "bold"),
    legend.position = "right"
  )

ggsave(file.path(cd_dir, "dotplot_paired_temporal_cell_death.pdf"),
       p_paired, width = 11, height = 7)

# =============================================================================
# PART 7: LEADING EDGE ANALYSIS
# =============================================================================
# Extract leading edge genes from significant pathways in the temporal
# and interaction contrasts, and annotate with their functional descriptions

message("\nExtracting and annotating leading edge genes...")

# Focus on temporal_pos and interaction (Myc-relevant changes)
contrasts_of_interest <- c("temporal_pos", "temporal_neg", "interaction")

le_results <- list()

for (ctr in contrasts_of_interest) {
  res <- fgsea_results[[ctr]]

  # All pathways (not just significant — for comprehensive annotation)
  for (i in seq_len(nrow(res))) {
    pw <- res$pathway[i]
    pw_clean <- str_replace(pw, "^CD_", "")
    le_genes <- unlist(res$leadingEdge[i])

    if (length(le_genes) == 0) next

    # Look up descriptions from the original CSV data
    le_annot <- cd_descriptions %>%
      filter(deathtype == pw_clean, mouse_symbol %in% le_genes) %>%
      dplyr::select(mouse_symbol, human_gene, description, comment) %>%
      distinct(mouse_symbol, .keep_all = TRUE)

    # Genes without annotation (mapped but no description match)
    missing <- setdiff(le_genes, le_annot$mouse_symbol)
    if (length(missing) > 0) {
      le_annot <- bind_rows(
        le_annot,
        tibble(mouse_symbol = missing, human_gene = NA,
               description = NA, comment = NA)
      )
    }

    le_annot$contrast  <- ctr
    le_annot$pathway   <- pw
    le_annot$NES       <- res$NES[i]
    le_annot$padj      <- res$padj[i]
    le_annot$le_size   <- length(le_genes)
    le_annot$set_size  <- res$size[i]

    le_results[[paste(ctr, pw, sep = "::")]] <- le_annot
  }
}

le_combined <- bind_rows(le_results)

message(sprintf("  Leading edge entries: %d across %d pathway-contrast combinations",
                nrow(le_combined), length(le_results)))

# =============================================================================
# PART 8: CLASSIFY LEADING EDGE GENES AS PRO- OR ANTI-DEATH
# =============================================================================
# Use the "comment" field from the CSVs to classify direction.
# Genes whose silencing/inhibition promotes death → anti-death (protective)
# Genes whose activation/overexpression promotes death → pro-death
#
# This is done through keyword heuristics on the comment field.

message("Classifying leading edge genes by pro-/anti-death role...")

classify_death_role <- function(comment_text) {
  if (is.na(comment_text) || comment_text == "") return("unclassified")
  # Sanitise encoding: some CSVs are cp1252, tolower() needs valid UTF-8
  txt <- iconv(comment_text, from = "", to = "UTF-8", sub = "byte")
  txt <- tolower(txt)

  # Pro-death indicators
  pro_keywords <- c("promot.*death", "promot.*apoptosis", "promot.*pyroptosis",
                     "promot.*ferroptosis", "promot.*necroptosis", "promot.*autophagy",
                     "induc.*death", "induc.*apoptosis", "induc.*pyroptosis",
                     "induc.*ferroptosis", "trigger", "activat.*caspase",
                     "enhanc.*death", "enhanc.*apoptosis", "execut",
                     "essential for.*death", "required.*death", "required.*ferroptosis",
                     "mediat.*death", "pro-apoptotic", "pro-death",
                     "positiv.*regulat.*death", "positiv.*regulat.*apoptosis",
                     "upregulat.*promot", "overexpression.*induc.*apoptosis",
                     "overexpression.*induc.*death")

  # Anti-death indicators
  anti_keywords <- c("inhibit.*death", "inhibit.*apoptosis", "inhibit.*ferroptosis",
                      "inhibit.*pyroptosis", "inhibit.*necroptosis",
                      "suppress.*death", "suppress.*apoptosis",
                      "protect.*against", "protect.*from", "resist",
                      "anti-apoptotic", "anti-death", "prevent.*death",
                      "negativ.*regulat.*death", "negativ.*regulat.*apoptosis",
                      "block.*death", "confer.*against", "survival",
                      "silencing.*confer", "knockdown.*protect")

  pro_score  <- sum(sapply(pro_keywords, function(k) grepl(k, txt)))
  anti_score <- sum(sapply(anti_keywords, function(k) grepl(k, txt)))

  if (pro_score > anti_score) return("pro-death")
  if (anti_score > pro_score) return("anti-death")
  return("ambiguous")
}

le_combined <- le_combined %>%
  mutate(death_role = sapply(comment, classify_death_role))

# Summary
role_summary <- le_combined %>%
  filter(contrast %in% c("temporal_pos", "interaction")) %>%
  group_by(contrast, pathway, death_role) %>%
  summarise(n = n(), .groups = "drop")

message("\nLeading edge death role classification (temporal_pos & interaction):")
role_summary %>%
  pivot_wider(names_from = death_role, values_from = n, values_fill = 0) %>%
  print()

# =============================================================================
# PART 9: BIOLOGICAL INTERPRETATION — RESISTANCE VS SUSCEPTIBILITY
# =============================================================================
# For each death modality in the temporal contrasts:
# - Positive NES (upregulated at 12W):
#   - If dominated by pro-death LE genes → MORE SUSCEPTIBLE at 12W
#   - If dominated by anti-death LE genes → MORE RESISTANT at 12W
# - Negative NES (downregulated at 12W):
#   - If dominated by pro-death LE genes → MORE RESISTANT at 12W
#   - If dominated by anti-death LE genes → MORE SUSCEPTIBLE at 12W

message("\n", strrep("=", 70))
message("BIOLOGICAL INTERPRETATION: CELL DEATH RESISTANCE AT 12W")
message(strrep("=", 70))

interpretation <- list()

for (ctr in c("temporal_pos", "temporal_neg")) {
  ctr_label <- ifelse(ctr == "temporal_pos", "Myc+", "Myc-")
  message(sprintf("\n--- %s (12W vs 6W) ---", ctr_label))

  res <- fgsea_results[[ctr]] %>% arrange(pval)

  for (i in seq_len(nrow(res))) {
    pw <- res$pathway[i]
    pw_clean <- str_replace(pw, "^CD_", "")
    nes <- res$NES[i]
    p <- res$padj[i]

    # Get leading edge role breakdown
    le_this <- le_combined %>%
      filter(contrast == ctr, pathway == pw)

    n_pro  <- sum(le_this$death_role == "pro-death")
    n_anti <- sum(le_this$death_role == "anti-death")
    n_amb  <- sum(le_this$death_role %in% c("ambiguous", "unclassified"))
    n_total <- nrow(le_this)

    # Determine direction
    if (n_pro + n_anti == 0) {
      verdict <- "Indeterminate (no classified LE genes)"
    } else if (nes > 0 && n_pro > n_anti) {
      verdict <- "MORE SUSCEPTIBLE at 12W (pro-death genes upregulated)"
    } else if (nes > 0 && n_anti > n_pro) {
      verdict <- "MORE RESISTANT at 12W (anti-death genes upregulated)"
    } else if (nes > 0 && n_pro == n_anti) {
      verdict <- "MIXED (equal pro/anti upregulated)"
    } else if (nes < 0 && n_pro > n_anti) {
      verdict <- "MORE RESISTANT at 12W (pro-death genes downregulated)"
    } else if (nes < 0 && n_anti > n_pro) {
      verdict <- "MORE SUSCEPTIBLE at 12W (anti-death genes downregulated)"
    } else if (nes < 0 && n_pro == n_anti) {
      verdict <- "MIXED (equal pro/anti downregulated)"
    } else {
      verdict <- "Indeterminate"
    }

    sig_flag <- if (p < 0.05) "***" else if (p < 0.10) "*" else ""
    message(sprintf("  %-30s NES=%6.2f padj=%.3f %s | LE: %d pro, %d anti, %d other | %s",
                    pw_clean, nes, p, sig_flag, n_pro, n_anti, n_amb, verdict))

    interpretation[[paste(ctr, pw, sep = "::")]] <- tibble(
      contrast   = ctr,
      genotype   = ctr_label,
      pathway    = pw_clean,
      NES        = nes,
      padj       = p,
      n_le_total = n_total,
      n_pro      = n_pro,
      n_anti     = n_anti,
      n_ambiguous = n_amb,
      verdict    = verdict
    )
  }
}

interpretation_df <- bind_rows(interpretation)

# =============================================================================
# PART 10: SUMMARY VISUALISATION — VERDICT HEATMAP
# =============================================================================

message("\nGenerating interpretation summary plot...")

# Create a numeric score: +1 = more susceptible, -1 = more resistant, 0 = mixed
interpretation_df <- interpretation_df %>%
  mutate(
    direction_score = case_when(
      grepl("MORE SUSCEPTIBLE", verdict) ~  1,
      grepl("MORE RESISTANT", verdict)   ~ -1,
      TRUE                               ~  0
    ),
    sig_label = case_when(
      padj < 0.01  ~ "**",
      padj < 0.05  ~ "*",
      padj < 0.10  ~ ".",
      TRUE         ~ ""
    )
  )

# Wide format for heatmap
verdict_wide <- interpretation_df %>%
  dplyr::select(pathway, genotype, direction_score) %>%
  pivot_wider(names_from = genotype, values_from = direction_score) %>%
  column_to_rownames("pathway")

sig_wide <- interpretation_df %>%
  dplyr::select(pathway, genotype, sig_label) %>%
  pivot_wider(names_from = genotype, values_from = sig_label) %>%
  column_to_rownames("pathway")

verdict_mat <- as.matrix(verdict_wide)
sig_label_mat <- as.matrix(sig_wide)

col_verdict <- circlize::colorRamp2(
  c(-1, 0, 1),
  c("#2166AC", "grey95", "#B2182B")  # blue = resistant, red = susceptible
)

pdf(file.path(cd_dir, "heatmap_cell_death_verdict.pdf"),
    width = 7, height = 7)
ht_verdict <- ComplexHeatmap::Heatmap(
  verdict_mat,
  name = "Direction",
  col  = col_verdict,
  cluster_columns = FALSE,
  cluster_rows    = TRUE,
  row_names_side  = "left",
  row_names_gp    = gpar(fontsize = 9),
  column_names_gp = gpar(fontsize = 11),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sig_label_mat[i, j], x, y, gp = gpar(fontsize = 9))
  },
  column_title = "Cell Death Pathway Direction at 12W vs 6W",
  column_title_gp = gpar(fontsize = 12, fontface = "bold"),
  heatmap_legend_param = list(
    title = "Direction",
    at = c(-1, 0, 1),
    labels = c("More resistant", "Mixed/unclear", "More susceptible"),
    title_gp = gpar(fontsize = 10),
    labels_gp = gpar(fontsize = 9)
  )
)
draw(ht_verdict)
dev.off()

# =============================================================================
# PART 11: EXPORT LEADING EDGE DETAILS FOR SIGNIFICANT PATHWAYS
# =============================================================================

message("Exporting results...")

# --- 11a. Full fGSEA results table ---
fgsea_export <- fgsea_combined %>%
  dplyr::select(pathway, pathway_label, contrast, NES, pval, padj, size) %>%
  arrange(contrast, pval)

write.csv(fgsea_export,
          file.path(cd_dir, "fgsea_cell_death_all_results.csv"),
          row.names = FALSE)

# --- 11b. Leading edge with annotations (significant pathways only) ---
le_export <- le_combined %>%
  filter(padj < 0.10) %>%
  dplyr::select(contrast, pathway, NES, padj, mouse_symbol, human_gene,
                description, comment, death_role, le_size, set_size) %>%
  arrange(contrast, pathway, death_role, mouse_symbol)

write.csv(le_export,
          file.path(cd_dir, "leading_edge_annotated_sig_pathways.csv"),
          row.names = FALSE)

# --- 11c. Interpretation summary ---
write.csv(interpretation_df,
          file.path(cd_dir, "cell_death_interpretation_summary.csv"),
          row.names = FALSE)

# =============================================================================
# PART 12: SAVE RDS
# =============================================================================

message("Saving RDS...")

cell_death_results <- list(
  # Gene sets
  cd_human_sets  = cd_human_sets,
  cd_mouse_sets  = cd_mouse_sets,
  cd_descriptions = cd_descriptions,

  # fGSEA results
  fgsea_results  = fgsea_results,
  fgsea_combined = fgsea_combined,

  # NES matrices
  nes_wide  = nes_wide,
  padj_wide = padj_wide,

  # Leading edge
  le_combined = le_combined,

  # Interpretation
  interpretation = interpretation_df,

  # Metadata
  analysis_date = Sys.Date(),
  reference = "Tang et al. (2024) Comput Struct Biotechnol J. doi:10.1016/j.csbj.2024.08.012",
  description = paste(
    "fGSEA analysis of 15 regulated cell death modalities across 5 DESeq2",
    "contrasts (Myc effect at 6W/12W, temporal in Myc-/Myc+, interaction).",
    "Gene sets from Tang et al. (2024), converted from human to mouse orthologs.",
    "Leading edge genes annotated with pro-/anti-death roles from source",
    "descriptions. Interpretation: whether 12W samples are more resistant or",
    "susceptible to each death modality."
  )
)

saveRDS(cell_death_results, here("results", "cell_death_fgsea.rds"))

# =============================================================================
# PART 13: CONSOLE SUMMARY
# =============================================================================

message("\n", strrep("=", 70))
message("CELL DEATH PATHWAY ANALYSIS COMPLETE")
message(strrep("=", 70))

message("\n--- Significant pathways per contrast (padj < 0.05) ---")
for (nm in names(fgsea_results)) {
  sig <- fgsea_results[[nm]] %>% filter(padj < 0.05) %>% arrange(pval)
  if (nrow(sig) > 0) {
    message(sprintf("\n  %s:", nm))
    for (j in seq_len(nrow(sig))) {
      message(sprintf("    %-35s NES=%6.2f  padj=%.4f",
                      str_replace(sig$pathway[j], "^CD_", ""),
                      sig$NES[j], sig$padj[j]))
    }
  } else {
    message(sprintf("\n  %s: (none significant)", nm))
  }
}

message("\n--- Temporal interpretation (significant at padj < 0.10) ---")
interpretation_df %>%
  filter(padj < 0.10) %>%
  arrange(genotype, padj) %>%
  mutate(summary = sprintf("%-8s %-25s NES=%5.2f padj=%.3f → %s",
                           genotype, pathway, NES, padj, verdict)) %>%
  pull(summary) %>%
  { message(paste("  ", ., collapse = "\n")) }

message(sprintf("\nResults saved to: results/cell_death_fgsea.rds"))
message(sprintf("Figures saved to: %s", cd_dir))
message(strrep("=", 70))
