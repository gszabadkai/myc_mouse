# scripts/07_heatmap_oxphos_annotated.R

source("scripts/00_setup_packages.R")
source("functions/generate_heatmap.R")

# === Load inputs ===
dds_int <- readRDS("results/dds_int_run.rds")
gene_sets_list <- readRDS("results/gene_sets_list.rds")

set_name <- "MC_OXPHOS > OXPHOS subunits"
gene_symbols <- gene_sets_list[[set_name]]

# === make OXPHOS annotation ===
oxphos_df <- read.csv("data/OXPHOS_subunits.csv", header = FALSE, col.names = c("Complex", "Genes"))
oxphos_long <- oxphos_df %>%
  mutate(Genes = strsplit(as.character(Genes), ",")) %>%
  tidyr::unnest(Genes) %>%
  mutate(Genes = str_trim(Genes))

oxphos_long$mtDNA <- ifelse(str_starts(oxphos_long$Genes, "mt-"), "mtDNA", "nuclear")

# === Generate heatmap (returns heatmap object) ===

for (lfc_variant in c("shrunk", "raw")) {
  # Load correct LFC version
  lfc_df <- switch(
    lfc_variant,
    shrunk = readRDS("results/combined_df_annotated.rds"),
    raw    = readRDS("results/combined_df_annotated_raw.rds")
  )
  ht_result <- generate_heatmaps_for_gene_set(
    set_name = set_name,
    gene_symbols = gene_symbols,
    dds_obj = dds_int,
    lfc_variant = lfc_variant,
    show_plot = FALSE
  )
  # === Extract row order and heatmap object ===
  ht_main <- ht_result$heatmap_object
  lfc_mat <- ht_result$lfc_matrix
  rownames_lfc <- rownames(lfc_mat)
  # Reorder annotation to match heatmap gene order
  oxphos_anno_df <- oxphos_long %>%
    distinct(Genes, Complex, mtDNA) %>%
    filter(Genes %in% rownames_lfc) %>%
    column_to_rownames("Genes") %>%
    .[rownames_lfc, , drop = FALSE]
  # === Custom row annotations (to avoid overwriting the main function's annotation need to add it again as raw_ha) ===
  row_ha_oxphos <- rowAnnotation(
    Complex = oxphos_anno_df$Complex,
    mtDNA = oxphos_anno_df$mtDNA,
    col = list(
      Complex = structure(
        RColorBrewer::brewer.pal(length(unique(oxphos_anno_df$Complex)), "Blues"),
        names = unique(oxphos_anno_df$Complex)
      ),
      mtDNA = c(mtDNA = "black", nuclear = "grey90")
    ),
    show_annotation_name = TRUE
  )
  
  row_meta <- lfc_df %>%
    filter(mgi_symbol %in% rownames(lfc_mat)) %>%
    dplyr::select(mgi_symbol, starts_with("time_effect_category"), group_sig_status) %>%
    distinct() %>%
    column_to_rownames("mgi_symbol") %>%
    .[rownames(lfc_mat), , drop = FALSE]
  
  # Define colors
  time_colors <- c(
    "direct_Myc_reduction" = "firebrick",
    "baseline_driven_reduction" = "purple",
    "baseline_driven_increase" = "steelblue",
    "direct_Myc_increase" = "darkgreen",
    "no_change" = "grey90"
  )
  
  sig_colors <- c(
    up = "blue",
    down = "orange",
    ns = "grey90"
  )
  
  # Build row annotation
  ra_list <- list()
  time_col <- grep("time_effect_category", colnames(row_meta), value = TRUE)
  if (length(time_col) > 0) {
    ra_list$Time_Effect <- row_meta[[time_col]]
  }
  if ("group_sig_status" %in% colnames(row_meta)) {
    ra_list$Group_Change <- row_meta$group_sig_status
  }
  
  row_ha <- rowAnnotation(
    df = as.data.frame(ra_list),
    col = list(
      Time_Effect = time_colors,
      Group_Change = sig_colors
    ),
    show_annotation_name = TRUE
  )
  
  
  # === Save combined heatmap ===
  pdf_file <- sprintf("outputs/heatmaps_int/OXPHOS_%s_with_extra_annotation2.pdf", lfc_variant)
  pdf(pdf_file, width = 12, height = 15)
  grid.newpage()
  draw(row_ha + row_ha_oxphos + ht_main, merge_legend = TRUE)
  dev.off()
}

