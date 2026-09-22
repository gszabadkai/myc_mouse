# scripts/05_heatmap_utils.R

source("scripts/00_setup_packages.R")
source("functions/generate_heatmap.R")

# === Load data ===
gene_sets_list <- readRDS("results/gene_sets_list.rds")
dds_int <- readRDS("results/dds_int_run.rds")

# === Create output folder ===
output_dir <- "outputs/heatmaps_int/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# === Generate heatmaps for both LFC variants ===
heatmap_outputs <- list()

for (variant in c("shrunk", "raw")) {
  message("Generating heatmaps for: ", variant)
  heatmap_outputs[[variant]] <- lapply(names(gene_sets_list), function(set_name) {
    message(" → ", set_name)
    generate_heatmaps_for_gene_set(
      set_name = set_name,
      gene_symbols = gene_sets_list[[set_name]],
      dds_obj = dds_int,
      lfc_variant = variant,
      output_dir = output_dir,
      show_plot = TRUE
    )
  })
  names(heatmap_outputs[[variant]]) <- names(gene_sets_list)
}


# Save file paths for reference
saveRDS(heatmap_outputs, file = "results/heatmap_paths_by_variant.rds")
