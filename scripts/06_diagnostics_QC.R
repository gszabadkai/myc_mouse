# scripts/06_diagnostics_QC.R

source("scripts/00_setup_packages.R")

# === Load inputs ===
dds <- readRDS("results/dds_int_run.rds")
cts <- readRDS("results/count_matrix.rds")
coldata <- readRDS("results/coldata.rds")

# Ensure output directory exists
qc_dir <- "outputs/qc"
dir.create(qc_dir, showWarnings = FALSE, recursive = TRUE)

# === Match samples between coldata and cts ===
stopifnot(all(rownames(coldata) %in% colnames(cts)))

# === Total counts per sample ===
coldata$tcts <- round(colSums(cts)[rownames(coldata)] / 1e6, 1)

coldata$Exp_groups <- factor(
  paste(coldata$myc_status, coldata$timepoint, sep = "_"),
  levels = c("neg_6W", "pos_6W", "neg_12W", "pos_12W")
)

# === Plot total count distribution ===
ggsave(
  file.path(qc_dir, "total_counts_by_group.pdf"),
  ggbetweenstats(
    data = coldata,
    x = Exp_groups,
    y = tcts,
    type = "parametric",
    pairwise.comparisons = TRUE,
    pairwise.display = "significant",
    title = "Total counts per sample",
    xlab = "MYC status ~ timepoint",
    ylab = "Total counts (M)"
  ),
  width = 6, height = 4
)

# === Size factors ===
dds <- estimateSizeFactors(dds)
coldata$size_factors <- sizeFactors(dds)

ggsave(
  file.path(qc_dir, "size_factors_by_group.pdf"),
  ggbetweenstats(
    data = coldata,
    x = Exp_groups,
    y = size_factors,
    type = "parametric",
    pairwise.comparisons = TRUE,
    pairwise.display = "significant",
    title = "Size factors",
    xlab = "MYC status ~ timepoint",
    ylab = "Size factor"
  ),
  width = 6, height = 4
)

# === VST and rlog transforms ===
vsd <- vst(dds, blind = FALSE)
rld <- rlog(dds, blind = FALSE)

# === Compare variance transforms ===
df_var <- bind_rows(
  as.data.frame(log2(counts(dds, normalized = TRUE)[, 1:2] + 1)) %>%
    mutate(transformation = "log2(x + 1)"),
  as.data.frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"),
  as.data.frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog")
)

colnames(df_var)[1:2] <- c("x", "y")
df_var$transformation <- factor(df_var$transformation, levels = c("log2(x + 1)", "vst", "rlog"))

ggsave(
  file.path(qc_dir, "variance_transforms.pdf"),
  ggplot(df_var, aes(x = x, y = y)) +
    geom_hex(bins = 80) +
    coord_fixed() +
    facet_grid(. ~ transformation) +
    labs(title = "Variance stabilization across transformations"),
  width = 9, height = 3
)

# === Sample distance heatmaps ===
colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Blues")))(255)

sample_dists <- function(mat, label_prefix, outname) {
  dmat <- dist(t(mat))
  dmat_mtx <- as.matrix(dmat)
  rownames(dmat_mtx) <- paste(coldata$group, rownames(coldata), sep = " / ")
  colnames(dmat_mtx) <- NULL
  
  pdf(file.path(qc_dir, outname), width = 6, height = 6)
  pheatmap(dmat_mtx,
           clustering_distance_rows = dmat,
           clustering_distance_cols = dmat,
           col = colors,
           main = paste("Sample distances by", label_prefix))
  dev.off()
}

sample_dists(assay(vsd), "VST", "sample_distance_vsd.pdf")
sample_dists(assay(rld), "rlog", "sample_distance_rlog.pdf")

# === PCA plot ===
ggsave(
  file.path(qc_dir, "PCA_vsd.pdf"),
  plotPCA(vsd, intgroup = c("group")),
  width = 5, height = 4
)

# === Optional: Poisson distance (if installed) ===
if (requireNamespace("PoiClaClu", quietly = TRUE)) {
  poisd <- PoiClaClu::PoissonDistance(t(counts(dds)))
  mtx <- as.matrix(poisd$dd)
  rownames(mtx) <- paste(coldata$group, rownames(coldata), sep = " / ")
  colnames(mtx) <- NULL
  
  pdf(file.path(qc_dir, "sample_distance_poisson.pdf"), width = 6, height = 6)
  pheatmap(mtx,
           clustering_distance_rows = poisd$dd,
           clustering_distance_cols = poisd$dd,
           col = colors,
           main = "Sample distances by Poisson")
  dev.off()
}
