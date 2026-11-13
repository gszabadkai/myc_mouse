# ----------------------------
# Setup
# ----------------------------
##### Remove all environmental variables

rm(list = ls())
setwd('/Users/kossifm/Documents/Bulkseq_New_Amy')

library(DESeq2)
library(pheatmap)
library(dplyr)
library(tibble)
library(ggplot2)
library(ggvenn)
library(RColorBrewer)
library(ggrepel)
library(msigdbr)
library(AnnotationHub)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(clusterProfiler)
library(enrichplot)
library(ggnewscale)
library(cowplot)
library(fgsea)
library(data.table)
library(ggtext)
library(decoupleR)
library(OmnipathR)


# ----------------------------
# 1) Read input counts
# ----------------------------
counts_all <- read.table("salmon.merged.gene_counts.tsv",
                         header = TRUE, row.names = 1, check.names = FALSE, sep = "\t")

if ("gene_name" %in% colnames(counts_all)) {
  countmat <- counts_all[ , !(colnames(counts_all) %in% "gene_name") ]
} else {
  countmat <- counts_all
}

sample_names <- colnames(countmat)

# ----------------------------
# 2) Extract conditions (robust)
# ----------------------------
tmp <- sub("^[^-]+-", "", sample_names)
tmp <- sub("BRNS.*$", "", tmp)
tmp <- sub("_.*$", "", tmp)
condition_raw <- sub("([A-Za-z]+).*", "\\1", tmp)

condition <- dplyr::case_when(
  condition_raw %in% c("BoMY", "BoMYBFGE") ~ "BoMY1",   # <- merged here
  condition_raw == "FaMY"    ~ "FaMY1",
  condition_raw == "Pgc"     ~ "Pgc1a",
  condition_raw == "KOPgc"   ~ "KOPgc1a",
  condition_raw == "NTPgc"   ~ "NTPgc1a",
  TRUE ~ condition_raw
)

expected_conditions <- c(
  "BclxL", "BoMY1", "EV", "FaMY1",     # <- "BoMYBFGE" removed, it can never appear now
  "KOEV", "KOPgc1a", "MYAZ", "NTEV",
  "NTPgc1a", "NTPgc", "Pgc1a"
)

if (!all(condition %in% expected_conditions)) {
  stop(
    "Unexpected condition labels found:\n",
    paste(setdiff(unique(condition), expected_conditions), collapse = ", ")
  )
}
conds <- factor(condition, levels = expected_conditions)
colData <- data.frame(row.names = sample_names, condition = conds)

# ----------------------------
# 3) Build DESeq2 object
# ----------------------------
dds <- DESeqDataSetFromMatrix(countData = round(countmat),
                              colData   = colData,
                              design    = ~ condition)

keep <- rowSums(counts(dds) >= 10) >= 2
dds <- dds[keep, ]
dds <- DESeq(dds)

# ----------------------------
# 4) Gene set collections
# ----------------------------
h_gene_sets  <- msigdbr(species = "mouse", category = "H")
c3_gene_sets  <- msigdbr(species = "mouse", category = "C3")
cpg_gene_sets <- msigdbr(species = "mouse", category = "C2")
c5_gene_sets  <- msigdbr(species = "mouse", category = "C5")
c6_gene_sets  <- msigdbr(species = "mouse", category = "C6")

msigdbr_h_list  <- split(h_gene_sets$entrez_gene,  f = h_gene_sets$gs_name)
msigdbr_c3_list  <- split(c3_gene_sets$entrez_gene,  f = c3_gene_sets$gs_name)
msigdbr_cpg_list <- split(cpg_gene_sets$entrez_gene, f = cpg_gene_sets$gs_name)
msigdbr_c5_list  <- split(c5_gene_sets$entrez_gene,  f = c5_gene_sets$gs_name)
msigdbr_c6_list  <- split(c6_gene_sets$entrez_gene,  f = c6_gene_sets$gs_name)

# ----------------------------
# 5) Function for DESeq + annotation + GSEA with deduplication
# ----------------------------
run_deseq_gsea <- function(dds, num, den, prefix) {
  cond_levels <- levels(dds$condition)
  if (!(num %in% cond_levels)) stop("Numerator condition '", num, "' not found.")
  if (!(den %in% cond_levels)) stop("Denominator condition '", den, "' not found.")
  
  message("Running DESeq2 + GSEA for: ", prefix)
  
  # DESeq2 results using contrast
  res <- results(dds, contrast = c("condition", num, den))
  res <- na.omit(res)
  res <- res[order(res$log2FoldChange, decreasing = TRUE), ]
  res <- as.data.frame(res) %>% rownames_to_column("ENSEMBL")
  
  # Annotation
  res_anno <- AnnotationDbi::select(org.Mm.eg.db, keys = res$ENSEMBL,
                                    columns = c("ENSEMBL","SYMBOL","GENENAME","ENTREZID"),
                                    keytype = "ENSEMBL") %>%
    filter(!duplicated(ENSEMBL))
  
  res_annotated <- left_join(res_anno, res, by="ENSEMBL") %>% na.omit()
  res_annotated <- res_annotated[order(-res_annotated$log2FoldChange), ]
  write.csv(res_annotated, paste0(prefix, "_annotated_genelist.csv"), row.names = FALSE)
  
  # ----------------------------
  # Deduplicate Entrez IDs for FGSEA
  # ----------------------------
  res_ENTREZ <- res_annotated %>%
    group_by(ENTREZID) %>%
    summarize(log2FC = log2FoldChange[which.max(abs(log2FoldChange))]) %>%
    ungroup()
  
  res_vect <- res_ENTREZ$log2FC
  res_vect <- res_vect[is.finite(res_vect)]
  print (res_vect)
  
  names(res_vect) <- res_ENTREZ$ENTREZID
  
  # Run fgsea for all collections
  gs_lists <- list(H = msigdbr_h_list, C2 = msigdbr_cpg_list, C5 = msigdbr_c5_list, C6 = msigdbr_c6_list, C3 = msigdbr_c3_list)
  
  for (gs in names(gs_lists)) {
    res_gsea <- fgseaMultilevel(pathways = gs_lists[[gs]],
                                stats = res_vect,
                                minSize = 15,
                                maxSize = 1000)
    
    fwrite(res_gsea, paste0(prefix, "_MSigDB_", gs, "_full.csv"))
    fwrite(res_gsea[res_gsea$padj < 0.01, ], paste0(prefix, "_MSigDB_", gs, "_significant.csv"))
    
    # Plot Hallmark only
    if (gs == "H") {
      res_filt <- res_gsea[res_gsea$padj < 0.05, ]
      if (nrow(res_filt) > 0) {
        p <- ggplot(res_filt, aes(x = reorder(pathway, NES), y = NES)) +
          geom_col(aes(fill = as.factor(sign(NES)))) +
          scale_fill_manual(values = c("red", "blue"),
                            labels = c("Underrepresented", "Overrepresented")) +
          coord_flip() +
          labs(x = "Pathway", y = "Normalized Enrichment Score", fill = "NES Sign") +
          ggtitle(paste0(num, " vs ", den, " — MSigDB Hallmark pathways")) +
          theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
                panel.grid.major.x = element_blank(),
                panel.grid.minor.x = element_blank(),
                panel.background = element_blank())
        
        ggsave(paste0(prefix, "_MSigDB_Hallmarks_Barplot.tiff"),
               plot = p, width = 12, height = 6, dpi = 800)
      }
    }
  }
  
  invisible(res_annotated)
}

# ----------------------------
# 6) Run target contrasts
# ----------------------------
contrasts_to_run <- list(
  list(num = "BoMY1",   den = "FaMY1",    prefix = "BoMY_vs_FaMY"),
  list(num = "BoMY1",   den = "MYAZ",     prefix = "BoMY_vs_MYAZ"),
  list(num = "FaMY1",   den = "MYAZ",     prefix = "FaMY_vs_MYAZ"),
  list(num = "Pgc1a",   den = "EV",       prefix = "Pgc_vs_EV"),
  list(num = "Pgc1a",   den = "BclxL",    prefix = "Pgc_vs_BclxL"),
  list(num = "KOEV",    den = "NTEV",     prefix = "KOEV_vs_NTEV"),
  list(num = "NTPgc1a", den = "NTEV",     prefix = "NTPgc_vs_NTEV"),
  list(num = "BclxL", den = "EV",     prefix = "BclxL_vs_EV"),
  list(num = "KOPgc1a", den = "NTPgc1a",     prefix = "KOPgc1a_vs_NTPgc1a"),
  list(num = "KOPgc1a", den = "NTEV",     prefix = "KOPgc1a_vs_NTEV")
)

for (c in contrasts_to_run) {
  run_deseq_gsea(dds, c$num, c$den, c$prefix)
}


##############################################################################
##### Running TF analysis from online vignette (looped, saved outputs) #######
##############################################################################

## ===============================
## Load CollecTRI TF network (once)
## ===============================
net <- get_collectri(
  organism = "mouse",
  split_complexes = FALSE
)

## ===============================
## Define comparisons (derived from contrasts_to_run so every
## DESeq2 contrast above automatically gets a TF analysis)
## ===============================
prefixes <- sapply(contrasts_to_run, `[[`, "prefix")
nums     <- sapply(contrasts_to_run, `[[`, "num")
dens     <- sapply(contrasts_to_run, `[[`, "den")

comparisons <- tibble::tibble(
  comparison = prefixes,
  file       = paste0(prefixes, "_annotated_genelist.csv"),
  ylab       = paste0("Activity score (", nums, " vs ", dens, ")")
)

n_tfs <- 25

## ===============================
## Create output directories
## ===============================
dir.create("TF_results", showWarnings = FALSE)
dir.create("TF_plots", showWarnings = FALSE)

## ===============================
## Loop over comparisons
## ===============================
for (i in seq_len(nrow(comparisons))) {
  
  message("Running TF analysis: ", comparisons$comparison[i])
  
  ## -------------------------------
  ## Read DESeq2 results
  ## -------------------------------
  deg_raw <- read.csv(
    comparisons$file[i],
    stringsAsFactors = FALSE
  )
  
  ## -------------------------------
  ## Build gene-level statistic
  ## -------------------------------
  deg <- deg_raw %>%
    select(SYMBOL, stat) %>%
    filter(!is.na(stat)) %>%
    group_by(SYMBOL) %>%
    slice_max(order_by = abs(stat), n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    column_to_rownames("SYMBOL") %>%
    as.matrix()
  
  ## -------------------------------
  ## Infer TF activities (ULM)
  ## -------------------------------
  tf_acts <- run_ulm(
    mat = deg,
    net = net,
    .source = "source",
    .target = "target",
    .mor = "mor",
    minsize = 5
  )
  
  ## -------------------------------
  ## Save full TF results
  ## -------------------------------
  write.csv(
    tf_acts,
    file = file.path(
      "TF_results",
      paste0(comparisons$comparison[i], "_TF_activity_all.csv")
    ),
    row.names = FALSE
  )
  
  ## -------------------------------
  ## Select top TFs
  ## -------------------------------
  tf_top <- tf_acts %>%
    mutate(
      rnk = ifelse(
        score > 0,
        rank(-score),
        rank(-abs(score))
      )
    ) %>%
    arrange(rnk) %>%
    head(n_tfs)
  
  ## -------------------------------
  ## Save top TFs
  ## -------------------------------
  write.csv(
    tf_top,
    file = file.path(
      "TF_results",
      paste0(comparisons$comparison[i], "_TF_activity_top25.csv")
    ),
    row.names = FALSE
  )
  
  ## -------------------------------
  ## Plot TF activity
  ## -------------------------------
  p <- ggplot(
    tf_top,
    aes(
      x = reorder(source, score),
      y = score,
      fill = score
    )
  ) +
    geom_bar(stat = "identity", color = "black") +
    scale_fill_gradient2(
      low = brewer.pal(11, "RdBu")[10],
      mid = "white",
      high = brewer.pal(11, "RdBu")[2],
      midpoint = 0
    ) +
    coord_flip() +
    theme_minimal(base_size = 12) +
    labs(
      title = paste(
        "Transcription factor activity (",
        gsub("_", " ", comparisons$comparison[i]),
        ")",
        sep = ""
      ),
      x = "Transcription Factor",
      y = comparisons$ylab[i]
    )
  
  ## -------------------------------
  ## Save plot
  ## -------------------------------
  ggsave(
    filename = file.path(
      "TF_plots",
      paste0(comparisons$comparison[i], "_TF_activity_top25.png")
    ),
    plot = p,
    width = 8,
    height = 6,
    dpi = 300
  )
  
  print(p)
}



##########################################################################
#########################################################################
##### Plot MMP plots and integrins

df <- read.csv("Pgc_vs_EV_annotated_genelist.csv")

mmp_df <- df[grepl("Mmp", df$SYMBOL), ]
mmp_df <- mmp_df[order(mmp_df$log2FoldChange), ]

ggplot(mmp_df, aes(x = log2FoldChange, y = SYMBOL)) +
  geom_col() +
  labs(title = "Differential Expression of MMPs",
       x = "log2 Fold Change",
       y = "MMP Genes") +
  theme_minimal()


# Mouse gene symbols (correct capitalization)
genes_of_interest <- c(
  "Itga3", "Itga4", "Itga10", "Itgb3", "Itgb5",
  "Cav1", "Cav2", "Actn2", "Lama4", "Col4a1",
  "Inhba", "Tgfbi", "Tgfbr3", "Tgfbr2",
  "Smad3", "Smad7", "Il8", "Il11",
  "Lef1", "Tcf7l2", "Dkk3", "Ppp3ca", "Sfrp1"
)
genes_of_interest <- c(
  "Itga5","Itga3", "Itga4", "Itga10", "Itgb3", "Itgb5","Itgb1",
  "Cav1", "Cav2", "Actn2", "Lama4", "Col4a1",
  "Inhba"
)

# Subset genes
sub_df <- df[df$SYMBOL %in% genes_of_interest, ]

# Order by log2 fold change
sub_df <- sub_df[order(sub_df$log2FoldChange), ]

# Dot plot
ggplot(sub_df, aes(x = log2FoldChange, y = SYMBOL)) +
  geom_point(size = 3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  labs(
    title = "Differential Expression of Selected Genes",
    x = "log2 Fold Change",
    y = "Gene"
  ) +
  theme_minimal()


# ==========================================
# PGC vs EV: TF activity & MMP-specific analysis (Mouse)
# ==========================================

# ------------------------------
# Load libraries
# ------------------------------
library(tidyverse)    # Data manipulation + ggplot2
library(dorothea)     # TF regulons
library(decoupleR)    # TF activity inference
library(org.Mm.eg.db) # Mouse gene annotation
library(AnnotationDbi)

# ------------------------------
# 1. Load Dorothea mouse regulons (A/B/C)
# ------------------------------
dorothea_mm <- dorothea::dorothea_mm
dorothea_df <- dorothea_mm %>%
  filter(confidence %in% c("A","B","C")) %>%
  dplyr::rename(tf = tf) %>%
  dplyr::select(target, tf, mor) %>%
  as.data.frame()
dorothea_df$likelihood <- 1

# ------------------------------
# 2. Load RNA-seq DE results
# ------------------------------
RNA_diff <- read_csv("Pgc_vs_BclxL_annotated_genelist.csv")

# Collapse duplicate SYMBOLs by mean log2FC
RNA_diff_collapsed <- RNA_diff %>%
  group_by(SYMBOL) %>%
  summarise(log2FoldChange = mean(log2FoldChange, na.rm = TRUE)) %>%
  ungroup()

# Convert to matrix with unique rownames
RNA_diff_mat <- RNA_diff_collapsed %>%
  tibble::column_to_rownames("SYMBOL") %>%
  dplyr::select(log2FoldChange) %>%
  as.matrix()

# ------------------------------
# 3. Map genes to Dorothea targets
# ------------------------------
mapped_genes <- mapIds(org.Mm.eg.db,
                       keys = rownames(RNA_diff_mat),
                       column = "SYMBOL",
                       keytype = "SYMBOL",
                       multiVals = "first")

RNA_diff_mat <- RNA_diff_mat[!is.na(mapped_genes), , drop=FALSE]
rownames(RNA_diff_mat) <- mapped_genes[!is.na(mapped_genes)]

# ------------------------------
# 4. Visualize top 10 DE genes
# ------------------------------
top_genes <- head(RNA_diff_mat[order(abs(RNA_diff_mat[,1]), decreasing = TRUE), , drop=FALSE], 10)
top_genes_df <- data.frame(
  gene = rownames(top_genes),
  log2FC = top_genes[,1]
)

ggplot(top_genes_df, aes(x = log2FC, y = reorder(gene, log2FC))) +
  geom_point(size = 3) +
  labs(title = "Top 10 Differentially Expressed Genes", x = "log2FC", y = "Gene") +
  theme_minimal()

# ------------------------------
# 5. Run TF activity inference using decoupleR
# ------------------------------
TF_activities <- run_wmean(
  mat = RNA_diff_mat,
  network = dorothea_df,
  .source = "tf",
  times = 1000
)

TF_activities <- as.data.frame(TF_activities) %>%
  filter(statistic == "norm_wmean") %>%
  select(source, score) %>%
  tibble::column_to_rownames("source")

# ------------------------------
# 6. Plot top 10 TFs overall
# ------------------------------
top_TFs <- head(TF_activities[order(abs(TF_activities$score), decreasing = TRUE), , drop=FALSE], 10)
top_TFs_df <- data.frame(TF = rownames(top_TFs), Activity = top_TFs$score)

ggplot(top_TFs_df, aes(x = Activity, y = reorder(TF, Activity))) +
  geom_point(size = 3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  labs(title = "Top 10 TF Activities (PGC vs EV)", x = "Activity score", y = "TF") +
  theme_minimal()

# ------------------------------
# 7. Focus on downregulated MMP genes
# ------------------------------
MMP_genes <- c("Mmp1a","Mmp2","Mmp3","Mmp8","Mmp9","Mmp10","Mmp12","Mmp13","Mmp14")
MMP_in_data <- intersect(rownames(RNA_diff_mat), MMP_genes)

MMP_targets <- dorothea_df %>%
  filter(target %in% MMP_in_data)

TF_MMP_links <- MMP_targets %>%
  group_by(tf) %>%
  summarize(n_targets = n(), .groups = "drop") %>%
  mutate(activity = TF_activities[tf, "score"]) %>%
  arrange(activity)

# Print candidate TFs regulating MMPs
print(TF_MMP_links)

# ------------------------------
# 8. Plot TFs regulating MMPs
# ------------------------------


ggplot(TF_MMP_links, aes(x = activity, y = reorder(tf, activity))) +
  geom_point(size = 3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  labs(title = "Candidate TFs regulating downregulated MMPs",
       x = "TF activity score",
       y = "TF") +
  theme_minimal()

# ------------------------------
# 9. Optional: MMP log2FC per TF heatmap
# ------------------------------
if(nrow(TF_MMP_links) > 0 & length(MMP_in_data) > 0){
  TF_MMP_mat <- sapply(TF_MMP_links$tf, function(tf){
    targets <- MMP_targets$target[MMP_targets$tf == tf]
    mean(RNA_diff_mat[targets, 1], na.rm = TRUE)
  })
  TF_MMP_mat <- matrix(TF_MMP_mat, nrow = length(TF_MMP_mat),
                       dimnames = list(TF_MMP_links$tf, "avg_MMP_log2FC"))
  
  ggplot(as.data.frame(TF_MMP_mat) %>% tibble::rownames_to_column("TF"),
         aes(x = TF, y = avg_MMP_log2FC, fill = avg_MMP_log2FC)) +
    geom_bar(stat = "identity") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Average log2FC of MMPs per TF")
}


TF_MMP_mat <- matrix(
  TF_MMP_mat,
  nrow = length(TF_MMP_mat),
  dimnames = list(TF_MMP_links$tf, "avg_MMP_log2FC")
)

ggplot(as.data.frame(TF_MMP_mat) %>% tibble::rownames_to_column("TF"),
       aes(x = TF, y = avg_MMP_log2FC, fill = avg_MMP_log2FC)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Average log2FC of MMPs per TF")


####. MMP-consistency score per TF (direction-aware)
# ------------------------------
TF_MMP_score <- sapply(TF_MMP_links$tf, function(tf) {
  targets <- MMP_targets %>% filter(tf == !!tf)
  mean(
    RNA_diff_mat[targets$target, 1] * targets$mor,
    na.rm = TRUE
  )
})

TF_MMP_score_df <- data.frame(
  TF = names(TF_MMP_score),
  MMP_directional_score = TF_MMP_score
)

ggplot(TF_MMP_score_df,
       aes(x = reorder(TF, MMP_directional_score),
           y = MMP_directional_score,
           fill = MMP_directional_score)) +
  geom_col() +
  coord_flip() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(
    title = "Direction-consistent TF influence on MMPs",
    x = "Transcription factor",
    y = "MMP directional score"
  ) +
  theme_minimal()



