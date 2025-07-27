# scripts/02_deseq_interaction_model.R

source("scripts/00_setup_packages.R")

# === Load processed data ===
dds_int <- readRDS("results/dds_int.rds")

# === Run DESeq2 ===
dds_int <- DESeq(dds_int)

# Save DESeq2 object
saveRDS(dds_int, "results/dds_int_run.rds")

# === Extract contrast results (IHW used for filtering) ===

# 1. Raw (unshrunken, IHW padj)
myc_6W_log2FC_raw         <- results(dds_int, name = "myc_status_pos_vs_neg", filterFun = ihw)
timepoint_neg_log2FC_raw <- results(dds_int, name = "timepoint_12W_vs_6W", filterFun = ihw)
myc_12W_log2FC_raw        <- results(dds_int, contrast = list(c("myc_status_pos_vs_neg", "timepoint12W.myc_statuspos")), filterFun = ihw)
timepoint_pos_log2FC_raw <- results(dds_int, contrast = list(c("timepoint_12W_vs_6W", "timepoint12W.myc_statuspos")), filterFun = ihw)

# 2. Shrunk (IHW + ashr)
myc_6W_log2FC <- lfcShrink(dds_int, coef = "myc_status_pos_vs_neg", type = "ashr", res = myc_6W_log2FC_raw)
timepoint_neg_log2FC <- lfcShrink(dds_int, coef = "timepoint_12W_vs_6W", type = "ashr", res = timepoint_neg_log2FC_raw)
myc_12W_log2FC <- lfcShrink(dds_int, contrast = list(c("myc_status_pos_vs_neg", "timepoint12W.myc_statuspos")), type = "ashr", res = myc_12W_log2FC_raw)
timepoint_pos_log2FC <- lfcShrink(dds_int, contrast = list(c("timepoint_12W_vs_6W", "timepoint12W.myc_statuspos")), type = "ashr", res = timepoint_pos_log2FC_raw)

# === Build annotated dataframes ===

## Shrunken
combined_df_annotated <- data.frame(
  gene = rownames(myc_6W_log2FC),
  myc_6W_log2FC = myc_6W_log2FC$log2FoldChange,
  timepoint_neg_log2FC = timepoint_neg_log2FC$log2FoldChange,
  myc_12W_log2FC = myc_12W_log2FC$log2FoldChange,
  timepoint_pos_log2FC = timepoint_pos_log2FC$log2FoldChange,
  baseMean = myc_6W_log2FC$baseMean,
  padj = myc_6W_log2FC$padj
)

## Raw
combined_df_annotated_raw <- data.frame(
  gene = rownames(myc_6W_log2FC_raw),
  myc_6W_log2FC_raw = myc_6W_log2FC_raw$log2FoldChange,
  timepoint_neg_log2FC_raw = timepoint_neg_log2FC_raw$log2FoldChange,
  myc_12W_log2FC_raw = myc_12W_log2FC_raw$log2FoldChange,
  timepoint_pos_log2FC_raw = timepoint_pos_log2FC_raw$log2FoldChange,
  baseMean = myc_6W_log2FC_raw$baseMean,
  padj = myc_6W_log2FC_raw$padj
)

# === Add gene symbols ===
gene_annotations <- readRDS("results/ortholog_table.rds") %>%
  dplyr::select(ensembl_gene_id, external_gene_name) %>%
  distinct()

gene_annotations <- gene_annotations %>%
  rename(gene = ensembl_gene_id, mgi_symbol = external_gene_name)

combined_df_annotated <- combined_df_annotated %>%
  left_join(gene_annotations, by = "gene")

combined_df_annotated_raw <- combined_df_annotated_raw %>%
  left_join(gene_annotations, by = "gene")

# === Save results ===
saveRDS(combined_df_annotated, "results/combined_df_annotated.rds")
saveRDS(combined_df_annotated_raw, "results/combined_df_annotated_raw.rds")

