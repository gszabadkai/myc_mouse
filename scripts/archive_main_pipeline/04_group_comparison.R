# scripts/04_group_comparison.R

source("scripts/00_setup_packages.R")

# === Load input ===
cts <- readRDS("results/count_matrix.rds")
coldata <- readRDS("results/coldata.rds")
combined_df_annotated <- readRDS("results/combined_df_annotated.rds")

# === Build group-based DESeq2 design ===
dds_group <- DESeqDataSetFromMatrix(countData = cts,
                                    colData = coldata,
                                    design = ~ group)

keep <- rowSums(counts(dds_group) >= 10) >= 4
dds_group <- dds_group[keep, ]
dds_group <- DESeq(dds_group)

# Save for reuse
saveRDS(dds_group, "results/dds_group_run.rds")

# === Extract 12W_pos vs 6W_pos comparison ===
res_group <- results(dds_group, contrast = c("group", "12W_pos", "6W_pos"), filterFun = ihw)
res_group_shrunk <- lfcShrink(dds_group,
                              contrast = c("group", "12W_pos", "6W_pos"),
                              type = "ashr",
                              res = res_group)

df_group <- as.data.frame(res_group_shrunk) %>%
  rownames_to_column("gene") %>%
  mutate(
    group_12Wpos_vs_6Wpos_log2FC = log2FoldChange,
    group_12Wpos_vs_6Wpos_padj = res_group[rownames(.), "padj"],
    group_sig_status = case_when(
      group_12Wpos_vs_6Wpos_padj < 0.1 & log2FoldChange > 0 ~ "up",
      group_12Wpos_vs_6Wpos_padj < 0.1 & log2FoldChange < 0 ~ "down",
      TRUE ~ "ns"
    )
  ) %>%
  dplyr::select(gene, group_12Wpos_vs_6Wpos_log2FC, group_sig_status)

df_group <- as.data.frame(res_group_shrunk) %>%
  rownames_to_column("gene") %>%
  mutate(
    group_12Wpos_vs_6Wpos_log2FC = log2FoldChange,
    group_12Wpos_vs_6Wpos_padj = res_group[rownames(res_group_shrunk), "padj"],
    group_sig_status = case_when(
      group_12Wpos_vs_6Wpos_padj < 0.1 & log2FoldChange > 0 ~ "up",
      group_12Wpos_vs_6Wpos_padj < 0.1 & log2FoldChange < 0 ~ "down",
      TRUE ~ "ns"
    )
  ) %>%
  dplyr::select(gene, group_12Wpos_vs_6Wpos_log2FC, group_sig_status)


# === Join into combined_df_annotated ===
combined_df_annotated <- combined_df_annotated %>%
  left_join(df_group, by = "gene")

combined_df_annotated_raw <- combined_df_annotated_raw %>%
  left_join(df_group, by = "gene")

# Save both
saveRDS(combined_df_annotated,     "results/combined_df_annotated.rds")
saveRDS(combined_df_annotated_raw, "results/combined_df_annotated_raw.rds")

