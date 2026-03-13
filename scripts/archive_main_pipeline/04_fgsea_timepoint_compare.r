# scripts/04_fgsea_timepoint_compare.R

source("scripts/00_setup_packages.R")

# Inputs
combined        <- readRDS("results/combined_df_annotated.rds")
genesets_symbol <- readRDS("results/gene_sets_list.rds")

# --- Safety checks -----------------------------------------------------------
stopifnot(all(c("gene","mgi_symbol","mu6_stat","mu12_stat") %in% colnames(combined)))

# If Ensembl IDs may have version suffixes, normalize them once
strip_ver <- function(x) sub("\\.\\d+$","", x)
combined$gene <- strip_ver(combined$gene)

# --- Map gene sets (mouse symbols) -> Ensembl IDs present in 'combined' -----
symbol2ens <- combined %>%
  dplyr::select(gene, mgi_symbol) %>%
  filter(!is.na(mgi_symbol), mgi_symbol != "", !is.na(gene), gene != "") %>%
  mutate(mgi_symbol = toupper(mgi_symbol)) %>%
  distinct(mgi_symbol, gene)

genesets_ens <- lapply(genesets_symbol, function(vec_sym){
  tibble(mgi_symbol = toupper(vec_sym)) %>%
    inner_join(symbol2ens, by = "mgi_symbol") %>%
    pull(gene) %>%
    unique()
})

# --- Rank vectors (Wald stats; shrinkage-agnostic) --------------------------
rank_6W  <- combined$mu6_stat;  names(rank_6W)  <- combined$gene
rank_12W <- combined$mu12_stat; names(rank_12W) <- combined$gene

# Drop NAs and zero-variance edge cases
rank_6W  <- rank_6W[!is.na(rank_6W)]
rank_12W <- rank_12W[!is.na(rank_12W)]

# Optional: tiny jitter if many exact ties (rare, but fgsea warns sometimes)
if (length(unique(rank_6W)) < length(rank_6W)/5)  rank_6W  <- rank_6W  + rnorm(length(rank_6W), 0, 1e-6)
if (length(unique(rank_12W)) < length(rank_12W)/5) rank_12W <- rank_12W + rnorm(length(rank_12W), 0, 1e-6)

# --- Run fgsea ---------------------------------------------------------------
fg_6W  <- fgseaMultilevel(pathways = genesets_ens, stats = rank_6W,  scoreType = "std")
fg_12W <- fgseaMultilevel(pathways = genesets_ens, stats = rank_12W, scoreType = "std")

# --- Compare NES across timepoints ------------------------------------------
gs_comp <- fg_6W %>%
  dplyr::select(pathway, NES, padj) %>%
  dplyr::rename(NES_6W = NES, padj_6W = padj) %>%
  inner_join(
    fg_12W %>% dplyr::select(pathway, NES, padj) %>%
      dplyr::rename(NES_12W = NES, padj_12W = padj),
    by = "pathway"
  ) %>%
  mutate(delta_NES = NES_12W - NES_6W)

# --- Outputs ----------------------------------------------------------------
dir.create("results/fgsea", showWarnings = FALSE)
write_csv(fg_6W,  "results/fgsea/fgsea_6W.csv")
write_csv(fg_12W, "results/fgsea/fgsea_12W.csv")
write_csv(gs_comp,"results/fgsea/fgsea_compare_6W_vs_12W.csv")

# NES scatter
# Add annotation columns
gs_comp <- gs_comp %>%
  mutate(
    category = case_when(
      grepl("^MC_", pathway)  ~ "Mitocarta",
      grepl("^MYC_", pathway) ~ "MYC_signature",
      TRUE ~ "Other"
    ),
    sig_any = (padj_6W < 0.1 | padj_12W < 0.1)
  )

# --- Plot 1: All pathways ---
p_all <- ggplot(gs_comp, aes(NES_6W, NES_12W, color=category)) +
  geom_abline(slope=1, intercept=0, linetype=2) +
  geom_point(aes(size = -log10(pmin(padj_6W, padj_12W))), alpha=0.7) +
  geom_text_repel(aes(label=pathway),
                  size=1.3, max.overlaps=25, box.padding=0.3) +
  labs(x="NES @ 6W (MYC+ vs MYC-)",
       y="NES @ 12W (MYC+ vs MYC-)",
       size="-log10(FDR)",
       color="Pathway category",
       title="All pathways") +
  theme_minimal()

ggsave("results/fgsea/NES_scatter_all.pdf", p_all, width=7, height=6)

# --- Plot 2: Significant pathways only ---
p_sig <- ggplot(filter(gs_comp, sig_any), aes(NES_6W, NES_12W, color=category)) +
  geom_abline(slope=1, intercept=0, linetype=2) +
  geom_point(aes(size = -log10(pmin(padj_6W, padj_12W))), alpha=0.8) +
  geom_text_repel(aes(label=pathway), size=1.3, max.overlaps=25, box.padding=0.3) +
  labs(x="NES @ 6W (MYC+ vs MYC-)",
       y="NES @ 12W (MYC+ vs MYC-)",
       size="-log10(FDR)",
       color="Pathway category",
       title="Significant pathways (FDR < 0.1 at 6W or 12W)") +
  theme_minimal()

ggsave("results/fgsea/NES_scatter_sig.pdf", p_sig, width=7, height=6)


