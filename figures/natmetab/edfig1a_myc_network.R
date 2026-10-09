# =============================================================================
# edfig1a_myc_network.R -- Extended Data Fig. 1a: the MYC/MAX/MXD network
# -----------------------------------------------------------------------------
# The text: "... without any fall in the MYC/MAX/MXD network transcripts (Fig. 1b,
# ED Fig. 1a)." A negative control on the alternative to a dose explanation. MYC
# works as a heterodimer with MAX and competes for the same E-boxes with the
# MXD/MNT repressors, so a fall in MAX or a rise in the repressors would weaken
# MYC's output at constant MYC. The panel shows each gene's MYC effect at both
# ages: if the network moved, the two ages would differ.
#
# FORM: ten genes x two effects, signed -> dots with 95% intervals, 0 marked
#   (data-to-viz, "The Moire effect": dots, not bars; "Order your data": genes
#   ordered by the 6-week effect within each role).
#
# STATISTICS: DESeq2 Wald, raw (unshrunken) ~ timepoint * myc_status fit. The
#   intervals are log2FC +/- 1.96 SE; the adjusted P values are IHW and are
#   printed for the legend, not drawn: nothing here is significant.
#
# Reads:  results/interaction_results.rds
# Output: outputs/natmetab/EDFig1/EDFig1a_myc_network.pdf
# =============================================================================

source(here::here("figures", "natmetab", "_style.R"))
source(here::here("functions", "reconcile_gene_symbols.R"))

ir <- readRDS(here::here("results", "interaction_results.rds"))
universe <- rownames(as.data.frame(ir[["myc_6W_raw"]]))

# the roster, by role. MXI1 is MXD2, so the MXD family is complete; MLX pairs with
# MLXIP (MondoA) and MLXIPL (ChREBP) and competes for MAX-family partners.
net <- data.frame(
  gene = c("Max", "Mxd1", "Mxd3", "Mxd4", "Mxi1", "Mnt", "Mga", "Mlx", "Mlxip", "Mlxipl"),
  role = c("Obligate\npartner", rep("MXD/MNT\nrepressors", 6), rep("MLX arm", 3)),
  stringsAsFactors = FALSE)
net$role <- factor(net$role, levels = unique(net$role))
net$ens  <- vapply(net$gene, function(g) recon_to_ensembl(g, universe), character(1))
stopifnot(!anyNA(net$ens))

grab <- function(slot) {
  r <- as.data.frame(ir[[slot]])[net$ens, ]
  data.frame(lfc = r$log2FoldChange, se = r$lfcSE, padj = r$padj, baseMean = r$baseMean)
}
g6 <- grab("myc_6W_raw"); g12 <- grab("myc_12W_raw"); gi <- grab("interaction_raw")
stopifnot(nrow(g6) == 10L)

ord <- order(net$role, g6$lfc)
net$gene <- factor(net$gene, levels = net$gene[ord])
dl <- rbind(data.frame(net, contrast = "myc_6W",  lfc = g6$lfc,  se = g6$se),
            data.frame(net, contrast = "myc_12W", lfc = g12$lfc, se = g12$se))
dl$contrast <- factor(dl$contrast, levels = names(contrast_cols))
dl$lo <- dl$lfc - 1.96 * dl$se
dl$hi <- dl$lfc + 1.96 * dl$se

dodge <- position_dodge(width = 0.6, orientation = "y", reverse = TRUE)
p <- ggplot(dl, aes(lfc, gene, colour = contrast)) +
  geom_vline(xintercept = 0, linewidth = NM_LINE, colour = box_line) +
  geom_linerange(aes(xmin = lo, xmax = hi), orientation = "y", position = dodge,
                 linewidth = 0.4) +
  geom_point(position = dodge, size = 1.1, shape = 16) +
  facet_grid(role ~ ., scales = "free_y", space = "free_y") +
  scale_colour_manual(values = contrast_cols, labels = contrast_labels, name = NULL) +
  scale_x_nm(step = 0.5, top = max(dl$hi), bottom = min(dl$lo), minor = TRUE,
             labels = lab_signed) +
  labs(x = "MYC+ vs WT (log2 fold change)", y = NULL) +
  theme_nm(legend = "bottom") +
  theme(axis.text.y  = element_text(face = "italic"),
        axis.line.y  = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.y = element_text(size = NM_TXT, angle = 0, hjust = 0),
        panel.spacing.y = unit(1.5, "mm"),
        legend.margin = margin(0, 0, 0, 0))

save_panel(p, fig = "EDFig1", panel = "a", name = "myc_network", width = 80, height = 64)

# --- numbers for the legend (printed, never drawn) --------------------------------
cat("\nED Fig. 1a -- for the legend\n")
cat("  DESeq2 Wald, raw ~ timepoint * genotype fit; bars are log2FC +/- 1.96 SE;",
    "n = 6 mice per group.\n")
cat(sprintf("  smallest IHW-adjusted P, MYC+ vs WT: 6W %s %.3f, 12W %s %.3f\n",
            net$gene[which.min(g6$padj)], min(g6$padj, na.rm = TRUE),
            net$gene[which.min(g12$padj)], min(g12$padj, na.rm = TRUE)))
cat(sprintf("  smallest IHW-adjusted interaction P: %s %.2f\n",
            net$gene[which.min(gi$padj)], min(gi$padj, na.rm = TRUE)))
cat(sprintf("  standard errors at the two ages differ by at most %.3f log2\n",
            max(abs(g6$se - g12$se))))
low <- net$gene[g6$baseMean < 100]
if (length(low)) cat("  low counts (baseMean < 100):", paste(low, collapse = ", "), "\n")

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {
  print(p)
  cbind(gene = as.character(net$gene), round(cbind(g6, lfc12 = g12$lfc, p12 = g12$padj,
                                                  int = gi$lfc, pint = gi$padj), 3))
}
