# =============================================================================
# PROPOSAL (alternative to p1h) -- Fig. 1h as 2G (alt)'s plane: the BCL-2 family
# and Foxo3 on both timelines
# -----------------------------------------------------------------------------
# x = 12W vs 6W in WT, y = 12W vs 6W in MYC+ (DESeq2 log2 fold changes). Dashed:
# the same change in both genotypes. Dotted: 0.20 either side, the declared
# MYC-specific threshold (script 54's rule). Distance below the dashed line is the
# interaction, which is within batch; the axes themselves each coincide with batch.
# Black: the four genes the text names; grey: the rest of the BCL-2 family.
# Shape: pro-apoptotic (circle), anti-apoptotic (square), Foxo3 (diamond).
#
# Reads: results/interaction_results.rds (script 03);
#        results/priming_arm_teb.rds (script 42, $machinery: cross-check)
# =============================================================================

source(here::here("figures", "natmetab", "_style.R"))
source(here::here("functions", "reconcile_gene_symbols.R"))

ir <- readRDS(here::here("results", "interaction_results.rds"))
U  <- rownames(as.data.frame(ir$interaction_raw))
PRO  <- c("Bax", "Bak1", "Bid", "Bbc3", "Bcl2l11", "Bmf", "Pmaip1")
ANTI <- c("Bcl2l1", "Mcl1", "Bcl2")
G    <- c(PRO, ANTI, "Foxo3")
NAMED <- c("Bbc3", "Bax", "Bcl2l1", "Foxo3")
ens <- vapply(G, function(g) recon_to_ensembl(g, U), character(1))
v   <- function(slot, col) as.data.frame(ir[[slot]])[ens, col]
d <- data.frame(gene = G, x = v("timepoint_neg_raw", "log2FoldChange"),
                y = v("timepoint_pos_raw", "log2FoldChange"),
                int = v("interaction_raw", "log2FoldChange"), int_p = v("interaction_raw", "pvalue"))
d$class <- factor(ifelse(d$gene %in% PRO, "pro", ifelse(d$gene %in% ANTI, "anti", "foxo")),
                  levels = c("pro", "anti", "foxo"))
d$named <- d$gene %in% NAMED

m <- as.data.frame(readRDS(here::here("results", "priming_arm_teb.rds"))$machinery)
k <- match(c(PRO, ANTI), m$gene)
stopifnot(!anyNA(k), max(abs(m$lfc_wt_time[k] - d$x[1:10])) < 1e-9,
          max(abs(m$lfc_myc_time[k] - d$y[1:10])) < 1e-9,
          max(abs((d$y - d$x) - d$int)) < 1e-9)

lo <- -0.6; hi <- 1.1
xs <- c(lo, hi)
p <- ggplot(d, aes(x, y)) +
  geom_hline(yintercept = 0, linewidth = NM_LINE, colour = "grey80") +
  geom_vline(xintercept = 0, linewidth = NM_LINE, colour = "grey80") +
  geom_abline(slope = 1, intercept = 0, linetype = "22", linewidth = NM_LINE, colour = "grey30") +
  geom_abline(slope = 1, intercept = c(-0.2, 0.2), linetype = "12", linewidth = NM_LINE, colour = "grey55") +
  geom_point(aes(shape = class, colour = named), size = 1.4) +
  ggrepel::geom_text_repel(aes(label = gene), size = pt(NM_TXT), fontface = "italic",
                           family = NM_FONT, colour = "black", seed = 2,
                           min.segment.length = 0.2, segment.size = 0.2, segment.colour = "grey50",
                           box.padding = 0.15, point.padding = 0.1, max.overlaps = Inf) +
  scale_shape_manual(values = c(pro = 16, anti = 15, foxo = 18), name = NULL,
                     labels = c(pro = "Pro\u00adapoptotic", anti = "Anti\u00adapoptotic", foxo = "PUMA activator")) +
  scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = "grey60"), guide = "none") +
  scale_x_nm(step = 0.5, top = hi, bottom = lo, minor = TRUE, labels = lab_signed) +
  scale_y_nm(step = 0.5, top = hi, bottom = lo, minor = TRUE, labels = lab_signed) +
  coord_equal() +
  labs(x = "12W vs 6W, WT (log2 fold change)", y = "12W vs 6W, MYC+ (log2 fold change)") +
  guides(shape = guide_legend(override.aes = list(colour = "black", size = 1.3))) +
  theme_nm(legend = "inside") +
  theme(legend.position.inside = c(0.02, 0.98), legend.justification.inside = c(0, 1),
        legend.key.size = unit(2.4, "mm"), legend.background = element_blank(),
        legend.key = element_blank())

save_panel(p, fig = "Fig1", panel = "h", name = "alt_bcl2_plane", width = 60, height = 60,
           dir = here::here("outputs", "natmetab", "proposals"))
print(d[order(d$int), c("gene", "x", "y", "int", "int_p")], digits = 3, row.names = FALSE)
