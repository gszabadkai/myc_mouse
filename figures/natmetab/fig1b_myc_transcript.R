# =============================================================================
# fig1b_myc_transcript.R -- Fig. 1b, RNA-seq part: Myc mRNA in every animal
# -----------------------------------------------------------------------------
# Sits beside the MYC western in Fig. 1b. The text: "MYC protein levels fell by
# ~50% during the pubertal-adult (6-12 weeks) transition, without any fall in the
# MYC/MAX/MXD network transcripts (Fig. 1b, ED Fig. 1a)." This panel is Myc
# itself: the message does not fall where the protein does. ED Fig. 1a is the
# rest of the network.
#
# FORM: one value per animal, four groups of six -> every animal over a box
#   (data-to-viz, "Do boxplots hide information?": the points are drawn). A linear
#   axis from 0, so a halving would read as half the height beside the blot.
#   Groups age-major: three of the four tests compare within an age or across the
#   window within a genotype, and age-major keeps every bracket short.
#
# STATISTICS: DESeq2 Wald tests on the raw (unshrunken) ~ timepoint * myc_status
#   fit, IHW-adjusted -- read from results/interaction_results.rds, never computed
#   here. The script asserts that the drawn group means give the same log2 ratios.
#
# Reads:  results/dds_int_run.rds          (normalized counts, sample table)
#         results/interaction_results.rds  (the four contrasts)
# Output: outputs/natmetab/Fig1/Fig1b_myc_transcript.pdf
# =============================================================================

source(here::here("figures", "natmetab", "_style.R"))
source(here::here("functions", "reconcile_gene_symbols.R"))
suppressPackageStartupMessages(library(DESeq2))

dds <- readRDS(here::here("results", "dds_int_run.rds"))
ir  <- readRDS(here::here("results", "interaction_results.rds"))

myc <- recon_to_ensembl("Myc", rownames(dds))
stopifnot(identical(unname(myc), "ENSMUSG00000022346"))

cd <- as.data.frame(colData(dds))
d  <- data.frame(group = factor(paste(cd$timepoint, cd$myc_status, sep = "_"),
                                levels = order_by_age),
                 y = counts(dds, normalized = TRUE)[myc, ])
stopifnot(!anyNA(d$group), all(table(d$group) == 6L))

# --- the four tests, as reported -----------------------------------------------
tests <- data.frame(
  group1 = c("6W_neg", "12W_neg", "6W_neg", "6W_pos"),
  group2 = c("6W_pos", "12W_pos", "12W_neg", "12W_pos"),
  slot   = c("myc_6W_raw", "myc_12W_raw", "timepoint_neg_raw", "timepoint_pos_raw"),
  level  = c(1, 1, 2, 3))
res <- do.call(rbind, lapply(tests$slot, function(s)
  as.data.frame(ir[[s]])[myc, c("log2FoldChange", "lfcSE", "pvalue", "padj")]))
tests <- cbind(tests, res)
tests$p <- tests$padj                                  # IHW-adjusted

# the drawn comparison and the reported statistic are the same comparison
gm  <- tapply(d$y, d$group, mean)
emp <- log2(gm[tests$group2] / gm[tests$group1])
stopifnot(max(abs(emp - tests$log2FoldChange)) < 0.1)

# --- the panel -------------------------------------------------------------------
br <- brackets(tests[, c("group1", "group2", "p", "level")],
               base = max(d$y) * 1.06, step = 1700)

p <- ggplot(d, aes(group, y)) +
  geom_boxplot(aes(fill = group), outlier.shape = NA, width = 0.6,
               colour = box_line, linewidth = NM_LINE) +
  ggbeeswarm::geom_quasirandom(aes(colour = group), width = 0.2, size = 0.9, shape = 16) +
  br +
  scale_colour_manual(values = group_cols) +
  scale_fill_manual(values = group_fill) +
  scale_x_discrete(labels = group_labels) +
  scale_y_nm(step = 4000, top = attr(br, "top"), minor = TRUE) +
  labs(x = NULL, y = expression(italic(Myc) ~ "mRNA (normalized counts)")) +
  theme_nm()

save_panel(p, fig = "Fig1", panel = "b", name = "myc_transcript", width = 48, height = 56)

# --- numbers for the legend (printed, never drawn) --------------------------------
cat("\nFig. 1b -- for the legend\n")
cat("  n = 6 mice per group; points are animals; boxes show median and interquartile",
    "range, whiskers 1.5 x IQR.\n")
cat("  DESeq2 Wald test, ~ timepoint * genotype, IHW-adjusted P:\n")
for (i in seq_len(nrow(tests))) {
  cat(sprintf("    %-8s vs %-8s log2FC %+.2f (SE %.2f)  P = %.2g\n",
              sub("\n", " ", group_labels[[tests$group2[i]]]), sub("\n", " ", group_labels[[tests$group1[i]]]),
              tests$log2FoldChange[i], tests$lfcSE[i], tests$padj[i]))
}

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {
  print(p)
  tests
  tapply(d$y, d$group, summary)
}
