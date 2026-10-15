# =============================================================================
# fig1b_myc_transcript.R -- Fig. 1b, RNA-seq part: Myc mRNA in every animal
# -----------------------------------------------------------------------------
# Sits beside the MYC western in Fig. 1b. The text: "MYC protein levels fell by
# ~50% during the pubertal-adult (6-12 weeks) transition, without any fall in the
# MYC/MAX/MXD network transcripts (Fig. 1b, ED Fig. 1a)." This panel is Myc
# itself: the message does not fall where the protein does. ED Fig. 1a is the
# rest of the network.
#
# FORM: one value per animal, four groups of six -> every animal over a box (median
#   and interquartile range; whiskers minimum to maximum since 2026-10-14, as Prism)
#   (data-to-viz, "Do boxplots hide information?": the points are drawn). A linear
#   axis from 0, so a halving would read as half the height beside the blot.
#   Groups age-major: three of the four tests compare within an age or across the
#   window within a genotype, and age-major keeps every bracket short.
#
# PRISM FORMAT (author, 2026-10-15): the panel sits in the bench row among Prism
#   graphs, so it takes their format rather than the R panels':
#   - the graph (the data area) is 20 x 20 mm, and the y axis ends at the tick just
#     above the top animal (12.5 thousand); labels in thousands, "x1,000" in the title;
#   - no x axis and no group labels: Illustrator sets them, as on the Prism graphs;
#   - axis title 6 pt, tick labels and P values 5 pt;
#   - every line 0.9 pt (lw_pt(0.9); 1 pt read heavier than the Prism graphs): axis,
#     ticks, boxes, whiskers, medians, brackets;
#   - brackets are straight lines, without the downward tips, each just above the
#     animals it compares (second round, to save space): the 12W genotype bracket
#     inside the graph over the 12W groups, the 6W one above the axis end, the MYC+
#     window above that. The WT window (P = 0.48) is not drawn; it goes in the legend.
#   - the title is two lines 2.4 mm apart, close to the axis: set as a column of its
#     own, because plotmath's atop() spaces its lines about twice as far apart.
#   The page is the title, the graph with its axis, and the brackets (natural_size()),
#   so the PDF has no margin to trim.
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
  level  = c(1, 2, 3, 4))       # the genotype pair cannot share a level (header, PRISM FORMAT)
res <- do.call(rbind, lapply(tests$slot, function(s)
  as.data.frame(ir[[s]])[myc, c("log2FoldChange", "lfcSE", "pvalue", "padj")]))
tests <- cbind(tests, res)
tests$p <- tests$padj                                  # IHW-adjusted

# the drawn comparison and the reported statistic are the same comparison
gm  <- tapply(d$y, d$group, mean)
emp <- log2(gm[tests$group2] / gm[tests$group1])
stopifnot(max(abs(emp - tests$log2FoldChange)) < 0.1)

# --- the panel, in the Prism format -------------------------------------------------
GRAPH <- 20           # mm; the data area, as Prism's graph size
LINE  <- lw_pt(0.9)   # every line 0.9 pt (author: 1 pt read heavier than the Prism graphs)
TXT   <- 5            # pt; tick labels and P values
TITLE <- 6            # pt; the axis title

# The y axis ends on the tick just above the top animal: labels every 5 (thousand) and
# minor ticks between, so it ends on 12.5 rather than 15.
STEP <- 5000
TOP  <- ceiling(max(d$y) / (STEP / 2)) * STEP / 2
stopifnot(TOP > max(d$y), (TOP - max(d$y)) / TOP < 0.05)
MM <- TOP / GRAPH                                      # data units per printed mm

# Each bracket sits just above the animals it compares (author, 2026-10-15): the 12W
# genotype bracket inside the graph over the 12W groups, the 6W one above the axis end,
# and the MYC+ window above that. The WT window is not drawn; its P goes in the legend.
top_of <- function(g) max(d$y[d$group %in% g])
drawn  <- tests[tests$slot != "timepoint_neg_raw", ]
drawn$base <- c(myc_6W_raw        = top_of(c("6W_neg", "6W_pos")) + 1.5 * MM,
                myc_12W_raw       = top_of(c("12W_neg", "12W_pos")) + 1.5 * MM,
                timepoint_pos_raw = top_of(c("6W_neg", "6W_pos")) + 4.8 * MM)[drawn$slot]
br <- brackets(drawn[, c("group1", "group2", "p")], base = drawn$base, step = 3.3 * MM,
               txt = TXT, line = LINE, tips = FALSE)
above <- (attr(br, "top") - TOP) / MM                  # mm the brackets need above the axis
stopifnot(drawn$base[drawn$slot == "myc_12W_raw"] + 3 * MM < TOP)   # its P label stays in the graph

p <- ggplot(d, aes(group, y)) +
  geom_boxplot(aes(fill = group), outlier.shape = NA, width = 0.6, coef = Inf,   # min-max whiskers
               colour = box_line, linewidth = LINE, median.linewidth = LINE) +
  ggbeeswarm::geom_quasirandom(aes(colour = group), width = 0.2, size = 0.9, shape = 16) +
  br +
  scale_colour_manual(values = group_cols) +
  scale_fill_manual(values = group_fill) +
  # the brackets lie above the axis end: kept (oob_keep) and drawn there (clip off)
  scale_y_nm(step = STEP, top = max(d$y), minor = TRUE, oob = scales::oob_keep,
             labels = function(x) lab_num(x / 1000)) +
  labs(x = NULL, y = NULL) +
  coord_cartesian(clip = "off") +
  theme_nm(txt = TXT, txt_title = TITLE, line = LINE) +
  # right margin 1.4 mm: the 12W genotype P label overhangs the graph by 0.9 mm
  theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.x = element_blank(), plot.margin = margin(above, 1.4, 0.8, 0, "mm"))
stopifnot(isTRUE(all.equal(ggplot_build(p)$layout$panel_params[[1]]$y$continuous_range, c(0, TOP))))

# The y title: two lines LEAD mm apart, baseline to baseline (6 pt at 1.15), closer than
# plotmath's atop() sets them, so it is drawn as a column of its own left of the axis,
# centred on the graph as an axis title is; x is each line's baseline, GAP mm from the
# tick labels. The break keeps both lines inside the 20 mm.
LEAD <- 2.4; GAP <- 0.6; COL <- 5.2                    # mm
title <- ggplot() +
  annotate("text", x = COL - GAP - c(LEAD, 0), y = 0.5, angle = 90, vjust = 0, parse = TRUE,
           label = c('italic(Myc)~"mRNA (×1,000,"', '"normalised counts)"'),
           size = pt(TITLE), family = NM_FONT, colour = "black") +
  scale_x_continuous(limits = c(0, COL), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme_void(base_family = NM_FONT) + theme(plot.margin = margin(0, 0, 0, 0))

# the page is the title, the graph with its axis, and the brackets, to the half millimetre
panel <- stack_parts(title, p, ncol = 2, widths = grid::unit(c(COL, GRAPH), "mm"),
                     heights = grid::unit(GRAPH, "mm"))
size  <- ceiling(natural_size(panel) * 2) / 2
save_panel(panel, fig = "Fig1", panel = "b", name = "myc_transcript",
           width = size[["width"]], height = size[["height"]])

# --- numbers for the legend (printed, never drawn) --------------------------------
cat(sprintf("\nFig. 1b -- for the legend (graph 20 x 20 mm, saved %g x %g mm)\n", size[["width"]], size[["height"]]))
cat("  n = 6 mice per group; points are animals; boxes show median and interquartile",
    "range, whiskers minimum to maximum.\n")
cat("  DESeq2 Wald test, ~ timepoint * genotype, IHW-adjusted P:\n")
for (i in seq_len(nrow(tests))) {
  cat(sprintf("    %-8s vs %-8s log2FC %+.2f (SE %.2f)  P = %.2g%s\n",
              sub("\n", " ", group_labels[[tests$group2[i]]]), sub("\n", " ", group_labels[[tests$group1[i]]]),
              tests$log2FoldChange[i], tests$lfcSE[i], tests$padj[i],
              if (tests$slot[i] %in% drawn$slot) "" else "   (not drawn: legend only)"))
}

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {
  print(p)
  tests
  tapply(d$y, d$group, summary)
}
