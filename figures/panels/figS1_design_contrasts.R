# =============================================================================
# figS1_design_contrasts.R -- the four comparisons, drawn once
# -----------------------------------------------------------------------------
# SLOT: Fig. S1A (was S1B until 2026-07-31; paragraph 1 as rewritten cites the
# design before the library, so S1A and S1B swapped). Filenames no longer carry
# the slot letter -- figures/panels/PANELS.md is the slug -> slot map.
#
# SUPPORTS (Results, "Myc drives early breast tumourigenesis by inducing OXPHOS
# and biosynthetic pathways", paragraph 1):
#   "we compared the transcriptome of the purified MEC population using both
#    longitudinal (6W versus 12W for WT or Myc+ genotypes) and cross-sectional
#    (WT versus Myc+ at 6W or 12W) comparisons (Fig. S1A)"
#
# A schematic, not a result: the 2x2 design and the five contrasts the DESeq2
# `~ timepoint * myc_status` model yields, drawn once so every later panel can
# be read against it. Contrast labels are the author's figure vocabulary
# (_panel_common.R: contrast_levels) and the legend maps them to the
# results/interaction_results.rds slot names.
#
# WHY THE TWO AXES ARE DRAWN DIFFERENTLY. The key says "genotype" and
# "development" -- the two contrast families -- and nothing more. The 6W and 12W
# cohorts were extracted as two separate batches, so batch is perfectly
# confounded with age and the development contrasts are not separable from it;
# that is a Methods statement (author, 2026-07-30), so it is in the legend block
# below and NOT on the page.
#
# Input:  none (geometry only; labels from figures/theme_myc.R)
# Output: outputs/figures/panels/figS1_design_contrasts.pdf
# =============================================================================

source(here::here("figures", "panels", "_panel_common.R"))

# --- node grid: x = age, y = genotype (WT on top, control first) -------------
# x spans ~1.29 units across ~76 mm and y spans ~1.58 across ~38 mm, so an x unit
# is about 2.4x a y unit on the page: HALF_W is set well below HALF_H. The node
# ends up a landscape rectangle roughly 13 x 10 mm, deliberately the same shape
# as a tile of Fig. 1B, so the two panels read as the same object. (Tightening
# the x limits to pull the WT / Myc+ labels in shrank the x range, so HALF_W came
# down from 0.13 with it -- the node's physical size on the page is unchanged.)
HALF_W <- 0.11          # node half-width  (x units)
HALF_H <- 0.19          # node half-height (y units)
GAP    <- 0.03          # clearance between a node edge and an arrow tail

# `ink` follows the palette's lightness, not the genotype: the 6W fills are the
# saturated Okabe-Ito pair and the 12W fills are the light pair, so "n = 6" is
# white at 6W and near-black at 12W.
nodes <- data.frame(
  group = c("6W_neg", "12W_neg", "6W_pos", "12W_pos"),
  x     = c(1,        2,         1,        2),
  y     = c(2,        2,         1,        1),
  ink   = c("white",  "grey10",  "white",  "grey10"),
  stringsAsFactors = FALSE)
nodes$fill <- unname(group_cols[nodes$group])
stopifnot(!any(is.na(nodes$fill)))

# --- contrasts ---------------------------------------------------------------
# `kind` drives the linetype and is the contrast FAMILY: horizontal = development
# (6->12W within one genotype), vertical = genotype (the Myc effect at one age).
# Labels are the author's figure vocabulary; the legend maps them to the DESeq2
# slot names (6>12W_wt = timepoint_neg, 6>12W_myc = timepoint_pos).
arrows <- data.frame(
  label = c(contrast_dev, contrast_geno),
  kind  = c("development", "development", "genotype", "genotype"),
  x     = c(1 + HALF_W + GAP, 1 + HALF_W + GAP, 1, 2),
  xend  = c(2 - HALF_W - GAP, 2 - HALF_W - GAP, 1, 2),
  y     = c(2, 1, 2 - HALF_H - GAP, 2 - HALF_H - GAP),
  yend  = c(2, 1, 1 + HALF_H + GAP, 1 + HALF_H + GAP),
  stringsAsFactors = FALSE)

arrow_lab <- data.frame(
  label = c(contrast_dev, contrast_geno),
  x     = c(1.5,  1.5,  1 - 0.105, 2 + 0.105),
  y     = c(2.14, 0.86, 1.5,       1.5),
  angle = c(0,    0,    90,        90),
  stringsAsFactors = FALSE)

# grid draws a closed arrowhead with the segment's own line type, so an arrow on
# a dashed segment gets a dashed HEAD -- it renders as a broken triangle. Split
# each arrow into a dashed shaft and a short solid stub that carries the head.
STUB <- 0.13                                   # fraction of the segment
arrows$hx <- arrows$xend - STUB * (arrows$xend - arrows$x)
arrows$hy <- arrows$yend - STUB * (arrows$yend - arrows$y)

# The interaction is the difference between the two genotype contrasts, so it is
# drawn as the connector between them rather than as a fifth arrow. Same stub
# treatment, at both ends.
INT_X0 <- 1 + 0.045; INT_X1 <- 2 - 0.045; INT_STUB <- 0.05
int_seg  <- data.frame(x = INT_X0 + INT_STUB, xend = INT_X1 - INT_STUB,
                       y = 1.5, yend = 1.5)
int_head <- data.frame(x    = c(INT_X0 + INT_STUB, INT_X1 - INT_STUB),
                       xend = c(INT_X0,            INT_X1),
                       y    = c(1.5, 1.5), yend = c(1.5, 1.5))

line_types <- c(genotype = "solid", development = "22")
line_lab   <- c(genotype = "genotype", development = "development")

# --- panel -------------------------------------------------------------------
p <- ggplot2::ggplot() +
  # the four groups
  ggplot2::geom_tile(data = nodes,
                     ggplot2::aes(x = x, y = y),
                     fill = nodes$fill, colour = NA,
                     width = 2 * HALF_W, height = 2 * HALF_H) +
  ggplot2::geom_text(data = nodes,
                     ggplot2::aes(x = x, y = y, label = "n = 6"),
                     colour = nodes$ink, size = 2.2) +
  # the interaction connector: dotted shaft, solid stubs carrying the heads
  ggplot2::geom_segment(data = int_seg,
                        ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
                        linetype = "12", linewidth = 0.3, colour = "grey45") +
  ggplot2::geom_segment(data = int_head,
                        ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
                        linewidth = 0.3, colour = "grey45",
                        arrow = grid::arrow(length = grid::unit(1.3, "mm"),
                                            type = "closed")) +
  ggplot2::annotate("text", x = 1.5, y = 1.60, label = "interaction",
                    size = 2.1, colour = "grey30") +
  # the four contrasts: shaft carries the linetype, stub carries the head
  ggplot2::geom_segment(data = arrows,
                        ggplot2::aes(x = x, xend = hx, y = y, yend = hy,
                                     linetype = kind),
                        linewidth = 0.4, colour = "grey20") +
  ggplot2::geom_segment(data = arrows,
                        ggplot2::aes(x = hx, xend = xend, y = hy, yend = yend),
                        linewidth = 0.4, colour = "grey20",
                        arrow = grid::arrow(length = grid::unit(1.6, "mm"),
                                            type = "closed")) +
  ggplot2::geom_text(data = arrow_lab,
                     ggplot2::aes(x = x, y = y, label = label, angle = angle),
                     size = 2.1, colour = "grey20") +
  # axes carry the factor levels, so the nodes stay uncluttered. The x limits are
  # tight against the genotype arrow labels at x = 1 +/- 0.105, which pulls the
  # WT / Myc+ labels in close to the scheme (author, 2026-07-30); the residual
  # gap is the axis-text margin, set to a hairline below.
  ggplot2::scale_x_continuous(breaks = c(1, 2), labels = c("6 weeks", "12 weeks"),
                              limits = c(0.855, 2.145)) +
  ggplot2::scale_y_continuous(breaks = c(1, 2),
                              labels = c(unname(geno_labels[["pos"]]),
                                         unname(geno_labels[["neg"]])),
                              limits = c(0.72, 2.30)) +
  ggplot2::scale_linetype_manual(values = line_types, labels = line_lab,
                                 name = NULL,
                                 breaks = c("genotype", "development")) +
  ggplot2::labs(x = NULL, y = NULL) +
  theme_panel() +
  ggplot2::theme(
    axis.line       = ggplot2::element_blank(),
    axis.ticks      = ggplot2::element_blank(),
    axis.text       = ggplot2::element_text(colour = "black", face = "bold"),
    axis.text.y     = ggplot2::element_text(colour = "black", face = "bold",
                                            margin = ggplot2::margin(r = 0.3,
                                                                     unit = "mm")),
    plot.margin     = ggplot2::margin(2, 2, 2, 0.5, "mm"),
    legend.position = "bottom",
    legend.margin   = ggplot2::margin(-2, 0, 0, 0),
    legend.key.width = ggplot2::unit(6, "mm"))

# --- the legend text (never drawn) -------------------------------------------
LEGEND <- panel_legend(
  slot = "Fig. S1A",
  what = paste0(
    "Design and the four comparisons. Purified mammary epithelial cells from ",
    "MMTV-Myc transgenic mice and wild-type littermates at 6 and 12 weeks, ",
    "six animals per group (n = 24)."),
  detail = c(
    "Model: DESeq2 on a `~ timepoint * myc_status` interaction design.",
    "Genotype contrasts (solid, vertical): myc_6W and myc_12W, the Myc effect at each age.",
    "Development contrasts (dashed, horizontal): 6>12W_wt and 6>12W_myc, the 6 to 12 week change within each genotype.",
    "Interaction (dotted): the difference between the two genotype contrasts (myc_12W minus myc_6W), i.e. how much of the Myc effect is lost by 12 weeks; equivalently the difference between the two development contrasts.",
    "The four drawn names map to the slots of results/interaction_results.rds as myc_6W, myc_12W, timepoint_neg (= 6>12W_wt) and timepoint_pos (= 6>12W_myc)."),
  bounds = c(
    "Batch = timepoint (Methods). The 6W and 12W cohorts were extracted as two separate batches, so batch is perfectly confounded with age and the two development contrasts are not separable from it. Genotype is balanced within each batch, so the design absorbs the batch effect into the timepoint main effect.",
    "Consequently every genotype gap and every interaction is clean, and so is any difference between two temporal contrasts, because the shared offset cancels in a difference. The pure between-age comparison - which is the developmental reading - cannot be separated from batch post hoc; no dissociation-batch, viability or RIN metadata exists.",
    "Purification is by enzymatic dissociation only, with no sorting step, so residual stromal, endothelial and immune signal (3-12%) is contamination rather than tissue composition, and the warm digest itself induces an immediate-early signature."),
  source = c(
    "Design and contrasts: scripts/03_deseq_results_qc.R; contrast objects in results/interaction_results.rds",
    "Batch = timepoint: CLAUDE.md, and paper/myc_mito.qmd @sec-overview"))

save_panel_p(p, "figS1_design_contrasts",
             width = fig_w[["single"]], height = 52)

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {

  print(p)
  print(LEGEND)

  ## the geometry, if a node or arrow needs nudging
  nodes  |> print()
  arrows |> print()
}
