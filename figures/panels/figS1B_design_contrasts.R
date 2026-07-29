# =============================================================================
# figS1B_design_contrasts.R -- the four comparisons, drawn once
# -----------------------------------------------------------------------------
# SUPPORTS (Results, "Myc drives early breast tumourigenesis by inducing OXPHOS
# and biosynthetic pathways", paragraph 1):
#   "to compare the transcriptome on the purified MEC population, using
#    timeline: 6W vs 12W of the WT or Myc+ genotypes, or cross sectional:
#    WT vs Myc+ at 6W or 12W comparisons"
#
# A schematic, not a result: the 2x2 design and the five contrasts the DESeq2
# `~ timepoint * myc_status` model yields, drawn once so every later panel can
# be read against it. Contrast labels are the names used throughout the analysis
# (results/interaction_results.rds slots), so a reader can go from the panel to
# the object without a translation step.
#
# WHY THE TWO AXES ARE DRAWN DIFFERENTLY. The 6W and 12W cohorts were extracted
# as two separate batches, so batch is perfectly confounded with age; genotype
# is balanced WITHIN each batch. The vertical (genotype) contrasts are therefore
# clean and the horizontal (age) contrasts are not separable from batch. That is
# a property of the design, which is what this panel is for, so it is encoded --
# solid vs dashed with a two-entry key -- rather than written on the page.
#
# Input:  none (geometry only; labels from figures/theme_myc.R)
# Output: outputs/figures/panels/figS1B_design_contrasts.pdf
# =============================================================================

source(here::here("figures", "panels", "_panel_common.R"))

# --- node grid: x = age, y = genotype (WT on top, control first) -------------
# x spans ~1.8 units across ~70 mm and y spans ~1.6 across ~34 mm, so an x unit
# is about 1.8x an y unit on the page: HALF_W is set smaller than HALF_H to make
# the nodes read square.
HALF_W <- 0.13          # node half-width  (x units)
HALF_H <- 0.19          # node half-height (y units)
GAP    <- 0.03          # clearance between a node edge and an arrow tail

nodes <- data.frame(
  group = c("6W_neg", "12W_neg", "6W_pos", "12W_pos"),
  x     = c(1,        2,         1,        2),
  y     = c(2,        2,         1,        1),
  ink   = c("black",  "white",   "black",  "white"),   # legibility on the fill
  stringsAsFactors = FALSE)
nodes$fill <- unname(group_cols[nodes$group])
stopifnot(!any(is.na(nodes$fill)))

# --- contrasts ---------------------------------------------------------------
# `kind` drives the linetype: "genotype" = within batch (clean),
# "age" = between batch (confounded). Names match interaction_results.rds.
arrows <- data.frame(
  label = c("timepoint_neg", "timepoint_pos", "myc_6W", "myc_12W"),
  kind  = c("age", "age", "genotype", "genotype"),
  x     = c(1 + HALF_W + GAP, 1 + HALF_W + GAP, 1, 2),
  xend  = c(2 - HALF_W - GAP, 2 - HALF_W - GAP, 1, 2),
  y     = c(2, 1, 2 - HALF_H - GAP, 2 - HALF_H - GAP),
  yend  = c(2, 1, 1 + HALF_H + GAP, 1 + HALF_H + GAP),
  stringsAsFactors = FALSE)

arrow_lab <- data.frame(
  label = c("timepoint_neg", "timepoint_pos", "myc_6W", "myc_12W"),
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

line_types <- c(genotype = "solid", age = "22")
line_lab   <- c(genotype = "genotype (within batch)", age = "age (= batch)")

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
  # axes carry the factor levels, so the nodes stay uncluttered
  ggplot2::scale_x_continuous(breaks = c(1, 2), labels = c("6 weeks", "12 weeks"),
                              limits = c(0.74, 2.26)) +
  ggplot2::scale_y_continuous(breaks = c(1, 2),
                              labels = c(unname(geno_labels[["pos"]]),
                                         unname(geno_labels[["neg"]])),
                              limits = c(0.72, 2.30)) +
  ggplot2::scale_linetype_manual(values = line_types, labels = line_lab,
                                 name = NULL, breaks = c("genotype", "age")) +
  ggplot2::labs(x = NULL, y = NULL) +
  theme_panel() +
  ggplot2::theme(
    axis.line       = ggplot2::element_blank(),
    axis.ticks      = ggplot2::element_blank(),
    axis.text       = ggplot2::element_text(colour = "black", face = "bold"),
    legend.position = "bottom",
    legend.margin   = ggplot2::margin(-2, 0, 0, 0),
    legend.key.width = ggplot2::unit(6, "mm"))

# --- the legend text (never drawn) -------------------------------------------
LEGEND <- panel_legend(
  slot = "Fig. S1B",
  what = paste0(
    "Design and the four comparisons. Purified mammary epithelial cells from ",
    "MMTV-Myc transgenic mice and wild-type littermates at 6 and 12 weeks, ",
    "six animals per group (n = 24)."),
  detail = c(
    "Model: DESeq2 on a `~ timepoint * myc_status` interaction design.",
    "Cross-sectional (vertical): myc_6W and myc_12W, the genotype gap at each age.",
    "Timeline (horizontal): timepoint_neg and timepoint_pos, the 6W to 12W change within each genotype.",
    "Interaction: the difference between the two genotype gaps (myc_12W minus myc_6W), i.e. how much of the Myc effect is lost by 12 weeks.",
    "Contrast names are the slot names in results/interaction_results.rds, so figure and analysis object use one vocabulary."),
  bounds = c(
    "Batch = timepoint. The 6W and 12W cohorts were extracted as two separate batches, so batch is perfectly confounded with age (dashed arrows). Genotype is balanced within each batch, so the design absorbs the batch effect into the timepoint main effect.",
    "Consequently every genotype gap and every interaction is clean, and so is any difference between two temporal contrasts, because the shared offset cancels in a difference. The pure between-age comparison - which is the developmental reading - cannot be separated from batch post hoc; no dissociation-batch, viability or RIN metadata exists.",
    "Purification is by enzymatic dissociation only, with no sorting step, so residual stromal, endothelial and immune signal (3-12%) is contamination rather than tissue composition, and the warm digest itself induces an immediate-early signature."),
  source = c(
    "Design and contrasts: scripts/03_deseq_results_qc.R; contrast objects in results/interaction_results.rds",
    "Batch = timepoint: CLAUDE.md, and paper/myc_mito.qmd @sec-overview"))

save_panel_p(p, "figS1B_design_contrasts",
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
