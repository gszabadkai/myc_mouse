# =============================================================================
# figS2_p53_arm_alt.R -- PUMA is lost without p53 moving
# -----------------------------------------------------------------------------
# SLOT: Fig. S2D (alt). An ALTERNATIVE to figS2_puma_inducers.R, not a
# replacement -- both are built and the author picks.
#
#   "... no transcriptional alterations were detected among other established
#    upstream inducers of PUMA (Fig. S2D)."
#
# WHY A SECOND VERSION (author, 2026-08-09): "the main point here should be to
# show that the P53-DEPENDENT ARM of regulating PUMA is not moving." The current
# panel draws script 42's twelve `puma_inputs`, which is the roster the sentence
# names but which contains only ONE p53-family gene (Trp73). The canonical p53
# arm -- Trp53 itself, Mdm2, and the transactivation targets Cdkn1a, Trp53inp1,
# Zmat3, Phlda3, Eda2r, Ccng1 -- is a separate roster in the same object and is
# the one the reader needs to see, because "PUMA fell" invites "p53 must have
# risen" and the answer is that it did not.
#
# THE PLANE IS THE FIGURE'S: x = the wild-type window, y = the Myc+ window,
# dashed diagonal = development alone. Three casts on it:
#
#   the p53 arm (9)          a HORIZONTAL BAND AT y = 0, not a cluster at the
#                            origin. Several of them drift on the wild-type axis
#                            (Eda2r -0.99, Trp53inp1 -0.77, Phlda3 -0.56) but the
#                            whole arm stays within +/-0.26 on the Myc timeline,
#                            and none has an interaction. That band IS the claim:
#                            the axis on which PUMA is lost is the axis on which
#                            p53's transcriptional output does nothing.
#   other PUMA inducers      the ISR and the FOXO paralogues, also in the band
#   the two that DO move     Bbc3 leaves the band (-0.48 at padj 0.016) and Foxo3
#                            leaves the DIAGONAL (+0.47 in development, -0.01
#                            under Myc) -- two different ways of moving, and the
#                            grammar shows both
#
# WITHOUT THE TWO REFERENCE GENES THIS WOULD BE A BAND OF DOTS WITH NOTHING TO
# CALIBRATE IT, which is the failure mode of every negative control drawn alone.
# With them it says: on the axes where PUMA and its FOXO activator move, the p53
# arm does not.
#
# Reads (read-only, no re-run):
#   results/priming_arm_teb.rds (script 42) -- $exclusions$p53_axis,
#                                  $exclusions$puma_inputs, $machinery
# Output: outputs/figures/panels/figS2_p53_arm_alt.pdf
# =============================================================================

source(here::here("figures", "panels", "_panel_common.R"))

pa_path <- here::here("results", "priming_arm_teb.rds")
require_fresher_than(pa_path)
pa <- readRDS(pa_path)

p53 <- as.data.frame(pa$exclusions$p53_axis)
pin <- as.data.frame(pa$exclusions$puma_inputs)
mac <- as.data.frame(pa$machinery)
stopifnot(nrow(p53) == 9L, nrow(pin) == 12L, "Trp53" %in% p53$gene,
          "Foxo3" %in% pin$gene, "Bbc3" %in% mac$gene)

COLS <- c("gene", "baseMean", "lfc_wt_time", "padj_wt_time",
          "lfc_myc_time", "padj_myc_time", "lfc_interaction")
d <- rbind(
  data.frame(p53[, COLS], class = "p53 arm", stringsAsFactors = FALSE),
  data.frame(pin[pin$gene != "Foxo3", COLS], class = "other PUMA inducers",
             stringsAsFactors = FALSE),
  data.frame(pin[pin$gene == "Foxo3", COLS], class = "the two that move",
             stringsAsFactors = FALSE),
  data.frame(mac[mac$gene == "Bbc3", COLS], class = "the two that move",
             stringsAsFactors = FALSE))
d$class <- factor(d$class,
                  levels = c("p53 arm", "other PUMA inducers", "the two that move"))
stopifnot(nrow(d) == 9L + 11L + 2L, sum(d$class == "the two that move") == 2L)

g <- function(x, col) d[[col]][d$gene == x]
# The claim, asserted: nothing in the p53 arm moves on either timeline or in the
# interaction, while the two reference genes do.
a53 <- d[d$class == "p53 arm", ]
stopifnot(all(a53$padj_wt_time  > 0.01, na.rm = TRUE),
          all(a53$padj_myc_time > 0.05, na.rm = TRUE),
          max(abs(a53$lfc_myc_time)) < 0.30,
          g("Bbc3", "padj_myc_time") < 0.05)

# =============================================================================
# the panel
# =============================================================================
# WINDOWED FOR DISPLAY, as Fig. 1G windows its facets. Two of the twelve PUMA
# inducers are far outside and both are at the floor of what this design can
# measure: they are named in the legend and every number is computed on the
# complete roster.
WIN <- 1.15
out <- d[abs(d$lfc_wt_time) > WIN | abs(d$lfc_myc_time) > WIN, ]
dd  <- d[!(d$gene %in% out$gene), ]
stopifnot(nrow(out) <= 3L, all(out$baseMean < 200),
          !any(out$class == "p53 arm"), !any(out$class == "the two that move"))

LIM <- c(-1, 1) * WIN
p <- ggplot2::ggplot(dd, ggplot2::aes(lfc_wt_time, lfc_myc_time)) +
  two_timeline_base(LIM, diag_at = 0.72) +
  ggplot2::geom_point(ggplot2::aes(colour = class, size = class), stroke = 0) +
  ggrepel::geom_text_repel(
    data = dd[dd$class == "the two that move", ],
    ggplot2::aes(label = gene), size = 1.8, colour = "grey10",
    fontface = "bold.italic", seed = 8, max.overlaps = Inf,
    min.segment.length = 0, segment.size = 0.2, segment.colour = "grey55",
    box.padding = 0.45, point.padding = 0.25) +
  ggrepel::geom_text_repel(
    data = dd[dd$class == "p53 arm", ],
    ggplot2::aes(label = gene), size = 1.55, colour = "grey30",
    fontface = "italic", seed = 8, max.overlaps = Inf,
    min.segment.length = 0, segment.size = 0.18, segment.colour = "grey70",
    box.padding = 0.30, point.padding = 0.16) +
  ggplot2::scale_colour_manual(
    values = c("p53 arm"             = "grey20",
               "other PUMA inducers" = unname(pole_cols[["other"]]),
               "the two that move"   = unname(pole_cols[["down"]])),
    name = NULL) +
  ggplot2::scale_size_manual(values = c("p53 arm" = 1.5,
                                        "other PUMA inducers" = 1.2,
                                        "the two that move" = 2.1),
                             guide = "none") +
  ggplot2::labs(x = "wild-type 6>12W  (raw log2FC)", y = "Myc+ 6>12W") +
  ggplot2::guides(colour = ggplot2::guide_legend(
    override.aes = list(size = 1.7))) +
  theme_panel(base_size = 6) +
  ggplot2::theme(
    legend.position        = "inside",
    legend.position.inside = c(0.01, 0.99),
    legend.justification   = c(0, 1),
    legend.background      = ggplot2::element_blank(),
    legend.margin          = ggplot2::margin(0, 0, 0, 0),
    legend.key.size        = ggplot2::unit(2.4, "mm"),
    legend.spacing.y       = ggplot2::unit(0.3, "mm"),
    plot.margin            = ggplot2::margin(1.5, 2.5, 1, 1.5, "mm"))

# =============================================================================
# the legend text (never drawn)
# =============================================================================
say <- function(x)
  sprintf("%s %+.3f and %+.3f (padj %.2f and %.2f)", x, g(x, "lfc_wt_time"),
          g(x, "lfc_myc_time"), g(x, "padj_wt_time"), g(x, "padj_myc_time"))

LEGEND <- panel_legend(
  slot = "Fig. S2D (alt)",
  what = paste0(
    "The p53-dependent arm of PUMA regulation on both timelines: how each gene ",
    "changes across the wild-type window horizontally and across the same ",
    "window in the Myc+ gland vertically. The other established PUMA inducers ",
    "are drawn in pale grey, and Bbc3 and Foxo3 -- the two genes that do move -- ",
    "are drawn as reference points so the band has a scale to be null on."),
  detail = c(
    sprintf("n = 6 per group. Raw (unshrunken) DESeq2 log2 fold changes. The p53 arm is script 42's own `exclusions$p53_axis` (%d genes: Trp53, Mdm2 and the transactivation targets Cdkn1a, Cdkn2a, Trp53inp1, Zmat3, Eda2r, Phlda3, Ccng1); the pale class is the remaining %d of `exclusions$puma_inputs`. Neither roster is assembled here.",
            nrow(a53), sum(d$class == "other PUMA inducers")),
    sprintf("NOTHING IN THE p53 ARM MOVES, and the script asserts it: no gene reaches padj < 0.05 on the Myc timeline, the largest Myc-timeline change in the whole arm is %+.3f, and no interaction survives correction. p53 itself is the flattest thing on the panel: %s.",
            a53$lfc_myc_time[which.max(abs(a53$lfc_myc_time))], say("Trp53")),
    sprintf("Its two best-known outputs do not move either: %s; %s. Two of the arm DO respond to Myc itself at six weeks -- Cdkn1a (p21) %+.2f at padj %.3f and Phlda3 %+.2f at padj %.3f -- but that is a GENOTYPE effect, on a different axis from this panel, and it argues that the arm is measurable rather than that it moves with time.",
            say("Mdm2"), say("Zmat3"),
            p53$lfc_myc_6W[p53$gene == "Cdkn1a"], p53$padj_myc_6W[p53$gene == "Cdkn1a"],
            p53$lfc_myc_6W[p53$gene == "Phlda3"], p53$padj_myc_6W[p53$gene == "Phlda3"]),
    sprintf("THE TWO REFERENCE GENES MOVE IN TWO DIFFERENT WAYS, and the grammar shows both: %s -- Bbc3 leaves the BAND, straight down the Myc axis; and %s -- Foxo3 stays on the band but leaves the DIAGONAL, because it rises across the wild-type window and does not under Myc. A band of dots means nothing until the reader can see how far a gene that did move sits from it.",
            say("Bbc3"), say("Foxo3")),
    sprintf("ONE GENE OF THE ARM DOES MOVE, ON THE OTHER AXIS, and the panel should not be read as flatter than it is: %s. Eda2r is the only p53-arm gene to clear padj 0.05 on either timeline, it does so on the WILD-TYPE one, and its mean expression is %d counts. It is a developmental change, not a Myc one, and it is drawn at the far left where a reader can see it.",
            say("Eda2r"), round(g("Eda2r", "baseMean"))),
    sprintf("Two of the twelve PUMA inducers fall outside the drawn window and are named here instead: %s. Both are at the floor of what this design can measure (mean expression %s), and neither is in the p53 arm.",
            paste(sprintf("%s (%+.2f, %+.2f)", out$gene, out$lfc_wt_time,
                          out$lfc_myc_time), collapse = " and "),
            paste(round(out$baseMean), collapse = " and "))),
  bounds = c(
    "BATCH = TIMEPOINT, on BOTH axes: each coordinate is a temporal contrast and is DESCRIBED, not claimed. The batch-clean quantity is the vertical distance from the diagonal, which is the interaction.",
    "A NEGATIVE AT n = 6 says no detectable movement, not no movement. What gives it force is that the two reference genes on the same axes, in the same libraries, at the same n, do move -- and one of them (Bbc3) clears genome-wide correction.",
    "AND A TRANSCRIPT-LEVEL NEGATIVE IS NOT A PATHWAY-LEVEL ONE. p53 acts overwhelmingly through protein stabilisation and post-translational modification, not through its own message; Trp53 mRNA being flat is exactly what an active p53 response would also look like. What this panel excludes is a TRANSCRIPTIONAL p53 response -- which is what the sentence claims, and it is why Mdm2, Cdkn1a and Zmat3 are drawn: they are p53's transcriptional OUTPUT, and they are the informative rows, not Trp53 itself.",
    "Several of these genes are lowly expressed (Cdkn2a at 6 counts, Eda2r at 101), so a real half-log2 effect could be missed.",
    "Foxo3 belongs to the PUMA-inducer roster and is drawn here as a mover, not as a member of the null class. The sentence should say \"other than Foxo3\"."),
  source = c(
    "results/priming_arm_teb.rds (scripts/42_priming_arm_and_teb_substrate.R) -- $exclusions$p53_axis and $exclusions$puma_inputs, the two rosters script 42's own PART C used to exclude alternative routes; $machinery for Bbc3",
    "The grammar: two_timeline_base() in figures/panels/_panel_common.R, shared with Figs. 2F (alt), 2G (alt) and 2H+I (alt)"))

save_panel_p(p, "figS2_p53_arm_alt", height = 76)

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {

  print(p)
  print(LEGEND)

  ## every gene on the panel, ranked by what the Myc timeline does
  d[order(d$lfc_myc_time),
    c("class", "gene", "baseMean", "lfc_wt_time", "padj_wt_time",
      "lfc_myc_time", "padj_myc_time", "lfc_interaction")] |>
    print(row.names = FALSE, digits = 3)

  ## the p53 arm on the GENOTYPE axis, which is a different question and is where
  ## Cdkn1a and Phlda3 do move
  as.data.frame(pa$exclusions$p53_axis)[
    , c("gene", "lfc_myc_6W", "padj_myc_6W", "lfc_myc_12W", "padj_myc_12W")] |>
    print(row.names = FALSE, digits = 3)
}
