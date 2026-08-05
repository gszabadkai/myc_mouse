# =============================================================================
# figS2_puma_inducers.R -- the other ways to switch PUMA on, and none of them
# is switched
# -----------------------------------------------------------------------------
# SLOT: Fig. S2D.
#
#   "... no changes in the transcriptome of other known PUMA inducers were
#    observed."
#
# THE JOB. Fig. 2H shows PUMA's Myc effect collapsing between six and twelve
# weeks, and Foxo3 -- its canonical p53-independent activator -- collapsing with
# it. The obvious alternative explanation is that some OTHER upstream input to
# PUMA moved instead, and that Foxo3 is a bystander. This panel is that
# alternative, tested: the twelve known inducers of PUMA, drawn as the Myc effect
# at each age, and eleven of them do nothing at either.
#
# THE ROSTER IS NOT ASSEMBLED HERE. It is script 42's `exclusions$puma_inputs`,
# the list that script's own PART C used to exclude alternative routes: the p53
# family (Trp73), the E2F arm (E2f1), the integrated stress response (Atf4,
# Ddit3, Trib3, Chac1, Eif2ak3, Nupr1, Sesn2) and the FOXO family (Foxo1, Foxo3,
# Foxo4).
#
# FOXO3 IS IN THE ROSTER AND IS DRAWN, because it is the one that is NOT "other".
# Leaving it out would make the panel look like a clean negative when what it
# actually shows is a negative WITH ONE EXCEPTION, and the exception is the
# sentence before it.
#
# WHY THE MYC EFFECT AT EACH AGE RATHER THAN THE INTERACTION. The claim is that
# these genes do not move, and a near-zero interaction can also be two large
# effects that happen to be equal. Drawing both ages shows the level and the
# change at once -- the idiom Fig. S1F uses for the MYC/MAX/MXD network, which is
# the same kind of negative.
#
# Reads (read-only, no re-run):
#   results/priming_arm_teb.rds (script 42) -- $exclusions$puma_inputs
# Output: outputs/figures/panels/figS2_puma_inducers.pdf
# =============================================================================

source(here::here("figures", "panels", "_panel_common.R"))

pa_path <- here::here("results", "priming_arm_teb.rds")
require_fresher_than(pa_path)
pa <- readRDS(pa_path)

d <- as.data.frame(pa$exclusions$puma_inputs)
stopifnot(nrow(d) == 12L,
          all(c("gene", "baseMean", "lfc_myc_6W", "padj_myc_6W", "lfc_myc_12W",
                "padj_myc_12W", "lfc_interaction", "padj_interaction") %in% names(d)),
          "Foxo3" %in% d$gene)

# the panel's whole content, asserted: nothing here is significant under Myc at
# either age, and no interaction survives correction
stopifnot(all(d$padj_myc_6W  >= 0.05, na.rm = TRUE),
          all(d$padj_myc_12W >= 0.05, na.rm = TRUE),
          all(d$padj_interaction >= 0.05, na.rm = TRUE))

# =============================================================================
# the panel
# =============================================================================
# Ordered by the six-week effect so the rows read as a ranking rather than as the
# roster's arbitrary order; Foxo3 is marked because it is the exception the
# previous panel is about.
d$is_foxo3 <- d$gene == "Foxo3"
d$row <- factor(sprintf("%s (%d)", d$gene, round(d$baseMean)),
                levels = sprintf("%s (%d)", d$gene, round(d$baseMean))[order(d$lfc_myc_6W)])

long <- rbind(
  data.frame(row = d$row, gene = d$gene, contrast = contrast_geno[1],
             lfc = d$lfc_myc_6W,  is_foxo3 = d$is_foxo3),
  data.frame(row = d$row, gene = d$gene, contrast = contrast_geno[2],
             lfc = d$lfc_myc_12W, is_foxo3 = d$is_foxo3))
long$contrast <- factor(long$contrast, levels = contrast_geno)

XR <- range(long$lfc) + c(-1, 1) * diff(range(long$lfc)) * 0.10

p <- ggplot2::ggplot(long, ggplot2::aes(lfc, row)) +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey45") +
  # the connector carries the change; the two points carry the levels
  ggplot2::geom_line(ggplot2::aes(group = row), linewidth = 0.3,
                     colour = "grey60") +
  ggplot2::geom_point(ggplot2::aes(colour = contrast), size = 1.4) +
  ggplot2::scale_colour_manual(values = contrast_cols[contrast_geno],
                               breaks = contrast_geno, name = NULL) +
  ggplot2::scale_x_continuous(limits = XR, labels = lab_signed,
                              expand = ggplot2::expansion(mult = 0)) +
  ggplot2::scale_y_discrete(expand = ggplot2::expansion(add = 0.7)) +
  ggplot2::labs(x = "Myc effect  (raw log2FC)", y = NULL) +
  theme_panel(base_size = 6) +
  # Key UNDER the plot, as Fig. S1F does and for the same reason: twelve rows
  # spanning the full width leave no empty corner, and inside it landed on the
  # bottom row. The one gene that is not "other" is marked on its own axis label
  # instead of in the key.
  ggplot2::theme(
    legend.position  = "bottom",
    legend.key.size  = ggplot2::unit(2.6, "mm"),
    legend.margin    = ggplot2::margin(-1.5, 0, 0, 0, "mm"),
    axis.text.y = ggplot2::element_text(
      size = 6, face = ifelse(levels(d$row) %in%
                                sprintf("%s (%d)", "Foxo3",
                                        round(d$baseMean[d$is_foxo3])),
                              "bold.italic", "italic")),
    axis.ticks.y = ggplot2::element_blank(),
    plot.margin  = ggplot2::margin(1.5, 2.5, 1, 1.5, "mm"))

# =============================================================================
# the legend text (never drawn)
# =============================================================================
g  <- function(x, col) d[[col]][d$gene == x]
mx <- d$gene[which.max(abs(d$lfc_myc_6W))]

LEGEND <- panel_legend(
  slot = "Fig. S2D",
  what = paste0(
    "The twelve known upstream inducers of PUMA, each drawn as the Myc genotype ",
    "effect at six weeks and at twelve, joined. The number after each name is ",
    "its mean expression. Foxo3 is set in bold because it is the one gene here ",
    "that is not \"other\" -- it is the subject of the preceding panel."),
  detail = c(
    sprintf("n = 6 per group. Raw (unshrunken) DESeq2 log2 fold changes; %d genes covering the p53 family (Trp73), the E2F arm (E2f1), the integrated stress response (Atf4, Ddit3, Trib3, Chac1, Eif2ak3, Nupr1, Sesn2) and the FOXO family (Foxo1, Foxo3, Foxo4). The roster is script 42's own `exclusions$puma_inputs`, not assembled here.",
            nrow(d)),
    sprintf("NOTHING MOVES, and the script asserts it: no gene reaches padj < 0.05 for the Myc effect at either age (the smallest adjusted p is %.2f at six weeks and %.2f at twelve) and no interaction survives correction (all padj = 1). The largest six-week effect on the panel is %s at %+.2f, on a mean expression of %d.",
            min(d$padj_myc_6W, na.rm = TRUE), min(d$padj_myc_12W, na.rm = TRUE),
            mx, g(mx, "lfc_myc_6W"), round(g(mx, "baseMean"))),
    sprintf("THE ONE EXCEPTION IS THE ONE THE PREVIOUS PANEL IS ABOUT: Foxo3 goes %+.3f at six weeks to %+.3f at twelve, an interaction of %+.3f -- the same sign change as Bbc3, and the reason it is drawn here rather than quietly dropped from a roster it belongs to.",
            g("Foxo3", "lfc_myc_6W"), g("Foxo3", "lfc_myc_12W"),
            g("Foxo3", "lfc_interaction")),
    sprintf("The other two FOXO paralogues do not do it: Foxo1 %+.3f to %+.3f and Foxo4 %+.3f to %+.3f. Neither is a substitute reading of the Foxo3 result, and neither is significant.",
            g("Foxo1", "lfc_myc_6W"), g("Foxo1", "lfc_myc_12W"),
            g("Foxo4", "lfc_myc_6W"), g("Foxo4", "lfc_myc_12W")),
    sprintf("ONE GENE DOES MOVE ON A DIFFERENT AXIS, and it is worth naming so the panel is not read as flatter than it is: Trp73 rises %+.2f across the WILD-TYPE window at padj %.4f. Its mean expression is %d counts, which is at the floor of what this design can measure, and the contrast is the batch-confounded one. It is not a Myc effect and it is not on this panel's axis.",
            g("Trp73", "lfc_wt_time"), g("Trp73", "padj_wt_time"),
            round(g("Trp73", "baseMean")))),
  bounds = c(
    "THIS IS A NEGATIVE AT n = 6, so it says these transcripts show no detectable movement, not that they are unchanged. Its force comes from the contrast with the panels either side of it: the same design, the same libraries and the same n do detect the Myc effects of Fig. 2G and the departures of Fig. 2H.",
    "A TRANSCRIPT-LEVEL NEGATIVE IS NOT A PATHWAY-LEVEL ONE. The integrated stress response and the p53 axis act substantially through protein stability, phosphorylation and localisation; ATF4 in particular is translationally controlled and can be fully active with an unchanged message. What this panel excludes is a transcriptional re-routing of PUMA's inputs, which is the alternative the text raises.",
    "Several of these genes are lowly expressed (Trp73 at 12 counts, Chac1 at 95, Trib3 at 118), so their standard errors are wide and a real half-log2 effect could be missed. The mean expression is on the face of the panel for exactly that reason.",
    "Foxo3 is in this roster and is drawn. The sentence says \"OTHER known PUMA inducers\", so it should name Foxo3 as the exception rather than let the reader infer that the roster excludes it."),
  source = c(
    "results/priming_arm_teb.rds (scripts/42_priming_arm_and_teb_substrate.R) -- $exclusions$puma_inputs, the twelve-gene roster with all four contrasts and the interaction"))

save_panel_p(p, "figS2_puma_inducers", height = 52)

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {

  print(p)
  print(LEGEND)

  ## the roster with every contrast, ranked by the six-week Myc effect
  d[order(d$lfc_myc_6W),
    c("gene", "baseMean", "lfc_myc_6W", "padj_myc_6W", "lfc_myc_12W",
      "padj_myc_12W", "lfc_wt_time", "padj_wt_time", "lfc_interaction")] |>
    print(row.names = FALSE, digits = 3)

  ## the two other rosters script 42 used to exclude alternative routes
  as.data.frame(pa$exclusions$p53_axis)[, c("gene", "baseMean", "lfc_myc_6W",
                                            "padj_myc_6W", "lfc_interaction")] |>
    print(row.names = FALSE, digits = 3)
}
