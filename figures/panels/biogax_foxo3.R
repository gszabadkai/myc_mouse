# =============================================================================
# biogax_foxo3.R -- the one readable TF signal, and the three things it is not
# -----------------------------------------------------------------------------
# DISCUSSION PANEL (see biogax_regulon_split.R for why the `biogax_` prefix).
#
# FOXO3 is the strongest READABLE transcription-factor signal in the wild-type
# timeline: top of the 69 context-free lanes, and the only one of them that is
# neither a mitochondrial-content artifact (its target set contains no OXPHOS
# subunit at all) nor a cell-state artifact (it carries no Gray context). That
# makes it the best-supported TF statement this dataset can make.
#
# It is still not the driver, and this panel is built to make both halves of that
# visible at once rather than letting the first half stand alone.
#
# LEFT -- the lane across the four contrasts, with FOXO1 as the specificity
# control. Two readings: the programme rises with age in BOTH genotypes (so no
# "Myc blocks the developmental rise"), and Myc suppresses it at BOTH ages (a
# stable genotype effect, which is the part worth keeping).
#
# RIGHT -- the programme gene by gene, grouped into the three arms FOXO3 is known
# to control, named in script 47 before the values were read. Only the
# atrophy/turnover arm moves. The arrest arm does not, and neither does the
# apoptotic arm -- including Bbc3, which is the gene the whole death story turns
# on. "The adult gland raises FOXO3 and arms PUMA" is the reading this panel
# exists to prevent.
#
# Reads : results/biogenesis_axis_developmental.rds (script 47 PART I)
# Output: outputs/figures/panels/biogax_foxo3.pdf
# =============================================================================

source(here::here("figures", "panels", "_panel_common.R"))

ba_path <- here::here("results", "biogenesis_axis_developmental.rds")
if (!file.exists(ba_path)) stop("run scripts/47_... first")
ba <- readRDS(ba_path)
if (is.null(ba$foxo_lanes))
  stop("this panel needs script 47 PART I -- re-source scripts/47_... and try again")

fl <- as.data.frame(ba$foxo_lanes)
ft <- as.data.frame(ba$foxo_targets)
fr <- as.data.frame(ba$foxo_reach)

# --- left: the lane across the contrasts, in the declared vocabulary ----------
RANK_MAP <- c(myc_6W = "myc_6W", myc_12W = "myc_12W",
              timepoint_neg = "6>12W_wt", timepoint_pos = "6>12W_myc")
d <- fl[fl$ranking %in% names(RANK_MAP) &
        fl$pathway %in% c("TFT_FOXO3_CHUNG", "TFT_FOXO1_CHUNG"), ]
d$contrast <- factor(unname(RANK_MAP[d$ranking]), levels = rev(contrast_levels))
d$factor   <- factor(ifelse(d$pathway == "TFT_FOXO3_CHUNG", "FOXO3", "FOXO1"),
                     levels = c("FOXO3", "FOXO1"))
d$sig <- d$padj_within_category < 0.05
stopifnot(nrow(d) == 8L)

# ASSERTIONS -- the three readings, so a re-run cannot flip them quietly.
f3 <- d[d$factor == "FOXO3", ]
stopifnot(
  all(f3$NES[f3$contrast %in% contrast_geno] < 0),      # Myc suppresses it
  all(f3$NES[f3$contrast %in% contrast_dev]  > 0),      # both genotypes rise
  fr$covers_ox[fr$set == "TFT_FOXO3_CHUNG"] == 0)       # and it reaches nothing

p1 <- ggplot2::ggplot(d, ggplot2::aes(x = NES, y = contrast)) +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey55") +
  ggplot2::geom_line(ggplot2::aes(group = contrast), linewidth = 0.2,
                     colour = "grey75") +
  ggplot2::geom_point(ggplot2::aes(shape = factor, fill = factor), size = 1.9,
                      stroke = 0.3, colour = "grey20") +
  ggplot2::geom_point(data = d[d$sig, ], shape = 21, size = 3.0, fill = NA,
                      stroke = 0.32, colour = unname(sig_cols[["sig"]])) +
  ggplot2::scale_shape_manual(values = c(FOXO3 = 21, FOXO1 = 24), name = NULL) +
  ggplot2::scale_fill_manual(values = c(FOXO3 = "grey15", FOXO1 = "white"),
                             name = NULL) +
  ggplot2::scale_x_continuous(name = "target-set NES", labels = lab_signed,
                              limits = c(-2.4, 2.4), breaks = seq(-2, 2, 1)) +
  ggplot2::scale_y_discrete(name = NULL) +
  theme_panel() +
  ggplot2::theme(legend.position = "top",
                 legend.margin = ggplot2::margin(0, 0, -3, 0),
                 panel.grid.major.y = ggplot2::element_blank())

# --- right: the programme gene by gene, by arm --------------------------------
ARM_LAB <- c(atrophy_turnover = "atrophy and\nturnover",
             antioxidant      = "antioxidant",
             arrest           = "cell-cycle\narrest",
             apoptotic        = "apoptotic")
ft <- ft[!is.na(ft$wt_lfc), ]
ft$arm <- factor(ft$arm, levels = names(ARM_LAB), labels = ARM_LAB)
ft$sig <- !is.na(ft$wt_padj) & ft$wt_padj < 0.05
ft <- ft[order(ft$arm, ft$wt_lfc), ]
ft$sym <- factor(ft$sym, levels = ft$sym)

stopifnot(ft$wt_lfc[ft$sym == "Fbxo32"] > 1,
          abs(ft$wt_lfc[ft$sym == "Bbc3"]) < 0.15)     # PUMA is flat: the point

p2 <- ggplot2::ggplot(ft, ggplot2::aes(x = wt_lfc, y = sym)) +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey55") +
  ggplot2::geom_segment(ggplot2::aes(x = 0, xend = wt_lfc, yend = sym),
                        linewidth = 0.35, colour = "grey70") +
  ggplot2::geom_point(ggplot2::aes(fill = sig), shape = 21, size = 1.5,
                      stroke = 0.25, colour = "white") +
  ggplot2::scale_fill_manual(
    values = c(`TRUE` = unname(verdict_cols[["rises"]]), `FALSE` = "grey45"),
    guide = "none") +
  ggplot2::facet_grid(arm ~ ., scales = "free_y", space = "free_y", switch = "y") +
  ggplot2::scale_x_continuous(name = "log2 fold change, wild type 6 to 12 weeks",
                              labels = lab_signed) +
  ggplot2::scale_y_discrete(name = NULL, position = "right") +
  theme_panel() +
  ggplot2::theme(
    strip.placement = "outside",
    strip.text.y.left = ggplot2::element_text(angle = 0, hjust = 1, size = 6,
                                              lineheight = 0.95),
    panel.grid.major.y = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_text(face = "italic", size = 5.4))

p <- patchwork::wrap_plots(p1, p2, widths = c(1, 1.12))

g1 <- function(s, col) ft[[col]][ft$sym == s]
LEGEND <- panel_legend(
  slot = "Discussion D5",
  what = paste(
    "Left: the FOXO3 target set across the four contrasts, with FOXO1 as the",
    "specificity control; red rings mark BH-adjusted p below 0.05 within the TF",
    "category. Right: the same programme gene by gene, grouped into the three arms",
    "FOXO3 controls, filled where the adjusted p is below 0.05."),
  detail = c(
    sprintf("The programme rises with age in BOTH genotypes -- %+.2f in wild type and %+.2f in Myc+ -- so Myc does not block a developmental rise. The interaction is not significant.",
            f3$NES[f3$contrast == "6>12W_wt"], f3$NES[f3$contrast == "6>12W_myc"]),
    sprintf("What Myc does is suppress the programme at BOTH ages, %+.2f at six weeks and %+.2f at twelve, and FOXO1 is null in every contrast.",
            f3$NES[f3$contrast == "myc_6W"], f3$NES[f3$contrast == "myc_12W"]),
    sprintf("Only the atrophy and turnover arm moves: Fbxo32 %+.2f (padj %.1e), Bnip3 %+.2f, Sirt1 %+.2f. The arrest arm does not move at all.",
            g1("Fbxo32", "wt_lfc"), g1("Fbxo32", "wt_padj"),
            g1("Bnip3", "wt_lfc"), g1("Sirt1", "wt_lfc")),
    sprintf("Bbc3 (PUMA) is flat in the wild-type gland: %+.3f, padj %.2f. The adult gland does not arm PUMA through FOXO3.",
            g1("Bbc3", "wt_lfc"), g1("Bbc3", "wt_padj"))),
  bounds = c(
    sprintf("FOXO3 REACHES NONE OF THE GENES THAT FELL: its target set contains %d of the 89 OXPHOS subunits, against 59 for the ERRa/NRF1/GABP regulon. It cannot be the proximal cause of the respiratory decline, however clean its own signal is -- it is the readable MARKER of the state, not the driver.",
            fr$covers_ox[fr$set == "TFT_FOXO3_CHUNG"]),
    "The two genotype rows must NOT be read as 'the FOXO3 programme does not attenuate'. NES is scale-free and cannot see an amplitude change; equal NES at both ages is not equal effect size.",
    "mRNA and target-set enrichment are not activity: FOXO3 is regulated by nuclear exclusion. This is concordant evidence, not a measurement of activity.",
    "The set is 39 genes from one curated resource (Chung), and its rise does not clear BH within the TF category (padj 0.10)."),
  source = "results/biogenesis_axis_developmental.rds (script 47 PART I)")

save_panel_p(p, "biogax_foxo3", width = fig_w[["onehalf"]], height = 82)

if (FALSE) {
  print(p)
  ba$foxo_verdict |> print()
  ba$foxo_reach |> print()
  ba$foxo_arm_summary |> print()
  ba$foxo_separability |> print()
}
