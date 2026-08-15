# =============================================================================
# fig2_net_state_cross.R -- the young normal gland against the gland at initial
# tumour expansion, and what each of the two steps contributed
# -----------------------------------------------------------------------------
# SLOT: not currently cited. This panel has no sentence yet -- it was built to
# help decide what the reduced Results section says, so its sentence comes after
# it rather than before. `figS1_mb_fork_specificity.R` carries the same string.
#
# THE CONTRAST IS NEW AND IT IS FREE. 6W_wt -> 12W_myc had never been computed:
# script 03 defines seven contrasts and neither group contrast crosses genotype.
# It needs no new model. `~ group` and `~ timepoint*myc_status` are both saturated
# over the same four groups, so gene by gene
#
#     6W_wt>12W_myc  ==  myc_12W + 6>12W_wt  ==  myc_6W + 6>12W_myc
#
# to optimiser precision (script 45 PART A: median deviation 1.7e-08 over 18,523
# genes, 2 genes over 0.01 where the MLE is unstable on a near-empty group). THAT
# IS WHY THE ROWS ARE DRAWN AS A HEAD-TO-TAIL SUM. The two segments are not a
# schematic of a decomposition; they are the decomposition, and their far end is
# the diagonal itself.
#
# WHAT THE PANEL SAYS, in the order the eye meets it:
#   * EVERY ARM IS UP except the lineage signature -- the compartment of the 12W
#     Myc+ gland is larger than the young normal gland's on every mitochondrial
#     arm drawn (bottom strip: 96.5% of 143 MitoPathways above zero, median
#     +0.247).
#   * THE RESPIRATORY CHAIN GAINS LEAST, and it gains least because it is the one
#     arm the developing gland takes back: its grey segment runs to -0.255 and its
#     orange segment returns it to +0.061. Every other arm's two segments point
#     the same way or the grey one barely moves.
#   * AND IT IS THE ONE ARM THAT LOSES ITS ENRICHMENT. In the strip on the right
#     the OXPHOS-subunit row is deeply coloured under myc_6W and myc_12W and pale
#     under the diagonal (NES 2.73 and 3.13 -> 1.27, and the only drawn arm whose
#     diagonal NES misses BH < 0.05).
#
# THE ORDER OF THE TWO SEGMENTS IS THE SUBSTRATE FRAME, AND IT IS A CHOICE.
# Development is drawn first, from zero, and the Myc effect at 12W is drawn on top
# of it -- so the gland's own trajectory is the baseline and the oncogene's effect
# is the increment. The other route (myc_6W then 6>12W_myc) reaches the identical
# endpoint and is reported in the legend block. Drawn the other way round the
# OXPHOS row would read "Myc builds it and development takes it back", which is
# equally true and is the wrong emphasis for a manuscript whose claim is that the
# gland's withdrawal is what permits transformation.
#
# WHAT THIS PANEL MUST NOT BE READ AS SAYING. "The respiratory chain returns to
# baseline" is true only RELATIVE TO THE REST OF THE COMPARTMENT. Against
# expression-matched random gene sets the OXPHOS subunits still rise, at the
# 95.4th percentile (script 45's matched null, 2000 draws) -- they are the lowest
# of the mitochondrial arms, not at chance. The legend block says so in those
# words and the panel does not draw a "no change" mark of any kind.
#
# Reads (read-only, no re-run):
#   results/state_readings.rds  (script 45) -- $arms (the four contrasts + the
#                                  diagonal on the content ruler), $arm_null,
#                                  $ruler (all 144 MitoPathways), $fgsea (six
#                                  rankings), $mitopps, $direction, $identity
# Output: outputs/figures/panels/fig2_net_state_cross.pdf
# =============================================================================

source(here::here("figures", "panels", "_panel_common.R"))

sr_path <- here::here("results", "state_readings.rds")
require_fresher_than(sr_path)
sr <- readRDS(sr_path)

arms <- as.data.frame(sr$arms)
anul <- as.data.frame(sr$arm_null)
rul  <- as.data.frame(sr$ruler)
fg   <- as.data.frame(sr$fgsea)
dir_ <- as.data.frame(sr$direction)
idn  <- as.data.frame(sr$identity)

# --- guards: a stale or half-run object must fail here, not at review ---------
# (1) the identity the head-to-tail drawing rests on
stopifnot(nrow(idn) == 2L, all(idn$median_dev < 1e-5), all(idn$n_over_001 <= 10L))
# (2) the arm-level identity, re-derived rather than trusted
stopifnot(max(abs(arms$c_cross - (arms$c_myc_12W + arms$c_wt_time))) < 1e-6)
# (3) script 43's positive control still reads, so the arms are the published ones
stopifnot(abs(arms$c_wt_time[arms$arm == "OXPHOS subunits"] + 0.2548) < 0.02)
# (4) the diagonal ranking exists among the six
stopifnot("cross" %in% fg$ranking, dplyr::n_distinct(fg$ranking) == 6L)

# =============================================================================
# LEFT -- the decomposition on the ten named arms
# =============================================================================
# The roster is script 43's own (`$defs$arms` plus the pooled proliferation set),
# so this panel and Fig. 2F name the same things. Shortened only where a name will
# not fit a 30 mm label column.
SHORT <- c("TEB vs ductal (HS)"    = "TEB vs ductal",
           "PROLIF_* pooled"       = "proliferation, pooled",
           "amino-acid metabolism" = "amino acid",
           "nucleotide metabolism" = "nucleotide",
           "lipid metabolism"      = "lipid")
arms$lab <- ifelse(arms$arm %in% names(SHORT), SHORT[arms$arm], arms$arm)
arms$lab <- sprintf("%s (%d)", arms$lab, arms$n_genes)
arms <- arms[order(arms$c_cross), ]
arms$y <- seq_len(nrow(arms))
stopifnot(nrow(arms) == 10L)

# The two segments, head to tail. `x0` of the second is the far end of the first,
# which is what makes the pair a sum rather than two independent bars.
seg <- rbind(
  data.frame(y = arms$y, x = 0,               xend = arms$c_wt_time,
             part = "6>12W_wt"),
  data.frame(y = arms$y, x = arms$c_wt_time,  xend = arms$c_cross,
             part = "myc_12W"))
seg$part <- factor(seg$part, levels = c("6>12W_wt", "myc_12W"))

XR <- range(c(0, arms$c_cross, arms$c_wt_time))
XR <- XR + c(-1, 1) * diff(XR) * 0.06

pL <- ggplot2::ggplot() +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey80") +
  ggplot2::geom_segment(
    data = seg, ggplot2::aes(x = x, xend = xend, y = y, yend = y, colour = part),
    linewidth = 1.5, lineend = "butt") +
  # the diagonal itself: the far end of the pair, in the declared net ink
  ggplot2::geom_point(
    data = arms, ggplot2::aes(x = c_cross, y = y),
    shape = 23, size = 1.5, stroke = 0.3,
    fill = unname(contrast_cols[[contrast_net]]),
    colour = unname(contrast_cols[[contrast_net]])) +
  ggplot2::scale_colour_manual(
    values = contrast_cols[c("6>12W_wt", "myc_12W")], name = NULL,
    labels = c("6>12W_wt" = "development, wild type", "myc_12W" = "Myc at 12W")) +
  ggplot2::scale_y_continuous(breaks = arms$y, labels = arms$lab,
                              expand = ggplot2::expansion(add = 0.6)) +
  # ONE X SCALE FOR BOTH PARTS, with this part's tick labels dropped so the axis
  # at the foot serves both -- Fig. 1F's idiom, and here it is what makes the
  # alignment usable: a vertical from any row lands on the strip's own axis.
  ggplot2::scale_x_continuous(limits = XR, labels = NULL,
                              expand = ggplot2::expansion(mult = 0)) +
  ggplot2::labs(x = NULL, y = NULL) +
  theme_panel(base_size = 6) +
  ggplot2::theme(
    axis.ticks.x      = ggplot2::element_blank(),
    legend.position   = "top",
    legend.key.width  = ggplot2::unit(3.2, "mm"),
    legend.key.height = ggplot2::unit(2.2, "mm"),
    legend.margin     = ggplot2::margin(0, 0, -1, 0, "mm"),
    panel.grid.major.y = ggplot2::element_line(linewidth = 0.15, colour = "grey93"),
    plot.margin = ggplot2::margin(1.5, 0.5, 1, 1.5, "mm"))

# =============================================================================
# RIGHT -- the same rows on the enrichment ruler
# =============================================================================
# Three columns, not four: the fourth contrast (6>12W_wt) is already drawn, as the
# grey segment. What the strip adds is the thing a set-average cannot show -- that
# the respiratory arm was a top-ranked enrichment at BOTH ages and is not one on
# the diagonal.
NES_RANKS <- c("myc_6W", "myc_12W", "cross")
NES_LABS  <- c(myc_6W = "myc_6W", myc_12W = "myc_12W", cross = contrast_net)

set_of <- c(
  "OXPHOS subunits"       = "MITOCARTA_OXPHOS_SUBUNITS",
  "OXPHOS (all)"          = "MITOCARTA_OXPHOS",
  "OXPHOS assembly"       = "MITOCARTA_OXPHOS_ASSEMBLY_FACTORS",
  "nucleotide metabolism" = "MITOCARTA_NUCLEOTIDE_METABOLISM",
  "mitoribosome"          = "MITOCARTA_MITOCHONDRIAL_RIBOSOME",
  "TCA cycle"             = "MITOCARTA_TCA_CYCLE",
  "amino-acid metabolism" = "MITOCARTA_AMINO_ACID_METABOLISM",
  "lipid metabolism"      = "MITOCARTA_LIPID_METABOLISM",
  "TEB vs ductal (HS)"    = "MG_TEB_VS_DUCTAL_HS_GRAY_UP")
# PROLIF_* pooled is a UNION of several sets and has no single NES. Its cells are
# left empty rather than filled with a representative set's value.

nes <- do.call(rbind, lapply(arms$arm, function(a) {
  s <- unname(set_of[a])
  do.call(rbind, lapply(NES_RANKS, function(rk) {
    v <- if (is.na(s)) NULL else fg[fg$ranking == rk & fg$pathway == s, ]
    if (is.null(v) || !nrow(v)) {
      data.frame(arm = a, ranking = rk, NES = NA_real_, padj = NA_real_)
    } else {
      v <- v[which.min(v$padj_within_category), ]
      data.frame(arm = a, ranking = rk, NES = v$NES[1],
                 padj = v$padj_within_category[1])
    }
  }))
}))
nes$y <- arms$y[match(nes$arm, arms$arm)]
nes$ranking <- factor(nes$ranking, levels = NES_RANKS)
stopifnot(sum(is.na(nes$NES)) == length(NES_RANKS))   # exactly the pooled row

NLIM <- c(-max(abs(nes$NES), na.rm = TRUE), max(abs(nes$NES), na.rm = TRUE))

# THE NUMBER IS IN THE CELL AND THE FILL REPEATS IT. Same "encoded twice" idiom as
# Fig. 1B and Fig. 2F, where the fill repeats an axis: here it repeats the printed
# value, so no colour bar is drawn and the strip needs no key. Three columns at
# 89 mm cannot afford a bar as well, and a reader comparing 2.73 with 1.27 is
# better served by the numbers than by two shades of mint.
#
# SIGNIFICANCE IS THE FONT WEIGHT, not a glyph: an asterisk beside a number in a
# 4 mm cell collides with it, and bold reads at 5 pt through the base pdf device.
nes$lab  <- ifelse(is.na(nes$NES), "", sprintf("%.2f", nes$NES))
nes$face <- ifelse(!is.na(nes$padj) & nes$padj < 0.05, "bold", "plain")

pR <- ggplot2::ggplot(nes[!is.na(nes$NES), ],
                      ggplot2::aes(x = ranking, y = y, fill = NES)) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
  ggplot2::geom_text(ggplot2::aes(label = lab, fontface = face,
                                  colour = ink_on_fill(NES, NLIM[2])),
                     size = 1.55, show.legend = FALSE) +
  ggplot2::scale_colour_identity() +
  heat_fill(NLIM) +
  ggplot2::guides(fill = "none") +
  ggplot2::scale_x_discrete(labels = NES_LABS, position = "top",
                            expand = ggplot2::expansion(add = 0)) +
  ggplot2::scale_y_continuous(breaks = arms$y, labels = NULL,
                              limits = range(arms$y) + c(-0.6, 0.6),
                              expand = ggplot2::expansion(add = 0)) +
  ggplot2::labs(x = NULL, y = NULL) +
  theme_panel(base_size = 6) +
  # 90 degrees, not 45: a rotated label runs off the right edge of a strip this
  # narrow, and the third column's name is the longest of the three.
  ggplot2::theme(
    axis.text.x  = ggplot2::element_text(angle = 90, hjust = 0, vjust = 0.5,
                                         size = 5),
    axis.text.y  = ggplot2::element_blank(),
    axis.ticks   = ggplot2::element_blank(),
    panel.grid   = ggplot2::element_blank(),
    plot.margin  = ggplot2::margin(1.5, 1.5, 1, 0.5, "mm"))

# =============================================================================
# BOTTOM -- the compartment the ten arms came out of
# =============================================================================
# Fig. 1E's lower-strip construction exactly, on the diagonal instead of on the
# Myc effect at 6W, so the two strips are directly comparable: one point per
# MitoPathway over a kernel density, coloured by sign, zero and the median marked.
# The synthetic mtDNA-encoded pathway is excluded by the same rule scripts 40/43
# and Figs. 1E/1F use; its value is in the legend block.
r143 <- rul[!rul$is_mtdna, ]
stopifnot(nrow(r143) == 143L)
med143 <- stats::median(r143$c_cross)

# WINDOWED FOR DISPLAY, because the strip shares the x axis with the rows above
# and that sharing is the point: a vertical dropped from any arm row lands where
# that pathway sits in the whole compartment. Fig. 2F and Fig. 1G window the same
# way. What falls outside is the membership-loose caveat made visible -- every one
# is a small set, where a set mean's noise scales as one over the square root of n
# -- and all six are on the HIGH side, so the window UNDERSTATES the panel's own
# claim rather than flattering it. Every number in the legend block is computed on
# the complete 143.
outside <- r143[r143$c_cross < XR[1] | r143$c_cross > XR[2], ]
stopifnot(nrow(outside) == 6L, all(outside$c_cross > 0), all(outside$n_genes <= 13L),
          min(r143$c_cross) > XR[1])
drawn143 <- r143[r143$c_cross >= XR[1] & r143$c_cross <= XR[2], ]

set.seed(45)                       # the jitter must not move between rebuilds
dens <- stats::density(r143$c_cross, adjust = 1.1)
dd   <- data.frame(x = dens$x, y = dens$y / max(dens$y))
dd   <- dd[dd$x >= XR[1] & dd$x <= XR[2], ]      # clipped, so the area cannot warn
drawn143$jit <- stats::runif(nrow(drawn143), -0.30, -0.06)
ox_row <- drawn143[drawn143$pathway == "OXPHOS subunits", ]
stopifnot(nrow(ox_row) == 1L)

pB <- ggplot2::ggplot() +
  ggplot2::geom_area(data = dd, ggplot2::aes(x, y), fill = "grey92") +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey55") +
  ggplot2::geom_vline(xintercept = med143, linewidth = 0.35, linetype = "22",
                      colour = "grey25") +
  ggplot2::geom_point(data = drawn143,
                      ggplot2::aes(x = c_cross, y = jit, fill = c_cross),
                      shape = 21, size = 0.85, stroke = 0.12, colour = "grey45") +
  # the arm the panel is about, ringed where it sits in its own compartment
  ggplot2::geom_point(data = ox_row,
                      ggplot2::aes(x = c_cross, y = jit),
                      shape = 21, size = 2.2, stroke = 0.45,
                      fill = NA, colour = unname(contrast_cols[[contrast_net]])) +
  ggplot2::annotate("text", x = med143, y = 1.02, vjust = 0, hjust = 0.5,
                    size = 1.65, colour = "grey25",
                    label = sprintf("median %s", lab_signed(round(med143, 3)))) +
  ggplot2::annotate("text", x = ox_row$c_cross, y = -0.40, vjust = 1, hjust = 0.5,
                    size = 1.65, colour = unname(contrast_cols[[contrast_net]]),
                    label = "OXPHOS subunits") +
  heat_fill(range(r143$c_cross)) +
  ggplot2::guides(fill = "none") +
  ggplot2::scale_x_continuous(limits = XR, labels = lab_signed,
                              expand = ggplot2::expansion(mult = 0)) +
  ggplot2::scale_y_continuous(limits = c(-0.62, 1.18), breaks = NULL,
                              expand = ggplot2::expansion(mult = 0)) +
  ggplot2::labs(x = sprintf("%s, all %d MitoPathways  (set-average log2FC)",
                            contrast_net, nrow(r143)),
                y = NULL) +
  theme_panel(base_size = 6) +
  ggplot2::theme(panel.grid = ggplot2::element_blank(),
                 plot.margin = ggplot2::margin(1, 1.5, 1, 1.5, "mm"))

# =============================================================================
# assembly
# =============================================================================
# The bottom strip carries a SPACER of the strip's own width, so its x axis lines
# up with the rows above. That alignment is load-bearing, not cosmetic: a vertical
# dropped from any arm row has to land where that pathway sits in the compartment,
# and it cannot if the two panels are scaled to different widths.
WIDTHS <- c(1, 0.26)
top <- patchwork::wrap_plots(pL, pR, nrow = 1, widths = WIDTHS)
bot <- patchwork::wrap_plots(pB, patchwork::plot_spacer(), nrow = 1,
                             widths = WIDTHS)
p   <- patchwork::wrap_plots(top, bot, ncol = 1, heights = c(2.9, 1))

# =============================================================================
# the legend text (never drawn)
# =============================================================================
A  <- function(a, col) arms[[col]][arms$arm == a]
NP <- function(a) anul$percentile[anul$arm == a]
NN <- function(a) anul$null_median[anul$arm == a]
NE <- function(a, rk) {
  v <- nes[nes$arm == a & nes$ranking == rk, ]
  if (!nrow(v) || is.na(v$NES)) NA_real_ else v$NES[1]
}
NPADJ <- function(a, rk) {
  v <- nes[nes$arm == a & nes$ranking == rk, ]
  if (!nrow(v) || is.na(v$padj)) NA_real_ else v$padj[1]
}
mpp   <- as.data.frame(sr$mitopps)
MP    <- function(pw) mpp$diff[match(pw, mpp$pathway)]
mtrow <- rul[rul$is_mtdna, ]
n_sig_cross <- sum(fg$ranking == "cross" & fg$padj_within_category < 0.05,
                   na.rm = TRUE)

LEGEND <- panel_legend(
  slot = "not currently cited",
  what = paste0(
    "The young normal gland against the gland at initial tumour expansion ",
    "(6W_wt to 12W_myc), decomposed. Each row is one gene set; the grey segment ",
    "is the wild-type 6 to 12 week change, the orange segment is the Myc effect ",
    "measured at 12 weeks, drawn head to tail from zero, and the black diamond at ",
    "their far end is the diagonal itself. Right, the same rows on the enrichment ",
    "ruler at both ages and on the diagonal: the normalised enrichment score is ",
    "printed in each cell and the fill repeats it, with the number in bold where ",
    "BH-adjusted p < 0.05 within category. Below, the whole mitochondrial compartment on the ",
    "diagonal, one point per MitoPathway over a kernel density."),
  detail = c(
    "n = 6 animals per group. Values are set-average RAW (unshrunken) log2 fold changes, as CLAUDE.md requires of an averaged-LFC visual.",
    sprintf("THE TWO SEGMENTS ARE A SUM, NOT A SCHEMATIC. `~ group` and `~ timepoint*myc_status` are both saturated over the same four groups, so the diagonal is a re-reading of a fit that already exists: gene by gene it equals myc_12W + 6>12W_wt to a median deviation of %.1e over %d genes (%d genes exceed 0.01, where the MLE is unstable on a near-empty group). At set level the same identity holds to %.0e. The contrast required no new model and no new fit.",
            idn$median_dev[1], idn$n_genes[1], idn$n_over_001[1],
            max(abs(arms$c_cross - (arms$c_myc_12W + arms$c_wt_time)))),
    sprintf("THE OTHER ROUTE REACHES THE SAME POINT and is not drawn: myc_6W then 6>12W_myc, identical to %.1e. Development is drawn first because the manuscript's frame is the substrate frame -- the gland's own trajectory is the baseline and the oncogene's effect is the increment on it. Drawn the other way the OXPHOS row would read \"Myc builds it and development takes it back\", which is equally true of the same two numbers.",
            idn$median_dev[2]),
    sprintf("EVERY MITOCHONDRIAL ARM IS UP ON THE DIAGONAL, and the strip below is the scope of that: %.1f%% of the %d MitoPathways sit above zero, median %+.3f. The one negative row is the lineage signature (%s %+.3f), which is a cell-identity readout and not a mitochondrial arm.",
            100 * mean(r143$c_cross > 0), nrow(r143), med143,
            "TEB vs ductal", A("TEB vs ductal (HS)", "c_cross")),
    sprintf("THE RESPIRATORY CHAIN GAINS LEAST, AND THE PANEL SHOWS WHY. OXPHOS subunits: development %+.3f, Myc at 12W %+.3f, diagonal %+.3f. The assembly factors of the same complexes, in the same libraries: development %+.3f, diagonal %+.3f. The difference between those two rows is the whole specificity claim, and it is the same pair Fig. 2F draws on the wild-type window alone.",
            A("OXPHOS subunits", "c_wt_time"), A("OXPHOS subunits", "c_myc_12W"),
            A("OXPHOS subunits", "c_cross"),
            A("OXPHOS assembly", "c_wt_time"), A("OXPHOS assembly", "c_cross")),
    sprintf("AND IT IS THE ONLY ARM THAT LOSES ITS ENRICHMENT. OXPHOS subunits score NES %+.2f on myc_6W and %+.2f on myc_12W, both far beyond BH < 0.05, and %+.2f on the diagonal at BH-adjusted p = %.3f -- the only drawn arm that misses. For contrast, amino-acid metabolism goes %+.2f, %+.2f, %+.2f and stays significant throughout. Across the whole library the diagonal is a strong contrast, not a weak one: %d of %d sets clear BH < 0.05 within category.",
            NE("OXPHOS subunits", "myc_6W"), NE("OXPHOS subunits", "myc_12W"),
            NE("OXPHOS subunits", "cross"), NPADJ("OXPHOS subunits", "cross"),
            NE("amino-acid metabolism", "myc_6W"), NE("amino-acid metabolism", "myc_12W"),
            NE("amino-acid metabolism", "cross"),
            n_sig_cross, sum(fg$ranking == "cross")),
    sprintf("A THIRD RULER AGREES, and it is content-blind. On mitoPPS -- the pairwise-ratio score, which reads how the compartment divides its budget and cancels a uniform scaling -- the diagonal DEMOTES the respiratory arm and promotes the biosynthetic ones: OXPHOS subunits %+.3f, OXPHOS assembly %+.3f, against amino acid %+.3f, lipid %+.3f, mitoribosome %+.3f. So the compartment grows while the respiratory share of its budget falls.",
            MP("OXPHOS subunits"), MP("OXPHOS assembly factors"),
            MP("Amino acid metabolism"), MP("Lipid metabolism"),
            MP("Mitochondrial ribosome")),
    sprintf("EXCLUDED FROM THE STRIP, AND IT IS THE LARGEST SINGLE MOVEMENT: the 13 mtDNA-encoded OXPHOS subunits rise %+.3f on the diagonal. They are dropped by the rule scripts 40 and 43 and Figs. 1E, 1F and 2F all use -- the mtDNA-encoded read fraction is confounded three ways (real content, the proliferation denominator, dissociation leak) and is time-associated rather than genotype-associated, so a contrast that spans time is exactly where it cannot be read.",
            as.numeric(unname(mtrow$c_cross)))),
  bounds = c(
    "BATCH = TIMEPOINT, AND THIS CONTRAST SPANS IT. The 6W and 12W cohorts were extracted as two separate batches, so the diagonal carries the batch offset in full -- the same status as 6>12W_wt, and unlike the genotype-within-age contrasts, which are clean. Every value here is DESCRIBED, not claimed.",
    sprintf("\"RETURNS TO BASELINE\" IS RELATIVE TO THE COMPARTMENT, NOT TO THE TRANSCRIPTOME, and the panel must not be read the other way. Against %d expression-matched random gene sets the OXPHOS subunits still RISE on the diagonal, at the %.1f percentile (null median %+.3f). What is true is that they are the lowest of the mitochondrial arms and the only one that does not clear its compartment: every other drawn arm sits at the %.1f percentile or above. No \"no change\" mark is drawn anywhere on this panel.",
            2000L, NP("OXPHOS subunits"), NN("OXPHOS subunits"),
            min(anul$percentile[!anul$arm %in% c("OXPHOS subunits", "TEB vs ductal (HS)")])),
    "A PERCENTILE IS NOT AN EFFECT SIZE. The pooled proliferation set reaches the 100th percentile of its null on a diagonal value of only +0.057, because 731 genes give a very tight null. Read the percentile as a rank against chance and the segment length as the size.",
    "A NET NEAR ZERO IS NOT AN ABSENCE OF CHANGE. It is two large opposite changes that cancel, which is exactly why both components are drawn and why the diamond is never shown without them.",
    "THE ABUNDANCE-WEIGHTED RULER DISAGREES ON THIS ARM. Every value drawn here is an unweighted per-gene mean. Weighting the same 87 OXPHOS subunits by expression gives +0.201, and the summed normalised counts give +0.226, against +0.061 unweighted -- the highly expressed subunits move up more. The difference is the largest of any arm and it is drawn in its own panel; a sentence quoting one of these numbers must name which ruler it is on.",
    "PROLIF_* pooled is a UNION of several library sets and therefore has no single enrichment score; its three cells in the strip are empty rather than filled with a representative. The lineage row is a Gray signature, not a MitoCarta arm, and is included as the negative direction on the axis.",
    "MitoCarta sets are membership-loose and two drawn arms are small (TCA cycle 20 genes, nucleotide metabolism 34). A large effect on a small set can be one gene, which is why every row label carries its gene count.",
    sprintf("THE LOWER STRIP IS WINDOWED, because it shares its x axis with the rows above and that sharing is the point -- a vertical dropped from any row lands where that pathway sits in the whole compartment. %d of the %d MitoPathways fall outside the drawn range and every one is a small set: %s. All %d are on the HIGH side, so the window UNDERSTATES the compartment's gain rather than flattering it, and nothing falls off the left. Every number in this legend is computed on the complete %d.",
            nrow(outside), nrow(r143),
            paste(sprintf("%s %+.2f (%d genes)", outside$pathway[order(outside$c_cross)],
                          outside$c_cross[order(outside$c_cross)],
                          outside$n_genes[order(outside$c_cross)]), collapse = "; "),
            nrow(outside), nrow(r143)),
    "n = 6 per cell. This panel RANKS arms; no arm-level difference on it is a confirmatory test."),
  source = c(
    "results/state_readings.rds (scripts/45_state_readings.R) -- $arms (four published contrasts plus the diagonal on the content ruler), $arm_null (2000 expression-matched draws), $ruler (all 144 MitoPathways), $fgsea (script 20's five rankings plus the diagonal as a sixth, same recipe and same per-category BH), $mitopps (the diagonal on the mitoPPS ruler), $identity",
    "the arm roster and its set names are script 43's own ($defs$arms in results/substrate_specificity_tradeoff.rds), so this panel and Fig. 2F name the same things",
    "the diagonal contrast: results(dds_group, contrast = c(\"group\", \"12W_pos\", \"6W_neg\"), filterFun = ihw) on results/dds_group_run.rds -- no refit",
    "mitoPPS: Monzel et al. 2025, external/mitotyping/; pairwise-ratio scores on linear-scale DESeq2 normalised counts"))

save_panel_p(p, "fig2_net_state_cross", height = 88)

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {

  print(p)
  print(LEGEND)

  ## the ten arms, both components and the diagonal, with the matched null
  merge(arms[, c("arm", "n_genes", "c_wt_time", "c_myc_12W", "c_cross")],
        anul[, c("arm", "null_median", "percentile")], by = "arm") |>
    (\(x) x[order(x$c_cross), ])() |> print(row.names = FALSE, digits = 3)

  ## the enrichment strip as a table -- the OXPHOS row is the one that fades
  stats::reshape(nes[, c("arm", "ranking", "NES")], idvar = "arm",
                 timevar = "ranking", direction = "wide") |>
    print(row.names = FALSE, digits = 3)

  ## where each arm sits in its own compartment on the diagonal
  q <- function(a) 100 * mean(r143$c_cross < arms$c_cross[arms$arm == a])
  data.frame(arm = arms$arm, c_cross = arms$c_cross,
             compartment_pct = vapply(arms$arm, q, numeric(1))) |>
    print(row.names = FALSE, digits = 3)

  ## the tier medians behind the strip
  as.data.frame(sr$ruler_tiers) |> print(row.names = FALSE, digits = 3)

  ## and the two rulers side by side -- the disagreement this panel flags
  as.data.frame(sr$ruler_compare) |> print(row.names = FALSE, digits = 3)
}
