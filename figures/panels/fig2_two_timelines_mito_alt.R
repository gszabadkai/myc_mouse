# =============================================================================
# fig2_two_timelines_mito_alt.R -- development strips the respiratory chain, and
# Myc strips it again
# -----------------------------------------------------------------------------
# SLOT: Fig. 2F (alt). An ALTERNATIVE to fig2_wt_mito_contraction.R, not a
# replacement -- both are built and the author picks.
#
#   "... the normal gland development showed a specific mitochondrial
#    transcriptomic decline pattern. All OXPHOS subunit logFC and MitoPPS dropped
#    significantly, while the overall mitochondrial transcriptome abundance
#    remained relatively stable ... the maturing gland upregulated amino-acid and
#    lipid catabolism ... While OXPHOS subunits were downregulated in the MYC
#    6->12W timeline in accordance with the reduced Myc effect, the normal
#    maturation of the mammary gland added A FURTHER REDUCTION in OXPHOS in the WT
#    6->12W timeline."
#
# WHY THE CURRENT PANEL IS NOT ENOUGH. fig2_wt_mito_contraction.R draws the
# wild-type window on two rulers, which is correct and answers the first half of
# the sentence. It cannot show the second half at all, because the Myc+ timeline
# is not on it -- and the second half is the one the section turns on.
#
# THE PLANE. x = what the wild-type gland does over the window, y = what the Myc+
# gland does. The dashed diagonal is DEVELOPMENT ALONE. The vertical drop from it
# to a point is the Myc-specific term, and that is exact rather than approximate:
# script 40's ruler satisfies c_tp = c_tn + c_int to the last bit (asserted).
#
#   lower left            down in BOTH timelines -- and on the content ruler this
#                         is where the respiratory chain sits, nearly alone
#   right of zero,        up in development, down under Myc -- the catabolic arms
#   below zero
#   below the diagonal    93% of the 143 pathways: the global Myc fade
#
# SO THE PANEL SAYS, WITHOUT A SENTENCE: development RAISES the compartment
# (median +0.041, 72% above zero) while STRIPPING OXPHOS; the Myc+ gland lowers
# everything; and OXPHOS ends deepest because the two terms ADD. OXPHOS subunits
# -0.255 in development, -0.150 more from Myc, -0.405 in the Myc+ gland.
#
# ONE RULER, DELIBERATELY. The content ruler is the one the sentence's "logFC"
# names. Drawing both rulers AND both timelines in one 89 mm panel is what made
# the first version unreadable; the priority ruler gives the same picture and is
# quoted in the legend (OXPHOS subunits -0.109 -> -0.153, compartment median
# +0.003 -> -0.023), and the sandbox redraws the panel on it with one constant.
#
# BATCH = TIMEPOINT (CLAUDE.md): BOTH axes are the batch-confounded contrast, so
# every value is DESCRIBED, not claimed. What is batch-CLEAN is the vertical
# distance from the diagonal -- it is the interaction, and genotype is balanced
# within each batch. Say it that way round.
#
# Reads (read-only, no re-run):
#   results/background_vs_myc.rds (script 40) -- $ruler, both rulers x four
#                                    contrasts x 144 MitoPathways
# Output: outputs/figures/panels/fig2_two_timelines_mito_alt.pdf
# =============================================================================

source(here::here("figures", "panels", "_panel_common.R"))

bv_path <- here::here("results", "background_vs_myc.rds")
require_fresher_than(bv_path)
r <- as.data.frame(readRDS(bv_path)$ruler)

# Script 40 saves a tibble and its content columns carry names; strip them or a
# row lookup returns a named vector and sprintf prints the name too (Fig. 1E).
num <- function(x) as.numeric(unname(x))
r$n <- as.integer(num(r$n_genes))
for (v in c("c_tn", "c_tp", "c_int", "p_tn", "p_tp", "p_int", "c_m6"))
  r[[v]] <- num(r[[v]])

# 143, not 144: the synthetic mtDNA-encoded pathway is dropped exactly as script
# 40 drops it, so this panel and every regression in the corpus describe the same
# compartment.
d <- r[!r$is_mtdna, ]
stopifnot(nrow(d) == 143L)

# THE DEVICE, PROVED: a vertical drop from the diagonal IS the interaction.
# Exact to floating point -- the priority ruler is bit-identical and the content
# ruler differs by 1.2e-16, which is one representable step, not an approximation.
stopifnot(max(abs(d$c_tn + d$c_int - d$c_tp)) < 1e-12,
          max(abs(d$p_tn + d$p_int - d$p_tp)) < 1e-12)

# =============================================================================
# the three classes and the labelled arms
# =============================================================================
# The respiratory chain and the catabolic arms are named here because the
# sentence names them; everything else is the compartment they are read against.
# `programme_group()` is not the grouping (it is built for the library's fGSEA
# categories, not for MitoCarta tiers), and `tier_cols` is not used: four of its
# seven hues are the sample palette, and this panel has no samples on it.
RESP  <- c("OXPHOS subunits", "CI subunits", "CII subunits", "CIII subunits",
           "CIV subunits", "CV subunits", "Complex I", "Complex III",
           "Complex IV", "Complex V", "OXPHOS")
CATAB <- c("Fatty acid oxidation", "Carnitine shuttle",
           "Carnitine synthesis and transport", "Glycine cleavage system",
           "Amino acid metabolism", "Lipid metabolism", "Pyruvate metabolism",
           "Branched-chain amino acid metabolism")
d$class <- ifelse(d$pathway %in% RESP, "respiratory chain",
                  ifelse(d$pathway %in% CATAB, "catabolic arms", "other"))
d$class <- factor(d$class, levels = c("respiratory chain", "catabolic arms", "other"))

LAB <- c("CIV subunits", "CI subunits", "CIII subunits", "CV subunits",
         "OXPHOS assembly factors", "Mitochondrial ribosome",
         "Mitochondrial central dogma", "Fatty acid oxidation",
         "Carnitine shuttle", "Glycine cleavage system", "Amino acid metabolism",
         "Lipid metabolism")
w <- d[match(LAB, d$pathway), ]
stopifnot(!anyNA(w$c_tn), nrow(w) == length(LAB),
          # the claim the panel is built to show
          w$c_tp[w$pathway == "CI subunits"] < w$c_tn[w$pathway == "CI subunits"],
          sum(d$class == "respiratory chain") >= 10L)

# Shorter names where a MitoCarta name will not fit beside a point; the
# central-dogma shortening is taken from the declared `tier_labels`.
SHORT <- c("OXPHOS assembly factors"     = "OXPHOS assembly",
           "Mitochondrial ribosome"      = "mitoribosome",
           "Mitochondrial central dogma" = "central dogma",
           "Glycine cleavage system"     = "glycine cleavage",
           "Amino acid metabolism"       = "amino acid",
           "Lipid metabolism"            = "lipid",
           "Fatty acid oxidation"        = "fatty acid ox.",
           "Carnitine shuttle"           = "carnitine")
w$lab <- ifelse(w$pathway %in% names(SHORT), SHORT[w$pathway], w$pathway)
# A labelled arm must be READABLE even when its class ink is the pale "other"
# grey -- and the three arms in that class here are the controls (the assembly
# factors, the mitoribosome, the central dogma), which are the point of the
# panel, not background. Their TEXT therefore takes a dark neutral while their
# POINT keeps the class colour. Passed as a per-row constant rather than an
# aesthetic, so the panel keeps one colour scale.
w$lab_col <- ifelse(w$class == "other", "grey25",
                    ifelse(w$class == "respiratory chain",
                           unname(pole_cols[["down"]]), unname(pole_cols[["up"]])))
stopifnot(identical(unname(tier_labels[["Mitochondrial central dogma"]]),
                    "Central dogma"))       # the declared vector, lower-cased here

# =============================================================================
# the panel
# =============================================================================
LIM <- c(-1, 1) * max(abs(c(d$c_tn, d$c_tp))) * 1.04

# The two in-panel annotations sit where nothing can: the diagonal's label on the
# empty lower-left run of the line (nothing falls that far in BOTH timelines
# except CIV, which is well above it), and the quadrant note in the corner
# beneath. Checked against the data rather than eyeballed.
at <- function(f) LIM[1] + diff(LIM) * f
stopifnot(
  !any(abs(d$c_tn - at(0.13)) < diff(LIM) * 0.10 &
       abs(d$c_tp - at(0.13)) < diff(LIM) * 0.06),
  !any(d$c_tn > at(0.34) & d$c_tn < at(0.66) & d$c_tp < at(0.10)))

p <- ggplot2::ggplot(d, ggplot2::aes(c_tn, c_tp)) +
  # the diagonal's label goes to the empty lower-left run of the line, and the
  # quadrant note to the corner beneath it; both regions are checked below.
  two_timeline_base(LIM, diag_at = 0.13,
                    quadrant = "below the line = lost under Myc as well",
                    quadrant_at = c(0.50, 0.02), quadrant_hjust = 0.5) +
  # the Myc-specific term, drawn: from where development alone would have put the
  # arm, down to where it actually is. Only for the named arms -- 143 arrows is a
  # hedgehog, and the 93% figure belongs in the legend.
  ggplot2::geom_segment(data = w, inherit.aes = FALSE,
                        ggplot2::aes(x = c_tn, xend = c_tn, y = c_tn, yend = c_tp),
                        arrow = ggplot2::arrow(length = ggplot2::unit(1, "mm"),
                                               type = "closed"),
                        linewidth = 0.22, colour = "grey55") +
  ggplot2::geom_point(ggplot2::aes(colour = class), size = 1.15, stroke = 0,
                      alpha = 0.9) +
  ggplot2::geom_point(data = w, ggplot2::aes(colour = class), size = 1.9,
                      stroke = 0) +
  ggrepel::geom_text_repel(data = w, ggplot2::aes(label = lab),
                           colour = w$lab_col,
                           size = 1.65, seed = 4, max.overlaps = Inf,
                           min.segment.length = 0, segment.size = 0.2,
                           segment.colour = "grey60", box.padding = 0.32,
                           point.padding = 0.18) +
  ggplot2::scale_colour_manual(
    values = c("respiratory chain" = unname(pole_cols[["down"]]),
               "catabolic arms"    = unname(pole_cols[["up"]]),
               "other"             = unname(pole_cols[["other"]])),
    breaks = c("respiratory chain", "catabolic arms"), name = NULL) +
  ggplot2::labs(x = "wild-type 6>12W  (set-average log2FC)",
                y = "Myc+ 6>12W") +
  ggplot2::guides(colour = ggplot2::guide_legend(
    override.aes = list(size = 1.8, alpha = 1))) +
  theme_panel(base_size = 6) +
  # Key inside, top left: nothing gains under Myc while losing in development, so
  # the wedge above the diagonal on the left is provably empty.
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
g   <- function(pw, col) d[[col]][d$pathway == pw]
arm <- function(pw, nm = pw)
  sprintf("%s %+.3f in development and %+.3f under Myc, a Myc-specific %+.3f",
          nm, g(pw, "c_tn"), g(pw, "c_tp"), g(pw, "c_int"))

LEGEND <- panel_legend(
  slot = "Fig. 2F (alt)",
  what = paste0(
    "Every mitochondrial pathway on both timelines at once: the change across ",
    "the wild-type window on the horizontal axis, the change across the same ",
    "window in the Myc+ gland on the vertical. The dashed diagonal is what the ",
    "pathway would have done under development alone, so the arrow from it down ",
    "to a point is the Myc-specific part of the loss."),
  detail = c(
    sprintf("n = 6 per group; %d MitoPathways (the synthetic mtDNA-encoded pathway excluded, as in Figs. 1E, 1F and 2F). Set-average raw (unshrunken) log2 fold changes on `6>12W_wt` and `6>12W_myc`.",
            nrow(d)),
    sprintf("THE ARROW IS THE INTERACTION, EXACTLY. Script 40's ruler satisfies (Myc+ change) = (wild-type change) + (interaction) to the last bit, asserted in the script, so the vertical drop from the diagonal is not an approximation of the Myc-specific term -- it is that term.",
            NULL),
    sprintf("THE RESPIRATORY ARM IS THE ONLY THING IN THE LOWER-LEFT QUADRANT: %s; %s; %s; %s. Development takes about a quarter of a log2 out of the chain and Myc takes a further sixth.",
            arm("CIV subunits"), arm("CI subunits"), arm("CIII subunits"),
            arm("CV subunits")),
    sprintf("AND THE CONTROL IS ON THE PANEL: the assembly factors of those same complexes sit at the origin in development (%+.3f) and fall only by the global Myc term (%+.3f). What development strips is the stoichiometry of the chain, not the machinery that builds it -- and the Myc+ gland does not add anything specific to them either.",
            g("OXPHOS assembly factors", "c_tn"),
            g("OXPHOS assembly factors", "c_tp")),
    sprintf("THE COMPARTMENT AS A WHOLE RISES IN DEVELOPMENT AND FALLS UNDER MYC, which is what makes the respiratory arm's position mean something: median %+.3f across the wild-type window with %.0f%% above zero, against median %+.3f and %.0f%% under Myc. %.0f%% of the %d pathways sit below the diagonal -- that is the global fade of the Myc effect, not a mitochondrial result.",
            stats::median(d$c_tn), 100 * mean(d$c_tn > 0),
            stats::median(d$c_tp), 100 * mean(d$c_tp > 0),
            100 * mean(d$c_int < 0), nrow(d)),
    sprintf("CATABOLISM IS THE OTHER CORNER: %s; %s; %s; %s. They rise in development and are pulled below zero only by the Myc term, which is a different thing from being repressed.",
            arm("Fatty acid oxidation", "fatty acid oxidation"),
            arm("Carnitine shuttle", "the carnitine shuttle"),
            arm("Glycine cleavage system", "glycine cleavage"),
            arm("Amino acid metabolism", "amino acid metabolism")),
    sprintf("THE PRIORITY RULER GIVES THE SAME PICTURE and is not drawn, to keep the labels legible: OXPHOS subunits %+.3f to %+.3f, compartment median %+.3f to %+.3f, and the same additive identity holds. fig2_wt_mito_contraction.R is the panel that draws both rulers; this one draws both timelines.",
            g("OXPHOS subunits", "p_tn"), g("OXPHOS subunits", "p_tp"),
            stats::median(d$p_tn), stats::median(d$p_tp)),
    sprintf("For scale, the Myc genotype effect at six weeks on the same arm is %+.3f -- the chain is something Myc BUILDS at six weeks and both timelines then take apart.",
            g("OXPHOS subunits", "c_m6"))),
  bounds = c(
    "BATCH = TIMEPOINT, and here it bites BOTH axes: the 6W and 12W cohorts were extracted as two batches, so each coordinate is DESCRIBED, not claimed. What is batch-CLEAN is the VERTICAL DISTANCE from the diagonal -- that is the interaction, and genotype is balanced within each batch. The panel should be read down from the line rather than along the axes.",
    "No per-pathway significance is drawn and none exists for these contrasts in the saved object; script 40 carries adjusted p-values for the genotype contrasts only. The arm-level empirical nulls (script 43, 2000 expression-matched sets) exist for the wild-type axis only and are quoted in fig2_wt_mito_contraction.R's legend.",
    "MitoCarta sets are membership-loose and several named arms are small (glycine cleavage 4 genes, the carnitine shuttle 5). A large move on a four-gene set is one gene, not a module.",
    "The two poles of the manuscript diverging ramp are spent here as CATEGORICAL colours -- espresso for the respiratory arm, mint for the catabolic one -- as Figs. 1G and 1H spend them. Neither is a sample colour.",
    "\"Respiratory chain\" and \"catabolic arms\" are rosters named in this script from the sentence, not a saved classification; every other pathway is drawn in grey and none is hidden."),
  source = c(
    "results/background_vs_myc.rds (scripts/40_background_vs_myc_decomposition.R) -- $ruler: c_tn / c_tp / c_int (content) and p_tn / p_tp / p_int (mitoPPS priority) for 144 MitoPathways",
    "The grammar: two_timeline_base() in figures/panels/_panel_common.R, shared with Figs. 2G (alt), 2H+I (alt) and S2D (alt)"))

save_panel_p(p, "fig2_two_timelines_mito_alt", height = 76)

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {

  print(p)
  print(LEGEND)

  ## the named arms, both timelines and the Myc-specific term
  w[, c("pathway", "n", "c_tn", "c_tp", "c_int", "p_tn", "p_tp", "p_int")] |>
    (\(x) x[order(x$c_tp), ])() |> print(row.names = FALSE, digits = 3)

  ## the same panel on the PRIORITY ruler -- one substitution
  d2 <- d; d2$c_tn <- d$p_tn; d2$c_tp <- d$p_tp
  w2 <- w; w2$c_tn <- w$p_tn; w2$c_tp <- w$p_tp
  L2 <- c(-1, 1) * max(abs(c(d2$c_tn, d2$c_tp))) * 1.04
  ggplot2::ggplot(d2, ggplot2::aes(c_tn, c_tp)) +
    two_timeline_base(L2) +
    ggplot2::geom_point(ggplot2::aes(colour = class), size = 1.2, stroke = 0) +
    ggrepel::geom_text_repel(data = w2, ggplot2::aes(label = lab), size = 1.8,
                             seed = 4, max.overlaps = Inf) +
    theme_panel(base_size = 6)

  ## who else is in the lower-left quadrant besides the chain
  d[d$c_tn < 0 & d$c_tp < 0, c("pathway", "tier", "n", "c_tn", "c_tp")] |>
    (\(x) x[order(x$c_tp), ])() |> print(row.names = FALSE, digits = 3)

  ## and the few pathways ABOVE the diagonal -- where Myc adds rather than fades
  d[d$c_int > 0, c("pathway", "tier", "n", "c_tn", "c_tp", "c_int")] |>
    (\(x) x[order(-x$c_int), ])() |> print(row.names = FALSE, digits = 3)
}
