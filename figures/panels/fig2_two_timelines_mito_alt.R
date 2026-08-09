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
# THE FILL IS THE HORIZONTAL AXIS (author, 2026-08-09): one continuous manuscript
# diverging ramp on the wild-type change, replacing the two hand-named rosters the
# first version coloured by. It is a REDUNDANT encoding by construction -- fill
# and x are the same number -- and that is the point: it sorts the plane by what
# development did without a roster deciding for the reader, so the respiratory arm
# arrives espresso in the lower left and the catabolic arms mint on the right. The
# rosters survive only as the twelve labels, which the sentence names anyway.
# Points are drawn as outlined circles, not solid ones, because zero on this ramp
# is #FAFAFA -- a solid point at the median would be invisible on the page.
#
# NO ARROWS (author, 2026-08-09). The first version drew the drop from the line to
# each labelled arm. The distance is still the interaction and still exact; twelve
# arrowheads over a 143-point cloud just cost more ink than they returned.
#
# ONE RULER, DELIBERATELY. The content ruler is the one the sentence's "logFC"
# names. Drawing both rulers AND both timelines in one 89 mm panel is what made
# the first version unreadable; the mitoPPS ruler gives the same picture and is
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
# Exact to floating point -- the mitoPPS ruler is bit-identical and the content
# ruler differs by 1.2e-16, which is one representable step, not an approximation.
stopifnot(max(abs(d$c_tn + d$c_int - d$c_tp)) < 1e-12,
          max(abs(d$p_tn + d$p_int - d$p_tp)) < 1e-12)

# =============================================================================
# the labelled arms
# =============================================================================
# The respiratory chain, its assembly and translation controls, and the catabolic
# arms are labelled because the sentence names them; everything else is the
# compartment they are read against, and none of it is hidden. There is no
# categorical colour any more, so this roster decides only what gets a name.
LAB <- c("CIV subunits", "CI subunits", "CIII subunits", "CV subunits",
         "OXPHOS assembly factors", "Mitochondrial ribosome",
         "Mitochondrial central dogma", "Fatty acid oxidation",
         "Carnitine shuttle", "Glycine cleavage system", "Amino acid metabolism",
         "Lipid metabolism")
w <- d[match(LAB, d$pathway), ]
stopifnot(!anyNA(w$c_tn), nrow(w) == length(LAB),
          # the claim the panel is built to show
          w$c_tp[w$pathway == "CI subunits"] < w$c_tn[w$pathway == "CI subunits"])

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
stopifnot(identical(unname(tier_labels[["Mitochondrial central dogma"]]),
                    "Central dogma"))       # the declared vector, lower-cased here

# =============================================================================
# the panel
# =============================================================================
LIM <- c(-1, 1) * max(abs(c(d$c_tn, d$c_tp))) * 1.04
at  <- function(f) LIM[1] + diff(LIM) * f

# Both in-panel annotations sit where nothing can, and both are CHECKED rather
# than eyeballed -- in panel fractions, because that is the space a collision
# happens in.
#
# (1) the diagonal's label now runs up the TOP end of the line (author,
#     2026-08-09: the lower-left end sat on the respiratory cluster). Its glyph
#     run is a segment offset just above the line, so the test is a real
#     point-to-segment distance rather than a bounding box: the nearest point is
#     catechol metabolism, which sits BELOW the line and left of where the text
#     starts.
# (2) the quadrant note is right-aligned at mid-width along the bottom edge, so
#     its glyph run reaches back to about a seventh of the width; the lowest
#     thing in that strip is the CIII/CV cluster, a good 8 mm above it.
fx <- (d$c_tn - LIM[1]) / diff(LIM)
fy <- (d$c_tp - LIM[1]) / diff(LIM)
seg_dist <- function(s0, s1) {
  tt <- pmax(0, pmin(1, ((fx - s0[1]) * (s1[1] - s0[1]) +
                         (fy - s0[2]) * (s1[2] - s0[2])) / sum((s1 - s0)^2)))
  sqrt((fx - (s0[1] + tt * (s1[1] - s0[1])))^2 +
       (fy - (s0[2] + tt * (s1[2] - s0[2])))^2)
}
DIAG_AT <- 0.79
stopifnot(
  # the label's glyph run: 14 mm of text lying along the line, offset about one
  # line-height ABOVE it. The nearest thing to it is catechol metabolism, which
  # sits BELOW the line; 0.025 of the panel is 1.7 mm at this size, and the
  # measured gap is 1.9 mm. A tighter bar cannot be met at this font size --
  # the run needs 0.15 of the panel and would clip the corner if pushed further
  # up the line.
  min(seg_dist(c(DIAG_AT, DIAG_AT + 0.012), c(DIAG_AT + 0.15, DIAG_AT + 0.162))) > 0.025,
  # and the bottom strip the quadrant note occupies
  !any(fx > 0.13 & fx < 0.51 & fy < 0.10),
  # and the top-left block the colour key sits in -- generous, because the whole
  # upper-left quadrant is empty for a reason: a pathway there would be one the
  # normal gland strips and Myc restores, and there is none
  !any(fx < 0.40 & fy > 0.60))

# The fill is the horizontal axis. Symmetric about zero so equal ink is equal
# magnitude on both sides -- the asymmetric form of heat_fill() exists for
# lopsided data and this is only mildly lopsided (-0.41 / +0.56).
CLIM <- c(-1, 1) * max(abs(d$c_tn))

p <- ggplot2::ggplot(d, ggplot2::aes(c_tn, c_tp)) +
  two_timeline_base(LIM, diag_at = DIAG_AT,
                    quadrant = "below the line = lost under MYC",
                    quadrant_at = c(0.50, 0.02), quadrant_hjust = 1) +
  # Outlined circles, not solid ones: zero on the manuscript ramp is #FAFAFA, so
  # a filled point at the compartment median would vanish into the page.
  ggplot2::geom_point(ggplot2::aes(fill = c_tn), shape = 21, size = 1.35,
                      colour = "grey45", stroke = 0.12) +
  ggplot2::geom_point(data = w, ggplot2::aes(fill = c_tn), shape = 21,
                      size = 2.1, colour = "grey20", stroke = 0.3) +
  ggrepel::geom_text_repel(data = w, ggplot2::aes(label = lab),
                           colour = "grey15",
                           size = 1.65, seed = 4, max.overlaps = Inf,
                           min.segment.length = 0, segment.size = 0.2,
                           segment.colour = "grey60", box.padding = 0.32,
                           point.padding = 0.18) +
  heat_fill(CLIM, name = "6>12W_wt", breaks = c(-0.4, 0, 0.4)) +
  ggplot2::labs(x = "wild-type 6>12W  (set average log2FC)",
                y = "Myc+ 6>12W  (set average log2FC)") +
  ggplot2::guides(fill = ggplot2::guide_colourbar(
    barwidth = ggplot2::unit(16, "mm"), barheight = ggplot2::unit(2, "mm"),
    ticks.colour = "grey30", frame.colour = "grey30",
    frame.linewidth = 0.2, title.position = "top", title.hjust = 0)) +
  theme_panel(base_size = 6) +
  # THE BAR IS HORIZONTAL AND INSIDE, TOP LEFT (author, 2026-08-09): laid the same
  # way as the axis it repeats, so the eye reads the two as one scale. The quadrant
  # it sits in is the empty one -- x below zero with y above it means a pathway the
  # normal gland strips and Myc restores, and no pathway does that -- which the
  # assertion below states as a fact about the data rather than a hope.
  ggplot2::theme(
    legend.position        = "inside",
    legend.position.inside = c(0.01, 0.98),
    legend.justification   = c(0, 1),
    legend.direction       = "horizontal",
    legend.background      = ggplot2::element_blank(),
    # theme_classic sets axis.text to rel(0.8) of base_size, so the key's own
    # numbers have to be set to 4.8 explicitly or they print LARGER than the axis
    # they repeat -- which is backwards for a key.
    legend.title           = ggplot2::element_text(size = 5.2,
                               margin = ggplot2::margin(b = 0.5, unit = "mm")),
    legend.text            = ggplot2::element_text(size = 4.8),
    legend.margin          = ggplot2::margin(0, 0, 0, 0),
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
    "window in the Myc+ gland on the vertical, and the same wild-type change ",
    "again as the fill. The dashed diagonal is what the pathway would have done ",
    "under development alone, so the vertical distance below it is the ",
    "Myc-specific part of the loss."),
  detail = c(
    sprintf("n = 6 per group; %d MitoPathways (the synthetic mtDNA-encoded pathway excluded, as in Figs. 1E, 1F and 2F). Set-average raw (unshrunken) log2 fold changes on `6>12W_wt` and `6>12W_myc`.",
            nrow(d)),
    "THE DROP FROM THE LINE IS THE INTERACTION, EXACTLY. Script 40's ruler satisfies (Myc+ change) = (wild-type change) + (interaction) to the last bit, asserted in the script, so the vertical distance from the diagonal is not an approximation of the Myc-specific term -- it is that term.",
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
    sprintf("THE mitoPPS RULER GIVES THE SAME PICTURE and is not drawn, to keep the labels legible: OXPHOS subunits %+.3f to %+.3f, compartment median %+.3f to %+.3f, and the same additive identity holds. fig2_wt_mito_contraction.R is the panel that draws both rulers; this one draws both timelines.",
            g("OXPHOS subunits", "p_tn"), g("OXPHOS subunits", "p_tp"),
            stats::median(d$p_tn), stats::median(d$p_tp)),
    sprintf("For scale, the Myc genotype effect at six weeks on the same arm is %+.3f -- the chain is something Myc BUILDS at six weeks and both timelines then take apart.",
            g("OXPHOS subunits", "c_m6"))),
  bounds = c(
    "BATCH = TIMEPOINT, and here it bites BOTH axes: the 6W and 12W cohorts were extracted as two batches, so each coordinate is DESCRIBED, not claimed. What is batch-CLEAN is the VERTICAL DISTANCE from the diagonal -- that is the interaction, and genotype is balanced within each batch. The panel should be read down from the line rather than along the axes.",
    "No per-pathway significance is drawn and none exists for these contrasts in the saved object; script 40 carries adjusted p-values for the genotype contrasts only. The arm-level empirical nulls (script 43, 2000 expression-matched sets) exist for the wild-type axis only and are quoted in fig2_wt_mito_contraction.R's legend.",
    "MitoCarta sets are membership-loose and several named arms are small (glycine cleavage 4 genes, the carnitine shuttle 5). A large move on a four-gene set is one gene, not a module.",
    "THE FILL IS THE HORIZONTAL AXIS -- the same wild-type log2 fold change, on the manuscript diverging ramp with white pinned to zero. It is redundant by construction and carries no second variable: it is there so the plane sorts by what development did without a roster deciding it. Read magnitude off the axes; the ramp is symmetric so equal ink is equal magnitude on both sides.",
    "The twelve named arms are a roster written in this script from the sentence, not a saved classification. They decide what is LABELLED and nothing else -- every one of the 143 pathways is drawn, on the same ramp, and none is hidden."),
  source = c(
    "results/background_vs_myc.rds (scripts/40_background_vs_myc_decomposition.R) -- $ruler: c_tn / c_tp / c_int (content) and p_tn / p_tp / p_int (mitoPPS) for 144 MitoPathways",
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
    ggplot2::geom_point(ggplot2::aes(fill = c_tn), shape = 21, size = 1.35,
                        colour = "grey45", stroke = 0.12) +
    heat_fill(c(-1, 1) * max(abs(d2$c_tn))) +
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
