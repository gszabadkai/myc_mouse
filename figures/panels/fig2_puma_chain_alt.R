# =============================================================================
# fig2_puma_chain_alt.R -- the chain, group by group, and the coupling that holds
# it together
# -----------------------------------------------------------------------------
# SLOT: Fig. 2H+I (alt). An ALTERNATIVE that MERGES fig2_departure_from_dose.R
# and fig2_oxphos_puma_coupling.R (author's call, 2026-08-09). All three are
# built; the author picks.
#
#   "... the expression of Bbc3 closely mirrored that of Foxo3, its recognized
#    p53-independent activator ... Indeed, the interaction between MYC and OXPHOS
#    subunit coupling to the PUMA/Bcl-XL ratio is highly significant (p = 0.0052)
#    ... intact respiratory state is required for MYC-PUMA mediated mammary
#    epithelial cell death."
#
# WHY MERGE. The section's closing claim is a CHAIN -- respiration, then Foxo3,
# then PUMA -- and the two current panels each show one link of it on a different
# instrument: a percentile among 8,774 genes, and a regression coefficient.
# Neither shows the chain, and a reader who does not read the text cannot assemble
# one from the other.
#
# WHY THE PER-ANIMAL HEATMAP IS GONE (author, 2026-08-09: "the heatmap does not
# work, the individual data are too much variable"). It was right about the data
# and wrong as a display: at n = 6 the per-animal tiles are noisy enough that the
# group pattern has to be averaged by eye in four places at once. The same three
# quantities are now DISTRIBUTIONS in the idiom of Fig. 1E -- box, every animal,
# and the contrast that matters bracketed -- so the noise is still visible but the
# comparison is drawn rather than inferred.
#
# PART A, LEFT: THE CHAIN IN ITS OWN ORDER. Respiration, then Foxo3, then the
# priming ratio -- left to right, so the panel reads the way the mechanism runs.
# Each facet carries the two WITHIN-GENOTYPE temporal contrasts AND, above them,
# the difference between the two -- which is the whole argument in nine brackets:
#
#                    wild type 6>12W    Myc+ 6>12W       difference
#   OXPHOS mitoPPS   -0.11  p 0.108     -0.15  p 0.027   p 0.615
#   Foxo3            +0.46  p 0.0044    -0.00  p 0.982   p 0.013
#   PUMA:Bcl-xL      +0.17  p 0.517     -0.50  p 0.019   p 0.042
#
# i.e. the respiratory arm falls on BOTH timelines; Foxo3 rises ONLY in the
# wild-type gland; the priming ratio collapses ONLY under Myc. That is the text's
# own sentence, and it is a set of interactions rather than a block that falls
# together -- which is what the first version of this panel got wrong until the
# group means were checked.
#
# THE UPPER BRACKET IS THE ONE BATCH = TIMEPOINT LEAVES CLEAN (author, 2026-08-09).
# Both timelines run across the same two extraction batches, so each lower bracket
# is DESCRIBED rather than claimed and what survives the confound is how much the
# two differ. It is drawn from 1.5 to 3.5 -- midpoint of one lower bracket to
# midpoint of the other -- because that is what a difference of differences is.
#
# AND THE GENOTYPE CONTRASTS ARE NOT DRAWN, ON PURPOSE. Myc RAISES Foxo3 at six
# weeks (+0.243, p 0.057) and LOWERS it at twelve (-0.242, p 0.059) by almost
# exactly as much: neither half clears 0.05 and their difference is the strongest
# term in the whole roster (-0.486, p 0.0074, rank 1 of script 44's 26 mechanism
# genes). Labelling either half would print a failed test beside the claim.
#
# PART A, RIGHT: THE TWO MEMBERS OF THE RATIO, as group medians. Bbc3 tracks the
# ratio (+0.74 -> -0.81 under Myc, flat in the wild type); Bcl2l1 follows neither
# timeline (p 0.53 both). So the reversal is entirely the numerator, which is what
# makes the ratio quotable as a PUMA result rather than a balance result.
#
# PART B IS THE STATISTIC. The same OXPHOS mitoPPS score against the same ratio.
# ONE fitted line, through the Myc+ animals, dashed, with its R2 (author,
# 2026-08-09). The wild-type animals trend the OTHER way just as strongly
# (R2 0.33, p 0.052), and it is the DIFFERENCE between the two slopes that script
# 43 tests at p = 0.0052 -- so that p is PRINTED beside the R2 and the ink says
# whose number each is: the R2 in the Myc+ colour, the interaction in neutral.
#
# WHAT IS RECONSTRUCTED AND WHAT IS QUOTED. The per-animal ratio is rebuilt from
# the VST matrix and reproduces script 43's interaction to 1.2% (asserted, the
# same check fig2_oxphos_puma_coupling.R makes). The number for the text is
# script 43's, not the panel's.
#
# READ IT AS A LEAD. The permutation null in the same object puts that
# interaction at the 91.7th percentile of 5,000 draws (p = 0.083) and a timepoint
# term moves p from 0.0052 to 0.088. It is the one axis-by-genotype interaction
# in the corpus that reaches nominal significance.
#
# Reads (read-only, no re-run):
#   results/priming_arm_teb.rds                (script 42) -- $axis_scores, $purity
#   results/gsva_scores.rds                    (script 15) -- $expr_mat
#   results/substrate_specificity_tradeoff.rds (script 43) -- $tradeoff, $tradeoff_perm
#   results/interaction_results.rds            (script 03) -- the DESeq2 record the
#                                                 drawn temporal tests are checked against
#   results/collapse_module_ownership.rds      (script 44) -- $mech_genes, for the
#                                                 ranking quoted in the legend
# Output: outputs/figures/panels/fig2_puma_chain_alt.pdf
# =============================================================================

source(here::here("figures", "panels", "_panel_common.R"))
if (!requireNamespace("patchwork", quietly = TRUE)) stop("fig2H+I alt needs patchwork")

pa_path <- here::here("results", "priming_arm_teb.rds")
ss_path <- here::here("results", "substrate_specificity_tradeoff.rds")
require_fresher_than(pa_path)
require_fresher_than(ss_path)
pa  <- readRDS(pa_path)
ss  <- readRDS(ss_path)
gs  <- readRDS(here::here("results", "gsva_scores.rds"))
cmo <- readRDS(here::here("results", "collapse_module_ownership.rds"))

ax <- as.data.frame(pa$axis_scores)
pu <- as.data.frame(pa$purity)
E  <- gs$expr_mat
sm <- as.data.frame(gs$sample_meta)
tr <- as.data.frame(ss$tradeoff)
tp <- as.data.frame(ss$tradeoff_perm)

OUT <- "Bbc3:Bcl2l1 (PUMA priming)"
stopifnot(identical(as.character(ax$sample), colnames(E)),
          identical(as.character(pu$sample), colnames(E)),
          all(c("Foxo3", "Bbc3", "Bcl2l1") %in% rownames(E)),
          OUT %in% tr$outcome)

# =============================================================================
# the five quantities
# =============================================================================
# "OXPHOS mitoPPS" is the author's name for this axis (2026-08-09). The rename is
# now project-wide across figures/panels/: the quantity is called mitoPPS
# wherever it is drawn OR named in a legend block, including "the mitoPPS ruler"
# where the text used to say "the priority ruler". `priority` survives only as
# the English word for the PHENOMENON -- reprioritisation, a pathway losing
# priority -- because renaming that would make the sentence untrue. Column names
# in the saved objects (`ruler == "priority"`, `$comparator_priority`) are DATA
# and are untouched.
ratio <- as.numeric(E["Bbc3", ] - E["Bcl2l1", ])
CH <- list(
  "OXPHOS mitoPPS" = ax$oxphos_ppd,
  "Foxo3"          = as.numeric(E["Foxo3", ]),
  "PUMA:Bcl-xL"    = ratio,
  "Bbc3"           = as.numeric(E["Bbc3", ]),
  "Bcl2l1"         = as.numeric(E["Bcl2l1", ]))
DIST <- names(CH)[1:3]      # drawn as distributions, in the order the chain runs
HEAT <- names(CH)[4:5]      # drawn as group medians beside them

# Standardised ACROSS THE 24 ANIMALS, which is the only way quantities in three
# different units (a mitoPPS score, a VST level, a log2 ratio) can share one axis
# and one fill. A z-score is a linear transform, so every p-value below is the
# p-value of the untransformed quantity.
z  <- function(x) as.numeric(scale(x))
gl <- names(group_cols)
L  <- do.call(rbind, lapply(names(CH), function(nm)
  data.frame(measure = nm, sample = colnames(E),
             group = factor(as.character(sm$group), levels = gl),
             tp = factor(as.character(sm$timepoint), levels = c("6W", "12W")),
             myc = factor(as.character(sm$myc_status), levels = c("neg", "pos")),
             v = z(CH[[nm]]), stringsAsFactors = FALSE)))
stopifnot(nrow(L) == length(CH) * 24L, !anyNA(L$v))

# =============================================================================
# the two within-genotype temporal contrasts, per quantity
# =============================================================================
# ONE instrument for all three facets: an ordinary least squares fit of the drawn
# values on timepoint, within each genotype -- the same idiom as Fig. 1E's
# fit_simple(). It is the only test available for the two composites (neither the
# mitoPPS axis nor the ratio has a DESeq2 test), so drawing a Wald p on the gene
# facet and a t on the other two would put two instruments on one panel.
# The DESeq2 record is checked against it below and quoted in the legend.
tt <- do.call(rbind, lapply(DIST, function(nm) do.call(rbind, lapply(
  c("neg", "pos"), function(gt) {
    k  <- L$measure == nm & L$myc == gt
    co <- summary(stats::lm(v ~ tp, data = L[k, ]))$coefficients
    data.frame(measure = nm, myc = gt, beta = co[2, 1], p = co[2, 4],
               stringsAsFactors = FALSE)
  }))))
pv <- function(nm, gt) tt$p[tt$measure == nm & tt$myc == gt]

# AND THE DIFFERENCE BETWEEN THE TWO, per quantity (author, 2026-08-09). It is the
# quantity the panel is really about and the only one BATCH = TIMEPOINT leaves
# clean: both timelines run across the same two extraction batches, so what
# survives the confound is how much they differ. Same instrument as the six
# within-genotype tests, so the three numbers in a facet are comparable.
ti <- do.call(rbind, lapply(DIST, function(nm) {
  co <- summary(stats::lm(v ~ tp * myc, data = L[L$measure == nm, ]))$coefficients
  data.frame(measure = nm, beta = co["tp12W:mycpos", 1],
             p = co["tp12W:mycpos", 4], stringsAsFactors = FALSE)
}))
pi_ <- function(nm) ti$p[ti$measure == nm]

# THE READING, ASSERTED. Foxo3 rises across the wild-type window and does not
# under Myc; the ratio does the opposite. If a re-run inverted either, the panel
# would still draw and the legend would be wrong.
bt <- function(nm, gt) tt$beta[tt$measure == nm & tt$myc == gt]
stopifnot(bt("Foxo3", "neg") > 0.3, pv("Foxo3", "neg") < 0.01,
          abs(bt("Foxo3", "pos")) < 0.05, pv("Foxo3", "pos") > 0.5,
          bt("PUMA:Bcl-xL", "pos") < -0.3, pv("PUMA:Bcl-xL", "pos") < 0.05,
          pv("PUMA:Bcl-xL", "neg") > 0.3)

# AND THE DRAWN TEST AGREES WITH THE ANALYSIS OF RECORD. The two composites have
# no DESeq2 test, but Foxo3 does, and so do both members of the ratio: same sign,
# and the same side of 0.05, on both timelines. That is what licenses one
# instrument across the panel.
ir  <- readRDS(here::here("results", "interaction_results.rds"))
ann <- as.data.frame(readRDS(here::here("results", "combined_df_annotated_raw.rds")))
deseq_time <- function(g) {
  ens <- ann$gene[match(g, ann$mgi_symbol)]
  wt  <- as.data.frame(ir$timepoint_neg_raw)[ens, ]
  my  <- as.data.frame(ir$timepoint_pos_raw)[ens, ]
  c(wt_lfc = wt$log2FoldChange, wt_padj = wt$padj,
    myc_lfc = my$log2FoldChange, myc_padj = my$padj)
}
DE <- vapply(c("Foxo3", "Bbc3", "Bcl2l1"), deseq_time, numeric(4))
own <- vapply(c("Foxo3", "Bbc3", "Bcl2l1"), function(g) {
  k <- L$measure == g
  c(wt = stats::coef(stats::lm(v ~ tp, L[k & L$myc == "neg", ]))[2],
    myc = stats::coef(stats::lm(v ~ tp, L[k & L$myc == "pos", ]))[2])
}, numeric(2))
int_de <- vapply(c("Foxo3", "Bbc3", "Bcl2l1"), function(g) {
  r <- as.data.frame(ir$interaction_raw)[ann$gene[match(g, ann$mgi_symbol)], ]
  c(lfc = r$log2FoldChange, p = r$pvalue)
}, numeric(2))
stopifnot(all(sign(own["wt.tp12W", ])  == sign(DE["wt_lfc", ])),
          all(sign(own["myc.tp12W", ]) == sign(DE["myc_lfc", ])),
          identical(DE["wt_padj", ] < 0.05,  c(Foxo3 = TRUE, Bbc3 = FALSE, Bcl2l1 = FALSE)),
          identical(DE["myc_padj", ] < 0.05, c(Foxo3 = FALSE, Bbc3 = TRUE, Bcl2l1 = FALSE)),
          # and the drawn interaction agrees with the record where there is one
          sign(ti$beta[ti$measure == "Foxo3"]) == sign(int_de["lfc", "Foxo3"]),
          ti$p[ti$measure == "Foxo3"] < 0.05, int_de["p", "Foxo3"] < 0.05)

# =============================================================================
# PART A left -- the chain as distributions, Fig. 1E's idiom
# =============================================================================
D <- L[L$measure %in% DIST, ]
D$measure <- factor(D$measure, levels = DIST)

fmt_p <- function(p) if (p < 0.001) "<0.001" else formatC(p, format = "g", digits = 2)
# THREE brackets per facet, in two tiers. The lower two span an AGE PAIR WITHIN A
# GENOTYPE -- groups are drawn genotype-major, so x 1-2 is the wild-type timeline
# and x 3-4 the Myc+ one. The upper one runs from 1.5 to 3.5, i.e. from the
# midpoint of one lower bracket to the midpoint of the other: it brackets the two
# brackets, which is what a difference of differences is. The declared
# bracket_frame()/bracket_layers() draw all three.
cmp <- do.call(rbind, lapply(DIST, function(nm) data.frame(
  measure = factor(nm, levels = DIST),
  x1 = c(1, 3, 1.5), x2 = c(2, 4, 3.5), level = c(1L, 1L, 2L),
  p = c(pv(nm, "neg"), pv(nm, "pos"), pi_(nm)),
  lab = c(fmt_p(pv(nm, "neg")), fmt_p(pv(nm, "pos")),
          paste("interaction", fmt_p(pi_(nm)))),
  stringsAsFactors = FALSE)))
cmp$col <- ifelse(cmp$p < 0.05, "sig", "ns")
brk <- bracket_frame(cmp, range(D$v), pad = 0.035, step = 0.105, tick = 0.020)

pts_layer <- function(dat, size) {
  if (requireNamespace("ggbeeswarm", quietly = TRUE)) {
    ggbeeswarm::geom_quasirandom(data = dat, width = 0.22, size = size, alpha = 0.95)
  } else {
    ggplot2::geom_jitter(data = dat, width = 0.15, height = 0, size = size, alpha = 0.95)
  }
}

pA <- ggplot2::ggplot(D, ggplot2::aes(group, v, colour = group, fill = group)) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey88") +
  ggplot2::geom_boxplot(outlier.shape = NA, width = 0.66, alpha = 0.28,
                        colour = "grey35", linewidth = 0.25) +
  pts_layer(D, size = 0.7) +
  bracket_layers(brk, size = 1.7, linewidth = 0.22) +
  ggplot2::geom_blank(ggplot2::aes(y = y),
                      data = data.frame(y = attr(brk, "headroom"),
                                        group = factor(gl[1], levels = gl),
                                        measure = factor(DIST[1], levels = DIST)),
                      inherit.aes = FALSE) +
  ggplot2::facet_wrap(~ measure, nrow = 1) +
  ggplot2::scale_colour_manual(values = c(group_cols, sig_cols),
                               breaks = names(group_cols), labels = group_labels,
                               name = NULL) +
  ggplot2::scale_fill_manual(values = group_cols, guide = "none") +
  ggplot2::scale_y_continuous(labels = lab_signed) +
  ggplot2::labs(x = NULL, y = "per animal  (z)") +
  ggplot2::guides(colour = ggplot2::guide_legend(
    nrow = 1, override.aes = list(size = 1.5, alpha = 1, shape = 16, label = ""))) +
  theme_panel(base_size = 6) +
  ggplot2::theme(
    axis.text.x     = ggplot2::element_blank(),
    axis.ticks.x    = ggplot2::element_blank(),
    axis.line.x     = ggplot2::element_blank(),
    axis.text.y     = ggplot2::element_text(size = 5.2),
    axis.title.y    = ggplot2::element_text(
                        margin = ggplot2::margin(r = 0.6, unit = "mm")),
    panel.spacing.x = ggplot2::unit(1.4, "mm"),
    strip.clip      = "off",
    strip.text      = ggplot2::element_text(face = "plain", size = 5.6,
                        margin = ggplot2::margin(0, 0, 0.6, 0, "mm")),
    legend.position = "bottom",
    legend.key.size = ggplot2::unit(2.4, "mm"),
    legend.margin   = ggplot2::margin(-1.5, 0, 0, 0, "mm"),
    plot.margin     = ggplot2::margin(1, 1, 0.5, 1.5, "mm"))

# =============================================================================
# PART A right -- the two members of the ratio, as group medians
# =============================================================================
# The MEDIAN, as the author asked, and it is also the right summary here: these
# are six-animal groups and the mean of six z-scores is pulled by exactly the kind
# of single animal that made the per-animal version unreadable. The value is
# PRINTED in the cell, because a group summary on a diverging ramp is necessarily
# paler than the animals it summarises and would otherwise be the weakest mark on
# the panel exactly where the claim lives (the idiom of scripts/24 and 39).
M <- do.call(rbind, lapply(HEAT, function(nm) do.call(rbind, lapply(gl, function(gp) {
  k <- L$measure == nm & L$group == gp
  data.frame(measure = nm, group = gp, med = stats::median(L$v[k]),
             stringsAsFactors = FALSE)
}))))
M$measure <- factor(M$measure, levels = rev(HEAT))     # Bbc3 on top
M$group   <- factor(M$group, levels = gl)
MLIM <- max(abs(M$med)) * 1.02
stopifnot(nrow(M) == 8L, MLIM > 0.5)

# The columns carry NO text. Written out ("12W_myc" rotated under a 6 mm column)
# they cost 10 mm of height, and patchwork then pads the boxplots' blank x axis to
# match, which is where the first version's gap between the two halves came from.
# Instead the four groups are named by a strip of the SAME sample colours the
# boxes beside them use, in the same left-to-right order, keyed by the one legend
# under part A. `colour` is free because the tiles use `fill`, so this needs no
# second fill scale.
pH <- ggplot2::ggplot(M, ggplot2::aes(group, measure, fill = med)) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.4) +
  ggplot2::geom_text(ggplot2::aes(label = sprintf("%+.1f", med)),
                     colour = ink_on_fill(M$med, MLIM), size = 1.6) +
  ggplot2::geom_point(ggplot2::aes(x = group, y = 2.62, colour = group),
                      shape = 15, size = 1.5, inherit.aes = FALSE) +
  heat_fill(c(-MLIM, MLIM), name = NULL, breaks = c(-1, 0, 1)) +
  ggplot2::scale_colour_manual(values = group_cols, guide = "none") +
  ggplot2::scale_x_discrete(expand = ggplot2::expansion(mult = 0)) +
  ggplot2::scale_y_discrete(expand = ggplot2::expansion(add = c(0, 0.72))) +
  # Square cells: without it the two rows stretch to whatever height the row
  # beside them needs, and a 6 x 20 mm tile stops reading as a heat map.
  ggplot2::coord_fixed() +
  ggplot2::labs(x = NULL, y = NULL) +
  ggplot2::guides(fill = ggplot2::guide_colourbar(
    barwidth = ggplot2::unit(13, "mm"), barheight = ggplot2::unit(1.5, "mm"),
    ticks.colour = NA, direction = "horizontal", label.position = "bottom")) +
  theme_panel(base_size = 6) +
  ggplot2::theme(
    axis.text.x     = ggplot2::element_blank(),
    axis.text.y     = ggplot2::element_text(size = 6, face = "italic"),
    axis.ticks      = ggplot2::element_blank(),
    axis.line       = ggplot2::element_blank(),
    legend.position = "bottom",
    legend.text     = ggplot2::element_text(size = 4.8),
    legend.margin   = ggplot2::margin(-1, 0, 0, 0, "mm"),
    plot.margin     = ggplot2::margin(1, 1.5, 0.5, 0.5, "mm"))

# =============================================================================
# PART B -- the coupling, and the reconstruction proved against the record
# =============================================================================
d <- data.frame(y = z(ratio),
                myc = factor(as.character(sm$myc_status), levels = c("neg", "pos")),
                group = factor(as.character(sm$group), levels = gl),
                a = ax$oxphos_ppd, epi = pu$epithelial, imm = pu$immune)
m   <- summary(stats::lm(y ~ myc * a + epi + imm, d))$coefficients
rec <- tr[tr$outcome == OUT & tr$axis == "oxphos_ppd", ]
per <- tp[tp$outcome == OUT & tp$axis == "oxphos_ppd", ]
stopifnot(nrow(rec) == 1L,
          abs(m["mycpos:a", 1] - rec$myc_x_axis) / abs(rec$myc_x_axis) < 0.05,
          rec$p < 0.01)

fit_g <- function(g) stats::lm(y ~ a, data = d[d$myc == g, ])
r2  <- summary(fit_g("pos"))$r.squared
r2w <- summary(fit_g("neg"))$r.squared

# TWO labels, bottom right, right-aligned, and the INK SAYS WHOSE NUMBER IT IS:
# the R2 in the Myc+ colour because it belongs to the dashed line, the interaction
# p in neutral ink because it belongs to both genotypes and is the statistic the
# text quotes. Without the second line the panel would show one slope and let a
# reader take it for the result (author, 2026-08-09). Checked empty rather than
# eyeballed -- the strip covers both lines.
lab_at <- c(x = max(d$a), y = min(d$y) + 0.035 * diff(range(d$y)))
LINE   <- 0.085 * diff(range(d$y))
stopifnot(!any(d$a > lab_at["x"] - 0.32 * diff(range(d$a)) &
               d$y < lab_at["y"] + LINE + 0.10 * diff(range(d$y))))

pB <- ggplot2::ggplot(d, ggplot2::aes(a, y)) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey85") +
  # ONE line, through the Myc+ animals (author, 2026-08-09). The wild-type
  # animals keep their points; their slope is in the legend block, and so is the
  # reason it matters -- see the bounds.
  ggplot2::geom_smooth(data = d[d$myc == "pos", ], method = "lm", formula = y ~ x,
                       se = FALSE, linewidth = 0.5, linetype = "22",
                       colour = unname(geno_cols[["pos"]])) +
  ggplot2::geom_point(ggplot2::aes(fill = group), shape = 21, size = 1.7,
                      stroke = 0.25, colour = "grey25") +
  ggplot2::annotate("text", x = lab_at["x"], y = lab_at["y"] + LINE,
                    hjust = 1, vjust = 0, size = 1.9, parse = TRUE,
                    colour = unname(geno_cols[["pos"]]),
                    label = as.character(as.expression(
                      bquote(Myc * "+ " * R^2 == .(sprintf("%.2f", r2)))))) +
  ggplot2::annotate("text", x = lab_at["x"], y = lab_at["y"], hjust = 1, vjust = 0,
                    size = 1.9, colour = "grey20",
                    label = sprintf("interaction p = %.4f", rec$p)) +
  # No key of its own: part A names the four groups directly above, in the same
  # palette, and a second identical legend would cost 6 mm to say it twice.
  ggplot2::scale_fill_manual(values = group_cols, guide = "none") +
  ggplot2::scale_x_continuous(labels = function(x) sprintf("%.1f", x)) +
  ggplot2::scale_y_continuous(labels = lab_signed) +
  ggplot2::labs(x = "OXPHOS mitoPPS, per animal", y = "PUMA:Bcl-xL  (z)") +
  theme_panel(base_size = 6) +
  ggplot2::theme(
    axis.title.y = ggplot2::element_text(
                     margin = ggplot2::margin(r = 0.6, unit = "mm")),
    plot.margin  = ggplot2::margin(1, 2.5, 0.5, 1.5, "mm"))

p <- patchwork::wrap_plots(
  patchwork::wrap_plots(pA, pH, nrow = 1, widths = c(2.55, 1)),
  pB, ncol = 1, heights = c(1.15, 1))

# =============================================================================
# the legend text (never drawn)
# =============================================================================
gmed <- function(nm, gp) M$med[M$measure == nm & M$group == gp]
slope <- function(g) stats::coef(fit_g(g))[2]
mg <- as.data.frame(cmo$mech_genes)
mg$int_p <- as.data.frame(ir$interaction_raw)$pvalue[
  match(ann$gene[match(mg$gene, ann$mgi_symbol)],
        rownames(as.data.frame(ir$interaction_raw)))]
mg <- mg[order(mg$int_p), ]
stopifnot(identical(mg$gene[1:2], c("Foxo3", "Bbc3")))

say <- function(nm) sprintf(
  "%s %+.2f (p = %s) across the wild-type window and %+.2f (p = %s) across the Myc+ one",
  nm, bt(nm, "neg"), fmt_p(pv(nm, "neg")), bt(nm, "pos"), fmt_p(pv(nm, "pos")))

LEGEND <- panel_legend(
  slot = "Fig. 2H+I (alt)",
  what = paste0(
    "A. The chain in the order it runs: the animal's OXPHOS-subunit mitoPPS ",
    "score, its Foxo3 level and its PUMA:Bcl-xL ratio, one point per mouse, ",
    "each quantity standardised across the 24 animals. The lower brackets are ",
    "the two within-genotype 6-to-12-week contrasts and the upper one, spanning ",
    "their midpoints, is the difference between them. Beside them, the two ",
    "members of the ratio as group medians. B. The same mitoPPS score against ",
    "the same ratio, with a line fitted through the Myc+ animals."),
  detail = c(
    "n = 24 animals, 6 per group. Every quantity is z-scored ACROSS THE 24, which is the only way a mitoPPS score, a VST expression level and a log2 ratio can share one axis and one fill -- so a value says high or low FOR THAT QUANTITY and never compares one facet with another. A z-score is a linear transform, so each p-value is the p-value of the untransformed quantity.",
    sprintf("THE THREE FACETS ARE THREE DIFFERENT ANSWERS, WHICH IS THE POINT. OXPHOS mitoPPS falls on BOTH timelines: %s. Foxo3 rises ONLY in the wild-type gland: %s. And the priming ratio collapses ONLY under Myc: %s.",
            say("OXPHOS mitoPPS"), say("Foxo3"), say("PUMA:Bcl-xL")),
    sprintf("THE UPPER BRACKET IS THE ONE BATCH = TIMEPOINT LEAVES CLEAN, and it is the one the section's claim lives in: Foxo3 %+.2f, p = %s; PUMA:Bcl-xL %+.2f, p = %s; OXPHOS mitoPPS %+.2f, p = %s -- so the two arms of the chain differ between the genotypes and the respiratory arm falls the same way in both. For the two transcripts the DESeq2 record puts the same interaction at Foxo3 %+.3f (p = %.4f) and Bbc3 %+.3f (p = %.4f); Bbc3 is not a facet here but is the numerator of the ratio that is. Neither clears genome-wide BH -- what licenses the test is pre-specification, and Foxo3 and Bbc3 are ranks 1 and 2 of script 44's 26 mechanism genes.",
            ti$beta[ti$measure == "Foxo3"], fmt_p(pi_("Foxo3")),
            ti$beta[ti$measure == "PUMA:Bcl-xL"], fmt_p(pi_("PUMA:Bcl-xL")),
            ti$beta[ti$measure == "OXPHOS mitoPPS"], fmt_p(pi_("OXPHOS mitoPPS")),
            int_de["lfc", "Foxo3"], int_de["p", "Foxo3"],
            int_de["lfc", "Bbc3"], int_de["p", "Bbc3"]),
    sprintf("NEITHER GENOTYPE CONTRAST IS SIGNIFICANT ON ITS OWN AND THAT IS THE RESULT, not a gap in it. Myc RAISES Foxo3 at six weeks (%+.3f, p = %.3f) and LOWERS it at twelve (%+.3f, p = %.3f) by almost exactly as much; each half misses 0.05 and their difference is the strongest term in the roster. A panel that labelled either half would print a failed test beside the claim, so neither is drawn.",
            as.data.frame(ir$myc_6W_raw)[ann$gene[match("Foxo3", ann$mgi_symbol)], "log2FoldChange"],
            as.data.frame(ir$myc_6W_raw)[ann$gene[match("Foxo3", ann$mgi_symbol)], "pvalue"],
            as.data.frame(ir$myc_12W_raw)[ann$gene[match("Foxo3", ann$mgi_symbol)], "log2FoldChange"],
            as.data.frame(ir$myc_12W_raw)[ann$gene[match("Foxo3", ann$mgi_symbol)], "pvalue"]),
    sprintf("THE DRAWN TEST IS ONE INSTRUMENT ACROSS THE PANEL and it agrees with the analysis of record where one exists. Neither the mitoPPS axis nor the ratio has a DESeq2 test, so all six brackets are ordinary least squares on the drawn values (Fig. 1E's fit_simple idiom). For Foxo3 the DESeq2 record is %+.3f, padj %.4f across the wild-type window and %+.3f, padj %.3f across the Myc+ one; for Bbc3 %+.3f, padj %.3f and %+.3f, padj %.3f. Same signs, same side of 0.05, asserted in the script.",
            DE["wt_lfc", "Foxo3"], DE["wt_padj", "Foxo3"],
            DE["myc_lfc", "Foxo3"], DE["myc_padj", "Foxo3"],
            DE["wt_lfc", "Bbc3"], DE["wt_padj", "Bbc3"],
            DE["myc_lfc", "Bbc3"], DE["myc_padj", "Bbc3"]),
    sprintf("THE REVERSAL IS ALL NUMERATOR. Bbc3 group medians run %+.2f, %+.2f in the wild type and %+.2f, %+.2f under Myc, while Bcl2l1 runs %+.2f, %+.2f and %+.2f, %+.2f and moves on NEITHER timeline (DESeq2 padj %.2f and %.2f). That is what makes the ratio quotable as a PUMA result rather than a balance result.",
            gmed("Bbc3", "6W_neg"), gmed("Bbc3", "12W_neg"),
            gmed("Bbc3", "6W_pos"), gmed("Bbc3", "12W_pos"),
            gmed("Bcl2l1", "6W_neg"), gmed("Bcl2l1", "12W_neg"),
            gmed("Bcl2l1", "6W_pos"), gmed("Bcl2l1", "12W_pos"),
            DE["wt_padj", "Bcl2l1"], DE["myc_padj", "Bcl2l1"]),
    sprintf("B IS ONE HALF OF AN INTERACTION. The drawn line is the Myc+ animals, slope %+.2f, R2 %.2f, p = %.3f. The wild-type animals are drawn as points and trend the OTHER way just as strongly (slope %+.2f, R2 %.2f, p = %.3f). Script 43's model (ratio ~ genotype * axis + epithelial + immune) tests the DIFFERENCE between the two slopes: %+.2f, p = %.4f.",
            slope("pos"), r2, summary(fit_g("pos"))$coefficients[2, 4],
            slope("neg"), r2w, summary(fit_g("neg"))$coefficients[2, 4],
            rec$myc_x_axis, rec$p),
    sprintf("AND FOXO3 IS WHY THE MIDDLE FACET IS THERE: among the %d biogenesis- and cell-death-related genes script 44 curates, Foxo3 has the LOWEST interaction p (%.4f) and Bbc3 the second (%.4f), with a clear gap to the third (%s, %.3f). That is the ranking the text quotes, and it is a far smaller and more meaningful universe than the 8,774-gene scan Fig. 2H uses.",
            nrow(mg), mg$int_p[1], mg$int_p[2], mg$gene[3], mg$int_p[3]),
    sprintf("PART B IS A RECONSTRUCTION AND SAYS SO. Script 43 fits on its own log matrix; this rebuilds the ratio from the VST matrix and reproduces the recorded interaction to %.1f%% (asserted). The number for the text is script 43's.",
            100 * abs(m["mycpos:a", 1] - rec$myc_x_axis) / abs(rec$myc_x_axis)),
    sprintf("Redox is the control and it is not drawn here (fig2_oxphos_puma_coupling.R has it): both genotypes couple to redox strongly and in the SAME direction, so the interaction there is %+.2f at p = %.2f.",
            tr$myc_x_axis[tr$outcome == OUT & tr$axis == "redox_ppd"],
            tr$p[tr$outcome == OUT & tr$axis == "redox_ppd"])),
  bounds = c(
    sprintf("PART B DRAWS ONE OF THE TWO SLOPES AND PRINTS THE STATISTIC FOR BOTH. The wild-type animals are points only, and they trend the opposite way at R2 %.2f (p = %.3f) -- as strong as the drawn line. The two numbers on the panel are therefore about different things and the ink says so: the R2 in the Myc+ colour is the dashed fit, the neutral p = %.4f is the DIFFERENCE between the two slopes, which is what the text quotes. Without reading the second number a viewer sees a Myc-specific correlation where the result is a genotype-by-slope interaction.",
            r2w, summary(fit_g("neg"))$coefficients[2, 4], rec$p),
    sprintf("A LEAD, NOT A RESULT. The permutation null in script 43 puts that interaction at the %.1fth percentile of 5,000 draws (empirical p = %.3f), and adding a timepoint term moves p from %.4f to %.4f. At n = 6 per cell an interaction between a genotype and a slope is the least powered thing this design can be asked for.",
            per$percentile, per$p_emp, rec$p, rec$p_with_tp),
    "BATCH = TIMEPOINT: the six brackets in part A are all WITHIN-genotype temporal contrasts, and the 6W and 12W cohorts were extracted as two batches. Each bracket is therefore DESCRIBED, not claimed; what is batch-clean is the CONTRAST BETWEEN the two brackets in a facet, because genotype is balanced within each batch.",
    "THE THREE FACETS ARE NOT INDEPENDENT MEASUREMENTS: PUMA:Bcl-xL is built from the two transcripts in the heatmap beside it, so its pattern must follow theirs. It is drawn because it is the quantity the text names, not as separate evidence.",
    "mitoPPS is a RELATIVE score -- a high OXPHOS-subunit value means the compartment spends more of its budget there, not that the cell respires more. Neither axis of part B is a rate.",
    "PUMA:Bcl-xL is a transcript ratio, not priming; and the direction of the chain is not established by any part of this panel. PUMA restrains the mitochondrial pyruvate carrier (Kim, Cancer Cell 2019), so respiration sits upstream and FOXO3 -> BBC3 closes the circuit rather than starting it."),
  source = c(
    "results/priming_arm_teb.rds (scripts/42) -- $axis_scores (per-animal mitoPPS axes), $purity (the composition covariates the model adjusts for)",
    "results/gsva_scores.rds (scripts/15) -- $expr_mat, the VST matrix the per-animal values and the ratio are built from",
    "results/substrate_specificity_tradeoff.rds (scripts/43) -- $tradeoff and $tradeoff_perm, the interaction and its permutation null",
    "results/interaction_results.rds (scripts/03) -- the raw DESeq2 temporal contrasts the drawn tests are checked against",
    "results/collapse_module_ownership.rds (scripts/44) -- $mech_genes, the 26 biogenesis/cell-death genes the Foxo3 ranking is computed over"))

save_panel_p(p, "fig2_puma_chain_alt", height = 100)

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {

  print(p); print(pA); print(pH); print(pB)
  print(LEGEND)

  ## the six drawn contrasts
  tt |> print(row.names = FALSE, digits = 3)

  ## the drawn test against the DESeq2 record, for the three genes that have one
  data.frame(gene = colnames(DE), own_wt = own["wt.tp12W", ],
             deseq_wt = DE["wt_lfc", ], deseq_wt_padj = DE["wt_padj", ],
             own_myc = own["myc.tp12W", ], deseq_myc = DE["myc_lfc", ],
             deseq_myc_padj = DE["myc_padj", ]) |>
    print(row.names = FALSE, digits = 3)

  ## the group medians drawn in the heatmap, and the means for comparison
  stats::aggregate(v ~ measure + group, data = L,
                   FUN = function(x) c(median = stats::median(x), mean = mean(x))) |>
    print(row.names = FALSE, digits = 2)

  ## the 26 mechanism genes ranked by interaction p -- Foxo3 and Bbc3 on top
  mg[, c("group", "gene", "lfc_6W", "lfc_12W", "int_p", "percentile")] |>
    head(10) |> print(row.names = FALSE, digits = 3)

  ## both slopes, the reconstruction against the record, and the redox control
  data.frame(slope_wt = slope("neg"), slope_myc = slope("pos"),
             refit_int = m["mycpos:a", 1], recorded_int = rec$myc_x_axis,
             p_recorded = rec$p) |> print(row.names = FALSE, digits = 4)
  tr[tr$outcome == OUT, ] |> print(row.names = FALSE, digits = 3)
}
