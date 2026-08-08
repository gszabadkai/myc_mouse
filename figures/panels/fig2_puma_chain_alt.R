# =============================================================================
# fig2_puma_chain_alt.R -- the chain, animal by animal, and the coupling that
# holds it together
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
# then PUMA, then the priming ratio -- and the two current panels each show one
# link of it on a different instrument: a percentile among 8,774 genes, and a
# regression coefficient. Neither shows the chain, and a reader who does not read
# the text cannot assemble one from the other.
#
# PART A PUTS THE WHOLE CHAIN ON ONE PICTURE, ANIMAL BY ANIMAL. Five rows, 24
# columns, each column one mouse, grouped by age and genotype; every row
# standardised across the 24, so a colour says high or low FOR THAT ROW.
#
# WHAT THE FOUR BLOCKS ACTUALLY SAY -- and it is the interaction, not a simple
# co-fall. Group means of the standardised rows:
#
#                 OXPHOS   Foxo3    Bbc3   Bcl2l1   PUMA:Bcl-xL
#   6W_wt          +0.07   -0.90   -0.12    +0.74      -0.55
#   12W_wt         -0.78   +0.89   +0.18    +0.39      -0.16
#   6W_myc         +0.95   +0.01   +0.84    -0.71      +0.94
#   12W_myc        -0.25   -0.00   -0.90    -0.42      -0.23
#
# Read DOWN the two timelines rather than looking for one dark block. In the
# wild-type gland OXPHOS falls (+0.07 -> -0.78), Foxo3 RISES (-0.90 -> +0.89) and
# PUMA is untouched (-0.12 -> +0.18). In the Myc+ gland OXPHOS falls from a much
# higher start (+0.95 -> -0.25), Foxo3 FAILS TO RISE (+0.01 -> -0.00) and PUMA
# COLLAPSES (+0.84 -> -0.90). That is exactly the text -- "although Foxo3 was
# upregulated in the maturing WT gland, its expression was reduced along the
# 6>12W_myc timeline" -- and Bcl2l1 does not follow either row, which is why the
# ratio moves with its numerator.
#
# No p-value is drawn on part A and none is needed: it is the mechanism as
# measured, and every inferential number is in part B or the legend block.
#
# PART B IS THE STATISTIC. The same OXPHOS priority against the same PUMA:Bcl-xL
# ratio, with a line per genotype: the two genotypes couple in OPPOSITE
# directions (wild type -4.80, Myc+ +4.35) and script 43's adjusted interaction
# is +6.09 at p = 0.0052. Drawing it directly beneath the heatmap means the
# reader sees what is being correlated before seeing the correlation.
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
# the per-animal chain
# =============================================================================
ratio <- as.numeric(E["Bbc3", ] - E["Bcl2l1", ])
CH <- list(
  "OXPHOS priority" = ax$oxphos_ppd,
  "Foxo3"           = as.numeric(E["Foxo3", ]),
  "Bbc3"            = as.numeric(E["Bbc3", ]),
  "Bcl2l1"          = as.numeric(E["Bcl2l1", ]),
  "PUMA:Bcl-xL"     = ratio)
ROWS <- names(CH)

# Each row is standardised ACROSS THE 24 ANIMALS, which is the only way five
# quantities in five different units can share a fill scale. It also means the
# fill says "high or low FOR THIS QUANTITY", never "more than that other row".
z <- function(x) as.numeric(scale(x))
H <- do.call(rbind, lapply(ROWS, function(nm)
  data.frame(row = nm, sample = colnames(E), group = as.character(sm$group),
             z = z(CH[[nm]]), stringsAsFactors = FALSE)))
H$row   <- factor(H$row, levels = rev(ROWS))          # first row at the top
H$group <- factor(H$group, levels = names(group_cols))
# Columns keep their on-disk order within a group. Sorting them by one of the
# rows would make the block look more coherent than it is.
H$sample <- factor(H$sample, levels = colnames(E)[order(
  match(as.character(sm$group), names(group_cols)), colnames(E))])
stopifnot(nrow(H) == length(ROWS) * 24L, !anyNA(H$z))

# The reading the panel is built on, asserted: Foxo3 RISES across the wild-type
# window and does not under Myc, while Bbc3 does the opposite. If a re-run
# inverted either, the panel would still draw and the legend would be wrong.
gmz <- function(nm, grp) mean(H$z[H$row == nm & H$group == grp])
stopifnot(gmz("Foxo3", "12W_neg") - gmz("Foxo3", "6W_neg") > 1,
          abs(gmz("Foxo3", "12W_pos") - gmz("Foxo3", "6W_pos")) < 0.2,
          gmz("Bbc3", "12W_pos") - gmz("Bbc3", "6W_pos") < -1,
          abs(gmz("Bbc3", "12W_neg") - gmz("Bbc3", "6W_neg")) < 0.5)

# The fill is CLIPPED AT TWO STANDARD DEVIATIONS rather than stretched to the
# most extreme animal. One 6W_wt mouse sits at |z| = 2.7 on the ratio, and
# scaling to it leaves every other tile pale -- the pattern the panel exists to
# show becomes invisible so that one outlier can be drawn accurately. heat_fill()
# squishes out-of-range values, so nothing is hidden, only saturated; the count
# clipped is reported in the legend block.
# EACH GROUP ALSO GETS ITS MEAN, drawn as a seventh column with a dark border AND
# ITS VALUE PRINTED IN THE CELL. At n = 6 the per-animal tiles are genuinely
# noisy -- that is the honest state of the data and it stays on the panel -- but
# the claim is a group-level interaction, and without the summary the reader has
# to average six tiles by eye in four places.
#
# THE NUMBER IS NOT DECORATION. A mean of six z-scores is necessarily PALER than
# the animals it averages, so on a shared fill scale the summary column is the
# weakest-looking thing on the panel exactly where the claim lives. Printing the
# value fixes that without giving the means a second scale, which would have
# broken the comparison the panel is for. Same idiom as scripts/24's
# A1_imbalance_group_map and scripts/39's B_fgsea_contrasts; the ink is the
# declared ink_on_fill().
MEAN_LAB <- "mean"
Hm <- do.call(rbind, lapply(levels(H$group), function(gp)
  do.call(rbind, lapply(ROWS, function(nm) data.frame(
    row = nm, sample = paste(gp, MEAN_LAB), group = gp,
    z = mean(H$z[H$row == nm & H$group == gp]), stringsAsFactors = FALSE)))))
Hm$row <- factor(Hm$row, levels = levels(H$row))
Hm$group <- factor(Hm$group, levels = levels(H$group))
H$is_mean  <- FALSE
Hm$is_mean <- TRUE
lev <- unlist(lapply(levels(H$group), function(gp) c(
  levels(H$sample)[levels(H$sample) %in% H$sample[H$group == gp]],
  paste(gp, MEAN_LAB))))
H <- rbind(H, Hm)
H$sample <- factor(as.character(H$sample), levels = lev)
stopifnot(!anyNA(H$sample), nlevels(H$sample) == 28L)

ZL <- 2
n_clipped <- sum(abs(H$z) > ZL)
pA <- ggplot2::ggplot(H, ggplot2::aes(sample, row, fill = z)) +
  ggplot2::geom_tile(ggplot2::aes(colour = is_mean, linewidth = is_mean)) +
  ggplot2::scale_colour_manual(values = c(`FALSE` = "white", `TRUE` = "grey15"),
                               guide = "none") +
  ggplot2::scale_linewidth_manual(values = c(`FALSE` = 0.25, `TRUE` = 0.45),
                                  guide = "none") +
  ggplot2::geom_text(data = H[H$is_mean, ],
                     ggplot2::aes(label = sprintf("%+.1f", z)),
                     colour = ink_on_fill(H$z[H$is_mean], ZL), size = 1.5,
                     angle = 90) +   # a 28-column tile is 2.5 mm wide; "+0.95"
                                     # is not, so the number runs up the cell
  ggplot2::facet_grid(~ group, scales = "free_x", space = "free_x",
                      labeller = ggplot2::as_labeller(group_labels)) +
  heat_fill(c(-ZL, ZL), name = NULL, breaks = c(-2, 0, 2)) +
  ggplot2::scale_x_discrete(expand = ggplot2::expansion(mult = 0)) +
  ggplot2::scale_y_discrete(expand = ggplot2::expansion(mult = 0)) +
  ggplot2::labs(x = NULL, y = NULL) +
  # A horizontal bar under the tiles. In a right-hand legend slot a 14 x 1.6 mm
  # bar collapses its own tick labels on top of each other, which is what the
  # first version did.
  ggplot2::guides(fill = ggplot2::guide_colourbar(
    barwidth = ggplot2::unit(16, "mm"), barheight = ggplot2::unit(1.5, "mm"),
    ticks.colour = NA, title = NULL, direction = "horizontal",
    label.position = "bottom")) +
  theme_panel(base_size = 6) +
  ggplot2::theme(
    axis.text.x  = ggplot2::element_blank(),
    axis.ticks   = ggplot2::element_blank(),
    axis.line    = ggplot2::element_blank(),
    axis.text.y  = ggplot2::element_text(size = 6, face = "italic"),
    strip.text   = ggplot2::element_text(face = "plain", size = 6,
                                         margin = ggplot2::margin(0, 0, 0.8, 0, "mm")),
    strip.clip   = "off",
    panel.spacing.x = ggplot2::unit(1.2, "mm"),
    legend.position = "bottom",
    legend.justification = "right",
    legend.margin   = ggplot2::margin(-2, 0, 0, 0, "mm"),
    plot.margin     = ggplot2::margin(1.5, 2.5, 0, 1.5, "mm"))

# =============================================================================
# the coupling, and the reconstruction proved against the record
# =============================================================================
d <- data.frame(y = z(ratio),
                myc = factor(as.character(sm$myc_status), levels = c("neg", "pos")),
                group = factor(as.character(sm$group), levels = names(group_cols)),
                a = ax$oxphos_ppd, epi = pu$epithelial, imm = pu$immune)
m <- summary(stats::lm(y ~ myc * a + epi + imm, d))$coefficients
rec <- tr[tr$outcome == OUT & tr$axis == "oxphos_ppd", ]
per <- tp[tp$outcome == OUT & tp$axis == "oxphos_ppd", ]
stopifnot(nrow(rec) == 1L,
          abs(m["mycpos:a", 1] - rec$myc_x_axis) / abs(rec$myc_x_axis) < 0.05,
          rec$p < 0.01)

pB <- ggplot2::ggplot(d, ggplot2::aes(a, y)) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey85") +
  ggplot2::geom_smooth(ggplot2::aes(colour = myc), method = "lm", formula = y ~ x,
                       se = FALSE, linewidth = 0.5) +
  ggplot2::geom_point(ggplot2::aes(fill = group), shape = 21, size = 1.7,
                      stroke = 0.25, colour = "grey25") +
  ggplot2::scale_colour_manual(values = geno_cols, guide = "none") +
  ggplot2::scale_fill_manual(values = group_cols, labels = group_labels,
                             breaks = names(group_cols), name = NULL) +
  ggplot2::scale_x_continuous(labels = function(x) sprintf("%.1f", x)) +
  ggplot2::scale_y_continuous(labels = lab_signed) +
  ggplot2::labs(x = "OXPHOS priority, per animal  (mitoPPS)",
                y = "PUMA:Bcl-xL  (z)") +
  ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1,
                                               override.aes = list(size = 1.8))) +
  theme_panel(base_size = 6) +
  ggplot2::theme(
    legend.position = "bottom",
    legend.key.size = ggplot2::unit(2.6, "mm"),
    legend.margin   = ggplot2::margin(-1.5, 0, 0, 0, "mm"),
    plot.margin     = ggplot2::margin(1, 2.5, 0.5, 1.5, "mm"))

p <- patchwork::wrap_plots(pA, pB, ncol = 1, heights = c(1, 1.55))

# =============================================================================
# the legend text (never drawn)
# =============================================================================
gm <- function(nm, grp) mean(H$z[H$row == nm & H$group == grp])
slope <- function(g) {
  k <- d$myc == g
  stats::coef(stats::lm(d$y[k] ~ d$a[k]))[2]
}
mg <- as.data.frame(cmo$mech_genes)
ir <- readRDS(here::here("results", "interaction_results.rds"))
ann <- as.data.frame(readRDS(here::here("results", "combined_df_annotated_raw.rds")))
mg$int_p <- as.data.frame(ir$interaction_raw)$pvalue[
  match(ann$gene[match(mg$gene, ann$mgi_symbol)], rownames(as.data.frame(ir$interaction_raw)))]
mg <- mg[order(mg$int_p), ]
stopifnot(identical(mg$gene[1:2], c("Foxo3", "Bbc3")))

LEGEND <- panel_legend(
  slot = "Fig. 2H+I (alt)",
  what = paste0(
    "A. The chain, one column per mouse: the animal's OXPHOS-subunit priority ",
    "score, its Foxo3, Bbc3 and Bcl2l1 levels and the PUMA:Bcl-xL ratio, each ",
    "row standardised across the 24 animals and grouped by age and genotype. ",
    "B. The same OXPHOS priority against the same ratio, with a line fitted ",
    "within each genotype; the statistic in the text is the difference between ",
    "those two slopes."),
  detail = c(
    sprintf("The fill is clipped at two standard deviations: %d of the %d tiles lie beyond it and are drawn saturated rather than shrinking the whole scale to the single most extreme animal. Nothing is hidden.",
            n_clipped, nrow(H)),
    "n = 24 animals, 6 per group. Rows are z-scored ACROSS THE 24, which is the only way five quantities in five different units can share one fill scale -- so a colour says high or low FOR THAT ROW and never compares one row with another. Columns keep their on-disk order within each group; sorting them by one of the rows would make the block look more coherent than it is.",
    sprintf("READ THE BLOCKS DOWN THE TWO TIMELINES, NOT AS ONE DARK CORNER -- what part A shows is the interaction. Group means of the standardised rows. In the WILD-TYPE gland OXPHOS priority falls (%+.2f to %+.2f), Foxo3 RISES (%+.2f to %+.2f) and Bbc3 is untouched (%+.2f to %+.2f). In the Myc+ gland OXPHOS falls from a much higher start (%+.2f to %+.2f), Foxo3 FAILS TO RISE (%+.2f to %+.2f) and Bbc3 COLLAPSES (%+.2f to %+.2f).",
            gm("OXPHOS priority", "6W_neg"), gm("OXPHOS priority", "12W_neg"),
            gm("Foxo3", "6W_neg"), gm("Foxo3", "12W_neg"),
            gm("Bbc3", "6W_neg"), gm("Bbc3", "12W_neg"),
            gm("OXPHOS priority", "6W_pos"), gm("OXPHOS priority", "12W_pos"),
            gm("Foxo3", "6W_pos"), gm("Foxo3", "12W_pos"),
            gm("Bbc3", "6W_pos"), gm("Bbc3", "12W_pos")),
    sprintf("AND Bcl2l1 FOLLOWS NEITHER (%+.2f to %+.2f in the wild type, %+.2f to %+.2f under Myc), which is why the PUMA:Bcl-xL row moves with its numerator: %+.2f to %+.2f under Myc against %+.2f to %+.2f in the wild type.",
            gm("Bcl2l1", "6W_neg"), gm("Bcl2l1", "12W_neg"),
            gm("Bcl2l1", "6W_pos"), gm("Bcl2l1", "12W_pos"),
            gm("PUMA:Bcl-xL", "6W_pos"), gm("PUMA:Bcl-xL", "12W_pos"),
            gm("PUMA:Bcl-xL", "6W_neg"), gm("PUMA:Bcl-xL", "12W_neg")),
    sprintf("B IS AN INTERACTION, WHICH IS TWO SLOPES THAT DIFFER: within the wild-type animals the ratio falls as OXPHOS priority rises (%+.2f) and within the Myc+ animals it rises (%+.2f). Script 43's model (ratio ~ genotype * axis + epithelial + immune) puts the interaction at %+.2f, p = %.4f.",
            slope("neg"), slope("pos"), rec$myc_x_axis, rec$p),
    sprintf("AND FOXO3 IS WHY THE MIDDLE ROW IS THERE: among the %d biogenesis- and cell-death-related genes script 44 curates, Foxo3 has the LOWEST interaction p (%.4f) and Bbc3 the second (%.4f), with a clear gap to the third (%s, %.3f). That is the ranking the text quotes, and it is a far smaller and more meaningful universe than the 8,774-gene scan Fig. 2H uses.",
            nrow(mg), mg$int_p[1], mg$int_p[2], mg$gene[3], mg$int_p[3]),
    sprintf("PART B IS A RECONSTRUCTION AND SAYS SO. Script 43 fits on its own log matrix; this rebuilds the ratio from the VST matrix and reproduces the recorded interaction to %.1f%% (asserted). The number for the text is script 43's.",
            100 * abs(m["mycpos:a", 1] - rec$myc_x_axis) / abs(rec$myc_x_axis)),
    sprintf("Redox is the control and it is not drawn here (fig2_oxphos_puma_coupling.R has it): both genotypes couple to redox strongly and in the SAME direction, so the interaction there is %+.2f at p = %.2f.",
            tr$myc_x_axis[tr$outcome == OUT & tr$axis == "redox_ppd"],
            tr$p[tr$outcome == OUT & tr$axis == "redox_ppd"])),
  bounds = c(
    sprintf("A LEAD, NOT A RESULT. The permutation null in script 43 puts this interaction at the %.1fth percentile of 5,000 draws (empirical p = %.3f), and adding a timepoint term moves p from %.4f to %.4f. At n = 6 per cell an interaction between a genotype and a slope is the least powered thing this design can be asked for.",
            per$percentile, per$p_emp, rec$p, rec$p_with_tp),
    "PART A IS DESCRIPTIVE AND CARRIES NO TEST. It is a display of the measured values, chosen because the chain is what the section claims and no single statistic shows a chain. Every inferential number is in part B or in this block.",
    "THE ROWS ARE NOT INDEPENDENT: PUMA:Bcl-xL is built from the two rows above it, so its column pattern must follow theirs. It is drawn because it is the quantity the text names, not as separate evidence.",
    "mitoPPS is a RELATIVE score -- a high OXPHOS-subunit value means the compartment spends more of its budget there, not that the cell respires more. Neither axis of part B is a rate.",
    "PUMA:Bcl-xL is a transcript ratio, not priming; and the direction of the chain is not established by either part. PUMA restrains the mitochondrial pyruvate carrier (Kim, Cancer Cell 2019), so respiration sits upstream and FOXO3 -> BBC3 closes the circuit rather than starting it."),
  source = c(
    "results/priming_arm_teb.rds (scripts/42) -- $axis_scores (per-animal mitoPPS axes), $purity (the composition covariates the model adjusts for)",
    "results/gsva_scores.rds (scripts/15) -- $expr_mat, the VST matrix the per-animal rows and the ratio are built from",
    "results/substrate_specificity_tradeoff.rds (scripts/43) -- $tradeoff and $tradeoff_perm, the interaction and its permutation null",
    "results/collapse_module_ownership.rds (scripts/44) -- $mech_genes, the 26 biogenesis/cell-death genes the Foxo3 ranking is computed over"))

save_panel_p(p, "fig2_puma_chain_alt", height = 88)

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {

  print(p); print(pA); print(pB)
  print(LEGEND)

  ## the group means of every row -- the block structure as numbers
  stats::aggregate(z ~ row + group, data = H, FUN = mean) |>
    (\(x) stats::reshape(x, idvar = "row", timevar = "group", direction = "wide"))() |>
    print(row.names = FALSE, digits = 2)

  ## the 26 mechanism genes ranked by interaction p -- Foxo3 and Bbc3 on top
  mg[, c("group", "gene", "lfc_6W", "lfc_12W", "int_p", "percentile")] |>
    head(10) |> print(row.names = FALSE, digits = 3)

  ## the reconstruction against the record, and the redox control
  data.frame(refit = m["mycpos:a", 1], recorded = rec$myc_x_axis,
             p_refit = m["mycpos:a", 4], p_recorded = rec$p) |>
    print(row.names = FALSE, digits = 4)
  tr[tr$outcome == OUT, ] |> print(row.names = FALSE, digits = 3)
}
