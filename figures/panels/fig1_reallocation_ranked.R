# =============================================================================
# fig1_reallocation_ranked.R -- what Myc promotes and demotes INSIDE the
# mitochondrial compartment, and how much of it survives to twelve weeks
# -----------------------------------------------------------------------------
# SLOT: Fig. 1F. Cited twice, which is what fixes its form.
#
#   paragraph 3: "Second, a resource reallocation altered intra-compartmental
#     priorities, as quantified by MitoPPS analysis (Monzel et al. 2025)
#     (Fig. 1F). MYC significantly promoted protein import/homeostasis,
#     translation, OXPHOS subunits (all complexes), and biosynthetic pathways
#     including glycine cleavage, pyruvate, and serine metabolism. Conversely, it
#     demoted dynamics and surveillance, including fission, mitophagy, and
#     apoptosis alongside calcium signalling."
#   paragraph 4: "a similar effect on the mitopathway priority scores was
#     observed (see Fig. 1F)."
#
# So the panel must carry (a) the named winners and losers, (b) their statistics,
# and (c) the 6W -> 12W fade in the SAME panel, because paragraph 4 points back
# at it rather than at a panel of its own.
#
#   TOP     the thirteen pathways the sentence names, ranked, each with its 6W and
#           12W Myc effect joined by a connector. Filled = padj < 0.05 at 6W.
#   BOTTOM  all 143 non-mtDNA MitoPathways at 6W, the same construction as Fig.
#           1E's lower strip and DELIBERATELY so: on the content ruler that
#           distribution is one-sided (95% above zero), here it straddles zero
#           (57% above). Two panels, one pair of rulers, and the difference
#           between the two strips IS the two-mechanism sentence.
#
# WHY THE ROWS ARE NAMED RATHER THAN ALL 144 RANKED. The sentence names twelve
# things; a 144-point ranked chart with twelve leader labels is unreadable at 89
# mm, and unlabelled it cannot support the sentence at all. The cherry-picking
# objection is answered by the lower strip, which shows the full distribution the
# thirteen are drawn from. figures/figS4_reallocation_full.R has the labelled
# all-pathway version if a supplementary ever wants it.
#
# WHY n IS PRINTED IN EVERY ROW LABEL. MitoCarta sets are membership-loose and
# several of the largest effects here rest on three or four genes (glycine
# cleavage n = 4, serine n = 4). A big effect on a four-gene set is one gene, not
# a module -- so the reader is given n on the face of the panel rather than in the
# legend.
#
# THE RULER. mitoPPS (Monzel 2025) is pairwise-ratio and therefore CONTENT-BLIND:
# each pathway is read against the rest of the compartment, so a uniform rise in
# mitochondrial content cancels. "Demoted" here means "rose more slowly than the
# compartment", never "fell" -- all thirteen rows rise in absolute content, the
# smallest by +0.117 log2. That is the whole reason Fig. 1E has to be read first.
#
# Reads (read-only, no re-run):
#   results/background_vs_myc.rds (script 40) -- $ruler, which carries BOTH rulers
#       for all 144 MitoPathways: p_m6 / p_m12 (priority) and c_m6 (content), plus
#       tier, n_genes and p_m6_padj. Fig. 1E's lower strip reads the same object.
#   results/mitopps_scores.rds (script 08) -- read for a provenance assertion only
#       (see below), never gated.
# Output: outputs/figures/panels/fig1_reallocation_ranked.pdf
# =============================================================================

source(here::here("figures", "panels", "_panel_common.R"))
if (!requireNamespace("patchwork", quietly = TRUE)) stop("Fig. 1F needs patchwork")

bg    <- readRDS(here::here("results", "background_vs_myc.rds"))
ruler <- as.data.frame(bg$ruler)                # script 40 saves tibbles

# --- provenance: the ruler's priority column IS script 08's mitoPPS contrast ---
# Asserted rather than assumed, because the two objects are written by different
# scripts three weeks apart and a silent divergence would be invisible.
#
# NOT gated by require_fresher_than(). results/mitopps_scores.rds is timestamped
# ONE MINUTE before results/gsva_scores.rds although both came out of the same
# post-reconciliation re-run of 2026-07-24, so the freshness guard would stop this
# panel for a reason that is an artefact of the ordering inside that session. The
# guard is right in general and wrong here; the identity check below is stronger.
mp <- readRDS(here::here("results", "mitopps_scores.rds"))
e6 <- mp$mitopps_pairwise[mp$mitopps_pairwise$contrast == "Myc_effect_6W", ]
i  <- match(ruler$pathway, e6$pathway)
stopifnot(!anyNA(i),
          max(abs(ruler$p_m6      - e6$diff[i])) == 0,
          max(abs(ruler$p_m6_padj - e6$padj[i])) == 0)

# =============================================================================
# TOP -- the thirteen the sentence names
# =============================================================================
# Keys are the MitoCarta pathway names on disk; the second column is what is
# drawn, shortened only where a single column cannot hold the full name. The
# three tier-level rows take their shortening from `tier_labels`.
rows <- rbind(
  data.frame(pathway = "Glycine cleavage system",                 label = "Glycine cleavage"),
  data.frame(pathway = "Pyruvate metabolism",                     label = "Pyruvate metabolism"),
  data.frame(pathway = "Serine metabolism",                       label = "Serine metabolism"),
  data.frame(pathway = "Translation",                             label = "Mitochondrial translation"),
  data.frame(pathway = "Protein import, sorting and homeostasis",
             label = unname(tier_labels[["Protein import, sorting and homeostasis"]])),
  data.frame(pathway = "OXPHOS subunits",                         label = "OXPHOS subunits"),
  data.frame(pathway = "OXPHOS assembly factors",                 label = "OXPHOS assembly factors"),
  data.frame(pathway = "Calcium homeostasis",                     label = "Calcium homeostasis"),
  data.frame(pathway = "Small molecule transport",
             label = unname(tier_labels[["Small molecule transport"]])),
  data.frame(pathway = "Mitochondrial dynamics and surveillance",
             label = unname(tier_labels[["Mitochondrial dynamics and surveillance"]])),
  data.frame(pathway = "Mitophagy",                               label = "Mitophagy"),
  data.frame(pathway = "Fission",                                 label = "Fission"),
  data.frame(pathway = "Apoptosis",                               label = "Apoptosis"),
  stringsAsFactors = FALSE)

d <- ruler[match(rows$pathway, ruler$pathway), ]
stopifnot(nrow(d) == nrow(rows), !anyNA(d$pathway), !any(d$is_mtdna))
d$label <- sprintf("%s (%d)", rows$label, d$n_genes)
d <- d[order(d$p_m6), ]                          # ascending -> promoted at top
d$label <- factor(d$label, levels = d$label)
d$sig   <- ifelse(!is.na(d$p_m6_padj) & d$p_m6_padj < 0.05, "yes", "no")

# The connector needs the two effects in one long frame so a single colour scale
# carries the contrast vocabulary.
dl <- rbind(
  data.frame(label = d$label, eff = d$p_m6,  contrast = "myc_6W",  sig = d$sig),
  data.frame(label = d$label, eff = d$p_m12, contrast = "myc_12W", sig = "12W"))
dl$contrast <- factor(dl$contrast, levels = contrast_geno)

# ONE x SCALE FOR BOTH PARTS. The two parts plot the same quantity, so separate
# scales would let a reader compare a row against the distribution and be wrong by
# whatever the two ranges happened to differ by. Fixed limits over everything
# drawn, patchwork aligns the panel regions, and the top part therefore drops its
# tick labels -- the axis at the foot serves both, and a vertical dropped from any
# row lands where that pathway sits in the compartment.
XLIM <- range(c(ruler$p_m6[!ruler$is_mtdna], ruler$p_m12[!ruler$is_mtdna],
                d$p_m6, d$p_m12))
XLIM <- XLIM + c(-1, 1) * diff(XLIM) * 0.05

p_top <- ggplot2::ggplot(dl, ggplot2::aes(x = eff, y = label)) +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey55") +
  ggplot2::geom_segment(data = d, inherit.aes = FALSE,
                        ggplot2::aes(x = p_m6, xend = p_m12, y = label, yend = label),
                        colour = "grey72", linewidth = 0.3) +
  # 12W first so the 6W point, which carries the significance, sits on top
  ggplot2::geom_point(data = dl[dl$contrast == "myc_12W", ],
                      ggplot2::aes(colour = contrast), size = 0.9, shape = 16) +
  ggplot2::geom_point(data = dl[dl$contrast == "myc_6W", ],
                      ggplot2::aes(colour = contrast, shape = sig),
                      size = 1.5, fill = "white", stroke = 0.4) +
  ggplot2::scale_colour_manual(values = contrast_cols[contrast_geno],
                               breaks = contrast_geno, name = NULL) +
  ggplot2::scale_shape_manual(values = c(yes = 16, no = 21),
                              breaks = c("yes", "no"),
                              labels = c("padj < 0.05", "n.s."), name = NULL) +
  ggplot2::scale_x_continuous(labels = lab_signed, limits = XLIM,
                              expand = ggplot2::expansion(mult = c(0, 0))) +
  ggplot2::labs(x = NULL, y = NULL) +
  ggplot2::guides(
    colour = ggplot2::guide_legend(order = 1, override.aes = list(size = 1.5, shape = 16)),
    shape  = ggplot2::guide_legend(order = 2, override.aes = list(size = 1.5,
                                                                  colour = "grey25"))) +
  theme_panel(base_size = 6) +
  ggplot2::theme(
    axis.text.y  = ggplot2::element_text(size = 5.6),
    axis.line.y  = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    axis.text.x  = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank(),
    plot.margin  = ggplot2::margin(1, 1.5, 0.5, 1.5, "mm"))

# =============================================================================
# BOTTOM -- the compartment the thirteen came out of
# =============================================================================
# Identical construction to Fig. 1E's lower strip, on the other ruler. The
# synthetic mtDNA-encoded pathway carries is_mtdna and is excluded exactly as
# script 40 excludes it from every fit and null.
r143 <- ruler[!ruler$is_mtdna, ]
rs   <- as.data.frame(bg$ruler_summary)
rs6  <- rs[rs$ruler == "priority" & rs$metric == "m6", ]

pct_up_143 <- 100 * mean(r143$p_m6 > 0)
med_143    <- stats::median(r143$p_m6)

stopifnot(nrow(r143) == 143L, nrow(ruler) == 144L,
          abs(100 * mean(ruler$p_m6 > 0) - rs6$pct_up) < 1e-6,
          abs(stats::median(ruler$p_m6)  - rs6$median) < 1e-6)

r143$dir <- ifelse(r143$p_m6 > 0, "up", "down")
dn <- stats::density(r143$p_m6, adjust = 0.9)
PT <- -max(dn$y) * 0.14

p_bot <- ggplot2::ggplot(r143, ggplot2::aes(x = p_m6)) +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey55") +
  ggplot2::geom_density(adjust = 0.9, fill = "grey88", colour = "grey40",
                        linewidth = 0.3) +
  ggplot2::geom_segment(x = med_143, xend = med_143, y = 0, yend = max(dn$y) * 1.02,
                        linewidth = 0.3, colour = "grey20", linetype = "22") +
  ggplot2::geom_jitter(ggplot2::aes(y = PT, colour = dir), height = max(dn$y) * 0.06,
                       width = 0, size = 0.55, alpha = 0.85, show.legend = FALSE) +
  ggplot2::scale_colour_manual(values = direction_cols) +
  ggplot2::scale_x_continuous(labels = lab_signed, limits = XLIM,
                              expand = ggplot2::expansion(mult = c(0, 0))) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.16, 0.06))) +
  ggplot2::labs(x = "Myc priority effect at 6W, per MitoPathway (mitoPPS)", y = NULL) +
  theme_panel(base_size = 6) +
  ggplot2::theme(
    axis.text.y  = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    axis.line.y  = ggplot2::element_blank(),
    plot.margin  = ggplot2::margin(0.5, 1.5, 1, 1.5, "mm"))

p <- patchwork::wrap_plots(p_top, p_bot, ncol = 1, heights = c(3, 1) ) +
  patchwork::plot_layout(guides = "collect") &
  ggplot2::theme(legend.position = "bottom",
                 legend.margin   = ggplot2::margin(-2, 0, 0, 0),
                 legend.box.spacing = ggplot2::unit(1, "mm"),
                 legend.key.size = ggplot2::unit(2.4, "mm"))

# =============================================================================
# the legend text (never drawn)
# =============================================================================
f3 <- function(x) sprintf("%+.3f", x)
row_line <- function(pw) {
  r <- d[d$pathway == pw, ]
  sprintf("%s: %s at 6W (padj %.3f), %s at 12W",
          gsub(" \\(\\d+\\)$", "", as.character(r$label)),
          f3(r$p_m6), r$p_m6_padj, f3(r$p_m12))
}
n_sig <- sum(!is.na(ruler$p_m6_padj) & ruler$p_m6_padj < 0.05)

LEGEND <- panel_legend(
  slot = "Fig. 1F",
  what = paste0(
    "Myc reprioritises the mitochondrial compartment. TOP: the Myc genotype ",
    "effect on mitoPPS priority for the thirteen MitoPathways named in the text, ",
    "ranked, with the twelve-week effect joined to the six-week one so the fade ",
    "reads in the same panel. BOTTOM: the same effect across all 143 ",
    "nuclear-encoded MitoPathways -- the distribution the thirteen are drawn from, ",
    "and the counterpart of Fig. 1E's lower strip on the other ruler."),
  detail = c(
    "mitoPPS (Monzel et al. 2025) is a pairwise-ratio score: every pathway is normalised against the rest of the mitochondrial compartment, so the quantity is RELATIVE priority and a uniform change in mitochondrial content cancels. The contrast is Myc+ minus wild type within one age, which is the clean axis (genotype is balanced within each extraction batch). n = 6 animals per group.",
    "Filled points are padj < 0.05 at six weeks, open points are not; padj is script 08's within-contrast Benjamini-Hochberg adjustment across all 144 pathways. Row labels carry the number of genes in the pathway.",
    sprintf("Of the 144 pathways, %d clear padj < 0.05 at six weeks. Among the thirteen drawn: %s.",
            n_sig, paste(gsub(" \\(\\d+\\)$", "", as.character(d$label[d$sig == "yes"])),
                         collapse = ", ")),
    row_line("Glycine cleavage system"), row_line("Pyruvate metabolism"),
    row_line("Serine metabolism"), row_line("Translation"),
    row_line("Protein import, sorting and homeostasis"),
    row_line("OXPHOS subunits"), row_line("OXPHOS assembly factors"),
    row_line("Calcium homeostasis"), row_line("Small molecule transport"),
    row_line("Mitochondrial dynamics and surveillance"), row_line("Mitophagy"),
    row_line("Fission"), row_line("Apoptosis"),
    sprintf("BOTTOM: one point per MitoPathway coloured by sign, over a kernel density; the dashed line is the median (%+.4f) and the solid line is zero. %.1f%% of the 143 are above zero. Set against Fig. 1E's lower strip, where 95.1%% are above zero on the content ruler, this is the two-mechanism sentence in two pictures: Myc raises almost the whole compartment and reorders it at the same time.",
            med_143, pct_up_143),
    sprintf("The compartment-wide spread narrows with age as everything else does: the 6W priority effects have SD %.3f and the 12W effects %.3f, and the regression of one on the other has slope 0.644 with R2 0.789 (Fig. 1G).",
            rs$sd[rs$ruler == "priority" & rs$metric == "m6"],
            rs$sd[rs$ruler == "priority" & rs$metric == "m12"])),
  bounds = c(
    "\"DEMOTED\" NEVER MEANS \"FELL\". The ruler is content-blind, so a negative value means the pathway rose more slowly than the compartment around it. Every one of the thirteen rows rises in absolute content at six weeks, the smallest of them by +0.117 log2 (Apoptosis) and the largest by +1.95 (glycine cleavage). Fig. 1E has to be read first for this panel to mean what it says.",
    "THE SENTENCE'S \"SIGNIFICANTLY\" DOES NOT COVER EVERY ITEM IT LISTS. Supported at padj < 0.05: translation, pyruvate, serine, dynamics and surveillance, fission, apoptosis, calcium homeostasis, small-molecule transport. NOT supported: OXPHOS subunits (padj 0.27) and every individual complex (CI 0.37, CII 0.54, CIII 0.16, CIV 0.86, CV 0.15) -- no part of the OXPHOS arm is significant on this ruler; glycine cleavage, the largest single promotion (padj 0.12); mitophagy (padj 0.084); and the protein-import tier itself, which is marginal at padj 0.053 although three of its sub-pathways clear it (SAM 0.016, chaperones 0.039, protein homeostasis 0.044). The one OXPHOS row that IS significant is CIII assembly factors, and it is DEMOTED (-0.141, padj 0.022).",
    "MitoCarta sets are membership-loose and several rows are small: glycine cleavage and serine are four genes each, mitophagy 14, fission 15. On a four-gene set a large effect is one or two genes and not a module, which is why n is printed on the panel. Resolve gene by gene before making any mechanistic claim from a small row.",
    "The 144 MitoPathways NEST inside one another -- OXPHOS subunits is contained in OXPHOS, fission and mitophagy in dynamics and surveillance -- so neither the thirteen rows nor the 143 points are independent, and the Benjamini-Hochberg adjustment across them is anticonservative in an unknown direction. Read the pattern.",
    "This is a genotype contrast at each age and both are clean. The 6W-to-12W connector is a comparison BETWEEN two clean contrasts, not a temporal contrast, so it is not exposed to batch = timepoint. The temporal mitoPPS contrasts, which are, are not drawn.",
    "Descriptive and exploratory at n = 6 per group. The compartment-wide claim (a two-sided priority distribution over a one-sided content distribution) is the statement; no single row is a result."),
  source = c(
    "results/background_vs_myc.rds (scripts/40_background_vs_myc_decomposition.R) -- $ruler, $ruler_summary",
    "results/mitopps_scores.rds (scripts/08_mitoPPS_analysis.R) -- $mitopps_pairwise, read for the identity assertion only",
    "Monzel et al. 2025, mitoPPS; reference implementation in external/mitotyping/",
    "Earlier double-column forms: figures/fig02_reallocation_ranked.R (all 144 ranked), figures/figS4_reallocation_full.R (all 144 labelled), figures/figS3_reallocation_two_rulers.R (both rulers side by side)"))

save_panel_p(p, "fig1_reallocation_ranked", height = 74)

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {

  print(p)
  print(LEGEND)

  ## the thirteen rows, both rulers, with padj
  d[, c("pathway", "tier", "n_genes", "p_m6", "p_m6_padj", "p_m12", "c_m6", "c_m12")] |>
    print(row.names = FALSE, digits = 3)

  ## every OXPHOS row -- the sentence says "all complexes" and none is significant
  ruler[ruler$tier == "OXPHOS",
        c("pathway", "n_genes", "p_m6", "p_m6_padj", "p_m12")] |>
    print(row.names = FALSE, digits = 3)

  ## everything that IS significant at 6W, ranked -- the honest version of the list
  s <- ruler[!is.na(ruler$p_m6_padj) & ruler$p_m6_padj < 0.05, ]
  s[order(-s$p_m6), c("pathway", "tier", "n_genes", "p_m6", "p_m6_padj")] |>
    print(row.names = FALSE, digits = 3)

  ## the tier medians the narrative doc quotes, recomputed
  stats::aggregate(cbind(p_m6, p_m12, c_m6) ~ tier, data = ruler[!ruler$is_mtdna, ],
                   FUN = stats::median) |> print(digits = 3)

  ## alternative form kept for comparison: all 144 ranked, only the 7 tiers labelled
  ## source(here::here("figures", "fig02_reallocation_ranked.R"))
}
