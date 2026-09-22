# =============================================================================
# figureS2_controls.R -- SUPPLEMENTARY FIGURE 2
# "The controls behind the window" -- what Figure 2 rests on.
# -----------------------------------------------------------------------------
# Figure 2 reads PUMA's collapse against a global x0.55 rescaling of the MYC
# programme. Three questions have to be answered before that reading is allowed,
# and a fourth makes the result independent of the ratio it is measured with.
#
#   A  IS THE x0.55 JUST A FALLING DOSE? No. Everything that would register a
#      falling dose WIDENS over the window -- Myc itself, Myc:MAX, Myc relative
#      to the MXD/MNT repressor sum -- while MYC's OUTPUT, the genotype gap in
#      Hallmark MYC targets V2, NARROWS. Driver up, output down: the attenuation
#      is DOWNSTREAM of dose. (Transcript is not protein; the settling
#      measurement is the MYC blot, and it agrees.)
#   B  IS THE x0.55 REAL, AND IS IT UNIFORM? The 12W MYC effect IS the 6W effect
#      at ~55% of the amplitude: same pathways, same order, less of it. Its
#      coherence sits at the 100th percentile of gene-label shuffles matched on
#      set size, expression and pathway overlap. That is what makes "the rest of
#      the programme is merely scaled down" a measurement rather than a metaphor.
#   C  WHAT DOES MYC ADD OVER THE WINDOW, ON TOP OF WHAT DEVELOPMENT DOES?
#      The wild-type gland drifts UP across the window (median +0.04; 103 of 143
#      pathways rise) while the Myc gland drifts DOWN (median -0.13; 22 of 143
#      rise). Per pathway, Myc+ minus wild-type is NEGATIVE in 133 of 143 cases,
#      median -0.18 -- and that difference IS the Myc-specific temporal change, by
#      script 40's exact identity. It is BATCH-CLEAN: both temporal contrasts carry
#      the same batch offset and it cancels in their difference. Ranked, a flat band
#      shows the near-uniformity directly.
#        The inset decomposes the same thing as a regression. TWO PHRASINGS TO
#      AVOID. (i) "A UNIFORM OFFSET" is exact only at slope 1; at 0.75 there is also
#      a component scaling with how far the wild-type moved, so the defensible claim
#      is that the INTERCEPT is the extreme, batch-clean term -- not that the offset
#      is literally constant. (ii) "The Myc change is NOT DUE TO the wild-type
#      change" is too strong the other way: slope 0.75 says much of the movement IS
#      shared. The panel's point is that the shared part is UNREMARKABLE (82nd pct)
#      while an ADDITIONAL downward component is not.
#        WHAT THE NULL IS A NULL OF: a gene-label shuffle preserving set size,
#      expression and pathway overlap -- a test against the data's own internal
#      structure, NOT against sampling new animals. 500 draws, so the intercept's
#      "0th percentile" means p < 1/500 = 0.002 and no smaller. The bootstrap CI on
#      the slope (0.57-0.94) resamples PATHWAYS, not mice.
#   D  DOES PUMA'S COLLAPSE SURVIVE WITHOUT THE BCL-xL DENOMINATOR? Yes. Ranking
#      all Myc-responsive genes by their residual from the global rate -- a
#      statistic that never touches Bcl2l1 -- puts Bbc3 at the 0.5th percentile.
#      Bcl2l11 (BIM), the other half of the western prediction, sits at the 91.6th:
#      the pair splits in vivo, and both halves are reported. Myc's own transcript
#      sits at the top (retention 1.09), which is panel A's premise recovered from
#      inside the RNA.
#
# Reads (read-only; the author runs scripts 40, 41 and 44 first):
#   results/dds_int_run.rds               -- normalised counts (A)
#   results/combined_df_annotated.rds     -- symbol -> ensembl (A)
#   results/background_vs_myc.rds         -- $ruler, $regressions, $regression_null,
#                                            $regression_boot, $defs (B, C)
#   results/collapse_module_ownership.rds -- $collapse_genes, $defs (D)
#   data/genesets_from_library/mammary_mito_myc_metab_v1_mouse.gmt (A)
# =============================================================================

source(here::here("figures", "theme_myc.R"))
source(here::here("functions", "reconcile_gene_symbols.R"))
if (!requireNamespace("patchwork", quietly = TRUE)) stop("figureS2 needs patchwork")
if (!requireNamespace("ggrepel", quietly = TRUE))   stop("figureS2 needs ggrepel")
if (!requireNamespace("DESeq2", quietly = TRUE))    stop("figureS2 needs DESeq2")
if (!requireNamespace("fgsea", quietly = TRUE))     stop("figureS2 needs fgsea (gmtPathways)")

out_dir <- here::here("outputs", "figures")

dds <- readRDS(here::here("results", "dds_int_run.rds"))
cdf <- readRDS(here::here("results", "combined_df_annotated.rds"))
bg  <- readRDS(here::here("results", "background_vs_myc.rds"))
cmo <- readRDS(here::here("results", "collapse_module_ownership.rds"))
gmt <- fgsea::gmtPathways(
  here::here("data", "genesets_from_library", "mammary_mito_myc_metab_v1_mouse.gmt"))

need <- function(obj, fields, what) {
  miss <- fields[!fields %in% names(obj)]
  if (length(miss)) stop("figureS2: ", what, " is missing -> ", paste(miss, collapse = ", "))
}
need(bg,  c("ruler", "regressions", "regression_null", "regression_boot", "defs"),
     "background_vs_myc.rds")
need(cmo, c("collapse_genes", "defs"), "collapse_module_ownership.rds")

# =============================================================================
# PANEL A -- driver against output: the attenuation is downstream of dose
# =============================================================================
nc <- DESeq2::counts(dds, normalized = TRUE)
sm <- as.data.frame(SummarizedExperiment::colData(dds))
sm$group <- factor(sm$group, levels = names(group_labels))

sym2ens <- cdf[!is.na(cdf$mgi_symbol) & !duplicated(cdf$mgi_symbol), c("mgi_symbol", "gene")]
ens_of  <- function(s) sym2ens$gene[match(s, sym2ens$mgi_symbol)]
expr_of <- function(s) {
  e <- ens_of(s)
  if (is.na(e) || !e %in% rownames(nc)) NULL else nc[e, ]
}
if (is.null(expr_of("Myc")))
  stop("figureS2 A: Myc is not in the count matrix -- re-run script 03")

myc      <- expr_of("Myc")
rep_syms <- c("Mxd1", "Mxd3", "Mxd4", "Mxi1", "Mnt", "Mga")
rep_syms <- rep_syms[!vapply(rep_syms, function(g) is.null(expr_of(g)), logical(1))]
if (length(rep_syms) < 4)
  stop("figureS2 A: fewer than 4 MXD/MNT repressors resolved -- check the annotation")
rep_sum  <- Reduce(`+`, lapply(rep_syms, expr_of))

tgt_ens <- recon_to_ensembl(gmt[["MYC_HALLMARK_MYC_TARGETS_V2"]], rownames(nc))
tgt_ens <- tgt_ens[!is.na(tgt_ens)]
tgt_z   <- t(scale(t(log2(nc[tgt_ens, , drop = FALSE] + 1))))
tgt_z   <- tgt_z[is.finite(rowSums(tgt_z)), , drop = FALSE]

measures <- list(`Myc`                 = log2(myc),
                 `Myc : MAX`           = log2(myc / expr_of("Max")),
                 `Myc : repressor sum` = log2(myc / rep_sum),
                 `MYC targets (V2)`    = as.numeric(colMeans(tgt_z)))
kind_of  <- c(`Myc` = "driver", `Myc : MAX` = "driver",
              `Myc : repressor sum` = "driver", `MYC targets (V2)` = "output")
gap_at <- function(y, tp)
  mean(y[sm$group == paste0(tp, "_pos")]) - mean(y[sm$group == paste0(tp, "_neg")])

A <- data.frame(
  measure = factor(names(measures), levels = names(measures)),
  kind    = unname(kind_of[names(measures)]),
  g6      = vapply(measures, gap_at, numeric(1), tp = "6W"),
  g12     = vapply(measures, gap_at, numeric(1), tp = "12W"),
  stringsAsFactors = FALSE)
A$delta <- A$g12 - A$g6

# The panel's whole point, asserted: the three drivers must widen and the output
# must narrow. If that flips, the "downstream of dose" reading is gone.
if (!(all(A$delta[A$kind == "driver"] > 0) && A$delta[A$kind == "output"] < 0))
  warning("figureS2 A: driver/output directions have changed -- re-read the panel text")

kind_col <- c(driver = "#D73027", output = "#4575B4")
pA <- ggplot2::ggplot(A, ggplot2::aes(y = measure)) +
  ggplot2::geom_vline(xintercept = 0, colour = "grey80", linewidth = 0.3) +
  ggplot2::geom_segment(ggplot2::aes(x = g6, xend = g12, yend = measure, colour = kind),
                        linewidth = 0.5,
                        arrow = ggplot2::arrow(length = ggplot2::unit(1.6, "mm"),
                                               type = "closed")) +
  ggplot2::geom_point(ggplot2::aes(x = g6, colour = kind), shape = 21, fill = "white",
                      size = 2, stroke = 0.5) +
  ggplot2::geom_text(ggplot2::aes(x = pmax(g6, g12), label = sprintf("%+.2f", delta),
                                  colour = kind),
                     hjust = -0.3, size = 2.2, fontface = "bold") +
  ggplot2::scale_colour_manual(values = kind_col, name = NULL) +
  ggplot2::scale_y_discrete(limits = rev(levels(A$measure))) +
  ggplot2::coord_cartesian(xlim = c(0, max(c(A$g6, A$g12)) * 1.28)) +
  ggplot2::labs(
    x = "genotype gap (Myc+ - WT), log2", y = NULL,
    title = "A   The driver does not fall while its output does",
    subtitle = "open = 6W, arrow to 12W, label = the change in the gap.\nRatios are read between groups only, so transcript length cancels.") +
  theme_myc(base_size = 8) +
  ggplot2::theme(legend.position = "bottom",
                 legend.key.size = ggplot2::unit(3, "mm"),
                 legend.text     = ggplot2::element_text(size = 6.4),
                 axis.text.y     = ggplot2::element_text(size = 6.8),
                 plot.title      = ggplot2::element_text(face = "bold", size = 8.5),
                 plot.subtitle   = ggplot2::element_text(size = 6.2, colour = "grey25",
                                                         lineheight = 1.15))

# =============================================================================
# PANELS B and C -- the two regressions, on the CONTENT ruler
# -----------------------------------------------------------------------------
# The synthetic mtDNA-encoded pathway is a 4-6.5 SD outlier on every temporal
# contrast and alone moves the shared-vector slope 0.75 -> 0.84, so it is excluded
# from both fits (script 40 flags it as `is_mtdna`).
# =============================================================================
d   <- bg$ruler[!bg$ruler$is_mtdna, ]
reg <- bg$regressions
nul <- bg$regression_null
gv  <- function(model, field)
  reg[[field]][grepl(model, reg$model) & reg$scope == "content, all"]
nv  <- function(stat, field) nul[[field]][nul$statistic == stat]
bv  <- function(model, field) bg$regression_boot[[field]][bg$regression_boot$model == model]

tier_lv  <- bg$defs$tier_levels
tier_lab <- c("Protein import / homeostasis", "Central dogma", "OXPHOS", "Metabolism",
              "Signaling", "Dynamics & surveillance", "SM transport")
names(tier_lab) <- tier_lv
tier_col <- c("Protein import / homeostasis" = "#D55E00", "Central dogma" = "#E69F00",
              "OXPHOS"                  = "#009E73", "Metabolism"    = "grey72",
              "Signaling"               = "#56B4E9",
              "Dynamics & surveillance" = "#0072B2", "SM transport"  = "#CC79A7")
d$Tier <- factor(unname(tier_lab[d$tier]), levels = unname(tier_lab[tier_lv]))
anc <- d[d$pathway %in% tier_lv, ]

scatter_base <- function(dat, xv, yv, xlab, ylab, ttl, sub, note, slope, intercept,
                         legend = "none") {
  rng <- range(c(dat[[xv]], dat[[yv]]), na.rm = TRUE)
  pad <- diff(rng) * 0.05
  ggplot2::ggplot(dat, ggplot2::aes(.data[[xv]], .data[[yv]])) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey88", linewidth = 0.3) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey88", linewidth = 0.3) +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = "grey60",
                         linetype = 2, linewidth = 0.35) +
    ggplot2::geom_point(ggplot2::aes(colour = Tier), size = 1.2, alpha = 0.85) +
    ggplot2::geom_point(data = anc, ggplot2::aes(fill = Tier), shape = 21, size = 2.2,
                        colour = "black", stroke = 0.4, show.legend = FALSE) +
    ggplot2::geom_abline(slope = slope, intercept = intercept, colour = "black",
                         linewidth = 0.5) +
    ggplot2::annotate("text", x = rng[1], y = rng[2], hjust = 0, vjust = 1,
                      size = 1.9, colour = "grey20", lineheight = 1.2, label = note) +
    ggplot2::scale_colour_manual(values = tier_col, name = NULL, drop = FALSE) +
    ggplot2::scale_fill_manual(values = tier_col, guide = "none") +
    ggplot2::coord_cartesian(xlim = rng + c(-pad, pad), ylim = rng + c(-pad, pad)) +
    ggplot2::labs(x = xlab, y = ylab, title = ttl, subtitle = sub) +
    theme_myc(base_size = 8) +
    ggplot2::theme(legend.position = legend,
                   plot.title    = ggplot2::element_text(face = "bold", size = 8.5),
                   plot.subtitle = ggplot2::element_text(size = 6.2, colour = "grey25",
                                                         lineheight = 1.15))
}

pB <- scatter_base(
  d, "c_m6", "c_m12",
  "MYC effect at 6W  (set log2FC)", "MYC effect at 12W  (set log2FC)",
  "B   Rescaled, not reshaped",
  "one point per MitoPathway; dashed = no attenuation",
  sprintf("slope %.2f (boot %.2f-%.2f)\nR2 %.2f, rho %.2f\ncoherence: %.0fth pct of the null\n\nsame pathways, same order,\n%.0f%% of the amplitude",
          gv("^rescale", "slope"), bv("rescale content", "lo"), bv("rescale content", "hi"),
          gv("^rescale", "r2"), gv("^rescale", "rho"), nv("rescale_r2", "percentile"),
          100 * gv("^rescale", "slope")),
  gv("^rescale", "slope"), gv("^rescale", "intercept"))

# --- PANEL C, main: the Myc-specific temporal change, one number per pathway ---
# `c_tp - c_tn` IS the Myc-specific part of the temporal move, by script 40's exact
# identity (change in the genotype gap = Myc+ temporal - WT temporal). It is the
# BATCH-CLEAN quantity for the same reason the regression intercept is: both
# temporal contrasts carry the same batch offset, so it cancels in their difference.
#
# ARTIFACT LEDGER (script 40 PART G) -- DO NOT plot or correlate this difference
# AGAINST either temporal contrast. It IS their difference, so that relationship is
# algebraically forced and carries no information whatsoever.
#
# Ranking it also answers "is the offset uniform?" BY EYE, which the intercept can
# only answer indirectly: a flat band across the ranking is the uniformity claim.
Cd <- d[is.finite(d$c_tn) & is.finite(d$c_tp), ]
Cd$delta <- Cd$c_tp - Cd$c_tn
Cd <- Cd[order(Cd$delta), ]
Cd$rank <- seq_len(nrow(Cd))
ancC    <- Cd[Cd$pathway %in% tier_lv, ]
n_down  <- sum(Cd$delta < 0)
med_d   <- stats::median(Cd$delta)

# The direction is the message, so assert it rather than trusting the caption.
if (!(n_down > 0.8 * nrow(Cd) && med_d < 0))
  warning("figureS2 C: the Myc-specific temporal shift is no longer predominantly ",
          "downward -- re-read the panel text before using it")

pC_main <- ggplot2::ggplot(Cd, ggplot2::aes(rank, delta)) +
  ggplot2::geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.35) +
  ggplot2::geom_hline(yintercept = med_d, colour = "#B2182B", linewidth = 0.4,
                      linetype = 2) +
  ggplot2::geom_point(ggplot2::aes(colour = Tier), size = 1.1, alpha = 0.85) +
  ggplot2::geom_point(data = ancC, ggplot2::aes(fill = Tier), shape = 21, size = 2.1,
                      colour = "black", stroke = 0.4, show.legend = FALSE) +
  ggplot2::annotate("text", x = nrow(Cd) * 0.99, y = med_d, hjust = 1, vjust = 1.6,
                    size = 1.95, colour = "#B2182B",
                    label = sprintf("median %+.2f", med_d)) +
  ggplot2::geom_label(
    data = data.frame(rank = nrow(Cd) * 0.02, delta = max(Cd$delta),
                      lab = sprintf("%d of %d pathways move\nFURTHER DOWN with MYC",
                                    n_down, nrow(Cd))),
    ggplot2::aes(label = lab), inherit.aes = TRUE, hjust = 0, vjust = 1, size = 2.1,
    colour = "grey20", lineheight = 1.2, fill = "white", linewidth = 0,
    label.padding = ggplot2::unit(0.6, "mm")) +
  ggplot2::scale_colour_manual(values = tier_col, guide = "none", drop = FALSE) +
  ggplot2::scale_fill_manual(values = tier_col, guide = "none") +
  ggplot2::labs(
    x = sprintf("%d MitoPathways, ranked", nrow(Cd)),
    y = "Myc+ minus wild-type  (set log2FC)",
    title = "C   What MYC adds over the window",
    subtitle = "The wild-type gland drifts UP over the window and the MYC gland drifts DOWN.\nThe shared batch offset cancels in a difference, which is what makes this clean.") +
  theme_myc(base_size = 8) +
  ggplot2::theme(plot.title    = ggplot2::element_text(face = "bold", size = 8.5),
                 plot.subtitle = ggplot2::element_text(size = 6.2, colour = "grey25",
                                                       lineheight = 1.15))

# --- PANEL C, inset: the same thing as a regression ---------------------------
# Kept because it decomposes the difference into a SHARED slope and a MYC-SPECIFIC
# intercept, which the one-dimensional view cannot. Orientation is deliberate and
# must not be flipped: regression is NOT symmetric (c_tp ~ c_tn gives 0.749/-0.172,
# c_tn ~ c_tp gives 0.559/+0.118, and the geometric mirror would be 1.335), and the
# wild-type gland is the substrate, hence the predictor.
rngC <- range(c(Cd$c_tn, Cd$c_tp), na.rm = TRUE)
pC_inset <- ggplot2::ggplot(Cd, ggplot2::aes(c_tn, c_tp)) +
  ggplot2::geom_abline(slope = 1, intercept = 0, colour = "grey65", linetype = 2,
                       linewidth = 0.3) +
  ggplot2::geom_point(colour = "grey45", size = 0.45, alpha = 0.7) +
  ggplot2::geom_abline(slope = gv("^shared", "slope"),
                       intercept = gv("^shared", "intercept"),
                       colour = "black", linewidth = 0.45) +
  ggplot2::coord_cartesian(xlim = rngC, ylim = rngC) +
  ggplot2::labs(x = "WT", y = "Myc+") +
  theme_myc(base_size = 6) +
  ggplot2::theme(
    axis.title = ggplot2::element_text(size = 5.2),
    axis.text  = ggplot2::element_text(size = 4.2),
    axis.ticks = ggplot2::element_line(linewidth = 0.2),
    axis.line  = ggplot2::element_line(linewidth = 0.2),
    plot.background  = ggplot2::element_rect(fill = "white", colour = "grey75",
                                             linewidth = 0.25),
    plot.margin = ggplot2::margin(1, 2, 1, 1))

# The inset sits BOTTOM-RIGHT. Ranked ascending, the most negative values are at
# LOW rank, so the high-rank / low-value corner is the empty one -- top-left is NOT
# empty, because the flat band sits high in the y range.
# annotation_custom rather than patchwork::inset_element: this figure nests
# wrap_plots, and inset_element interacts with nesting in ways annotation_custom
# does not.
inset_note <- sprintf(
  "the same, decomposed:\nslope %.2f = the SHARED move\n  null %.2f, %.0fth pct -> NOT beyond it\nintercept %.2f = MYC-SPECIFIC\n  beyond all %d shuffles, p<%.3f",
  gv("^shared", "slope"), nv("shared_slope", "null_median"),
  nv("shared_slope", "percentile"), gv("^shared", "intercept"),
  bg$defs$n_perm, 1 / bg$defs$n_perm)

yr <- range(Cd$delta)
pC <- pC_main +
  ggplot2::annotation_custom(
    ggplot2::ggplotGrob(pC_inset),
    xmin = nrow(Cd) * 0.66, xmax = nrow(Cd) * 1.00,
    ymin = yr[1] - diff(yr) * 0.02, ymax = yr[1] + diff(yr) * 0.43) +
  ggplot2::annotate("text", x = nrow(Cd) * 0.63, y = yr[1] + diff(yr) * 0.21,
                    hjust = 1, vjust = 0.5, size = 1.8, colour = "grey25",
                    lineheight = 1.3, label = inset_note)

# =============================================================================
# PANEL D -- the collapse, without the Bcl-xL denominator
# =============================================================================
cg <- as.data.frame(cmo$collapse_genes)
MARK <- c("Bbc3", "Bcl2l11", "Myc")
if (!all(MARK %in% cg$gene))
  stop("figureS2 D: collapse_genes has lost Bbc3 / Bcl2l11 / Myc -- re-run script 44")
mk <- cg[match(MARK, cg$gene), c("gene", "z_resid", "pct_z", "retention")]
mk$lab <- sprintf("%s\n%.2f pct", mk$gene, mk$pct_z)

zr    <- c(min(cg$z_resid) - 0.2, 6.2)
n_off <- sum(cg$z_resid > zr[2])
top   <- max(graphics::hist(cg$z_resid, breaks = 80, plot = FALSE)$counts)

pD <- ggplot2::ggplot(cg, ggplot2::aes(z_resid)) +
  ggplot2::geom_histogram(bins = 80, fill = "grey85", colour = NA) +
  ggplot2::geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.3) +
  ggplot2::geom_segment(data = mk,
                        ggplot2::aes(x = z_resid, xend = z_resid, y = 0, yend = top * 0.6),
                        colour = "#D73027", linewidth = 0.45) +
  ggrepel::geom_text_repel(data = mk,
                           ggplot2::aes(x = z_resid, y = top * 0.6, label = lab),
                           size = 2, seed = 3, colour = "#D73027", direction = "y",
                           nudge_y = top * 0.22, segment.size = 0.25,
                           min.segment.length = 0, box.padding = 0.3,
                           lineheight = 0.95) +
  ggplot2::coord_cartesian(xlim = zr, expand = FALSE) +
  ggplot2::labs(
    x = "collapse statistic  (residual from the global rate)", y = "genes",
    title = "D   PUMA collapses without BCL-xL",
    subtitle = sprintf("%s MYC-responsive genes. Bbc3 and Bcl2l11 were\nPRE-SPECIFIED from the PGC-1a westerns, and they SPLIT.",
                       format(nrow(cg), big.mark = ","))) +
  theme_myc(base_size = 8) +
  ggplot2::theme(plot.title    = ggplot2::element_text(face = "bold", size = 8.5),
                 plot.subtitle = ggplot2::element_text(size = 6.2, colour = "grey25",
                                                       lineheight = 1.15))

# =============================================================================
# ASSEMBLY
# =============================================================================
p <- patchwork::wrap_plots(pA, pB, pC, pD, nrow = 2, byrow = TRUE,
                           widths = c(1, 1), heights = c(1, 1)) +
  patchwork::plot_annotation(
    caption = paste(
      sprintf("A: normalised counts, n=6 per group; ratios read between groups only. Transcript is not protein -- the settling measurement is the MYC blot, and it agrees (protein -50%%, i.e. the x%.2f rate).",
              gv("^rescale", "slope")),
      "B/C: one point per MitoPathway on the content ruler (set-average RAW log2FC), coloured by MitoPathway tier as in Figure 1B; ringed points are the seven top-level arms. The synthetic",
      "   mtDNA-encoded pathway is excluded as a 4-6.5 SD outlier on every temporal contrast.",
      "   Nulls are gene-label shuffles preserving set size, expression AND pathway overlap -- MitoPathways nest, so overlap alone would manufacture cross-pathway correlation. They test against the",
      sprintf("   data's own structure, NOT against sampling new animals; %d draws, so a value beyond all of them means p < %.3f and no smaller. The bootstrap CI resamples PATHWAYS, not mice.",
              bg$defs$n_perm, 1 / bg$defs$n_perm),
      "C: Myc+ minus wild-type per pathway IS the Myc-specific temporal change (script 40's identity), and it is batch-clean because the shared offset cancels in a difference. It must never be plotted",
      "   AGAINST either temporal contrast -- it is their difference, so that relationship is algebraically forced (script 40's artifact ledger). Inset: the same decomposed; orientation is deliberate,",
      "   since regression is not symmetric (0.75/-0.17 one way, 0.56/+0.12 the other) and the wild-type gland is the predictor.",
      sprintf("D: rank on sign(LFC_6W) x (LFC_12W - rate x LFC_6W), standardised by the two contrasts' independent SEs; fitted rate %.2f, axis clipped at %.0f (%d gene beyond).",
              cmo$defs$global_rate_fitted, zr[2], n_off),
      "BATCH = TIMEPOINT: it cancels in panel C's difference (and in the inset's intercept), which is why those are the quantities read, but not in the temporal contrasts themselves.",
      sep = "\n"),
    theme = ggplot2::theme(
      plot.caption = ggplot2::element_text(size = 5.6, hjust = 0, colour = "grey30",
                                           lineheight = 1.15)))

if (!isTRUE(getOption("myc.fig.nosave"))) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  save_panel(p, file.path(out_dir, "figureS2_controls.pdf"),
             width = fig_w[["double"]], height = 162)
}

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {
  ## A: drivers must widen, the output must narrow
  A[, c("measure", "kind", "g6", "g12", "delta")] |> print()

  ## B/C: the two regressions and their nulls
  reg[reg$scope == "content, all", c("model", "slope", "intercept", "r2", "rho")] |> print()
  nul |> print()

  ## D: the three marked genes, and the fitted rate against the assumed one
  print(mk)
  cmo$defs[c("global_rate_assumed", "global_rate_fitted", "n_reported_genes")] |> str()

  print(pA); print(pB); print(pC); print(pD)
  print(p)
}
