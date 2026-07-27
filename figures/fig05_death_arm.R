# =============================================================================
# fig05_death_arm.R -- Myc builds a death-ready mitochondrion, and loses the trigger.
# MANUSCRIPT FIGURE 5 (build number 05; see the manifest in paper/myc_mito.qmd).
# -----------------------------------------------------------------------------
# On the substrate of fig04, what does Myc do -- and what does it stop doing?
#
#   A  IT BUILDS THE ORGANELLE AND THE MEANS OF ITS OWN EXECUTION, TOGETHER.
#      Genotype effect at 6W across the hand-curated roster: import and cristae up,
#      OXPHOS subunits up, HTRA2 up, BAX up, BCL-xL down. Asset and liability
#      constructed in one move -- the trade-off, stated as molecules.
#   B  WHAT FADES BY 12W FADES AT THE GLOBAL RATE -- EXCEPT THE TRIGGER. Retention
#      (12W effect / 6W effect) against the x0.55 global rescaling that scripts 40
#      and 41 established. Nearly every gene sits on the line. `Bbc3` reverses sign.
#      OPEN POINTS ARE GENES WHOSE 6W EFFECT IS NOT SIGNIFICANT: their retention is
#      a ratio of noise and is shown, not leaned on.
#   C  THE SAME THING AT THE PAIR LEVEL, WITH THE NULL THAT MATTERS. Each pro:anti
#      ratio against a matched-pair null CONDITIONED on a 6W effect at least as
#      large -- i.e. "among pairs that start this high, how many retain less?".
#      The null's own median retention reconstructs x0.55 from random gene pairs,
#      which is the built-in positive control. Read only the pairs that HAVE a 6W
#      effect (filled): among those, `Bbc3:Bcl2l1` is the one that collapses.
#   D  AND IT REPRODUCES GENOME-WIDE ON A STATISTIC THAT NEVER TOUCHES Bcl2l1.
#      Residual from the global rate over 8774 Myc-responsive genes. `Bbc3` at the
#      0.5th percentile; `Bcl2l11` (BIM) at the 91.6th -- the westerns' PUMA+BIM
#      pair does NOT transfer as a pair, and both halves are reported. `Myc` itself
#      at the 100th (retention 1.09): the transgene message does not attenuate, so
#      the downstream fade is not transcript-level dose. That control was unplanned.
#
# EPISTEMIC CONTRACT (docs/2026-07-25_death_narrative_dose_vs_competence.md sec 0):
# the in-vivo transcriptome GENERATES the hypothesis; the cell perturbations PROVE
# it. `Bbc3` and `Bcl2l11` were PRE-SPECIFIED from the PGC1a westerns, not found by
# scanning -- which is why they are reported despite failing a padj gate at 6W.
# Nothing here survives BH across the pair panel (min int_p_bh 0.33); this is
# ranking plus pre-specification, not a confirmatory test at n=24.
#
# Reads (read-only; the author runs scripts 42 and 44 first):
#   results/priming_arm_teb.rds            -- $machinery (A, B), $priming +
#       $pair_null (C), $params$GLOBAL_RATE
#   results/collapse_module_ownership.rds  -- $collapse_genes + $defs (D)
# =============================================================================

source(here::here("figures", "theme_myc.R"))
if (!requireNamespace("patchwork", quietly = TRUE)) stop("fig05 needs patchwork")
if (!requireNamespace("ggrepel", quietly = TRUE))   stop("fig05 needs ggrepel")

out_dir <- here::here("outputs", "figures")

pa  <- readRDS(here::here("results", "priming_arm_teb.rds"))
cmo <- readRDS(here::here("results", "collapse_module_ownership.rds"))

need <- function(obj, fields, what) {
  miss <- fields[!fields %in% names(obj)]
  if (length(miss)) stop("fig05: ", what, " is missing -> ", paste(miss, collapse = ", "))
}
need(pa,  c("machinery", "priming", "pair_null", "params"), "priming_arm_teb.rds")
need(cmo, c("collapse_genes", "defs"), "collapse_module_ownership.rds")

RATE <- pa$params$GLOBAL_RATE
stopifnot(is.finite(RATE), RATE > 0, RATE < 1)

# --- guard: the pre-specified genes must still be locatable -------------------
cg <- as.data.frame(cmo$collapse_genes)
if (!all(c("Bbc3", "Bcl2l11", "Myc") %in% cg$gene))
  stop("fig05: collapse_genes has lost Bbc3 / Bcl2l11 / Myc -- re-run script 44")

# --- roster arms collapsed to the six functional blocks the panels read -------
block_of <- function(arm) {
  ifelse(grepl("^execution", arm),               "execution",
  ifelse(grepl("^effector", arm),                "effector (Bax/Bak)",
  ifelse(grepl("^BH3-only", arm),                "BH3-only trigger",
  ifelse(grepl("^brake", arm),                   "brake",
  ifelse(grepl("^OXPHOS subunit", arm),          "OXPHOS subunit",
  ifelse(grepl("^biogenesis TF|coactivator|mtDNA machinery", arm), "biogenesis TF",
                                                 "import / cristae"))))))
}
BLOCKS <- c("biogenesis TF", "import / cristae", "OXPHOS subunit",
            "brake", "BH3-only trigger", "effector (Bax/Bak)", "execution")

M <- as.data.frame(pa$machinery)
M$block <- factor(block_of(M$arm), levels = BLOCKS)
M$sig6  <- !is.na(M$padj_myc_6W) & M$padj_myc_6W < 0.05

sig_col <- c(`TRUE` = "#D73027", `FALSE` = "grey80")   # project Myc+ red
jit <- ggplot2::position_jitter(width = 0, height = 0.16, seed = 11)

# =============================================================================
# PANEL A -- the 6W genotype effect across the roster
# =============================================================================
labA <- M$gene %in% c("Htra2", "Bax", "Bcl2l1", "Bbc3", "Hspd1", "Tomm22",
                      "Atp5f1a", "Esrra", "Birc5", "mt-Co1")

pA <- ggplot2::ggplot(M, ggplot2::aes(lfc_myc_6W, block)) +
  ggplot2::geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.35) +
  ggplot2::geom_point(ggplot2::aes(fill = sig6), shape = 21, size = 1.9, stroke = 0.3,
                      colour = "grey35", position = jit) +
  ggrepel::geom_text_repel(data = M[labA, ], ggplot2::aes(label = gene), size = 2,
                           seed = 11, position = jit, min.segment.length = 0.15,
                           segment.size = 0.25, box.padding = 0.42,
                           point.padding = 0.2, max.overlaps = Inf) +
  ggplot2::scale_fill_manual(values = sig_col, guide = "none") +
  ggplot2::labs(
    x = "Myc effect at 6W  (raw log2FC)", y = NULL,
    title = "A   Myc builds a death-ready mitochondrion",
    subtitle = "filled = padj<0.05. Import, cristae and OXPHOS up;\nHTRA2 and BAX up, BCL-xL down -- asset and liability together.") +
  theme_myc(base_size = 8) +
  ggplot2::theme(axis.text.y   = ggplot2::element_text(size = 6.4),
                 plot.title    = ggplot2::element_text(face = "bold", size = 8.5),
                 plot.subtitle = ggplot2::element_text(size = 6.2, colour = "grey25",
                                                       lineheight = 1.15))

# =============================================================================
# PANEL B -- retention against the x0.55 global rate
# =============================================================================
# `retention` is NA where |6W effect| < 0.2: script 42 refuses to divide by a
# non-effect, and so does this panel.
R <- M[is.finite(M$retention), ]
labB <- R$gene %in% c("Bbc3", "Bcl2l11", "Bax", "Bcl2l1", "Htra2", "Bak1", "Birc5")

pB <- ggplot2::ggplot(R, ggplot2::aes(retention, block)) +
  ggplot2::geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3) +
  ggplot2::geom_vline(xintercept = RATE, colour = "grey30", linewidth = 0.4,
                      linetype = 2) +
  ggplot2::geom_point(ggplot2::aes(fill = sig6), shape = 21, size = 1.9, stroke = 0.3,
                      colour = "grey35", position = jit) +
  ggrepel::geom_text_repel(data = R[labB, ], ggplot2::aes(label = gene), size = 2,
                           seed = 11, position = jit, min.segment.length = 0.15,
                           segment.size = 0.25, box.padding = 0.42,
                           point.padding = 0.2, max.overlaps = Inf) +
  ggplot2::annotate("text", x = RATE, y = 0.52,
                    label = sprintf(" global rate x%.2f", RATE), size = 1.95,
                    colour = "grey30", hjust = 0) +
  ggplot2::scale_fill_manual(values = sig_col, guide = "none") +
  ggplot2::coord_cartesian(ylim = c(0.35, length(BLOCKS) + 0.6)) +
  ggplot2::labs(
    x = "retention  (12W Myc effect / 6W Myc effect)", y = NULL,
    title = "B   Everything fades at that rate -- except PUMA",
    subtitle = "open = 6W effect not significant, so its retention is a ratio\nof noise: shown, not leaned on. Bbc3 reverses sign.") +
  theme_myc(base_size = 8) +
  ggplot2::theme(axis.text.y   = ggplot2::element_text(size = 6.4),
                 plot.title    = ggplot2::element_text(face = "bold", size = 8.5),
                 plot.subtitle = ggplot2::element_text(size = 6.2, colour = "grey25",
                                                       lineheight = 1.15))

# =============================================================================
# PANEL C -- the pro:anti ratios against the conditional matched-pair null
# =============================================================================
P <- merge(as.data.frame(pa$pair_null),
           as.data.frame(pa$priming)[, c("pair", "d6", "p6", "int_p", "int_p_adj")],
           by = "pair")
P$real6 <- P$p6 < 0.05
P <- P[order(P$retention), ]
P$pair <- factor(P$pair, levels = P$pair)
P$note <- sprintf("%.0f   %.3f", P$pct_retention_cond, P$p_emp_cond)

xrC   <- range(c(P$retention, P$null_retention_median, 0))
x_txC <- xrC[2] + diff(xrC) * 0.40

pC <- ggplot2::ggplot(P, ggplot2::aes(y = pair)) +
  ggplot2::geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3) +
  ggplot2::geom_vline(xintercept = RATE, colour = "grey30", linewidth = 0.4,
                      linetype = 2) +
  ggplot2::geom_segment(ggplot2::aes(x = null_retention_median, xend = retention,
                                     yend = pair), colour = "grey60", linewidth = 0.35) +
  ggplot2::geom_point(ggplot2::aes(x = null_retention_median), shape = 124, size = 2,
                      colour = "grey35") +
  ggplot2::geom_point(ggplot2::aes(x = retention, fill = real6), shape = 21, size = 2.2,
                      stroke = 0.35, colour = "grey25") +
  ggplot2::geom_text(ggplot2::aes(x = x_txC, label = note), hjust = 1, size = 1.9,
                     colour = "grey25") +
  ggplot2::annotate("text", x = x_txC, y = nrow(P) + 0.85, label = "pct      p",
                    hjust = 1, size = 1.9, colour = "grey40", fontface = "italic") +
  ggplot2::scale_fill_manual(values = sig_col, guide = "none") +
  ggplot2::coord_cartesian(xlim = c(xrC[1] - diff(xrC) * 0.06, x_txC + 0.02),
                           ylim = c(0.4, nrow(P) + 1.2), expand = FALSE) +
  ggplot2::labs(
    x = "retention of the pro:anti ratio", y = NULL,
    title = "C   The priming ratios, against a conditional null",
    subtitle = sprintf(paste("tick = median retention of matched pairs conditioned on the",
                             "same 6W effect\n(it reconstructs x%.2f from random pairs --",
                             "the positive control).\nFilled = the pair HAS a 6W effect",
                             "(p6<0.05); open pairs are unreadable."),
                       stats::median(P$null_retention_median))) +
  theme_myc(base_size = 8) +
  ggplot2::theme(axis.text.y   = ggplot2::element_text(size = 6.2),
                 plot.title    = ggplot2::element_text(face = "bold", size = 8.5),
                 plot.subtitle = ggplot2::element_text(size = 6.2, colour = "grey25",
                                                       lineheight = 1.15))

# =============================================================================
# PANEL D -- the genome-wide reproduction, independent of the Bcl2l1 denominator
# =============================================================================
MARK <- c("Bbc3", "Bcl2l11", "Myc")
mk <- cg[match(MARK, cg$gene), c("gene", "z_resid", "pct_z", "retention")]
mk$lab <- sprintf("%s\n%.2f pct", mk$gene, mk$pct_z)

# one gene (Firrm, z 9.66) sits far beyond the rest; clipping the axis keeps the
# marked genes readable and is stated in the caption.
zr  <- c(min(cg$z_resid) - 0.2, 6.2)
n_off <- sum(cg$z_resid > zr[2])
top <- max(graphics::hist(cg$z_resid, breaks = 80, plot = FALSE)$counts)

pD <- ggplot2::ggplot(cg, ggplot2::aes(z_resid)) +
  ggplot2::geom_histogram(bins = 80, fill = "grey85", colour = NA) +
  ggplot2::geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.3) +
  ggplot2::geom_segment(data = mk,
                        ggplot2::aes(x = z_resid, xend = z_resid, y = 0,
                                     yend = top * 0.62),
                        colour = "#D73027", linewidth = 0.45) +
  ggrepel::geom_text_repel(data = mk,
                           ggplot2::aes(x = z_resid, y = top * 0.62, label = lab),
                           size = 2, seed = 3, colour = "#D73027",
                           direction = "y", nudge_y = top * 0.2,
                           segment.size = 0.25, min.segment.length = 0,
                           box.padding = 0.3, lineheight = 0.95) +
  ggplot2::coord_cartesian(xlim = zr, expand = FALSE) +
  ggplot2::labs(
    x = "collapse statistic  (residual from the global rate)",
    y = "genes",
    title = "D   Reproduced genome-wide, without Bcl-xL",
    subtitle = sprintf(paste("%s Myc-responsive genes. Bbc3 and Bcl2l11 were",
                             "PRE-SPECIFIED\nfrom the PGC1a westerns -- and they split."),
                       format(nrow(cg), big.mark = ","))) +
  theme_myc(base_size = 8) +
  ggplot2::theme(plot.title    = ggplot2::element_text(face = "bold", size = 8.5),
                 plot.subtitle = ggplot2::element_text(size = 6.2, colour = "grey25",
                                                       lineheight = 1.15))

# =============================================================================
# ASSEMBLY
# =============================================================================
p <- patchwork::wrap_plots(pA, pB, pC, pD, nrow = 2, byrow = TRUE,
                           widths = c(1, 1), heights = c(1, 1.05)) +
  patchwork::plot_annotation(
    caption = paste(
      "Raw (unshrunken) DESeq2 log2FC throughout; n=6/group. A/B: the curated death-and-mitochondrion roster of script 42 PART A, grouped by function.",
      "B/C: retention is read against the x0.55 global rescaling of the whole Myc programme (scripts 40, 41) -- the null here is 'fades like everything else', not 'no change'.",
      "C: 4000 baseMean-matched random (pro-like, anti-like) pairs, conditioned on a 6W effect at least as large as the observed one. No pair survives BH (min adjusted interaction p 0.33).",
      sprintf(paste("D: rank on sign(LFC_6W) x (LFC_12W - rate x LFC_6W), standardised by the two contrasts' independent SEs;",
                    "fitted rate %.2f, axis clipped at %.0f (%d gene beyond)."),
              cmo$defs$global_rate_fitted, zr[2], n_off),
      "Ranking plus pre-specification, not a confirmatory test at n=24.",
      sep = "\n"),
    theme = ggplot2::theme(
      plot.caption = ggplot2::element_text(size = 5.8, hjust = 0, colour = "grey30",
                                           lineheight = 1.15)))

if (!isTRUE(getOption("myc.fig.nosave"))) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  save_panel(p, file.path(out_dir, "fig05_death_arm.pdf"),
             width = fig_w[["double"]], height = 158)
}

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {
  ## panel B in table form: what retains at the global rate and what does not
  R[order(R$retention), c("gene", "block", "lfc_myc_6W", "padj_myc_6W",
                          "lfc_myc_12W", "retention", "vs_global")] |> print()

  ## panel C, and the honesty check: Bmf sits further out than Bbc3 on the
  ## conditional null (97.5th, p 0.025) but has NO 6W effect (d6 0.24, p6 0.40),
  ## so it fails the precondition the whole statistic rests on. Filled vs open.
  P[, c("pair", "d6", "p6", "retention", "null_retention_median",
        "pct_retention_cond", "p_emp_cond", "int_p", "int_p_adj")] |> print()

  ## panel D: the three marked genes, and the fitted rate next to the assumed one
  print(mk)
  cmo$defs[c("global_rate_assumed", "global_rate_fitted", "rank_agreement_spearman",
             "n_reported_genes", "n_ranking_genes")] |> str()

  ## the single-gene null on the same three, from script 42 PART B
  pa$single_gene_null |> print()

  print(pA); print(pB); print(pC); print(pD)
  print(p)
}
