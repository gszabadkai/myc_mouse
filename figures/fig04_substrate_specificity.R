# =============================================================================
# fig04_substrate_specificity.R -- what the normal gland does to itself.
# MANUSCRIPT FIGURE 1 (build number 04; see the manifest in paper/myc_mito.qmd).
# -----------------------------------------------------------------------------
# The argument order is SUBSTRATE FIRST: before any statement about Myc, what
# does the wild-type mammary epithelium do between 6W and 12W? Two things, and
# the paper turns on their conjunction.
#
#   A  IT WITHDRAWS FROM RESPIRATION AND FROM NOTHING ELSE THAT WOULD EXPLAIN IT.
#      Set-mean raw log2FC on `timepoint_neg`, each arm against its own
#      expression-matched random-set null. OXPHOS subunits at the extreme; the
#      ASSEMBLY FACTORS OF THE SAME COMPLEXES sit at the 50th percentile, which
#      is the internal control that makes the withdrawal subunit-specific.
#   B  AND THE CLAIM IS COMPARATIVE, SO IT IS TESTED COMPARATIVELY. `PROLIF_*`
#      pooled is itself beyond its own null (1.3rd percentile) -- it is small,
#      not immobile -- so a one-set null cannot carry "de-respires without
#      de-proliferating". The paired null redraws both sets together and tests
#      the DIFFERENCE. That is the actual result.
#   C  ON BOTH RULERS. Content (set-average log2FC) against priority (mitoPPS,
#      content-blind). The respiratory arms fall on both; the growth-coupled arms
#      on neither -- so this is not a normalisation artifact.
#   D  IT DOES NOT DISMANTLE THE DEATH MACHINERY, AND IT DOES NOT BUFFER.
#      Of 32 MitoCarta pro-/anti-apoptotic transcripts plus five non-MitoCarta
#      brakes, exactly one moves at padj<0.05 -- and it moves UP. What the gland
#      withdraws is the STATE the machinery depends on, not the machinery.
#
# BATCH = TIMEPOINT (CLAUDE.md): the 6W and 12W cohorts were extracted as two
# batches, so every wild-type temporal statement in this figure is DESCRIBED, not
# claimed. Panel D's flatness is "no detectable movement at n=6 on a confounded
# axis", not "no movement". Exploratory.
#
# Reads (read-only; the author runs scripts 43 and 44 first):
#   results/substrate_specificity_tradeoff.rds -- $comparator + $wt_null (A),
#       $paired_null (B), $comparator_priority (C), $buffer (D), $defs
#   results/collapse_module_ownership.rds      -- $wt_genes (D)
# =============================================================================

source(here::here("figures", "theme_myc.R"))
if (!requireNamespace("patchwork", quietly = TRUE)) stop("fig04 needs patchwork")
if (!requireNamespace("ggrepel", quietly = TRUE))   stop("fig04 needs ggrepel")

out_dir <- here::here("outputs", "figures")

ss  <- readRDS(here::here("results", "substrate_specificity_tradeoff.rds"))
cmo <- readRDS(here::here("results", "collapse_module_ownership.rds"))

# --- guards: the script-42 failure mode was a silently missing field ----------
need <- function(obj, fields, what) {
  miss <- fields[!fields %in% names(obj)]
  if (length(miss)) stop("fig04: ", what, " is missing -> ", paste(miss, collapse = ", "))
}
need(ss,  c("comparator", "wt_null", "paired_null", "buffer", "defs"),
     "substrate_specificity_tradeoff.rds")
need(cmo, "wt_genes", "collapse_module_ownership.rds")

# --- guard: a stale .rds must fail here, not at review ------------------------
# Script 43 PART A's built-in positive control is that the wild-type OXPHOS-subunit
# value reproduces Issue #4's -0.2548. If membership resolution has drifted, every
# number in this panel is incomparable with the committed ones.
ox_obs <- ss$comparator$c_wt_time[ss$comparator$arm == "OXPHOS subunits"]
ox_ref <- ss$defs$oxphos_wt_reference
if (!isTRUE(abs(ox_obs - ox_ref) < 0.02))
  stop(sprintf(paste("fig04: stale results object -- OXPHOS-subunit wild-type %.4f",
                     "against the reference %.4f. Re-run script 43."), ox_obs, ox_ref))

# --- shared verdict palette (same semantics as fig03 panel D) -----------------
verdict_col <- c("withdraws" = "#762A83", "at chance" = "grey55", "rises" = "#1B7837")
verdict_of  <- function(pct) ifelse(pct < 5, "withdraws",
                             ifelse(pct > 95, "rises", "at chance"))

# =============================================================================
# PANEL A -- the dissociation, each arm against its own matched null
# =============================================================================
A <- merge(ss$comparator[, c("arm", "n_genes", "c_wt_time")],
           ss$wt_null[, c("arm", "percentile", "p_emp_lower", "null_median")],
           by = "arm")
A$verdict <- verdict_of(A$percentile)
A <- A[order(A$c_wt_time), ]
A$label <- sprintf("%s  (%d)", A$arm, A$n_genes)
# level 1 sits at the bottom in a discrete y scale -> reverse so the most
# negative arm is at the TOP and the panel reads as a ranked list.
A$label <- factor(A$label, levels = rev(A$label))
A$pct_lab <- sprintf("%.1f", A$percentile)

xr    <- range(c(A$c_wt_time, 0))
x_txt <- xr[2] + diff(xr) * 0.42

pA <- ggplot2::ggplot(A, ggplot2::aes(y = label)) +
  ggplot2::geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.35) +
  ggplot2::geom_segment(ggplot2::aes(x = 0, xend = c_wt_time, yend = label,
                                     colour = verdict), linewidth = 0.5) +
  ggplot2::geom_point(ggplot2::aes(x = c_wt_time, colour = verdict), size = 2.2) +
  ggplot2::geom_text(ggplot2::aes(x = x_txt, label = pct_lab), hjust = 1,
                     size = 1.95, colour = "grey25") +
  ggplot2::annotate("text", x = x_txt, y = nrow(A) + 0.85, label = "null pct",
                    hjust = 1, size = 1.95, colour = "grey40", fontface = "italic") +
  ggplot2::scale_colour_manual(values = verdict_col, name = NULL) +
  ggplot2::coord_cartesian(xlim = c(xr[1] - 0.03, x_txt + 0.01),
                           ylim = c(0.4, nrow(A) + 1.3), expand = FALSE) +
  ggplot2::labs(
    x = "wild-type 6->12W  (set-average raw log2FC)", y = NULL,
    title = "A   The gland de-respires without de-proliferating",
    subtitle = paste("percentile against 2000 expression-matched random sets",
                     "(every null median is within 0.015 of zero).\nOXPHOS assembly =",
                     "the same complexes as the subunits: the internal control.")) +
  theme_myc(base_size = 8) +
  ggplot2::theme(
    legend.position = "bottom",
    legend.key.size = ggplot2::unit(3, "mm"),
    legend.text     = ggplot2::element_text(size = 6.2),
    axis.text.y     = ggplot2::element_text(size = 6.4),
    plot.title      = ggplot2::element_text(face = "bold", size = 8.5),
    plot.subtitle   = ggplot2::element_text(size = 6.2, colour = "grey25",
                                            lineheight = 1.15))

# =============================================================================
# PANEL B -- the paired null: the claim is a DIFFERENCE, so the null is too
# =============================================================================
B <- ss$paired_null
B$lab <- sprintf("%s\nminus %s", B$arm_a, B$arm_b)
B <- B[order(B$observed_diff), ]
B$lab <- factor(B$lab, levels = rev(B$lab))
B$note <- sprintf("%.2f   %s", B$percentile,
                  ifelse(B$p_emp_lower < 1 / ss$defs$n_set_draws,
                         sprintf("<%.4f", 1 / ss$defs$n_set_draws),
                         sprintf("%.4f", B$p_emp_lower)))

xrB   <- range(c(B$observed_diff, B$null_median, 0))
x_txB <- xrB[2] + diff(xrB) * 0.46

pB <- ggplot2::ggplot(B, ggplot2::aes(y = lab)) +
  ggplot2::geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.35) +
  ggplot2::geom_segment(ggplot2::aes(x = null_median, xend = observed_diff, yend = lab),
                        colour = "grey55", linewidth = 0.4,
                        arrow = ggplot2::arrow(length = ggplot2::unit(1.4, "mm"),
                                               type = "closed")) +
  ggplot2::geom_point(ggplot2::aes(x = null_median), shape = 21, size = 1.8,
                      fill = "white", colour = "grey40", stroke = 0.4) +
  ggplot2::geom_point(ggplot2::aes(x = observed_diff), size = 2.2,
                      colour = verdict_col[["withdraws"]]) +
  ggplot2::geom_text(ggplot2::aes(x = x_txB, label = note), hjust = 1,
                     size = 1.9, colour = "grey25") +
  ggplot2::annotate("text", x = x_txB, y = nrow(B) + 0.72, label = "pct     p",
                    hjust = 1, size = 1.9, colour = "grey40", fontface = "italic") +
  ggplot2::coord_cartesian(xlim = c(xrB[1] - diff(xrB) * 0.06, x_txB + 0.005),
                           ylim = c(0.4, nrow(B) + 1.05), expand = FALSE) +
  ggplot2::labs(
    x = "difference in wild-type 6->12W set-average log2FC", y = NULL,
    title = "B   Tested as a comparison",
    subtitle = "open = paired-null median (both sets redrawn), filled = observed") +
  theme_myc(base_size = 8) +
  ggplot2::theme(axis.text.y   = ggplot2::element_text(size = 6.2, lineheight = 0.95),
                 plot.title    = ggplot2::element_text(face = "bold", size = 8.5),
                 plot.subtitle = ggplot2::element_text(size = 6.2, colour = "grey25",
                                                       lineheight = 1.15))

# =============================================================================
# PANEL C -- both rulers (content is not the whole story, and neither is priority)
# =============================================================================
# comparator_priority is NULL if script 40's output was absent when 43 was run.
if (!is.null(ss$comparator_priority)) {
  Cd <- merge(ss$comparator[, c("arm", "c_wt_time")],
              ss$comparator_priority[, c("arm", "prio_wt_time")], by = "arm")
  Cd <- Cd[is.finite(Cd$prio_wt_time), ]
  Cd$verdict <- verdict_of(ss$wt_null$percentile[match(Cd$arm, ss$wt_null$arm)])
  pC <- ggplot2::ggplot(Cd, ggplot2::aes(c_wt_time, prio_wt_time)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3) +
    ggplot2::geom_point(ggplot2::aes(colour = verdict), size = 2.1) +
    ggrepel::geom_text_repel(ggplot2::aes(label = arm, colour = verdict), size = 2,
                             seed = 1, min.segment.length = 0.2,
                             segment.size = 0.25, box.padding = 0.3,
                             show.legend = FALSE) +
    ggplot2::scale_colour_manual(values = verdict_col, guide = "none") +
    ggplot2::labs(
      x = "content  (set-average log2FC)", y = "priority  (mitoPPS)",
      title = "C   Both rulers agree",
      subtitle = "mitoPPS is content-blind by construction") +
    theme_myc(base_size = 8) +
    ggplot2::theme(plot.title    = ggplot2::element_text(face = "bold", size = 8.5),
                   plot.subtitle = ggplot2::element_text(size = 6.2, colour = "grey25",
                                                         lineheight = 1.15))
} else {
  warning("fig04: comparator_priority is NULL -- panel C omitted (run script 40, then 43)")
  pC <- patchwork::plot_spacer()
}

# =============================================================================
# PANEL D -- the death machinery is untouched, and so is the buffer
# =============================================================================
wg <- cmo$wt_genes
Dd <- data.frame(
  gene  = wg$gene,
  class = ifelse(wg$arm == "PRO", "pro-apoptotic", "anti-apoptotic"),
  lfc   = wg$wt_time,
  padj  = wg$padj_wt,
  stringsAsFactors = FALSE)

# the non-MitoCarta brakes (IAPs and the Bcl2a1 paralog) are not in the apoptosis
# sets, so they are added explicitly -- "does the gland buffer?" is a separate
# question from "does it move its apoptotic transcripts?".
bf <- ss$buffer[!ss$buffer$gene %in% Dd$gene, ]
if (nrow(bf))
  Dd <- rbind(Dd, data.frame(gene = bf$gene, class = "brake (non-MitoCarta)",
                             lfc = bf$lfc_wt_time, padj = bf$padj_wt_time,
                             stringsAsFactors = FALSE))
Dd <- Dd[is.finite(Dd$lfc), ]
# two of the 32 apoptosis-set members have no current MGI symbol in the annotation
# table; they are kept in the distribution and simply go unlabelled.
Dd$gene[is.na(Dd$gene)] <- ""
Dd$sig <- !is.na(Dd$padj) & Dd$padj < 0.05
n_class <- table(Dd$class)
Dd$class <- factor(sprintf("%s  (%d)", Dd$class, n_class[Dd$class]),
                   levels = sprintf("%s  (%d)", c("brake (non-MitoCarta)",
                                                  "anti-apoptotic", "pro-apoptotic"),
                                    n_class[c("brake (non-MitoCarta)",
                                              "anti-apoptotic", "pro-apoptotic")]))
lab_these <- Dd$gene %in% c("Bnip3", "Bbc3", "Bcl2l11", "Bcl2l1", "Bcl2", "Mcl1", "Bax")
jit <- ggplot2::position_jitter(width = 0, height = 0.17, seed = 7)

pD <- ggplot2::ggplot(Dd, ggplot2::aes(lfc, class)) +
  ggplot2::geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.35) +
  ggplot2::geom_point(ggplot2::aes(fill = sig), shape = 21, size = 1.8, stroke = 0.3,
                      colour = "grey35", position = jit) +
  ggrepel::geom_text_repel(data = Dd[lab_these, ],
                           ggplot2::aes(label = gene), size = 2, seed = 7,
                           position = jit, min.segment.length = 0.15,
                           segment.size = 0.25, box.padding = 0.45,
                           point.padding = 0.2, max.overlaps = Inf) +
  ggplot2::scale_fill_manual(values = c(`FALSE` = "grey85", `TRUE` = verdict_col[["rises"]]),
                             guide = "none") +
  ggplot2::labs(
    x = "wild-type 6->12W  (raw log2FC)", y = NULL,
    title = "D   The machinery is left intact",
    subtitle = sprintf(paste("%d of %d transcripts reach padj<0.05 (green)",
                             "--\nand it RISES (Bnip3 %+.2f). No brake moves."),
                       sum(Dd$sig), nrow(Dd), Dd$lfc[which(Dd$gene == "Bnip3")])) +
  theme_myc(base_size = 8) +
  ggplot2::theme(axis.text.y   = ggplot2::element_text(size = 6.4),
                 plot.title    = ggplot2::element_text(face = "bold", size = 8.5),
                 plot.subtitle = ggplot2::element_text(size = 6.2, colour = "grey25",
                                                       lineheight = 1.15))

# =============================================================================
# ASSEMBLY
# =============================================================================
p <- patchwork::wrap_plots(pA, pC, pB, pD, nrow = 2, byrow = TRUE,
                           widths = c(1.45, 1), heights = c(1.35, 1)) +
  patchwork::plot_annotation(
    caption = paste(
      "Wild-type 6->12W only (n=6/timepoint). Content = DESeq2 set-average RAW log2FC on `timepoint_neg`; priority = mitoPPS (Monzel 2025), pairwise-ratio, content-blind.",
      "A and B: 2000 expression-matched random sets drawn within baseMean ventiles. B redraws BOTH sets of a contrast together -- the null a comparative claim needs.",
      "D: the 25 pro- and 7 anti-apoptotic MitoCarta transcripts plus 5 non-MitoCarta brakes (Bcl2 padj 0.45, Bcl2l1 0.79, Mcl1 0.44, Xiap 0.79, Birc2/3/5 n.s.).",
      "BATCH = TIMEPOINT: the two cohorts were extracted as separate batches, so every value here is DESCRIBED, not claimed; D reads as 'no detectable movement at n=6'.",
      sep = "\n"),
    theme = ggplot2::theme(
      plot.caption = ggplot2::element_text(size = 5.8, hjust = 0, colour = "grey30",
                                           lineheight = 1.15)))

if (!isTRUE(getOption("myc.fig.nosave"))) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  save_panel(p, file.path(out_dir, "fig04_substrate_specificity.pdf"),
             width = fig_w[["double"]], height = 150)
}

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {
  ## PART A's ranking, which is the panel: OXPHOS subunits and the TEB signature
  ## at the extreme, the assembly factors of the same complexes at the 50th.
  A[, c("arm", "n_genes", "c_wt_time", "null_median", "percentile")] |> print()

  ## the comparative test -- this is the one that makes A a result
  B[, c("arm_a", "arm_b", "observed_diff", "null_median", "percentile", "p_emp_lower")] |>
    print()

  ## the caveat the prose must carry, and it is NOT in any panel: within the 12
  ## wild-type mice OXPHOS and proliferation are correlated at 0.82, and the
  ## timepoint term survives adjustment for it. Descriptive only -- batch = timepoint.
  ss$wt_within |> print()

  ## the one mover, and the buffer
  Dd[Dd$sig, ] |> print()
  ss$buffer |> print()

  print(pA); print(pB); print(pC); print(pD)
  print(p)
}
