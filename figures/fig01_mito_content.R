# =============================================================================
# fig01_mito_content.R -- "Myc raises mitochondrial content"
# -----------------------------------------------------------------------------
# The single strongest new result of Block B, taken to publication quality.
# Genotype (Myc+ vs WT) WITHIN a timepoint is the CLEAN axis (depth-balanced,
# litter-controlled). The 6W->12W time axis is cohort/batch-confounded
# (batch = timepoint) and is drawn only as context, in grey.
#
# Reads (read-only, already on disk):
#   results/mito_content_proxies.rds  (script 32) -- $shares (per-sample % of the
#       non-mtDNA transcriptome, n=6/group; denominator drops ONLY the 13 mt-*
#       genes, so nuclear MitoCarta stays IN it), $share_stats (genotype main
#       effect per arm). Bracket p-values are recomputed here by script 32's OWN
#       method (log2-share simple-effect lm within a subset) so they reconcile.
#
# Five claim-bearing arms: nuclear arms rise (mass markers, nuclear MitoCarta,
# nuclear OXPHOS, mito biogenesis) while mtDNA-encoded output stays flat -- content
# AND where the imbalance is NOT, in one panel. "Mito biogenesis" = BIOGENESIS_FULL
# (script 32 roster) = MitoCarta central dogma + protein import/sorting/homeostasis,
# 314 genes of which only 7 are in nuclear OXPHOS -- so OXPHOS vs general biogenesis
# is a genuinely disjoint comparison (the old mitoribosome facet was a sub-arm).
#
# Design (skills: distribution-plots, ggplot2-fundamentals, figure-export):
#   n=6/group -> every point (quasirandom) over a light box; never bar-of-mean.
#   x-axis dropped; the four groups carry a colour legend (hue=genotype,
#   lightness=age). p-values sit on brackets between conditions.
# =============================================================================

source(here::here("figures", "theme_myc.R"))

out_dir <- here::here("outputs", "figures")

content <- readRDS(here::here("results", "mito_content_proxies.rds"))

# --- the five claim-bearing arms, ordered nuclear-up first, mtDNA (flat) last -
arms <- c("MASS_MARKERS_NOCHAP",
          "MITOCARTA_NUCLEAR_ENCODED",
          "MITOCARTA_OXPHOS_NU",
          "BIOGENESIS_FULL",
          "MITOCARTA_MTDNA_ENCODED")
arm_name <- c(MASS_MARKERS_NOCHAP       = "Mass markers",
              MITOCARTA_NUCLEAR_ENCODED = "Nuclear MitoCarta",
              MITOCARTA_OXPHOS_NU       = "Nuclear OXPHOS",
              BIOGENESIS_FULL           = "Mito biogenesis",
              MITOCARTA_MTDNA_ENCODED   = "mtDNA-encoded")

# --- per-sample data for the five arms (title carries the claim; strips carry
#     only the arm name -- keep text light) -------------------------------------
df <- content$shares
df <- df[df$panel %in% arms, ]
df$panel <- factor(df$panel, levels = arms, labels = arm_name[arms])
df$group <- factor(df$group, levels = names(group_labels))

n_per <- min(table(df$group[df$panel == arm_name[[arms[1]]]]))

# --- off-scale handling: one 12W WT mtDNA sample sits at ~67%; every other point
#     is <40. Cap the DISPLAY at 40 and mark the capped point with its true value.
CAP <- 40
df$y_disp <- pmin(df$share_nomt, CAP)
df$capped <- df$share_nomt > CAP

# --- bracket p-values, by script 32's method (simple-effect lm on log2 share) --
#   Four comparisons the author specified:
#     C1 genotype @ 6W   (WT vs Myc+, x1-x2)  -- CLEAN  (black)
#     C2 time in WT      (6W vs 12W, x1-x3)    -- CONFOUNDED (grey; batch=timepoint)
#     C3 genotype @ 12W  (WT vs Myc+, x3-x4)   -- CLEAN  (black)
#     C4 time in Myc+    (6W vs 12W, x2-x4)    -- CONFOUNDED (grey)
#   x positions follow the group factor: 6W_neg=1, 6W_pos=2, 12W_neg=3, 12W_pos=4.
comp <- data.frame(
  cid   = c("C1", "C2", "C3", "C4"),
  x1    = c(1, 1, 3, 2),
  x2    = c(2, 3, 4, 4),
  kind  = c("geno", "time", "geno", "time"),
  key   = c("6W", "neg", "12W", "pos"),   # geno -> timepoint subset; time -> myc subset
  level = c(1, 3, 1, 2),                  # vertical stacking (wider brackets higher)
  clean = c(TRUE, FALSE, TRUE, FALSE),
  stringsAsFactors = FALSE)

pval_for <- function(d, kind, key) {
  d$yv <- log2(d$share_nomt)
  m <- if (kind == "geno") {
    stats::lm(yv ~ myc_status, data = d[d$timepoint == key, ])
  } else {
    stats::lm(yv ~ timepoint, data = d[d$myc_status == key, ])
  }
  summary(m)$coefficients[2, "Pr(>|t|)"]
}

fmt_p <- function(p) if (p < 0.001) "<0.001" else formatC(p, format = "g", digits = 2)

brk <- do.call(rbind, lapply(arms, function(a) {
  d      <- content$shares[content$shares$panel == a, ]
  facmax <- max(pmin(d$share_nomt, CAP))
  b <- comp
  b$panel <- factor(arm_name[[a]], levels = arm_name[arms])
  b$p     <- vapply(seq_len(nrow(b)),
                    function(i) pval_for(d, b$kind[i], b$key[i]), numeric(1))
  b$lab   <- vapply(b$p, fmt_p, character(1))
  b$y     <- facmax * (1 + 0.075 * b$level)  # tight stacking (levels 1..3)
  b$tick  <- facmax * 0.012                  # short end-ticks
  b$xmid  <- (b$x1 + b$x2) / 2
  b
}))
# significant p-values in red so they stand out; the rest neutral grey.
# (Batch-confounding of the 6W-vs-12W brackets is handled in the text, not here.)
brk$col <- ifelse(brk$p < 0.05, "sig", "ns")

# per-facet headroom so the top bracket + its label fit while the data keep
# ~2/3 of the axis (brackets compressed into the top strip).
hr <- do.call(rbind, lapply(arms, function(a) {
  d <- content$shares[content$shares$panel == a, ]
  data.frame(panel = factor(arm_name[[a]], levels = arm_name[arms]),
             group = factor("6W_neg", levels = names(group_labels)),
             y     = max(pmin(d$share_nomt, CAP)) * 1.27)
}))

# --- point geom: quasirandom if available, else jitter (n=6 -> show all) -------
pts_layer <- function(dat) {
  if (requireNamespace("ggbeeswarm", quietly = TRUE)) {
    ggbeeswarm::geom_quasirandom(data = dat, width = 0.22, size = 1.2, alpha = 0.95)
  } else {
    ggplot2::geom_jitter(data = dat, width = 0.16, height = 0, size = 1.2, alpha = 0.95)
  }
}

p <- ggplot2::ggplot(df, ggplot2::aes(group, y_disp, colour = group, fill = group)) +
  ggplot2::geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.28,
                        colour = "grey35", linewidth = 0.3) +
  pts_layer(df[!df$capped, ]) +
  # capped (off-scale) points: triangle at the cap + true value
  ggplot2::geom_point(data = df[df$capped, ], shape = 17, size = 1.6) +
  ggplot2::geom_text(data = df[df$capped, ],
                     ggplot2::aes(label = sprintf("%.0f", share_nomt)),
                     vjust = -0.7, size = 2.1, colour = "grey30",
                     show.legend = FALSE) +
  # brackets: horizontal bar + two end-ticks + p label (red if p<0.05)
  ggplot2::geom_segment(data = brk, inherit.aes = FALSE,
                        ggplot2::aes(x = x1, xend = x2, y = y, yend = y, colour = col),
                        linewidth = 0.25) +
  ggplot2::geom_segment(data = brk, inherit.aes = FALSE,
                        ggplot2::aes(x = x1, xend = x1, y = y, yend = y - tick, colour = col),
                        linewidth = 0.25) +
  ggplot2::geom_segment(data = brk, inherit.aes = FALSE,
                        ggplot2::aes(x = x2, xend = x2, y = y, yend = y - tick, colour = col),
                        linewidth = 0.25) +
  ggplot2::geom_text(data = brk, inherit.aes = FALSE,
                     ggplot2::aes(x = xmid, y = y, label = lab, colour = col),
                     vjust = -0.25, size = 1.9) +
  ggplot2::geom_blank(data = hr, ggplot2::aes(x = group, y = y), inherit.aes = FALSE) +
  ggplot2::facet_wrap(~ panel, nrow = 1, scales = "free_y") +
  ggplot2::scale_colour_manual(values = c(group_cols, sig = "#E41A1C", ns = "grey45"),
                               breaks = names(group_cols), labels = group_labels) +
  ggplot2::scale_fill_manual(values = group_cols, guide = "none") +
  ggplot2::labs(
    x = NULL, y = "fraction of the nuclear transcriptome (%)", colour = NULL,
    title = "Myc raises the nuclear mitochondrial arms; mtDNA-encoded output stays flat",
    caption = paste(
      sprintf("Points, n=%d/group; box = median/IQR. Brackets: genotype (WT vs Myc+) and 6W-vs-12W p (log2-share simple-effect models, script 32); red = p<0.05.", n_per),
      "Mito biogenesis = MitoCarta central dogma + protein import/sorting/homeostasis (314 genes; 7 shared with nuclear OXPHOS). One 12W WT mtDNA sample (67%) off-scale.",
      "Fraction of the nuclear transcriptome (13 mtDNA-encoded genes excluded); transcript share, not per-cell content, so a LOWER BOUND under global RNA amplification.",
      sep = "\n")) +
  theme_myc(base_size = 9) +
  ggplot2::theme(
    axis.text.x     = ggplot2::element_blank(),
    axis.ticks.x    = ggplot2::element_blank(),
    axis.line.x     = ggplot2::element_blank(),
    legend.position = "bottom",
    legend.key.size = ggplot2::unit(3.5, "mm")) +
  ggplot2::guides(colour = ggplot2::guide_legend(
    nrow = 1, override.aes = list(size = 2.4, alpha = 1, shape = 16)))

# Guard: when this script is sourced only to obtain `p` (e.g. from the Quarto doc),
# set options(myc.fig.nosave = TRUE) to display without re-exporting the canonical PDF.
if (!isTRUE(getOption("myc.fig.nosave"))) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  save_panel(p, file.path(out_dir, "fig01_mito_content.pdf"),
             width = fig_w[["double"]], height = 95)
}

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {
  print(brk[, c("panel", "cid", "kind", "p", "lab", "clean")])
  print(hr)
  print(n_per)
  print(p)
}
