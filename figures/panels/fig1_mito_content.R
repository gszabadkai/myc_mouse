# =============================================================================
# fig1_mito_content.R -- Myc builds MORE mitochondrion, and it does so across the
# whole compartment
# -----------------------------------------------------------------------------
# SLOT: Fig. 1E. The slot was freed when the partial-correlation paragraph was
# superseded (author, 2026-08-02); figures/panels/PANELS.md is the slug -> slot map.
#
# SUPPORTS one written sentence, which carries TWO claims and therefore needs two
# rulers:
#   "First, it drove a quantitative expansion by systemically upregulating 94% of
#    143 nuclear-encoded mitochondrial pathways, increasing the mitochondrial
#    transcriptomic fraction by 21-27% (padj < 0.001), demonstrating robust
#    mitochondrial biogenesis (Fig. 1E)."
#
#   TOP     the ABSOLUTE claim -- what fraction of the transcriptome is
#           mitochondrial, per animal, five arms, four groups.
#   BOTTOM  the SYSTEMIC claim -- the Myc content effect across all 143
#           nuclear-encoded MitoPathways, so the reader can see that almost the
#           whole compartment moves and not a favoured corner of it.
#
# Fig. 1F is the SAME 143 pathways on the content-blind priority ruler. The pair
# is the two-mechanism sentence: a one-sided distribution here, a two-sided one
# there. That is why the tier colour key lives on 1F and this panel stays neutral.
#
# WHICH TEST THE BRACKET DRAWS, AND WHY IT IS THE POOLED ONE. The number the
# sentence quotes is the genotype MAIN effect, lm(log2 share ~ timepoint +
# myc_status), which is script 32's own model and is reproduced here to 1e-8. The
# groups are drawn genotype-major (both WT boxes, then both Myc+ boxes), so one
# bracket spanning the two halves IS that contrast. The additive model is the
# right one: every arm's genotype x timepoint interaction is far from significant
# (p 0.35-0.81), i.e. the content effect does not differ between the two ages.
# Split by age it is present at both and clears p<0.05 only at 12 weeks (6W p
# 0.03-0.11, 12W p 0.0015-0.016) -- at n=6 per cell, which is why the pooled test
# is the one that answers the sentence. Both splits are in the legend block.
#
# Reads (read-only, no re-run):
#   results/mito_content_proxies.rds (script 32) -- $shares, per-sample percentage
#       of the NUCLEAR transcriptome (denominator drops only the 13 mt-* genes, so
#       nuclear MitoCarta stays in it); $share_stats for the assertion.
#   results/background_vs_myc.rds (script 40) -- $ruler, set-average RAW log2FC per
#       MitoPathway (CLAUDE.md: averaged-LFC visuals run on unshrunken LFCs);
#       $ruler_summary for the assertion.
# Output: outputs/figures/panels/fig1_mito_content.pdf
# =============================================================================

source(here::here("figures", "panels", "_panel_common.R"))
if (!requireNamespace("patchwork", quietly = TRUE)) stop("Fig. 1E needs patchwork")

content_path <- here::here("results", "mito_content_proxies.rds")
require_fresher_than(content_path)
content <- readRDS(content_path)
bg      <- readRDS(here::here("results", "background_vs_myc.rds"))

# =============================================================================
# TOP -- the absolute claim: five arms, four groups, every animal
# =============================================================================

# Ordered nuclear-up first, the flat mtDNA-encoded arm last. It is drawn because
# it is the internal control for "is this just more RNA?": a global amplification
# would lift it too, and it does not move. It is NOT the mitonuclear imbalance --
# that claim was retracted (script 32, 2026-07-16) and the bounds say so.
arms <- c("MASS_MARKERS_NOCHAP",
          "MITOCARTA_NUCLEAR_ENCODED",
          "MITOCARTA_OXPHOS_NU",
          "BIOGENESIS_FULL",
          "MITOCARTA_MTDNA_ENCODED")
# Two lines each: at 89 mm a facet is ~15 mm wide and a one-line "Nuclear
# MitoCarta" is wider than that, so it clips. The break is where the qualifier is.
arm_name <- c(MASS_MARKERS_NOCHAP       = "Mass\nmarkers",
              MITOCARTA_NUCLEAR_ENCODED = "Nuclear\nMitoCarta",
              MITOCARTA_OXPHOS_NU       = "Nuclear\nOXPHOS",
              BIOGENESIS_FULL           = "Mito\nbiogenesis",
              MITOCARTA_MTDNA_ENCODED   = "mtDNA-\nencoded")

shares <- as.data.frame(content$shares)
# as.data.frame first: script 32 saves tibbles, and a tibble's [ returns a tibble
# rather than a scalar, which would silently poison every sprintf below.
stats  <- as.data.frame(content$share_stats)
stats  <- stats[stats$denominator == "share_nomt", ]
rownames(stats) <- stats$panel

df <- shares[shares$panel %in% arms, ]
df$panel <- factor(df$panel, levels = arms, labels = arm_name[arms])
# genotype-major, which is names(group_cols): 6W_wt, 12W_wt | 6W_myc, 12W_myc
df$group <- factor(df$group, levels = names(group_cols))
stopifnot(nrow(df) == length(arms) * 24, !any(is.na(df$group)),
          all(table(df$panel, df$group) == 6L))

# --- script 32's model, reproduced ------------------------------------------
# lm(log2 share ~ timepoint + myc_status). `simple` gives the four simple effects
# for the legend block; `main` gives the drawn contrast.
fit_main <- function(d) {
  d$yv <- log2(d$share_nomt)
  d$tp <- factor(d$timepoint,  levels = c("6W", "12W"))
  d$mc <- factor(d$myc_status, levels = c("neg", "pos"))
  co <- summary(stats::lm(yv ~ tp + mc, data = d))$coefficients
  c(beta = unname(co["mcpos", 1]), p = unname(co["mcpos", 4]))
}
fit_simple <- function(d, kind, key) {
  d$yv <- log2(d$share_nomt)
  m <- if (kind == "geno") {
    stats::lm(yv ~ myc_status, data = d[d$timepoint == key, ])
  } else {
    stats::lm(yv ~ timepoint,  data = d[d$myc_status == key, ])
  }
  co <- summary(m)$coefficients
  c(beta = unname(co[2, 1]), p = unname(co[2, 4]))
}
arm_fit <- vapply(arms, function(a) fit_main(shares[shares$panel == a, ]),
                  numeric(2))

# ASSERTION: the drawn contrast IS the analysis of record, not a lookalike.
stopifnot(
  max(abs(arm_fit["beta", ] - stats[arms, "geno_beta"])) < 1e-8,
  max(abs(arm_fit["p",    ] - stats[arms, "geno_p"]))    < 1e-8)

# --- off-scale handling ------------------------------------------------------
# One 12W wild-type animal sits at ~67% mtDNA-encoded while every other point in
# that arm is under 40. Capping the DISPLAY and printing the true value keeps the
# other 23 points readable without deleting an animal.
CAP <- 40
df$y_disp <- pmin(df$share_nomt, CAP)
df$capped <- df$share_nomt > CAP

fmt_p <- function(p) if (p < 0.001) "<0.001" else formatC(p, format = "g", digits = 2)

# One bracket per facet: the two WT boxes (x 1,2) against the two Myc+ boxes
# (x 3,4), which is the genotype main effect the sentence quotes.
brk <- do.call(rbind, lapply(arms, function(a) {
  d      <- df[df$panel == arm_name[[a]], ]
  facmax <- max(d$y_disp)
  data.frame(panel = factor(arm_name[[a]], levels = arm_name[arms]),
             x = 1.5, xend = 3.5, xmid = 2.5,
             y     = facmax * 1.07,
             y_lab = facmax * 1.10,       # label sits ABOVE the bar, not on it
             tick  = facmax * 0.028,
             p    = unname(arm_fit["p", a]),
             lab  = fmt_p(unname(arm_fit["p", a])),
             stringsAsFactors = FALSE)
}))
brk$col <- ifelse(brk$p < 0.05, "sig", "ns")

hr <- do.call(rbind, lapply(arms, function(a) {
  d <- df[df$panel == arm_name[[a]], ]
  data.frame(panel = factor(arm_name[[a]], levels = arm_name[arms]),
             group = factor("6W_neg", levels = names(group_cols)),
             y     = max(d$y_disp) * 1.20)
}))

pts_layer <- function(dat, size) {
  if (requireNamespace("ggbeeswarm", quietly = TRUE)) {
    ggbeeswarm::geom_quasirandom(data = dat, width = 0.24, size = size, alpha = 0.95)
  } else {
    ggplot2::geom_jitter(data = dat, width = 0.16, height = 0, size = size, alpha = 0.95)
  }
}

p_top <- ggplot2::ggplot(df, ggplot2::aes(group, y_disp,
                                          colour = group, fill = group)) +
  ggplot2::geom_boxplot(outlier.shape = NA, width = 0.68, alpha = 0.28,
                        colour = "grey35", linewidth = 0.25) +
  pts_layer(df[!df$capped, ], size = 0.7) +
  ggplot2::geom_point(data = df[df$capped, ], shape = 17, size = 1.1) +
  ggplot2::geom_text(data = df[df$capped, ],
                     ggplot2::aes(label = sprintf("%.0f", share_nomt)),
                     hjust = -0.45, vjust = 0.5, size = 1.6, colour = "grey30",
                     show.legend = FALSE) +
  ggplot2::geom_segment(data = brk, inherit.aes = FALSE,
                        ggplot2::aes(x = x, xend = xend, y = y, yend = y,
                                     colour = col), linewidth = 0.22) +
  ggplot2::geom_segment(data = brk, inherit.aes = FALSE,
                        ggplot2::aes(x = x, xend = x, y = y, yend = y - tick,
                                     colour = col), linewidth = 0.22) +
  ggplot2::geom_segment(data = brk, inherit.aes = FALSE,
                        ggplot2::aes(x = xend, xend = xend, y = y, yend = y - tick,
                                     colour = col), linewidth = 0.22) +
  ggplot2::geom_text(data = brk, inherit.aes = FALSE,
                     ggplot2::aes(x = xmid, y = y_lab, label = lab, colour = col),
                     vjust = 0, size = 1.75) +
  ggplot2::geom_blank(data = hr, ggplot2::aes(x = group, y = y), inherit.aes = FALSE) +
  ggplot2::facet_wrap(~ panel, nrow = 1, scales = "free_y") +
  ggplot2::scale_colour_manual(values = c(group_cols, sig_cols),
                               breaks = names(group_cols), labels = group_labels) +
  ggplot2::scale_fill_manual(values = group_cols, guide = "none") +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.08, 0.02))) +
  ggplot2::labs(x = NULL, y = "% of nuclear transcriptome", colour = NULL) +
  ggplot2::guides(colour = ggplot2::guide_legend(
    nrow = 1, override.aes = list(size = 1.5, alpha = 1, shape = 16, label = ""))) +
  theme_panel(base_size = 6) +
  ggplot2::theme(
    axis.text.x     = ggplot2::element_blank(),
    axis.ticks.x    = ggplot2::element_blank(),
    axis.line.x     = ggplot2::element_blank(),
    axis.text.y     = ggplot2::element_text(size = 5.2),
    panel.spacing.x = ggplot2::unit(1.2, "mm"),
    strip.clip      = "off",
    strip.text      = ggplot2::element_text(face = "plain", size = 5.4, lineheight = 0.95,
                                            margin = ggplot2::margin(0, 0, 0.6, 0, "mm")),
    plot.margin     = ggplot2::margin(1, 1.5, 0.5, 1.5, "mm"))

# =============================================================================
# BOTTOM -- the systemic claim: the whole compartment, one point per pathway
# =============================================================================
# `c_m6` is the set-average RAW (unshrunken) log2 fold change of the Myc genotype
# effect at 6 weeks. The synthetic mtDNA-encoded pathway carries `is_mtdna` and is
# excluded from every fit in script 40, so it is excluded here too and the panel
# and the regressions describe the same 143 pathways.
ruler <- as.data.frame(bg$ruler)                    # tibbles, as above
r143  <- ruler[!ruler$is_mtdna, ]
rs    <- as.data.frame(bg$ruler_summary)
rs6   <- rs[rs$ruler == "content" & rs$metric == "m6", ]

pct_up_143 <- 100 * mean(r143$c_m6 > 0)
med_143    <- stats::median(r143$c_m6)

# ASSERTION: the drawn 143 are script 40's 143, and the all-144 summary this
# reproduces is the one the manuscript's "94%" currently quotes.
stopifnot(
  nrow(r143) == 143L, nrow(ruler) == 144L,
  abs(100 * mean(ruler$c_m6 > 0) - rs6$pct_up)      < 1e-6,
  abs(stats::median(ruler$c_m6)  - rs6$median)      < 1e-6)

r143$dir <- ifelse(r143$c_m6 > 0, "up", "down")
dn <- stats::density(r143$c_m6, adjust = 0.9)
PT <- -max(dn$y) * 0.14        # the point strip sits under the curve's baseline

p_bot <- ggplot2::ggplot(r143, ggplot2::aes(x = c_m6)) +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey55") +
  ggplot2::geom_density(adjust = 0.9, fill = "grey88", colour = "grey40",
                        linewidth = 0.3) +
  ggplot2::geom_segment(x = med_143, xend = med_143, y = 0,
                        yend = max(dn$y) * 1.02, linewidth = 0.3,
                        colour = "grey20", linetype = "22") +
  ggplot2::geom_jitter(ggplot2::aes(y = PT, colour = dir), height = max(dn$y) * 0.06,
                       width = 0, size = 0.55, alpha = 0.85, show.legend = FALSE) +
  ggplot2::scale_colour_manual(values = direction_cols) +
  ggplot2::scale_x_continuous(labels = lab_signed,
                              expand = ggplot2::expansion(mult = c(0.03, 0.03))) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.16, 0.06))) +
  ggplot2::labs(x = "Myc content effect at 6W, per MitoPathway (log2)", y = NULL) +
  theme_panel(base_size = 6) +
  ggplot2::theme(
    axis.text.y  = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    axis.line.y  = ggplot2::element_blank(),
    plot.margin  = ggplot2::margin(0.5, 1.5, 1, 1.5, "mm"))

p <- patchwork::wrap_plots(p_top, p_bot, ncol = 1, heights = c(2.5, 1)) +
  patchwork::plot_layout(guides = "collect") &
  ggplot2::theme(legend.position = "bottom",
                 legend.margin   = ggplot2::margin(-2, 0, 0, 0),
                 legend.key.size = ggplot2::unit(2.4, "mm"))

# =============================================================================
# the legend text (never drawn)
# =============================================================================
pc  <- function(b) sprintf("%+.0f%%", 100 * (2^b - 1))
# "p = <0.001" reads badly; make the operator part of the string.
p_txt <- function(p) if (p < 0.001) "p < 0.001" else sprintf("p = %s", fmt_p(p))
sim <- function(a, kind, key) fit_simple(shares[shares$panel == a, ], kind, key)
arm_line <- function(a) {
  g6 <- sim(a, "geno", "6W"); g12 <- sim(a, "geno", "12W")
  sprintf("%s (%d genes): %s overall, %s; %s at 6 weeks (p = %.3f) and %s at 12 (p = %.4f)",
          gsub("\n", " ", arm_name[[a]]), stats[a, "n_genes"],
          pc(arm_fit["beta", a]), p_txt(arm_fit["p", a]),
          pc(g6[["beta"]]), g6[["p"]], pc(g12[["beta"]]), g12[["p"]])
}

LEGEND <- panel_legend(
  slot = "Fig. 1E",
  what = paste0(
    "Myc raises mitochondrial content. TOP: the percentage of the nuclear ",
    "transcriptome contributed by five mitochondrial arms, one point per animal. ",
    "BOTTOM: the Myc genotype effect at six weeks on each of the 143 ",
    "nuclear-encoded MitoPathways, as a distribution -- the same 143 pathways ",
    "Fig. 1F ranks on the content-blind priority ruler."),
  detail = c(
    "TOP: n = 6 animals per group, n = 24 per arm; every animal is drawn. Boxes are median and quartiles, whiskers 1.5x the interquartile range. Groups run genotype-major (both wild-type boxes, then both Myc+ boxes), so the single bracket spans the contrast it labels. The x axis is blank because the colour key carries the group names.",
    "The bracket is the genotype MAIN effect from script 32's own model, lm(log2 share ~ timepoint + myc_status), reproduced here and asserted equal to results/mito_content_proxies.rds$share_stats to 1e-8. Red marks p < 0.05. Per arm, overall and split by age:",
    arm_line("MASS_MARKERS_NOCHAP"),
    arm_line("MITOCARTA_NUCLEAR_ENCODED"),
    arm_line("MITOCARTA_OXPHOS_NU"),
    arm_line("BIOGENESIS_FULL"),
    arm_line("MITOCARTA_MTDNA_ENCODED"),
    sprintf("The additive model is the correct one here: the genotype x timepoint interaction is far from significant on every arm (p = %s), so the content effect does not differ between the two ages.",
            paste(sprintf("%.2f", stats[arms, "int_p"]), collapse = ", ")),
    "Mass markers are 12 canonical abundance proteins with the chaperones removed; Mito biogenesis is the MitoCarta central dogma plus protein import, sorting and homeostasis (314 genes, of which 7 are also in nuclear OXPHOS, so the two arms are effectively disjoint). One 12-week wild-type animal sits at 67% mtDNA-encoded and is drawn as a triangle at the axis cap of 40% with its true value printed.",
    sprintf("BOTTOM: one point per MitoPathway, coloured by sign, over a kernel density; the dashed line is the median (%+.3f) and the solid line is zero. %.1f%% of the %d pathways are above zero, %d below. Values are set-average RAW, unshrunken log2 fold changes, which is what CLAUDE.md requires of an averaged-LFC visual.",
            med_143, pct_up_143, nrow(r143), sum(r143$c_m6 <= 0)),
    sprintf("The manuscript sentence currently reads \"94%% of 143\", which mixes two numbers: %.1f%% is the fraction over ALL 144 rows of script 40's ruler (the synthetic mtDNA-encoded pathway included) and %.1f%% is the fraction over the 143 the regressions and this panel use. Pick one.",
            rs6$pct_up, pct_up_143)),
  bounds = c(
    "These are transcript SHARES, not per-cell content. Myc amplifies global transcription, which inflates the denominator, so every share here is a LOWER BOUND on the content increase. Blot, qPCR or EM settles the absolute question; this panel bounds it from below.",
    "A share can in principle rise because other genes fall. It does not here: the same arms rise in absolute per-gene log2 fold change, 13 of the 14 mass-marker genes move up with Myc and 11 individually clear p < 0.05, and the effect survives adjustment for dissociation stress and for residual non-epithelial contamination (script 33).",
    "The quoted p-value is the pooled genotype main effect. Split by age at n = 6 per cell it clears p < 0.05 at twelve weeks on all four nuclear arms and only on Mito biogenesis at six (p = 0.034), the others sitting at p = 0.08-0.11. The effect size is essentially the same at both ages; what differs is power.",
    "THE mtDNA ARM IS A CONTROL, NOT A MITONUCLEAR IMBALANCE. Script 32 withdrew that reading on 2026-07-16: tested as a difference (nuclear minus mtDNA, per sample) it gives p = 0.57 raw and p = 0.82 composition-adjusted, and the arm is underpowered (13 genes, share 3.4-40% across animals). It is drawn because a global RNA amplification would lift it too and it does not move with genotype; it is not evidence for or against an imbalance.",
    "Every genotype comparison here is clean -- genotype is balanced within each extraction batch. The two 6W-versus-12W comparisons are not drawn because batch is aligned with timepoint by design; they are in the sandbox and in script 32.",
    "The bottom distribution is descriptive. 33 of the 143 pathways clear padj < 0.05 individually at six weeks; MitoPathways nest inside one another, so the 143 are not independent and the shape of the distribution is the statement, not any single point."),
  source = c(
    "results/mito_content_proxies.rds (scripts/32_mito_content_proxies.R) -- $shares, $share_stats",
    "results/background_vs_myc.rds (scripts/40_background_vs_myc_decomposition.R) -- $ruler, $ruler_summary",
    "Earlier double-column form of the top row: figures/fig01_mito_content.R; its per-gene log2FC companion, figures/fig01b_mito_content_lfc.R"))

save_panel_p(p, "fig1_mito_content", height = 78)

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {

  print(p)
  print(LEGEND)

  ## every arm, all four simple effects plus the main effect, in percent
  do.call(rbind, lapply(arms, function(a) data.frame(
    arm      = arm_name[[a]],
    n_genes  = stats[a, "n_genes"],
    main_pct = pc(arm_fit["beta", a]),  main_p = signif(arm_fit["p", a], 3),
    myc_6W   = pc(sim(a, "geno", "6W")[["beta"]]),
    p_6W     = signif(sim(a, "geno", "6W")[["p"]], 3),
    myc_12W  = pc(sim(a, "geno", "12W")[["beta"]]),
    p_12W    = signif(sim(a, "geno", "12W")[["p"]], 3),
    wt_time  = pc(sim(a, "time", "neg")[["beta"]]),
    p_wt     = signif(sim(a, "time", "neg")[["p"]], 3),
    myc_time = pc(sim(a, "time", "pos")[["beta"]]),
    p_myc    = signif(sim(a, "time", "pos")[["p"]], 3),
    int_p    = signif(stats[a, "int_p"], 3),
    row.names = NULL))) |> print()

  ## the bottom distribution, by tier -- the 7 pathways that go DOWN
  r143[r143$c_m6 <= 0, c("pathway", "tier", "n_genes", "c_m6", "p_m6_padj")] |>
    print()

  ## leave-one-out sensitivity of the drawn contrast (script 32's own check)
  content$loo_sensitivity[content$loo_sensitivity$panel %in% arms, ] |> print()

  ## the four-arm fallback if the five-facet row is too tight at 89 mm:
  ## nuclear OXPHOS is the one arm Figs. 1F and 2F both carry.
  ## arms <- setdiff(arms, "MITOCARTA_OXPHOS_NU")
}
