# =============================================================================
# explainer_nes_sign.R -- why a rank-STABLE Myc effect gives a NEGATIVE NES on
# the Myc+ timeline
# -----------------------------------------------------------------------------
# SLOT: none, and deliberately so. The filename does not start with "fig", so it
# falls outside the `^fig.*\.R$` glob that rebuild_panels.R:44 and
# panels_to_pdf.R:65 both use. It therefore does not enter the 28-panel count,
# does not need a slot letter, and does not need a chapter in
# paper/analysis_record.qmd. It is a TEACHING figure, built 2026-08-16 to answer
# the author's question, and it is the only figure in this layer whose top three
# quarters are SIMULATED. That is marked on the page, not only here.
#
# THE QUESTION. Fig. 1D/1H say the Myc programme keeps its ranking between six
# and twelve weeks (Spearman 0.933 over 866 sets) while the effect size halves
# (Fig. 1G). scripts/04 -> outputs/fgsea/ says most pathways have a strongly
# NEGATIVE NES on the two timelines. If the ranking does not move, why does the
# timeline move at all?
#
# THE ANSWER, and it is one sentence: fGSEA never compares two rankings. It takes
# ONE ranked gene list and asks where a set sits in it -- and there are FOUR
# lists here, not two.
#
#     myc_6W      the Myc effect at 6W        Myc+ minus WT, at six weeks
#     myc_12W     the Myc effect at 12W       Myc+ minus WT, at twelve weeks
#     6>12W_wt    the wild-type timeline      12W minus 6W, WITHIN WT animals
#     6>12W_myc   the Myc+ timeline           12W minus 6W, WITHIN Myc+ animals
#
# "The ranking is preserved" is a statement about the first two. The last two are
# different lists with their own gene order, and the Myc+ one is built from the
# DIFFERENCE between the first two -- the part that was lost. A difference of two
# proportional things is not zero. Halving is a change: a target 2-fold over wild
# type at six weeks and 1.4-fold at twelve went DOWN inside the Myc+ animals, and
# so did every other target, in the same order. So the Myc+ timeline is the
# NEGATIVE of the Myc effect, and the set that tops the first two lists sits at
# the BOTTOM of it.
#
# So the negative timeline NES is not evidence against stability. It is what
# stability PRODUCES. A rank-UNSTABLE fade would give a weak, incoherent timeline
# NES; that the timeline is a clean mirror is further evidence the shape held.
#
# REVISED 2026-08-16, same day, on the author's reading of the first version:
# "it explains the preserved gene ranking of the Myc effects at the different
# timepoints, but does not show the ranking in the 6>12W_myc (nor the 6>12W WT)".
# Correct, and the first version could not have: its simulation had a STATIC wild
# type, so there was no 6>12W_wt list to draw. The simulation now carries a
# wild-type trajectory, and REGION B draws the gene order of every list against
# the gene order of myc_6W -- which is the thing the question was actually about.
#
# WHAT THE PANEL DRAWS
#   A   the premise. Every gene's Myc effect at each age, halved exactly and with
#       no reordering. The LINE connecting the two ages IS the Myc+ timeline when
#       the wild type is static, and every line moves toward zero. Ranks cannot
#       show amplitude, which is why this region exists at all.
#   B   THE ANSWER. Gene rank in myc_6W against gene rank in each of the other
#       three lists, all 2,000 genes, set members picked out. A perfect diagonal,
#       a cloud, and an anti-diagonal -- one picture per list, and each list's own
#       gene order is on the page.
#   C   what fGSEA computes from those four orders: the running enrichment score
#       of the same set in each, with the NES it produces.
#   D   the real data, for two named sets, on the same four contrasts -- so the
#       toy can be checked against the thing it claims to explain, and so the one
#       way the real data DIFFERS from the toy is visible: for the OXPHOS subunits
#       the wild-type gland moves as much as the fade does.
#
# Reads (read-only, no re-run):
#   results/fgsea_percategory.rds (script 20) -- $fgsea, REGION D and two numbers
# Output: outputs/figures/panels/explainer_nes_sign.pdf
# =============================================================================

source(here::here("figures", "panels", "_panel_common.R"))

stopifnot(requireNamespace("fgsea", quietly = TRUE),
          requireNamespace("patchwork", quietly = TRUE))

# =============================================================================
# THE SIMULATION
# =============================================================================
# The simplest world in which all four lists exist:
#   * the twelve-week Myc effect is EXACTLY half the six-week one, gene for gene,
#     so the rank correlation between those two lists is 1 by construction and no
#     reordering can be smuggled in;
#   * the wild-type gland HAS a trajectory, and it is independent of the Myc
#     effect. That is the measured situation, not a simplification: over the 866
#     library sets Spearman(myc_6W, 6>12W_wt) = -0.068. The real wild-type
#     timeline has strong structure of its own -- it is simply structure on
#     DIFFERENT genes, which is what independence means here;
#   * the ranking statistic is the effect itself, i.e. constant standard error.
#     That is also measured: over the 1,967 genes Myc moves at six weeks the
#     median lfcSE is 0.189 at six weeks and 0.193 at twelve, so the fade is in
#     the numerator alone and is not a change in power.

set.seed(42)

N_GENE  <- 2000L
N_SET   <- 120L
FADE    <- 0.5     # the fraction of the Myc effect surviving to twelve weeks
DEV_SD  <- 0.25    # the wild-type gland's own trajectory, per gene

genes   <- sprintf("g%04d", seq_len(N_GENE))
set_ids <- genes[seq_len(N_SET)]   # the first N_SET genes are the "Myc target" set

# The Myc effect at six weeks: the set is induced, the background is not.
e6 <- stats::rnorm(N_GENE, mean = 0, sd = 0.40)
e6[seq_len(N_SET)] <- stats::rnorm(N_SET, mean = 1.20, sd = 0.40)
names(e6) <- genes

# The Myc effect at twelve weeks: the SAME numbers, times a constant.
e12 <- FADE * e6

# The wild-type gland's trajectory. Independent of the Myc effect and of set
# membership -- so this list is the control: it is a real, structured list, and
# the Myc-target set has no business being anywhere in particular within it.
dev <- stats::rnorm(N_GENE, mean = 0, sd = DEV_SD)
names(dev) <- genes

# The two timelines, and the interaction.
tl_wt  <- dev                    # 12W_wt  - 6W_wt
tl_myc <- dev + (e12 - e6)       # 12W_myc - 6W_myc = development + the fade
inter  <- e12 - e6               # myc_12W - myc_6W = -(1 - FADE) * e6

# The identities the whole explanation rests on, asserted rather than claimed.
stopifnot(
  # the two genotype lists are in the SAME order -- no reordering anywhere
  isTRUE(all.equal(stats::cor(e6, e12, method = "spearman"), 1)),
  # the INTERACTION is the exact mirror of the Myc effect
  isTRUE(all.equal(stats::cor(e6, inter, method = "spearman"), -1)),
  isTRUE(all.equal(unname(inter), unname(-(1 - FADE) * e6))),
  # the Myc+ timeline is that mirror PLUS development, so it is a strong but not
  # perfect mirror -- which is exactly the real situation
  isTRUE(all.equal(unname(tl_myc), unname(dev + inter))),
  stats::cor(e6, tl_myc, method = "spearman") < -0.5,
  # and the wild-type timeline has nothing to do with the Myc effect
  abs(stats::cor(e6, tl_wt, method = "spearman")) < 0.10)

LISTS <- list("myc_6W"    = e6,     "myc_12W"   = e12,
              "6>12W_wt"  = tl_wt,  "6>12W_myc" = tl_myc)
stopifnot(identical(names(LISTS), contrast_levels))   # the declared vocabulary

# =============================================================================
# FOUR REAL fGSEA RUNS -- every NES on the panel is computed, never asserted
# =============================================================================
set.seed(42)
sim_nes <- do.call(rbind, lapply(names(LISTS), function(nm) {
  r <- sort(LISTS[[nm]], decreasing = TRUE)
  f <- fgsea::fgsea(pathways = list(MYC_TARGET_SET = set_ids), stats = r,
                    minSize = 10, maxSize = 500, eps = 0, nPermSimple = 10000)
  data.frame(list_name = nm, NES = f$NES[1], padj = f$padj[1],
             stringsAsFactors = FALSE)
}))
sim_nes$list_name <- factor(sim_nes$list_name, levels = contrast_levels)
nes_at <- function(nm) sim_nes$NES[sim_nes$list_name == nm]

# The result the panel exists to show: the two genotype lists agree because they
# ARE the same ranking, the wild-type timeline says nothing about this set, and
# the Myc+ timeline flips sign.
stopifnot(nes_at("myc_6W") > 2, nes_at("myc_12W") > 2,
          abs(nes_at("myc_6W") - nes_at("myc_12W")) < 0.15,
          abs(nes_at("6>12W_wt")) < 1.4,
          nes_at("6>12W_myc") < -2)

# --- the running enrichment score, the classic weighted-KS walk ---------------
# Implemented here rather than taken from fgsea::plotEnrichment() so the four
# facets share one geometry and one set of aesthetics.
running_es <- function(stats_vec, set_genes) {
  r   <- sort(stats_vec, decreasing = TRUE)
  hit <- names(r) %in% set_genes
  inc <- abs(r) * hit / sum(abs(r[hit]))
  dec <- (!hit) / sum(!hit)
  data.frame(rank = seq_along(r), es = cumsum(inc - dec), hit = hit)
}

es_df <- do.call(rbind, lapply(names(LISTS), function(nm) {
  d <- running_es(LISTS[[nm]], set_ids); d$list_name <- nm; d
}))
es_df$list_name <- factor(es_df$list_name, levels = contrast_levels)

FAM_COLS <- c("Myc-target set" = unname(contrast_cols[["myc_6W"]]),
              "background"     = "grey78")
fam_of <- function(g) factor(ifelse(g %in% set_ids, "Myc-target set", "background"),
                             levels = c("Myc-target set", "background"))

# =============================================================================
# REGION A -- the premise: same order, half the size
# =============================================================================
# A SAMPLE of genes is drawn, because 2,000 lines at 26 mm is a black rectangle.
# The sample is stratified so the set's whole range is represented; the rank
# correlation quoted in the legend is over all 2,000.
set.seed(7)
show <- c(set_ids[round(seq(1, N_SET, length.out = 11))],
          sample(genes[(N_SET + 1L):N_GENE], 15L))

slope_df <- data.frame(
  gene = rep(show, 2L),
  age  = factor(rep(contrast_geno, each = length(show)), levels = contrast_geno),
  lfc  = c(e6[show], e12[show]),
  fam  = rep(fam_of(show), 2L))

pA <- ggplot2::ggplot(slope_df,
                      ggplot2::aes(age, lfc, group = gene, colour = fam)) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey55") +
  ggplot2::geom_line(linewidth = 0.34) +
  ggplot2::geom_point(size = 0.6) +
  ggplot2::scale_colour_manual(values = FAM_COLS, name = NULL) +
  ggplot2::scale_x_discrete(expand = ggplot2::expansion(add = c(0.30, 0.30))) +
  ggplot2::scale_y_continuous(labels = lab_signed) +
  ggplot2::labs(x = NULL, y = "log2 fold change\nvs wild type") +
  theme_panel() +
  ggplot2::theme(legend.position = "top",
                 legend.margin = ggplot2::margin(0, 0, -2, 0),
                 axis.text.x = ggplot2::element_text(face = "bold"),
                 panel.grid.major.x = ggplot2::element_blank())

# =============================================================================
# REGION B -- THE ANSWER: each list's own gene order, against myc_6W's
# =============================================================================
# Rank 1 is the top of the list. Drawn for ALL 2,000 genes, so this is the whole
# ranking and not a sample: a perfect diagonal means the order is identical, a
# cloud means the two lists are unrelated, and an anti-diagonal means one is the
# reverse of the other.
rk <- function(v) rank(-v, ties.method = "first")

rank_df <- do.call(rbind, lapply(setdiff(contrast_levels, "myc_6W"), function(nm) {
  data.frame(x = rk(e6), y = rk(LISTS[[nm]]),
             fam = fam_of(genes), list_name = nm)
}))
rank_df$list_name <- factor(rank_df$list_name, levels = contrast_levels)
# set members last, so they are not buried under 1,880 grey points
rank_df <- rank_df[order(rank_df$list_name, rank_df$fam == "Myc-target set"), ]

rho_df <- data.frame(
  list_name = factor(setdiff(contrast_levels, "myc_6W"), levels = contrast_levels),
  rho = vapply(setdiff(contrast_levels, "myc_6W"),
               function(nm) stats::cor(e6, LISTS[[nm]], method = "spearman"),
               numeric(1)))
stopifnot(rho_df$rho[1] > 0.99, abs(rho_df$rho[2]) < 0.10, rho_df$rho[3] < -0.5)

# rho goes in the STRIP, not inside the panel. No corner is free in all three
# facets -- the diagonal fills one pair, the anti-diagonal the other and the cloud
# fills every one -- so a label placed inside would sit on data in at least one.
strip_B <- stats::setNames(
  sprintf("%s\nrho %+.2f", as.character(rho_df$list_name), rho_df$rho),
  as.character(rho_df$list_name))

# BOTH AXES RUN 1 -> 2000 UPWARD AND RIGHTWARD, so agreement reads as the familiar
# "/" and reversal as "\". Reversing the y axis to put rank 1 at the top drew the
# PERFECT agreement facet as a "\", which reads as anti-correlation to anyone who
# has ever seen a scatter plot. Rank 1 is the top of the list; the axis titles say
# so, and that is cheaper than fighting the reader's eye.
pB <- ggplot2::ggplot(rank_df, ggplot2::aes(x, y, colour = fam)) +
  ggplot2::geom_point(size = 0.09, alpha = 0.55, show.legend = FALSE) +
  ggplot2::facet_wrap(~ list_name, nrow = 1,
                      labeller = ggplot2::labeller(list_name = strip_B)) +
  ggplot2::scale_colour_manual(values = FAM_COLS) +
  ggplot2::scale_x_continuous(breaks = c(1, N_GENE), labels = c("1", "2000")) +
  ggplot2::scale_y_continuous(breaks = c(1, N_GENE), labels = c("1", "2000")) +
  ggplot2::labs(x = "gene rank in myc_6W   (1 = top)",
                y = "gene rank in\nthat list  (1 = top)") +
  ggplot2::coord_fixed() +
  theme_panel() +
  ggplot2::theme(strip.text = ggplot2::element_text(face = "bold", size = 6.0,
                                                    lineheight = 0.95),
                 panel.spacing.x = ggplot2::unit(2.0, "mm"))

# =============================================================================
# REGION C -- what fGSEA computes from those four orders
# =============================================================================
TICK_LO <- -0.10
TICK_HI <- -0.02

tick_df <- es_df[es_df$hit, c("rank", "list_name")]

nes_lab <- sim_nes
# Significance is on the label because the wild-type list's +0.95 is NOT a weak
# positive result, it is nothing (padj 0.56), and a reader who takes it for a
# small enrichment has learned the wrong thing from this figure.
nes_lab$lab  <- sprintf("NES\n%+.2f%s", nes_lab$NES,
                        ifelse(nes_lab$padj < 0.05, "", "\n(ns)"))
nes_lab$xpos <- N_GENE * 0.96
nes_lab$ypos <- 0.98
# NEUTRAL INK, deliberately. The strip names the list and the curve carries its
# declared colour, so the number is unambiguous already; spending a third
# encoding on it would either repeat the strip or invent a sign colour that the
# curve's own direction already shows.

pC <- ggplot2::ggplot(es_df, ggplot2::aes(rank, es)) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey55") +
  ggplot2::geom_line(ggplot2::aes(colour = list_name), linewidth = 0.42) +
  ggplot2::geom_segment(data = tick_df,
                        ggplot2::aes(x = rank, xend = rank,
                                     y = TICK_LO, yend = TICK_HI),
                        inherit.aes = FALSE, linewidth = 0.14, colour = "grey30") +
  ggplot2::geom_text(data = nes_lab,
                     ggplot2::aes(xpos, ypos, label = lab), inherit.aes = FALSE,
                     hjust = 1, vjust = 1, size = 1.95, lineheight = 0.9,
                     fontface = "bold", colour = "grey15") +
  # Name the tick row once, in the first facet, rather than leaving the reader to
  # infer it from the GSEA idiom.
  ggplot2::geom_text(data = data.frame(list_name = factor("myc_6W",
                                                          levels = contrast_levels)),
                     ggplot2::aes(x = 60, y = -0.26, label = "set members"),
                     inherit.aes = FALSE, hjust = 0, size = 1.75, colour = "grey30") +
  ggplot2::facet_wrap(~ list_name, nrow = 1) +
  ggplot2::scale_colour_manual(values = contrast_cols, guide = "none") +
  ggplot2::scale_x_continuous(breaks = c(1, N_GENE), labels = c("1", "2000")) +
  ggplot2::scale_y_continuous(labels = lab_signed, limits = c(-1.02, 1.02)) +
  ggplot2::labs(x = "gene rank in that list (high to low)",
                y = "running\nenrichment") +
  theme_panel() +
  ggplot2::theme(strip.text = ggplot2::element_text(face = "bold", size = 6.2),
                 panel.spacing.x = ggplot2::unit(2.0, "mm"))

# =============================================================================
# REGION D -- the real data, so the toy can be checked against it
# =============================================================================
fg_path <- here::here("results", "fgsea_percategory.rds")
fg <- readRDS(fg_path)$fgsea

RANK_MAP <- c(myc_6W = "myc_6W", myc_12W = "myc_12W",
              timepoint_neg = "6>12W_wt", timepoint_pos = "6>12W_myc")

# The compression is TOWARD ZERO, not downward, so the fade lifts the sets Myc had
# REPRESSED. Computed here rather than quoted, because it is the one prediction of
# region A that a reader can check against the whole library.
fam_tbl <- fg |>
  dplyr::filter(ranking %in% c("myc_6W", "timepoint_pos")) |>
  dplyr::select(ranking, pathway, NES) |>
  tidyr::pivot_wider(names_from = ranking, values_from = NES) |>
  dplyr::filter(!is.na(myc_6W), !is.na(timepoint_pos)) |>
  dplyr::mutate(fam = ifelse(myc_6W > 0, "induced", "repressed")) |>
  dplyr::group_by(fam) |>
  dplyr::summarise(n = dplyr::n(), med_tl = stats::median(timepoint_pos),
                   .groups = "drop")
med_tl <- function(f) fam_tbl$med_tl[fam_tbl$fam == f]
n_fam  <- function(f) fam_tbl$n[fam_tbl$fam == f]
stopifnot(nrow(fam_tbl) == 2L, sum(fam_tbl$n) == 866L,
          med_tl("induced") < 0, med_tl("repressed") > 0)

# Two sets, chosen because they differ in the one way the toy cannot show. For the
# MYC targets the Myc+ timeline is much more negative than the wild-type one, so
# the fade dominates -- the toy's case. For the OXPHOS subunits the two timelines
# are the SAME, so what moves them is the gland's own trajectory, not the fade.
REAL <- c("HALLMARK_MYC_TARGETS_V1"   = "Hallmark MYC targets V1",
          "MITOCARTA_OXPHOS_SUBUNITS" = "MitoCarta OXPHOS subunits")

real <- as.data.frame(fg[fg$pathway %in% names(REAL) &
                           fg$ranking %in% names(RANK_MAP),
                         c("pathway", "ranking", "NES")])
real$contrast <- factor(unname(RANK_MAP[real$ranking]), levels = contrast_levels)
real$set      <- factor(unname(REAL[real$pathway]), levels = unname(REAL))
stopifnot(nrow(real) == 8L, !anyNA(real$contrast), !anyNA(real$set))

nes_of <- function(s, c) real$NES[real$set == unname(REAL[s]) & real$contrast == c]
stopifnot(
  nes_of("HALLMARK_MYC_TARGETS_V1", "myc_6W")    > 2,
  nes_of("HALLMARK_MYC_TARGETS_V1", "myc_12W")   > 2,
  nes_of("HALLMARK_MYC_TARGETS_V1", "6>12W_myc") < -2,
  # the fade dominates for MYC targets: the Myc+ timeline is the more negative
  nes_of("HALLMARK_MYC_TARGETS_V1", "6>12W_myc") <
    nes_of("HALLMARK_MYC_TARGETS_V1", "6>12W_wt") - 0.5,
  # and for OXPHOS the two timelines are within a quarter of a unit of each other
  abs(nes_of("MITOCARTA_OXPHOS_SUBUNITS", "6>12W_myc") -
      nes_of("MITOCARTA_OXPHOS_SUBUNITS", "6>12W_wt")) < 0.25)

pD <- ggplot2::ggplot(real, ggplot2::aes(contrast, NES, fill = contrast)) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey55") +
  ggplot2::geom_col(width = 0.66) +
  ggplot2::geom_text(ggplot2::aes(label = sprintf("%+.2f", NES),
                                  vjust = ifelse(NES > 0, -0.35, 1.25)),
                     size = 1.85, colour = "grey15") +
  ggplot2::facet_wrap(~ set, nrow = 1) +
  ggplot2::scale_fill_manual(values = contrast_cols, guide = "none") +
  ggplot2::scale_y_continuous(labels = lab_signed,
                              expand = ggplot2::expansion(mult = c(0.17, 0.17))) +
  ggplot2::labs(x = NULL, y = "NES, real data") +
  theme_panel() +
  ggplot2::theme(strip.text = ggplot2::element_text(face = "bold", size = 6.2),
                 axis.text.x = ggplot2::element_text(angle = 40, hjust = 1),
                 panel.grid.major.x = ggplot2::element_blank(),
                 panel.spacing.x = ggplot2::unit(1.6, "mm"))

p <- patchwork::wrap_plots(pA, pB, pC, pD, ncol = 1,
                           heights = c(0.86, 1.30, 0.86, 0.98))

LEGEND <- panel_legend(
  slot = "not currently cited -- teaching figure",
  what = paste0(
    "Why a Myc effect that keeps its ranking gives a NEGATIVE enrichment score ",
    "on the Myc+ timeline. A to C are a simulation; D is the real data on the ",
    "same four contrasts."),
  detail = c(
    sprintf("THE SIMULATION: %d genes, of which %d are a 'Myc-target' set induced at six weeks. The twelve-week Myc effect is EXACTLY %.2f times the six-week one, gene for gene, so those two lists are in identical order by construction. The wild-type gland has a trajectory of its own, independent of the Myc effect and of set membership. All four correlations below are computed, and the three structural ones are asserted in the script rather than claimed here.",
            N_GENE, N_SET, FADE),
    "A, THE PREMISE: each line joins one gene's Myc effect at the two ages. No line crosses another -- that is the preserved ranking. Every line moves toward zero -- that is the fade. Ranks cannot show amplitude, which is the only reason this region exists; everything else on the panel is about order.",
    sprintf("TOWARD ZERO, NOT DOWNWARD, and the difference is testable. Halving compresses both directions, so the genes Myc had pushed DOWN come back UP -- the pale lines below zero rise. The library says exactly that: over the %d sets Myc induces the median Myc+ timeline NES is %+.2f, and over the %d sets Myc represses it is %+.2f. A fade produces a SIGN FLIP relative to the Myc effect, not a uniform decline.",
            n_fam("induced"), med_tl("induced"),
            n_fam("repressed"), med_tl("repressed")),
    sprintf("B, THE ANSWER TO THE QUESTION: the gene order of each list against the gene order of myc_6W, all %d genes drawn, rank 1 being the top of a list on both axes. myc_12W is a PERFECT DIAGONAL (rho %+.2f) -- the same list, which is what 'the ranking is preserved' means and the whole of what it means. 6>12W_wt is a CLOUD (rho %+.2f) -- a real, structured list in which this set has no business being anywhere in particular. 6>12W_myc is an ANTI-DIAGONAL (rho %+.2f) -- the reverse of myc_6W, because it is built from what the Myc effect LOST. FOLLOW THE ORANGE: the set sits at the TOP of myc_6W and myc_12W, nowhere in particular in 6>12W_wt, and at the BOTTOM of 6>12W_myc. Same genes, four lists.",
            N_GENE, rho_df$rho[1], rho_df$rho[2], rho_df$rho[3]),
    sprintf("C, WHAT fGSEA DOES WITH THOSE FOUR ORDERS: the same set through fgsea::fgsea() once per list, drawn as the running enrichment score with set members ticked. myc_6W %+.2f and myc_12W %+.2f -- the same ranking, so the same score, differing only by permutation noise, and the halved amplitude is invisible because NES is scale-free. 6>12W_wt %+.2f, nothing. 6>12W_myc %+.2f. Same genes, same set, four lists.",
            nes_at("myc_6W"), nes_at("myc_12W"),
            nes_at("6>12W_wt"), nes_at("6>12W_myc")),
    "THE ONE SENTENCE: fGSEA never compares two rankings. It takes one ranked list and asks where a set sits in it, and there are four lists here. 'The ranking is preserved' is about the two genotype lists; the Myc+ timeline is built from the DIFFERENCE between them, which is the part that was lost. Halving is a change: a target 2-fold over wild type at six weeks and 1.4-fold at twelve went DOWN inside the Myc+ animals, and so did every other target, in the same order.",
    "SO THE NEGATIVE TIMELINE IS NOT EVIDENCE AGAINST STABILITY -- IT IS WHAT STABILITY PRODUCES. A rank-unstable fade would give a weak, incoherent timeline score; a clean mirror image is what a preserved shape looks like from the other side. In the real data Spearman(myc_6W, interaction) = -0.878 over 866 sets and Spearman(myc_6W, 6>12W_myc) = -0.838, against -0.068 for the wild-type timeline.",
    sprintf("D, THE REAL NES for two named sets on all four contrasts, from script 20. Hallmark MYC targets V1 is the toy's case -- the Myc+ timeline (%+.2f) is far more negative than the wild-type one (%+.2f), so the fade dominates. The OXPHOS subunits are NOT: the two timelines are the same (%+.2f against %+.2f), so what moves them is the gland's own trajectory, not the fade. That difference is the substrate frame, and it is the one thing the simulation is deliberately built too simply to show.",
            nes_of("HALLMARK_MYC_TARGETS_V1", "6>12W_myc"),
            nes_of("HALLMARK_MYC_TARGETS_V1", "6>12W_wt"),
            nes_of("MITOCARTA_OXPHOS_SUBUNITS", "6>12W_myc"),
            nes_of("MITOCARTA_OXPHOS_SUBUNITS", "6>12W_wt"))),
  bounds = c(
    "REGIONS A TO C ARE SIMULATED and are the only simulated drawing in this layer. They carry no result and support no claim about the mouse; they show what the arithmetic does. Only region D is data.",
    "THE SIMULATED WILD-TYPE TIMELINE IS INDEPENDENT OF THE MYC EFFECT, WHICH IS THE MEASURED SITUATION BUT NOT THE WHOLE OF IT. Spearman(myc_6W, 6>12W_wt) = -0.068 over the 866 library sets, so 'unrelated' is right on average -- but that near-zero is an average over two opposite-signed families (Myc-induced sets -0.314, Myc-repressed +0.345, both FALLING). The real wild-type timeline is strongly structured; its structure simply sits on other genes. The simulation draws the average, not the structure.",
    "AND THE REAL 6>12W_myc IS 6>12W_wt PLUS THE FADE, so for some arms -- the OXPHOS subunits in region D -- the wild-type half is the whole of it. Never read a negative Myc+ timeline as the fade without checking the wild-type timeline beside it. The simulation's Myc+ timeline is likewise the mirror plus development, which is why its rho is well short of -1 while the INTERACTION's is exactly -1.",
    "The constant-standard-error assumption is the measured situation and not a convenience: over the 1,967 genes Myc moves at six weeks the median lfcSE is 0.189 at six weeks against 0.193 at twelve, so the amplitude change is in the numerator alone.",
    "A TEMPORAL NES AND A GENOTYPE NES ARE NOT ON THE SAME RULER, which is why region D's four bars must be read for sign and ordering rather than as comparable magnitudes. fGSEA normalises each score against a permutation null built from that same ranked list, so the normaliser differs between the four.",
    "The lines drawn in REGION A are a stratified sample of 26 of the 2,000 genes, because 2,000 lines at that height is a filled rectangle. Region B draws all 2,000, and every correlation quoted is over all 2,000.",
    "The two timeline contrasts in region D span the two extraction batches (batch = timepoint), so they are descriptive. The two genotype contrasts are clean. Nothing on this figure rests on the temporal bars being a clean effect -- they are drawn to show a SIGN."),
  source = c(
    "REGIONS A to C: simulated in this script, set.seed(42); the enrichment scores are real fgsea::fgsea() calls on the simulated lists",
    "REGION D: results/fgsea_percategory.rds (scripts/20_fgsea_percategory.R) -- NES per ranking x category on the unshrunken Wald statistic",
    "The lfcSE figures quoted in the bounds: results/interaction_results.rds, myc_6W_raw and myc_12W_raw",
    "The 866-set correlations quoted in the detail and bounds: results/fgsea_percategory.rds, rankings myc_6W, interaction, timepoint_neg and timepoint_pos"))

save_panel_p(p, "explainer_nes_sign", width = fig_w[["single"]], height = 152)

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {

  print(p)
  print(LEGEND)

  ## the four simulated scores, and the rank correlations behind them
  print(sim_nes, row.names = FALSE, digits = 3)
  vapply(names(LISTS), function(nm) stats::cor(e6, LISTS[[nm]], method = "spearman"),
         numeric(1)) |> round(3) |> print()
  cat("Spearman(myc_6W, interaction) =",
      round(stats::cor(e6, inter, method = "spearman"), 3), "\n")

  ## the sign does not depend on the fade being 0.5, only on it being below 1.
  ## A fade of 1 (no change) leaves the Myc+ timeline equal to the wild-type one;
  ## a fade above 1 (a GROWING effect) flips the timeline positive.
  for (f in c(0.25, 0.5, 0.75, 1.0, 1.5)) {
    r <- sort(dev + (f * e6 - e6), decreasing = TRUE)
    set.seed(1)
    ff <- fgsea::fgsea(list(S = set_ids), r, minSize = 10, maxSize = 500,
                       eps = 0, nPermSimple = 10000)
    cat(sprintf("fade %.2f  Myc+ timeline NES %+.2f\n", f, ff$NES[1]))
  }

  ## the real 866-set correlations the legend quotes
  fgw <- fg |>
    dplyr::filter(ranking %in% c("myc_6W", "myc_12W", "timepoint_neg",
                                 "timepoint_pos", "interaction")) |>
    dplyr::select(ranking, pathway, NES) |>
    tidyr::pivot_wider(names_from = ranking, values_from = NES)
  vapply(c("myc_12W", "timepoint_neg", "timepoint_pos", "interaction"),
         function(k) stats::cor(fgw$myc_6W, fgw[[k]], method = "spearman"),
         numeric(1)) |> round(3) |> print()
}
