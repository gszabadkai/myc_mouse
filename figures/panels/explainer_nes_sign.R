# =============================================================================
# explainer_nes_sign.R -- why a rank-STABLE Myc effect gives a NEGATIVE NES on
# the Myc+ timeline
# -----------------------------------------------------------------------------
# SLOT: none, and deliberately so. The filename does not start with "fig", so it
# falls outside the `^fig.*\.R$` glob that rebuild_panels.R:44 and
# panels_to_pdf.R:65 both use. It therefore does not enter the 28-panel count,
# does not need a slot letter, and does not need a chapter in
# paper/analysis_record.qmd. It is a TEACHING figure, built 2026-08-16 to answer
# the author's question, and it is the only figure in this layer whose left two
# thirds are SIMULATED. That is marked on the page, not only here.
#
# THE QUESTION. Fig. 1D/1H say the Myc programme keeps its ranking between six
# and twelve weeks (Spearman 0.933 over 866 sets) while the effect size halves
# (Fig. 1G). scripts/04 -> outputs/fgsea/ says most pathways have a strongly
# NEGATIVE NES on the two timelines. If the ranking does not move, why does the
# timeline move at all?
#
# THE ANSWER, and it is one sentence: fGSEA never compares two rankings. It takes
# ONE ranked gene list and asks where a set sits in it -- and there are THREE
# lists here, not two.
#
#     list 1   the Myc effect at 6W       Myc+ minus WT, at six weeks
#     list 2   the Myc effect at 12W      Myc+ minus WT, at twelve weeks
#     list 3   the Myc+ timeline          12W minus 6W, WITHIN Myc+ animals
#
# "The ranking is preserved" is a statement about lists 1 and 2. List 3 is built
# from the DIFFERENCE between them, and a difference of two proportional things is
# not zero -- it is the part that was lost. Halving is a change. A Myc target that
# was 2-fold over wild type at six weeks and 1.4-fold at twelve went DOWN inside
# the Myc+ animals, and so did every other Myc target, in the same order. List 3
# is therefore the NEGATIVE of list 1, and the set that topped lists 1 and 2 sits
# at the BOTTOM of list 3.
#
# So the negative timeline NES is not evidence against stability. It is what
# stability PRODUCES. A rank-UNSTABLE fade would give a weak, incoherent timeline
# NES; that the timeline is a clean mirror is further evidence the shape held.
#
# WHAT THE PANEL DRAWS
#   TOP     the simulation, as a slopegraph. Every gene's Myc effect at each age,
#           halved exactly and with no reordering. The LINE connecting the two
#           ages IS the timeline contrast, and every line descends.
#   MIDDLE  the same simulation put through real fgsea::fgsea() three times, one
#           per list, drawn as the running enrichment score. Same set, same genes,
#           two positive scores and one negative one.
#   BOTTOM  the real data, for two named sets, on all four contrasts -- so the toy
#           can be checked against the thing it claims to explain, and so the one
#           way the real data DIFFERS from the toy is visible: the wild-type gland
#           moves too, and for OXPHOS it moves as much as the fade does.
#
# Reads (read-only, no re-run):
#   results/fgsea_percategory.rds (script 20) -- $fgsea, the BOTTOM region only
# Output: outputs/figures/panels/explainer_nes_sign.pdf
# =============================================================================

source(here::here("figures", "panels", "_panel_common.R"))

stopifnot(requireNamespace("fgsea", quietly = TRUE),
          requireNamespace("patchwork", quietly = TRUE))

# =============================================================================
# THE SIMULATION
# =============================================================================
# Deliberately the simplest world in which the question can be asked:
#   * the wild-type gland does NOT change with age, so the Myc+ timeline is the
#     fade and nothing else (the real gland does change -- that is the BOTTOM
#     region's job, and the difference is the point of drawing it);
#   * the twelve-week effect is EXACTLY half the six-week one, gene for gene, so
#     the rank correlation between lists 1 and 2 is 1 by construction and no
#     reordering can be smuggled in;
#   * the ranking statistic is the effect itself, i.e. constant standard error.
#     That is the measured situation and not a convenience: over the 1,967 genes
#     Myc moves at six weeks the median lfcSE is 0.189 at six weeks and 0.193 at
#     twelve, so the fade is in the numerator alone.

set.seed(42)

N_GENE  <- 2000L
N_SET   <- 120L
FADE    <- 0.5                      # the amplitude that survives to twelve weeks

genes   <- sprintf("g%04d", seq_len(N_GENE))
set_ids <- genes[seq_len(N_SET)]    # the first N_SET genes are the "Myc target" set

# The Myc effect at six weeks: the set is induced, the background is not.
e6 <- stats::rnorm(N_GENE, mean = 0, sd = 0.40)
e6[seq_len(N_SET)] <- stats::rnorm(N_SET, mean = 1.20, sd = 0.40)
names(e6) <- genes

# The Myc effect at twelve weeks: the SAME numbers, times a constant.
e12 <- FADE * e6

# The Myc+ timeline: what happened inside the Myc+ animals between the two ages.
# With a static wild type this is exactly e12 - e6 = -(1 - FADE) * e6.
tl <- e12 - e6

# The two identities the whole explanation rests on, asserted rather than claimed.
stopifnot(
  # lists 1 and 2 are in the SAME order -- no reordering anywhere
  isTRUE(all.equal(stats::cor(e6, e12, method = "spearman"), 1)),
  # list 3 is the exact negative of list 1, up to the scale factor
  isTRUE(all.equal(stats::cor(e6, tl, method = "spearman"), -1)),
  isTRUE(all.equal(unname(tl), unname(-(1 - FADE) * e6))))

# =============================================================================
# THREE REAL fGSEA RUNS -- the NES on the panel is computed, never asserted
# =============================================================================
LISTS <- list("myc_6W" = e6, "myc_12W" = e12, "6>12W_myc" = tl)

set.seed(42)
sim_nes <- do.call(rbind, lapply(names(LISTS), function(nm) {
  r <- sort(LISTS[[nm]], decreasing = TRUE)
  f <- fgsea::fgsea(pathways = list(MYC_TARGET_SET = set_ids), stats = r,
                    minSize = 10, maxSize = 500, eps = 0, nPermSimple = 10000)
  data.frame(list_name = nm, NES = f$NES[1], padj = f$padj[1],
             stringsAsFactors = FALSE)
}))
sim_nes$list_name <- factor(sim_nes$list_name, levels = names(LISTS))

# The result the panel exists to show: the two genotype lists agree (they are the
# same ranking), and the timeline flips sign.
stopifnot(sim_nes$NES[1] > 2, sim_nes$NES[2] > 2, sim_nes$NES[3] < -2,
          abs(sim_nes$NES[1] - sim_nes$NES[2]) < 0.15)

# --- the running enrichment score, the classic weighted-KS walk ---------------
# Implemented here rather than taken from fgsea::plotEnrichment() so the three
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
es_df$list_name <- factor(es_df$list_name, levels = names(LISTS))

# =============================================================================
# REGION 1 -- the slopegraph: same order, half the size, every line descending
# =============================================================================
# A SAMPLE of genes is drawn, because 2,000 lines at 34 mm is a black rectangle.
# The sample is stratified so the set's whole range is represented; the rank
# correlation quoted in the legend is over all 2,000.
set.seed(7)
show <- c(set_ids[round(seq(1, N_SET, length.out = 11))],
          sample(genes[(N_SET + 1L):N_GENE], 15L))

slope_df <- data.frame(
  gene = rep(show, 2L),
  age  = factor(rep(c("myc_6W", "myc_12W"), each = length(show)),
                levels = c("myc_6W", "myc_12W")),
  lfc  = c(e6[show], e12[show]),
  fam  = rep(ifelse(show %in% set_ids,
                    "Myc-target set", "background"), 2L),
  stringsAsFactors = FALSE)
slope_df$fam <- factor(slope_df$fam, levels = c("Myc-target set", "background"))

FAM_COLS <- c("Myc-target set" = unname(contrast_cols[["myc_6W"]]),
              "background"     = "grey78")

p1 <- ggplot2::ggplot(slope_df,
                      ggplot2::aes(age, lfc, group = gene, colour = fam)) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey55") +
  ggplot2::geom_line(linewidth = 0.34) +
  ggplot2::geom_point(size = 0.62) +
  ggplot2::scale_colour_manual(values = FAM_COLS, name = NULL) +
  ggplot2::scale_x_discrete(expand = ggplot2::expansion(add = c(0.30, 0.62))) +
  ggplot2::scale_y_continuous(labels = lab_signed) +
  # The line IS the timeline contrast. Naming it on the page is the whole point
  # of the region, and it is a NAME of a drawn element, not prose.
  ggplot2::annotate("text", x = 2.46, y = 0.92,
                    label = "each line =\n6>12W_myc", hjust = 0.5, vjust = 0.5,
                    size = 1.95, lineheight = 0.95, colour = "grey20") +
  ggplot2::annotate("segment", x = 2.26, xend = 2.055, y = 0.86, yend = 0.64,
                    linewidth = 0.22, colour = "grey40",
                    arrow = grid::arrow(length = grid::unit(0.9, "mm"),
                                        type = "closed")) +
  ggplot2::labs(x = NULL, y = "log2 fold change vs wild type") +
  theme_panel() +
  ggplot2::theme(legend.position = "top",
                 legend.margin = ggplot2::margin(0, 0, -2, 0),
                 axis.text.x = ggplot2::element_text(face = "bold"),
                 panel.grid.major.x = ggplot2::element_blank())

# =============================================================================
# REGION 2 -- the same set through fgsea three times
# =============================================================================
TICK_LO <- -0.09   # where the member ticks sit, below the curve
TICK_HI <- -0.02

tick_df <- es_df[es_df$hit, c("rank", "list_name")]

nes_lab <- sim_nes
nes_lab$lab  <- sprintf("NES %+.2f", nes_lab$NES)
nes_lab$xpos <- N_GENE * 0.97
nes_lab$ypos <- 0.72
# NEUTRAL INK, deliberately. The strip names the list and the curve carries its
# declared colour, so the number is unambiguous already; spending a third encoding
# on it would either repeat the strip or invent a sign colour the curve's own
# direction already shows. It also keeps 6>12W_myc's light grey off 2 mm type.
nes_lab$col  <- "grey15"

p2 <- ggplot2::ggplot(es_df, ggplot2::aes(rank, es)) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey55") +
  ggplot2::geom_line(ggplot2::aes(colour = list_name), linewidth = 0.45) +
  ggplot2::geom_segment(data = tick_df,
                        ggplot2::aes(x = rank, xend = rank,
                                     y = TICK_LO, yend = TICK_HI),
                        inherit.aes = FALSE, linewidth = 0.16, colour = "grey30") +
  ggplot2::geom_text(data = nes_lab,
                     ggplot2::aes(xpos, ypos, label = lab, colour = NULL),
                     inherit.aes = FALSE, hjust = 1, vjust = 1, size = 2.1,
                     fontface = "bold", colour = nes_lab$col) +
  # Name the tick row once, in the first facet, rather than leaving the reader to
  # infer it from the GSEA idiom.
  ggplot2::geom_text(data = data.frame(list_name = factor("myc_6W",
                                                          levels = names(LISTS))),
                     ggplot2::aes(x = 60, y = -0.24, label = "set members"),
                     inherit.aes = FALSE, hjust = 0, size = 1.85,
                     colour = "grey30") +
  ggplot2::facet_wrap(~ list_name, nrow = 1) +
  ggplot2::scale_colour_manual(values = contrast_cols, guide = "none") +
  ggplot2::scale_x_continuous(breaks = c(1, 2000), labels = c("1", "2000")) +
  ggplot2::scale_y_continuous(labels = lab_signed) +
  ggplot2::labs(x = "gene rank in that list (high to low)",
                y = "running enrichment") +
  theme_panel() +
  ggplot2::theme(strip.text = ggplot2::element_text(face = "bold", size = 6.4),
                 panel.spacing.x = ggplot2::unit(2.8, "mm"))

# =============================================================================
# REGION 3 -- the real data, so the toy can be checked against it
# =============================================================================
fg_path <- here::here("results", "fgsea_percategory.rds")
fg <- readRDS(fg_path)$fgsea

# Two sets, chosen because they differ in the one way the toy cannot show. For the
# MYC targets the Myc+ timeline is much more negative than the wild-type one, so
# the fade dominates -- the toy's case. For the OXPHOS subunits the two timelines
# are the SAME, so what moves them is the gland's own trajectory, not the fade.
REAL <- c("HALLMARK_MYC_TARGETS_V1"   = "Hallmark MYC targets V1",
          "MITOCARTA_OXPHOS_SUBUNITS" = "MitoCarta OXPHOS subunits")
RANK_MAP <- c(myc_6W = "myc_6W", myc_12W = "myc_12W",
              timepoint_neg = "6>12W_wt", timepoint_pos = "6>12W_myc")

# The compression is TOWARD ZERO, not downward, so the fade lifts the sets Myc had
# REPRESSED. Computed here rather than quoted, because it is the one prediction of
# the top region that a reader can check against the whole library.
fam <- fg |>
  dplyr::filter(ranking %in% c("myc_6W", "timepoint_pos")) |>
  dplyr::select(ranking, pathway, NES) |>
  tidyr::pivot_wider(names_from = ranking, values_from = NES) |>
  dplyr::filter(!is.na(myc_6W), !is.na(timepoint_pos)) |>
  dplyr::mutate(fam = ifelse(myc_6W > 0, "induced", "repressed")) |>
  dplyr::group_by(fam) |>
  dplyr::summarise(n = dplyr::n(), med_tl = stats::median(timepoint_pos),
                   .groups = "drop")
med_tl <- function(f) fam$med_tl[fam$fam == f]
n_fam  <- function(f) fam$n[fam$fam == f]
stopifnot(nrow(fam) == 2L, sum(fam$n) == 866L,
          med_tl("induced") < 0, med_tl("repressed") > 0)

real <- fg[fg$pathway %in% names(REAL) & fg$ranking %in% names(RANK_MAP),
           c("pathway", "ranking", "NES")]
real <- as.data.frame(real)
real$contrast <- factor(unname(RANK_MAP[real$ranking]), levels = contrast_levels)
real$set      <- factor(unname(REAL[real$pathway]), levels = unname(REAL))
stopifnot(nrow(real) == 8L, !anyNA(real$contrast), !anyNA(real$set))

# The two facts the region is drawn to make visible, asserted before drawing.
nes_of <- function(s, c) real$NES[real$set == unname(REAL[s]) & real$contrast == c]
stopifnot(
  nes_of("HALLMARK_MYC_TARGETS_V1", "myc_6W")  > 2,
  nes_of("HALLMARK_MYC_TARGETS_V1", "myc_12W") > 2,
  nes_of("HALLMARK_MYC_TARGETS_V1", "6>12W_myc") < -2,
  # the fade dominates for MYC targets: the Myc+ timeline is the more negative
  nes_of("HALLMARK_MYC_TARGETS_V1", "6>12W_myc") <
    nes_of("HALLMARK_MYC_TARGETS_V1", "6>12W_wt") - 0.5,
  # and for OXPHOS the two timelines are within a quarter of a unit of each other
  abs(nes_of("MITOCARTA_OXPHOS_SUBUNITS", "6>12W_myc") -
      nes_of("MITOCARTA_OXPHOS_SUBUNITS", "6>12W_wt")) < 0.25)

p3 <- ggplot2::ggplot(real, ggplot2::aes(contrast, NES, fill = contrast)) +
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
  ggplot2::theme(strip.text = ggplot2::element_text(face = "bold", size = 6.4),
                 axis.text.x = ggplot2::element_text(angle = 40, hjust = 1),
                 panel.grid.major.x = ggplot2::element_blank(),
                 panel.spacing.x = ggplot2::unit(1.6, "mm"))

p <- patchwork::wrap_plots(p1, p2, p3, ncol = 1, heights = c(1.05, 1, 1.05))

LEGEND <- panel_legend(
  slot = "not currently cited -- teaching figure",
  what = paste0(
    "Why a Myc effect that keeps its ranking gives a NEGATIVE enrichment score ",
    "on the Myc+ timeline. TOP and MIDDLE are a simulation; BOTTOM is the real ",
    "data on the same four contrasts."),
  detail = c(
    sprintf("THE SIMULATION: %d genes, of which %d are a 'Myc-target' set induced at six weeks. The twelve-week effect is EXACTLY %.2f times the six-week effect, gene for gene, and the wild-type gland does not change. Spearman between the six- and twelve-week lists is therefore 1 by construction, and between the six-week list and the timeline it is -1; both are asserted in the script, not claimed here.",
            N_GENE, N_SET, FADE),
    "TOP: each line joins one gene's Myc effect at the two ages. No line crosses another -- that is the preserved ranking. Every line moves toward zero -- that is the fade. The two statements are about the same picture, and THE LINE ITSELF IS THE TIMELINE CONTRAST, because with a static wild type 6>12W_myc = myc_12W - myc_6W.",
    sprintf("TOWARD ZERO, NOT DOWNWARD, and the difference is testable. Halving compresses both directions, so the genes Myc had pushed DOWN come back UP -- the pale lines below zero rise. The library says exactly that: over the %d sets Myc induces the median Myc+ timeline NES is %+.2f, and over the %d sets Myc represses it is %+.2f. A fade produces a SIGN FLIP relative to the Myc effect, not a uniform decline.",
            n_fam("induced"), med_tl("induced"),
            n_fam("repressed"), med_tl("repressed")),
    sprintf("MIDDLE: the same set through fgsea::fgsea() three times, once per ranked list, drawn as the running enrichment score with set members ticked. The two genotype lists give %+.2f and %+.2f -- the same ranking, so the same score, differing only by permutation noise. The timeline gives %+.2f. The set that sits at the LEFT of the first two lists sits at the RIGHT of the third.",
            sim_nes$NES[1], sim_nes$NES[2], sim_nes$NES[3]),
    "THE ONE SENTENCE: fGSEA never compares two rankings. It takes one ranked list and asks where a set sits in it, and there are three lists here. 'The ranking is preserved' is about lists 1 and 2; list 3 is built from the DIFFERENCE between them, which is the part that was lost. Halving is a change: a target 2-fold over wild type at six weeks and 1.4-fold at twelve went DOWN inside the Myc+ animals, and so did every other target, in the same order.",
    "SO THE NEGATIVE TIMELINE IS NOT EVIDENCE AGAINST STABILITY -- IT IS WHAT STABILITY PRODUCES. A rank-unstable fade would give a weak, incoherent timeline score; a clean mirror image is what a preserved shape looks like from the other side. In the real data Spearman(myc_6W, interaction) = -0.878 over 866 sets and Spearman(myc_6W, 6>12W_myc) = -0.838.",
    sprintf("BOTTOM: the real NES for two named sets on all four contrasts, from script 20. Hallmark MYC targets V1 is the toy's case -- the Myc+ timeline (%+.2f) is far more negative than the wild-type one (%+.2f), so the fade dominates. The OXPHOS subunits are NOT: the two timelines are the same (%+.2f against %+.2f), so what moves them is the gland's own trajectory, not the fade. That difference is the substrate frame, and it is the one thing the simulation is built too simply to show.",
            nes_of("HALLMARK_MYC_TARGETS_V1", "6>12W_myc"),
            nes_of("HALLMARK_MYC_TARGETS_V1", "6>12W_wt"),
            nes_of("MITOCARTA_OXPHOS_SUBUNITS", "6>12W_myc"),
            nes_of("MITOCARTA_OXPHOS_SUBUNITS", "6>12W_wt"))),
  bounds = c(
    "THE TOP TWO REGIONS ARE SIMULATED and are the only simulated drawing in this layer. They carry no result and support no claim about the mouse; they show what the arithmetic does. Only the bottom region is data.",
    "THE SIMULATION IS DELIBERATELY SIMPLER THAN THE GLAND in one way that matters: its wild type does not change with age. The real 6>12W_myc is 6>12W_wt PLUS the fade, and for some arms -- the OXPHOS subunits above -- the wild-type half is the whole of it. Never read a negative Myc+ timeline as the fade without checking the wild-type timeline beside it.",
    "The constant-standard-error assumption is the measured situation and not a convenience: over the 1,967 genes Myc moves at six weeks the median lfcSE is 0.189 at six weeks against 0.193 at twelve, so the amplitude change is in the numerator alone.",
    "A TEMPORAL NES AND A GENOTYPE NES ARE NOT ON THE SAME RULER, which is why the bottom region's four bars must be read for sign and ordering rather than as comparable magnitudes. fGSEA normalises each score against a permutation null built from that same ranked list, so the normaliser differs between the four.",
    "The lines drawn in the TOP region are a stratified sample of 26 of the 2,000 genes, because 2,000 lines at this height is a filled rectangle. The rank correlations quoted are over all 2,000.",
    "The two timeline contrasts in the bottom region span the two extraction batches (batch = timepoint), so they are descriptive. The two genotype contrasts are clean. Nothing on this figure rests on the temporal bars being a clean effect -- they are drawn to show a SIGN."),
  source = c(
    "TOP and MIDDLE: simulated in this script, set.seed(42); the enrichment scores are real fgsea::fgsea() calls on the simulated lists",
    "BOTTOM: results/fgsea_percategory.rds (scripts/20_fgsea_percategory.R) -- NES per ranking x category on the unshrunken Wald statistic",
    "The lfcSE figures quoted in the bounds: results/interaction_results.rds, myc_6W_raw and myc_12W_raw",
    "The 866-set correlations quoted in the detail: results/fgsea_percategory.rds, rankings myc_6W, interaction and timepoint_pos"))

save_panel_p(p, "explainer_nes_sign", width = fig_w[["single"]], height = 112)

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {

  print(p)
  print(LEGEND)

  ## the three simulated scores, and the two identities behind them
  print(sim_nes, row.names = FALSE, digits = 3)
  cat("Spearman(myc_6W, myc_12W)   =", stats::cor(e6, e12, method = "spearman"), "\n")
  cat("Spearman(myc_6W, timeline)  =", stats::cor(e6, tl,  method = "spearman"), "\n")

  ## the same three scores at other fade factors -- the sign does not depend on 0.5,
  ## only on the fade being less than 1. A fade of 1 (no change) gives no timeline
  ## at all, and a fade above 1 (a growing effect) flips the timeline positive.
  for (f in c(0.25, 0.5, 0.75, 1.0, 1.5)) {
    t2 <- f * e6 - e6
    if (all(t2 == 0)) { cat(sprintf("fade %.2f  timeline is identically zero\n", f)); next }
    r  <- sort(t2, decreasing = TRUE)
    set.seed(1)
    ff <- fgsea::fgsea(list(S = set_ids), r, minSize = 10, maxSize = 500,
                       eps = 0, nPermSimple = 10000)
    cat(sprintf("fade %.2f  timeline NES %+.2f\n", f, ff$NES[1]))
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
