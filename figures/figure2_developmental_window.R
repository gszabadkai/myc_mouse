# =============================================================================
# figure2_developmental_window.R -- MANUSCRIPT FIGURE 2
# "Development opens the window"
# -----------------------------------------------------------------------------
# Beats three to five of the Results paragraph, then the hand-off to the cells.
#
#   A  THE MATURING GLAND WITHDRAWS FROM RESPIRATION, AND FROM NOTHING ELSE THAT
#      WOULD EXPLAIN IT. Set-average change across the wild-type window, each arm
#      against its own expression-matched null. The ASSEMBLY FACTORS OF THE SAME
#      COMPLEXES do not move: that is the internal control that makes this
#      specific to the respiratory subunits rather than to mitochondria at large.
#   B  AND THE CLAIM IS COMPARATIVE, SO IT IS TESTED AS A COMPARISON. Proliferation
#      is small but not immobile, so a one-set null cannot carry "de-respires
#      without de-proliferating". Redrawing both sets together tests the DIFFERENCE.
#   C  IT IS NOT AN ARTEFACT OF THE CONTENT RISE. The withdrawal persists on
#      mitoPPS, a within-compartment ratio that is blind to how many mitochondria
#      a cell has -- so it cannot be produced by scaling the compartment.
#   D  THE GLAND DOES NOT DISMANTLE THE DEATH MACHINERY, AND DOES NOT BUFFER IT.
#      One transcript of 37 moves, and it moves UP. DEVELOPMENT REMOVES THE INPUT,
#      NOT THE MACHINE.
#   E  SO WHAT MYC LOSES IS THE TRIGGER. Every arm of the MYC programme is rescaled
#      by ~0.55 over the window; BAX:BCL-xL retains exactly that, at the 50th
#      percentile of matched pairs. PUMA:BCL-xL does not fade -- it REVERSES.
#   F  WHAT THE MOUSE HANDS TO THE CELL EXPERIMENTS. Three propositions the same
#      data licenses, each with its control and the experiment that adjudicates it.
#
# EPISTEMIC STATUS OF PANEL F, on its face and not only in the legend: these are
# RANKINGS, not tests. Nothing beats its permutation null; the respiratory axis is
# collinear with the timepoint it would explain (the term goes p 0.005 -> 0.088 once
# tp*myc is in the model); and H2 is not separable from H1 at n=24. The mouse
# generates; the perturbations prove.
#
# BATCH = TIMEPOINT: every wild-type temporal value in A-D is DESCRIBED, not
# claimed. Three things keep it usable: the comparison is internal to one contrast
# (a shared batch offset cannot make subunits fall while assembly factors of the
# same complexes do not), it reproduces on two rulers with different denominators,
# and the protein blots move the same way off the RNA batch entirely.
#
# Reads (read-only; the author runs scripts 42, 43 and 44 first):
#   results/substrate_specificity_tradeoff.rds -- $comparator/$wt_null (A),
#       $paired_null (B), $comparator_priority (C), $buffer (D), $tradeoff_perm (F)
#   results/collapse_module_ownership.rds      -- $wt_genes (D)
#   results/priming_arm_teb.rds                -- $machinery/$priming/$pair_null (E),
#       $perm_axis (F), $params$GLOBAL_RATE
# =============================================================================

source(here::here("figures", "theme_myc.R"))
if (!requireNamespace("patchwork", quietly = TRUE)) stop("figure2 needs patchwork")
if (!requireNamespace("ggrepel", quietly = TRUE))   stop("figure2 needs ggrepel")

out_dir <- here::here("outputs", "figures")

ss  <- readRDS(here::here("results", "substrate_specificity_tradeoff.rds"))
cmo <- readRDS(here::here("results", "collapse_module_ownership.rds"))
pa  <- readRDS(here::here("results", "priming_arm_teb.rds"))

need <- function(obj, fields, what) {
  miss <- fields[!fields %in% names(obj)]
  if (length(miss)) stop("figure2: ", what, " is missing -> ", paste(miss, collapse = ", "))
}
need(ss,  c("comparator", "wt_null", "paired_null", "buffer", "tradeoff_perm", "defs"),
     "substrate_specificity_tradeoff.rds")
need(cmo, "wt_genes", "collapse_module_ownership.rds")
need(pa,  c("machinery", "priming", "pair_null", "perm_axis", "params"),
     "priming_arm_teb.rds")

# --- guard: a stale object must fail here, not at review ----------------------
ox_obs <- ss$comparator$c_wt_time[ss$comparator$arm == "OXPHOS subunits"]
if (!isTRUE(abs(ox_obs - ss$defs$oxphos_wt_reference) < 0.02))
  stop(sprintf("figure2: stale rds -- OXPHOS-subunit wild-type %.4f vs reference %.4f",
               ox_obs, ss$defs$oxphos_wt_reference))

RATE <- pa$params$GLOBAL_RATE
stopifnot(is.finite(RATE), RATE > 0, RATE < 1)

verdict_col <- c("withdraws" = "#762A83", "at chance" = "grey55", "rises" = "#1B7837")
verdict_of  <- function(pct) ifelse(pct < 5, "withdraws",
                             ifelse(pct > 95, "rises", "at chance"))

# =============================================================================
# PANEL A -- the withdrawal, and its internal control
# =============================================================================
A <- merge(ss$comparator[, c("arm", "n_genes", "c_wt_time")],
           ss$wt_null[, c("arm", "percentile")], by = "arm")
A$verdict <- verdict_of(A$percentile)
A <- A[order(A$c_wt_time), ]
A$label <- factor(sprintf("%s  (%d)", A$arm, A$n_genes),
                  levels = rev(sprintf("%s  (%d)", A$arm, A$n_genes)))
xrA   <- range(c(A$c_wt_time, 0))
x_txA <- xrA[2] + diff(xrA) * 0.34

pA <- ggplot2::ggplot(A, ggplot2::aes(y = label)) +
  ggplot2::geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.35) +
  ggplot2::geom_segment(ggplot2::aes(x = 0, xend = c_wt_time, yend = label,
                                     colour = verdict), linewidth = 0.5) +
  ggplot2::geom_point(ggplot2::aes(x = c_wt_time, colour = verdict), size = 2.1) +
  ggplot2::geom_text(ggplot2::aes(x = x_txA, label = sprintf("%.1f", percentile)),
                     hjust = 1, size = 1.9, colour = "grey25") +
  ggplot2::annotate("text", x = x_txA, y = nrow(A) + 0.8, label = "null pct",
                    hjust = 1, size = 1.9, colour = "grey40", fontface = "italic") +
  ggplot2::scale_colour_manual(values = verdict_col, name = NULL) +
  ggplot2::coord_cartesian(xlim = c(xrA[1] - 0.03, x_txA + 0.01),
                           ylim = c(0.4, nrow(A) + 1.25), expand = FALSE) +
  ggplot2::labs(
    x = "wild-type 6 -> 12 weeks  (set-average log2FC)", y = NULL,
    title = "A   The gland withdraws from respiration",
    subtitle = "OXPHOS assembly = the same complexes as the subunits: the internal control.\nPercentile against 2000 expression-matched random sets.") +
  theme_myc(base_size = 8) +
  ggplot2::theme(legend.position = "bottom",
                 legend.key.size = ggplot2::unit(3, "mm"),
                 legend.text     = ggplot2::element_text(size = 6.2),
                 axis.text.y     = ggplot2::element_text(size = 6.4),
                 plot.title      = ggplot2::element_text(face = "bold", size = 8.5),
                 plot.subtitle   = ggplot2::element_text(size = 6.2, colour = "grey25",
                                                         lineheight = 1.15))

# =============================================================================
# PANEL B -- the comparison, tested as a comparison
# =============================================================================
B <- ss$paired_null
B <- B[order(B$observed_diff), ]
B$lab <- factor(sprintf("%s\nminus %s", B$arm_a, B$arm_b),
                levels = rev(sprintf("%s\nminus %s", B$arm_a, B$arm_b)))
B$note <- sprintf("%.2f   %s", B$percentile,
                  ifelse(B$p_emp_lower < 1 / ss$defs$n_set_draws,
                         sprintf("<%.4f", 1 / ss$defs$n_set_draws),
                         sprintf("%.4f", B$p_emp_lower)))
xrB   <- range(c(B$observed_diff, B$null_median, 0))
x_txB <- xrB[2] + diff(xrB) * 0.5

pB <- ggplot2::ggplot(B, ggplot2::aes(y = lab)) +
  ggplot2::geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.35) +
  ggplot2::geom_segment(ggplot2::aes(x = null_median, xend = observed_diff, yend = lab),
                        colour = "grey55", linewidth = 0.4,
                        arrow = ggplot2::arrow(length = ggplot2::unit(1.3, "mm"),
                                               type = "closed")) +
  ggplot2::geom_point(ggplot2::aes(x = null_median), shape = 21, size = 1.7,
                      fill = "white", colour = "grey40", stroke = 0.4) +
  ggplot2::geom_point(ggplot2::aes(x = observed_diff), size = 2.1,
                      colour = verdict_col[["withdraws"]]) +
  ggplot2::geom_text(ggplot2::aes(x = x_txB, label = note), hjust = 1, size = 1.85,
                     colour = "grey25") +
  ggplot2::annotate("text", x = x_txB, y = nrow(B) + 0.68, label = "pct     p",
                    hjust = 1, size = 1.85, colour = "grey40", fontface = "italic") +
  ggplot2::coord_cartesian(xlim = c(xrB[1] - diff(xrB) * 0.06, x_txB + 0.005),
                           ylim = c(0.4, nrow(B) + 1.0), expand = FALSE) +
  ggplot2::labs(
    x = "difference in wild-type 6 -> 12W log2FC", y = NULL,
    title = "B   Tested as a comparison",
    subtitle = "open = paired-null median (both sets redrawn together)") +
  theme_myc(base_size = 8) +
  ggplot2::theme(axis.text.y   = ggplot2::element_text(size = 6, lineheight = 0.95),
                 plot.title    = ggplot2::element_text(face = "bold", size = 8.5),
                 plot.subtitle = ggplot2::element_text(size = 6.2, colour = "grey25"))

# =============================================================================
# PANEL C -- and not an artefact of the content rise
# =============================================================================
if (is.null(ss$comparator_priority))
  stop("figure2 C: comparator_priority is NULL -- run script 40, then re-run 43")
Cd <- merge(ss$comparator[, c("arm", "c_wt_time")],
            ss$comparator_priority[, c("arm", "prio_wt_time")], by = "arm")
Cd <- Cd[is.finite(Cd$prio_wt_time), ]
Cd$verdict <- verdict_of(ss$wt_null$percentile[match(Cd$arm, ss$wt_null$arm)])

pC <- ggplot2::ggplot(Cd, ggplot2::aes(c_wt_time, prio_wt_time)) +
  ggplot2::geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +
  ggplot2::geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3) +
  ggplot2::geom_point(ggplot2::aes(colour = verdict), size = 2) +
  ggrepel::geom_text_repel(ggplot2::aes(label = arm, colour = verdict), size = 1.95,
                           seed = 1, min.segment.length = 0.2, segment.size = 0.25,
                           box.padding = 0.35, max.overlaps = Inf, show.legend = FALSE) +
  ggplot2::scale_colour_manual(values = verdict_col, guide = "none") +
  ggplot2::labs(
    x = "content  (set-average log2FC)", y = "priority  (mitoPPS)",
    title = "C   Not an artefact of the content rise",
    subtitle = "mitoPPS is a within-compartment ratio: blind to how\nmany mitochondria a cell has.") +
  theme_myc(base_size = 8) +
  ggplot2::theme(plot.title    = ggplot2::element_text(face = "bold", size = 8.5),
                 plot.subtitle = ggplot2::element_text(size = 6.2, colour = "grey25",
                                                       lineheight = 1.15))

# =============================================================================
# PANEL D -- the machinery is left intact
# =============================================================================
wg <- cmo$wt_genes
Dd <- data.frame(gene = wg$gene,
                 class = ifelse(wg$arm == "PRO", "pro-apoptotic", "anti-apoptotic"),
                 lfc = wg$wt_time, padj = wg$padj_wt, stringsAsFactors = FALSE)
bf <- ss$buffer[!ss$buffer$gene %in% Dd$gene, ]
if (nrow(bf))
  Dd <- rbind(Dd, data.frame(gene = bf$gene, class = "brake (IAP etc.)",
                             lfc = bf$lfc_wt_time, padj = bf$padj_wt_time,
                             stringsAsFactors = FALSE))
Dd <- Dd[is.finite(Dd$lfc), ]
Dd$gene[is.na(Dd$gene)] <- ""
Dd$sig <- !is.na(Dd$padj) & Dd$padj < 0.05
n_class <- table(Dd$class)
lv <- c("brake (IAP etc.)", "anti-apoptotic", "pro-apoptotic")
Dd$class <- factor(sprintf("%s  (%d)", Dd$class, n_class[Dd$class]),
                   levels = sprintf("%s  (%d)", lv, n_class[lv]))
labD <- Dd$gene %in% c("Bnip3", "Bbc3", "Bcl2l1", "Bcl2", "Mcl1", "Bax")
jitD <- ggplot2::position_jitter(width = 0, height = 0.17, seed = 7)

pD <- ggplot2::ggplot(Dd, ggplot2::aes(lfc, class)) +
  ggplot2::geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.35) +
  ggplot2::geom_point(ggplot2::aes(fill = sig), shape = 21, size = 1.8, stroke = 0.3,
                      colour = "grey35", position = jitD) +
  ggrepel::geom_text_repel(data = Dd[labD, ], ggplot2::aes(label = gene), size = 1.95,
                           seed = 7, position = jitD, min.segment.length = 0.15,
                           segment.size = 0.25, box.padding = 0.45,
                           point.padding = 0.2, max.overlaps = Inf) +
  ggplot2::scale_fill_manual(values = c(`FALSE` = "grey85",
                                        `TRUE` = verdict_col[["rises"]]), guide = "none") +
  ggplot2::labs(
    x = "wild-type 6 -> 12 weeks  (raw log2FC)", y = NULL,
    title = "D   The machine is left intact",
    subtitle = sprintf("%d of %d transcripts move (padj<0.05)\n-- and it RISES. No brake moves either.", sum(Dd$sig), nrow(Dd))) +
  theme_myc(base_size = 8) +
  ggplot2::theme(axis.text.y   = ggplot2::element_text(size = 6.2),
                 plot.title    = ggplot2::element_text(face = "bold", size = 8.5),
                 plot.subtitle = ggplot2::element_text(size = 6.2, colour = "grey25",
                                                       lineheight = 1.15))

# =============================================================================
# PANEL E -- the whole programme halves; PUMA:BCL-xL reverses
# =============================================================================
E <- merge(as.data.frame(pa$pair_null),
           as.data.frame(pa$priming)[, c("pair", "pro", "anti", "d6", "p6")], by = "pair")
E$real6 <- E$p6 < 0.05
E <- E[order(E$retention), ]
# build the label from the pro/anti COLUMNS, never by substituting into the pair
# string: sub("Bcl2l1", ...) turns Bcl2l11 (BIM) into "BCL-xL1".
nice <- c(Bbc3 = "PUMA", Bcl2l11 = "BIM", Bax = "BAX", Bak1 = "BAK", Bid = "BID",
          Bmf = "BMF", Pmaip1 = "NOXA", Bcl2l1 = "BCL-xL", Mcl1 = "MCL1")
lab_of <- function(g) ifelse(g %in% names(nice), nice[g], g)
E$pair_lab <- paste0(lab_of(E$pro), ":", lab_of(E$anti))
E$pair_lab <- factor(E$pair_lab, levels = rev(E$pair_lab))
E$note <- sprintf("%.0f", E$pct_retention_cond)

# the roster's own retentions, as the backdrop: the rest of the programme sits on
# the global rate. Restricted to genes with a real 6W effect, as script 42 does.
Rg <- as.data.frame(pa$machinery)
Rg <- Rg[is.finite(Rg$retention) & !is.na(Rg$padj_myc_6W) & Rg$padj_myc_6W < 0.05, ]

xrE   <- c(min(c(E$retention, Rg$retention)) - 0.1, 1.45)
x_txE <- xrE[2] - 0.02
E_off <- sum(E$retention > xrE[2])

pE <- ggplot2::ggplot(E, ggplot2::aes(y = pair_lab)) +
  ggplot2::geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3) +
  ggplot2::geom_vline(xintercept = RATE, colour = "grey30", linewidth = 0.4,
                      linetype = 2) +
  ggplot2::geom_segment(ggplot2::aes(x = null_retention_median, xend = retention,
                                     yend = pair_lab), colour = "grey60",
                        linewidth = 0.35) +
  ggplot2::geom_point(ggplot2::aes(x = null_retention_median), shape = 124, size = 1.9,
                      colour = "grey35") +
  ggplot2::geom_point(ggplot2::aes(x = retention, fill = real6), shape = 21, size = 2.1,
                      stroke = 0.35, colour = "grey25") +
  ggplot2::geom_text(ggplot2::aes(x = x_txE, label = note), hjust = 1, size = 1.85,
                     colour = "grey25") +
  # the rest of the MYC programme, as a rug on the same axis
  ggplot2::geom_point(data = Rg, ggplot2::aes(x = retention, y = nrow(E) + 0.85),
                      inherit.aes = FALSE, shape = 124, size = 1.6, colour = "grey45",
                      alpha = 0.8) +
  ggplot2::annotate("text", x = xrE[1] + 0.02, y = nrow(E) + 1.35, hjust = 0, size = 1.85,
                    colour = "grey30",
                    label = sprintf("the rest of the MYC programme (%d genes)", nrow(Rg))) +
  ggplot2::annotate("text", x = RATE, y = 0.55, hjust = -0.06, size = 1.9,
                    colour = "grey30", label = sprintf(" x%.2f", RATE)) +
  ggplot2::scale_fill_manual(values = c(`TRUE` = "#D73027", `FALSE` = "grey82"),
                             guide = "none") +
  ggplot2::coord_cartesian(xlim = xrE, ylim = c(0.4, nrow(E) + 1.75), expand = FALSE) +
  ggplot2::labs(
    x = "retention  (12W / 6W MYC effect)", y = NULL,
    title = "E   PUMA reverses",
    subtitle = sprintf("filled = the pair has a 6W effect;\ntick = the matched-pair null, numbers = its percentile.%s",
                       if (E_off) sprintf("\n%d pair sits beyond the axis.", E_off) else "")) +
  theme_myc(base_size = 8) +
  ggplot2::theme(axis.text.y   = ggplot2::element_text(size = 6.2),
                 plot.title    = ggplot2::element_text(face = "bold", size = 8.5),
                 plot.subtitle = ggplot2::element_text(size = 6.2, colour = "grey25",
                                                       lineheight = 1.15))

# =============================================================================
# PANEL F -- the three propositions, and the hand-off to the cell experiments
# -----------------------------------------------------------------------------
# H1/H2 come from script 42 PART H (`perm_axis`) and H3 from script 43 PART B
# (`tradeoff_perm`). The two models are on DIFFERENT coefficient scales, so the
# panel plots the PERMUTATION PERCENTILE, which is the one currency they share,
# and the coefficient rides as text. Never compare the coefficients across blocks.
# =============================================================================
px <- as.data.frame(pa$perm_axis)
tp <- as.data.frame(ss$tradeoff_perm)

grab_px <- function(ratio, axis) px[px$ratio == ratio & px$axis == axis, ][1, ]
grab_tp <- function(outcome, axis) tp[tp$outcome == outcome & tp$axis == axis, ][1, ]

Frow <- function(h, lab, obs, pct, ctrl, readable = TRUE)
  data.frame(h = h, lab = lab, observed = obs, percentile = pct, control = ctrl,
             readable = readable, stringsAsFactors = FALSE)

r_ox_puma  <- grab_px("Bbc3:Bcl2l1", "oxphos_ppd")
r_lv_puma  <- grab_px("Bbc3:Bcl2l1", "oxphos_lvl")
r_ox_bax   <- grab_px("Bax:Bcl2l1",  "oxphos_ppd")
r_rx_puma  <- grab_px("Bbc3:Bcl2l1", "redox_ppd")
r_rx_bax   <- grab_px("Bax:Bcl2l1",  "redox_ppd")
r_teb_puma <- grab_px("Bbc3:Bcl2l1", "teb")
r_teb_bax  <- grab_px("Bax:Bcl2l1",  "teb")
r_pr_mark  <- grab_tp("proliferation (markers)",      "oxphos_ppd")
r_pr_e2f   <- grab_tp("proliferation (E2F hallmark)", "oxphos_ppd")
c_pr_mark  <- grab_tp("proliferation (markers)",      "redox_ppd")
c_pr_e2f   <- grab_tp("proliferation (E2F hallmark)", "redox_ppd")

# CROSS-SCRIPT CHECK: the H1 test is implemented twice, independently. If the two
# implementations have drifted apart, the panel is not drawable.
h1_alt <- grab_tp("Bbc3:Bcl2l1 (PUMA priming)", "oxphos_ppd")
if (abs(r_ox_puma$percentile - h1_alt$percentile) > 1)
  stop(sprintf(paste("figure2 F: the two implementations of the H1 test disagree",
                     "(%.2f vs %.2f) -- re-run scripts 42 and 43 before drawing this."),
               r_ox_puma$percentile, h1_alt$percentile))

Fd <- rbind(
  Frow("H1", "respiratory priority -> PUMA:BCL-xL",  r_ox_puma$observed,
       r_ox_puma$percentile,  r_rx_puma$percentile),
  Frow("H1", "respiratory level -> PUMA:BCL-xL",     r_lv_puma$observed,
       r_lv_puma$percentile,  NA_real_),
  Frow("H1", "respiratory priority -> BAX:BCL-xL",   r_ox_bax$observed,
       r_ox_bax$percentile,   r_rx_bax$percentile, readable = FALSE),
  Frow("H2", "TEB signature -> PUMA:BCL-xL",         r_teb_puma$observed,
       r_teb_puma$percentile, NA_real_),
  Frow("H2", "TEB signature -> BAX:BCL-xL",          r_teb_bax$observed,
       r_teb_bax$percentile,  NA_real_, readable = FALSE),
  Frow("H3", "respiratory priority -> proliferation", r_pr_mark$observed,
       r_pr_mark$percentile,  c_pr_mark$percentile),
  Frow("H3", "respiratory priority -> E2F programme", r_pr_e2f$observed,
       r_pr_e2f$percentile,   c_pr_e2f$percentile))
Fd$lab <- factor(Fd$lab, levels = rev(Fd$lab))
Fd$note <- sprintf("%+.2f", Fd$observed)
h_col <- c(H1 = "#D55E00", H2 = "#0072B2", H3 = "#762A83")

pF_left <- ggplot2::ggplot(Fd, ggplot2::aes(y = lab)) +
  ggplot2::annotate("rect", xmin = 40, xmax = 60, ymin = -Inf, ymax = Inf,
                    fill = "grey92") +
  ggplot2::geom_vline(xintercept = 50, colour = "grey45", linewidth = 0.35,
                      linetype = 2) +
  ggplot2::geom_segment(ggplot2::aes(x = 50, xend = percentile, yend = lab, colour = h),
                        linewidth = 0.45) +
  ggplot2::geom_point(ggplot2::aes(x = control), shape = 124, size = 2.3,
                      colour = "grey30", na.rm = TRUE) +
  ggplot2::geom_point(ggplot2::aes(x = percentile, colour = h, alpha = readable),
                      size = 2.4) +
  ggplot2::geom_text(ggplot2::aes(x = 101, label = note), hjust = 1, size = 1.85,
                     colour = "grey25") +
  ggplot2::scale_colour_manual(values = h_col, guide = "none") +
  ggplot2::scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.28), guide = "none") +
  ggplot2::annotate("text", x = 2, y = nrow(Fd) + 0.75, hjust = 0, size = 1.85,
                    colour = "grey35", label = "MYC does LESS where the axis is high") +
  ggplot2::annotate("text", x = 98, y = nrow(Fd) + 0.75, hjust = 1, size = 1.85,
                    colour = "grey35", label = "MYC does MORE") +
  ggplot2::coord_cartesian(xlim = c(0, 102), ylim = c(0.4, nrow(Fd) + 1.15),
                           expand = FALSE) +
  ggplot2::labs(
    x = "percentile of the within-timepoint permutation null", y = NULL,
    title = "F   What the mouse hands to the cell experiments",
    subtitle = "Tick = the same test on `redox`, a matched axis MYC does not drive.\nFaded rows have a control that behaves like the axis, so they are unreadable.") +
  theme_myc(base_size = 8) +
  ggplot2::theme(axis.text.y   = ggplot2::element_text(size = 6.2),
                 plot.title    = ggplot2::element_text(face = "bold", size = 8.5),
                 plot.subtitle = ggplot2::element_text(size = 6.2, colour = "grey25",
                                                       lineheight = 1.15))

# --- the propositions themselves, as a text block with the H3 sign flip drawn --
txt <- data.frame(
  y    = c(9.6, 8.7, 8.0,   6.7, 5.8, 5.1,   3.8, 2.9, 2.2),
  h    = rep(c("H1", "H2", "H3"), each = 3),
  kind = rep(c("head", "body", "test"), 3),
  lab  = c(
    "H1   Respiratory state gates PUMA priming",
    "the effect on PUMA is not shared with BAX, and the control is dead centre",
    "TEST: raise PGC-1a in MYC-driven cells                                 DONE",
    "H2   The TEB-to-ductal transition moves both together",
    "TEB falls hardest of any arm and moderates the same ratio",
    "TEST: single-cell / deconvolution, or TEB- vs ductal-derived cells     OPEN",
    "H3   Respiration costs proliferation -- unless death is blocked",
    "in a mouse, death is intact, so the liability limb is the one expressed",
    "TEST: PGC-1a in MYAZ, with and without Bcl-xL                          DONE"),
  stringsAsFactors = FALSE)
txt$size  <- ifelse(txt$kind == "head", 2.35, 1.95)
txt$face  <- ifelse(txt$kind == "head", "bold", "plain")
txt$col   <- ifelse(txt$kind == "head", h_col[txt$h], "grey30")

flip <- data.frame(x = c(0.06, 0.06), xend = c(0.42, 0.42), y = c(0.85, 0.30),
                   lab = c("mouse, death INTACT:  more respiration -> LESS proliferation",
                           "MYAZ + Bcl-xL, death BLOCKED:  more PGC-1a -> MORE growth"),
                   col = c("#762A83", "#1B7837"), stringsAsFactors = FALSE)

pF_right <- ggplot2::ggplot() +
  ggplot2::geom_text(data = txt,
                     ggplot2::aes(x = 0.02, y = y, label = lab, colour = I(col),
                                  size = I(size), fontface = I(face)), hjust = 0) +
  ggplot2::annotate("segment", x = 0.02, xend = 0.98, y = c(7.4, 4.5), yend = c(7.4, 4.5),
                    colour = "grey85", linewidth = 0.3) +
  ggplot2::annotate("rect", xmin = 0.01, xmax = 0.99, ymin = -0.25, ymax = 1.35,
                    fill = "grey95", colour = NA) +
  ggplot2::geom_text(data = flip, ggplot2::aes(x = 0.03, y = y, label = lab,
                                               colour = I(col)),
                     hjust = 0, size = 1.95, fontface = "bold") +
  ggplot2::annotate("text", x = 0.03, y = 1.62, hjust = 0, size = 1.9, colour = "grey25",
                    label = "The same intervention, opposite sign -- the switch is whether the cell can die.") +
  ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(-0.5, 10.1), expand = FALSE) +
  ggplot2::theme_void()

pF <- patchwork::wrap_plots(pF_left, pF_right, nrow = 1, widths = c(1, 1.25))

# =============================================================================
# ASSEMBLY
# =============================================================================
row1 <- patchwork::wrap_plots(pA, pC, nrow = 1, widths = c(1.5, 1))
row2 <- patchwork::wrap_plots(pB, pD, pE, nrow = 1, widths = c(0.95, 0.95, 1.35))

p <- patchwork::wrap_plots(row1, row2, pF, ncol = 1, heights = c(1, 0.92, 1.02)) +
  patchwork::plot_annotation(
    caption = paste(
      "A-D: wild-type only, n=6 per timepoint. Nulls are 2000 expression-matched random sets; B redraws BOTH sets of a contrast together, which is the null a comparative claim needs.",
      "E: retention is read against the x0.55 rescaling of the whole MYC programme -- the null is 'fades like everything else', not 'no change'. 4000 baseMean-matched pairs.",
      "F: RANKINGS, NOT TESTS. Nothing beats its permutation null (best p 0.081); the respiratory axis is collinear with the timepoint it would explain; H1 and H2 are not separable at n=24.",
      "   The two source models are on different coefficient scales, so only the percentile is comparable across blocks. The H1 test is implemented twice independently: 91.9th and 91.7th percentile.",
      "BATCH = TIMEPOINT: the two cohorts were extracted as separate batches, so every wild-type temporal value is DESCRIBED, not claimed. The mouse generates; the perturbations prove.",
      sep = "\n"),
    theme = ggplot2::theme(
      plot.caption = ggplot2::element_text(size = 5.6, hjust = 0, colour = "grey30",
                                           lineheight = 1.15)))

if (!isTRUE(getOption("myc.fig.nosave"))) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  save_panel(p, file.path(out_dir, "figure2_developmental_window.pdf"),
             width = fig_w[["double"]], height = 215)
}

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {
  ## A/B: the withdrawal and the comparative test
  A[, c("arm", "c_wt_time", "percentile")] |> print()
  B[, c("arm_a", "arm_b", "observed_diff", "percentile", "p_emp_lower")] |> print()

  ## D: the one mover
  Dd[Dd$sig, ] |> print()

  ## E: and the honesty check the panel encodes -- Bmf sits further out than PUMA
  ## on the conditional null but has NO 6W effect to lose.
  E[, c("pair", "d6", "p6", "retention", "null_retention_median",
        "pct_retention_cond", "p_emp_cond")] |> print()

  ## F: the propositions table as drawn, plus the cross-script agreement
  Fd |> print()
  cat(sprintf("H1 implemented twice: %.2f (script 42) vs %.2f (script 43)\n",
              r_ox_puma$percentile, h1_alt$percentile))

  ## and the limit that keeps F a ranking: the axis is collinear with timepoint
  ss$tradeoff[, c("outcome", "axis", "myc_x_axis", "p",
                  "myc_x_axis_with_tp", "p_with_tp")] |> print()

  print(pA); print(pB); print(pC); print(pD); print(pE); print(pF)
  print(p)
}
