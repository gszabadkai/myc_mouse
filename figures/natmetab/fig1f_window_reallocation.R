# =============================================================================
# fig1f_window_reallocation.R -- Fig. 1f: across the pubertal-adult window the
# gland withdraws from respiration, not from proliferation (ED Fig. 2f's tests folded in)
# -----------------------------------------------------------------------------
# PROMOTED from proposals/ on 2026-10-14 (author: "promote them to final panels").
# The 2026-10-11 Fig. 1f is retired/fig1f_window_reallocation.R, and the old ED 2f
# is retired/edfig2f_window_matched_controls.R.
#
# The author (2026-10-14): "Fig. 1f and ED Fig 2f could be combined into one panel
# going on the main panel. The ED figure adds the statistics, but does not add more
# information ... The tick for the median are difficult to grasp -- why there are
# two and which one is which? Probably could be left out."
#
# Fig. 1f's rows and two rulers (fig1f_window_reallocation.R, the panel of record,
# unchanged), with ED 2f's statistics drawn on them and its median ticks dropped:
#   EXPRESSION  filled = P < 0.05 against 2,000 random gene sets of the same size
#               from the same expression ventiles (beyond their central 95%,
#               two-sided); open = not. Script 54's draws, both windows.
#   PRIORITY    filled = BH-adjusted P < 0.05, script 08's Welch test -- Fig. 1e's
#               grammar on the same ruler. No window change of these six passes
#               (0 of 144 in WT), so every priority point is open; the legend gives
#               the ranks among the 143 the text leans on.
#   BRACKET     OXPHOS subunits minus proliferation, WT window, against paired
#               matched sets redrawn together (script 43): the text's comparison.
#               The rows sit next to each other so that it can be drawn; the other
#               two paired comparisons (mitoribosome, nucleotide metabolism) go to
#               the legend.
# The two groups share one panel (no facets) so that the bracket can cross from a
# development programme to a mitochondrial one; the group names sit at the right.
#
# THE TWO EXTRA ROWS: script 54's matched null has no arm for the carnitine shuttle
# or fatty acid oxidation, which Fig. 1f draws for priority. Script 56 PART A drew
# theirs by script 54's recipe (proved there on three of script 54's own arms);
# this panel reads it.
#
# BATCH: 6W and 12W were extracted as separate batches, so both windows coincide
# with batch; matched sets and paired sets share it (legend).
#
# TYPE (2026-10-15, as Fig. 1d and 1e): row and tick labels, the key and the P value
# 5 pt, the axis titles and the group names 6 pt.
#
# Reads:  results/background_vs_myc.rds            (script 40: $ruler)
#         results/substrate_specificity_tradeoff.rds (script 43: $comparator, $paired_null)
#         results/two_timeline_verification.rds    (script 54: $arm_null)
#         results/window_and_chain_statistics.rds  (script 56: $window_null, the two extra rows)
#         results/mitopps_scores.rds               (script 08: $mitopps_pairwise)
# Output: outputs/natmetab/Fig1/Fig1f_window_reallocation.pdf
# =============================================================================

source(here::here("figures", "natmetab", "_style.R"))

ruler <- as.data.frame(readRDS(here::here("results", "background_vs_myc.rds"))$ruler)
s43   <- readRDS(here::here("results", "substrate_specificity_tradeoff.rds"))
comp  <- as.data.frame(s43$comparator)
pn    <- as.data.frame(s43$paired_null)
tv    <- readRDS(here::here("results", "two_timeline_verification.rds"))
an    <- as.data.frame(tv$arm_null)
wn56  <- as.data.frame(readRDS(here::here("results", "window_and_chain_statistics.rds"))$window_null)
pw    <- as.data.frame(readRDS(here::here("results", "mitopps_scores.rds"))$mitopps_pairwise)

TXT   <- 5      # pt; row and tick labels, the key and the P value
TITLE <- 6      # pt; the axis titles and the group names

GROUPS <- c(dev = "Development\nprogrammes", mito = "Mitochondrial\npathways")

# --- the rows: Fig. 1f's, with the arm each has in script 54 ------------------------
dev <- data.frame(label = c("TEB vs ductal", "Proliferation"),
                  arm = c("TEB vs ductal (HS)", "PROLIF_* pooled"), group = "dev",
                  stringsAsFactors = FALSE)
k <- match(dev$arm, comp$arm); stopifnot(!anyNA(k))
dev$n <- comp$n_genes[k]
dev$expr_neg <- comp$c_wt_time[k]; dev$expr_pos <- comp$c_myc_time[k]
dev$prio_neg <- NA_real_;          dev$prio_pos <- NA_real_

mito <- data.frame(label = c("OXPHOS subunits", "Mitochondrial ribosome", "OXPHOS assembly factors",
                             "TCA cycle", "Carnitine shuttle", "Fatty acid oxidation"),
                   arm = c("OXPHOS subunits", "mitoribosome", "OXPHOS assembly", "TCA cycle", NA, NA),
                   set = c(NA, NA, NA, NA, "MITOCARTA_CARNITINE_SHUTTLE", "MITOCARTA_FATTY_ACID_OXIDATION"),
                   group = "mito", stringsAsFactors = FALSE)
k <- match(mito$label, ruler$pathway); stopifnot(!anyNA(k))
mito$n <- ruler$n_genes[k]
mito$expr_neg <- ruler$c_tn[k]; mito$expr_pos <- ruler$c_tp[k]
mito$prio_neg <- ruler$p_tn[k]; mito$prio_pos <- ruler$p_tp[k]
dev$set <- NA_character_

rows <- rbind(dev[, names(mito)], mito)

# --- expression: script 54's matched null where it has the arm --------------------
k <- match(rows$arm, an$arm)
stopifnot(identical(is.na(k), is.na(rows$arm)))
rows$pct_neg <- an$pct_wt[k]; rows$pct_pos <- an$pct_myc[k]
# the observed values of record agree with the drawn ones
stopifnot(max(abs(an$obs_wt[k] - rows$expr_neg), na.rm = TRUE) < 1e-8,
          max(abs(an$obs_myc[k] - rows$expr_pos), na.rm = TRUE) < 1e-8)

# --- the two rows script 54 has no arm for: script 56 PART A ----------------------------
k56 <- match(rows$label[!is.na(rows$set)], wn56$arm)
stopifnot(!anyNA(k56), identical(wn56$set[k56], rows$set[!is.na(rows$set)]),
          max(abs(wn56$obs_wt[k56]  - rows$expr_neg[!is.na(rows$set)])) < 1e-8,    # the drawn values
          max(abs(wn56$obs_myc[k56] - rows$expr_pos[!is.na(rows$set)])) < 1e-8)
rows$pct_neg[!is.na(rows$set)] <- wn56$pct_wt[k56]
rows$pct_pos[!is.na(rows$set)] <- wn56$pct_myc[k56]
stopifnot(!anyNA(rows$pct_neg), !anyNA(rows$pct_pos))

# --- priority: script 08's Welch test, BH across the MitoPathways --------------------
wt08 <- pw[pw$contrast == "Temporal_Myc-", ]; my08 <- pw[pw$contrast == "Temporal_Myc+", ]
rows$padj_neg <- wt08$padj[match(rows$label, wt08$pathway)]
rows$padj_pos <- my08$padj[match(rows$label, my08$pathway)]
stopifnot(identical(is.na(rows$padj_neg), rows$group == "dev"),
          max(abs(wt08$diff[match(rows$label, wt08$pathway)] - rows$prio_neg), na.rm = TRUE) < 1e-6)

# --- the paired comparison the bracket draws (WT window, script 43) ---------------------
pp <- pn[pn$arm_a == "OXPHOS subunits" & pn$arm_b == "PROLIF_* pooled", ]
stopifnot(nrow(pp) == 1L, pp$observed_diff < 0)
p_lab <- if (pp$p_emp_lower == 0) "italic(P)<'0.0005'" else fmt_p(pp$p_emp_lower, exponent = "reduced")

# --- rows top to bottom, one panel, a gap between the groups -------------------------------
rows <- rows[order(match(rows$group, names(GROUPS)), rows$expr_neg), ]
stopifnot(identical(rows$label[2:3], c("Proliferation", "OXPHOS subunits")))
rows$label_n <- sprintf("%s (%d)", rows$label, rows$n)
# discrete rows, bottom to top, with one empty level between the groups
LEV <- c(rev(rows$label_n[rows$group == "mito"]), " ", rev(rows$label_n[rows$group == "dev"]))
rows$y <- match(rows$label_n, LEV)
grp_y <- tapply(rows$y, rows$group, mean)

sig_state <- function(sig) factor(ifelse(sig, "sig", "ns"), levels = c("sig", "ns"))
two_sided <- function(pct) pct <= 2.5 | pct >= 97.5
longer <- function(a, b, sa, sb) {
  d <- rbind(data.frame(rows[, c("label_n", "y")], geno = "neg", val = rows[[a]], sig = sa),
             data.frame(rows[, c("label_n", "y")], geno = "pos", val = rows[[b]], sig = sb))
  d$label_n <- factor(d$label_n, levels = LEV)
  d$geno <- factor(d$geno, levels = c("neg", "pos"))
  d$sig  <- sig_state(d$sig)
  d
}
ex <- longer("expr_neg", "expr_pos", two_sided(rows$pct_neg), two_sided(rows$pct_pos))
pr <- longer("prio_neg", "prio_pos", rows$padj_neg < 0.05, rows$padj_pos < 0.05)
pr <- pr[!is.na(pr$val), ]

# --- the text's numbers, asserted (as on the panel of record) --------------------------
r143 <- ruler[!ruler$is_mtdna, ]
ox <- r143$pathway == "OXPHOS subunits"
stopifnot(all(rows$pct_neg[rows$label == "OXPHOS subunits"] == 0, rows$pct_pos[rows$label == "OXPHOS subunits"] == 0),
          rows$pct_neg[rows$label == "TEB vs ductal"] == 0,
          mean(r143$p_tn > r143$p_tn[ox]) >= 0.9, mean(r143$p_tp > r143$p_tp[ox]) >= 0.9,
          !any(pr$sig == "sig"))                      # every priority point is open, as the header says

# --- drawing ------------------------------------------------------------------------------
PT_SIZE <- c(neg = 1.6, pos = 1.1)
pad <- function(v) diff(range(v, na.rm = TRUE)) * 0.04

one_ruler <- function(d, x0, x1, xtitle, show_rows) {
  seg <- stats::reshape(d[, c("label_n", "geno", "val")], idvar = "label_n",
                        timevar = "geno", direction = "wide")
  g <- ggplot(d, aes(val, label_n)) +
    geom_vline(xintercept = 0, linewidth = NM_LINE, colour = box_line) +
    geom_segment(data = seg, aes(x = val.neg, xend = val.pos, y = label_n, yend = label_n),
                 inherit.aes = FALSE, linewidth = 0.4, colour = "grey70") +
    geom_point(aes(colour = geno, size = geno, shape = sig), stroke = 0.5) +
    scale_colour_manual(values = geno_cols, labels = dev_labels, name = NULL) +
    scale_size_manual(values = PT_SIZE, guide = "none") +
    scale_shape_manual(values = c(sig = 16, ns = 1), drop = FALSE, name = NULL,
                       labels = c(sig = expression(italic(P) < "0.05"), ns = "Not significant")) +
    # oob_keep: the group names sit right of the priority ruler's last tick
    scale_x_nm(step = 0.2, top = x1, bottom = x0, minor = TRUE, labels = lab_signed,
               oob = scales::oob_keep) +
    scale_y_discrete(limits = LEV, drop = FALSE) +
    labs(x = xtitle, y = NULL) +
    guides(colour = guide_legend(order = 1, override.aes = list(size = unname(PT_SIZE), shape = 16)),
           shape  = guide_legend(order = 2, override.aes = list(size = 1.4, colour = "black"))) +
    coord_cartesian(clip = "off") +
    theme_nm(legend = "bottom", txt = TXT, txt_title = TITLE) +
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank())
  if (show_rows) g else g + theme(axis.text.y = element_blank())
}

# the bracket: right of both rows' points, tips toward them, P beside it
XB <- 0.05; TIP <- 0.02
yb <- rows$y[rows$label %in% c("Proliferation", "OXPHOS subunits")]
stopifnot(XB > max(ex$val[ex$label_n %in% rows$label_n[rows$label %in% c("Proliferation", "OXPHOS subunits")]]))
br <- data.frame(x = c(XB - TIP, XB, XB), xend = c(XB, XB, XB - TIP),
                 y = c(yb[1], yb[1], yb[2]), yend = c(yb[1], yb[2], yb[2]))

p_ex <- one_ruler(ex, min(ex$val) - pad(ex$val), max(ex$val) + pad(ex$val),
                  "Expression\n(log2 fold change)", TRUE) +
  geom_segment(data = br, aes(x = x, xend = xend, y = y, yend = yend), inherit.aes = FALSE,
               linewidth = NM_LINE, colour = "black") +
  annotate("text", x = XB + 0.02, y = mean(yb), label = p_lab, parse = TRUE, hjust = 0,
           size = pt(TXT), family = NM_FONT)
p_pr <- one_ruler(pr, min(pr$val) - pad(pr$val), max(pr$val) + pad(pr$val),
                  "Priority\n(mitoPPS)", FALSE) +
  guides(shape = "none")              # one shape key, from the expression ruler, for both
x_right <- layer_scales(p_pr)$x$get_limits()[2]
p_pr <- p_pr +
  annotate("text", x = x_right + 0.03, y = grp_y[names(GROUPS)], label = GROUPS, hjust = 0,
           vjust = 0.5, lineheight = 0.95, size = pt(TITLE), family = NM_FONT) +
  theme(plot.margin = margin(1, 15, 1, 1, "mm"))

save_panel(stack_parts(p_ex, p_pr, ncol = 2, widths = c(1.15, 1), collect_guides = TRUE),
           fig = "Fig1", panel = "f", name = "window_reallocation", width = 100, height = 56)

# --- numbers for the legend (printed, never drawn) --------------------------------
cat("\n1f (with ED 2f) -- for the legend\n")
cat("  12W vs 6W within each genotype; n = 6 mice per group. Expression: set-average raw log2 fold change;\n")
cat("  filled = beyond the central 95% of 2,000 random sets of the same size from the same expression\n")
cat("  ventiles (two-sided empirical P < 0.05; script 54, the two marked * from script 56 by its recipe)\n")
for (i in seq_len(nrow(rows))) cat(sprintf("    %-28s%s expr WT %+.3f (pct %5.1f)  MYC+ %+.3f (pct %5.1f)%s\n",
  rows$label_n[i], if (!is.na(rows$set[i])) "*" else " ", rows$expr_neg[i], rows$pct_neg[i],
  rows$expr_pos[i], rows$pct_pos[i],
  if (is.na(rows$prio_neg[i])) "" else sprintf("   priority WT %+.3f (adj %.2f, %s)  MYC+ %+.3f (adj %.2f, %s)",
    rows$prio_neg[i], rows$padj_neg[i],
    if (rows$prio_neg[i] < 0) sprintf("%d lowest", sum(r143$p_tn <= rows$prio_neg[i])) else sprintf("%d highest", sum(r143$p_tn >= rows$prio_neg[i])),
    rows$prio_pos[i], rows$padj_pos[i],
    if (rows$prio_pos[i] < 0) sprintf("%d lowest", sum(r143$p_tp <= rows$prio_pos[i])) else sprintf("%d highest", sum(r143$p_tp >= rows$prio_pos[i])))))
cat("  bracket and legend: WT window, OXPHOS subunits minus comparator against 2,000 paired matched draws (script 43):\n")
for (i in seq_len(nrow(pn))) cat(sprintf("    OXPHOS - %-24s %+.3f, P %s\n", pn$arm_b[i], pn$observed_diff[i],
                                         if (pn$p_emp_lower[i] == 0) "< 0.0005" else formatC(pn$p_emp_lower[i], format = "fg", digits = 2)))
cat(sprintf("  priority: script 08's Welch test, BH across MitoPathways; window changes passing: WT %d, MYC+ %d of %d\n",
            sum(wt08$padj < 0.05, na.rm = TRUE), sum(my08$padj < 0.05, na.rm = TRUE), nrow(wt08)))
cat("  BATCH: each window coincides with an extraction batch; matched and paired sets share it\n")
