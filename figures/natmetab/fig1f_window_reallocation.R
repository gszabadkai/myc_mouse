# =============================================================================
# fig1f_window_reallocation.R -- Fig. 1f: across the pubertal-adult window the
# gland withdraws from respiration, not from proliferation (ED Fig. 2f's tests folded in)
# -----------------------------------------------------------------------------
# 2026-10-16 (author: "the actual OXPHOS deprioritisation is weak"): the RIGHT SIDE
# is now the window scatter, the old ED Fig. 2f (retired/edfig2f_window_two_rulers.R)
# with its two facets overlaid: every nuclear MitoPathway's change in expression
# against its change in priority, WT and MYC+ in their colours, "WT r" and "MYC+ r"
# on the panel. OXPHOS subunits, the carnitine shuttle and fatty acid oxidation (the
# text's three) have their two points joined by a thin dashed line and one black
# label with a leader to the WT point (author, later on 2026-10-16: the joins thinner
# and dashed, one leader each). The priority ruler and the group names ("Development programmes",
# "Mitochondrial pathways") are retired; the priority numbers stay in the legend.
# The left side, the expression ruler with its tests and bracket, is unchanged.
#
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
# development programme to a mitochondrial one; a blank row parts them (their names,
# at the right of the priority ruler, were retired with it on 2026-10-16).
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
# 5 pt, the axis titles 6 pt.
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
TITLE <- 6      # pt; the axis titles

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

one_ruler <- function(d, x0, x1, xtitle) {
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
    scale_x_nm(step = 0.2, top = x1, bottom = x0, minor = TRUE, labels = lab_signed) +
    scale_y_discrete(limits = LEV, drop = FALSE) +
    labs(x = xtitle, y = NULL) +
    guides(colour = guide_legend(order = 1, override.aes = list(size = unname(PT_SIZE), shape = 16)),
           shape  = guide_legend(order = 2, override.aes = list(size = 1.4, colour = "black"))) +
    coord_cartesian(clip = "off") +
    theme_nm(legend = "bottom", txt = TXT, txt_title = TITLE) +
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank())
  g
}

# the bracket: right of both rows' points, tips toward them, P beside it
XB <- 0.05; TIP <- 0.02
yb <- rows$y[rows$label %in% c("Proliferation", "OXPHOS subunits")]
stopifnot(XB > max(ex$val[ex$label_n %in% rows$label_n[rows$label %in% c("Proliferation", "OXPHOS subunits")]]))
br <- data.frame(x = c(XB - TIP, XB, XB), xend = c(XB, XB, XB - TIP),
                 y = c(yb[1], yb[1], yb[2]), yend = c(yb[1], yb[2], yb[2]))

p_ex <- one_ruler(ex, min(ex$val) - pad(ex$val), max(ex$val) + pad(ex$val),
                  "Expression\n(log2 fold change)") +
  geom_segment(data = br, aes(x = x, xend = xend, y = y, yend = yend), inherit.aes = FALSE,
               linewidth = NM_LINE, colour = "black") +
  annotate("text", x = XB + 0.02, y = mean(yb), label = p_lab, parse = TRUE, hjust = 0,
           size = pt(TXT), family = NM_FONT)
# --- the right side: the window scatter (author, 2026-10-16) ------------------------------
# The old ED 2f, its two facets overlaid: every nuclear MitoPathway's change in
# expression against its change in priority, WT and MYC+ in their colours. The three
# pathways the text names have their two points joined (thin, dashed), and one black
# label with a leader to the WT point; WT r and MYC+ r are written on the panel. It replaces the
# priority ruler and the group names, retired.
r143 <- ruler[!ruler$is_mtdna, ]
stopifnot(nrow(r143) == 143L)
sc <- rbind(data.frame(pathway = r143$pathway, geno = "neg", x = r143$c_tn, y = r143$p_tn),
            data.frame(pathway = r143$pathway, geno = "pos", x = r143$c_tp, y = r143$p_tp))
sc$geno <- factor(sc$geno, levels = c("neg", "pos"))
rr <- c(neg = stats::cor(r143$c_tn, r143$p_tn), pos = stats::cor(r143$c_tp, r143$p_tp))

NAMED <- c("OXPHOS subunits", "Carnitine shuttle", "Fatty acid oxidation")
hl   <- sc[sc$pathway %in% NAMED, ]
link <- data.frame(pathway = NAMED,
                   x  = r143$c_tn[match(NAMED, r143$pathway)], y  = r143$p_tn[match(NAMED, r143$pathway)],
                   xe = r143$c_tp[match(NAMED, r143$pathway)], ye = r143$p_tp[match(NAMED, r143$pathway)])
# Each label is placed by hand in empty space (checked on the drawn panel), with a
# leader from its anchor to the pathway's WT point.
lab <- data.frame(pathway = NAMED,
                  text = c("OXPHOS\nsubunits", "Carnitine\nshuttle", "Fatty acid\noxidation"),
                  lx = c(-0.10, 0.26, 0.30), ly = c(-0.34, 0.34, -0.10), hj = c(0, 0, 0))
lab$x0 <- lab$lx + ifelse(lab$hj == 0, -0.015, 0.015)          # the leaders start just off the text
# one leader per label, to the WT point (author, 2026-10-16)
lead <- merge(lab[, c("pathway", "x0", "ly")], hl[hl$geno == "neg", c("pathway", "x", "y")], by = "pathway")
stopifnot(nrow(lead) == length(NAMED))
xl <- range(sc$x) + c(-1, 1) * 0.03 * diff(range(sc$x))
yl <- range(sc$y) + c(-1, 1) * 0.04 * diff(range(sc$y))
# WT r and MYC+ r, two lines in the empty upper-left corner (atop() would set them far apart)
r_lab <- data.frame(y = c(0.38, 0.38 - 0.045),
                    label = sprintf('"%s"~italic(r)=="%.2f"', c("WT", "MYC+"), rr))

p_sc <- ggplot(sc, aes(x, y)) +
  geom_hline(yintercept = 0, linewidth = NM_LINE, colour = "grey80") +
  geom_vline(xintercept = 0, linewidth = NM_LINE, colour = "grey80") +
  geom_point(aes(colour = geno), size = 0.7, alpha = 0.5, stroke = 0) +
  geom_segment(data = lead, aes(x = x0, y = ly, xend = x, yend = y), inherit.aes = FALSE,
               linewidth = 0.25, colour = "grey50") +
  geom_segment(data = link, aes(x = x, y = y, xend = xe, yend = ye), inherit.aes = FALSE,
               linewidth = NM_LINE, linetype = "22", colour = "grey30") +
  geom_point(data = hl, aes(fill = geno), shape = 21, size = 1.4, colour = "black", stroke = 0.35) +
  geom_text(data = lab, aes(x = lx, y = ly, label = text, hjust = hj), inherit.aes = FALSE,
            size = pt(TXT), family = NM_FONT, lineheight = 0.9, colour = "black") +
  geom_text(data = r_lab, aes(x = xl[1] + 0.03, y = y, label = label), inherit.aes = FALSE,
            parse = TRUE, hjust = 0, vjust = 1, size = pt(TXT), family = NM_FONT, colour = "black") +
  scale_colour_manual(values = geno_cols, guide = "none") +       # the key is the ruler's
  scale_fill_manual(values = geno_cols, guide = "none") +
  scale_x_nm(step = 0.4, top = xl[2], bottom = xl[1], minor = TRUE, labels = lab_signed) +
  scale_y_nm(step = 0.2, top = yl[2], bottom = yl[1], minor = TRUE, labels = lab_signed) +
  labs(x = "Expression\n(log2 fold change)", y = "Priority (mitoPPS)") +
  coord_cartesian(clip = "off") +
  theme_nm(txt = TXT, txt_title = TITLE)

save_panel(stack_parts(p_ex, p_sc, ncol = 2, widths = c(1, 1.15), collect_guides = TRUE),
           fig = "Fig1", panel = "f", name = "window_reallocation", width = 100, height = 56)

# --- numbers for the legend (printed, never drawn) --------------------------------
cat("\nFig. 1f -- for the legend\n")
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
cat(sprintf("  RIGHT: %d nuclear MitoPathways, window change in expression (set-average raw log2 fold change)\n", nrow(r143)))
cat(sprintf("  against priority (mitoPPS); r(expression, priority): WT %.2f, MYC+ %.2f; named and joined:\n", rr[["neg"]], rr[["pos"]]))
for (k in NAMED) cat(sprintf("    %-22s WT (%+.3f, %+.3f)  MYC+ (%+.3f, %+.3f)\n", k,
  r143$c_tn[r143$pathway == k], r143$p_tn[r143$pathway == k], r143$c_tp[r143$pathway == k], r143$p_tp[r143$pathway == k]))
cat("  BATCH: each window coincides with an extraction batch; matched and paired sets share it\n")
