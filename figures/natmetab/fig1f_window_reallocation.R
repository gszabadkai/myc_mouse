# =============================================================================
# fig1f_window_reallocation.R -- Fig. 1f: across the pubertal-adult window the
# gland withdraws from respiration, not from proliferation
# -----------------------------------------------------------------------------
# The text: "between 6 and 12 weeks, both the WT and the MYC+ gland demoted
# OXPHOS in priority and in expression, while fatty acid oxidation and the
# carnitine shuttle, which supply NADH to the respiratory chain, gained priority
# (Fig. 1f, ED Fig. 2g). In the WT gland, OXPHOS subunits fell further than all
# 2,000 expression-matched gene sets, whereas the proliferative signature fell
# only slightly ... (Fig. 1d,f; ED Fig. 2f)".
#
# ONE PANEL, TWO RULERS SHARING ROWS (one PDF, side by side):
#   left   EXPRESSION: the set-average raw log2 fold change, 12W vs 6W
#   right  PRIORITY: the mitoPPS change, 12W vs 6W -- defined only for
#          mitochondrial pathways, so the two development programmes (TEB,
#          proliferation) have no point there
# Each row has the WT window and the MYC+ window, in the genotypes' colours, the
# MYC+ point smaller; the key says "12W vs 6W", so the two axis titles name only
# the ruler and its unit.
#
# THE ROWS: the two development programmes the sentence compares against, then
# the mitochondrial arms ordered by the WT expression change, so that
# proliferation and OXPHOS subunits sit next to each other.
#
# FORM: signed changes, two conditions -> dot plot, the pair joined by a line
#   (data-to-viz, "Order your data"). No significance encoding: the
#   matched-control test of expression is ED Fig. 2f; priority has no matched null,
#   and no window change passes script 08's per-pathway test in WT (0 of 144), so
#   the legend gives priority as ranks among the 143 (printed below).
#
# BATCH: 6W and 12W were extracted as separate batches, so both windows coincide
# with batch. The relative readings (against expression-matched sets, OXPHOS
# against proliferation, mitoPPS within the compartment) cancel a batch effect
# that acts on all genes or by expression level; the OXPHOS protein decline is
# Fig. 1b. Say both in the legend.
#
# Reads:  results/background_vs_myc.rds          (script 40: $ruler)
#         results/substrate_specificity_tradeoff.rds (script 43: $comparator,
#             $wt_null; legend numbers)
#         results/two_timeline_verification.rds  (script 54: $arm_null; legend)
#         results/mitopps_scores.rds             (script 08: $mitopps_pairwise; legend)
# Output: outputs/natmetab/Fig1/Fig1f_window_reallocation.pdf
# =============================================================================

source(here::here("figures", "natmetab", "_style.R"))

ruler <- as.data.frame(readRDS(here::here("results", "background_vs_myc.rds"))$ruler)
s43   <- readRDS(here::here("results", "substrate_specificity_tradeoff.rds"))
comp  <- as.data.frame(s43$comparator)
tv    <- readRDS(here::here("results", "two_timeline_verification.rds"))

GROUPS <- c(dev = "Development\nprogrammes", mito = "Mitochondrial\npathways")

dev <- data.frame(label = c("TEB vs ductal", "Proliferation"),
                  arm = c("TEB vs ductal (HS)", "PROLIF_* pooled"), group = "dev",
                  stringsAsFactors = FALSE)
k <- match(dev$arm, comp$arm); stopifnot(!anyNA(k))
dev$n <- comp$n_genes[k]
dev$expr_neg <- comp$c_wt_time[k]; dev$expr_pos <- comp$c_myc_time[k]
dev$prio_neg <- NA_real_;          dev$prio_pos <- NA_real_

mito <- data.frame(pathway = c("OXPHOS subunits", "Mitochondrial ribosome",
                               "OXPHOS assembly factors", "TCA cycle",
                               "Carnitine shuttle", "Fatty acid oxidation"),
                   group = "mito", stringsAsFactors = FALSE)
k <- match(mito$pathway, ruler$pathway); stopifnot(!anyNA(k))
mito$label <- mito$pathway
mito$n <- ruler$n_genes[k]
mito$expr_neg <- ruler$c_tn[k]; mito$expr_pos <- ruler$c_tp[k]
mito$prio_neg <- ruler$p_tn[k]; mito$prio_pos <- ruler$p_tp[k]

# the two objects compute the same set-average the same way: cross-check them
stopifnot(abs(mito$expr_neg[1] - comp$c_wt_time[comp$arm == "OXPHOS subunits"]) < 1e-8,
          abs(mito$expr_pos[1] - comp$c_myc_time[comp$arm == "OXPHOS subunits"]) < 1e-8)

rows <- rbind(dev[, c("label", "group", "n", "expr_neg", "expr_pos", "prio_neg", "prio_pos")],
              mito[, c("label", "group", "n", "expr_neg", "expr_pos", "prio_neg", "prio_pos")])
rows <- rows[order(match(rows$group, names(GROUPS)), rows$expr_neg), ]   # top row first
rows$label_n <- factor(sprintf("%s (%d)", rows$label, rows$n),
                       levels = rev(sprintf("%s (%d)", rows$label, rows$n)))
rows$group <- factor(GROUPS[rows$group], levels = GROUPS)
stopifnot(identical(as.character(rows$label[2:3]), c("Proliferation", "OXPHOS subunits")))

longer <- function(a, b) {
  d <- rbind(data.frame(rows[, c("label_n", "group")], geno = "neg", val = rows[[a]]),
             data.frame(rows[, c("label_n", "group")], geno = "pos", val = rows[[b]]))
  d$geno <- factor(d$geno, levels = c("neg", "pos"))
  d
}
# the WT window is the gland's own programme, the one the text leads with; the MYC+
# window is drawn smaller (author, 2026-10-11)
PT_SIZE <- c(neg = 1.6, pos = 1.1)
ex <- longer("expr_neg", "expr_pos")
pr <- longer("prio_neg", "prio_pos")

one_ruler <- function(d, x0, x1, xtitle, step, show_rows) {
  pts <- d[!is.na(d$val), ]
  seg <- stats::reshape(pts[, c("label_n", "group", "geno", "val")], idvar = c("label_n", "group"),
                        timevar = "geno", direction = "wide")
  g <- ggplot(d, aes(val, label_n)) +
    geom_blank(aes(x = 0)) +                         # every row exists in both rulers
    geom_vline(xintercept = 0, linewidth = NM_LINE, colour = box_line) +
    geom_segment(data = seg, aes(x = val.neg, xend = val.pos, y = label_n, yend = label_n),
                 inherit.aes = FALSE, linewidth = 0.4, colour = "grey70") +
    geom_point(data = pts, aes(colour = geno, size = geno)) +
    facet_grid(group ~ ., scales = "free_y", space = "free_y") +
    scale_colour_manual(values = geno_cols, labels = dev_labels, name = NULL) +
    scale_size_manual(values = PT_SIZE, guide = "none") +
    scale_x_nm(step = step, top = x1, bottom = x0, minor = TRUE, labels = lab_signed) +
    labs(x = xtitle, y = NULL) +
    guides(colour = guide_legend(override.aes = list(size = unname(PT_SIZE)))) +
    theme_nm(legend = "bottom") +
    theme(axis.line.y  = element_blank(),
          axis.ticks.y = element_blank(),
          panel.spacing.y = unit(1.5, "mm"))
  if (show_rows) {
    g + theme(strip.text.y = element_blank())
  } else {
    g + theme(axis.text.y = element_blank(),
              strip.text.y = element_text(size = NM_TXT, angle = 0, hjust = 0, lineheight = 0.95))
  }
}

# =============================================================================
# THE TEXT'S NUMBERS (author, 2026-10-11), asserted where they are computed:
#   "OXPHOS subunit expression fell in both the WT and the MYC+ gland, further than
#    all 2,000 expression-matched gene sets"   -> percentile 0 in both windows
#   "OXPHOS subunits lost more priority than nine in ten mitochondrial pathways"
#   "the carnitine shuttle and fatty acid oxidation ... ranked among the largest
#    gains"                                      -> top fifth of the 143, both windows
#   "OXPHOS subunits fell further than the proliferative signature (mean log2 fold
#    change -0.25 versus -0.05)"                -> WT; paired test below 0.05
#   "although the terminal end bud programme regressed" -> WT, percentile 0
# =============================================================================
wn <- as.data.frame(s43$wt_null); an <- as.data.frame(tv$arm_null)
pn <- as.data.frame(s43$paired_null)
r143 <- ruler[!ruler$is_mtdna, ]
stopifnot(nrow(r143) == 143L)
stopifnot(s43$defs$n_set_draws == 2000L, tv$params$NSET == 2000L,
          wn$percentile[wn$arm == "OXPHOS subunits"] == 0,
          an$pct_wt[an$arm == "OXPHOS subunits"] == 0,
          an$pct_myc[an$arm == "OXPHOS subunits"] == 0)
ox <- r143$pathway == "OXPHOS subunits"
stopifnot(mean(r143$p_tn > r143$p_tn[ox]) >= 0.9, mean(r143$p_tp > r143$p_tp[ox]) >= 0.9)
for (nm in c("Carnitine shuttle", "Fatty acid oxidation")) {
  k <- r143$pathway == nm
  stopifnot(mean(r143$p_tn >= r143$p_tn[k]) <= 0.2, mean(r143$p_tp >= r143$p_tp[k]) <= 0.2)
}
TEXT_WT_EXPR <- c(oxphos = -0.25, prolif = -0.05)
stopifnot(round(rows$expr_neg[rows$label == "OXPHOS subunits"], 2) == TEXT_WT_EXPR[["oxphos"]],
          round(rows$expr_neg[rows$label == "Proliferation"], 2) == TEXT_WT_EXPR[["prolif"]],
          pn$observed_diff[pn$arm_b == "PROLIF_* pooled"] < 0,
          pn$p_emp_lower[pn$arm_b == "PROLIF_* pooled"] < 0.05,
          rows$expr_neg[rows$label == "TEB vs ductal"] < 0,
          wn$percentile[wn$arm == "TEB vs ductal (HS)"] == 0)

# room for a point's radius before the axis rounds up to its end tick
pad  <- function(v) diff(range(v, na.rm = TRUE)) * 0.04
p_ex <- one_ruler(ex, min(ex$val) - pad(ex$val), max(ex$val) + pad(ex$val),
                  "Expression\n(log2 fold change)", 0.2, show_rows = TRUE)
p_pr <- one_ruler(pr, min(pr$val, na.rm = TRUE) - pad(pr$val), max(pr$val, na.rm = TRUE) + pad(pr$val),
                  "Priority\n(mitoPPS)", 0.2, show_rows = FALSE)

save_panel(stack_parts(p_ex, p_pr, ncol = 2, widths = c(1.15, 1), collect_guides = TRUE),
           fig = "Fig1", panel = "f", name = "window_reallocation", width = 100, height = 56)

# --- numbers for the legend (printed, never drawn) --------------------------------
pc <- function(arm) sprintf("WT %.1f (script 43), MYC+ %.1f (script 54)",
                            wn$percentile[wn$arm == arm], an$pct_myc[an$arm == arm])
cat("\nFig. 1f -- for the legend\n")
cat("  12W vs 6W within each genotype; n = 6 mice per group; expression = set-average\n")
cat("  raw (unshrunken) DESeq2 log2 fold change, unweighted; priority = mitoPPS change\n")
for (i in seq_len(nrow(rows))) {
  cat(sprintf("    %-32s expression WT %+.3f / MYC+ %+.3f   priority %s\n", rows$label_n[i],
              rows$expr_neg[i], rows$expr_pos[i],
              if (is.na(rows$prio_neg[i])) "not defined (not mitochondrial)"
              else sprintf("WT %+.3f / MYC+ %+.3f", rows$prio_neg[i], rows$prio_pos[i])))
}
cat("  percentile among 2,000 expression-matched sets:\n")
cat("    OXPHOS subunits:", pc("OXPHOS subunits"), "\n")
cat("    Proliferation:  ", pc("PROLIF_* pooled"), "\n")
cat("    TEB vs ductal:  ", pc("TEB vs ductal (HS)"), "\n")
# priority has no matched null: the per-pathway test of record is script 08's Welch
# test (BH across the MitoPathways), and the rank among the 143 is descriptive
pw   <- as.data.frame(readRDS(here::here("results", "mitopps_scores.rds"))$mitopps_pairwise)
tn08 <- pw[pw$contrast == "Temporal_Myc-", ]       # same quantity as script 40's p_tn
stopifnot(stats::cor(r143$p_tn, tn08$diff[match(r143$pathway, tn08$pathway)]) > 0.999)
cat("  priority, Welch test (script 08; BH across MitoPathways) and rank among 143:\n")
for (nm in mito$pathway) {
  one <- function(ct, v) {
    t <- pw[pw$contrast == ct & pw$pathway == nm, ]
    x <- v[r143$pathway == nm]
    sprintf("%+.3f, P %.2g, adj %.2g, %s", t$diff, t$p_value, t$padj,
            if (x < 0) sprintf("%d lowest", sum(v <= x)) else sprintf("%d highest", sum(v >= x)))
  }
  cat(sprintf("    %-26s WT %s | MYC+ %s\n", nm, one("Temporal_Myc-", r143$p_tn),
              one("Temporal_Myc+", r143$p_tp)))
}
cat(sprintf("  pathways passing BH on the window: WT %d, MYC+ %d of %d (script 08)\n",
            sum(pw$padj[pw$contrast == "Temporal_Myc-"] < 0.05, na.rm = TRUE),
            sum(pw$padj[pw$contrast == "Temporal_Myc+"] < 0.05, na.rm = TRUE),
            sum(pw$contrast == "Temporal_Myc-")))
cat("  BATCH: each window coincides with an extraction batch (see header)\n")

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {
  print(.nm_drawable(stack_parts(p_ex, p_pr, ncol = 2, widths = c(1.15, 1), collect_guides = TRUE)))
  rows
}
