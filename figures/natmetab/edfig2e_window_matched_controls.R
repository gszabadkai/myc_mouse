# =============================================================================
# edfig2e_window_matched_controls.R -- Extended Data Fig. 2e: across the window,
# OXPHOS falls further than comparably expressed genes, and further than
# proliferation
# -----------------------------------------------------------------------------
# The text: "In the WT gland, OXPHOS subunits fell further than all 2,000
# expression-matched gene sets, whereas the proliferative signature fell only
# slightly ... (ED Fig. 2e)." Two stacked parts, one PDF:
#
#   TOP     each programme's 12W vs 6W change in both genotypes, against 2,000
#           random gene sets of the same size drawn from the same expression
#           ventiles (script 54). Filled: beyond the central 95% of those sets;
#           open: within it. The grey tick is the sets' median, which is where a
#           comparably expressed programme would sit.
#   BOTTOM  the comparison itself, WT window (script 43): OXPHOS subunits minus
#           each comparator, against the same difference between paired matched
#           sets redrawn together (grey tick: the paired-null median).
#
# One x axis for both parts. SIZE: 95 x 66 mm, beside ED 2f (88 mm) in one row.
#
# WHY TWO SCRIPTS' NULLS. Script 54 ran the matched draws for both windows on the
# same seed, so the top part uses it throughout; the paired comparison exists only
# in script 43 (WT window). The two WT percentiles agree within 3.2 points
# (verification note, section 1.7); the text quotes neither.
#
# FORM: signed changes per programme, two windows -> dot plot, vertically offset
#   (data-to-viz, "Order your data").
#
# Reads:  results/two_timeline_verification.rds    (script 54: $arm_null)
#         results/substrate_specificity_tradeoff.rds (script 43: $paired_null)
# Output: outputs/natmetab/EDFig2/EDFig2e_window_matched_controls.pdf
# =============================================================================

source(here::here("figures", "natmetab", "_style.R"))

an <- as.data.frame(readRDS(here::here("results", "two_timeline_verification.rds"))$arm_null)
pn <- as.data.frame(readRDS(here::here("results", "substrate_specificity_tradeoff.rds"))$paired_null)

# "OXPHOS (all)" is subunits plus assembly factors, both drawn: left out as redundant
lab <- c("OXPHOS subunits" = "OXPHOS subunits", "OXPHOS assembly" = "OXPHOS assembly factors",
         "nucleotide metabolism" = "Nucleotide metabolism", "mitoribosome" = "Mitochondrial ribosome",
         "TCA cycle" = "TCA cycle", "amino-acid metabolism" = "Amino acid metabolism",
         "lipid metabolism" = "Lipid metabolism", "TEB vs ductal (HS)" = "TEB vs ductal",
         "PROLIF_* pooled" = "Proliferation")
stopifnot(setequal(c(names(lab), "OXPHOS (all)"), an$arm))
an <- an[an$arm %in% names(lab), ]
an$label <- sprintf("%s (%d)", unname(lab[an$arm]), an$n_matched)
an$label <- factor(an$label, levels = an$label[order(-an$obs_wt)])   # most negative at top

beyond <- function(pct) ifelse(pct <= 2.5 | pct >= 97.5, "beyond", "within")
top <- rbind(
  data.frame(label = an$label, geno = "neg", obs = an$obs_wt,  null = an$null_med_wt,
             pos = beyond(an$pct_wt)),
  data.frame(label = an$label, geno = "pos", obs = an$obs_myc, null = an$null_med_myc,
             pos = beyond(an$pct_myc)))
top$geno <- factor(top$geno, levels = c("neg", "pos"))
top$pos  <- factor(top$pos, levels = c("beyond", "within", "median"))

# --- the paired comparison, WT window ------------------------------------------------
# Row label "OXPHOS - proliferation": base pdf() draws the ASCII hyphen as a minus
# sign, which is what these labels mean (FIGURE_RULES.md, Type).
cmp <- c("PROLIF_* pooled" = "proliferation", "mitoribosome" = "mitochondrial ribosome",
         "nucleotide metabolism" = "nucleotide metabolism")
stopifnot(all(pn$arm_a == "OXPHOS subunits"), setequal(names(cmp), pn$arm_b))
pn$label <- factor(sprintf("OXPHOS - %s", cmp[pn$arm_b]),
                   levels = rev(sprintf("OXPHOS - %s", cmp)))
# P from the empirical lower tail of 2,000 paired draws; 0 of 2,000 -> P < 0.0005
pn$p_lab <- ifelse(pn$p_emp_lower == 0, "italic(P)<'0.0005'",
                   sprintf("italic(P)=='%s'", formatC(pn$p_emp_lower, format = "fg", digits = 2)))

# ONE x axis for both parts, so their zero lines align (patchwork aligns panel edges)
xr <- range(c(top$obs, top$null, pn$observed_diff, pn$null_median, 0)); xpad <- 0.04 * diff(xr)
x0 <- xr[1] - xpad; x1 <- xr[2] + xpad
p_lab_x <- 0.06
x_axis  <- function() scale_x_nm(step = 0.2, top = x1, bottom = x0, minor = TRUE, labels = lab_signed)

dodge <- position_dodge(width = 0.64, orientation = "y", reverse = TRUE)
PT_SIZE <- c(neg = 1.3, pos = 0.95)    # MYC+ smaller, as on Fig. 1f (author, 2026-10-11)
p_top <- ggplot(top, aes(obs, label, colour = geno)) +
  geom_vline(xintercept = 0, linewidth = NM_LINE, colour = box_line) +
  geom_point(aes(x = null, shape = "median"), position = dodge, size = 2, colour = "grey55") +
  geom_point(aes(shape = pos, size = geno), position = dodge, stroke = 0.45) +
  scale_colour_manual(values = geno_cols, labels = dev_labels, name = NULL) +
  scale_size_manual(values = PT_SIZE, guide = "none") +
  scale_shape_manual(values = c(beyond = 16, within = 1, median = 124), name = NULL, drop = FALSE,
                     labels = c(beyond = "beyond 95% of matched sets",
                                within = "within 95% of matched sets",
                                median = "median of matched sets")) +
  x_axis() +
  labs(x = "12W vs 6W (log2 fold change)", y = NULL) +
  guides(colour = guide_legend(order = 1, ncol = 1, override.aes = list(size = unname(PT_SIZE),
                                                                         shape = 16)),
         shape  = guide_legend(order = 2, ncol = 1,
                               override.aes = list(size = c(1.5, 1.5, 2), stroke = 0.45,
                                                   colour = c("black", "black", "grey55")))) +
  theme_nm(legend = "bottom") +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        legend.key.height = unit(2.6, "mm"))

p_bot <- ggplot(pn, aes(observed_diff, label)) +
  geom_vline(xintercept = 0, linewidth = NM_LINE, colour = box_line) +
  geom_segment(aes(x = null_median, xend = observed_diff, yend = label),
               linewidth = 0.4, colour = "grey70") +
  geom_point(aes(x = null_median), shape = 124, size = 2, colour = "grey55") +
  geom_point(size = 1.3, colour = geno_cols[["neg"]]) +
  geom_text(aes(x = p_lab_x, label = p_lab), parse = TRUE, hjust = 0, size = pt(NM_TXT)) +
  x_axis() +
  labs(x = "Difference, 12W vs 6W, WT (log2 fold change)", y = NULL) +
  theme_nm() +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank())

save_panel(stack_parts(p_top, p_bot, heights = c(nlevels(top$label), nrow(pn)), collect_guides = TRUE),
           fig = "EDFig2", panel = "e", name = "window_matched_controls", width = 95, height = 66)

# --- numbers for the legend (printed, never drawn) --------------------------------
cat("\nED Fig. 2e -- for the legend\n")
cat("  TOP: set-average raw log2 fold change, 12W vs 6W; 2,000 random sets of the same size\n")
cat("  from the same expression ventiles per programme (script 54); percentile of the observed value:\n")
for (i in order(an$obs_wt)) cat(sprintf("    %-32s WT %+.3f (pct %.1f)   MYC+ %+.3f (pct %.1f)\n",
                                        an$label[i], an$obs_wt[i], an$pct_wt[i], an$obs_myc[i], an$pct_myc[i]))
cat("  BOTTOM: WT window, OXPHOS subunits minus comparator, against 2,000 paired matched draws (script 43):\n")
for (i in seq_len(nrow(pn))) cat(sprintf("    %-40s %+.3f (null median %+.4f), P %s\n", pn$label[i],
                                         pn$observed_diff[i], pn$null_median[i],
                                         if (pn$p_emp_lower[i] == 0) "< 0.0005"
                                         else formatC(pn$p_emp_lower[i], format = "fg", digits = 2)))
cat("  BATCH: each window coincides with an extraction batch; matched sets share it\n")

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {
  print(.nm_drawable(stack_parts(p_top, p_bot, heights = c(nlevels(top$label), nrow(pn)),
                                 collect_guides = TRUE)))
  an
}
