# =============================================================================
# edfig2d_myc_reallocation_ages.R -- Extended Data Fig. 2d: MYC re-weights the
# compartment in the same pattern at both ages, over all 143 MitoPathways
# -----------------------------------------------------------------------------
# The text: "MYC promoted biogenesis and biosynthetic pathways and demoted
# mitochondrial dynamics, surveillance and signalling, in the same pattern at both
# ages (Fig. 1e, ED Fig. 2d)." Fig. 1e draws eleven rows; this is all of them.
# It is ED Fig. 1b's grammar on the priority ruler: the 12-week effect against the
# 6-week effect, one point per pathway.
#
# THE LINE: script 40's fit, lm(12W ~ 6W) over the 143 nuclear MitoPathways (the
#   synthetic mtDNA-encoded pathway excluded, as in every script-40 fit); asserted
#   against a refit. Dashed: equal effect at both ages.
#
# FORM: two signed effects per pathway -> scatter on equal axes; Fig. 1e's eleven
#   rows drawn in black, the rest grey. Not named: see the note at fig1e_rows.
#
# SIZE: 62 x 62 mm, so that ED 2d and ED 2g's two scatters fill one 183-mm row.
#
# Reads:  results/background_vs_myc.rds (script 40: $ruler, $regressions)
# Output: outputs/natmetab/EDFig2/EDFig2d_myc_reallocation_ages.pdf
# =============================================================================

source(here::here("figures", "natmetab", "_style.R"))

bg    <- readRDS(here::here("results", "background_vs_myc.rds"))
ruler <- as.data.frame(bg$ruler)
reg   <- as.data.frame(bg$regressions)
r143  <- ruler[!ruler$is_mtdna, ]
stopifnot(nrow(r143) == 143L)

fit  <- stats::lm(p_m12 ~ p_m6, data = r143)
row2 <- reg[reg$model == "rescale (Myc@12W ~ Myc@6W)" & reg$scope == "priority, all", ]
stopifnot(nrow(row2) == 1L,
          abs(stats::coef(fit)[2] - row2$slope) < 1e-8,
          abs(summary(fit)$r.squared - row2$r2) < 1e-8)

# Fig. 1e's rows, drawn black and not named: at this size eleven names cannot be
# placed off the cloud (tried 2026-10-11), and Fig. 1e names them at both ages.
fig1e_rows <- c("Translation", "Transcription", "Chaperones", "Serine metabolism",
                "Folate and 1-C metabolism", "OXPHOS subunits",
                "Mitochondrial dynamics and surveillance", "Fission", "Apoptosis",
                "Calcium homeostasis", "SLC25A family")
stopifnot(all(fig1e_rows %in% r143$pathway))
hl <- r143[r143$pathway %in% fig1e_rows, ]

pad <- 0.04 * diff(range(c(r143$p_m6, r143$p_m12)))
lo  <- min(c(r143$p_m6, r143$p_m12)) - pad
hi  <- max(c(r143$p_m6, r143$p_m12)) + pad
ann <- data.frame(x = lo + 0.02, y = hi - c(0.02, 0.075),
                  lab = c(sprintf('slope~"%.2f"', stats::coef(fit)[2]),
                          sprintf('R^displaystyle(2)~"%.2f"', summary(fit)$r.squared)))
xs  <- seq(lo, hi, length.out = 2)

p <- ggplot(r143, aes(p_m6, p_m12)) +
  geom_hline(yintercept = 0, linewidth = NM_LINE, colour = "grey80") +
  geom_vline(xintercept = 0, linewidth = NM_LINE, colour = "grey80") +
  geom_abline(slope = 1, intercept = 0, linetype = "22", linewidth = NM_LINE, colour = "grey45") +
  geom_point(size = 0.8, colour = "grey60", stroke = 0) +
  geom_line(data = data.frame(x = xs, y = stats::coef(fit)[1] + stats::coef(fit)[2] * xs),
            aes(x, y), linewidth = 0.45, colour = "black") +
  geom_point(data = hl, size = 1.1, colour = "black", stroke = 0) +
  geom_text(data = ann, aes(x, y, label = lab), parse = TRUE, hjust = 0, size = pt(NM_TXT)) +
  scale_x_nm(step = 0.2, top = hi, bottom = lo, minor = TRUE, labels = lab_signed) +
  scale_y_nm(step = 0.2, top = hi, bottom = lo, minor = TRUE, labels = lab_signed) +
  coord_equal() +
  labs(x = "MYC+ vs WT, 6W (mitoPPS)", y = "MYC+ vs WT, 12W (mitoPPS)") +
  theme_nm()

save_panel(p, fig = "EDFig2", panel = "d", name = "myc_reallocation_ages", width = 62, height = 62)

# --- numbers for the legend (printed, never drawn) --------------------------------
cat("\nED Fig. 2d -- for the legend\n")
cat(sprintf("  %d nuclear MitoPathways; mitoPPS, MYC+ vs WT; n = 6 mice per group\n", nrow(r143)))
cat(sprintf("  solid line: lm(12W ~ 6W), slope %.3f, intercept %+.3f, R2 %.3f; r = %.2f; dashed: equal effect\n",
            stats::coef(fit)[2], stats::coef(fit)[1], summary(fit)$r.squared, stats::cor(r143$p_m6, r143$p_m12)))
cat(sprintf("  same sign at both ages: %.0f%% of pathways; black: the %d rows of Fig. 1e:\n",
            100 * mean(sign(r143$p_m6) == sign(r143$p_m12)), nrow(hl)))
cat("   ", paste(hl$pathway[order(-hl$p_m6)], collapse = ", "), "\n")
top <- r143[which.max(r143$p_m12), ]
cat(sprintf("  top point: %s (%d genes)\n", top$pathway, top$n_genes))

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {
  print(p)
  r143[order(-abs(r143$p_m12 - r143$p_m6)), c("pathway", "p_m6", "p_m12")][1:10, ]
}
