# =============================================================================
# edfig1b_myc_effect_rescaled.R -- Extended Data Fig. 1b: the MYC effect at 12
# weeks is the 6-week effect at about half amplitude, gene by gene
# -----------------------------------------------------------------------------
# The text: "The amplitude of the MYC-driven transcriptional program was reduced
# to around half, but it retained its overall composition ..." Half amplitude is
# the slope of the line; the composition is how closely the genes follow it (R2).
# ED Fig. 1c carries the gene-set half: "... without selective repression of
# specific target gene sets".
#
# THE GENES: the 2,648 MYC-responsive genes of script 44's rate of record
#   (IHW-adjusted P < 0.1 at 6W, |log2FC| >= 0.2, baseMean >= 20), reproduced
#   here from its own definition and asserted.
# THE LINE: through the origin, because the quantity is a multiplier ("the 12-week
#   effect is 0.49 times the 6-week effect"); asserted equal to script 44's rate.
# R2: the squared Pearson correlation, not a through-origin lm's uncentred R2.
#
# FORM: two signed effects per gene, 2,648 points -> scatter on equal axes, with
#   the identity line for "no change" (data-to-viz, "Overplotting": small,
#   translucent points).
# WINDOW: drawn over +/-3 log2 so the cloud that carries the fit is readable. The
#   genes outside are not drawn; every number is computed on all of them, and the
#   count left out is printed for the legend.
#
# Reads:  results/collapse_module_ownership.rds (script 44: $collapse_genes, $defs)
# Output: outputs/natmetab/EDFig1/EDFig1b_myc_effect_rescaled.pdf
# =============================================================================

source(here::here("figures", "natmetab", "_style.R"))

cm  <- readRDS(here::here("results", "collapse_module_ownership.rds"))
cg  <- as.data.frame(cm$collapse_genes)
dfs <- cm$defs

in_fit <- !is.na(cg$padj_6W) & cg$padj_6W < dfs$padj6 &
  abs(cg$lfc_6W) >= dfs$lfc6_floor & cg$baseMean >= dfs$basemean_floor
stopifnot(identical(as.logical(in_fit), as.logical(cg$in_ranking)),
          sum(in_fit) == dfs$n_ranking_genes, nrow(cg) == dfs$n_reported_genes)
g <- cg[in_fit, c("gene", "lfc_6W", "lfc_12W")]

fit   <- stats::lm(lfc_12W ~ 0 + lfc_6W, data = g)
slope <- unname(stats::coef(fit))
r2    <- stats::cor(g$lfc_6W, g$lfc_12W)^2
stopifnot(abs(slope - dfs$global_rate_fitted) < 1e-6)

WIN <- 3
shown <- g[pmax(abs(g$lfc_6W), abs(g$lfc_12W)) <= WIN, ]
n_out <- nrow(g) - nrow(shown)

ann <- data.frame(x = -WIN + 0.15, y = WIN - c(0.25, 0.75),
                  lab = c(sprintf('slope~"%.2f"', slope),
                          sprintf('R^displaystyle(2)~"%.2f"', r2)))

p <- ggplot(shown, aes(lfc_6W, lfc_12W)) +
  geom_hline(yintercept = 0, linewidth = NM_LINE, colour = "grey80") +
  geom_vline(xintercept = 0, linewidth = NM_LINE, colour = "grey80") +
  geom_abline(slope = 1, intercept = 0, linetype = "22", linewidth = NM_LINE,
              colour = "grey45") +
  geom_point(size = 0.35, alpha = 0.35, stroke = 0, colour = "grey20") +
  geom_abline(slope = slope, intercept = 0, linewidth = 0.45, colour = "black") +
  geom_text(data = ann, aes(x, y, label = lab), parse = TRUE, hjust = 0,
            size = pt(NM_TXT)) +
  scale_x_nm(step = 1, top = WIN, bottom = -WIN, labels = lab_signed) +
  scale_y_nm(step = 1, top = WIN, bottom = -WIN, labels = lab_signed) +
  coord_equal() +
  labs(x = "MYC+ vs WT, 6W (log2 fold change)",
       y = "MYC+ vs WT, 12W (log2 fold change)") +
  theme_nm()

save_panel(p, fig = "EDFig1", panel = "b", name = "myc_effect_rescaled",
           width = 58, height = 58)

# --- numbers for the legend (printed, never drawn) --------------------------------
free <- stats::lm(lfc_12W ~ lfc_6W, data = g)
allg <- stats::lm(lfc_12W ~ 0 + lfc_6W, data = cg)
cat("\nED Fig. 1b -- for the legend\n")
cat(sprintf("  %s MYC-responsive genes: IHW-adjusted P < %g at 6W, |log2FC| >= %g, baseMean >= %g\n",
            format(nrow(g), big.mark = ","), dfs$padj6, dfs$lfc6_floor, dfs$basemean_floor))
cat(sprintf("  solid line: least squares through the origin, slope %.3f; dashed line: equal effect at both ages\n", slope))
cat(sprintf("  R2 (squared Pearson) %.3f\n", r2))
cat(sprintf("  with a free intercept: slope %.3f, intercept %+.3f\n",
            stats::coef(free)[2], stats::coef(free)[1]))
cat(sprintf("  all %s reported genes: slope %.3f, R2 %.3f\n", format(nrow(cg), big.mark = ","),
            stats::coef(allg), stats::cor(cg$lfc_6W, cg$lfc_12W, use = "complete.obs")^2))
cat(sprintf("  %d gene(s) beyond +/-%g log2 not drawn (included in every number)\n", n_out, WIN))

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {
  print(p)
  g[pmax(abs(g$lfc_6W), abs(g$lfc_12W)) > WIN, ]
}
