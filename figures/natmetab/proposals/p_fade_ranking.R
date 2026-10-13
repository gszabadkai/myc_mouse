# =============================================================================
# PROPOSAL (letter open) -- which apoptotic genes depart from MYC's overall fade?
# -----------------------------------------------------------------------------
# The author's argument (2026-10-13): if every MYC effect on the apoptotic machinery
# fell with MYC by the same factor, the balance would keep its composition and its
# tilt would shrink by that factor, the dose effect. The composition changes only
# where a gene changes more (or less) than the fade predicts. So rank the pro- and
# anti-apoptotic genes by that departure, and show where Bbc3 sits.
#
# THE MEASURE is script 44's (PART C, the old Fig. 2H scan): the 12W MYC effect
# minus what MYC's overall fade predicts (0.487 x the 6W effect, the rate fitted over
# the genes MYC moves), over its standard error, sqrt(SE12^2 + 0.487^2 SE6^2).
# Script 44 signs it by the direction of the 6W effect (negative = MYC's effect
# shrank faster than the fade). Here it is signed by FUNCTION instead, so that one
# axis reads the same for every gene:
#   positive (toward survival)  a pro-apoptotic transcript LOWER, or an
#                               anti-apoptotic transcript HIGHER, at 12W than the
#                               fade predicts
#   negative (toward death)     the reverse
# For genes in script 44's reporting set this is exactly sign(LFC6) x z_resid
# (asserted). Genes MYC did not move at 6W (|log2FC| < 0.2) are outside that set;
# for them the same formula is computed HERE (their expected 12W effect is ~0, so
# the departure is simply a 12W effect) and they are drawn grey. Before this
# becomes a panel of record, that computation belongs in a numbered script.
#
# THE BAND: the central 95% of the same quantity over script 44's 8,774
# MYC-responsive genes (|departure z| below its 95th percentile) -- what an ordinary
# gene does.
#
# THE ROSTERS, declared before plotting:
#   BCL-2 family  the anti-apoptotic members, the effectors and the BH3-only
#                 proteins, expressed at baseMean >= 20 (script 44's floor);
#                 Hrk, Bcl2a1d and Bcl2l10 fall below it
#   Downstream    script 42's curated execution and IAP genes (the 2G (alt) roster
#                 minus its BCL-2 members)
#
# Reads: results/collapse_module_ownership.rds (script 44: $collapse_genes, $defs)
#        results/interaction_results.rds (script 03: the 6W and 12W MYC effects)
# =============================================================================

source(here::here("figures", "natmetab", "_style.R"))
source(here::here("functions", "reconcile_gene_symbols.R"))
suppressPackageStartupMessages(library(DESeq2))

cmo  <- readRDS(here::here("results", "collapse_module_ownership.rds"))
cg   <- as.data.frame(cmo$collapse_genes)
RATE <- cmo$defs$global_rate_fitted
ir   <- readRDS(here::here("results", "interaction_results.rds"))
d6   <- as.data.frame(ir$myc_6W_raw); d12 <- as.data.frame(ir$myc_12W_raw)
stopifnot(nrow(cg) == 8774L, abs(RATE - 0.4872) < 1e-4)

BCL2_PRO  <- c("Bax", "Bak1", "Bok", "Bid", "Bbc3", "Bcl2l11", "Pmaip1", "Bad", "Bik", "Bmf", "Hrk")
BCL2_ANTI <- c("Bcl2", "Bcl2l1", "Bcl2l2", "Mcl1", "Bcl2a1a", "Bcl2a1b", "Bcl2a1d", "Bcl2l10")
DOWN_PRO  <- c("Htra2", "Diablo", "Cycs", "Apaf1", "Casp9", "Casp3")
DOWN_ANTI <- c("Xiap", "Birc5")
m42 <- as.data.frame(readRDS(here::here("results", "priming_arm_teb.rds"))$machinery)
stopifnot(setequal(c(DOWN_PRO, DOWN_ANTI),
                   m42$gene[grepl("execution|apoptosome|IAP", m42$arm)]))

ens_or_na <- function(g) {
  e <- tryCatch(recon_to_ensembl(g, rownames(d6)), error = function(e) NA_character_)
  if (length(e) != 1L) NA_character_ else unname(e)
}
r <- data.frame(gene  = c(BCL2_PRO, BCL2_ANTI, DOWN_PRO, DOWN_ANTI),
                cls   = c(rep("pro", length(BCL2_PRO)), rep("anti", length(BCL2_ANTI)),
                          rep("pro", length(DOWN_PRO)), rep("anti", length(DOWN_ANTI))),
                group = c(rep("BCL-2 family", length(BCL2_PRO) + length(BCL2_ANTI)),
                          rep("Downstream", length(DOWN_PRO) + length(DOWN_ANTI))),
                stringsAsFactors = FALSE)
r$ens <- vapply(r$gene, ens_or_na, character(1))
r$bm  <- d6[r$ens, "baseMean"]
r <- r[!is.na(r$ens) & !is.na(r$bm) & r$bm >= 20, ]
r$l6  <- d6[r$ens, "log2FoldChange"];  r$s6  <- d6[r$ens, "lfcSE"]
r$l12 <- d12[r$ens, "log2FoldChange"]; r$s12 <- d12[r$ens, "lfcSE"]
r$dz  <- (r$l12 - RATE * r$l6) / sqrt(r$s12^2 + RATE^2 * r$s6^2)   # unsigned departure
r$z   <- ifelse(r$cls == "anti", 1, -1) * r$dz                     # + = toward survival
k     <- match(r$ens, cg$ens)
r$responsive <- !is.na(k)                                          # in script 44's set
r$pct_z <- cg$pct_z[k]
stopifnot(all(abs(sign(r$l6[r$responsive]) * cg$z_resid[k[r$responsive]] -
                  r$dz[r$responsive]) < 1e-9),
          all(abs(r$l6[!r$responsive]) < 0.2))

# the band: what an ordinary MYC-responsive gene does
BAND <- unname(stats::quantile(abs(sign(cg$lfc_6W) * cg$z_resid), 0.95))

# --- the claims, asserted --------------------------------------------------------
fam <- r[r$group == "BCL-2 family", ]
dn  <- r[r$group == "Downstream", ]
stopifnot(fam$gene[which.max(fam$z)] == "Bbc3",              # Bbc3 at the survival end
          identical(fam$gene[abs(fam$z) > BAND], "Bbc3"),    # the only family member outside
          identical(dn$gene[abs(dn$z) > BAND], "Apaf1"))     # downstream: Apaf1 alone outside

# --- the panel -------------------------------------------------------------------
build <- function(d, strips) {
  d$row <- factor(paste(d$group, d$gene), levels = paste(d$group, d$gene)[order(d$group, d$z)])
  lab   <- stats::setNames(d$gene, paste(d$group, d$gene))
  lo <- min(d$z, -BAND); hi <- max(d$z, BAND)
  pad <- 0.04 * (hi - lo)
  ggplot(d, aes(z, row)) +
    annotate("rect", xmin = -BAND, xmax = BAND, ymin = -Inf, ymax = Inf, fill = "grey92") +
    geom_vline(xintercept = 0, linewidth = NM_LINE, colour = box_line) +
    geom_point(aes(shape = cls, colour = responsive), size = 1.5) +
    (if (strips) facet_grid(group ~ ., scales = "free_y", space = "free_y",
                            labeller = as_labeller(c(`BCL-2 family` = "BCL\u00ad2 family",
                                                     Downstream = "Downstream")))) +
    scale_y_discrete(labels = lab) +
    scale_shape_manual(values = c(pro = 16, anti = 15), name = NULL,
                       breaks = c("pro", "anti"),
                       labels = c(pro = "Pro­apoptotic", anti = "Anti­apoptotic")) +
    scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = "grey60"), name = NULL,
                        breaks = c("TRUE", "FALSE"),
                        labels = c(`TRUE` = "MYC effect at 6W", `FALSE` = "No MYC effect at 6W")) +
    scale_x_nm(step = 1, top = hi + pad, bottom = lo - pad, minor = TRUE, labels = lab_signed) +
    labs(x = expression(atop("Departure from MYC's fade at 12W (" * italic(z) * ")",
                             "" %<-% "toward death" ~~~~~~ "toward survival" %->% "")),
         y = NULL) +
    guides(shape = guide_legend(order = 1, override.aes = list(size = 1.4)),
           colour = guide_legend(order = 2, override.aes = list(shape = 16, size = 1.4))) +
    theme_nm(legend = "bottom") +
    theme(axis.text.y = element_text(face = "italic"),
          axis.line.y = element_blank(), axis.ticks.y = element_blank(),
          strip.text.y = element_text(size = NM_TXT, angle = 0, hjust = 0),
          panel.spacing.y = unit(1.5, "mm"),
          legend.box = "vertical", legend.spacing.y = unit(0, "mm"),
          legend.key.size = unit(2.4, "mm"))
}
p_fam <- build(fam, strips = FALSE)
p_all <- build(r, strips = TRUE)

OUTD <- here::here("outputs", "natmetab", "proposals")
save_panel(p_fam, fig = "Fig1", panel = "x", name = "fade_ranking_bcl2", width = 62, height = 60, dir = OUTD)
save_panel(p_all, fig = "Fig1", panel = "x", name = "fade_ranking_all",  width = 66, height = 86, dir = OUTD)

# --- one sheet with both, at 100% ----------------------------------------------------
if (!isTRUE(getOption("natmetab.dry_run"))) {
  mm <- function(v) grid::unit(v, "mm")
  hy <- function(x) gsub("(?<=[A-Za-z0-9])-(?=[A-Za-z0-9])", "­", x, perl = TRUE)
  grDevices::pdf(file.path(OUTD, "Fig1_fade_ranking_sheet.pdf"), width = 210 / 25.4,
                 height = 297 / 25.4, family = NM_FONT, useDingbats = FALSE)
  grid::grid.newpage()
  txt <- function(s, x, y, size = 7, face = "plain")
    grid::grid.text(hy(s), x = mm(x), y = mm(297 - y), just = c("left", "top"),
                    gp = grid::gpar(fontsize = size, fontface = face, lineheight = 1.15))
  put <- function(p, x, y, w, h) print(p, vp = grid::viewport(x = mm(x), y = mm(297 - y),
                                       width = mm(w), height = mm(h), just = c("left", "top")),
                                       newpage = FALSE)
  txt("Ranking the apoptotic genes by their departure from MYC's fade (proposal, at 100%)", 13, 10, 10, "bold")
  txt("A  BCL-2 family (62 x 60 mm)", 13, 22); put(p_fam, 13, 28, 62, 60)
  txt("B  BCL-2 family and downstream genes (66 x 86 mm)", 100, 22); put(p_all, 100, 28, 66, 86)
  txt(paste(c(
    sprintf("Departure: the 12W MYC effect minus %.3f x the 6W effect (the fade over the genes MYC moves), over its", RATE),
    "standard error (script 44). Signed by function: positive = a pro-apoptotic transcript lower, or an anti-apoptotic",
    "transcript higher, at 12W than the fade predicts. Grey band: the central 95% of 8,774 MYC-responsive genes",
    sprintf("(|z| < %.2f). Black: a MYC effect at 6W (|log2FC| >= 0.2); grey: none, so nothing to fade.", BAND),
    "",
    sprintf("BCL-2 family: Bbc3 alone lies outside the band (z %+.2f; 0.5th percentile of the scan).", fam$z[fam$gene == "Bbc3"]),
    sprintf("Downstream: Apaf1 also lies outside it (z %+.2f): its 6W repression deepened at 12W instead of fading.", dn$z[dn$gene == "Apaf1"]),
    sprintf("Htra2 sits on the death side (z %+.2f): it kept more of its strong 6W induction than the fade predicts.", dn$z[dn$gene == "Htra2"])),
    collapse = "\n"), 13, 122)
  invisible(grDevices::dev.off())
}

# --- numbers for the legend (printed, never drawn) ----------------------------------
cat(sprintf("\nfade rate %.4f (script 44, fitted over the genes MYC moves); band |z| < %.2f (95th pct, 8,774 genes)\n", RATE, BAND))
o <- r[order(r$group, -r$z), ]
for (i in seq_len(nrow(o))) cat(sprintf("  %-12s %-6s %-4s 6W %+.2f  12W %+.2f  z %+.2f  %s\n",
  o$group[i], o$gene[i], o$cls[i], o$l6[i], o$l12[i], o$z[i],
  if (o$responsive[i]) sprintf("(scan percentile %.1f)", o$pct_z[i]) else "(no 6W effect; computed here)"))
