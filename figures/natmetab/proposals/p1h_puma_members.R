# =============================================================================
# PROPOSAL -- Fig. 1h: PUMA, BCL-XL, BAX and FOXO3 transcripts in every gland, with
# BMF as the developmental control
# -----------------------------------------------------------------------------
# 2H+I (alt) part A's per-animal form, with 2G (alt)'s content for the four genes
# the text names. Every animal over a box (as Fig. 1b), groups genotype-major so
# each genotype's 6W -> 12W pair sits side by side.
#   lower brackets  12W vs 6W within each genotype: DESeq2 Wald, IHW-adjusted P
#   upper bracket   the difference between the two changes (age x genotype
#                   interaction): raw Wald P. Bbc3 is pre-specified; Bax and
#                   Foxo3 are exploratory (shown and labelled, never quoted);
#                   Bcl2l1 is the fixed denominator of the pre-specified pair.
# Each age was one extraction batch, so each lower bracket coincides with batch;
# the upper bracket is within batch.
# Not bracketed on purpose: the 6W genotype difference for Bbc3 (P 0.064; 3a.3).
#
# BMF IS THE DEVELOPMENTAL CONTROL (author, 2026-10-13). Its ratio to BCL-XL falls
# nearly as far below the eight-ratio line as PUMA's, but the gene rises about
# two-fold across the window in BOTH genotypes and its interaction is inside the
# declared 0.20 (on the diagonal, script 54's rule): a developmental change MYC
# barely modifies, where Bbc3 is flat in WT and falls only under MYC.
#
# Reads: results/dds_int_run.rds, results/interaction_results.rds (script 03);
#        results/two_timeline_verification.rds (script 54, $four_genes: cross-check);
#        results/priming_arm_teb.rds (script 42, $machinery: the BCL-2 family);
#        results/collapse_module_ownership.rds (script 44, $mech_genes: the 26)
# =============================================================================

source(here::here("figures", "natmetab", "_style.R"))
source(here::here("functions", "reconcile_gene_symbols.R"))
suppressPackageStartupMessages(library(DESeq2))

dds <- readRDS(here::here("results", "dds_int_run.rds"))
ir  <- readRDS(here::here("results", "interaction_results.rds"))
fg  <- as.data.frame(readRDS(here::here("results", "two_timeline_verification.rds"))$four_genes)

GENES <- c("Bbc3", "Bcl2l1", "Bax", "Bmf", "Foxo3")
FOUR  <- c("Bbc3", "Bcl2l1", "Bax", "Foxo3")             # script 54's $four_genes
# "\u00ad" (Latin-1 soft hyphen) prints as a true hyphen from base pdf(); the ASCII
# hyphen prints as a minus sign (FIGURE_RULES.md)
PROT  <- c(Bbc3 = "PUMA", Bcl2l1 = "BCL\u00adXL", Bax = "BAX", Bmf = "BMF", Foxo3 = "FOXO3")
ens   <- vapply(GENES, function(g) recon_to_ensembl(g, rownames(dds)), character(1))
stopifnot(identical(unname(ens[FOUR]), fg$ensembl[match(FOUR, fg$gene)]))

cd  <- as.data.frame(colData(dds))
grp <- factor(paste(cd$timepoint, cd$myc_status, sep = "_"), levels = order_by_geno)
stopifnot(!anyNA(grp), all(table(grp) == 6L))
Y <- counts(dds, normalized = TRUE)[ens, , drop = FALSE]

# --- the statistics, as reported (never computed here) ---------------------------
st <- function(slot, col) vapply(ens, function(e) as.data.frame(ir[[slot]])[e, col], numeric(1))
S <- data.frame(gene = GENES,
                wt_lfc = st("timepoint_neg_raw", "log2FoldChange"), wt_padj = st("timepoint_neg_raw", "padj"),
                my_lfc = st("timepoint_pos_raw", "log2FoldChange"), my_padj = st("timepoint_pos_raw", "padj"),
                in_lfc = st("interaction_raw", "log2FoldChange"),  in_p   = st("interaction_raw", "pvalue"),
                in_se  = st("interaction_raw", "lfcSE"))
S4 <- S[match(FOUR, S$gene), ]; k <- match(FOUR, fg$gene)
stopifnot(max(abs(S4$wt_lfc - fg$wt_lfc[k])) < 1e-9, max(abs(S4$my_lfc - fg$myc_lfc[k])) < 1e-9,
          max(abs(S4$in_lfc - fg$int_lfc[k])) < 1e-9, max(abs(S4$in_p - fg$int_p[k])) < 1e-12,
          max(abs(S4$wt_padj - fg$wt_padj_ihw[k])) < 1e-9, max(abs(S4$my_padj - fg$myc_padj_ihw[k])) < 1e-9)

# the drawn group means give the reported changes (as Fig. 1b)
gm <- t(apply(Y, 1, function(v) tapply(v, grp, mean)))
stopifnot(max(abs(log2(gm[, "12W_neg"] / gm[, "6W_neg"]) - S$wt_lfc)) < 0.1,
          max(abs(log2(gm[, "12W_pos"] / gm[, "6W_pos"]) - S$my_lfc)) < 0.1)

# --- the sentences, asserted ---------------------------------------------------------
g <- function(x, col) S[[col]][S$gene == x]
stopifnot(g("Bbc3", "my_padj") < 0.05, g("Bbc3", "wt_padj") > 0.5, g("Bbc3", "in_p") < 0.01,
          g("Bcl2l1", "wt_padj") > 0.5, g("Bcl2l1", "my_padj") > 0.5,
          g("Bax", "my_lfc") < 0, g("Bax", "my_padj") < 0.05, g("Bax", "my_lfc") > g("Bbc3", "my_lfc"),
          g("Foxo3", "wt_lfc") > 0, g("Foxo3", "wt_padj") < 0.05, g("Foxo3", "my_padj") > 0.5)
# "of the BCL-2 family, Bbc3 fell furthest in the MYC+ gland" (script 42's table)
m  <- as.data.frame(readRDS(here::here("results", "priming_arm_teb.rds"))$machinery)
FAMILY <- c("Bax", "Bak1", "Bid", "Bbc3", "Bcl2l11", "Bmf", "Pmaip1", "Bcl2l1", "Mcl1", "Bcl2")
fam <- m[match(FAMILY, m$gene), ]
stopifnot(!anyNA(fam$gene), fam$gene[which.min(fam$lfc_myc_time)] == "Bbc3")
# Bmf, the developmental control: the same values as script 42's table, up in both
# genotypes (IHW P < 0.05 each) and on the diagonal by script 54's declared rule
# (|interaction| < 0.20 and raw P > 0.05)
mb <- m[m$gene == "Bmf", ]
stopifnot(abs(mb$lfc_wt_time - g("Bmf", "wt_lfc")) < 1e-9, abs(mb$lfc_myc_time - g("Bmf", "my_lfc")) < 1e-9,
          g("Bmf", "wt_lfc") > 0, g("Bmf", "wt_padj") < 0.05, g("Bmf", "my_lfc") > 0, g("Bmf", "my_padj") < 0.05,
          abs(g("Bmf", "in_lfc")) < 0.20, g("Bmf", "in_p") > 0.05)
# "Foxo3 ... separating the genotypes more sharply than any of the 26" (script 44's set)
mg  <- as.data.frame(readRDS(here::here("results", "collapse_module_ownership.rds"))$mech_genes)
e26 <- vapply(mg$gene, function(x) recon_to_ensembl(x, rownames(dds)), character(1))
z26 <- as.data.frame(ir$interaction_raw)[e26, "stat"]
stopifnot(nrow(mg) == 26L, !anyNA(z26), mg$gene[which.max(abs(z26))] == "Foxo3")

# --- one part per gene ------------------------------------------------------------------
nice_step <- function(top) {
  s <- c(1, 2, 2.5, 5) * 10^floor(log10(top))
  s <- c(s, s * 10)
  s[which(top / s <= 3.2)[1]]
}
# Bracket spacing is set in millimetres of the printed part, so a P label clears the
# bracket above it whatever the gene's scale: PLOT_MM is the data area's height.
PANEL_W <- 96; PANEL_H <- 64
PLOT_MM <- c(top = 19, bottom = 19)
STEP_MM <- 3.4
part <- function(gene, left, bottom) {
  d   <- data.frame(group = grp, y = Y[ens[[gene]], ])
  d$lab <- sprintf('italic("%s")~"(%s)"', gene, PROT[[gene]])
  top <- max(d$y)
  # A P below 0.001 prints with an exponent and is too wide to sit beside its
  # neighbour, so the two lower brackets are staggered for that gene (Bmf).
  pw   <- c(g(gene, "wt_padj"), g(gene, "my_padj"))
  wide <- any(pw < 0.001)
  nlev <- if (wide) 3 else 2
  fr  <- STEP_MM / PLOT_MM[[if (bottom) "bottom" else "top"]]
  T0  <- 1.04 * top / (1 - nlev * fr - 0.03)
  stp <- fr * T0
  win <- data.frame(group1 = c("6W_neg", "6W_pos"), group2 = c("12W_neg", "12W_pos"),
                    p = pw, level = if (wide) c(1, 2) else c(1, 1))
  br1 <- brackets(win, base = 1.04 * top, step = stp)
  br2 <- brackets(data.frame(group1 = 1.5, group2 = 3.5, p = g(gene, "in_p"), level = 1),
                  base = 1.04 * top + (nlev - 1) * stp, step = stp)
  ggplot(d, aes(group, y)) +
    geom_boxplot(aes(fill = group), outlier.shape = NA, width = 0.62,
                 colour = box_line, linewidth = NM_LINE) +
    ggbeeswarm::geom_quasirandom(aes(colour = group), width = 0.18, size = 0.7, shape = 16) +
    br1 + br2 +
    facet_wrap(~ lab, labeller = label_parsed) +
    scale_colour_manual(values = group_cols, guide = "none") +
    scale_fill_manual(values = group_fill, guide = "none") +
    scale_x_discrete(labels = group_labels) +
    scale_y_nm(step = nice_step(attr(br2, "top")), top = attr(br2, "top"), minor = TRUE) +
    labs(x = NULL, y = if (left) "Normalized counts" else NULL) +
    coord_cartesian(clip = "off") +          # a wide P label may overhang the y axis
    theme_nm() +
    # group labels at the floor (5 pt): at 6 pt "MYC+ 6W" and "MYC+ 12W" touch at this width
    theme(axis.text.x = if (bottom) element_text(size = NM_TXT_MIN, lineheight = 0.85) else element_blank(),
          strip.text = element_text(size = NM_TXT_TITLE, margin = margin(0, 0, 0.5, 0, "mm")),
          plot.margin = margin(0.5, 1.5, 0.5, 0.5, "mm"))
}
# row 1: the PUMA:BCL-XL pair and BAX; row 2: BMF (the developmental control) and FOXO3
p <- stack_parts(part("Bbc3", TRUE, TRUE), part("Bcl2l1", FALSE, TRUE), part("Bax", FALSE, TRUE),
                 part("Bmf", TRUE, TRUE),  part("Foxo3", FALSE, TRUE), ncol = 3)

save_panel(p, fig = "Fig1", panel = "h", name = "puma_members", width = PANEL_W, height = PANEL_H,
           dir = here::here("outputs", "natmetab", "proposals"))

cat("\n1h -- for the legend (DESeq2 Wald, raw MLE log2FC; lower brackets IHW-adjusted, upper raw)\n")
for (i in seq_len(nrow(S))) cat(sprintf(
  "  %-6s 12W vs 6W: WT %+.3f (IHW P %.2g)  MYC+ %+.3f (IHW P %.2g)  interaction %+.3f (SE %.3f, raw P %.2g)\n",
  S$gene[i], S$wt_lfc[i], S$wt_padj[i], S$my_lfc[i], S$my_padj[i], S$in_lfc[i], S$in_se[i], S$in_p[i]))
cat("  BCL-2 family, MYC+ 12W vs 6W:", paste(sprintf("%s %+.2f", fam$gene, fam$lfc_myc_time)[order(fam$lfc_myc_time)], collapse = ", "), "\n")
cat(sprintf("  Foxo3: largest |interaction z| of script 44's 26 (%.2f; next %s %.2f)\n",
            max(abs(z26)), mg$gene[order(-abs(z26))][2], sort(abs(z26), decreasing = TRUE)[2]))
