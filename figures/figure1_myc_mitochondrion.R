# =============================================================================
# figure1_myc_mitochondrion.R -- MANUSCRIPT FIGURE 1
# "MYC builds a death-competent mitochondrion"
# -----------------------------------------------------------------------------
# The first two beats of the Results paragraph, one claim per panel:
#
#   A  MYC RAISES MITOCHONDRIAL CONTENT. Compartment share across the four groups
#      for the mass markers, the whole nuclear-encoded MitoCarta compartment and
#      nuclear OXPHOS -- with the 13 mtDNA-encoded transcripts as the reference
#      arm that does NOT move. Genotype brackets only: genotype is the clean axis
#      (batch = timepoint), so the time comparison is not drawn here at all.
#   B  BUT NOT UNIFORMLY. Every MitoPathway ranked by the Myc effect on mitoPPS --
#      a within-compartment PRIORITY score. Some arms are promoted and others
#      demoted on top of a broadly positive content change, and the ordering is
#      structured by function, not by size.
#   C  THAT IS A REALLOCATION, AND IT HAS A DIRECTION. Content against priority,
#      one point per pathway. The content-DOWN / priority-UP quadrant is EMPTY:
#      demoted pathways are not switched off, they rise more slowly than the rest
#      and lose share. This is what "reprioritisation" means, drawn.
#   D  AND THE ORGANELLE MYC BUILDS IS DEATH-COMPETENT. The same genotype contrast
#      across a curated roster: import, cristae and OXPHOS up -- and with them
#      HTRA2 and BAX up while BCL-xL goes down. Asset and liability, built in one
#      move.
#
# SCOPE. Every panel is the GENOTYPE contrast, which is the batch-clean axis of
# this design. The developmental axis is Figure 2. Raw (unshrunken) DESeq2 log2FC
# throughout, per CLAUDE.md. n = 6 per group.
#
# Reads (read-only; the author runs scripts 32, 08 and 42 first):
#   results/mito_content_proxies.rds  -- $shares (A)
#   results/mitopps_scores.rds        -- $mitopps_pairwise, $pathway_tier1_map,
#                                        $pathway_levels, $gene_to_pathway (B, C)
#   results/interaction_results.rds   -- $myc_6W_raw (C, D)
#   results/priming_arm_teb.rds       -- $machinery (D)
# =============================================================================

source(here::here("figures", "theme_myc.R"))
source(here::here("functions", "reconcile_gene_symbols.R"))
if (!requireNamespace("patchwork", quietly = TRUE)) stop("figure1 needs patchwork")
if (!requireNamespace("ggrepel", quietly = TRUE))   stop("figure1 needs ggrepel")
if (!requireNamespace("DESeq2", quietly = TRUE))    stop("figure1 needs DESeq2")

out_dir <- here::here("outputs", "figures")

content <- readRDS(here::here("results", "mito_content_proxies.rds"))
mp      <- readRDS(here::here("results", "mitopps_scores.rds"))
ir      <- readRDS(here::here("results", "interaction_results.rds"))
pa      <- readRDS(here::here("results", "priming_arm_teb.rds"))

need <- function(obj, fields, what) {
  miss <- fields[!fields %in% names(obj)]
  if (length(miss)) stop("figure1: ", what, " is missing -> ", paste(miss, collapse = ", "))
}
need(content, "shares", "mito_content_proxies.rds")
need(mp, c("mitopps_pairwise", "pathway_tier1_map", "pathway_levels", "gene_to_pathway"),
     "mitopps_scores.rds")
need(pa, "machinery", "priming_arm_teb.rds")

# --- shared MitoPathway tier palette (identical to fig02 / figS5 / fig03) -----
tier_lv  <- c("Protein import, sorting and homeostasis", "Mitochondrial central dogma",
              "OXPHOS", "Metabolism", "Signaling",
              "Mitochondrial dynamics and surveillance", "Small molecule transport")
tier_lab <- c("Protein import / homeostasis", "Central dogma", "OXPHOS", "Metabolism",
              "Signaling", "Dynamics & surveillance", "SM transport")
names(tier_lab) <- tier_lv
tier_col <- c("Protein import / homeostasis" = "#D55E00", "Central dogma" = "#E69F00",
              "OXPHOS"                  = "#009E73", "Metabolism"    = "grey72",
              "Signaling"               = "#56B4E9",
              "Dynamics & surveillance" = "#0072B2", "SM transport"  = "#CC79A7")

# =============================================================================
# PANEL A -- content rises, and the mtDNA-encoded arm does not
# =============================================================================
arms <- c("MASS_MARKERS_NOCHAP", "MITOCARTA_NUCLEAR_ENCODED",
          "MITOCARTA_OXPHOS_NU", "MITOCARTA_MTDNA_ENCODED")
arm_name <- c(MASS_MARKERS_NOCHAP       = "Mass markers",
              MITOCARTA_NUCLEAR_ENCODED = "Nuclear MitoCarta",
              MITOCARTA_OXPHOS_NU       = "Nuclear OXPHOS",
              MITOCARTA_MTDNA_ENCODED   = "mtDNA-encoded")
if (!all(arms %in% content$shares$panel))
  stop("figure1 A: arms absent from mito_content_proxies$shares -- re-run script 32")

sh <- content$shares[content$shares$panel %in% arms, ]
sh$panel <- factor(sh$panel, levels = arms, labels = arm_name[arms])
sh$group <- factor(sh$group, levels = names(group_labels))

# one 12W WT mtDNA sample sits near 67%; every other point is below 40. Cap the
# DISPLAY so the other three facets stay readable, and mark the capped point.
CAP <- 40
sh$y_disp <- pmin(sh$share_nomt, CAP)
sh$capped <- sh$share_nomt > CAP

# genotype p-values by script 32's own method (simple-effect lm on log2 share),
# so the brackets reconcile with the numbers in the text. ONLY the genotype
# comparison is drawn: batch = timepoint makes the temporal one confounded.
pval_geno <- function(d, tp) {
  d <- d[d$timepoint == tp, ]
  summary(stats::lm(log2(share_nomt) ~ myc_status, data = d))$coefficients[2, "Pr(>|t|)"]
}
fmt_p <- function(p) if (p < 0.001) "<0.001" else formatC(p, format = "g", digits = 2)

brk <- do.call(rbind, lapply(arms, function(a) {
  d  <- content$shares[content$shares$panel == a, ]
  hi <- max(pmin(d$share_nomt, CAP))
  data.frame(panel = factor(arm_name[[a]], levels = arm_name[arms]),
             x1 = c(1, 3), x2 = c(2, 4), tp = c("6W", "12W"),
             p = c(pval_geno(d, "6W"), pval_geno(d, "12W")),
             y = hi * c(1.06, 1.06), tick = hi * 0.014, stringsAsFactors = FALSE)
}))
brk$lab  <- vapply(brk$p, fmt_p, character(1))
brk$xmid <- (brk$x1 + brk$x2) / 2
brk$col  <- ifelse(brk$p < 0.05, "sig", "ns")

hr <- do.call(rbind, lapply(levels(sh$panel), function(f) {
  hi <- max(sh$y_disp[sh$panel == f])
  data.frame(panel = factor(f, levels = levels(sh$panel)), y = hi * 1.14)
}))

pA <- ggplot2::ggplot(sh, ggplot2::aes(group, y_disp)) +
  ggplot2::geom_boxplot(ggplot2::aes(fill = group), outlier.shape = NA,
                        linewidth = 0.3, width = 0.62, colour = "grey25") +
  ggplot2::geom_jitter(width = 0.13, height = 0, size = 0.8, colour = "grey20",
                       alpha = 0.85) +
  ggplot2::geom_point(data = sh[sh$capped, ], shape = 24, size = 1.5, fill = "white",
                      colour = "grey20", stroke = 0.35) +
  ggplot2::geom_segment(data = brk, ggplot2::aes(x = x1, xend = x2, y = y, yend = y),
                        inherit.aes = FALSE, linewidth = 0.3, colour = "grey30") +
  ggplot2::geom_segment(data = brk, ggplot2::aes(x = x1, xend = x1, y = y, yend = y - tick),
                        inherit.aes = FALSE, linewidth = 0.3, colour = "grey30") +
  ggplot2::geom_segment(data = brk, ggplot2::aes(x = x2, xend = x2, y = y, yend = y - tick),
                        inherit.aes = FALSE, linewidth = 0.3, colour = "grey30") +
  ggplot2::geom_text(data = brk, ggplot2::aes(x = xmid, y = y, label = lab, colour = col),
                     inherit.aes = FALSE, vjust = -0.35, size = 1.95) +
  ggplot2::geom_blank(data = hr, ggplot2::aes(x = 1, y = y), inherit.aes = FALSE) +
  ggplot2::facet_wrap(~ panel, nrow = 1, scales = "free_y") +
  ggplot2::scale_fill_manual(values = group_cols, labels = group_labels, name = NULL) +
  ggplot2::scale_colour_manual(values = c(sig = "#B2182B", ns = "grey45"), guide = "none") +
  ggplot2::scale_x_discrete(labels = NULL) +
  ggplot2::labs(
    x = NULL, y = "share of the nuclear transcriptome  (%)",
    title = "A   MYC raises mitochondrial content",
    subtitle = "brackets = the genotype comparison at each age, the batch-clean axis (red = p<0.05).\nThe 13 mtDNA-encoded transcripts are the reference arm and do not move.") +
  theme_myc(base_size = 8) +
  ggplot2::theme(
    legend.position = "bottom",
    legend.key.size = ggplot2::unit(3.2, "mm"),
    legend.text     = ggplot2::element_text(size = 6.4),
    axis.ticks.x    = ggplot2::element_blank(),
    strip.text      = ggplot2::element_text(size = 7, face = "bold"),
    plot.title      = ggplot2::element_text(face = "bold", size = 8.5),
    plot.subtitle   = ggplot2::element_text(size = 6.2, colour = "grey25",
                                            lineheight = 1.15))

# =============================================================================
# PANEL B -- the same effect is NOT uniform across the compartment
# =============================================================================
pw6  <- mp$mitopps_pairwise[mp$mitopps_pairwise$contrast == "Myc_effect_6W",
                            c("pathway", "diff", "padj")]
names(pw6) <- c("pathway", "eff6", "padj6")
B <- pw6
B$tier <- unname(mp$pathway_tier1_map[B$pathway])
B <- B[!is.na(B$tier) & !is.na(B$eff6), ]
B$Tier <- factor(unname(tier_lab[B$tier]), levels = unname(tier_lab[tier_lv]))
B <- B[order(-B$eff6), ]
B$rank <- seq_len(nrow(B))
B$sig  <- !is.na(B$padj6) & B$padj6 < 0.05

l1    <- mp$pathway_levels$Pathway[mp$pathway_levels$Level == "Pathway_Level1"]
ancB  <- B[B$pathway %in% l1, ]
ancB$lab <- unname(tier_lab[unname(mp$pathway_tier1_map[ancB$pathway])])
ancB <- ancB[order(-ancB$eff6), ]

# Seven leader labels do not fit in a third of the page. The tier COLOUR legend
# (collected across the band) names them; the panel annotates only the two
# top-level arms at the extremes, which is the sentence the panel is making.
top_arm <- ancB$lab[1]; bot_arm <- ancB$lab[nrow(ancB)]
PAD <- diff(range(B$eff6)) * 0.12

pB <- ggplot2::ggplot(B, ggplot2::aes(rank, eff6)) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey55") +
  ggplot2::geom_point(ggplot2::aes(colour = Tier), size = 1.2, alpha = 0.9) +
  ggplot2::geom_text(data = B[B$sig, ], ggplot2::aes(y = eff6 + 0.025), label = "*",
                     size = 2.4, colour = "grey25") +
  ggplot2::geom_point(data = ancB, ggplot2::aes(fill = Tier), shape = 21, size = 2.4,
                      colour = "black", stroke = 0.45, show.legend = FALSE) +
  ggplot2::annotate("text", x = 1, y = max(B$eff6) + PAD * 0.75,
                    label = paste("promoted:", top_arm), hjust = 0, size = 2.05,
                    fontface = "bold", colour = tier_col[[top_arm]]) +
  ggplot2::annotate("text", x = nrow(B), y = min(B$eff6) - PAD * 0.75,
                    label = paste("demoted:", bot_arm), hjust = 1, size = 2.05,
                    fontface = "bold", colour = tier_col[[bot_arm]]) +
  ggplot2::scale_colour_manual(values = tier_col, name = NULL, drop = FALSE) +
  ggplot2::scale_fill_manual(values = tier_col, guide = "none") +
  # the band's ONE tier legend lives here. It must be set on this plot and not
  # with patchwork's `&`, which would apply guides() to every panel in the band
  # and resurrect the guides that panels C and D deliberately suppress.
  ggplot2::guides(colour = ggplot2::guide_legend(nrow = 2,
                                                 override.aes = list(size = 2))) +
  ggplot2::coord_cartesian(ylim = range(B$eff6) + c(-PAD, PAD)) +
  ggplot2::labs(
    x = sprintf("all %d MitoPathways, ranked", nrow(B)),
    y = "MYC effect on pathway priority",
    title = "B   Reprioritised, not just enlarged",
    subtitle = sprintf("* = padj<0.05 (%d of %d). Ringed = the seven top-level arms.",
                       sum(B$sig), nrow(B))) +
  theme_myc(base_size = 8) +
  ggplot2::theme(legend.position = "bottom",
                 plot.title    = ggplot2::element_text(face = "bold", size = 8.5),
                 plot.subtitle = ggplot2::element_text(size = 6.2, colour = "grey25",
                                                       lineheight = 1.15))

# =============================================================================
# PANEL C -- content against priority: the reallocation has a direction
# =============================================================================
gp   <- mp$gene_to_pathway
pcol <- names(gp)[1]; gcol <- names(gp)[2]
L6   <- { r <- as.data.frame(ir$myc_6W_raw)
          stats::setNames(r$log2FoldChange, rownames(r)) }
set_lfc <- function(p) {
  e <- recon_to_ensembl(unique(gp[[gcol]][gp[[pcol]] == p]), names(L6))
  v <- L6[e]
  if (length(v)) mean(v, na.rm = TRUE) else NA_real_
}
Cp <- B[, c("pathway", "eff6", "padj6", "Tier")]
Cp$content <- vapply(Cp$pathway, set_lfc, numeric(1))
Cp <- Cp[is.finite(Cp$content), ]
names(Cp)[names(Cp) == "eff6"] <- "priority"
Cp$sig <- !is.na(Cp$padj6) & Cp$padj6 < 0.05
ancC <- Cp[Cp$pathway %in% l1, ]
ancC$lab <- unname(tier_lab[unname(mp$pathway_tier1_map[ancC$pathway])])

# the claim the panel makes: the content-down / priority-up quadrant is empty.
n_q_empty <- sum(Cp$content < 0 & Cp$priority > 0)
n_q_full  <- sum(Cp$content > 0 & Cp$priority < 0)
pct_up    <- round(100 * mean(Cp$content > 0))

pC <- ggplot2::ggplot(Cp, ggplot2::aes(priority, content)) +
  ggplot2::annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = 0,
                    fill = "grey93") +
  ggplot2::geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
  ggplot2::geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
  ggplot2::geom_point(data = Cp[!Cp$sig, ], ggplot2::aes(colour = Tier), size = 1.3,
                      alpha = 0.8) +
  ggplot2::geom_point(data = Cp[Cp$sig, ], ggplot2::aes(fill = Tier), shape = 21,
                      size = 1.7, colour = "black", stroke = 0.35) +
  ggplot2::geom_point(data = ancC, ggplot2::aes(fill = Tier), shape = 21, size = 2.4,
                      colour = "black", stroke = 0.45, show.legend = FALSE) +
  ggplot2::annotate("text", x = Inf, y = -Inf,
                    label = sprintf("  content DOWN and\n  priority UP: %d  ", n_q_empty),
                    hjust = 1, vjust = -0.4, size = 2, colour = "grey35",
                    lineheight = 1.1) +
  ggplot2::scale_colour_manual(values = tier_col, guide = "none", drop = FALSE) +
  ggplot2::scale_fill_manual(values = tier_col, guide = "none") +
  ggplot2::labs(
    x = "priority:  MYC effect on mitoPPS",
    y = "content:  MYC set-average log2FC",
    title = "C   Demoted arms rise more slowly",
    subtitle = sprintf("%d%% rise in content; %d of those lose priority.\nThe reverse quadrant (shaded) is empty.",
                       pct_up, n_q_full)) +
  theme_myc(base_size = 8) +
  ggplot2::theme(legend.position = "none",
                 plot.title    = ggplot2::element_text(face = "bold", size = 8.5),
                 plot.subtitle = ggplot2::element_text(size = 6.2, colour = "grey25",
                                                       lineheight = 1.15))

# =============================================================================
# PANEL D -- the organelle MYC builds is death-competent
# =============================================================================
block_of <- function(arm) {
  ifelse(grepl("^execution", arm),               "execution",
  ifelse(grepl("^effector", arm),                "effector (BAX/BAK)",
  ifelse(grepl("^BH3-only", arm),                "BH3-only trigger",
  ifelse(grepl("^brake", arm),                   "brake",
  ifelse(grepl("^OXPHOS subunit", arm),          "OXPHOS subunit",
  ifelse(grepl("^biogenesis TF|coactivator|mtDNA machinery", arm), "biogenesis TF",
                                                 "import / cristae"))))))
}
BLOCKS <- c("biogenesis TF", "import / cristae", "OXPHOS subunit",
            "brake", "BH3-only trigger", "effector (BAX/BAK)", "execution")
M <- as.data.frame(pa$machinery)
M$block <- factor(block_of(M$arm), levels = BLOCKS)
M$sig6  <- !is.na(M$padj_myc_6W) & M$padj_myc_6W < 0.05
labD <- M$gene %in% c("Htra2", "Bax", "Bcl2l1", "Cycs", "Hspd1", "Tomm22", "Atp5f1a")
jitD <- ggplot2::position_jitter(width = 0, height = 0.16, seed = 11)

pD <- ggplot2::ggplot(M, ggplot2::aes(lfc_myc_6W, block)) +
  ggplot2::geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.35) +
  ggplot2::geom_point(ggplot2::aes(fill = sig6), shape = 21, size = 1.9, stroke = 0.3,
                      colour = "grey35", position = jitD) +
  ggrepel::geom_text_repel(data = M[labD, ], ggplot2::aes(label = gene), size = 2.1,
                           seed = 11, position = jitD, min.segment.length = 0.15,
                           segment.size = 0.25, box.padding = 0.45,
                           point.padding = 0.2, max.overlaps = Inf) +
  ggplot2::scale_fill_manual(values = c(`TRUE` = "#D73027", `FALSE` = "grey82"),
                             guide = "none") +
  ggplot2::labs(
    x = "MYC effect at 6 weeks  (raw log2FC)", y = NULL,
    title = "D   Built to execute death",
    subtitle = "filled = padj<0.05. HTRA2 and BAX rise with\nthe organelle; BCL-xL falls.") +
  theme_myc(base_size = 8) +
  ggplot2::theme(axis.text.y   = ggplot2::element_text(size = 6.4),
                 plot.title    = ggplot2::element_text(face = "bold", size = 8.5),
                 plot.subtitle = ggplot2::element_text(size = 6.2, colour = "grey25",
                                                       lineheight = 1.15))

# =============================================================================
# ASSEMBLY -- A across the top, B/C/D as the lower band
# =============================================================================
bottom <- patchwork::wrap_plots(pB, pC, pD, nrow = 1, widths = c(1.1, 1, 1)) +
  patchwork::plot_layout(guides = "collect") &
  ggplot2::theme(legend.position = "bottom",
                 legend.key.size = ggplot2::unit(3, "mm"),
                 legend.text     = ggplot2::element_text(size = 6),
                 legend.margin   = ggplot2::margin(0, 0, 0, 0))

p <- patchwork::wrap_plots(pA, bottom, ncol = 1, heights = c(1, 1.12)) +
  patchwork::plot_annotation(
    caption = paste(
      "Genotype contrast throughout (Myc+ vs wild-type), which is the batch-clean axis of this design; n=6 per group. Raw, unshrunken DESeq2 log2FC.",
      "Content = share of the nuclear transcriptome, excluding the 13 mtDNA-encoded transcripts from numerator and denominator. Priority = mitoPPS (Monzel 2025), a pairwise",
      "ratio within the mitochondrial compartment and therefore blind to total content: a pathway can rise in content and still lose priority, which is what panel C plots.",
      sep = "\n"),
    theme = ggplot2::theme(
      plot.caption = ggplot2::element_text(size = 5.8, hjust = 0, colour = "grey30",
                                           lineheight = 1.15)))

if (!isTRUE(getOption("myc.fig.nosave"))) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  save_panel(p, file.path(out_dir, "figure1_myc_mitochondrion.pdf"),
             width = fig_w[["double"]], height = 155)
}

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {
  ## panel A: the genotype p-values the brackets carry
  brk[, c("panel", "tp", "p", "lab")] |> print()

  ## panel B: what sits at each end of the ranking
  head(B[, c("pathway", "eff6", "padj6", "Tier")], 8) |> print()
  tail(B[, c("pathway", "eff6", "padj6", "Tier")], 8) |> print()

  ## panel C: the claim, as counts
  cat("content up:", pct_up, "%   content-up/priority-down:", n_q_full,
      "  content-down/priority-up:", n_q_empty, "\n")

  ## panel D: the three genes the caption names
  M[M$gene %in% c("Htra2", "Bax", "Bcl2l1"),
    c("gene", "arm", "lfc_myc_6W", "padj_myc_6W")] |> print()

  print(pA); print(pB); print(pC); print(pD)
  print(p)
}
