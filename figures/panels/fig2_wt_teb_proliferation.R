# =============================================================================
# fig2_wt_teb_proliferation.R -- what the normal gland loses between six and
# twelve weeks, and what it keeps
# -----------------------------------------------------------------------------
# SLOT: Fig. 2E.
#
#   "The 6>12W_wt comparison showed that while the TEB signature was lost as
#    expected in the puberty-adult transition, the overall proliferation
#    signalling remained relatively stable (Fig. 2E)."
#
# THE JOB. This is the first panel of section 2 and it establishes the substrate:
# what the wild-type mammary epithelium withdraws from over the pubertal-to-adult
# window is a MORPHOGENETIC programme, not a proliferative one. Everything after
# it depends on that -- the MYC-ER animals never saw Myc, so the change that
# licenses the whole of Figure 2 has to be developmental. It is also the negative
# limb that keeps Fig. 2F honest: had proliferation collapsed here, the
# respiratory withdrawal of 2F would be its shadow rather than a finding.
#
# THE RULER IS THE fGSEA NES on `timepoint_neg` (= 6>12W_wt), which is the same
# instrument Figs. 1D and 1H use, so "enrichment" means one thing across the
# figure set. ONE RANKED LIST, SO THE NORMALISER IS SHARED: the comparison this
# panel makes is WITHIN a single ranking, which is safe in exactly the way Fig.
# 1H's comparison ACROSS two rankings was not (there the flatter twelve-week list
# shrinks the normaliser and lifts every NES). Said again in the legend block.
#
# THREE ROWS, and the third is not decoration:
#   TEB programmes (38)          median NES -1.60, 29 significant, ALL depleted
#   ductal, TEB-down (3)         median NES +1.06 -- the direction control
#   TEB x MitoCarta lanes (12)   median NES +0.49, none significant
#   proliferation  (14)          median NES -0.89, NONE significant
#
# TWO OF THOSE ROWS EXIST BECAUSE A FAMILY WITH A SIGN BUILT INTO IT CANNOT SHARE
# A ROW WITH ONE THAT HAS NOT. The `*_VS_DUCTAL_*_DN` sets are the genes lower in
# the end bud than in the duct, so pooling them with the TEB programmes would
# cancel part of the effect being measured; on their own row they are the control
# a directional claim needs, and they rise. The `*_TEB_MITO` lanes are MitoCarta
# subsets by construction (17-59 genes each) and cannot be evidence about the TEB
# programme -- Fig. 1D's rule. Deleting either row would flatter the result and
# folding either into the TEB row would blur it.
#
# WHAT THE PANEL DELIBERATELY DOES NOT DRAW, because two other rulers say the same
# thing and belong in the text: the gene-level set-average raw log2FC with its
# expression-matched null (TEB -0.4220 at the 0th percentile of 2000 draws;
# PROLIF_* pooled -0.0472 at the 1.3rd), and the purity-adjusted per-sample
# composite (-0.4200, p 0.0053). All three are computed and asserted below.
#
# Reads (read-only, no re-run):
#   results/fgsea_percategory.rds              (script 20) -- the drawn ruler
#   results/substrate_specificity_tradeoff.rds (script 43) -- $comparator, $wt_null,
#                                                             $defs (rosters, guard)
#   results/priming_arm_teb.rds                (script 42) -- $teb_signatures
# Output: outputs/figures/panels/fig2_wt_teb_proliferation.pdf
# =============================================================================

source(here::here("figures", "panels", "_panel_common.R"))
if (!requireNamespace("ggrepel", quietly = TRUE)) stop("fig2E needs ggrepel")

fg_path <- here::here("results", "fgsea_percategory.rds")
ss_path <- here::here("results", "substrate_specificity_tradeoff.rds")
pa_path <- here::here("results", "priming_arm_teb.rds")
require_fresher_than(fg_path)      # 12:37, newer than gsva_scores.rds at 12:36
require_fresher_than(ss_path)
require_fresher_than(pa_path)

fg <- as.data.frame(readRDS(fg_path)$fgsea)
ss <- readRDS(ss_path)
pa <- readRDS(pa_path)
stopifnot(all(c("ranking", "category", "pathway", "NES", "pval",
                "padj_within_category", "size") %in% names(fg)))

# Scripts 42 and 43 save TIBBLES, and a tibble's `[` returns a tibble rather than
# a scalar -- which silently poisons every sprintf downstream (Fig. 1E's trap).
cmp <- as.data.frame(ss$comparator)
wtn <- as.data.frame(ss$wt_null)
tsg <- as.data.frame(pa$teb_signatures)

# --- guard: a stale results object must fail here, not at review --------------
# Script 43's own positive control, reused verbatim from figures/fig04: the
# wild-type OXPHOS-subunit value reproduces Issue #4's -0.2548. If membership
# resolution has drifted, nothing in this panel is comparable with the committed
# numbers.
ox_obs <- cmp$c_wt_time[cmp$arm == "OXPHOS subunits"]
if (!isTRUE(abs(ox_obs - ss$defs$oxphos_wt_reference) < 0.02))
  stop(sprintf("fig2E: stale results object -- OXPHOS-subunit wild-type %.4f against the reference %.4f. Re-run script 43.",
               ox_obs, ss$defs$oxphos_wt_reference))

# =============================================================================
# the three families
# =============================================================================
# `programme_group()` is deliberately NOT the grouping here: it would scatter the
# TEB lanes across five of its rows (mammary development, TF target sets, curated
# mitochondrial ...) and the sentence is about the TEB context, which cuts across
# them. The rosters are therefore named here -- and the proliferation one is the
# ANALYSIS OF RECORD's, asserted identical to script 43's `defs$prolif_sets`, so
# the panel and the arm-level statistic quoted in the legend describe one set.
CONTRAST <- "timepoint_neg"          # = 6>12W_wt, the declared contrast vocabulary
ANCHOR   <- "MG_TEB_VS_DUCTAL_HS_GRAY_UP"

tn <- fg[fg$ranking == CONTRAST, ]
stopifnot(nrow(tn) == 866L)

# FOUR rows, and two of them exist because a family with a sign built into it
# cannot share a row with one that does not.
#
#   *_VS_DUCTAL_*_GRAY_DN holds the genes LOWER in the terminal end bud than in
#   the mature duct -- the ductal genes. Their sign is inverted BY CONSTRUCTION,
#   so pooling them with the TEB programmes would cancel part of the very effect
#   the panel measures. Given their own row they become the internal control the
#   claim needs: "lost as expected in the puberty-adult transition" is a
#   DIRECTIONAL statement, and a panel showing only sets going down could be
#   showing a ranking that drifts down. These rise. (Same role the OXPHOS
#   assembly factors play in figures/fig04_substrate_specificity.R panel A.)
#
#   *_MITO lanes are MitoCarta subsets by construction -- Fig. 1D's rule.
fam_levels <- c("proliferation", "TEB x MitoCarta", "ductal (TEB-down)",
                "TEB programmes")
d <- tn[grepl("TEB", tn$pathway) | tn$pathway %in% ss$defs$prolif_sets, ]
d$fam <- ifelse(d$pathway %in% ss$defs$prolif_sets, "proliferation",
                ifelse(grepl("_MITO$", d$pathway), "TEB x MitoCarta",
                       ifelse(grepl("_DN$", d$pathway), "ductal (TEB-down)",
                              "TEB programmes")))
d$fam <- factor(d$fam, levels = fam_levels)
d$sig <- d$padj_within_category < 0.05

n_fam <- table(d$fam)
stopifnot(
  # the proliferation roster IS script 43's, not a list re-typed here
  setequal(d$pathway[d$fam == "proliferation"], ss$defs$prolif_sets),
  n_fam[["TEB programmes"]]    == 38L,
  n_fam[["ductal (TEB-down)"]] ==  3L,
  n_fam[["TEB x MitoCarta"]]   == 12L,
  n_fam[["proliferation"]]     == 14L,
  # every TEB programme is depleted, which is what puts the key top right
  max(d$NES[d$fam == "TEB programmes"]) < -1,
  ANCHOR %in% d$pathway,
  !anyNA(d$NES))

# Row labels carry n. Significance is a SHAPE (Fig. 1F's declared encoding) and
# never a fill or a point colour, so the counts stay off the page and go to the
# legend block.
lab_of <- vapply(fam_levels, function(f) sprintf("%s (%d)", f, n_fam[[f]]),
                 character(1))
d$row <- factor(unname(lab_of[as.character(d$fam)]),
                levels = unname(lab_of[fam_levels]))

med <- data.frame(row = levels(d$row),
                  m   = tapply(d$NES, d$row, stats::median)[levels(d$row)],
                  stringsAsFactors = FALSE)
med$row <- factor(med$row, levels = levels(d$row))

# =============================================================================
# the panel
# =============================================================================
XR  <- range(d$NES) + c(-1, 1) * diff(range(d$NES)) * 0.06
JIT <- ggplot2::position_jitter(width = 0, height = 0.27, seed = 3)
HW  <- 0.36                                   # half-width of a median rule

p <- ggplot2::ggplot(d, ggplot2::aes(NES, row)) +
  # +/-1 is the scale of an UNMOVED set: fGSEA divides the enrichment score by the
  # mean of the same-signed permutation null, so a set that is nowhere in
  # particular lands near |NES| = 1. It is the reference the proliferation row is
  # read against.
  ggplot2::geom_vline(xintercept = c(-1, 1), linetype = "22",
                      linewidth = 0.25, colour = "grey78") +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey45") +
  ggplot2::geom_point(ggplot2::aes(shape = sig), position = JIT,
                      size = 1.15, stroke = 0.3, colour = "grey25",
                      fill = "grey25") +
  ggplot2::geom_segment(data = med, inherit.aes = FALSE,
                        ggplot2::aes(x = m, xend = m,
                                     y = as.integer(row) - HW,
                                     yend = as.integer(row) + HW),
                        linewidth = 0.5, colour = "black") +
  # The anchor set, labelled because it is the one set all three rulers share.
  # `position` and `nudge_y` are mutually exclusive in ggrepel, so the label is
  # pushed up by constraining the repulsion to the strip above the top row
  # instead -- which is empty, the TEB programmes all being strongly depleted.
  ggrepel::geom_text_repel(
    data = d[d$pathway == ANCHOR, ],
    ggplot2::aes(label = "TEB vs ductal (HS)"), position = JIT,
    size = 1.75, colour = "grey20", seed = 3,
    direction = "y", ylim = c(4.34, 4.56), min.segment.length = 0,
    segment.size = 0.2, segment.colour = "grey55",
    box.padding = 0.12, point.padding = 0.12) +
  ggplot2::scale_shape_manual(values = c(`TRUE` = 21, `FALSE` = 1),
                              breaks = c(TRUE, FALSE),
                              labels = c("padj < 0.05", "n.s."), name = NULL) +
  ggplot2::scale_x_continuous(limits = XR, labels = lab_signed,
                              expand = ggplot2::expansion(mult = 0)) +
  ggplot2::scale_y_discrete(expand = ggplot2::expansion(add = 0.58)) +
  ggplot2::labs(x = "NES, wild-type 6>12W", y = NULL) +
  ggplot2::guides(shape = ggplot2::guide_legend(override.aes = list(size = 1.4))) +
  theme_panel(base_size = 6) +
  # Key inside, TOP RIGHT -- the wedge nothing occupies, because every TEB
  # programme is depleted, so the top row ends well left of zero. (Top LEFT is
  # where the anchor label goes, and bottom right is the only place the
  # proliferation row reaches.)
  ggplot2::theme(
    legend.position        = "inside",
    legend.position.inside = c(0.995, 0.99),
    legend.justification   = c(1, 1),
    legend.background      = ggplot2::element_blank(),
    legend.margin          = ggplot2::margin(0, 0, 0, 0),
    legend.key.size        = ggplot2::unit(2.4, "mm"),
    legend.spacing.y       = ggplot2::unit(0.3, "mm"),
    axis.text.y            = ggplot2::element_text(size = 6),
    plot.margin            = ggplot2::margin(1.5, 2.5, 1, 1.5, "mm"))

# =============================================================================
# the other two rulers -- computed and asserted, reported, NOT drawn
# =============================================================================
teb_gene  <- cmp[cmp$arm == "TEB vs ductal (HS)", ]
prol_gene <- cmp[cmp$arm == "PROLIF_* pooled", ]
teb_null  <- wtn[wtn$arm == "TEB vs ductal (HS)", ]
prol_null <- wtn[wtn$arm == "PROLIF_* pooled", ]
teb_samp  <- tsg[tsg$set == ANCHOR, ]
stopifnot(nrow(teb_gene) == 1L, nrow(prol_gene) == 1L, nrow(teb_samp) == 1L,
          # the two script-43 tables must describe the same arms
          abs(teb_null$observed_wt  - teb_gene$c_wt_time)  < 1e-8,
          abs(prol_null$observed_wt - prol_gene$c_wt_time) < 1e-8,
          # gene-level and purity-adjusted per-sample agree on the TEB arm
          abs(teb_samp$wt_time - teb_gene$c_wt_time) < 0.01)

anchor_row <- d[d$pathway == ANCHOR, ]
fam_line <- function(f) {
  s <- d[d$fam == f, ]
  sprintf("%s: %d sets, median NES %+.2f (range %+.2f to %+.2f), %d significant at padj < 0.05.",
          f, nrow(s), stats::median(s$NES), min(s$NES), max(s$NES), sum(s$sig))
}
# a self-contained multiplicity scope over the sets actually drawn
bh_panel <- stats::p.adjust(d$pval, "BH")
bh_prol  <- min(bh_panel[d$fam == "proliferation"])
bh_teb   <- sum(bh_panel < 0.05 & d$fam == "TEB programmes")

LEGEND <- panel_legend(
  slot = "Fig. 2E",
  what = paste0(
    "Gene-set enrichment across the wild-type 6 to 12 week transition. One point ",
    "per gene set, positioned by its fGSEA normalised enrichment score on the ",
    "`6>12W_wt` contrast, in four families: the terminal-end-bud programmes, the ",
    "ductal sets that are their reciprocal by construction, the TEB sets whose ",
    "membership is intersected with MitoCarta, and the curated proliferation ",
    "sets. Filled points are significant; the vertical rule is the family median."),
  detail = c(
    sprintf("n = %d gene sets on one ranked list (the unshrunken Wald statistic of the wild-type temporal contrast, `timepoint_neg`). %s %s %s %s",
            nrow(d), fam_line("TEB programmes"), fam_line("ductal (TEB-down)"),
            fam_line("TEB x MitoCarta"), fam_line("proliferation")),
    sprintf("THE STRONGEST SINGLE SET IS THE ONE THE OTHER RULERS ALSO USE: %s at NES %+.2f (p %.1e), labelled on the panel.",
            ANCHOR, anchor_row$NES, anchor_row$pval),
    sprintf("THE DIRECTION IS CONTROLLED, WHICH IS WHY THE DUCTAL ROW IS DRAWN. The `_DN` sets hold the genes lower in the terminal end bud than in the mature duct, so they are the reciprocal of the `_UP` sets by construction and a ranking that merely drifted downward could not move them the other way. Two of the three lineages show the reciprocal pattern outright -- basal %+.2f against its UP set's %+.2f, hormone-sensing %+.2f against %+.2f -- and the alveolar pair does not (%+.2f against %+.2f, neither significant). Only the basal DN set clears padj < 0.05 (%.3f).",
            d$NES[d$pathway == "MG_TEB_VS_DUCTAL_BA_GRAY_DN"],
            d$NES[d$pathway == "MG_TEB_VS_DUCTAL_BA_GRAY_UP"],
            d$NES[d$pathway == "MG_TEB_VS_DUCTAL_HS_GRAY_DN"],
            d$NES[d$pathway == ANCHOR],
            d$NES[d$pathway == "MG_TEB_VS_DUCTAL_AP_GRAY_DN"],
            d$NES[d$pathway == "MG_TEB_VS_DUCTAL_AP_GRAY_UP"],
            d$padj_within_category[d$pathway == "MG_TEB_VS_DUCTAL_BA_GRAY_DN"]),
    sprintf("THREE INDEPENDENT RULERS AGREE ON THE TEB ARM, and two of them are independent of cell composition. Gene-level set-average raw log2FC %+.4f over %d genes, the %.1fth percentile of %d expression-matched random sets (empirical p %.4f); and a purity-adjusted per-sample composite (score ~ timepoint * genotype + epithelial + immune) %+.4f, p %.4f. The two agree to %.3f, so residual stromal or immune contamination does not explain the loss.",
            teb_gene$c_wt_time, teb_gene$n_genes, teb_null$percentile,
            ss$defs$n_set_draws, teb_null$p_emp_lower, teb_samp$wt_time,
            teb_samp$wt_p, abs(teb_samp$wt_time - teb_gene$c_wt_time)),
    sprintf("PROLIFERATION IS SMALL, NOT IMMOBILE, and the sentence should carry the qualifier. Not one of the %d proliferation sets is significant on this ruler, and the family median (%+.2f) sits inside the +/-1 band an unmoved set occupies. But on the gene-level ruler the pooled %d proliferation genes move %+.4f log2, which is the %.1fth percentile of their own matched null (p = %.3f). The honest form of \"relatively stable\" is the ratio: the TEB arm moves %.1f times as far.",
            n_fam[["proliferation"]], stats::median(d$NES[d$fam == "proliferation"]),
            prol_gene$n_genes, prol_gene$c_wt_time, prol_null$percentile,
            prol_null$p_emp_lower, abs(teb_gene$c_wt_time / prol_gene$c_wt_time)),
    "THE GENERIC PROLIFERATION REGULONS AGREE WITH THE CURATED SETS -- TFT_E2F1_DOROTHEA_ABC -1.18 and TFT_E2F1_CHIPATLAS +0.86, neither significant -- but the Gray lineage-context lanes of the same transcription factors fall hard (FOXM1_HS_LE -2.10, MYBL2_HS_LE -2.04, TFDP1_HS_LE -2.08, E2F1_HS_LE -1.72, all padj < 1e-3). Those are lineage-identity lanes, not cell-cycle sets. What the normal gland withdraws from is lineage and morphogenesis; the cell cycle itself does not move.",
    sprintf("The %d TEB x MitoCarta lanes are drawn rather than deleted because deleting them would flatter the result and folding them into the TEB row would dilute it (the rule Fig. 1D applies to the same construction class). They are 17 to 59 genes each and their median sits at %+.2f, so they do not track the TEB programme they are named for.",
            n_fam[["TEB x MitoCarta"]],
            stats::median(d$NES[d$fam == "TEB x MitoCarta"])),
    sprintf("Multiplicity, both scopes: drawn significance is script 20's Benjamini-Hochberg WITHIN each library category, the analysis of record. A panel-scope BH over the %d sets drawn here changes nothing -- the same %d TEB programmes survive and the best proliferation set reaches %.3f.",
            nrow(d), bh_teb, bh_prol)),
  bounds = c(
    "BATCH = TIMEPOINT. The 6W and 12W cohorts were extracted as two batches, so every value on this panel is DESCRIBED, not claimed. The mitigation is that the panel makes a COMPARISON, not a measurement: both families ride the same two batches, so a batch artefact would have to be TEB-specific to produce this picture.",
    "The comparison is WITHIN one ranked list, so the fGSEA normaliser is shared. That is what makes it safe here and is exactly what was not true of Fig. 1H, where the flatter twelve-week ranking shrinks the normaliser and lifts every NES. No statement about the SIZE of an NES is made across contrasts anywhere in this panel.",
    "fGSEA's gene permutation is anti-conservative for sets of correlated genes, which is most of this library; and the TEB family spans two categories (03_mammary_development, 06_tf_targets) whose BH scopes are separate by script 20's design.",
    "\"The TEB signature was lost as expected\" -- the expectation is external (Gray et al.'s pubertal-to-adult series); nothing on this panel tests it. What the panel shows is that this cohort reproduces it.",
    "The proliferation family is the 14 curated PROLIF_* sets, which is script 43's own roster and the same genes as the pooled arm quoted above. It is not a claim about every proliferative programme in the library."),
  source = c(
    "results/fgsea_percategory.rds (scripts/20_fgsea_percategory.R) -- NES per ranking x category on the unshrunken Wald statistic; the `timepoint_neg` ranking is the contrast drawn",
    "results/substrate_specificity_tradeoff.rds (scripts/43_substrate_specificity_and_tradeoff.R) -- $comparator and $wt_null for the gene-level arms and their expression-matched nulls; $defs$prolif_sets is the drawn proliferation roster",
    "results/priming_arm_teb.rds (scripts/42_priming_arm_and_teb_substrate.R) -- $teb_signatures, the purity-adjusted per-sample composite"))

save_panel_p(p, "fig2_wt_teb_proliferation", height = 44)

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {

  print(p)
  print(LEGEND)

  ## every drawn set, by family and rank
  d[order(d$fam, d$NES), c("fam", "pathway", "NES", "padj_within_category", "size")] |>
    print(row.names = FALSE, digits = 3)

  ## the TEB programmes by lineage context -- the HS (hormone-sensing) context is
  ## the apoptosis-loaded one (TEB lumen clearance) and it falls furthest
  tb  <- d[d$fam == "TEB programmes", ]
  ctx <- ifelse(grepl("_HS_", tb$pathway), "HS",
         ifelse(grepl("_BA_", tb$pathway), "BA",
         ifelse(grepl("_AP_", tb$pathway), "AP", "other")))
  tapply(tb$NES, ctx, function(v) c(n = length(v), median = stats::median(v))) |> print()

  ## the two rulers the panel does not draw
  cmp[cmp$arm %in% c("TEB vs ductal (HS)", "PROLIF_* pooled"), ] |> print(digits = 4)
  wtn[wtn$arm %in% c("TEB vs ductal (HS)", "PROLIF_* pooled"), ] |> print(digits = 4)
  tsg[1:8, ] |> print(digits = 3)

  ## the generic vs lineage-context regulons of the proliferation TFs
  tn[grepl("^TFT_(E2F1|E2F7|FOXM1|MYBL2|CENPA|TFDP1)_", tn$pathway) &
     !grepl("TEB", tn$pathway),
     c("pathway", "NES", "padj_within_category", "size")] |>
    (\(x) x[order(x$NES), ])() |> print(row.names = FALSE, digits = 3)

  ## the mitochondrial arms, which are Fig. 2F's subject, on this same ruler
  tn[tn$category == "01_mitocarta" & grepl("OXPHOS|COMPLEX", tn$pathway),
     c("pathway", "NES", "padj_within_category", "size")] |>
    (\(x) x[order(x$NES), ])() |> print(row.names = FALSE, digits = 3)
}
