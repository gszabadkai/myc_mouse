# =============================================================================
# fig2_wt_mito_contraction.R -- the normal gland withdraws from the respiratory
# chain, on both rulers, while the rest of the compartment rises
# -----------------------------------------------------------------------------
# SLOT: Fig. 2F.
#
#   "MECs withdraw from the respiratory chain; OXPHOS subunit LFC AND MitoPPS
#    drop across all complexes while biogenesis pathways remain relatively
#    stable ... the maturing gland upregulates amino-acid and lipid catabolism."
#
# THE FORM IS figures/fig04_substrate_specificity.R PANEL C (author's review,
# 2026-08-04): the two rulers plotted AGAINST EACH OTHER rather than side by side.
# The first version was eleven named rows drawn twice; this one puts every
# MitoPathway on the plane and lets the named ones sit in it.
#
#   x  CONTENT   set-average raw log2FC. What the compartment HAS.
#   y  PRIORITY  mitoPPS (Monzel 2025), pairwise-ratio, content-blind. What the
#                compartment SPENDS ITS BUDGET ON -- a uniform scaling cancels.
#
# WHAT THE PLANE SHOWS THAT TWO RANKED LISTS DID NOT:
#   - THE TWO RULERS AGREE (Spearman 0.82 over the 143), so the cloud runs along
#     a diagonal. A content fall alone could be normalisation; a content-blind
#     ratio fall alone could be a reshuffle inside a growing compartment. The
#     lower-left quadrant is the conjunction, and that is where the respiratory
#     chain sits.
#   - THE CLOUD'S CENTRE IS UP AND TO THE RIGHT OF THE ORIGIN. The median
#     MitoPathway GAINS content over this window (+0.041, 72% above zero), so the
#     respiratory arm is not falling with the compartment, it is falling against
#     it. That is what makes "withdraw" the right verb, and on the plane it is
#     visible rather than asserted.
#   - THE ASSEMBLY FACTORS OF THE SAME COMPLEXES SIT AT THE ORIGIN, a few
#     millimetres from their own subunits. That is the internal control, and the
#     distance between the two points is the whole specificity claim.
#
# COLOUR IS THE MANUSCRIPT DIVERGING RAMP ON THE CONTENT AXIS (author's request):
# espresso down, white at zero, mint up. It repeats x rather than adding a
# variable -- the same "encoded twice" idiom as Fig. 1B -- so no colour bar is
# drawn: the x axis IS the key, and a reader can check the asymmetric arm scaling
# against it directly.
#
# SIGNIFICANCE: NONE IS DRAWN, AND THE LEGEND BLOCK SAYS WHAT EXISTS AND WHY NOT.
# Short form: script 08's own mitoPPS test on this contrast gives 0 of 143 at
# BH < 0.05 (23 at raw p < 0.05 against 7.2 expected, binomial p = 8e-7 -- the
# compartment moves coherently, no single pathway carries it), and the tests that
# DO reach significance are at the arm level, where they are decisive: the OXPHOS
# subunits sit at percentile 0.0 of 2000 expression-matched sets while their
# assembly factors sit at 50.2, and the PAIRED null on the difference -- the test
# this sentence's comparison actually needs -- gives OXPHOS subunits minus
# mitoribosome -0.230 at p = 0.0005. A glyph would have to mark 8 tested pathways
# among 135 untested ones, which reads as 135 negatives.
#
# THE mtDNA-ENCODED OXPHOS SUBUNITS ARE EXCLUDED, exactly as Figs. 1E and 1F
# exclude them and as script 40 does (`is_mtdna`). They are the largest single
# movement here (+0.63 content, +0.53 priority) and they are the one quantity in
# this project that cannot be read on a temporal axis -- see the legend block.
#
# BATCH = TIMEPOINT (CLAUDE.md): the 6W and 12W cohorts were extracted as two
# batches, so every value on this panel is DESCRIBED, not claimed.
#
# Reads (read-only, no re-run):
#   results/background_vs_myc.rds              (script 40) -- $ruler, both rulers
#                                                 x four contrasts x 144 pathways
#   results/substrate_specificity_tradeoff.rds (script 43) -- $comparator,
#                                                 $comparator_priority, $wt_null,
#                                                 $paired_null, $defs
#   results/mitopps_scores.rds                 (script 08) -- $mitopps_pairwise,
#                                                 the per-pathway test on this
#                                                 contrast (`Temporal_Myc-`)
# Output: outputs/figures/panels/fig2_wt_mito_contraction.pdf
# =============================================================================

source(here::here("figures", "panels", "_panel_common.R"))

bv_path <- here::here("results", "background_vs_myc.rds")
ss_path <- here::here("results", "substrate_specificity_tradeoff.rds")
mp_path <- here::here("results", "mitopps_scores.rds")
require_fresher_than(bv_path)
require_fresher_than(ss_path)
# NOT require_fresher_than(mp_path): mitopps_scores.rds is timestamped one minute
# BEFORE gsva_scores.rds although both came out of the same post-reconciliation
# re-run, so the guard would stop this panel for an artefact of the ordering
# inside that session -- the same exemption Fig. 1F takes, and for the same file.
# The identity check below is the stronger guard anyway.

bv <- readRDS(bv_path)
ss <- readRDS(ss_path)
mp <- readRDS(mp_path)

# Scripts 40 and 43 save TIBBLES, and a tibble's `[` returns a tibble rather than
# a scalar -- which silently poisons every sprintf downstream (Fig. 1E's trap).
r   <- as.data.frame(bv$ruler)
cmp <- as.data.frame(ss$comparator)
cpr <- as.data.frame(ss$comparator_priority)
wtn <- as.data.frame(ss$wt_null)
pnl <- as.data.frame(ss$paired_null)
arm <- as.data.frame(ss$defs$arms)

# `n_genes` and the content columns carry names in the saved object; strip them or
# a data.frame row lookup returns a named vector and sprintf prints the name too.
r$n <- as.integer(unname(r$n_genes))
r$c <- as.numeric(unname(r$c_tn))       # content, wild-type 6->12W
r$p <- as.numeric(unname(r$p_tn))       # priority, wild-type 6->12W
stopifnot(all(c("pathway", "tier", "is_mtdna") %in% names(r)), nrow(r) == 144L)

# --- guard: a stale results object must fail here, not at review --------------
# Script 43's own positive control, as figures/fig04 uses it: the wild-type
# OXPHOS-subunit value reproduces Issue #4's -0.2548.
ox_obs <- cmp$c_wt_time[cmp$arm == "OXPHOS subunits"]
if (!isTRUE(abs(ox_obs - ss$defs$oxphos_wt_reference) < 0.02))
  stop(sprintf("fig2F: stale results object -- OXPHOS-subunit wild-type %.4f against the reference %.4f. Re-run script 43.",
               ox_obs, ss$defs$oxphos_wt_reference))

# --- the compartment, and three scripts agreeing on it ------------------------
# 143, not 144: the synthetic mtDNA-encoded pathway is dropped exactly as script
# 40 drops it, so this panel, Fig. 1E's lower strip and every regression in the
# corpus describe the same set of pathways.
r143 <- r[!r$is_mtdna, ]
stopifnot(nrow(r143) == 143L)

# (1) scripts 40 and 43 computed the wild-type arms independently and must agree.
# `defs$arms` is script 43's own arm -> ruler-pathway map, so this is an identity
# check and not a re-derivation.
chk <- merge(cmp[, c("arm", "c_wt_time")], arm[, c("arm", "ruler_pathway")], by = "arm")
chk <- chk[!is.na(chk$ruler_pathway), ]
chk$c_ruler <- r$c[match(chk$ruler_pathway, r$pathway)]
cpr$p_ruler <- r$p[match(cpr$ruler_pathway, r$pathway)]
stopifnot(nrow(chk) >= 7L,
          max(abs(chk$c_wt_time - chk$c_ruler)) < 1e-6,
          max(abs(cpr$prio_wt_time - cpr$p_ruler)) < 1e-6)

# (2) the priority ruler IS script 08's, so its per-pathway test can be attached.
# `Temporal_Myc-` is the wild-type 6->12W contrast in script 08's naming.
mpp <- as.data.frame(mp$mitopps_pairwise)
mpp <- mpp[mpp$contrast == "Temporal_Myc-", c("pathway", "diff", "p_value", "padj")]
r143$padj_prio <- mpp$padj[match(r143$pathway, mpp$pathway)]
r143$p_prio    <- mpp$p_value[match(r143$pathway, mpp$pathway)]
stopifnot(!anyNA(r143$padj_prio),
          max(abs(r143$p - mpp$diff[match(r143$pathway, mpp$pathway)])) == 0)

# =============================================================================
# the pathways the sentence names
# =============================================================================
# One per clause, plus the control the OXPHOS clause needs. Named here rather than
# derived, because the sentence names them -- and each is a MitoCarta pathway of
# the ruler, so nothing is re-aggregated in the figure layer.
#
#   respiratory  the five complexes' SUBUNITS (the sentence says "subunit"), and
#                the assembly factors of the same complexes as the control
#   biogenesis   the central-dogma tier and the mitoribosome
#   catabolic    fatty acid oxidation, amino acid, lipid
NAMED <- c("CIV subunits", "CI subunits", "CIII subunits", "CV subunits",
           "CII subunits", "OXPHOS assembly factors",
           "Mitochondrial ribosome", "Mitochondrial central dogma",
           "Lipid metabolism", "Amino acid metabolism", "Fatty acid oxidation")
# Shortened where a MitoCarta name will not fit beside a point. The central-dogma
# shortening is taken from the declared `tier_labels` rather than typed again.
SHORT <- c("OXPHOS assembly factors"     = "OXPHOS assembly",
           "Mitochondrial ribosome"      = "mitoribosome",
           "Mitochondrial central dogma" = "central dogma",   # tier_labels, lower-cased
           "Amino acid metabolism"       = "amino acid",
           "Lipid metabolism"            = "lipid",
           "Fatty acid oxidation"        = "fatty acid ox.")

w <- r143[match(NAMED, r143$pathway), ]
stopifnot(!anyNA(w$c), !anyNA(w$p), nrow(w) == length(NAMED))
w$lab <- sprintf("%s (%d)",
                 ifelse(w$pathway %in% names(SHORT), SHORT[w$pathway], w$pathway), w$n)

# The four proton-circuit complexes fall on BOTH rulers and CII does not. Asserted
# so that a re-run cannot flip the exception without stopping the panel.
FOUR <- c("CI subunits", "CIII subunits", "CIV subunits", "CV subunits")
stopifnot(all(w$c[w$pathway %in% FOUR] < 0), all(w$p[w$pathway %in% FOUR] < 0),
          w$c[w$pathway == "CII subunits"] > 0, w$p[w$pathway == "CII subunits"] > 0)

# =============================================================================
# the panel
# =============================================================================
# WINDOWED FOR DISPLAY, and what falls outside is the membership-loose caveat
# made visible: EVERY pathway at the periphery of this cloud is a set of three to
# six genes (Vitamin D metabolism 4, Selenoproteins 5, GABA metabolism 6, catechol
# 3, molybdenum cofactor 5, the carnitine pair 5 and 6). A set mean's noise scales
# as 1/sqrt(n), so the spread out there is sampling, not biology. Two points sit
# below the window; they are named in the legend block and every number in it is
# computed on the complete 143. Fig. 1G windows its facets the same way.
YLO <- -0.26
XR <- range(r143$c) + c(-1, 1) * diff(range(r143$c)) * 0.05
YR <- c(YLO, max(r143$p) + diff(range(r143$p)) * 0.05)
outside <- r143[r143$p < YLO, ]
drawn   <- r143[r143$p >= YLO, ]
stopifnot(nrow(outside) == 2L, all(outside$n <= 6L),
          !any(outside$pathway %in% NAMED))

# --- where the labels go, and why it is arithmetic rather than nudging ---------
# EVERY label is placed; none is repelled (author's review, 2026-08-04: the leader
# must leave the label ON THE LABEL'S OWN LINE, which ggrepel cannot promise --
# it draws from wherever the box edge happens to be). Seven of the eleven also sit
# inside the knot at the origin or on the crowded upper diagonal, where repel had
# nowhere local to put them and stacked them on each other twice.
#
# THREE LANES, each in a region the assertions below prove is EMPTY of data:
#   upper left   nothing gains priority while losing content     -> the knot
#   bottom band  nothing sits below -0.135 right of -0.10        -> the complexes
#   right band   past +0.285 content, nothing sits near zero     -> the catabolic
#
# The knot's order is the author's: central dogma ABOVE OXPHOS assembly. Ordering
# lanes by the points' own y would put assembly above, and its leader -- which
# ends further LEFT -- then crosses the central-dogma leader, which starts lower
# and ends further right. Rather than encode a rule of thumb, the crossing test
# below checks every pair of leaders directly.
KNOT    <- c("CII subunits", "Mitochondrial central dogma",
             "OXPHOS assembly factors", "Mitochondrial ribosome")
COMPLEX <- c("CV subunits", "CIII subunits", "CI subunits", "CIV subunits")
CATAB   <- c("Fatty acid oxidation", "Amino acid metabolism", "Lipid metabolism")
stopifnot(setequal(c(KNOT, COMPLEX, CATAB), NAMED), setequal(COMPLEX, FOUR))

placed <- rbind(
  data.frame(pathway = KNOT,    x_lab = -0.070, h = 1,
             y_lab = c(0.270, 0.215, 0.160, 0.105)),
  data.frame(pathway = COMPLEX, x_lab = -0.100, h = 0,
             y_lab = c(-0.161, -0.189, -0.217, -0.245)),
  data.frame(pathway = CATAB,   x_lab =  0.285, h = 0,
             y_lab = c(0.100, 0.050, 0.000)))
placed$c   <- w$c[match(placed$pathway, w$pathway)]
placed$p   <- w$p[match(placed$pathway, w$pathway)]
placed$lab <- w$lab[match(placed$pathway, w$pathway)]
# the leader leaves the text on the text's own line, just past its end
placed$x0  <- placed$x_lab + ifelse(placed$h == 1, 0.012, -0.012)

# --- the two things that could go wrong, both checked against the data ---------
in_box <- function(x1, x2, y1, y2)
  sum(r143$c >= x1 & r143$c <= x2 & r143$p >= y1 & r143$p <= y2)

# (a) no lane sits on top of a data point. Boxes cover the text extent, not just
# the anchor, and are re-derived here rather than trusted from a render.
stopifnot(in_box(XR[1],  -0.070, 0.085,  0.290) == 0L,     # knot, right-aligned
          in_box(-0.110,  0.200, -0.258, -0.148) == 0L,    # complexes, left-aligned
          in_box(0.273,   XR[2], -0.020,  0.120) == 0L)    # catabolic, left-aligned

# (b) no two leaders cross. Standard orientation test, all 55 pairs.
ccw <- function(ax, ay, bx, by, cx, cy) (cy - ay) * (bx - ax) > (by - ay) * (cx - ax)
crosses <- function(i, j) {
  a <- c(placed$x0[i], placed$y_lab[i]); b <- c(placed$c[i], placed$p[i])
  d <- c(placed$x0[j], placed$y_lab[j]); e <- c(placed$c[j], placed$p[j])
  ccw(a[1], a[2], d[1], d[2], e[1], e[2]) != ccw(b[1], b[2], d[1], d[2], e[1], e[2]) &&
  ccw(a[1], a[2], b[1], b[2], d[1], d[2]) != ccw(a[1], a[2], b[1], b[2], e[1], e[2])
}
pairs_ij <- utils::combn(nrow(placed), 2)
stopifnot(!any(apply(pairs_ij, 2, function(k) crosses(k[1], k[2]))))

p <- ggplot2::ggplot(drawn, ggplot2::aes(c, p)) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey80") +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey80") +
  ggplot2::geom_point(ggplot2::aes(fill = c), shape = 21, size = 1.1,
                      stroke = 0.15, colour = "grey50") +
  # the named pathways, drawn again a size larger so a leader line has something
  # to land on
  ggplot2::geom_point(data = w, ggplot2::aes(fill = c), shape = 21, size = 2.1,
                      stroke = 0.35, colour = "grey15") +
  # LABELS ARE STEERED, NOT LEFT TO REPEL. Four of the eleven sit inside the dense
  # knot at the origin, where repel has nowhere local to put them and drops the
  # label on its neighbours. The plane has two provably empty regions -- upper
  # left (nothing gains priority while losing content) and the far right -- so the
  # knot's labels are constrained into the upper left and the catabolic three
  # upward, each with a leader. Same principle as Fig. S1C's lanes: place them
  # where the data cannot be, rather than nudging until it looks right.
  ggplot2::geom_segment(data = placed, inherit.aes = FALSE,
                        ggplot2::aes(x = x0, y = y_lab, xend = c, yend = p),
                        linewidth = 0.2, colour = "grey55") +
  ggplot2::geom_text(data = placed, inherit.aes = FALSE,
                     ggplot2::aes(x = x_lab, y = y_lab, label = lab, hjust = h),
                     size = 1.65, colour = "grey15") +
  # White pinned to zero and the two arms scaled independently -- the declared
  # helper. No colour bar: the fill repeats the x axis, so the axis is the key and
  # the asymmetry is inspectable there.
  heat_fill(range(r143$c)) +
  ggplot2::guides(fill = "none") +
  ggplot2::scale_x_continuous(limits = XR, labels = lab_signed,
                              expand = ggplot2::expansion(mult = 0)) +
  ggplot2::scale_y_continuous(limits = YR, labels = lab_signed,
                              expand = ggplot2::expansion(mult = 0)) +
  ggplot2::labs(x = "content, 6>12W_wt  (set-average log2FC)",
                y = "priority, 6>12W_wt  (mitoPPS)") +
  theme_panel(base_size = 6) +
  ggplot2::theme(plot.margin = ggplot2::margin(1.5, 2.5, 1, 1.5, "mm"))

# =============================================================================
# the legend text (never drawn)
# =============================================================================
val <- function(pw, col) w[[col]][w$pathway == pw]
pct <- function(pw, col) 100 * mean(r143[[col]] < val(pw, col))
# "1th" and "3th" are what %.0f plus a hard-coded "th" produces; the suffix has to
# follow the number.
ordinal <- function(x) {
  n <- round(x)
  suf <- if (n %% 100 %in% 11:13) "th" else
    switch(as.character(n %% 10), "1" = "st", "2" = "nd", "3" = "rd", "th")
  paste0(n, suf)
}
line <- function(pw, label = pw)
  sprintf("%s %+.3f content (%s percentile of the compartment) and %+.3f priority (%s)",
          label, val(pw, "c"), ordinal(pct(pw, "c")),
          val(pw, "p"), ordinal(pct(pw, "p")))
null_of <- function(a) wtn$percentile[wtn$arm == a]
pair_of <- function(b) pnl[pnl$arm_b == b, ]
n_raw   <- sum(r143$p_prio < 0.05)
binom_p <- stats::binom.test(n_raw, nrow(r143), 0.05, alternative = "greater")$p.value
mt <- r[r$is_mtdna, ]

LEGEND <- panel_legend(
  slot = "Fig. 2F",
  what = paste0(
    "Every mitochondrial pathway over the wild-type 6 to 12 week window, on two ",
    "rulers at once. Horizontal, content: the set-average raw log2 fold change, ",
    "what the compartment has. Vertical, priority: the mitoPPS pairwise-ratio ",
    "score, which is blind to content and reads how the compartment divides its ",
    "budget. Fill repeats the horizontal axis on the manuscript diverging scale. ",
    "The pathways the text names are drawn larger and labelled, with their gene ",
    "counts."),
  detail = c(
    sprintf("n = 6 wild-type animals per timepoint; %d MitoPathways (the synthetic mtDNA-encoded pathway excluded, as in script 40 and Figs. 1E and 1F). Percentiles quoted below are positions within that compartment.",
            nrow(r143)),
    sprintf("THE TWO RULERS AGREE, which is why the cloud runs along a diagonal and why both are drawn: across all %d pathways they correlate at Spearman %.2f (Pearson %.2f). A drop in content alone could be a normalisation effect; a drop in a content-blind ratio alone could be a reshuffle inside a growing compartment. The lower-left quadrant is the conjunction.",
            nrow(r143), stats::cor(r143$c, r143$p, method = "spearman"),
            stats::cor(r143$c, r143$p)),
    sprintf("THE CLOUD SITS UP AND RIGHT OF THE ORIGIN, which is what makes \"withdraw\" the right verb: the median MitoPathway GAINS content over this window (%+.3f, %.0f%% of the %d above zero), while on the content-blind ruler the compartment is by construction near zero (median %+.3f). The respiratory arm is not falling with the compartment, it is falling against it.",
            stats::median(r143$c), 100 * mean(r143$c > 0), nrow(r143),
            stats::median(r143$p)),
    sprintf("THE RESPIRATORY ARM, in the lower-left corner. %s; %s; %s; %s.",
            line("CIV subunits"), line("CI subunits"), line("CIII subunits"),
            line("CV subunits")),
    sprintf("THE INTERNAL CONTROL SITS AT THE ORIGIN, a few millimetres from its own subunits: %s. These are the assembly factors OF THE SAME COMPLEXES, measured in the same libraries on the same two batches. What the gland withdraws is the structural stoichiometry of the chain, not the machinery that builds it.",
            line("OXPHOS assembly factors")),
    sprintf("BIOGENESIS IS STABLE: %s; %s. This is the clause that separates the result from a general shrinkage of the mitochondrial programme.",
            line("Mitochondrial central dogma", "central dogma"),
            line("Mitochondrial ribosome", "mitoribosome")),
    sprintf("CATABOLISM RISES: %s; %s; %s.", line("Fatty acid oxidation"),
            line("Amino acid metabolism"), line("Lipid metabolism")),
    sprintf("SIGNIFICANCE, AND WHY NONE IS DRAWN. Script 08's own test on this contrast (`Temporal_Myc-`) gives %d of %d pathways at BH-adjusted p < 0.05 -- none -- with the smallest adjusted p at %.2f. What it does show is coherence: %d pathways clear an unadjusted p of 0.05 against %.1f expected (binomial p = %.0e), so the compartment moves together and no single pathway carries it. The tests that DO reach significance are at the arm level, and they are in the next two items. A glyph on this panel would have to mark 8 tested pathways among %d untested ones, which reads as %d negatives.",
            sum(r143$padj_prio < 0.05), nrow(r143), min(r143$padj_prio), n_raw,
            0.05 * nrow(r143), binom_p, nrow(r143) - 8L, nrow(r143) - 8L),
    sprintf("THE ARM-LEVEL TEST IS AN EXPRESSION-MATCHED NULL (script 43, %d random gene sets drawn within baseMean ventiles). Pooled over all 87 nuclear-encoded OXPHOS subunits the content effect is %+.4f at percentile %.1f (p < %.4f); their assembly factors sit at percentile %.1f, the middle of the distribution; the mitoribosome at %.1f; and amino-acid and lipid metabolism beat the null in the UPWARD direction, at percentiles %.1f and %.1f.",
            ss$defs$n_set_draws, cmp$c_wt_time[cmp$arm == "OXPHOS subunits"],
            null_of("OXPHOS subunits"), 1 / ss$defs$n_set_draws,
            null_of("OXPHOS assembly"), null_of("mitoribosome"),
            null_of("amino-acid metabolism"), null_of("lipid metabolism")),
    sprintf("AND THE COMPARISON THE SENTENCE MAKES HAS ITS OWN NULL, which is the number to quote: script 43 redraws BOTH sets of a contrast together and tests the DIFFERENCE. OXPHOS subunits minus mitoribosome = %+.3f, percentile %.2f, p = %.4f; minus nucleotide metabolism %+.3f, p = %.4f; minus the pooled proliferation set %+.3f, p < %.4f. \"The respiratory chain falls while biogenesis does not\" is a comparative claim, and this is the test of it.",
            pair_of("mitoribosome")$observed_diff, pair_of("mitoribosome")$percentile,
            pair_of("mitoribosome")$p_emp_lower,
            pair_of("nucleotide metabolism")$observed_diff,
            pair_of("nucleotide metabolism")$p_emp_lower,
            pair_of("PROLIF_* pooled")$observed_diff, 1 / ss$defs$n_set_draws),
    sprintf("\"ACROSS ALL COMPLEXES\" HAS ONE EXCEPTION AND IT IS INFORMATIVE: %s. Complex II is the only respiratory complex with no mtDNA-encoded subunit -- it is not in the proton circuit and it is also a TCA enzyme -- and it is the only one that does not fall. But the set is FOUR genes, and MitoCarta sets are membership-loose, so this is a direction to note and not a mechanism to claim.",
            line("CII subunits")),
    sprintf("EXCLUDED, AND IT IS THE LARGEST MOVEMENT IN THE COMPARTMENT: the 13 mtDNA-encoded OXPHOS subunits rise %+.3f on content and %+.3f on priority. They are dropped by the same rule Figs. 1E and 1F use. The reason is not tidiness: the mtDNA-encoded read fraction is confounded three ways in these data -- real content, the proliferation denominator, and dissociation leak -- and it is the one quantity that is time-associated rather than genotype-associated, so a temporal contrast is exactly where it cannot be read. The nuclear-encoded arm falling while the mtDNA arm rises is a mitonuclear discordance if it is real; on this axis it is not adjudicable.",
            as.numeric(unname(mt$c_tn)), as.numeric(unname(mt$p_tn)))),
  bounds = c(
    "BATCH = TIMEPOINT. The 6W and 12W cohorts were extracted as two separate batches, so every value on this panel is DESCRIBED, not claimed. Two things mitigate it and neither removes it: the withdrawal is SPECIFIC within the compartment (the assembly factors of the same complexes, in the same libraries on the same batches, sit at the origin), and it appears on a content-blind ratio ruler as well as on the content one.",
    "NO PATHWAY ON THIS PANEL IS INDIVIDUALLY SIGNIFICANT after correction, on either ruler. The claim is a pattern -- a quadrant, and a distance between two points that share their complexes -- supported by arm-level nulls, and it should be written that way.",
    "There is no per-pathway test on the CONTENT ruler for this contrast in the saved objects. One could be computed as a one-sample t-test over member-gene fold changes, but that null ignores inter-gene correlation and would be anti-conservative for exactly the coherently-regulated sets this panel is about; it is deliberately not done.",
    "mitoPPS is RELATIVE BY CONSTRUCTION: a pathway can be demoted while its absolute expression rises. On this panel both rulers point the same way for the respiratory arm, so no such reading is needed -- but a priority value must never be reported as a fall in expression.",
    "MitoCarta sets are membership-loose and several labelled pathways are small (CII subunits 4 genes, CIII subunits 9). A large effect on a four-gene set is one gene, not a module, which is why the gene count is on the face of the panel.",
    "The fill repeats the horizontal axis; it is not a second variable. White is pinned to zero and the two arms of the ramp are scaled independently, so equal ink does not mean equal magnitude across the sign change -- the x axis carries the magnitude and is the key.",
    "Within-compartment percentiles are a RANK AMONG MITOPATHWAYS, not a significance statement.",
    sprintf("WINDOWED FOR DISPLAY, as Fig. 1G windows its facets: %d pathways sit below the drawn priority range and are not shown (%s). Every number in this legend is computed on the complete %d. What falls outside is also the membership-loose caveat made visible -- both are sets of five genes or fewer, and every pathway at the periphery of this cloud is a set of three to six genes, where a set mean's noise scales as one over the square root of n.",
            nrow(outside), paste(sprintf("%s %+.3f, %d genes", outside$pathway,
                                         outside$p, outside$n), collapse = "; "),
            nrow(r143))),
  source = c(
    "results/background_vs_myc.rds (scripts/40_background_vs_myc_decomposition.R) -- $ruler, content (c_tn) and priority (p_tn) on the wild-type temporal contrast for 144 MitoPathways",
    "results/substrate_specificity_tradeoff.rds (scripts/43_substrate_specificity_and_tradeoff.R) -- $comparator, $comparator_priority, $wt_null and $paired_null for the arm-level values and their expression-matched nulls; $defs$arms is the arm-to-pathway map the identity check uses",
    "results/mitopps_scores.rds (scripts/08_mitoPPS_analysis.R) -- $mitopps_pairwise, contrast `Temporal_Myc-`, the per-pathway test on the priority ruler; its `diff` is asserted identical to the ruler's p_tn",
    "mitoPPS: Monzel et al. 2025, external/mitotyping/; pairwise-ratio scores on linear-scale DESeq2 normalised counts"))

save_panel_p(p, "fig2_wt_mito_contraction", height = 70)

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {

  print(p)
  print(LEGEND)

  ## the labelled pathways, both rulers, their compartment percentiles and the
  ## per-pathway mitoPPS test -- this is the first version of the panel, which
  ## drew these eleven as two ranked columns
  data.frame(pathway = w$pathway, n = w$n, content = w$c, priority = w$p,
             c_pct = vapply(w$pathway, pct, numeric(1), "c"),
             p_pct = vapply(w$pathway, pct, numeric(1), "p"),
             prio_p = w$p_prio, prio_padj = w$padj_prio) |>
    (\(x) x[order(x$content), ])() |> print(row.names = FALSE, digits = 3)

  ## the whole OXPHOS tier, which is where the subunit / assembly split shows
  r[r$tier == "OXPHOS", c("pathway", "n", "c", "p")] |>
    (\(x) x[order(x$c), ])() |> print(row.names = FALSE, digits = 3)

  ## the central-dogma tier -- nothing in it moves
  r[r$tier == "Mitochondrial central dogma", c("pathway", "n", "c", "p")] |>
    (\(x) x[order(x$c), ])() |> print(row.names = FALSE, digits = 3)

  ## the corners of the plane, labelled or not
  r143[order(r143$c), c("pathway", "tier", "n", "c", "p")] |> head(10) |>
    print(row.names = FALSE, digits = 3)
  r143[order(-r143$c), c("pathway", "tier", "n", "c", "p")] |> head(10) |>
    print(row.names = FALSE, digits = 3)

  ## the 23 pathways at unadjusted p < 0.05 on the priority ruler -- none survives
  r143[r143$p_prio < 0.05,
       c("pathway", "tier", "n", "c", "p", "p_prio", "padj_prio")] |>
    (\(x) x[order(x$p), ])() |> print(row.names = FALSE, digits = 3)

  ## the arm-level view with the expression-matched nulls, and the paired null,
  ## which figures/fig04 panels A and B draw
  merge(cmp[, c("arm", "n_genes", "c_wt_time")],
        wtn[, c("arm", "percentile", "p_emp_lower")], by = "arm") |>
    (\(x) x[order(x$c_wt_time), ])() |> print(row.names = FALSE, digits = 3)
  pnl |> print(digits = 3)
}
