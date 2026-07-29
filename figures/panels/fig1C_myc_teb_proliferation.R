# =============================================================================
# fig1C_myc_teb_proliferation.R -- endogenous Myc, the transgene, and the
# programmes each one moves
# -----------------------------------------------------------------------------
# SUPPORTS (Results, "Myc drives early breast tumourigenesis by inducing OXPHOS
# and biosynthetic pathways", paragraph 1):
#   "an endogenous Myc program contributed to the pubertal TEB state in WT
#    animals, promoting proliferation (Fig. 1C), evidenced by the reduction of
#    the canonical Myc, TEB and proliferation signatures in the adult stage
#    (WT 6W -> 12W). The addition of the Myc transgene amplified the TEB and
#    proliferation effects and suppressed the BMYO lineage in favor of
#    differentiation to LHS. However, the Myc+ 6W -> 12W TEB-proliferation
#    trajectory also showed a reduction similar to the WT trend"
#
# The four columns are the four contrasts of the design (Fig. S1B), which is the
# decomposition script 27 was built to make: ENDOGENOUS = the wild-type 6->12W
# shift, TRANSGENE = the genotype main effect. Everything is divided by each
# programme's own pooled within-group SD so the two halves are comparable -- and
# the comparison is the point, because the transgene effect is roughly three
# times the endogenous one.
#
# FOUR PLACES THE SENTENCE AND THE OBJECT DISAGREE, all visible in the panel:
#  1. "the transgene amplified the TEB ... effects" -- TEB-ductal has NO genotype
#     main effect (+0.63 SD, p = 0.14). Proliferation does (+1.01 SD, p = 0.019),
#     and so does MYC-in-TEB (+1.61 SD, p = 6.2e-04) -- but that set is MYC
#     TARGETS scored in the TEB context, not the TEB phenotype.
#  2. "in favor of differentiation to LHS" -- LHS has no genotype main effect
#     (+0.08 SD, p = 0.86). What it has is a time-dependent flip: below wild type
#     at 6W, above at 12W (interaction p = 0.097). So the claim holds at 12W, as
#     a trend, not as a pooled effect.
#  3. "suppressed the BMYO lineage" -- this one HOLDS and is the powered state
#     result (-0.98 SD, p = 0.023).
#  4. The endogenous half is DIRECTIONAL, not significant: every wild-type 6->12W
#     shift has p 0.21-0.65 at n=6. The six programmes all move the same way,
#     which is the evidence; no single one is a test.
#
# Input:  results/myc_endogenous_amplification.rds  (script 27 -- $prog_stats)
#         results/dev_program_myc_integration.rds   (script 26 -- $state_stats)
# Output: outputs/figures/panels/fig1C_myc_teb_proliferation.pdf
# =============================================================================

source(here::here("figures", "panels", "_panel_common.R"))

prog_path <- here::here("results", "myc_endogenous_amplification.rds")
dev_path  <- here::here("results", "dev_program_myc_integration.rds")
require_fresher_than(prog_path)
require_fresher_than(dev_path)

prog <- readRDS(prog_path)
dev  <- readRDS(dev_path)
stopifnot("prog_stats" %in% names(prog), "state_stats" %in% names(dev))

# --- one tidy table, both sources, everything in within-group SD units --------
# script 27 stores geno_d already; script 26 stores geno_beta on the raw GSVA
# scale. Dividing by within_sd reproduces geno_d exactly (asserted below), so the
# two objects can share an axis.
ps <- prog$prog_stats
ss <- dev$state_stats
stopifnot(max(abs(ps$geno_d - ps$geno_beta / ps$within_sd)) < 1e-8)

tidy_one <- function(df, label, block) {
  data.frame(
    programme = label,
    block     = block,
    `WT 6W->12W`   = df$wt_shift     / df$within_sd,
    `Myc+ 6W->12W` = df$mycpos_shift / df$within_sd,
    `Myc+ vs WT`   = df$geno_beta    / df$within_sd,
    `interaction`  = df$int_beta     / df$within_sd,
    p_1 = df$wt_p, p_2 = df$mycpos_p, p_3 = df$geno_p, p_4 = df$int_p,
    check.names = FALSE, stringsAsFactors = FALSE)
}

raw <- rbind(tidy_one(ps, ps$label,               "signatures"),
             tidy_one(ss, as.character(ss$state), "MEC states"))

contr_levels <- c("WT 6W->12W", "Myc+ 6W->12W", "Myc+ vs WT", "interaction")

eff <- raw |>
  tidyr::pivot_longer(dplyr::all_of(contr_levels),
                      names_to = "contrast", values_to = "effect")
pv  <- raw |>
  dplyr::select(programme, block, p_1, p_2, p_3, p_4) |>
  tidyr::pivot_longer(dplyr::starts_with("p_"), names_to = "slot", values_to = "p") |>
  dplyr::mutate(contrast = contr_levels[as.integer(sub("^p_", "", slot))]) |>
  dplyr::select(-slot)

dat <- dplyr::inner_join(eff, pv, by = c("programme", "block", "contrast")) |>
  dplyr::mutate(
    contrast = factor(contrast, levels = contr_levels),
    block    = factor(block, levels = c("signatures", "MEC states")),
    sig      = p < 0.05)
stopifnot(nrow(dat) == 36)

# top-to-bottom reading order; ggplot draws discrete y bottom-up, hence rev()
prog_order <- c("MYC signatures (17)", "Felsher", "Hallmark MYC V2",
                "MYC-in-TEB (Gray)", "Proliferation (14)", "TEB - ductal",
                "BMYO", "LASP", "LHS")
stopifnot(setequal(prog_order, unique(dat$programme)))
dat$programme <- factor(dat$programme, levels = rev(prog_order))

# --- panel -------------------------------------------------------------------
sig_fill <- c(`TRUE` = "#D73027", `FALSE` = "grey85")

p <- ggplot2::ggplot(dat, ggplot2::aes(x = effect, y = programme)) +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey80") +
  ggplot2::geom_segment(ggplot2::aes(x = 0, xend = effect,
                                     y = programme, yend = programme),
                        linewidth = 0.32, colour = "grey55") +
  ggplot2::geom_point(ggplot2::aes(fill = sig), shape = 21, size = 1.7,
                      colour = "grey25", stroke = 0.25) +
  ggplot2::facet_grid(block ~ contrast, scales = "free_y", space = "free_y") +
  ggplot2::scale_fill_manual(values = sig_fill, name = NULL,
                             breaks = "TRUE", labels = "p < 0.05") +
  ggplot2::labs(x = "effect (within-group SD)", y = NULL) +
  theme_panel() +
  ggplot2::theme(legend.position = "bottom",
                 legend.margin = ggplot2::margin(-3, 0, 0, 0),
                 axis.line.y   = ggplot2::element_blank(),
                 axis.ticks.y  = ggplot2::element_blank(),
                 strip.text.y  = ggplot2::element_text(angle = -90),
                 panel.spacing.x = ggplot2::unit(2.0, "mm"),
                 panel.spacing.y = ggplot2::unit(1.6, "mm"))

# --- the legend text (never drawn) -------------------------------------------
g   <- function(nm, col) ps[[col]][ps$program == nm]
f   <- function(x) sprintf("%+.2f", x)
st  <- function(nm, col) ss[[col]][as.character(ss$state) == nm]

LEGEND <- panel_legend(
  slot = "Fig. 1C",
  what = paste0(
    "Effect of each of the four contrasts on the canonical MYC, proliferation ",
    "and TEB signatures and on the three mammary epithelial states, in units of ",
    "each programme's own within-group standard deviation."),
  detail = c(
    "n = 6 per group, n = 24. Points filled where p < 0.05; every test is an ordinary least-squares contrast on the per-sample GSVA composite.",
    "Columns are the contrasts of Fig. S1B: the wild-type and Myc+ timelines, the genotype main effect pooled over age, and the interaction (how much of the genotype gap changes between 6 and 12 weeks).",
    sprintf("ENDOGENOUS (wild type, 6 to 12 weeks) - all six signatures fall together and none reaches significance: MYC signatures %s (p = %.2f), Felsher %s (p = %.2f), Hallmark MYC V2 %s (p = %.2f), MYC-in-TEB %s (p = %.2f), proliferation %s (p = %.2f), TEB-ductal %s (p = %.2f). The evidence is the shared direction, not any one test.",
            f(g("myc","wt_d")), g("myc","wt_p"), f(g("felsher","wt_d")), g("felsher","wt_p"),
            f(g("hallmark_v2","wt_d")), g("hallmark_v2","wt_p"),
            f(g("myc_in_teb","wt_d")), g("myc_in_teb","wt_p"),
            f(g("prolif","wt_d")), g("prolif","wt_p"),
            f(g("teb_ductal","wt_d")), g("teb_ductal","wt_p")),
    sprintf("TRANSGENE (genotype main effect) - powered and large on the MYC axis: MYC signatures %s (p = %.1e), Felsher %s, Hallmark MYC V2 %s; proliferation %s (p = %.3f); MYC-in-TEB %s (p = %.1e).",
            f(g("myc","geno_d")), g("myc","geno_p"), f(g("felsher","geno_d")),
            f(g("hallmark_v2","geno_d")), f(g("prolif","geno_d")), g("prolif","geno_p"),
            f(g("myc_in_teb","geno_d")), g("myc_in_teb","geno_p")),
    sprintf("The transgene effect is about three times the endogenous one on the same programmes (for example MYC signatures %s versus %s), which is what 'amplifies the same axis' means quantitatively.",
            f(g("myc","geno_d")), f(g("myc","wt_d"))),
    sprintf("Myc+ timeline - the Myc+ gland travels the same way as wild type, and for TEB-ductal significantly so: %s (p = %.4f).",
            f(g("teb_ductal","mycpos_shift") / g("teb_ductal","within_sd")),
            g("teb_ductal","mycpos_p")),
    sprintf("MEC states - BMYO is suppressed by the transgene, %s (p = %.3f), the only powered state effect. LASP %s (p = %.2f) and LHS %s (p = %.2f) have none; LHS instead flips with age (interaction p = %.2f), sitting below wild type at 6 weeks and above it at 12.",
            f(st("BMYO","geno_beta") / st("BMYO","within_sd")), st("BMYO","geno_p"),
            f(st("LASP","geno_beta") / st("LASP","within_sd")), st("LASP","geno_p"),
            f(st("LHS","geno_beta")  / st("LHS","within_sd")),  st("LHS","geno_p"),
            st("LHS","int_p"))),
  bounds = c(
    "The two halves are not equally powered and should not be reported as if they were. The genotype main effect is a 12-versus-12 comparison balanced within batch and is clean. The wild-type timeline is 6 versus 6, confounded with batch, and every one of its shifts is non-significant - it is a consistent direction across six programmes, described and not claimed.",
    sprintf("The TEB genotype effect does NOT reach significance (%s, p = %.2f). MYC-in-TEB does (p = %.1e) but is MYC targets scored in the TEB context, so it reports the transgene, not the TEB phenotype. A sentence saying the transgene amplified the TEB effect should rest on proliferation and on MYC-in-TEB, or be qualified.",
            f(g("teb_ductal","geno_d")), g("teb_ductal","geno_p"),
            g("myc_in_teb","geno_p")),
    "The signatures are heavily inter-correlated (MYC signatures with Felsher r = 0.99, with Hallmark V2 r = 0.98, with proliferation r = 0.77), so the six rows are not six independent observations and the shared direction is partly shared genes.",
    "Composite p-values are indicative: the sets inside each composite are correlated, so the composite is not independent replication.",
    "Endogenous Myc establishing the pubertal TEB phenotype is a literature claim; what these data show is wild-type co-variation consistent with it, plus that the transgene moves the same axis. Not causal."),
  source = c(
    "results/myc_endogenous_amplification.rds (scripts/27) -- $prog_stats; composites: MYC signatures 17 sets, proliferation 14, MYC-in-TEB 4, TEB-ductal = mean(UP) - mean(DN) over 3 + 3",
    "results/dev_program_myc_integration.rds (scripts/26, PART C2) -- $state_stats",
    "Both derived from results/gsva_scores.rds (scripts/15)"))

save_panel_p(p, "fig1C_myc_teb_proliferation",
             width = fig_w[["onehalf"]], height = 62)

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {

  print(p)
  print(LEGEND)

  ## the whole table behind the panel
  dat |> dplyr::arrange(contrast, programme) |> print(n = 40)

  ## endogenous vs transgene, side by side
  ps |> dplyr::select(label, wt_d, wt_p, geno_d, geno_p) |> print()
}
