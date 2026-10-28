# =============================================================================
# fig3i_prism_tables.R -- the Fig. 3i and Extended Data tables, for Prism
# -----------------------------------------------------------------------------
# The author is drawing these panels in Prism rather than R, to keep Figure 3
# uniform (2026-10-28). This writes what Prism pastes. It draws nothing.
#
# TWO KINDS OF TABLE, AND THEY ARE NOT INTERCHANGEABLE:
#
#   PER-ANIMAL (one column per arm, seven rows). Paste into a Prism Column or
#     Grouped table and let Prism compute whatever it computes. These are raw
#     values, so anything Prism derives from them is Prism's own statistic.
#
#   CONTRAST (one row per arm: value, and the asymmetric error offsets Prism
#     wants). These are script 52's difference of medians against each arm's own
#     matched control, with its 5,000-resample bootstrap interval. PRISM CANNOT
#     RECOMPUTE THEM from per-animal values -- a bootstrap of a difference of
#     medians between two arms is not a statistic Prism offers, and the two arms
#     are different tables. Enter them as summary values: in an XY or Grouped
#     table set "Enter and plot error values already calculated elsewhere", and
#     use the err_minus / err_plus columns, which are the DISTANCES from the
#     value to each end of the interval, not the ends themselves.
#
#   If a contrast panel is instead drawn from per-animal values and Prism's own
#   test, the legend must cite that test and not script 52's interval. That is
#   the same split Fig. 3g already has.
#
# ORDER is the escalation (author, 2026-10-28): no transgene, checkpoint
#   disrupted, pathway intact, effector blocked. Columns are written in that
#   order so a paste lands correctly without reordering.
#
# UNITS: CPM for Ppargc1a (the unit of record); mitoPPS and ox_rel in their own
#   units; MitoCarta share as a PERCENTAGE of the library.
#
# Reads : results/orthotopic_section3_panels.rds  (script 57, which reads script 52)
# Writes: outputs/natmetab/Fig3/Fig3i_prism_retained_dose.csv
#         outputs/natmetab/Fig3/Fig3i_prism_mitopps_contrast.csv
#         outputs/natmetab/Fig3/Fig3i_prism_mitopps_per_animal.csv
#         outputs/natmetab/EDFig5/EDFig5b_prism_ox_rel_contrast.csv
#         outputs/natmetab/EDFig5/EDFig5b_prism_ox_rel_per_animal.csv
#         outputs/natmetab/EDFig5/EDFig5_prism_mito_share_per_animal.csv
# =============================================================================

source(here::here("figures", "natmetab", "_style.R"))

o  <- readRDS(here::here("results", "orthotopic_section3_panels.rds"))
pa <- as.data.frame(o$per_animal)
sn <- as.data.frame(o$summary_numbers)

ROUTE <- c("BclxL", "KOPgc1a", "NTPgc1a", "Pgc1a_BclxL")          # the escalation
CTRL  <- c(BclxL = "EV", KOPgc1a = "KOEV", NTPgc1a = "NTEV", Pgc1a_BclxL = "BclxL")
stopifnot(setequal(ROUTE, sn$arm))

d3 <- here::here("outputs", "natmetab", "Fig3")
d5 <- here::here("outputs", "natmetab", "EDFig5")
for (d in c(d3, d5)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# --- per-animal, one column per arm -------------------------------------------
# Column names carry the route AND the arm, so a pasted block is self-identifying.
wide <- function(arms, col, digits = NULL) {
  m <- vapply(arms, function(a) {
    v <- pa[[col]][pa$arm == a]
    if (!is.null(digits)) v <- round(v, digits)
    v }, numeric(7))
  colnames(m) <- sprintf("%s (%s)", background_labels[arms], arms)
  as.data.frame(m)
}

# --- a contrast row, with the offsets Prism wants ------------------------------
contrast <- function(val, lo, hi) {
  data.frame(
    order        = seq_along(ROUTE),
    route        = unname(background_labels[ROUTE]),
    arm          = ROUTE,
    control      = unname(CTRL[ROUTE]),
    value        = round(val, 4),
    ci_low       = round(lo, 4),
    ci_high      = round(hi, 4),
    err_minus    = round(val - lo, 4),      # distance DOWN from value, for Prism
    err_plus     = round(hi - val, 4),      # distance UP   from value, for Prism
    excludes_zero = ifelse(hi < 0 | lo > 0, "yes", "no"),
    stringsAsFactors = FALSE) }

i <- match(ROUTE, sn$arm)
mp <- contrast(sn$mitopps_diff[i], sn$mitopps_lo[i], sn$mitopps_hi[i])
ox <- contrast(sn$ox_rel_diff[i],  sn$ox_rel_lo[i],  sn$ox_rel_hi[i])

# --- MitoCarta share: an ABSOLUTE per-animal quantity, so every arm is given ---
# This is the mass readout -- the one quantity here that reads as "overall
# mitochondrial level". Each PGC-1a arm is written beside its own control, in the
# escalation order, so the pairing survives the paste. Percentage of the library.
SHARE_ARMS <- as.vector(rbind(unname(CTRL[ROUTE]), ROUTE))      # control, arm, control, arm...
share <- vapply(SHARE_ARMS, function(a) round(100 * pa$mito_share[pa$arm == a], 3), numeric(7))
colnames(share) <- ifelse(SHARE_ARMS %in% ROUTE,
                          sprintf("%s (%s)", background_labels[SHARE_ARMS], SHARE_ARMS),
                          sprintf("control (%s)", SHARE_ARMS))
share <- as.data.frame(share)

w <- function(x, dir, file) {
  utils::write.csv(x, file.path(dir, file), row.names = FALSE)
  message("wrote ", file.path(dir, file)) }

w(wide(ROUTE, "Ppargc1a_cpm", 4),   d3, "Fig3i_prism_retained_dose.csv")
w(mp,                               d3, "Fig3i_prism_mitopps_contrast.csv")
w(wide(ROUTE, "mitopps_oxphos", 5), d3, "Fig3i_prism_mitopps_per_animal.csv")
w(ox,                               d5, "EDFig5b_prism_ox_rel_contrast.csv")
w(wide(ROUTE, "ox_rel", 5),         d5, "EDFig5b_prism_ox_rel_per_animal.csv")
w(share,                            d5, "EDFig5_prism_mito_share_per_animal.csv")

# --- printed, so the tables can be read without opening a file ----------------
op <- options(width = 130); on.exit(options(op))
show <- function(title, x, note = NULL) {
  cat("\n", strrep("=", 78), "\n", title, "\n", sep = "")
  if (!is.null(note)) cat(note, "\n")
  print(x, row.names = FALSE)
}
show("1. Fig. 3i LEFT -- Ppargc1a retained, CPM, one column per arm, seven animals",
     wide(ROUTE, "Ppargc1a_cpm", 4),
     "   per-animal values; Prism computes whatever it computes from these")
show("2. Fig. 3i RIGHT -- OXPHOS priority (mitoPPS) vs each arm's matched control",
     mp, "   SUMMARY values: enter already-calculated error, using err_minus / err_plus")
show("3. Extended Data b -- ox_rel vs each arm's matched control",
     ox, "   SUMMARY values, as above. Agrees with mitoPPS on every arm.")
show("4. MitoCarta share of the library, %, per animal -- each arm beside its control",
     share, "   the mass readout; absolute, not a difference")
cat("\n   MitoCarta share, median % per arm:\n")
for (a in SHARE_ARMS)
  cat(sprintf("     %-28s %.2f%%\n",
              if (a %in% ROUTE) sprintf("%s (%s)", background_labels[[a]], a) else paste("control", a),
              100 * stats::median(pa$mito_share[pa$arm == a])))
cat("\n   per-animal mitoPPS and ox_rel are written as files too, if a panel is to be\n",
    "  drawn from animals rather than from the saved contrast.\n")
