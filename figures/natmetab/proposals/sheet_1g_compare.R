# PROPOSALS (not panels of record) -- Fig. 1g, two forms of the same numbers, side by
# side at 100%: the ratio plane (p1g_bcl2_balance.R) and the three-step rows
# (p1g_alt_ratio_rows.R). Each panel is drawn by its own script (dry run, captured).
suppressPackageStartupMessages(library(grid))
PROP <- here::here("figures", "natmetab", "proposals")
OUT  <- here::here("outputs", "natmetab", "proposals", "Fig1g_compare_sheet.pdf")

grab <- function(script) {
  op <- options(natmetab.dry_run = TRUE, natmetab.capture = TRUE)
  on.exit(options(op))
  if (exists(".nm_captured", envir = globalenv())) rm(".nm_captured", envir = globalenv())
  invisible(capture.output(suppressMessages(source(script, local = new.env(parent = globalenv())))))
  get(".nm_captured", envir = globalenv())
}
plane <- grab(file.path(PROP, "p1g_bcl2_balance.R"))
rows  <- grab(file.path(PROP, "p1g_alt_ratio_rows.R"))

PH <- 297
mm <- function(v) unit(v, "mm")
hy <- function(x) gsub("(?<=[A-Za-z0-9])-(?=[A-Za-z0-9])", "­", x, perl = TRUE)
txt <- function(s, x, y, size = 7, face = "plain")
  grid.text(hy(s), x = mm(x), y = mm(PH - y), just = c("left", "top"),
            gp = gpar(fontsize = size, fontface = face, fontfamily = "Helvetica", lineheight = 1.15))
put <- function(p, x, y, w, h)
  print(p, vp = viewport(x = mm(x), y = mm(PH - y), width = mm(w), height = mm(h),
                         just = c("left", "top")), newpage = FALSE)

pdf(OUT, width = 210 / 25.4, height = PH / 25.4, family = "Helvetica", useDingbats = FALSE)
grid.newpage()
txt("Fig. 1g, two forms of the same nine ratios (proposals, at 100%)", 13, 10, 10, "bold")
txt("A  The plane (56 x 56 mm)", 13, 22);              put(plane, 13, 28, 56, 56)
txt("B  The three steps, one row per ratio (62 x 58 mm)", 95, 22); put(rows, 95, 28, 62, 58)
txt(paste(c(
  "Both read the same three steps, from script 54's object:",
  "  1. The pro-death tilt at 6W: in A, every point lies right of zero; in B, every orange dot.",
  "  2. The common shrink: in A, the line (slope 0.43, fitted to the eight ratios other than Bbc3:Bcl2l1);",
  "     in B, the black tick on each row (0.43 x that row's 6W tilt).",
  "  3. The added effect: in A, Bbc3:Bcl2l1's distance below the line; in B, its amber dot well left of",
  "     its tick (+0.674 at 6W, +0.289 expected, -0.061 observed: the largest departure of the nine).",
  "Grey: the line's 95% prediction interval (A, the band; B, a bar per row). Bbc3:Bcl2l1 lies inside it",
  "in both forms, so neither supports \"significantly\". Bbc3:Mcl1 (-0.26) and Bmf:Bcl2l1 (-0.31) also fall",
  "below their ticks, and Pmaip1:Bcl2l1 rises above its (+0.31)."), collapse = "\n"), 13, 96)
invisible(dev.off())
cat("wrote", OUT, "\n")
