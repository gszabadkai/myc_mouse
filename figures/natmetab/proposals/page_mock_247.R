# PROPOSALS (not panels of record) -- Fig. 1, ED Fig. 2 and ED Fig. 3 on the
# 183 x 247 mm page the author set for the initial submission (2026-10-14; ED 180
# wide). Every panel is drawn by its own script (dry run, captured through
# natmetab.capture) at the size it is placed; only Fig. 1's bench row is a
# placeholder. Writes outputs/natmetab/proposals/Fig1_ED2_ED3_247.pdf (A4 sheets).
#
# Second round (author's notes, 2026-10-14): 1f carries ED 2f's statistics
# (p1f_window_combined.R), so ED 2f leaves and ED 2g becomes 2f; ED 2c is lower;
# ED 2d and 2e share one data-area size; ED 3a is the per-animal form of
# Fig. 2H+I (alt) (ped3a_per_animal.R).
suppressPackageStartupMessages(library(grid))
NAT  <- here::here("figures", "natmetab")
PROP <- file.path(NAT, "proposals")
OUT  <- here::here("outputs", "natmetab", "proposals", "Fig1_ED2_ED3_247.pdf")

grab <- function(script) {
  op <- options(natmetab.dry_run = TRUE, natmetab.capture = TRUE)
  on.exit(options(op))
  if (exists(".nm_captured", envir = globalenv())) rm(".nm_captured", envir = globalenv())
  invisible(capture.output(suppressMessages(suppressWarnings(
    source(script, local = new.env(parent = globalenv()))))))
  get(".nm_captured", envir = globalenv())
}
P <- lapply(c(
  f1d = file.path(NAT, "fig1d_mito_content.R"),     f1e = file.path(NAT, "fig1e_myc_reallocation.R"),
  f1f = file.path(PROP, "p1f_window_combined.R"),   f1g = file.path(PROP, "p1g_heatmap.R"),
  f1h = file.path(PROP, "p1h_ratio_expected.R"),    f1i = file.path(PROP, "p1i_coupling_adjusted.R"),
  e2a = file.path(NAT, "edfig2a_design.R"),         e2b = file.path(NAT, "edfig2b_geneset_library.R"),
  e2c = file.path(PROP, "ped2c_pathway_axis_wide.R"), e2d = file.path(NAT, "edfig2d_window_independent.R"),
  e2e = file.path(NAT, "edfig2e_myc_reallocation_ages.R"),
  e2f = file.path(NAT, "edfig2g_window_two_rulers.R"),          # the old 2g
  e3a = file.path(PROP, "ped3a_per_animal.R"),      e3b = file.path(PROP, "ped3b_fade_ranking.R"),
  e3c = file.path(PROP, "ped3c_puma_regulators.R"), e3d = file.path(PROP, "p1i_coupling.R")), grab)

PW <- 210; PH <- 297; FY <- 22
mm <- function(v) unit(v, "mm")
hy <- function(x) gsub("(?<=[A-Za-z0-9])-(?=[A-Za-z0-9])", "­", x, perl = TRUE)
page <- function(title, W, H = 247) {
  grid.newpage()
  fx <- (PW - W) / 2
  grid.rect(x = mm(fx), y = mm(PH - FY), width = mm(W), height = mm(H), just = c("left", "top"),
            gp = gpar(fill = NA, col = "grey70", lwd = 0.4, lty = "dotted"))
  grid.text(hy(title), x = mm(fx), y = mm(PH - 8), just = c("left", "top"),
            gp = gpar(fontsize = 10, fontface = "bold", fontfamily = "Helvetica"))
  fx
}
place <- function(fx, p, x, y, w, h, letter) {
  print(p, vp = viewport(x = mm(fx + x), y = mm(PH - FY - y), width = mm(w), height = mm(h),
                         just = c("left", "top")), newpage = FALSE)
  grid.text(letter, x = mm(fx + x - 1.2), y = mm(PH - FY - y + 1.2), just = c("left", "top"),
            gp = gpar(fontsize = 8, fontface = "bold", fontfamily = "Helvetica"))
}

pdf(OUT, width = PW / 25.4, height = PH / 25.4, family = "Helvetica", useDingbats = FALSE)

# --- Fig. 1 -----------------------------------------------------------------------
fx <- page("Fig. 1 on 183 x 247 mm (proposal; the bench row a-c is a placeholder)", 183)
grid.rect(x = mm(fx), y = mm(PH - FY), width = mm(183), height = mm(56), just = c("left", "top"),
          gp = gpar(fill = "grey95", col = NA))
grid.text(hy("a-c  in vivo panels (bench): tumour onset, Ki-67 / cleaved caspase 3, MYC western + Myc mRNA, MYCER(T2)"),
          x = mm(fx + 91.5), y = mm(PH - FY - 28), gp = gpar(fontsize = 7, col = "grey40", fontfamily = "Helvetica"))
for (k in 0:2) grid.text(letters[k + 1], x = mm(fx - 1.2 + k * 61), y = mm(PH - FY + 1.2),
                         just = c("left", "top"), gp = gpar(fontsize = 8, fontface = "bold"))
place(fx, P$f1d, 0, 59, 100, 66, "d");    place(fx, P$f1e, 103, 59, 80, 66, "e")
place(fx, P$f1f, 0, 128, 100, 56, "f");   place(fx, P$f1g, 103, 128, 80, 56, "g")
place(fx, P$f1h, 0, 187, 100, 60, "h");   place(fx, P$f1i, 103, 187, 80, 60, "i")

# --- ED Fig. 2 --------------------------------------------------------------------
fx <- page("Extended Data Fig. 2 on 180 x 247 mm (proposal; old f in Fig. 1f, old g is f)", 180)
place(fx, P$e2a, 0, 0, 46, 40, "a");      place(fx, P$e2b, 50, 0, 58, 40, "b")
place(fx, P$e2c, 0, 43, 176, 46, "c")
place(fx, P$e2d, 0, 92, 114, 56, "d");    place(fx, P$e2e, 117, 92, 60, 56, "e")
place(fx, P$e2f, 0, 151, 84, 64, "f")

# --- ED Fig. 3 --------------------------------------------------------------------
fx <- page("Extended Data Fig. 3 (new) on 180 x 247 mm (proposal)", 180)
place(fx, P$e3a, 0, 0, 104, 60, "a");     place(fx, P$e3d, 107, 0, 73, 60, "d")
place(fx, P$e3b, 0, 63, 86, 100, "b");    place(fx, P$e3c, 92, 63, 86, 92, "c")

invisible(dev.off())
cat("wrote", OUT, "\n")
