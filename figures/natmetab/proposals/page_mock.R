# PROPOSALS (not panels of record) -- Fig. 1 page mock-ups: every panel drawn from its own script (dry run, nothing
# written by the panel scripts), placed at its printed size on a 183 x 170 mm page.
suppressPackageStartupMessages(library(grid))
SP   <- here::here("figures", "natmetab", "proposals")
NAT  <- here::here("figures", "natmetab")
OUT  <- here::here("outputs", "natmetab", "proposals", "Fig1_gh_proposals.pdf")

grab <- function(script, txt, ttl) {
  op <- options(natmetab.dry_run = TRUE, natmetab.capture = TRUE,
                natmetab.txt = txt, natmetab.txt_title = ttl)
  on.exit(options(op))
  if (exists(".nm_captured", envir = globalenv())) rm(".nm_captured", envir = globalenv())
  env <- new.env(parent = globalenv())
  invisible(capture.output(suppressMessages(source(script, local = env))))
  get(".nm_captured", envir = globalenv())
}
scripts <- c(d = file.path(NAT, "fig1d_mito_content.R"), e = file.path(NAT, "fig1e_myc_reallocation.R"),
             f = file.path(NAT, "fig1f_window_reallocation.R"), g = file.path(SP, "p1g_bcl2_balance.R"),
             h = file.path(SP, "p1h_puma_members.R"), halt = file.path(SP, "p1h_alt_bcl2_plane.R"),
             i = file.path(SP, "p1i_coupling.R"))
P6 <- lapply(scripts, grab, txt = 6, ttl = 7)
P5 <- lapply(scripts, grab, txt = 5, ttl = 6)

PW <- 210; PH <- 297                                # A4 sheets
FX <- (PW - 183) / 2; FY <- 22                      # figure area: 183 x 170 mm
mm <- function(v) unit(v, "mm")
hy <- function(x) gsub("(?<=[A-Za-z0-9])-(?=[A-Za-z0-9])", "\u00ad", x, perl = TRUE)
place <- function(p, x, y, w, h, letter = NULL) {
  vp <- viewport(x = mm(FX + x), y = mm(PH - FY - y), width = mm(w), height = mm(h),
                 just = c("left", "top"))
  print(p, vp = vp, newpage = FALSE)
  if (!is.null(letter))
    grid.text(letter, x = mm(FX + x - 1.2), y = mm(PH - FY - y + 1.2), just = c("left", "top"),
              gp = gpar(fontsize = 8, fontface = "bold", fontfamily = "Helvetica"))
}
bench <- function(h, label) {
  grid.rect(x = mm(FX), y = mm(PH - FY), width = mm(183), height = mm(h), just = c("left", "top"),
            gp = gpar(fill = "grey95", col = NA))
  grid.text(hy(label), x = mm(FX + 91.5), y = mm(PH - FY - h / 2),
            gp = gpar(fontsize = 7, col = "grey40", fontfamily = "Helvetica"))
  for (k in 0:2) grid.text(letters[k + 1], x = mm(FX - 1.2 + k * 61), y = mm(PH - FY + 1.2),
                           just = c("left", "top"),
                           gp = gpar(fontsize = 8, fontface = "bold", fontfamily = "Helvetica"))
}
frame <- function(title, lines) {
  grid.rect(x = mm(FX), y = mm(PH - FY), width = mm(183), height = mm(170), just = c("left", "top"),
            gp = gpar(fill = NA, col = "grey70", lwd = 0.4, lty = "dotted"))
  grid.text(hy(title), x = mm(FX), y = mm(PH - 8), just = c("left", "top"),
            gp = gpar(fontsize = 10, fontface = "bold", fontfamily = "Helvetica"))
  grid.text(hy(paste(lines, collapse = "\n")), x = mm(FX), y = mm(PH - FY - 174), just = c("left", "top"),
            gp = gpar(fontsize = 7.5, fontfamily = "Helvetica", lineheight = 1.15))
}

pdf(OUT, width = PW / 25.4, height = PH / 25.4, family = "Helvetica", useDingbats = FALSE)

# --- Layout A: house text (6 / 7 pt); Fig. 1e moves to Extended Data ------------
grid.newpage()
frame("Layout A -- house text size (6/7 pt); 1e moves to Extended Data",
      c("Dotted frame: the 183 x 170 mm Nature page. Grey band: the bench panels a-c (not drawn here; height assumed 52 mm).",
        "Row 2: content (d, 100 mm) beside the window panel (old 1f, narrowed to 80 mm). Row 3: the three new panels.",
        "Letters re-run after d: e = window, f = BCL-2 balance, g = four genes per animal, h = coupling.",
        "Every panel is the real script's output at the size shown; only the bench band is a placeholder."))
bench(52, "a-c  in vivo panels: tumour onset, Ki-67 / cleaved caspase 3, MYC western + Myc mRNA, MYCER(T2)")
place(P6$d, 0, 55, 100, 58, "d")
place(P6$f, 103, 55, 80, 58, "e")
place(P6$g, 0, 116, 50, 54, "f")
place(P6$h, 53, 116, 70, 54, "g")
place(P6$i, 126, 116, 57, 54, "h")

# --- Layout A2: house text; Fig. 1d moves to Extended Data -----------------------
grid.newpage()
frame("Layout A2 -- house text size (6/7 pt); 1d moves to Extended Data",
      c("Dotted frame: the 183 x 170 mm Nature page. Grey band: the bench panels a-c (height assumed 52 mm).",
        "Row 2: the MYC reallocation (old 1e, 80 mm) beside the window panel (old 1f) at its own 100 mm. Row 3: the three new panels.",
        "Letters re-run after c: d = MYC reallocation, e = window, f = BCL-2 balance, g = four genes, h = coupling.",
        "The content result (old 1d, 20-33%) would be cited from Extended Data."))
bench(52, "a-c  in vivo panels: tumour onset, Ki-67 / cleaved caspase 3, MYC western + Myc mRNA, MYCER(T2)")
place(P6$e, 0, 55, 80, 58, "d")
place(P6$f, 83, 55, 100, 58, "e")
place(P6$g, 0, 116, 50, 54, "f")
place(P6$h, 53, 116, 70, 54, "g")
place(P6$i, 126, 116, 57, 54, "h")

# --- Layout B: one point smaller (5 / 6 pt); every panel kept -------------------
grid.newpage()
frame("Layout B -- one point smaller (5/6 pt throughout); every panel kept",
      c("Dotted frame: the 183 x 170 mm Nature page. Grey band: the bench panels a-c (height assumed 50 mm).",
        "Row 2: d, e and f side by side (each narrowed from 100 / 83 / 100 mm). Row 3: the three new panels.",
        "Letters as in the Doc: g = BCL-2 balance, h = four genes per animal, i = coupling (new letter).",
        "5 pt is Nature's minimum, so nothing here needs redoing at acceptance."))
bench(50, "a-c  in vivo panels: tumour onset, Ki-67 / cleaved caspase 3, MYC western + Myc mRNA, MYCER(T2)")
place(P5$d, 0, 53, 66, 57, "d")
place(P5$e, 69, 53, 55, 57, "e")
place(P5$f, 127, 53, 56, 57, "f")
place(P5$g, 0, 113, 52, 57, "g")
place(P5$h, 55, 113, 70, 57, "h")
place(P5$i, 128, 113, 55, 57, "i")

# --- the new panels at 100%, house text, with the alternative 1h ----------------
grid.newpage()
grid.text("The new panels at 100% (house text, 6/7 pt)", x = mm(FX), y = mm(PH - 8), just = c("left", "top"),
          gp = gpar(fontsize = 10, fontface = "bold", fontfamily = "Helvetica"))
cap <- function(x, y, s) grid.text(hy(s), x = mm(FX + x), y = mm(PH - FY - y), just = c("left", "top"),
                                  gp = gpar(fontsize = 7, fontfamily = "Helvetica", lineheight = 1.1))
place(P6$g, 0, 8, 56, 56);   cap(0, 0, "1g  the BCL-2 balance (56 x 56 mm)")
place(P6$h, 70, 8, 72, 62);  cap(70, 0, "1h  four genes in every gland (72 x 62 mm)")
place(P6$i, 0, 86, 62, 56);  cap(0, 78, "1i  the coupling, both fits (62 x 56 mm)")
place(P6$halt, 70, 86, 60, 60); cap(70, 78, "1h, alternative: the BCL-2 family on both timelines (60 x 60 mm)")
cap(0, 156, paste(c(
  "1g  grey: the eight ratios the line is fitted to (slope 0.43, 95% prediction band); black: Bbc3:Bcl2l1, read against it.",
  "1h  lower brackets: 12W vs 6W in each genotype (DESeq2, IHW-adjusted P); upper bracket: the difference between the two",
  "      changes (interaction, raw P; Bbc3 pre-specified, Bax and Foxo3 exploratory). Each age was one extraction batch.",
  "1i  hypothesis-generating: the unadjusted difference clears its permutation null, the adjusted one does not.",
  "Alternative 1h: dashed, the same change in both genotypes; dotted, 0.20 either side (the declared threshold);",
  "      black, the four genes the text names; distance below the dashed line is the interaction."), collapse = "\n"))
invisible(dev.off())
cat("wrote", OUT, "\n")
