# =============================================================================
# panels_to_pdf.R -- one PDF, one panel per page, for reading alongside the text
# -----------------------------------------------------------------------------
# A circulation document, NOT a figure. Collaborators read the Results in the
# Google Doc and want to see the panel each sentence cites without assembling
# anything. So: every built panel of Figure 1 and Supplementary 1, in the order
# the Results cites them, one to a page, with the slot as the page title.
#
# USE
#   source(here::here("figures", "panels", "panels_to_pdf.R"))
#
# WHY IT RE-RENDERS RATHER THAN CONCATENATING THE EXISTING PDFs. There is no
# Ghostscript, qpdf, pdftk, pdfjam or LaTeX on this machine and neither the qpdf
# nor the pdftools R package is installed, so nothing here can place an existing
# PDF page inside another document. Re-rendering from the scripts is also the
# better answer: it stays vector, it cannot go stale against the scripts, and the
# titles come from the panel_legend() blocks rather than from a list kept by hand.
#
# EVERY PANEL IS DRAWN AT ITS DESIGNED PHYSICAL SIZE, centred on a page just large
# enough for the largest of them. Scaling a ggplot to fill a page does not
# magnify it -- the type stays at 6 pt while the plot area stretches, which is a
# different picture from the one that was signed off. At true size the page looks
# exactly like the manuscript figure will, and a reader magnifies with the
# viewer's zoom, which is what zoom is for.
#
# The size of each panel is READ FROM ITS OWN save_panel_p() CALL through the
# myc.fig.capture hook in _panel_common.R, so this file keeps no second copy of a
# number that lives in the panel script.
#
# Output: outputs/figures/panels/Figure1_and_S1_panels.pdf
# =============================================================================

if (!requireNamespace("here", quietly = TRUE)) stop("panels_to_pdf needs 'here'")

source(here::here("figures", "panels", "_panel_common.R"))

panels_dir <- here::here("figures", "panels")
out_file   <- file.path(here::here("outputs", "figures", "panels"),
                        "Figure1_and_S1_panels.pdf")

MARGIN    <- 8    # mm, around everything
TITLE_H   <- 11   # mm, the strip the page title sits in
TITLE_PT  <- 11

# --- collect every panel, its object, its slot and its designed size ----------
# One environment per script, as rebuild_panels.R does, so a panel object from one
# cannot quietly satisfy the next.
scripts <- sort(list.files(panels_dir, pattern = "^fig.*\\.R$", full.names = FALSE))
scripts <- setdiff(scripts, basename(c("panels_to_pdf.R")))
stopifnot(length(scripts) > 0L)

if (exists(".myc_panel_sizes", envir = globalenv())) {
  rm(".myc_panel_sizes", envir = globalenv())
}
# NO on.exit() ANYWHERE IN THIS FILE. At the top level of a sourced script
# on.exit() attaches to the frame evaluating that single expression, so it fires
# at once -- here it restored the option before a single panel had run, and every
# panel wrote its own PDF instead of reporting its size. Options and the device
# are restored explicitly at the foot.
old_opts <- options(myc.fig.capture = TRUE)

cat("\n", strrep("=", 74), "\nPANELS TO PDF\n", strrep("=", 74), "\n\n", sep = "")

panels <- list()
for (f in scripts) {
  slug <- sub("\\.R$", "", f)
  cat(sprintf("  %-34s ", f)); utils::flush.console()
  env <- new.env(parent = globalenv())
  ok <- tryCatch({ source(file.path(panels_dir, f), local = env); TRUE },
                 error = function(e) { cat("FAILED: ", conditionMessage(e), "\n"); FALSE })
  if (!ok) next
  if (!exists("p", envir = env, inherits = FALSE)) { cat("no `p`\n"); next }
  reg  <- get(".myc_panel_sizes", envir = globalenv())
  size <- get0(slug, envir = reg, ifnotfound = NULL)
  slot <- if (exists("LEGEND", envir = env, inherits = FALSE))
    get("LEGEND", envir = env)$slot else NA_character_
  if (is.null(size) || is.na(slot)) { cat("no size or no LEGEND slot\n"); next }
  stopifnot(identical(size$units, "mm"))
  panels[[slug]] <- list(p = get("p", envir = env), slot = slot,
                         w = size$width, h = size$height)
  cat(sprintf("%-9s %3.0f x %3.0f mm\n", slot, size$width, size$height))
}
options(old_opts)

# --- narrative order, and what is deliberately left out -----------------------
# slug_slot_order() puts main figures before supplementary, then number, then
# letter -- which for Figure 1 and Supplementary 1 IS the order the Results cites
# them, because the lettering was assigned from the written text.
keep <- vapply(panels, function(x) grepl("^Fig\\. S?1", x$slot), logical(1))
dropped <- names(panels)[!keep]
panels  <- panels[keep]
ord <- slug_slot_order(names(panels),
                       lapply(panels, function(x) list(slot = x$slot)))
panels <- panels[ord]

if (length(dropped)) {
  cat("\n  left out (no Figure 1 / S1 slot): ", paste(dropped, collapse = ", "), "\n", sep = "")
}
stopifnot(length(panels) > 0L)

# --- the page ----------------------------------------------------------------
PAGE_W <- max(vapply(panels, function(x) x$w, numeric(1))) + 2 * MARGIN
PAGE_H <- max(vapply(panels, function(x) x$h, numeric(1))) + TITLE_H + 2 * MARGIN

mm2in <- function(x) x / 25.4

grDevices::pdf(file = out_file, width = mm2in(PAGE_W), height = mm2in(PAGE_H),
               family = "Helvetica", useDingbats = FALSE, onefile = TRUE)

for (nm in names(panels)) {
  x <- panels[[nm]]
  grid::grid.newpage()

  # title strip, flush across the page
  grid::pushViewport(grid::viewport(
    y = grid::unit(PAGE_H - MARGIN, "mm"), height = grid::unit(TITLE_H, "mm"),
    just = "top", default.units = "mm"))
  grid::grid.text(x$slot, x = grid::unit(0.5, "npc"), y = grid::unit(0.5, "npc"),
                  gp = grid::gpar(fontsize = TITLE_PT, fontface = "bold"))
  grid::popViewport()

  # the panel itself, at its designed size, centred in what is left
  body_h <- PAGE_H - TITLE_H - 2 * MARGIN
  grid::pushViewport(grid::viewport(
    x = grid::unit(PAGE_W / 2, "mm"),
    y = grid::unit(MARGIN + body_h / 2, "mm"),
    width  = grid::unit(x$w, "mm"),
    height = grid::unit(x$h, "mm"),
    default.units = "mm"))
  print(x$p, vp = grid::viewport(width = grid::unit(1, "npc"),
                                 height = grid::unit(1, "npc")), newpage = FALSE)
  grid::popViewport()
}

grDevices::dev.off()

cat(sprintf("\n  wrote %s\n  %d page(s), %.0f x %.0f mm, panels at true size\n\n",
            out_file, length(panels), PAGE_W, PAGE_H))
cat(strrep("=", 74), "\n\n", sep = "")

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {

  ## what went in, in the order it went in
  data.frame(slot = vapply(panels, function(x) x$slot, character(1)),
             slug = names(panels),
             w = vapply(panels, function(x) x$w, numeric(1)),
             h = vapply(panels, function(x) x$h, numeric(1))) |>
    print(row.names = FALSE)

  ## the size registry the capture hook filled
  as.list(get(".myc_panel_sizes", envir = globalenv())) |> str()
}
