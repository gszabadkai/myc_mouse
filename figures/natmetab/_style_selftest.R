# =============================================================================
# _style_selftest.R -- proves that _style.R draws a good panel and refuses bad ones
# -----------------------------------------------------------------------------
# Synthetic data only; writes nothing but one temporary PDF, deleted (the 1-pt line
# is read from the PDF itself). Run after any change to _style.R:
#   Rscript figures/natmetab/_style_selftest.R
# Every case builds a FRESH plot: ggplot layers are shared by reference, so
# editing a copy of a plot also edits the original.
# =============================================================================

source(here::here("figures", "natmetab", "_style.R"))

d <- data.frame(group = factor(rep(order_by_age, each = 6), levels = order_by_age),
                y = c(2.6, 2.4, 2.9, 2.2, 2.7, 2.5, 7.9, 6.4, 9.1, 7.2, 8.3, 6.9,
                      2.1, 2.3, 1.9, 2.4, 2.2, 2.0, 7.6, 8.0, 7.7, 8.4, 7.4, 5.1))
comps <- data.frame(group1 = c("6W_neg", "12W_neg", "6W_neg", "6W_pos"),
                    group2 = c("6W_pos", "12W_pos", "12W_neg", "12W_pos"),
                    p = c(2.653779e-17, 2.143239e-20, 0.484, 0.764), level = c(1, 1, 2, 3))

good <- function(minor = FALSE) {
  br <- brackets(comps, base = 9.8, step = 1.2)
  ggplot(d, aes(group, y)) +
    geom_boxplot(aes(fill = group), outlier.shape = NA, width = 0.6, coef = Inf,
                 colour = box_line, linewidth = NM_LINE) +
    geom_point(aes(colour = group), position = position_jitter(width = 0.15, seed = 1),
               size = 0.9, shape = 16) +
    br +
    scale_colour_manual(values = group_cols) +
    scale_fill_manual(values = group_fill) +
    scale_x_discrete(labels = group_labels) +
    scale_y_nm(step = 4, top = attr(br, "top"), minor = minor) +
    labs(x = NULL, y = "Level (units)") +
    theme_nm()
}

n_fail <- 0L
expect_ok <- function(what, expr) {
  res <- tryCatch({ force(expr); "ok" }, error = function(e) conditionMessage(e))
  pass <- identical(res, "ok")
  cat(sprintf("[%s] %s%s\n", if (pass) "PASS" else "FAIL", what, if (pass) "" else paste(":", res)))
  if (!pass) n_fail <<- n_fail + 1L
}
expect_refused <- function(what, pattern, p, width = 50, height = 50) {
  res <- tryCatch({ check_panel(p, width, height); "not refused" },
                  error = function(e) conditionMessage(e))
  pass <- grepl(pattern, res)
  cat(sprintf("[%s] refuses %s\n", if (pass) "PASS" else "FAIL", what))
  if (!pass) { cat("       got:", res, "\n"); n_fail <<- n_fail + 1L }
}

# --- formatting -------------------------------------------------------------------
expect_ok("P labels format as specified", stopifnot(identical(
  fmt_p(c(0.484, 0.05, 1, 0.0074, 2.653779e-17, 0.000999)),
  c("italic(P)=='0.48'", "italic(P)=='0.050'", "italic(P)=='1.0'", "italic(P)=='0.0074'",
    "italic(P)=='2.7'%*%10^displaystyle('-17')", "italic(P)=='1'%*%10^displaystyle('-3')"))))
expect_ok("a reduced exponent is named scriptstyle", stopifnot(identical(
  fmt_p(c(2.653779e-17, 0.0074), exponent = "reduced"),
  c("italic(P)=='2.7'%*%10^scriptstyle('-17')", "italic(P)=='0.0074'"))))
expect_ok("tick labels: commas, '0', no trailing zeros", stopifnot(
  identical(lab_num(c(0, 4000, 16000)), c("0", "4,000", "16,000")),
  identical(lab_signed(c(-0.5, 0, 0.25)), c("-0.5", "0", "+0.25"))))

# --- a good panel passes every check and draws ---------------------------------------
expect_ok("a good panel passes and draws", {
  old <- options(natmetab.dry_run = TRUE)
  on.exit(options(old))
  suppressMessages(save_panel(good(), "Fig1", "z", "selftest", 50, 50))
})

# --- each rule refuses what it should ------------------------------------------------
p <- good(); p$layers[[3]]$aes_params$label.size <- pt(4.5)
expect_refused("text under 5 pt", "text at 4.50 pt", p)

p <- good() + annotate("text", x = 1, y = 1, label = "note", size = pt(6), colour = "red")
expect_refused("coloured text", "coloured text", p)

p <- good(); p$layers[[3]]$data$label <- "italic(P)==2.7%*%10^-17"
expect_refused("a superscript in no named style", "superscript", p)

p <- good() + theme(axis.text = element_text(size = 4.5))
expect_refused("theme text under 5 pt", "theme axis.text", p)

d10 <- transform(d, y = y + 10)        # every value above 10
p <- ggplot(d10, aes(group, y)) + geom_point() + scale_y_nm(step = 4, top = 20, bottom = 8) + theme_nm()
expect_refused("a magnitude axis not starting at 0", "starts at 8, not 0", p)
expect_ok("a box plot's y axis may start above 0, on the ticks box_axis() finds", {
  a <- box_axis(min(d10$y), max(d10$y), step = 4, graph = 40)
  stopifnot(isTRUE(all.equal(a$limits, c(10, 20))))
  check_panel(ggplot(d10, aes(group, y)) + geom_boxplot(coef = Inf) +
                scale_y_continuous(limits = a$limits, breaks = a$breaks, minor_breaks = a$minor_breaks,
                                   guide = a$guide, expand = c(0, 0)) + theme_nm(), 50, 50)
})
expect_ok("box_axis: the closest ticks outside the whiskers, its margins kept in mm", {
  a <- box_axis(lo = 0.31, hi = 0.94, step = 0.2, graph = 25, above = 3)
  stopifnot(isTRUE(all.equal(a$limits, c(0.3, 1.1))))
  b <- box_axis(lo = 0.31, hi = 0.94, step = 0.2, graph = 25, below = 0.5, above = 3)
  stopifnot(isTRUE(all.equal(b$limits, c(0.2, 1.1))))     # 0.3 leaves under 0.5 mm below 0.31
})

p <- good() + scale_y_continuous(limits = c(0, 15), breaks = seq(0, 12, 4), minor_breaks = NULL,
                                 expand = c(0, 0))
expect_refused("an axis not ending on a tick", "does not end on a tick", p)

p <- good() + scale_y_continuous()       # ggplot's default padding below 0
expect_refused("default padding", "does not end on a tick|must start at 0", p)

expect_ok("an axis may end on a drawn minor tick", check_panel(good(minor = TRUE), 50, 50))
expect_ok("a schematic with no drawn axes need not end on ticks", check_panel(
  ggplot(data.frame(x = c(0.3, 4.7), y = c(0.2, 2.9)), aes(x, y)) + geom_point() +
    theme_void(base_family = NM_FONT), 50, 40))
expect_refused("an Extended Data panel wider than 180 mm", "over the 180 mm page", {
  old <- options(natmetab.dry_run = TRUE); on.exit(options(old))
  suppressMessages(save_panel(good(), "EDFig2", "z", "selftest_wide", 182, 50))
})
expect_ok("tick-ending axes survive floating point (-0.2 + 3 * 0.2)", {
  q <- ggplot(data.frame(x = c(-0.17, 0.31), y = c(1, 2)), aes(x, y)) + geom_point() +
    scale_x_nm(step = 0.2, top = 0.31, bottom = -0.17, minor = TRUE) +
    scale_y_nm(step = 1, top = 2) + theme_nm()
  check_panel(q, 50, 50)
})

p <- good() + scale_y_continuous(limits = c(0, 14), breaks = seq(0, 12, 4), expand = c(0, 0))  # 14: an UNDRAWN minor break
expect_refused("an axis ending on an undrawn minor break", "does not end on a tick", p)

p <- good() + scale_y_log10()
expect_refused("a log axis under a drawn x axis", "has no 0", p)

expect_refused("a panel taller than the page", "over the 170 mm page", good(), height = 180)

p <- good(); p$layers[[1]]$stat_params$coef <- 1.5          # Tukey whiskers
expect_refused("box whiskers that do not run min to max", "whiskers must run min to max", p)

# --- the Prism format (Fig. 1b): 1-pt lines, brackets without tips, a fixed graph ----
tipless <- function(txt = 5) {
  br <- brackets(comps, base = 9.8, step = 1.2, txt = txt, line = lw_pt(1), tips = FALSE)
  ggplot(d, aes(group, y)) + geom_point() + br +
    scale_y_nm(step = 4, top = attr(br, "top")) + theme_nm(txt = 5, txt_title = 6, line = lw_pt(1))
}
expect_ok("lw_pt(1) prints a 1-pt line", {
  f <- tempfile(fileext = ".pdf"); on.exit(unlink(f))
  grDevices::pdf(f, compress = FALSE)
  grid::grid.lines(gp = grid::gpar(lwd = lw_pt(1) * ggplot2::.pt))
  invisible(grDevices::dev.off())
  stopifnot(any(grepl("^1.00 w$", readLines(f, warn = FALSE), useBytes = TRUE)))
})
expect_ok("a bracket without tips is its bar alone, with its P label", {
  g    <- ggplot2::layer_grob(tipless(), 2)[[1]]$children
  segs <- Filter(function(x) inherits(x, "segments"), g)
  labs <- Filter(function(x) inherits(x, "text"), g)
  stopifnot(length(segs) == nrow(comps), all(vapply(segs, function(s) length(s$x0), 1L) == 1L),
            length(labs) == nrow(comps), all(vapply(labs, function(t) is.expression(t$label), TRUE)))
})
expect_refused("a bracket without tips under 5 pt", "text at 4.50 pt", tipless(txt = 4.5))

# --- a reduced exponent (Fig. 1d): passes, and prints at 4 pt in a 5-pt label ----------
expect_ok("brackets with a reduced exponent pass", check_panel(
  ggplot(d, aes(group, y)) + geom_point() +
    brackets(comps, base = 9.8, step = 1.2, txt = 5, exponent = "reduced") +
    scale_y_nm(step = 4, top = 14) + theme_nm(txt = 5, txt_title = 6), 50, 50))
expect_ok("a reduced exponent prints at 4 pt in a 5-pt label (70%, rounded by pdf())", {
  f <- tempfile(fileext = ".pdf"); on.exit(unlink(f))
  grDevices::pdf(f, compress = FALSE, family = NM_FONT)
  grid::grid.text(parse(text = fmt_p(2.653779e-17, exponent = "reduced")), gp = grid::gpar(fontsize = 5))
  invisible(grDevices::dev.off())
  tm <- grep("Tm \\((2.7|-17)\\) Tj", readLines(f, warn = FALSE), value = TRUE, useBytes = TRUE)
  stopifnot(length(tm) == 2L, grepl("Tf 5.00 .*\\(2.7\\)", tm[1]), grepl("Tf 4.00 .*\\(-17\\)", tm[2]))
})
expect_ok("label_mm: 5-pt P labels measure 1.30 mm, 1.90 with a reduced exponent", {
  h <- label_mm(fmt_p(c(0.0030, 4.9e-4), exponent = "reduced"), txt = 5)
  stopifnot(abs(h - c(1.30, 1.90)) < 0.02)
})
expect_refused("a fixed-size panel saved smaller than it needs", "edges would be cut", {
  old <- options(natmetab.dry_run = TRUE); on.exit(options(old))
  suppressMessages(save_panel(fixed_size(good(), 30, 30), "Fig1", "z", "selftest_fixed", 30, 30))
})
expect_ok("natural_size measures a fixed stack with a collected key", {
  n <- natural_size(stack_parts(good(), good(), heights = grid::unit(c(30, 10), "mm"), collect_guides = TRUE))
  stopifnot(is.na(n[["width"]]), n[["height"]] > 40, n[["height"]] < 70)
})
expect_ok("a fixed-size panel saved at its natural size passes", {
  old <- options(natmetab.dry_run = TRUE); on.exit(options(old))
  q <- fixed_size(good(), 30, 30)
  n <- ceiling(natural_size(q) * 2) / 2
  suppressMessages(save_panel(q, "Fig1", "z", "selftest_fixed", n[["width"]], n[["height"]]))
})

# --- stacked parts: each part is checked on its own ----------------------------------
strip <- function(size = pt(6)) {
  ggplot(d, aes(y, group)) + geom_point(size = 0.6) +
    annotate("text", x = 1, y = 1, label = "x", size = size) +
    scale_x_nm(step = 4, top = 16) + theme_nm()
}
expect_ok("a stack of good parts passes and draws", {
  old <- options(natmetab.dry_run = TRUE)
  on.exit(options(old))
  suppressMessages(save_panel(stack_parts(good(), strip(), heights = c(3, 1),
                                          collect_guides = TRUE),
                              "Fig1", "z", "selftest_stack", 50, 60))
})
expect_ok("two parts side by side pass and draw", {
  old <- options(natmetab.dry_run = TRUE)
  on.exit(options(old))
  suppressMessages(save_panel(stack_parts(good(), good(), ncol = 2, widths = c(1, 1)),
                              "Fig1", "z", "selftest_side", 100, 50))
})
res <- tryCatch({ save_panel(stack_parts(good(), strip(size = pt(4)), heights = c(3, 1)),
                             "Fig1", "z", "selftest_stack", 50, 60); "not refused" },
                error = function(e) conditionMessage(e))
pass <- grepl("text at 4.00 pt", res)
cat(sprintf("[%s] refuses a stack with one bad part\n", if (pass) "PASS" else "FAIL"))
if (!pass) { cat("       got:", res, "\n"); n_fail <- n_fail + 1L }


cat(if (n_fail == 0L) "\nall checks behave as specified\n" else sprintf("\n%d FAILED\n", n_fail))
if (n_fail > 0L) stop("style self-test failed", call. = FALSE)
