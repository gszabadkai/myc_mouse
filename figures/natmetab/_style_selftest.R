# =============================================================================
# _style_selftest.R -- proves that _style.R draws a good panel and refuses bad ones
# -----------------------------------------------------------------------------
# Synthetic data only, writes nothing. Run after any change to _style.R:
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
    geom_boxplot(aes(fill = group), outlier.shape = NA, width = 0.6,
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
expect_refused("a shrunken superscript", "superscript", p)

p <- good() + theme(axis.text = element_text(size = 4.5))
expect_refused("theme text under 5 pt", "theme axis.text", p)

p <- good() + scale_y_continuous(limits = c(1, 16), breaks = seq(0, 16, 4), expand = c(0, 0))
expect_refused("a magnitude axis not starting at 0", "starts at 1, not 0|does not end on a tick", p)

p <- good() + scale_y_continuous(limits = c(0, 15), breaks = seq(0, 12, 4), minor_breaks = NULL,
                                 expand = c(0, 0))
expect_refused("an axis not ending on a tick", "does not end on a tick", p)

p <- good() + scale_y_continuous()       # ggplot's default padding below 0
expect_refused("default padding", "does not end on a tick|must start at 0", p)

expect_ok("an axis may end on a drawn minor tick", check_panel(good(minor = TRUE), 50, 50))

p <- good() + scale_y_continuous(limits = c(0, 14), breaks = seq(0, 12, 4), expand = c(0, 0))  # 14: an UNDRAWN minor break
expect_refused("an axis ending on an undrawn minor break", "does not end on a tick", p)

p <- good() + scale_y_log10()
expect_refused("a log axis under a drawn x axis", "has no 0", p)

expect_refused("a panel taller than the page", "over the 170 mm page", good(), height = 180)

cat(if (n_fail == 0L) "\nall checks behave as specified\n" else sprintf("\n%d FAILED\n", n_fail))
if (n_fail > 0L) stop("style self-test failed", call. = FALSE)
