# =============================================================================
# rebuild_panels.R -- render every panel under figures/panels/
# -----------------------------------------------------------------------------
# Sources each figures/panels/fig*.R in name order, writes its PDF to
# outputs/figures/panels/, and collects the panel_legend() blocks into
# outputs/figures/panels/legends.md -- the raw material for the figure legends,
# which is where all explanatory text lives under the publication rule.
#
# USE
#
#   source(here::here("figures", "panels", "rebuild_panels.R"))
#
#   options(myc.fig.nosave = TRUE)   # dry run: executes AND draws every panel
#                                    # to a null device, writes nothing at all
#                                    # (no PDFs, no legends.md).
#
# Each script runs in its own environment, so a panel object from one cannot
# quietly satisfy the next -- the same guard figures/rebuild_manuscript_figures.R
# uses, and the reason both locate the composite by the fixed name `p`.
#
# A dry run still DRAWS. save_panel_p() is where patchwork assembly, guides,
# clipping and cairo font embedding actually execute; a run that only builds
# plot objects would test the data paths and nothing else. grDevices::pdf(NULL)
# is a null device: it exercises the whole draw pipeline and writes no bytes.
#
# This is not a pipeline runner. Panels read results/*.rds and render; nothing
# here re-runs analysis. Where a panel needs an analysis object refreshed it says
# so and stops -- see require_fresher_than() in _panel_common.R and the re-run
# chain in PANELS.md.
# =============================================================================

if (!requireNamespace("here", quietly = TRUE)) stop("rebuild_panels needs 'here'")

# Each panel sources this into its OWN environment, so legend_md() would not be
# visible out here; the runner needs its own copy to render legends.md.
source(here::here("figures", "panels", "_panel_common.R"))

panels_dir <- here::here("figures", "panels")
out_dir    <- here::here("outputs", "figures", "panels")
dry_run    <- isTRUE(getOption("myc.fig.nosave", FALSE))

scripts <- sort(list.files(panels_dir, pattern = "^fig.*\\.R$", full.names = FALSE))
if (length(scripts) == 0) stop("no panel scripts found in ", panels_dir)

cat("\n", strrep("=", 78), "\n", sep = "")
cat("REBUILD PANELS", if (dry_run) "  [DRY RUN -- writes nothing]" else "", "\n", sep = "")
cat(strrep("=", 78), "\n\n", sep = "")

slugs  <- sub("\\.R$", "", scripts)
before <- vapply(file.path(out_dir, paste0(slugs, ".pdf")),
                 function(p) if (file.exists(p)) file.mtime(p) else NA_real_,
                 numeric(1))

status   <- character(length(scripts))
elapsed  <- numeric(length(scripts))
legends  <- vector("list", length(scripts))
names(legends) <- slugs

for (i in seq_along(scripts)) {
  cat(sprintf("[%d/%d] %s ... ", i, length(scripts), scripts[i]))
  utils::flush.console()
  t0 <- Sys.time()
  ok <- tryCatch({
    env <- new.env(parent = globalenv())
    source(file.path(panels_dir, scripts[i]), local = env)
    if (dry_run && exists("p", envir = env, inherits = FALSE)) {
      grDevices::pdf(NULL)
      on.exit(grDevices::dev.off(), add = TRUE)
      print(get("p", envir = env))
      grDevices::dev.off()
      on.exit()
    }
    if (exists("LEGEND", envir = env, inherits = FALSE)) {
      legends[[i]] <- get("LEGEND", envir = env)
    }
    TRUE
  }, error = function(e) {
    status[i] <<- conditionMessage(e)
    FALSE
  })
  elapsed[i] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (ok) {
    status[i] <- "ok"
    cat(sprintf("ok (%.1f s)%s\n", elapsed[i],
                if (is.null(legends[[i]])) "  [no LEGEND block]" else ""))
  } else {
    cat(sprintf("FAILED (%.1f s)\n         %s\n", elapsed[i], status[i]))
  }
}

# --- did the PDFs actually move? ---------------------------------------------
cat("\n", strrep("-", 78), "\n", sep = "")
after <- vapply(file.path(out_dir, paste0(slugs, ".pdf")),
                function(p) if (file.exists(p)) file.mtime(p) else NA_real_,
                numeric(1))
for (i in seq_along(slugs)) {
  written <- !is.na(after[i]) && (is.na(before[i]) || after[i] > before[i])
  cat(sprintf("  %-40s %s\n", paste0(slugs[i], ".pdf"),
              if (dry_run) "not written (dry run)"
              else if (written) paste0("written ",
                                       format(as.POSIXct(after[i], origin = "1970-01-01"),
                                              "%H:%M:%S"))
              else if (is.na(after[i])) "MISSING -- panel did not save"
              else "unchanged -- panel did not save"))
}

# --- collect the legends -----------------------------------------------------
# Written only on a real run, so a dry run leaves the tree byte-identical.
have <- !vapply(legends, is.null, logical(1))
if (!dry_run && any(have)) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  md <- c(
    "# Panel legend blocks",
    "",
    paste0("Generated by `figures/panels/rebuild_panels.R`, ",
           format(Sys.time(), "%Y-%m-%d %H:%M"), ". Do not edit by hand -- edit the",
           " `panel_legend()` call in the panel script."),
    "",
    "Raw material for the figure legends. Nothing here is drawn on the panel:",
    "under the publication rule the figure carries data, axes and a key, and every",
    "explanation lives in the legend.",
    "")
  for (s in slugs[have]) md <- c(md, legend_md(legends[[s]], slug = s), "")
  writeLines(md, file.path(out_dir, "legends.md"))
  cat("\n  legends.md   ", sum(have), " block(s) -> ",
      file.path("outputs", "figures", "panels", "legends.md"), "\n", sep = "")
} else if (dry_run) {
  cat("\n  legends.md    not written (dry run); ", sum(have),
      " block(s) collected\n", sep = "")
}

n_fail <- sum(status != "ok")
cat(sprintf("\n%d/%d panels ok, %.1f s total\n",
            length(status) - n_fail, length(status), sum(elapsed)))
if (n_fail > 0) {
  cat("\nFAILURES:\n")
  for (i in which(status != "ok")) cat("  ", scripts[i], ": ", status[i], "\n", sep = "")
}
cat(strrep("=", 78), "\n\n", sep = "")

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {

  ## dry run: every panel executes and draws, outputs/ is untouched
  options(myc.fig.nosave = TRUE)
  source(here::here("figures", "panels", "rebuild_panels.R"))

  ## real run
  options(myc.fig.nosave = NULL)
  source(here::here("figures", "panels", "rebuild_panels.R"))

  ## one panel only
  source(here::here("figures", "panels", "figS1A_geneset_library.R"))
  print(p); print(LEGEND)
}
