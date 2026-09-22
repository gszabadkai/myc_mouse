# =============================================================================
# rebuild_manuscript_figures.R -- re-render the four manuscript figures
# -----------------------------------------------------------------------------
# The autorun script. Sources the four manuscript figure scripts in order and
# writes their PDFs to outputs/figures/. Nothing is re-analysed: every panel is
# read from results/*.rds exactly as the figure scripts already do.
#
# USE
#
#   source(here::here("figures", "rebuild_manuscript_figures.R"))
#
#   options(myc.fig.nosave = TRUE)   # dry run: executes AND draws every panel
#                                    # to a null device, writes nothing. Same
#                                    # switch paper/myc_mito.qmd uses, so a check
#                                    # and a real pass differ by one line.
#
# WHAT THIS IS *NOT*
#
# It is not a pipeline runner. After a restart nothing needs re-running at all --
# results/*.rds are files, not session state. Run this only when a figure script
# has changed, or to confirm the panels still resolve against the current
# results/. Note that `quarto render paper/myc_mito.qmd` also re-executes all
# eighteen figure scripts, but its setup chunk sets myc.fig.nosave = TRUE, so it
# renders inline and writes no PDFs. This script is how outputs/figures/ is
# refreshed.
#
# The disaster-recovery chain (only if results/ is lost) is at the foot of this
# file, behind an explicit opt-in. It is recorded in executable form so the order
# is not re-derived under pressure -- not so it gets run by accident.
#
# Each script runs in its own environment, so a panel object from one cannot
# quietly satisfy the next. theme_myc.R and the reconciler still land in the
# global environment (they are sourced from inside each script), which is what
# you want for poking at the result afterwards.
# =============================================================================

if (!requireNamespace("here", quietly = TRUE)) stop("rebuild needs the 'here' package")

fig_scripts <- c("figure1_myc_mitochondrion.R",
                 "figure2_developmental_window.R",
                 "figureS1_compartment_detail.R",
                 "figureS2_controls.R")

fig_pdfs <- c("figure1_myc_mitochondrion.pdf",
              "figure2_developmental_window.pdf",
              "figureS1_compartment_detail.pdf",
              "figureS2_controls.pdf")

dry_run <- isTRUE(getOption("myc.fig.nosave", FALSE))
out_dir <- here::here("outputs", "figures")

cat("\n", strrep("=", 78), "\n", sep = "")
cat("REBUILD MANUSCRIPT FIGURES", if (dry_run) "  [DRY RUN -- writes nothing]" else "",
    "\n", sep = "")
cat(strrep("=", 78), "\n\n", sep = "")

before <- vapply(file.path(out_dir, fig_pdfs),
                 function(p) if (file.exists(p)) file.mtime(p) else NA_real_,
                 numeric(1))

status <- character(length(fig_scripts))
elapsed <- numeric(length(fig_scripts))

for (i in seq_along(fig_scripts)) {
  f <- fig_scripts[i]
  cat(sprintf("[%d/%d] %s ... ", i, length(fig_scripts), f))
  utils::flush.console()
  t0 <- Sys.time()
  ok <- tryCatch({
    env <- new.env(parent = globalenv())
    source(here::here("figures", f), local = env)
    # A dry run skips save_panel(), and save_panel() is where the plot is
    # actually DRAWN -- patchwork assembly, guides, clipping and font embedding
    # all happen at render, not at construction. Without this the dry run would
    # only prove the data code paths work. pdf(NULL) is a null device: it
    # exercises the whole draw pipeline and writes nothing.
    if (dry_run && exists("p", envir = env, inherits = FALSE)) {
      grDevices::pdf(NULL)
      on.exit(grDevices::dev.off(), add = TRUE)
      print(get("p", envir = env))
      grDevices::dev.off()
      on.exit()
    }
    TRUE
  }, error = function(e) {
    status[i] <<- conditionMessage(e)
    FALSE
  })
  elapsed[i] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (ok) {
    status[i] <- "ok"
    cat(sprintf("ok (%.1f s)\n", elapsed[i]))
  } else {
    cat(sprintf("FAILED (%.1f s)\n         %s\n", elapsed[i], status[i]))
  }
}

# --- did the PDFs actually move? ---------------------------------------------
cat("\n", strrep("-", 78), "\n", sep = "")
after <- vapply(file.path(out_dir, fig_pdfs),
                function(p) if (file.exists(p)) file.mtime(p) else NA_real_,
                numeric(1))
for (i in seq_along(fig_pdfs)) {
  written <- !is.na(after[i]) && (is.na(before[i]) || after[i] > before[i])
  cat(sprintf("  %-38s %s\n", fig_pdfs[i],
              if (dry_run) "not written (dry run)"
              else if (written) paste0("written ",
                                       format(as.POSIXct(after[i], origin = "1970-01-01"),
                                              "%H:%M:%S"))
              else if (is.na(after[i])) "MISSING -- script did not save"
              else "unchanged -- script did not save"))
}

n_fail <- sum(status != "ok")
cat(sprintf("\n%d/%d scripts ok, %.1f s total\n",
            length(status) - n_fail, length(status), sum(elapsed)))
if (n_fail > 0) {
  cat("\nFAILURES:\n")
  for (i in which(status != "ok")) cat("  ", fig_scripts[i], ": ", status[i], "\n", sep = "")
}
cat(strrep("=", 78), "\n\n", sep = "")

# =============================================================================
# DISASTER RECOVERY -- only if results/ is lost or corrupt
# -----------------------------------------------------------------------------
# Inert by default. This regenerates the eleven objects the whole figures/ layer
# reads, and nothing else. Do NOT run it to "refresh" anything:
#
#   * it is hours of compute for objects that already exist on disk;
#   * it overwrites objects whose numbers are quoted verbatim in
#     paper/myc_mito.qmd, in CLAUDE.md and in the dated docs;
#   * not every script seeds its randomness (21_ap6_permutation_null.R calls
#     sample() at line 125, BEFORE its set.seed(1) at line 135), and fgsea is
#     stochastic, so a re-run is not guaranteed bit-identical.
#
# Enable deliberately:  options(myc.rebuild.analysis = TRUE)
# =============================================================================

analysis_chain <- c(
  "scripts/00_setup_packages.R",
  "scripts/01_load_data.R",
  "scripts/03_deseq_results_qc.R",
  "scripts/archive_main_pipeline/02_deseq_interaction_model.R",  # cwd = repo root
  "scripts/08_mitoPPS_analysis.R",
  "scripts/32_mito_content_proxies.R",
  "scripts/40_background_vs_myc_decomposition.R",
  "scripts/42_priming_arm_and_teb_substrate.R",
  "scripts/43_substrate_specificity_and_tradeoff.R",
  "scripts/44_collapse_module_and_ownership.R")

if (isTRUE(getOption("myc.rebuild.analysis", FALSE))) {
  cat("myc.rebuild.analysis = TRUE. This rewrites results/*.rds. The order is:\n")
  for (s in analysis_chain) cat("  ", s, "\n", sep = "")
  cat("\nRun them BY HAND, one at a time, in Positron -- per the project's",
      "Option A workflow.\n  ",
      "archive_main_pipeline/02 uses bare relative paths: setwd(here::here()) first.\n")
}

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {

  ## dry run: every panel executes, outputs/ is untouched
  options(myc.fig.nosave = TRUE)
  source(here::here("figures", "rebuild_manuscript_figures.R"))

  ## real run
  options(myc.fig.nosave = NULL)
  source(here::here("figures", "rebuild_manuscript_figures.R"))

  ## one figure only
  source(here::here("figures", "figure2_developmental_window.R"))
}
