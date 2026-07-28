# =============================================================================
# 99_session_bootstrap.R -- come back to a working session after a restart
# -----------------------------------------------------------------------------
# NOT part of the analysis pipeline. Numbered 99 so it sorts clear of 00-44 and
# cannot be mistaken for a step. It is STRICTLY READ-ONLY: no saveRDS, no
# ggsave, no dir.create, nothing written anywhere.
#
# WHAT A RESTART ACTUALLY COSTS
#
# Nothing in results/. Those are files on disk, not session state, so a reboot
# cannot touch them and NO analysis script needs re-running. What is lost is the
# R session's memory and the terminal. This script restores the first in a few
# seconds and tells you the state of everything else.
#
# WHAT IT DOES
#
#   1. Checks the figure layer's packages are installed (reports, never installs).
#   2. Audits the 11 objects under results/ that the whole figures/ layer reads:
#      present? how big? how old? and is the producing SCRIPT newer than its own
#      output (i.e. edited but not re-run)?
#   3. Loads them into a single list S -- not loose globals, so nothing here can
#      shadow an object a figure script defines when you source one afterwards.
#   4. Prints the git and Google-Drive orientation you want before touching
#      anything: branch, HEAD, ahead/behind, dirty files, Icon stub counts.
#
# USE
#
#   source(here::here("scripts", "99_session_bootstrap.R"))
#
#   options(myc.bootstrap.full = TRUE)   # also source 00_setup_packages.R (slow,
#                                        # pulls DESeq2 etc.; only for analysis work)
#   options(myc.bootstrap.load = FALSE)  # audit only, load nothing
#
# The full pipeline stack is NOT loaded by default, on the same reasoning as
# figures/theme_myc.R: the figure layer only reads results/*.rds and renders.
#
# The producer map below is the closure of every readRDS() in figures/. If a
# figure script starts reading something new, add it here -- the audit is only
# as honest as this table.
# =============================================================================

if (!requireNamespace("here", quietly = TRUE)) stop("bootstrap needs the 'here' package")

boot_full <- isTRUE(getOption("myc.bootstrap.full", FALSE))
boot_load <- isTRUE(getOption("myc.bootstrap.load", TRUE))

root <- here::here()

# --- the 11 objects the figures/ layer reads, and what writes each ------------
# Two of these have NO writer in scripts/: combined_df_annotated{,_raw} come from
# the ARCHIVED main pipeline and are the one fragile link in the chain. They are
# gitignored, so they exist in exactly one place on disk.
boot_map <- data.frame(
  object = c("count_matrix",
             "dds_int_run",
             "interaction_results",
             "combined_df_annotated",
             "combined_df_annotated_raw",
             "mitopps_scores",
             "mito_content_proxies",
             "background_vs_myc",
             "priming_arm_teb",
             "substrate_specificity_tradeoff",
             "collapse_module_ownership"),
  producer = c("scripts/01_load_data.R",
               "scripts/03_deseq_results_qc.R",
               "scripts/03_deseq_results_qc.R",
               "scripts/archive_main_pipeline/02_deseq_interaction_model.R",
               "scripts/archive_main_pipeline/02_deseq_interaction_model.R",
               "scripts/08_mitoPPS_analysis.R",
               "scripts/32_mito_content_proxies.R",
               "scripts/40_background_vs_myc_decomposition.R",
               "scripts/42_priming_arm_and_teb_substrate.R",
               "scripts/43_substrate_specificity_and_tradeoff.R",
               "scripts/44_collapse_module_and_ownership.R"),
  s4 = c(FALSE, TRUE, TRUE, FALSE, FALSE,
         FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
  stringsAsFactors = FALSE
)

# --- packages the figure layer needs (checked, never installed) --------------
boot_pkgs <- c("here", "ggplot2", "dplyr", "tibble", "patchwork",
               "ggrepel", "ggbeeswarm", "DESeq2")

hr <- function(ch = "-") cat(strrep(ch, 78), "\n", sep = "")

cat("\n"); hr("=")
cat("SESSION BOOTSTRAP -- ", format(Sys.time(), "%Y-%m-%d %H:%M"), "\n", sep = "")
cat("repo: ", root, "\n", sep = "")
hr("=")

# --- 1. packages -------------------------------------------------------------
pkg_ok <- vapply(boot_pkgs, function(p) requireNamespace(p, quietly = TRUE), logical(1))
cat("\nPACKAGES (figure layer)\n")
if (all(pkg_ok)) {
  cat("  all ", length(boot_pkgs), " present\n", sep = "")
} else {
  cat("  MISSING: ", paste(boot_pkgs[!pkg_ok], collapse = ", "), "\n", sep = "")
  cat("  install with scripts/00_setup_packages.R\n")
}

if (boot_full) {
  cat("\n  myc.bootstrap.full = TRUE -> sourcing 00_setup_packages.R (slow)\n")
  source(here::here("scripts", "00_setup_packages.R"))
}

# --- 2. audit ----------------------------------------------------------------
boot_map$path       <- file.path(root, "results", paste0(boot_map$object, ".rds"))
boot_map$exists     <- file.exists(boot_map$path)
boot_map$bytes      <- ifelse(boot_map$exists, file.size(boot_map$path), NA_real_)
boot_map$mtime      <- as.POSIXct(ifelse(boot_map$exists,
                                         file.mtime(boot_map$path), NA),
                                  origin = "1970-01-01")
prod_path           <- file.path(root, boot_map$producer)
boot_map$prod_mtime <- as.POSIXct(ifelse(file.exists(prod_path),
                                         file.mtime(prod_path), NA),
                                  origin = "1970-01-01")
boot_map$stale      <- !is.na(boot_map$mtime) & !is.na(boot_map$prod_mtime) &
                       boot_map$prod_mtime > boot_map$mtime

fmt_size <- function(b) {
  if (is.na(b)) return("-")
  if (b >= 1024^2) sprintf("%.1f MB", b / 1024^2) else sprintf("%.0f KB", b / 1024)
}

cat("\nRESULTS AUDIT -- the 11 objects the figures/ layer reads\n\n")
cat(sprintf("  %-32s %8s  %-16s %-16s %s\n",
            "object", "size", "written", "producer edited", ""))
for (i in seq_len(nrow(boot_map))) {
  flag <- if (!boot_map$exists[i]) "MISSING" else if (boot_map$stale[i]) "script newer" else ""
  cat(sprintf("  %-32s %8s  %-16s %-16s %s\n",
              boot_map$object[i],
              fmt_size(boot_map$bytes[i]),
              if (boot_map$exists[i]) format(boot_map$mtime[i], "%Y-%m-%d %H:%M") else "-",
              if (is.na(boot_map$prod_mtime[i])) "?" else
                format(boot_map$prod_mtime[i], "%Y-%m-%d %H:%M"),
              flag))
}

if (any(!boot_map$exists)) {
  cat("\n  REBUILD the missing ones by sourcing, in this order:\n")
  for (p in unique(boot_map$producer[!boot_map$exists])) cat("    ", p, "\n", sep = "")
  cat("  (the archive_main_pipeline script uses bare relative paths:",
      "run it with cwd at the repo root)\n")
}
if (any(boot_map$stale)) {
  cat("\n  'script newer' means the producer was EDITED after its output was",
      "written.\n  That is not automatically wrong -- a header or comment edit",
      "moves the mtime\n  too -- but check the diff before trusting the object.\n")
}

# --- 3. load -----------------------------------------------------------------
S <- list()
if (boot_load) {
  have_deseq <- requireNamespace("DESeq2", quietly = TRUE)
  loadable <- boot_map$exists & (!boot_map$s4 | have_deseq)
  if (any(boot_map$s4 & !have_deseq)) {
    cat("\n  DESeq2 not installed: skipping ",
        paste(boot_map$object[boot_map$s4], collapse = ", "),
        " (S4 objects need it to deserialise)\n", sep = "")
  }
  t0 <- Sys.time()
  for (i in which(loadable)) S[[boot_map$object[i]]] <- readRDS(boot_map$path[i])
  cat(sprintf("\nLOADED %d objects into S in %.1f s -- S$mito_content_proxies$shares etc.\n",
              length(S), as.numeric(difftime(Sys.time(), t0, units = "secs"))))
} else {
  cat("\nmyc.bootstrap.load = FALSE -> nothing loaded (audit only)\n")
}

# --- 4. orientation ----------------------------------------------------------
git_say <- function(args) {
  out <- tryCatch(suppressWarnings(
    system2("git", c("-C", shQuote(root), args), stdout = TRUE, stderr = FALSE)),
    error = function(e) character(0))
  if (length(out) == 0) "?" else out
}

cat("\nGIT\n")
cat("  branch : ", git_say(c("rev-parse", "--abbrev-ref", "HEAD"))[1], "\n", sep = "")
# the format string holds a space, and system2() does not quote its arguments --
# unquoted, git reads the subject token as a pathspec and the call fails.
cat("  HEAD   : ", git_say(c("log", "-1", shQuote("--pretty=format:%h %s")))[1],
    "\n", sep = "")
ab <- strsplit(git_say(c("rev-list", "--left-right", "--count",
                         shQuote("HEAD...@{upstream}")))[1], "\t")[[1]]
cat("  origin : ", if (length(ab) == 2) paste0(ab[1], " ahead, ", ab[2], " behind")
                   else "no upstream", "\n", sep = "")
dirty <- git_say(c("status", "--porcelain"))
dirty <- dirty[nzchar(dirty)]
if (length(dirty) == 0 || identical(dirty, "?")) {
  cat("  tree   : clean\n")
} else {
  cat("  tree   : ", length(dirty), " entries\n", sep = "")
  for (d in utils::head(dirty, 10)) cat("           ", d, "\n", sep = "")
}

# Google Drive stamps a zero-byte "Icon\r" into every folder it syncs. Harmless
# in the tree (gitignored), NOT harmless inside .git -- one in refs/ reads as a
# ref and breaks fetch with "fatal: bad object refs/Icon?".
icon_name <- paste0("Icon", "\r")
n_tree <- length(list.files(root, pattern = paste0("^", icon_name, "$"),
                            recursive = TRUE, all.files = TRUE, full.names = TRUE))
n_git  <- length(list.files(file.path(root, ".git"),
                            pattern = paste0("^", icon_name, "$"),
                            recursive = TRUE, all.files = TRUE, full.names = TRUE))
cat("\nGOOGLE DRIVE Icon stubs\n")
cat("  working tree: ", n_tree, "   .git: ", n_git, "\n", sep = "")
if (n_git > 0) {
  cat("  >> STUBS INSIDE .git -- clear them BEFORE any git command:\n")
  cat("     ./paper/clean_icon_files.sh repo\n")
}

cat("\nNEXT: figures/rebuild_manuscript_figures.R rebuilds the four manuscript PDFs.\n")
cat("      quarto render paper/myc_mito.qmd rebuilds the document (writes no PDFs).\n")
hr("="); cat("\n")

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {

  ## what came back
  names(S) |> print()
  utils::str(S$mito_content_proxies, max.level = 1)

  ## the audit as a data frame, if you want to sort or filter it
  boot_map[, c("object", "bytes", "mtime", "producer", "stale")] |> print()

  ## the full stack, when you are doing analysis rather than figures
  options(myc.bootstrap.full = TRUE)
  source(here::here("scripts", "99_session_bootstrap.R"))
}
