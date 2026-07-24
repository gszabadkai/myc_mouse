# =============================================================================
# reconcile_gene_symbols.R -- gene-symbol vintage reconciliation
# -----------------------------------------------------------------------------
# ROOT CAUSE it fixes. The library GMTs and MitoCarta 3.0 (Sheet 4) carry the
# ORIGINAL MitoCarta gene symbols -- e.g. the ATP-synthase names Atp5a1, Atp5b,
# Atp5c1, Atp5o, Atpif1 -- while this project's expression annotation (the count
# matrix, combined_df$mgi_symbol and ortholog_table$external_gene_name) uses the
# CURRENT symbols (Atp5f1a, Atp5f1b, Atp5f1c, Atp5po, Atp5if1). A plain symbol
# match therefore SILENTLY DROPS the renamed genes: 17 of Complex V's 24 genes,
# ~12% of nuclear OXPHOS. The dropped genes behave like their set (indeed, the
# ATP-synthase genes are among the most strongly Myc-induced), so nothing reverses
# -- but set-level shares / NES / GSVA are conservatively UNDER-counted. This
# helper reconciles a set's symbols (any vintage) to the identifiers the data use.
#
# Two exported resolvers, both idempotent on already-current names:
#   recon_to_ensembl(symbols, universe_ensembl)
#       -> Ensembl IDs present in universe_ensembl. Use for Ensembl-keyed data
#          (count matrix, DESeqResults rownames). Three routes, union:
#          (1) current-symbol -> Ensembl via ortholog_table,
#          (2) MitoCarta Sheet 2 Symbol -> EnsemblGeneID (authoritative for old
#              MitoCarta names),
#          (3) org.Mm.eg.db ALIAS -> Ensembl (general, catches non-MitoCarta renames).
#   recon_to_current(symbols, universe_symbols)
#       -> CURRENT symbols present in universe_symbols. Use for symbol-keyed data
#          (e.g. mitoPPS expr_symbols, fGSEA ranks named by external_gene_name).
#
# Maps are built once and cached in a private environment. Needs: readxl,
# org.Mm.eg.db, AnnotationDbi, and results/ortholog_table.rds + Mouse.MitoCarta3.0.xls.
# =============================================================================

.recon_env <- new.env(parent = emptyenv())

.recon_build <- function() {
  if (!is.null(.recon_env$cur2ens)) return(invisible())
  for (p in c("readxl", "org.Mm.eg.db", "AnnotationDbi")) {
    if (!requireNamespace(p, quietly = TRUE)) {
      stop("reconcile_gene_symbols.R needs the '", p, "' package")
    }
  }
  ot <- readRDS(here::here("results", "ortholog_table.rds"))
  .recon_env$ens2cur <- stats::setNames(ot$external_gene_name, ot$ensembl_gene_id)
  c2e <- ot[!is.na(ot$external_gene_name) & !is.na(ot$ensembl_gene_id), ]
  .recon_env$cur2ens <- stats::setNames(c2e$ensembl_gene_id, c2e$external_gene_name)
  s2 <- as.data.frame(readxl::read_excel(
    here::here("data", "Mouse.MitoCarta3.0.xls"), sheet = 2))
  s2 <- s2[!is.na(s2$Symbol) & !is.na(s2$EnsemblGeneID) & s2$EnsemblGeneID != "", ]
  .recon_env$mc <- do.call(rbind, lapply(seq_len(nrow(s2)), function(i)
    data.frame(sym = s2$Symbol[i],
               ens = trimws(strsplit(s2$EnsemblGeneID[i], "[|]")[[1]]),
               stringsAsFactors = FALSE)))
  invisible()
}

recon_to_ensembl <- function(symbols, universe_ensembl) {
  .recon_build()
  symbols <- unique(symbols[!is.na(symbols)])
  out <- character(0)
  # Route 1 -- the symbol is already a CURRENT symbol (authoritative; one Ensembl each,
  # so an alias collision with an unrelated gene cannot occur for resolved symbols).
  hit1 <- symbols[symbols %in% names(.recon_env$cur2ens)]
  out  <- c(out, unname(.recon_env$cur2ens[hit1]))
  todo <- setdiff(symbols, hit1)
  # Route 2 -- unresolved names via MitoCarta's OWN Ensembl IDs (Sheet 2).
  if (length(todo)) {
    m2   <- .recon_env$mc[.recon_env$mc$sym %in% todo, ]
    out  <- c(out, m2$ens)
    todo <- setdiff(todo, m2$sym)
  }
  # Route 3 -- still unresolved via org.Mm.eg.db aliases (general renames).
  if (length(todo)) {
    al  <- suppressMessages(AnnotationDbi::mapIds(
      org.Mm.eg.db::org.Mm.eg.db, keys = todo, column = "ENSEMBL",
      keytype = "ALIAS", multiVals = "first"))
    out <- c(out, unname(al[!is.na(al)]))
  }
  intersect(unique(out), universe_ensembl)
}

recon_to_current <- function(symbols, universe_symbols) {
  .recon_build()
  ens <- recon_to_ensembl(symbols, names(.recon_env$ens2cur))
  cur <- unique(c(intersect(symbols, universe_symbols), unname(.recon_env$ens2cur[ens])))
  intersect(cur[!is.na(cur)], universe_symbols)
}
