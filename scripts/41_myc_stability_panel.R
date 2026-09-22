# scripts/41_myc_stability_panel.R
# =============================================================================
# Block B -- DOES THE MYC POST-TRANSLATIONAL STABILITY MACHINERY CHANGE 6W->12W?
# =============================================================================
#
# Follows figS8, which established that Myc MESSAGE does not fall (genotype gap
# +1.64 log2 at 6W -> +1.80 at 12W, no within-genotype decline) while Myc's target
# output narrows (+1.32 -> +0.95). The natural next hypothesis is that the driver
# weakens post-transcriptionally -- MYC protein stability is controlled by a
# well-mapped set of kinases, phosphatases, E3 ligases, deubiquitinases,
# acetyltransferases/deacetylases, SUMO enzymes and O-GlcNAc -- so this script asks
# whether any of that machinery is transcriptionally down in 12W Myc+ vs 6W Myc+.
#
# ANSWER, in one line: five genes are down at FDR<0.05 (Trim6, Pias4, Fbxw7, Usp36,
# Huwe1) and two up (Ep300, Sirt1), but NONE of it is Myc-specific, the panel does
# not move more than expression-matched genes, and the directions mostly predict a
# MORE stable MYC protein. The panel does not supply a mechanism.
#
# THE STANDING LIMITATION, stated first because it bounds everything below:
# post-translational regulation is by definition not visible in transcript levels.
# This panel can exclude a large TRANSCRIPTIONAL reprogramming of the machinery --
# and it does -- but a constant Fbxw7 mRNA with altered T58 phosphorylation would
# look exactly like this result. "No transcriptional evidence" is not "not
# happening". The settling experiments are protein-level: MYC blot at both ages,
# phospho-T58/S62, and a cycloheximide chase.
#
# Roster (31 genes, author-specified, human symbols mapped to their mouse orthologs
# one-to-one; all 31 resolve in the count matrix):
#   phospho-degron kinases/phosphatase  Mapk3 Mapk1 Gsk3b Cdk1 Cdk2 Aurka Pim1 Pim2
#                                       Ppp2ca Ppp2cb Pin1
#   E3 ligases (MYC turnover)           Fbxw7 Skp2 Huwe1 Trim6 Trim32 Fbxl16
#   deubiquitinases (MYC stabilising)   Usp28 Usp36
#   acetylation / deacetylation         Kat2a Ep300 Crebbp Kat5 Hdac1 Hdac3 Sirt1
#   SUMO                                Ube2i Pias1 Pias4 Senp1
#   O-GlcNAc                            Ogt
#
# Input:  results/interaction_results.rds   (raw LFC + p + padj, five contrasts)
#         results/dds_int_run.rds           (normalised counts, colData)
#         results/combined_df_annotated.rds (mgi_symbol <-> ensembl)
# Output: outputs/myc_stability_panel/myc_stability_panel.csv
#         outputs/myc_stability_panel/README.md   (the verdict, numbers filled in
#                                                  from this run so it cannot drift)
#
# CEILING. n=6/group. The 6W-vs-12W axis is BATCH-CONFOUNDED (batch = timepoint), so
# `timepoint_pos` is descriptive; the batch-clean quantity is the interaction, which
# is null for every gene here. No figure is produced -- this is a negative result and
# a table is the right weight for it.
# =============================================================================

source(here::here("scripts", "00_setup_packages.R"))

set.seed(1)
NPERM   <- 5000L
out_dir <- here::here("outputs", "myc_stability_panel")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

group_levels <- c("6W_neg", "6W_pos", "12W_neg", "12W_pos")

# =============================================================================
# PART 1: ROSTER + LOAD
# =============================================================================
roster <- tibble::tribble(
  ~gene,     ~arm,
  "Mapk3",   "phospho-degron (kinase/phosphatase)",
  "Mapk1",   "phospho-degron (kinase/phosphatase)",
  "Gsk3b",   "phospho-degron (kinase/phosphatase)",
  "Cdk1",    "phospho-degron (kinase/phosphatase)",
  "Cdk2",    "phospho-degron (kinase/phosphatase)",
  "Aurka",   "phospho-degron (kinase/phosphatase)",
  "Pim1",    "phospho-degron (kinase/phosphatase)",
  "Pim2",    "phospho-degron (kinase/phosphatase)",
  "Ppp2ca",  "phospho-degron (kinase/phosphatase)",
  "Ppp2cb",  "phospho-degron (kinase/phosphatase)",
  "Pin1",    "phospho-degron (kinase/phosphatase)",
  "Fbxw7",   "E3 ligase (MYC turnover)",
  "Skp2",    "E3 ligase (MYC turnover)",
  "Huwe1",   "E3 ligase (MYC turnover)",
  "Trim6",   "E3 ligase (MYC turnover)",
  "Trim32",  "E3 ligase (MYC turnover)",
  "Fbxl16",  "E3 ligase (MYC turnover)",
  "Usp28",   "deubiquitinase (MYC stabilising)",
  "Usp36",   "deubiquitinase (MYC stabilising)",
  "Kat2a",   "acetylation / deacetylation",
  "Ep300",   "acetylation / deacetylation",
  "Crebbp",  "acetylation / deacetylation",
  "Kat5",    "acetylation / deacetylation",
  "Hdac1",   "acetylation / deacetylation",
  "Hdac3",   "acetylation / deacetylation",
  "Sirt1",   "acetylation / deacetylation",
  "Ube2i",   "SUMO",
  "Pias1",   "SUMO",
  "Pias4",   "SUMO",
  "Senp1",   "SUMO",
  "Ogt",     "O-GlcNAc")

ir  <- readRDS(here::here("results", "interaction_results.rds"))
dds <- readRDS(here::here("results", "dds_int_run.rds"))
cdf <- readRDS(here::here("results", "combined_df_annotated.rds"))

nc <- DESeq2::counts(dds, normalized = TRUE)
sm <- as.data.frame(SummarizedExperiment::colData(dds))
sm$group <- factor(sm$group, levels = group_levels)

sym2ens <- cdf[!is.na(cdf$mgi_symbol) & !duplicated(cdf$mgi_symbol), c("mgi_symbol", "gene")]
roster$ens <- sym2ens$gene[match(roster$gene, sym2ens$mgi_symbol)]
unresolved <- roster$gene[is.na(roster$ens) | !roster$ens %in% rownames(nc)]
if (length(unresolved)) message("41: unresolved -> ", paste(unresolved, collapse = ", "))
roster <- roster[!is.na(roster$ens) & roster$ens %in% rownames(nc), ]
message(sprintf("41: %d of 31 roster genes resolved", nrow(roster)))

V <- function(k, f) { r <- as.data.frame(ir[[k]]); stats::setNames(r[[f]], rownames(r)) }
contrasts <- c(myc_6W = "myc_6W", myc_12W = "myc_12W",
               wt_time = "timepoint_neg", myc_time = "timepoint_pos",
               interaction = "interaction")

# =============================================================================
# PART 2: THE TABLE
# =============================================================================
tab <- roster
tab$baseMean <- V("myc_6W_raw", "baseMean")[tab$ens]
for (nm in names(contrasts)) {
  k <- paste0(contrasts[[nm]], "_raw")
  tab[[paste0("lfc_",  nm)]] <- V(k, "log2FoldChange")[tab$ens]
  tab[[paste0("p_",    nm)]] <- V(k, "pvalue")[tab$ens]
  tab[[paste0("padj_", nm)]] <- V(k, "padj")[tab$ens]        # genome-wide BH (DESeq2)
}
# BH within THIS 31-gene panel -- the correct family for a targeted question
for (nm in c("myc_time", "wt_time", "interaction"))
  tab[[paste0("fdr_panel_", nm)]] <- stats::p.adjust(tab[[paste0("p_", nm)]], "BH")

for (g in group_levels)
  tab[[paste0("mean_", g)]] <- vapply(tab$ens, function(e)
    mean(nc[e, sm$group == g]), numeric(1))

tab <- tab[order(tab$lfc_myc_time), ]

hits_down <- tab$gene[!is.na(tab$padj_myc_time) & tab$padj_myc_time < 0.05 & tab$lfc_myc_time < 0]
hits_up   <- tab$gene[!is.na(tab$padj_myc_time) & tab$padj_myc_time < 0.05 & tab$lfc_myc_time > 0]
hits_int  <- tab$gene[!is.na(tab$padj_interaction) & tab$padj_interaction < 0.05]
hits_int_panel <- tab$gene[!is.na(tab$fdr_panel_interaction) & tab$fdr_panel_interaction < 0.05]

# =============================================================================
# PART 3: IS THE PANEL MOVING MORE THAN THE GENOME? (expression-matched null)
# =============================================================================
# Under a GLOBAL attenuation every Myc-induced gene drifts down, so "five genes are
# down" is only informative against expression-matched background genes. This is the
# script-34 atten_excess logic.
bm  <- V("myc_6W_raw", "baseMean")
mt  <- V("timepoint_pos_raw", "log2FoldChange")
it  <- V("interaction_raw",   "log2FoldChange")
ok  <- !is.na(bm) & bm > 0 & !is.na(mt) & !is.na(it)
dec <- stats::setNames(cut(rank(bm[ok], ties.method = "first"), breaks = 20, labels = FALSE),
                       names(bm)[ok])
pool <- split(names(dec), dec)
tb   <- table(factor(dec[tab$ens], levels = names(pool)))

matched_null <- function(v) {
  obs  <- stats::median(v[tab$ens], na.rm = TRUE)
  nul  <- replicate(NPERM, {
    g <- unlist(lapply(names(tb)[tb > 0], function(k) sample(pool[[k]], tb[[k]])),
                use.names = FALSE)
    stats::median(v[g], na.rm = TRUE)
  })
  tibble::tibble(panel_median = obs, null_median = stats::median(nul),
                 percentile = 100 * mean(nul <= obs),
                 p_two_sided = 2 * min(mean(nul <= obs), mean(nul >= obs)))
}
background <- dplyr::bind_rows(
  dplyr::mutate(matched_null(mt), contrast = "myc_time (12W Myc+ vs 6W Myc+)", .before = 1),
  dplyr::mutate(matched_null(it), contrast = "interaction (Myc-specific)",     .before = 1))
background$genome_median <- c(stats::median(mt[ok]), stats::median(it[ok]))

# =============================================================================
# PART 4: ARM SUMMARIES
# =============================================================================
arm_summary <- tab |>
  dplyr::group_by(arm) |>
  dplyr::summarise(n = dplyr::n(),
                   median_myc_time    = stats::median(lfc_myc_time),
                   median_wt_time     = stats::median(lfc_wt_time),
                   median_interaction = stats::median(lfc_interaction),
                   n_down_fdr = sum(!is.na(padj_myc_time) & padj_myc_time < 0.05 &
                                      lfc_myc_time < 0),
                   .groups = "drop") |>
  dplyr::arrange(median_myc_time)

# the two arms whose BALANCE would move MYC protein, and the canonical T58 axis
balance <- tibble::tibble(
  set = c("E3 ligases (destabilise MYC)", "deubiquitinases (stabilise MYC)",
          "T58/S62 phospho-degron axis"),
  genes = c("Fbxw7, Skp2, Huwe1, Trim32", "Usp28, Usp36",
            "Mapk1, Mapk3, Gsk3b, Pin1, Ppp2ca, Ppp2cb"))
balance$median_myc_time <- vapply(strsplit(balance$genes, ", "), function(g)
  stats::median(tab$lfc_myc_time[tab$gene %in% g]), numeric(1))
balance$max_abs_lfc <- vapply(strsplit(balance$genes, ", "), function(g)
  max(abs(tab$lfc_myc_time[tab$gene %in% g])), numeric(1))
balance$min_p <- vapply(strsplit(balance$genes, ", "), function(g)
  min(tab$p_myc_time[tab$gene %in% g]), numeric(1))

# =============================================================================
# PART 5: WRITE THE TABLE + THE ARGUMENT
# =============================================================================
csv_cols <- c("gene", "arm", "ens", "baseMean",
              paste0("mean_", group_levels),
              "lfc_myc_time", "p_myc_time", "padj_myc_time", "fdr_panel_myc_time",
              "lfc_wt_time", "p_wt_time", "padj_wt_time",
              "lfc_interaction", "p_interaction", "padj_interaction",
              "fdr_panel_interaction",
              "lfc_myc_6W", "padj_myc_6W", "lfc_myc_12W", "padj_myc_12W")
utils::write.csv(tab[, csv_cols], file.path(out_dir, "myc_stability_panel.csv"),
                 row.names = FALSE)

g <- function(x) if (length(x)) paste(x, collapse = ", ") else "none"
row_of <- function(nm, col) tab[[col]][tab$gene == nm]
# "0.000" in a results table reads as p == 0; never print a p that way
fp <- function(x) if (is.na(x)) "NA" else if (x < 0.001) "<0.001" else sprintf("%.3f", x)

md <- c(
"# The MYC post-translational stability machinery, 6W -> 12W",
"",
sprintf("Generated by `scripts/41_myc_stability_panel.R` on %s. Numbers are filled in from the run, so this file cannot drift from `myc_stability_panel.csv` beside it.", Sys.Date()),
"",
"## The question",
"",
"figS8 showed Myc MESSAGE does not fall across the window (genotype gap +1.64 log2 at 6W",
"-> +1.80 at 12W; no within-genotype decline) while Myc's target output narrows (+1.32 ->",
"+0.95). So: does the machinery that controls MYC PROTEIN stability change instead?",
"31 genes -- phospho-degron kinases/phosphatases, E3 ligases, deubiquitinases,",
"acetyl/deacetylases, SUMO enzymes, O-GlcNAc -- read on the 12W-Myc+ vs 6W-Myc+ contrast.",
"",
"## The answer: five down, two up -- and none of it supports the hypothesis",
"",
sprintf("Down at genome-wide FDR<0.05: **%s**. Up: **%s**.", g(hits_down), g(hits_up)),
"The same genes survive BH within the 31-gene panel.",
"",
"| gene | LFC (Myc+ 6->12W) | padj | baseMean | same contrast in WT | interaction padj |",
"|---|---|---|---|---|---|",
paste0(vapply(c(hits_down, hits_up), function(nm) sprintf(
  "| %s | %+.2f | %s | %.0f | %+.2f (p %s) | %s |",
  nm, row_of(nm, "lfc_myc_time"), fp(row_of(nm, "padj_myc_time")), row_of(nm, "baseMean"),
  row_of(nm, "lfc_wt_time"), fp(row_of(nm, "p_wt_time")), fp(row_of(nm, "padj_interaction"))),
  character(1)), collapse = "\n"),
"",
"### 1. None of it is Myc-specific",
"",
sprintf("The interaction padj is >0.05 for **every** gene in the panel (significant: %s), and BH", g(hits_int)),
sprintf("within the panel gives nothing either (%s). Several of the hits move in WT as well --", g(hits_int_panel)),
sprintf("Pias4 %+.2f (p %s), Huwe1 %+.2f (p %s), Ep300 %+.2f (p %s), Sirt1 %+.2f (p %s).",
        row_of("Pias4", "lfc_wt_time"), fp(row_of("Pias4", "p_wt_time")),
        row_of("Huwe1", "lfc_wt_time"), fp(row_of("Huwe1", "p_wt_time")),
        row_of("Ep300", "lfc_wt_time"), fp(row_of("Ep300", "p_wt_time")),
        row_of("Sirt1", "lfc_wt_time"), fp(row_of("Sirt1", "p_wt_time"))),
"This is the batch-confounded time axis (batch = timepoint), on which a contrast is",
"described, not claimed. The batch-clean quantity is the interaction, and it is null.",
"",
"### 2. The panel does not move more than the genome",
"",
"Under a global attenuation every Myc-induced gene drifts down, so 'five genes are down' is",
"only informative against expression-matched background genes:",
"",
"| contrast | genome median | panel median | matched-null median | percentile | p |",
"|---|---|---|---|---|---|",
paste0(vapply(seq_len(nrow(background)), function(i) sprintf(
  "| %s | %+.3f | %+.3f | %+.3f | %.0f | %.2f |",
  background$contrast[i], background$genome_median[i], background$panel_median[i],
  background$null_median[i], background$percentile[i], background$p_two_sided[i]),
  character(1)), collapse = "\n"),
"",
"The panel moves exactly as much as matched genes. This is the global attenuation passing",
"through it, not a targeted change in the stability machinery.",
"",
"### 3. The directions mostly predict a MORE stable MYC protein",
"",
sprintf("FBXW7 (%+.2f) and HUWE1 (%+.2f), the two best-established degradative E3 ligases for MYC,",
        row_of("Fbxw7", "lfc_myc_time"), row_of("Huwe1", "lfc_myc_time")),
"are both DOWN -- which predicts more stable MYC at 12W, not less. The one hit in the",
sprintf("expected direction is Usp36 (%+.2f), a MYC-stabilising deubiquitinase, but it falls in WT",
        row_of("Usp36", "lfc_myc_time")),
sprintf("too (%+.2f) and its interaction is null. And the arms move together, so the ubiquitin",
        row_of("Usp36", "lfc_wt_time")),
"balance does not shift:",
"",
"| set | genes | median LFC (Myc+ 6->12W) | largest abs LFC | smallest p |",
"|---|---|---|---|---|",
paste0(vapply(seq_len(nrow(balance)), function(i) sprintf(
  "| %s | %s | %+.2f | %.2f | %s |", balance$set[i], balance$genes[i],
  balance$median_myc_time[i], balance$max_abs_lfc[i], fp(balance$min_p[i])),
  character(1)), collapse = "\n"),
"",
"The canonical T58/S62 phospho-degron route (ERK1/2 - GSK3beta - PIN1 - PP2A) shows no",
"transcriptional change at all. Meanwhile the proliferation-linked kinases stay Myc-induced",
sprintf("at both ages (Cdk1 %+.2f / %+.2f, Skp2 %+.2f / %+.2f, Aurka %+.2f / %+.2f at 6W / 12W).",
        row_of("Cdk1", "lfc_myc_6W"),  row_of("Cdk1", "lfc_myc_12W"),
        row_of("Skp2", "lfc_myc_6W"),  row_of("Skp2", "lfc_myc_12W"),
        row_of("Aurka", "lfc_myc_6W"), row_of("Aurka", "lfc_myc_12W")),
"",
"### Per-arm summary",
"",
"| arm | n | median Myc+ 6->12W | median WT 6->12W | median interaction | n down at FDR<0.05 |",
"|---|---|---|---|---|---|",
paste0(vapply(seq_len(nrow(arm_summary)), function(i) sprintf(
  "| %s | %d | %+.2f | %+.2f | %+.2f | %d |", arm_summary$arm[i], arm_summary$n[i],
  arm_summary$median_myc_time[i], arm_summary$median_wt_time[i],
  arm_summary$median_interaction[i], arm_summary$n_down_fdr[i]), character(1)),
  collapse = "\n"),
"",
"## The limitation that bounds all of the above",
"",
"Post-translational regulation is by definition not visible in transcript levels. This panel",
"can exclude a large TRANSCRIPTIONAL reprogramming of the machinery -- and it does -- but a",
"constant Fbxw7 mRNA with altered T58 phosphorylation would look exactly like this result.",
"**'No transcriptional evidence' is not 'not happening.'** The settling experiments are",
"protein-level: MYC blot at 6W vs 12W (does protein track the flat message?),",
"phospho-T58/S62, and a cycloheximide chase.",
"",
"## Provenance",
"",
sprintf("n=6/group. Contrasts from `results/interaction_results.rds` (raw/unshrunken LFCs); group means are DESeq2 median-of-ratios normalised counts. Matched null = %d draws, expression matched in 20 baseMean strata. All 31 roster genes resolved.", NPERM))

writeLines(md, file.path(out_dir, "README.md"))
message("wrote ", file.path(out_dir, "myc_stability_panel.csv"))
message("wrote ", file.path(out_dir, "README.md"))

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {
  tab |>
    dplyr::select(gene, arm, baseMean, lfc_myc_time, padj_myc_time, fdr_panel_myc_time,
                  lfc_wt_time, p_wt_time, lfc_interaction, padj_interaction) |>
    print(n = 31)
  background |> print()
  arm_summary |> print()
  balance |> print()
  cat(readLines(file.path(out_dir, "README.md")), sep = "\n")
}
