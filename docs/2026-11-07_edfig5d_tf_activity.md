---
date: 2026-11-07
tags: [project/myc_mouse, handoff, nat-metab, extended-data-5, ed-fig5d, orthotopic, tf-activity, decoupler, collectri, plan]
status: >
  DECIDED AND BUILT (evening update below). The morning's plan and findings follow it
  unchanged as the record. Script 58 is written and dry-run; the author runs it.
  Nothing committed.
relates-to:
  - scripts/external/Pgc bclxl analysis.R      (the copy that made ED 5d; untracked, the author's)
  - scripts/external/Pgc bclxl analysis 1.R    (a later edit of it; did not make the figure; untracked)
  - data/orthotopic_series/README.md            (the arms, the 2026-09-09 identity correction, the token trap)
  - docs/2026-09-09_orthotopic_identity_correction.md
  - docs/2026-10-28_section3_results_digest.md  (script 57: no DESeq2 Wald fit anywhere in this cohort)
  - figures/natmetab/FIGURE_RULES.md            (the construct colours; the Prism format)
  - ../myc_human_validation/data/collectri_human/README.md (the CollecTRI snapshot precedent; read-only)
---

# ED Fig. 5d, orthotopic TF activity: findings and plan, decisions pending

## Update, evening of 2026-11-07: the author's decisions and what was built

- **Analysis (b)**, without reproducing the original: 21 vector-series tumours, construct
  transcripts left out, the pinned network with resolved names (`data/collectri_mouse/`).
  `scripts/58_orthotopic_tf_activity.R` is written and dry-run; the author runs it.
- **The TF charts are dropped** (author: they do not test the bioenergetics-independent hypothesis;
  the mRNA does). TF activity stays in script 58 as the record of why the original charts went.
- **The mRNA hypothesis** (author: "we need a justification why are we looking for these genes"):
  35 genes fixed from the literature (Consensus and PubMed) before their read-out, except *Twist1*,
  *Cdh1*, *Mmp2* and *Mmp9*, already seen: Id proteins (Luo 2016, Oh 2021), the Tcf4-Twist1
  programme (Oh 2021, Fang 2022), MMPs (ED 5e), Luo 2016's PGC-1α-repressed melanoma signature (22 of
  23; the external script's own gene list) and Snail (Nam 2020). DOIs in `figures/natmetab/LEGENDS.md`.
- **ED 5d = heatmap H1** (author's choice over a compact heatmap and a forest plot,
  `figures/natmetab/proposals/p5d_expression_options.R`), with the author's changes: MMPs above the
  signature; the signature named "Adhesion, TGF-β and Wnt signature" and judged only as a set (no
  gene arrows, no bold; its fgsea test in the legend; 2026-11-08); bold names where a gene judged
  alone changes in the predicted direction; a star where Bcl-xL + Pgc-1α differs from at least one
  control; ±1.5. `figures/natmetab/edfig5d_invasion_genes.R`, 52.5 x 86 mm.
- **Dry-run reading** (the author's run is the record): supported -- *Id1* up, *Cdh1* up, *Cdh2*,
  *Tcf4*, *Twist1*, *Mmp2* down, *Mmp9* down against Bcl-xL + EV only; not -- *Id2*, *Vim*, *Pthlh*,
  *Snai1/2*, and the melanoma signature (fgsea NES +0.90 and +0.84, *P* 0.61 and 0.70; seven of its
  genes rise significantly). 25 of 35 genes carry a star; 7 of the 12 judged gene by gene are bold. Stromal markers are lower in
  the Pgc-1α arm (median log2 fold change -0.48), endothelial and adipose markers higher.
- Colour proposals for the TF charts (`p5d_tf_colours.R`) are moot.

## 0. Where things stand

- **The request (author, 2026-11-07).** ED 5d's two panels are the only ones in the paper not
  finalised. They come from the orthotopic analysis in `scripts/external/`. Read those scripts,
  reproduce the panels in the paper's colour code, and show the plan before writing the script.
- **What happened.** The plan was shown. Before writing it, Claude re-ran the original recipe in
  the session scratchpad. The re-run reproduces the panels, but the Twist1 bar, the one the text
  leans on, is carried partly by the PGC-1a construct's own transcript (section 4.4). The
  author: "great that the error is spotted"; they return with the decisions in section 5.
- **Nothing is built.** No numbered script, no panel script, no `results/` object, no PDF.
  The only repo changes of the session are this note and the dated paragraph at the top of
  `docs/handoff.md`. Nothing is committed.
- **Loose ends, closed by the author (2026-11-07):**
  - **The old `EDFig5` is superseded by the current structure.** That is the Bbc3-route ED of
    2026-10-28 (`figures/natmetab/edfig5a_bbc3_routes.R`, `edfig5b_ox_rel_routes.R`, and their
    PDFs and Prism tables in `outputs/natmetab/EDFig5/`). **The next session compares the old
    and current versions of the figures, which were finalised in the human_exploratory repo.
    Do nothing on them until then.**
  - **The Doc cites "ED Fig. 5a-f"; the assembled figure has a-e.** The author corrects the
    citation.

## 1. What ED 5d is

The Drive figure folder (`1qLbI6nQlLczmUBrmWyMvMmNPJoG6Fano`, file `NM_Extended_Data_Fig5.pdf`,
modified 2026-11-06) holds ED 5 as assembled: **a**, spontaneous lung metastatic burden (Prism;
EV + EV, Bcl-xL + EV, Bcl-xL + Pgc-1a); **b**, representative lungs; **c**, circulating tumour
cells; **d**, two TF-activity bar charts placed as PNGs, "Transcription factor activity (PGC vs
BCL-XL)" and "(PGC vs EV)", 25 TFs each, red-white-blue fill by score; **e**, westerns of Twist1,
MMP9 and MMP2.

The text it supports (Doc, "Nat metab version", the paragraph after the Fig. 4 legend):
"Notably, spontaneous metastatic dissemination from orthotopic lesions was attenuated, due to
diminished Twist, its transcriptional targets Mmp2 and Mmp9, and a reduced amount of circulating
tumour cells (ED Fig. 5a-f)." The Doc has an "Extended Data Figure 5" heading with no legend, and
no Methods paragraph for the TF analysis.

## 2. The original analysis

`scripts/external/Pgc bclxl analysis.R` made the figure: its TF section hard-codes the titles
"PGC_vs_BCL-XL" and "PGC_vs_EV" and the axis "Activity score (PGC − BCL-XL)" that the PNGs carry.
`Pgc bclxl analysis 1.R` is a later edit (working directory `Bulkseq_New_Amy`, BoMY and BoMYBFGE
merged into one group, two more contrasts, the TF loop driven by the contrast list, so its labels
read "Pgc1a vs BclxL"). It did not make the figure.

1. **Counts.** `salmon.merged.gene_counts.tsv`, the file pinned at
   `data/orthotopic_series/` (same name; the re-run reproduces the panels from the pinned copy).
2. **Groups.** Strips the plate prefix and the `BRNS...`/`_...` suffixes and keeps the leading
   letters, so every sample lands in its arm. The old copy keeps `BoMYBFGE` as its own group:
   11 groups, 67 samples (BclxL 7, BoMY1 3, BoMYBFGE 3, EV 7, FaMY1 6, KOEV 7, KOPgc1a 7, MYAZ 6,
   NTEV 7, NTPgc1a 7, Pgc1a 7).
3. **One DESeq2 model** over all 67, `~ condition`, rounded counts, genes with at least 10 counts in
   at least 2 samples (19,538 genes). MYAZ, FaMY, BoMY and the CRISPR arms share the dispersion fit.
4. **Contrasts** `Pgc1a` vs `BclxL` and `Pgc1a` vs `EV`, `results()` defaults. By the 2026-09-09
   identity correction `Pgc1a` is **Bcl-xL + PGC-1a**, so the panels compare Bcl-xL + Pgc-1a with
   Bcl-xL + EV (a clean PGC-1a contrast on a shared Bcl-xL pool) and with EV + EV (two
   manipulations). The old titles' "PGC" is the combined arm.
5. **Gene table.** `na.omit()` drops genes with no adjusted P (independent filtering, Cook's);
   symbols, Entrez IDs and names from `org.Mm.eg.db`, rows lacking any dropped; per symbol, the
   Wald `stat` with the largest |stat|.
6. **TF activity.** `get_collectri(organism = "mouse", split_complexes = FALSE)` (fetched live from
   OmniPath), `run_ulm(minsize = 5)`: per TF, the t-value of the regression of the Wald stats on
   its signed targets (`mor` +1/-1, 0 for non-targets).
7. **Plot.** The top 25 by `rnk = ifelse(score > 0, rank(-score), rank(-abs(score)))`. `ifelse()`
   ranks over all TFs, so repressed TFs are ranked by |score| among every TF, not among repressed
   ones as in decoupleR's vignette. The set it picks is close to, but not always the same as, the
   25 largest |score| (same set in 2 of the 8 re-runs, one to three TFs different in the other
   six). Fill:
   `scale_fill_gradient2` RdBu by score.
8. **Not on the figure:** fGSEA on MSigDB H/C2/C3/C5/C6 ranked by log2FC; MMP and integrin
   log2FC plots; a DoRothEA A-C `run_wmean` on log2FC whose plot is titled "PGC vs EV" but reads
   the `Pgc_vs_BclxL` table; a TF-to-MMP link table.

## 3. The scratch re-run (scratchpad only)

- **Run 1.** DESeq2 as the original (all 67, the old 11-group parse) and on the 21 vector-series
  samples alone (17,973 genes after the same filter; 15,459 genes into ULM, against 16,274 for the
  67). `get_collectri()` failed (section 4.2), so OmniPath's static mouse table was read from the
  OmnipathR cache and passed through `get_collectri()`'s own post-processing.
- **Run 2.** The mouse CollecTRI queried from the OmniPath web service directly, processed the same
  way (appendix B), ULM for each fit and contrast, once as the original and once with the two
  construct transcripts (`Ppargc1a`, `Bcl2l1`) left out of the gene statistics.
- Both runs compared against the old panels' scores, read off the figure (appendix A).
- Versions: R 4.6.1, DESeq2 1.52.0, decoupleR 2.17.0, OmnipathR 4.0.0, org.Mm.eg.db 3.23.0,
  dorothea 1.23.0, IHW 1.40.0 (installed, not yet used).

## 4. Findings

### 4.1 The panels reproduce

With the web-service network and the 67-sample fit, as the original:

| panel | TFs recovered | r with the old scores | differences |
|---|---|---|---|
| vs Bcl-xL + EV | 24 of 25 | 0.996 | Snai1 enters, Maz (A0A0U1RNL9) leaves; Twist1 -2.61 against -3.2 |
| vs EV + EV | 23 of 25 | 0.953 | Zbtb17 and Arid4a enter, Zfpm1 and Foxa3 leave; Foxa3 flips sign (-2.8 to +1.78); Twist1 -2.56 against -3.1 |

The 2023 static table recovers only 17 and 16 of 25 and has no Esrra at all, so the figure was
made with a live network close to today's (it shows the same accession names).

### 4.2 OmnipathR 4.0.0 cannot fetch CollecTRI on this machine

`OmnipathR::collectri()`, and so `decoupleR::get_collectri()`, fail for every organism, human
included: "Join columns in `y` must be present in the data. Problem with `ncbi_tax_id`". It
fails while building its organism table, which it parses from Ensembl's species page; the cached
copy of that page (2026-08-28) holds no species table (two `<tr>`), and the human repo found the
failure survives a cache wipe. The fallback, `static_table()`, then fails in `unnest_evidences()`. The human repo met the same bug on 2026-08-28 and snapshotted CollecTRI by
`curl` instead, with a README and a checksum (`../myc_human_validation/data/collectri_human/`).

The web service answers directly (HTTP 200 in 0.07 s):
`https://omnipathdb.org/interactions?resources=CollecTRI&genesymbols=yes&fields=sources,references,curation_effort&format=tsv&organisms=10090`
returned 39,849 rows on 2026-11-07, 39,521 edges and 1,106 sources after the post-processing. The
static fallback, `http://no-tls.static.omnipathdb.org/resources/interactions_collectri_10090.tsv.gz`,
is dated 2023-09-17 (38,389 rows).

### 4.3 The unreadable names, and the edges they cost

OmniPath's mouse table names 226 sources and 1,003 targets by UniProt accession, because no symbol
was found for them. **6,386 edges have an accession-named source** (kept, but the TF name is
unreadable). **5,870 have an accession-named target, which can never match a gene, so ULM drops
them without a message.** The six on the panels or in the re-run's top 25 (UniProt REST,
2026-11-07):

| on the panel | gene | targets |
|---|---|---|
| A0A979HLR9 | Sp1 | 1,339 |
| A0A1Y7VJR5 | Sp4 | 45 |
| A0A0U1RNL9 | Maz | 52 |
| A3KN90 | Nr3c2 | 54 |
| A0A0R4J132 | Mybl1 | |
| A0A0N4SV71 | Creb5 | |

None of the six symbols also exists as a symbol-named source, so for sources the fix is a relabel;
for targets it returns edges and will move scores. `org.Mm.eg.db` resolves 857 of the 1,077
accessions (not A0A979HLR9); 11 map to more than one symbol. The rest need UniProt.

### 4.4 The Twist1 bar depends partly on the transgene

CollecTRI lists **`Ppargc1a` as a target Twist1 represses** (`mor` -1). The construct raises
`Ppargc1a` by 10.0 log2 against Bcl-xL + EV and 11.8 against EV + EV (Wald stat +19 to +23), so
the construct's own transcript reads as low Twist1 activity. In the 21-sample contrast against
Bcl-xL + EV (run 1, 2023 network), `Ppargc1a` contributes -22.7 (stat x mor) to Twist1, the largest
of its 65 targets present; next come `Mmp2` -12.9, `Cdh1` -10.7 (a repressed target, up 2.9 log2)
and `Cdh5` -7.5. Fifty TFs have `Ppargc1a` as a target and sixty have `Bcl2l1`.

Neither construct transcript is a regulated target in these arms (`Bcl2l1` is uninterpretable in
the 14 construct-carrying samples; data/orthotopic_series/README.md). Leaving both out of the gene
statistics:

| fit | contrast | variant | Twist1: score (P; BH; rank) | Snai1: score (P; BH; rank) | TFs at BH < 0.05 |
|---|---|---|---|---|---|
| 67 | vs Bcl-xL + EV | as the original | -2.61 (0.0091; 0.29; 18) | -2.65 (0.0081; 0.29; 16) | 4, all up |
| 67 | vs Bcl-xL + EV | constructs out | -1.98 (0.047; 0.61; 45) | -2.65 (0.008; 0.30; 15) | Esrra, Ppara, Nrf1, Srebf1 |
| 67 | vs EV + EV | as the original | -2.56 (0.010; 0.33; 19) | -2.21 (0.027; 0.45; 36) | 4, all up |
| 67 | vs EV + EV | constructs out | -1.95 (0.051; 0.64; 47) | -2.21 (0.027; 0.50; 32) | Esrra, Ppara, Nrf1 |
| 21 | vs Bcl-xL + EV | as the original | -2.73 (0.0063; 0.24; 15) | -2.57 (0.010; 0.27; 22) | 5 (4 up, 1 down) |
| 21 | vs Bcl-xL + EV | constructs out | -1.98 (0.048; 0.65; 42) | -2.57 (0.010; 0.37; 16) | Esrra, Ppara, Nrf1, Srebf1 |
| 21 | vs EV + EV | as the original | -2.71 (0.0068; 0.28; 14) | -2.16 (0.031; 0.48; 38) | 4 (3 up, 1 down) |
| 21 | vs EV + EV | constructs out | -2.00 (0.046; 0.62; 42) | -2.16 (0.031; 0.55; 33) | Esrra, Ppara, Nrf1 |

Web-service network; BH across the 587-599 TFs scored; rank by |score|. **Twist1 leaves the top 25
in every variant once the constructs are out, and it never passed BH in any.** Snai1 is unmoved.
The 67- and 21-sample fits agree, so the sample set does not decide the result.

The top 25 with the constructs out (21 samples, before the accession targets are resolved, so
these will move):
- **vs Bcl-xL + EV** (25th |score| 2.30). Up: Esrra +7.32, Ppara +4.94, Nrf1 +4.10, Srebf1 +3.74,
  Ppargc1a +3.26, Arid4a +2.95, Sp3 +2.83, Nfya +2.80, Sp1 +2.75, Nfyb +2.74, Zfpm2 +2.72,
  Srebf2 +2.65, and eight more. Down: Nr3c2 -2.86, Pou5f1 -2.61, Snai1 -2.57, Pax6 -2.55,
  Fosl1 -2.37.
- **vs EV + EV** (25th |score| 2.24). Up: Esrra +7.46, Ppara +5.19, Nrf1 +3.99, Ppargc1a +3.47,
  Srebf1 +3.22, Sp4 +3.09, Zfpm2 +3.05, Zbtb14 +2.67, Arid4a +2.59, Btg2 +2.53, Nfya +2.48,
  Srebf2 +2.33, and five more. Down: Nr3c2 -3.28, E2f1 -3.01, Zbtb17 -2.67, Pou5f1 -2.45,
  Sfpq -2.40, Maz -2.39, Kat5 -2.25, Fosl1 -2.24.

What the corrected panel shows is PGC-1a's own programme, the only TFs that pass BH (ERRa, PPARa
and NRF1, with SREBP1 against Bcl-xL + EV), and a handful of nominal falls.

### 4.5 What the Twist sentence can still rest on

DESeq2 Wald, 21-sample fit, DESeq2's default BH (the script will report IHW, the paper's
convention):

| gene | vs Bcl-xL + EV | vs EV + EV |
|---|---|---|
| *Twist1* | -0.79 log2 (padj 8.6e-13) | -0.77 (2.9e-12) |
| *Mmp2* | -1.65 (1.6e-36) | -1.13 (1.7e-17) |
| *Mmp9* | -0.61 (0.043) | -0.28 (0.41) |
| *Ppargc1a* (construct) | +10.03 | +11.80 |
| *Esrra* | +1.79 (1.8e-32) | +1.75 (3e-31) |

The 67-sample fit gives the same within 0.01 log2. **In both CollecTRI and DoRothEA (all
confidence levels) `Mmp2` is a Twist1 target and `Mmp9` is not**, so "its transcriptional targets
Mmp2 and Mmp9" has network support for `Mmp2` only. The 5e westerns (TWIST1, MMP2, MMP9 down)
stand on their own. "Due to" is causal wording on what script 52 calls a selection experiment
("tumours that survived PGC1a expression show X").

### 4.6 The figure as assembled

ED 5 is an A4 page. Panel d occupies about 95 x 35 mm for both charts. The dominant stroke is
0.715 pt, the Fig. 3 Prism weight (`NM_LINE_PRISM[["Fig3"]]`), and type is 5-6 pt. In the
assembled CMYK file: EV + EV is near-white (cmyk 0.04/0.03/0.03/0), Bcl-xL + EV rich black,
Bcl-xL + Pgc-1a dark green (cmyk 0.89/0.30/1.00/0.19). ED 4's EV bars are near-white too, not
FIGURE_RULES' `#CCCCCC`, which lists EV + EV as open. The right-hand PNG carries a stray "E" panel
letter from an older figure.

## 5. Decisions for the author

1. **Which analysis.** Recommended: (b).
   - (a) The original as it was. Matches the figure; Twist1 stays on it, as a bar the transgene
     partly makes.
   - (b) Corrected: `Ppargc1a` and `Bcl2l1` left out of the gene statistics, accession names
     resolved, network pinned in the repo, and the 21 vector-series samples only (the repo
     excludes MYAZ, FaMY and BoMY, and both contrasts sit inside the vector series; the 67-sample
     fit gives the same Twist1). (a) is computed alongside as the record that the figure
     reproduces.
2. **If (b), what supports the Twist sentence?** Either the mRNA (4.5) and the 5e westerns, with an
   optional per-animal Prism table of *Twist1*, *Mmp2* and *Mmp9* (as script 57 exported for
   Fig. 3g), or Twist1 drawn on the TF chart below a break and labelled with its rank. Either way,
   "and Mmp9" and "due to" need the author's wording.
3. **Colour**, chosen from rendered options at printed size on an ED 5 page mock-up:
   - **A (Claude's pick):** every bar in the Bcl-xL + Pgc-1a green `#00573F`, filled where BH < 0.05
     and white where not, the Fig. 1e rule. Under (b) only the top three or four bars fill.
   - **B:** each bar in the colour of the arm where activity is higher: green, or black (Bcl-xL +
     EV), or EV + EV's colour (FIGURE_RULES `#CCCCCC`; near-white on the assembled ED 4-5).
   - **C:** the paper's espresso-white-mint diverging fill by score, the literal swap for red-blue.

## 6. The build, once decided

1. **Plumbing (Claude runs it).** Snapshot the web-service table (4.2) to
   `data/collectri_mouse/collectri_mouse_omnipath.tsv.gz` with a README (query URL verbatim, date,
   SHA-256 of the uncompressed TSV, the OmnipathR bug), the human repo's model. Pin the
   accession-to-symbol table beside it (`org.Mm.eg.db`, then UniProt for the rest, each row with
   its source). Add `decoupleR` to `00_setup_packages.R`; the script reads the snapshot, so
   OmnipathR is not needed. Point OmnipathR's log away from the repo root or gitignore
   `omnipathr-log/` (section 8).
2. **`scripts/58_orthotopic_tf_activity.R`** (Option A: Claude writes it and dry-runs it into the
   scratchpad; the author runs it) -> `results/orthotopic_tf_activity.rds`:
   - PART 0, the counts and the arms, `Pgc1a` mapped to `Pgc1a_BclxL`, 7/7/7 asserted;
   - PART 1, DESeq2 on the chosen samples, both contrasts, IHW beside BH, full tables saved. This
     is **the cohort's first DESeq2 Wald fit** (script 57 deliberately had none); the method needs
     a gene-level statistic, so the fit is scoped to these two contrasts and the header says so;
   - PART 2, the network: snapshot, `get_collectri()`'s post-processing (appendix B), names
     resolved, counts of what was resolved and lost;
   - PART 3, ULM as decided, BH across TFs;
   - PART 4, the reproduction record: the original recipe (67 samples, 11 groups, constructs in,
     names unresolved) against appendix A, with the recovery asserted;
   - PART 5, Twist1: its regulon, the driver table, its score with and without each construct;
     `Mmp2`/`Mmp9` membership in CollecTRI and DoRothEA;
   - PART 6, the gene read-out of 4.5 with 95% intervals; the per-animal export if chosen;
   - PART 7, the panel table (25 largest |score| per contrast, BH flags), and asserts on every
     number the text or legend will quote;
   - the `if (FALSE)` sandbox.
3. **Proposals** (`figures/natmetab/proposals/p5d_*.R` -> `outputs/natmetab/proposals/`), the three
   colours plus a page mock-up of ED 5; then **`figures/natmetab/edfig5d_tf_activity.R`** ->
   `outputs/natmetab/EDFig5/EDFig5d_tf_activity_vs_bclxl.pdf` and `..._vs_ev.pdf`. The Prism format
   of ED 5a and c: lines `lw_pt(0.72)`, 5-pt TF names and ticks, 6-pt titles; 25 rows at about
   2 mm, roughly 47 x 58 mm each, the same x scale on both, a zero line; arms named as on ED 5a
   (Bcl-xL + Pgc-1a, Bcl-xL + EV, EV + EV; alpha through plotmath, hyphens as the soft hyphen); TF
   names roman, since they name an inferred protein activity, not a transcript. Taller than the
   current 35-mm slot; ED 5 has room below row 2.
4. The ED 5d legend and a Methods paragraph drafted in `figures/natmetab/LEGENDS.md`;
   `FIGURE_RULES.md` updated with the decisions; one commit per phase.

## 7. For the Doc (the author's; Claude writes nothing to it)

- "ED Fig. 5a-f" -> a-e: the author corrects it.
- The Twist sentence: "due to", "and Mmp9", and what it cites (4.5, decision 2).
- No ED 5 legend and no TF-activity Methods yet: drafts come with the build.

## 8. Side effects and scratch files

- Scratch only, in this session's scratchpad (it may not survive): the two probe scripts, their
  reports, the DESeq2 and ULM tables, the web-service TSV, and the downloaded ED 4 and ED 5 PDFs.
  Nothing in them is needed: every number is in this note, and script 58 rebuilds all of it.
- The OmniPath web service and UniProt REST were queried from this machine (public, read-only).
- OmnipathR writes `omnipathr-log/` into the working directory. Claude's runs wrote two logs into
  the repo root's copy and deleted them; the third (00:50) predates those runs, came from another R
  session, and was left. The folder is untracked.
- `~/Library/Caches/OmnipathR` gained the 2023 static mouse table.
- `scripts/external/` is untracked; whether it is committed is the author's call.

## Appendix A. The old panels' scores, read off the figure

Read from `NM_Extended_Data_Fig5.pdf` rendered at 6,000 px; approximate to about 0.1. The
reproduction check in script 58 uses these.

**PGC vs BCL-XL** (Bcl-xL + Pgc-1a vs Bcl-xL + EV): Esrra 8.3, Ppara 5.75, Srebf1 4.45, Nrf1 4.2,
Esrrg 3.9, Ppargc1a 3.85, Fezf2 3.25, Nr1h4 3.15, Nfya 3.0, Nfyb 2.85, Sp3 2.85, Srebf2 2.75,
A0A1Y7VJR5 2.65, Yy1 2.6, Zfpm2 2.55, A0A979HLR9 2.45, Pou1f1 2.4, Arid4a 2.35, Zfpm1 2.35,
Irf7 2.35, A0A0U1RNL9 -2.4, Pax6 -2.6, A3KN90 -2.9, Sfpq -3.1, Twist1 -3.2.

**PGC vs EV** (Bcl-xL + Pgc-1a vs EV + EV): Esrra 8.5, Ppara 6.0, Ppargc1a 4.25, Nrf1 4.05,
Esrrg 3.7, Srebf1 3.5, Tfam 3.5, Fezf2 3.35, A0A1Y7VJR5 3.05, Ppard 3.0, Nr1h4 2.95, Zfpm2 2.85,
Nfya 2.7, Yy1 2.65, Six4 2.4, Sp3 2.4, Nfyb 2.4, Btg2 2.4, Zfpm1 2.35, A0A0U1RNL9 -2.8,
Foxa3 -2.8, E2f1 -2.9, Twist1 -3.1, Sfpq -3.4, A3KN90 -3.5.

## Appendix B. The network as `get_collectri(organism = "mouse", split_complexes = FALSE)` builds it

`decoupleR::get_collectri()`'s own post-processing (decoupleR 2.17.0), applied to the web-service
table `ws`. The mouse table has no `COMPLEX` rows, so the AP1/NFKB branch does nothing here.

```r
# the TSV writes True/False; decoupleR's case_when() expects 1/0
for (k in c("is_stimulation", "is_inhibition"))
  ws[[k]] <- as.integer(toupper(as.character(ws[[k]])) %in% c("TRUE", "1"))
cols <- c("source_genesymbol", "target_genesymbol", "is_stimulation", "is_inhibition")
ci <- ws[!grepl("COMPLEX", ws$source), cols]
cc <- ws[grepl("COMPLEX", ws$source), cols]
cc$source_genesymbol <- dplyr::case_when(
  grepl("JUN", cc$source_genesymbol) | grepl("FOS", cc$source_genesymbol) ~ "AP1",
  grepl("REL", cc$source_genesymbol) | grepl("NFKB", cc$source_genesymbol) ~ "NFKB")
net <- rbind(ci, cc) |>
  dplyr::distinct(source_genesymbol, target_genesymbol, .keep_all = TRUE) |>
  dplyr::mutate(mor = ifelse(is_stimulation == 1, 1, -1)) |>
  dplyr::transmute(source = source_genesymbol, target = target_genesymbol, mor)
```

An accession-named node matches
`^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$` (UniProt's accession
format).
