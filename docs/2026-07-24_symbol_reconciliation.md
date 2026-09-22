# Gene-symbol vintage reconciliation (2026-07-24)

## What was wrong

The library GMTs and MitoCarta 3.0 (Sheet 4) carry the **original** MitoCarta gene
symbols. This project's expression annotation -- the count matrix, `combined_df$mgi_symbol`
and `ortholog_table$external_gene_name` -- uses the **current** symbols. For most genes the
two agree, but a family of renamed genes does not, dominated by **ATP synthase**:

| old (MitoCarta) | current (data) |
|---|---|
| Atp5a1 / Atp5b / Atp5c1 / Atp5d / Atp5e | Atp5f1a / f1b / f1c / f1d / f1e |
| Atp5o / Atp5h / Atp5j / Atp5j2 / Atpif1 | Atp5po / pd / pf / mf / if1 |
| Atp5g1 / g2 / g3 / Atp5k / Atp5l | Atp5mc1 / mc2 / mc3 / me / mg |

A plain symbol match therefore **silently dropped** the renamed genes during *set
aggregation* (shares, mitoPPS pathway scores, fGSEA/GSVA set membership): **17 of Complex
V's 24 genes, ~12% of nuclear OXPHOS**. Scattered non-ATP-synthase renames cost other
pathways 1-6% (negligible).

**Nothing is reversed.** The DESeq2 per-gene results were always complete (every gene tested,
keyed by Ensembl). Only set-level numbers were affected, and the dropped ATP-synthase genes
are *more* strongly Myc-induced than the rest of OXPHOS (+0.57 vs +0.42 mean LFC at 6W, 100%
up), so every affected number was a **conservative under-count**, never an overstatement.
Cell-death fGSEA (script 12) was already immune (its sets and ranks both come from the
current-symbol ortholog map).

## The fix

`functions/reconcile_gene_symbols.R` -- resolves a set's symbols (any vintage) by priority
fallback: current-symbol -> MitoCarta Sheet 2 EnsemblGeneID -> org.Mm.eg.db ALIAS. Idempotent
on current names, never exceeds membership, mt-* and unknown names pass through safely.
- `recon_to_ensembl(symbols, universe)` for Ensembl-keyed data.
- `recon_to_current(symbols, universe)` / `recon_current_map(symbols)` for symbol-keyed data.

Wired in: **32** (`ens_of`), **08** (mitoPPS gene->pathway), **15** (GSVA sets), **20**
(fGSEA-per-category sets), **29 / 31** (`ens_of`), and the figure scripts **fig01b, figS1b,
figS2, figS2b**. Recovery: OXPHOS_NU 134->152, Complex V 7->24, Nuclear MitoCarta 1028->1073.

## Re-run order (Option A -- author runs; DESeq2 layer NOT re-run)

The DE contrasts (scripts 01-03) are unchanged, so **no DESeq re-fit**. Re-source only the
edited set-aggregation scripts, then the figures:

1. **Independent:** `08_mitoPPS_analysis.R`, `15_gsva_scoring.R`, `20_fgsea_percategory.R`,
   `29_attenuation_decomposition.R`.
2. `31_attenuation_mechanism.R` (reads 29's output).
3. `32_mito_content_proxies.R` (reads 08 + 29).
4. Re-render the figures: `fig01`, `fig01b`, `figS1`, `figS1b`, `figS2`, `figS2b`.

**Optional (full consistency of the exploratory chain):** scripts that consume 08/15/20/29
outputs -- 22, 24, 25, 30, 33, 34, 35, 36, 37 -- inherit the corrected sets. They are
exploratory and the change is conservative/conclusion-neutral, so re-run them only if you want
every derived exploratory number regenerated.

**What to check after re-running:** in `mito_content_proxies.rds`, `MITOCARTA_OXPHOS_NU`
resolves to ~152 genes (was 134); in mitoPPS, the Complex V pathway has ~24 member genes.
Every conclusion should hold, with OXPHOS-family effects slightly larger than before.

## Dependency note

The reconciler needs `org.Mm.eg.db` (add to `00_setup_packages.R` if not already installed):
it is a Bioconductor annotation package. `readxl` is already used by script 08.
