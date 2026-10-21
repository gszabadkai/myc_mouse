---
name: gene-symbol-vintage-reconciliation
description: MitoCarta/GMT old symbols vs current-symbol expression drop renamed genes (ATP synthase) during set aggregation; fixed by functions/reconcile_gene_symbols.R
metadata: 
  node_type: memory
  type: project
  originSessionId: 71ffe8cc-b22d-412b-9e40-f22881e23306
  modified: 2026-07-24T11:30:25.641Z
---

**The bug (found 2026-07-24 while building figS2).** Library GMTs + MitoCarta Sheet 4 carry
ORIGINAL gene symbols; the expression annotation (count matrix, `combined_df$mgi_symbol`,
`ortholog_table$external_gene_name`) uses CURRENT symbols. A plain symbol match SILENTLY DROPS
renamed genes during SET AGGREGATION (shares, mitoPPS, fGSEA/GSVA membership). Dominated by
**ATP synthase**: Atp5a1/b/c1/d/e -> Atp5f1a..e; Atp5o/h/j/j2/Atpif1 -> Atp5po/pd/pf/mf/if1;
Atp5g1-3/k/l -> Atp5mc1-3/me/mg. Cost: **17 of Complex V's 24 genes, ~12% of nuclear OXPHOS**;
other pathways 1-6% (negligible).

**NOTHING REVERSES.** Per-gene DESeq2 results were always complete (Ensembl-keyed, every gene
tested). Only set-level numbers were hit, and the dropped ATP-synthase genes are MORE
Myc-induced than the rest of OXPHOS (+0.57 vs +0.42 at 6W, 100% up) -> every affected number was
a CONSERVATIVE UNDER-count, never an overstatement. **Cell-death fGSEA (script 12) was immune**
(sets AND ranks both from the current-symbol ortholog map). The content headline (+21-27%) is
safe (mass markers = current symbols; nuclear MitoCarta loses only ~4%).

**The fix = `functions/reconcile_gene_symbols.R`.** Resolves any-vintage symbols by priority
fallback: current-symbol -> MitoCarta Sheet 2 EnsemblGeneID -> org.Mm.eg.db ALIAS. Idempotent on
current names, never exceeds membership, mt-* + unknowns pass through safely. Exposes
`recon_to_ensembl` (Ensembl-keyed data), `recon_to_current` / `recon_current_map` (symbol-keyed).
Needs org.Mm.eg.db + AnnotationDbi (added to `00_setup_packages.R`) + readxl.

**Wired in (author must RE-RUN -- Option A):** 32 (`ens_of`), 08 mitoPPS (gene->pathway
translate), 15 GSVA + 20 fGSEA (batch-translate library sets; Hallmark already current), 29/31
attenuation (`ens_of`), + figures fig01b/figS1b/figS2/figS2b. Recovery OXPHOS_NU 134->152, CV
7->24, Nuclear MitoCarta 1028->1073. **User chose "figures + core" scope.** Re-run order + what to
verify = `docs/2026-07-24_symbol_reconciliation.md`. NO DESeq re-fit (contrasts unchanged).
Downstream exploratory consumers (22/24/25/30/33-37) inherit the fix -> optional re-run.

**How to apply:** any NEW set-based analysis must resolve set symbols through the reconciler, not
a bare `mgi_symbol %in%` / `external_gene_name %in%` match. See [[figure-phase-architecture]]
(figS2 mapping fix) and [[mitocarta-sets-are-membership-loose]] (a related MitoCarta membership trap).
