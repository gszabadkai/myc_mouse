---
name: mitocarta-sets-are-membership-loose
description: MitoCarta process sets contain loosely-assigned multifunctional members that can carry a whole set-level effect; resolve any set-level MECHANISM claim gene by gene before believing it
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 71ffe8cc-b22d-412b-9e40-f22881e23306
---

**Never turn a MitoCarta set-level effect into a MECHANISM claim without resolving it gene by
gene first.** Learned the hard way 2026-07-16: I shipped a "commissioned but unbuilt" finding in
script 32 (Myc builds the mtDNA machinery, output does not follow) built on
`MITOCARTA_TRANSCRIPTION` +43% (p=8.5e-5) and `MITOCARTA_MTDNA_NUCLEOID` +22% (d=2.32, p=8.6e-6).
**Both were composition artifacts**, and the author had already sourced the script before I caught it.

Worked examples in this repo (share of non-mt transcriptome, Myc main effect):
- `MITOCARTA_TRANSCRIPTION` (9 genes) is carried by **Mrpl12 +116%** -- a *mitoribosomal* protein
  that moonlights with POLRMT -- plus Top1mt +76%, Mtres1 +50%. The actual core: Tfam +11% (ns),
  Tfb2m +16% (ns), Tefm +18% (ns), Mterf1a/1b ~0 (ns); only Polrmt +18% is nominal (p=0.036).
- `MITOCARTA_MTDNA_NUCLEOID` (11 genes) is carried by **Atad3a +113%** and **Poldip2 +58%** (both
  multifunctional; mito role secondary). The core: Polg +3% (ns), Twnk +1% (ns), Polg2 -9% (ns).

**Why:** MitoCarta process assignments are inclusive -- a gene lands in a process if it participates
at all. Small sets (9-15 genes) + a share/mean composite weighted toward high-expressed members means
**one or two multifunctional, strongly-regulated genes can carry the entire set-level effect**. The
set is fine for enrichment; it is NOT evidence about the pathway's dedicated machinery.

**How to apply:**
1. Any set-level claim of the form "Myc regulates process X's machinery" -> print the members with
   per-gene effect AND mean expression before writing a word of interpretation.
2. Watch for the tell: a small set with a large `geno_d` (>2) is usually one gene, not a module.
3. Having retracted a claim, **do not immediately substitute the opposite claim** -- the mirror story
   ("Myc SPARES the mtDNA machinery") was also unsupported here (39th percentile of the nuclear mito
   distribution, Wilcoxon p=0.12). Check it with the same rigour; an informative NEGATIVE is often the
   honest result.
4. Script 32 PART 3b (`mtdna_machinery_genes`, `mtdna_core_vs_compartment`) is the enforcement pattern
   -- reuse it.

Applies to the mito sets generally, and to any small process set in
`data/genesets_from_library/mammary_mito_myc_metab_v1_mouse.gmt`. See [[mito-content-bounded]].
