# Gray et al. (2023) developmental TF selection — Lane 2 of `07_build_tf_target_sets.R`

Memory note on how the **Gray-context CHEA lane** (Lane 2) of the TF target panel
was built. Lane 2 is the mammary-context, per-sample-scorable arm of Category 6
(`TF_targets`); it also runs the Decision-5 mito-overlap discovery that feeds the
Category-7/9 oncogenic-vs-developmental and apoptotic-sensitivity threads.

Source: Gray GK, et al. (Brugge lab) *A human breast atlas integrating single-cell
proteomics and transcriptomics*, Dev Cell 2022 — mouse mammary supplements
`Gray-et-al_2023_Suppl3.xlsx` and `Suppl5.xlsx` in `data/raw/mammary_dev_atlases/`.

## Inputs

Both supplements are **CHEA3 upstream-regulator enrichment tables** — for each
mammary cell-state expression signature, the TFs whose ChEA3 target sets are
enriched in that signature, ranked.

| File | Sheets (contexts) | Lifecycle |
|---|---|---|
| Suppl3 | `APhi/APlo/HShi/HSlo/BAhi/BAlo_CHEA` | adult major types, high/low-expressing |
| Suppl5 | `AP/HS/BA` x `TEB/Duct _CHEA` | pubertal terminal-end-bud vs duct |

Each sheet = 250 ranked TFs; columns `Rank, TF, Score, Library, Overlapping_Genes`.
833 distinct TFs across all 12 sheets. `Overlapping_Genes` is the comma-joined list
of signature genes overlapping that TF's ChEA3 targets, in **UPPERCASE human-style
symbols** — mapped to mouse via `R/03_ortholog_utilities.R::human_to_mouse`
(expand one-to-many; the whole union is mapped once for efficiency and clean QC).

Major types: **AP** = alveolar progenitor (secretory lineage), **HS** = hormone-
sensing luminal, **BA** = basal/myoepithelial. Context tokens normalise the
sheet's `hi/lo` to the project's `HE/LE` vocabulary, giving 12 contexts:
`AP_HE AP_LE HS_HE HS_LE BA_HE BA_LE AP_TEB AP_DUCT HS_TEB HS_DUCT BA_TEB BA_DUCT`.

## What Lane 2 emits

### 2a. Context TF-program sets — `TFT_<TF>_GRAY_<CONTEXT>` (gsva)
The TF's full overlapping-gene program (mouse symbols) per context, kept whole for
per-sample GSVA scoring. **TF selection** (to avoid 833 x 12): the generic roster
**MYC, E2F1, ESR1, ESRRA, NRF1, GABPA** wherever CHEA ranked them, PLUS the **top
10** ranked TFs per sheet; deduplicated. Result: **142 sets**, sizes 14-788.
`source = "GRAY_CHEA"`, `method = "gsva"`, `stage = <context>`. The roster TFs only
appear where their program is active (e.g. ESR1 in 8 contexts; ESRRA only AP_LE &
HS_HE; GABPA only AP_LE) — a sensible context-specificity signal, not a bug.

### 2b. Decision-5 mito-overlap search (discovery artifact, all 833 TFs)
For **every** (TF, context) row, the mouse-mapped overlapping genes are intersected
with **full MitoCarta** (`MITOCARTA_ALL`, 1037 genes within the ortholog-mappable
universe) and annotated by sub-pathway from script 02's splits: OXPHOS,
mitochondrial ribosome, apoptosis-PRO, apoptosis-ANTI, and metabolism-other.
Enrichment is a hypergeometric test vs the mappable mouse universe, BH-FDR adjusted
across the rows with any overlap. Written to
`outputs/qc/gray_chea_mito_tf_shortlist.csv` (2168 rows / 741 TFs).

**Validation (the point of the search):** the PGC1a-axis TFs recover strong
in-context mito enrichment — `ESRRA` AP_LE 47/165 (p = 1.1e-22), `GABPA` AP_LE
50/359 (p = 1.3e-10) — and novel mitochondrial TFs surface (`MTERF3`, `CHCHD3`).

### 2c. Promoted mito-overlap sets — `TFT_<TF>_GRAY_<CONTEXT>_MITO` (gsva)
Per the decision that **any** significant TF (not just the roster) becomes a set,
each shortlist row passing **BH FDR < 0.01 AND >= 15 mito genes** is emitted as a
descriptive set whose members are that TF's mito-overlap genes (TF program ∩
MitoCarta). Result: **218 sets / 142 TFs**, sizes 15-124. `source =
"GRAY_CHEA_MITO"`, `method = "gsva"`. Captures the PGC1a axis, MYC (BA_LE, 124 mito
genes), and the novel CHCHD3/MTERF3 programs.

**Lane 2 total: 142 (2a) + 218 (2c) = 360 sets**, plus the shortlist CSV.

## Apoptotic-sensitivity layer (for Category 9)

The MitoCarta apoptosis catalog is compact: **25 PRO** (Bax, Bak1, Bid, Bbc3/PUMA,
Bcl2l11/Bim, Pmaip1/NOXA, Casp3/8/9, Cycs, Diablo, ...) and **9 ANTI** (Bcl2,
Bcl2l1/Bcl-xL, Mcl1, Bcl2a1d, ...). Per-TF counts are small (1-3), so these are NOT
emitted as standalone sets here — instead the shortlist preserves, per (TF,context):
`n_apoptosis_pro`, `n_apoptosis_anti`, `apop_balance` (= n_pro - n_anti; >0 =
pro-apoptotic/sensitising lean), and the explicit `apoptosis_pro_genes` /
`apoptosis_anti_genes` symbols. These feed the pathway-level Category-9 intersections
(`MYC_MITO ∩ APOPTOSIS_PRO`, etc.) in script 10.

Observations: 478 TF x context rows (265 TFs) carry >=1 apoptosis-mito gene; the
`HS_TEB` context is apoptosis-loaded (TEB lumen clearance by apoptosis during
pubertal elongation). Roster TFs lean pro-apoptotic: MYC -> Bbc3/Bok/Bid; E2F1 ->
Bak1/Bid/Pmaip1; ESR1/ESRRA -> Aifm2; GABPA -> Bak1; NRF1 -> Bbc3. I.e. the
biogenesis TFs that build mito mass also touch pro-apoptotic effectors — a candidate
co-regulation of apoptotic sensitivity.

## Decisions log

- **2026-06-18** Granularity: roster + top-10 per sheet for context sets.
- **2026-06-18** Naming: `TFT_<TF>_GRAY_<CONTEXT>` and `..._MITO`; CHEA rank/score
  kept verbatim in the set description.
- **2026-06-18** Mito-set promotion gate: **BH FDR < 0.01 & n_mito >= 15**; applies
  to ANY CHEA TF, not just the PGC1a-axis roster.
- **2026-06-18** Apoptosis effectors recorded in the shortlist (genes + balance) for
  the Category-9 thread rather than emitted as tiny standalone sets.

## Caveats

- ChIP-Atlas regulons are the *generic* lane (Lane 1), not this one. Lane 2 is the
  mammary-context CHEA lane.
- `Overlapping_Genes` reflects expression-signature overlap with ChEA3 targets, not
  direct binding in mammary tissue — read enrichment as "TF-program associated".
- ER's broad cistrome makes ESR1 programs large/noisier.
- Generic transcription machinery (GTF3A, TBP, SNAPC5) appears in the mito shortlist
  by sheer target breadth; under "include if significant" these come along as
  descriptive sets — interpret high-overlap general factors with caution.
