# Orthotopic implant series — pinned data

Orthotopic mammary implants of MMTV-Myc tumour-derived (MYAZ) cells, 67 samples in ten
groups. Pinned for the `experimental-cohorts` branch. **This is not the 6W/12W purified-MEC
dataset** that scripts 00-48 analyse, and nothing in this folder feeds them.

## Provenance

- Source file: `/Users/gs/G/data/MK_myc_2022/orth_tumour_data/salmon.merged.gene_counts.tsv`
- Copied byte-identical on 2026-09-07; not edited here.
- MD5 `fba3f8fb41de4d4a2c2b76091f933c69`, 16,395,950 bytes, source mtime 2026-03-15.
- 78,334 genes x 67 samples. Columns 1-2 are `gene_id` (Ensembl `ENSMUSG...`) and
  `gene_name`; columns 3-69 are samples.
- The source directory also holds that project's own MSigDB/ORA contrast outputs and empty
  `TF_results/` `TF_plots/` folders. None of it is pinned and none has been used.

## Groups

**CORRECTED 2026-09-09.** The arm the count matrix labels `Pgc1a` is **PGC1a + Bcl-xL**, not
PGC1a alone. See `docs/2026-09-09_orthotopic_identity_correction.md`. It is written
`Pgc1a_BclxL` throughout this file; the raw column headers in the TSV are unchanged and still
say `Pgc1a`, so any script must map the label, not trust it.

| group | n | construct(s) | series |
|---|---|---|---|
| `EV` | 7 | empty vector | vector series |
| `BclxL` | 7 | Bcl-xL + EV | vector series |
| `Pgc1a_BclxL` | 7 | **Bcl-xL + PGC1a** | vector series |
| `NTEV` | 7 | non-targeting CRISPR + EV | CRISPR series |
| `NTPgc1a` | 7 | non-targeting CRISPR + **PGC1a alone** | CRISPR series |
| `KOEV` | 7 | p19^ARF KO + EV | CRISPR series |
| `KOPgc1a` | 7 | p19^ARF KO + PGC1a (transgene lost) | CRISPR series |
| `MYAZ` | 6 | none recorded | unidentified series — excluded |
| `FaMY` | 6 | none recorded | unidentified series — excluded |
| `BoMY` | 6 | none recorded | unidentified series — excluded |

### Derivation of the vector series

MYAZ was transduced with **Bcl-xL**, and **then split**: one half received EV, the other
PGC1a. Both were polyclonally selected. So `BclxL` and `Pgc1a_BclxL` are **the same starting
pool**, and **BCL-XL protein is equal between them by western**.

Two consequences, and they run in opposite directions:

- `Pgc1a_BclxL` vs `BclxL` **is a valid PGC1a contrast for every gene except `Bcl2l1`
  itself**, because the Bcl-xL background is shared and matched at the protein level.
- **`Bcl2l1` in the 14 construct-carrying samples is not interpretable.** Gene-level counts
  cannot separate construct from endogenous transcript in either arm, and equal protein means
  the transcript difference is not differential construct expression.

`EV` is the only construct-free arm of the three, so `Pgc1a_BclxL` vs `EV` differs by **two**
manipulations, not one. Every conclusion in scripts 50 and 51 that read that contrast as
"PGC1a alone" is void.

### The CRISPR series carries the only construct-free PGC1a arm

`NTEV` / `NTPgc1a` are a **non-targeting CRISPR clone** background, and `NTPgc1a` is a
**surviving PGC1a-only arm** — `Ppargc1a` 128.6 CPM, with some increase in OXPHOS subunits by
western. It is therefore the **only place in the cohort where `Bcl2l1` means endogenous
BCL-XL**.

`KOEV` / `KOPgc1a` are a **p19^ARF CRISPR clone**. `KOPgc1a` has **lost the transgene
entirely** (`Ppargc1a` 0.2 CPM), which is itself informative: it belongs to the escape series
rather than being a failed arm.

Note the two series differ in derivation — the CRISPR arms are clones, the vector arms a
polyclonal pool — so **clean contrasts are within series** (`NTPgc1a` vs `NTEV`,
`Pgc1a_BclxL` vs `BclxL`). Cross-series contrasts carry a derivation difference and are
descriptive only.

### The correct labels existed the whole time

The manuscript's **Figure 5 legend names these arms correctly** — *EV EV, Bcl-xL EV, Bcl-xL
Pgc1a*. Nothing had to be discovered to avoid this error; the sample labels simply were never
reconciled against the manuscript's own figure legends. That check is now a standing rule in
`docs/experimental_cohorts_branch_notes.md`.

### Scope

The original plan (`docs/2026-09-07_orthotopic_analysis_plan.md`) put the CRISPR series out of
scope. **That is reversed as of 2026-09-09, deliberately and with reason:** `NTEV`/`NTPgc1a`
carry the only construct-free PGC1a contrast in the cohort, and `KOEV`/`KOPgc1a` complete the
escape series. Both return.

The 21-sample vector-series scoring behind scripts 50 and 51 stays on the record as it is.
Anything involving the CRISPR arms re-scores **all 49 vector + CRISPR samples in one run** (21 vector + 28 CRISPR),
because `ox_rel`, GSVA and mitoPPS are cohort-relative. **Never quote a value from one scoring
run beside a value from the other.**

## Trap: the group tokens are inconsistent

Sample names are `<plate><n>-<GROUP><animal id>`, but the group token is not written the same
way twice:

- `A1-EV_BRNS336-5c` against `A1-EVBRNS364-4h` — the underscore is optional
- `Pgc1a` / `Pgc`, `BoMY1` / `BoMY`, `NTPgc1a` / `NTPgc`, `KOPgc1a` / `KOPgc`, `FaMY1`

**Match longest-token-first.** `EV` and `Pgc` must be matched *after* `NTEV`, `KOEV`, `NTPgc`
and `KOPgc`, or four groups silently collapse into two. Splitting on `_` or `-` does not work.
`scripts/50_orthotopic_vector_series_scoring.R` PART 0 does this and asserts zero unparsed
samples and group sizes 7/7/7/7/7/7/7/6/6/6.

## The three unidentified groups

**`MYAZ` is identifiable**: MMTV-Myc tumour-derived implant line, the parental background the
vector series is built on (`docs/2026-07-26_introduction_alignment_and_the_question.md:36-58`;
`cell_line_data/`). **`FaMY1` and `BoMY` are not determinable from files on disk** — no sample
sheet, no metadata, no coldata exists for this cohort, and the only text mention anywhere is
the analysis plan itself.

Three file-level facts put all three in a separate series, without inferring what they are:

1. **Different animal-ID namespace.** All 21 vector and 28 CRISPR samples are `BRNS###-#x`;
   these are bare (`-410-6h`, `-473-3f`, `-434-4g`) or `BFGE###-#x`.
2. **The originating analysis never crossed them.** The contrasts on disk are
   `BoMY_vs_FaMY`, `BoMY_vs_MYAZ`, `FaMY_vs_MYAZ` — and separately `BclxL_vs_EV`,
   `Pgc_vs_BclxL`, `NTPgc_vs_NTEV`, `KOEV_vs_NTEV`. **There is no MYAZ-vs-EV contrast.**
3. n = 6 each against n = 7 for every `BRNS` group, with `Ppargc1a` at 0.0 CPM and `Myc`
   1,280-1,370 against 1,440-2,368 in the `BRNS` groups.

They must not be pooled into the vector series. If provenance turns up later, that is a
decision to revisit explicitly, not silently.

## Scale, and what is NOT available

The matrix holds **salmon merged gene-level estimated counts**, which are **non-integer**.
Library sizes run 18.9-32.9 M.

**No transcript-level quantification exists anywhere on disk** — no `quant.sf`, no salmon
output directories, no `tx2gene`, no transcript count matrix, under
`/Users/gs/G/data/MK_myc_2022` or in this repo. **`tximport` is therefore unavailable**, and
the analysis goes through `DESeqDataSetFromMatrix` on rounded counts, which **silently
discards the transcript-length offset** that the `tximport` pathway would carry. This is a
stated limitation of every result from this cohort, not a detail.

Handling that follows from it:

- counts are **rounded** for `DESeqDataSetFromMatrix`, then `estimateSizeFactors`;
- z-score constructions (`ox_rel`, `ox_lvl`, `ox_mt`) take **log scale**, `log2(counts + 1)`;
- mitoPPS, if ever run, takes the **linear** normalised counts and **its own input object**.
  The two never share one.

## Trap: the filename

The analysis plan writes `salmon_merged_gene_counts.tsv`. The file is
**`salmon.merged.gene_counts.tsv`** — dot-separated, nf-core's convention. Load by the full
path under `data/orthotopic_series/`; never glob.

## Git

The TSV is gitignored, following this repo's convention for large data files
(`data/FULL.DAT.csv`, `data/fatpad_timeline/`, the MitoCarta workbooks). This README is
tracked. To reproduce the pin, copy the source file above and check the MD5.
