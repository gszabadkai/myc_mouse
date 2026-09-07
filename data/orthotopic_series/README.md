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

| group | n | series |
|---|---|---|
| `EV` | 7 | **vector series** (the scoring set) |
| `BclxL` | 7 | **vector series** |
| `Pgc1a` | 7 | **vector series** |
| `NTEV` | 7 | p21 series — out of scope |
| `NTPgc1a` | 7 | p21 series — out of scope |
| `KOEV` | 7 | p21 series — out of scope |
| `KOPgc1a` | 7 | p21 series — out of scope |
| `MYAZ` | 6 | unidentified series — excluded |
| `FaMY` | 6 | unidentified series — excluded |
| `BoMY` | 6 | unidentified series — excluded |

**The scoring set is the vector series alone, n = 21**, fixed before any number was computed
(`docs/2026-09-07_orthotopic_analysis_plan.md` Gate 1 Q3). `ox_rel` is a relative score and
GSVA is cohort-relative, so this choice sets every value; the wider set is a named sensitivity
for the composition table only.

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

1. **Different animal-ID namespace.** All 21 vector and 28 p21 samples are `BRNS###-#x`;
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
