# Fat-pad timeline cohort (Chandan) — pinned data

Whole mammary **fat pad** RNA-seq across an MMTV-Myc progression series, 30 samples.
Pinned here for the `experimental-cohorts` branch. **This is not the 6W/12W purified-MEC
dataset** that scripts 00-48 analyse, and nothing in this folder feeds them.

## Provenance

- Source file: `/Users/gs/G/data/MK_myc_2022/Chandan_data/Normalised_DESseq_counts.csv`
- Copied byte-identical on 2026-09-07; not edited here.
- MD5 `98db846e89bbe9881bcd813ff0d01485`, 13,761,651 bytes, source mtime 2023-08-19.
- 54,838 genes x 30 samples. Row names are **Ensembl mouse gene IDs** (`ENSMUSG...`).
- The source directory also holds that project's own DESeq2 contrast outputs and a
  `12_WK_POS_V_6WK_POS.R`; none of it is pinned, and none of it has been read.

## Design

| group | n |
|---|---|
| `6WK_NEG` | 5 |
| `6WK_POS` | 5 |
| `12WK_POS` | 5 |
| `SMALL_TUMOUR` | 6 |
| `LARGE_TUMOUR` | 9 |

Column names carry the group and a replicate suffix (`LARGE_TUMOUR_R7`), so the group
factor is recovered with `sub("_R[0-9]+$", "", colnames(x))` — there is no separate
coldata file.

**There is no 12W negative control.** The 6W arm is genotype-controlled
(`6WK_NEG` vs `6WK_POS`) and nothing after 6W is. Any 12W-or-later statement is a
within-genotype progression reading, never a genotype contrast, and the developmental
comparison that the MEC dataset supports at 12W cannot be made here at all.

## Tissue — read this before scoring anything

This is **whole fat pad, not enzymatically purified MECs**. The distinction is the
opposite of the main dataset's: in the MEC cohort a stromal or adipocyte signal is
*contamination* (~3-12%); here adipose is the *majority tissue* throughout and only
partly dilutes. Group medians of `Adipoq`: 59,342 (`6WK_NEG`) -> 26,415
(`LARGE_TUMOUR`). It never becomes tumour tissue; it stays tumour-bearing fat pad.

Adipocytes are mitochondria-rich, so a respiratory score computed here reads the adipose
fraction unless proven otherwise. It was not proven: see
`docs/2026-09-07_fatpad_timeline_oxphos_precheck.md`, which **fails** this dataset for
respiratory-axis work and records what it is cleared to contribute instead.

## Scale — handle explicitly, do not assume

The matrix is **DESeq-normalised and LINEAR** (non-integer values; normalisation was done
upstream by the originating analysis, not in this repo). There is no `dds` object, no
size factors on disk, and no raw counts, so:

- z-score constructions (`ox_rel`, `ox_lvl`, any composite) need **log scale**:
  `log2(counts + 1)`.
- mitoPPS needs the **linear** counts, per the project convention.
- The two never share an input object. Build each from `NC` separately.

Because normalisation is upstream, a DESeq2 re-analysis of these data is **not** possible
from what is pinned here — it would need the raw counts from the originating project.

## Trap: the filename collides

`Normalised_DESseq_counts.csv` is **also** the name of the 24-sample purified-MEC matrix at
`/Users/gs/G/data/MK_myc_2022/k693_rerun/Normalised_DESseq_counts.csv`
(MD5 `1dcda2b5208f431f875abf94321d710c`, 25 columns, sample-code headers `MYCF62_3f` etc.).
Different dataset, different tissue, same basename. Load by the full path under
`data/fatpad_timeline/`; **never glob `*Normalised*`**, exactly as with
`Mouse.MitoCarta3.0.xls` and `FULL.DAT.csv`.

## Git

The CSV is gitignored, following this repo's convention for large data files
(`data/FULL.DAT.csv`, the MitoCarta workbooks). This README is tracked. To reproduce the
pin, copy the source file above and check the MD5.
