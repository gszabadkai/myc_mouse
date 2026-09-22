# Claude Code kickoff — script 51, the specificity and normalisation pass

Paste as the first message of a fresh session in `myc_mouse`. Start in Plan Mode.

---

## Orientation

Repo `myc_mouse`, branch `experimental-cohorts`. Script 50 is run and signed off.

Read, in this order:

1. `docs/2026-09-08_orthotopic_vector_series.md` — the result. Authoritative.
2. `docs/2026-09-07_orthotopic_analysis_plan.md` — the spec 50 was built from.
3. `scripts/50_orthotopic_vector_series_scoring.R` — house pattern to follow.
4. `CLAUDE.md` — workflow, coding rules, data traps. Its branch section is stale
   and predates `paper-final`; everything else holds.

**Nothing here reopens a verdict.** C1 holds, C2 holds as an observation, C3 is
UNINFORMATIVE and stays that way. This pass adds specificity checks on results
already obtained, one normalisation robustness check, and one correction to how
C3 was constructed. It is not a fourth claim, and it must not be written as one.

---

## Mode

Plan Mode. **Gate 0 first — report and stop.** Option A throughout: you write
the numbered script, you do not run it. File inspection and git are yours.

---

## Gate 0 — a targeted file check, before anything is planned

§1 of the note records that no transcript-level quantification exists, from a
search for `quant.sf`, salmon output directories and `tx2gene`. **That search
would not have found the files below**, which nf-core writes alongside the gene
matrix under different names. Check for these specifically, by name, in
`/Users/gs/G/data/MK_myc_2022/orth_tumour_data/` and its siblings:

- `salmon.merged.gene_counts_length_scaled.tsv` (or `_scaled`)
- `salmon.merged.transcript_counts.tsv`
- any `.rds` SummarizedExperiment objects from the pipeline
- `salmon.merged.gene_lengths.tsv` or equivalent

Report what exists, with paths and MD5s. Then stop.

Why it matters, so you can see the stakes:

- **length-scaled counts** carry the transcript-length correction and are
  documented as safe for `DESeqDataSetFromMatrix` without an offset. If that
  file exists, the note's headline limitation largely dissolves.
- **transcript counts** would separate **Bcl-xL from Bcl-xS**, which arise from
  alternative 5' splice site usage at the same locus and are functionally
  opposite. Gene-level `Bcl2l1` pools them. For a paper resting a causal claim
  on `Bcl2l1` that is a sharper limitation than the length offset, and it would
  change what PART B below can claim.

If either exists, say so and wait — the build changes and I will decide. If
neither exists, say that plainly and we proceed on the current matrix.

---

## Script 51 — `scripts/51_orthotopic_specificity_checks.R`

House pattern of 48/49/50: long header stating the estimand and the fixed
reading rules, `set.seed(1)`, asserts, NOTES, `saveRDS`, `if (FALSE)` sandbox at
the end ordered gate-first.

Same fixed cohort as 50 — the 21 vector-series samples, EV / BclxL / Pgc1a.
Rebuild it from source and assert reproduction rather than trusting the rds.

**PART 0 — load control.** Reproduce script 50's C1 table and the §2 median-CPM
table against hard-coded reference values, `stopifnot` before anything new is
computed. If 51 cannot reproduce 50, stop; nothing below is readable.

**PART A — arm geometry.** Mean and per-sample `ox_rel`, `ox_lvl`, `Bcl2l1`,
`Mcl1` by arm. This exists to confirm or refute one specific explanation: that
the near-zero pooled C3 rho arises because the three arm-means form a
non-monotone triangle in the (`ox_rel`, `Bcl2l1`) plane — EV low-low, Pgc1a
high-mid, BclxL low-high — so the EV→Pgc1a limb and the BclxL arm cancel.
Report the arm-mean positions explicitly. Plot it.

**PART B — the C1 panel extended.** Same machinery as 50's C1 — Pgc1a vs EV,
n = 7 v 7, Wilcoxon plus difference in medians with bootstrap interval — applied
to `Mcl1`, `Bbc3` and the `Mcl1:Bcl2l1` group contrast, which 50 never ran.
Report the remaining BCL2-family twelve alongside as a panel. `Bcl2l1` repeated
as the anchor. **No per-gene FDR** across the twelve; estimates with intervals.
BclxL reported as the third arm and never pooled, exactly as 50 did.

**PART C — the asymmetry check.** `ox_lvl` and `ox_rel`, **BclxL vs EV**, same
machinery. Also decompose the PGC1α `ox_lvl` move into numerator and denominator
so it is visible how OXPHOS-selective the intervention was, rather than left as
an inference from the gap between +1.925 and +0.915.

**PART D — normalisation.** Three things:
1. Diagnostic — `sizeFactors(dds)` against `colSums(counts)/1e6` per sample,
   ratio summarised **by arm**; plus the share of each library sitting in
   MitoCarta genes, by arm. This measures the CPM deflation directly rather than
   arguing about it.
2. C2 recomputed on `counts(dds, normalized = TRUE)` instead of CPM — max and
   80th-percentile contrasts, same exact permutation null over all 3432
   relabellings.
3. Script 50's PART F composition-adjusted C3 fit recomputed on median-of-ratios
   rather than CPM. That is the one C3 reading where the normalisation choice
   bites, and it is the reading carrying the sign flip.

**PART E — threshold or dose-response.** Within the Pgc1a arm only, n = 7:
spread of `Ppargc1a`, and whether `Bcl2l1` tracks it. Report the correlation
with its interval and the per-sample values. This discriminates a viability
threshold from a dose-responsive induction, which is what decides how C1 is
worded.

*Optional external dose point, and hedged hard:* `NTPgc1a` carries `Ppargc1a` at
129 CPM against `NTEV`, with `Bcl2l1` 79.3 vs 81.6. Raw CPM is not
cohort-relative, so this is readable without rescoring and without touching the
p21 series' out-of-scope status. Report it **as a raw-CPM observation only** —
never scored, never tested, never in a claim, and flagged as a different genetic
background. It is in the deposited matrix and a reader will find it.

**PART F — the C3 construction correction.** Descriptive, no new verdict.
Recompute the pooled C3 readings on EV + Pgc1a alone and show that they
restate C1 rather than adding to it. The point being recorded: C1's spec said
BclxL is reported and **never pooled**; C3's spec pooled it, and seven of the
21 `Bcl2l1` values are set by an expression construct rather than by the tumour.
Same endpoint, opposite handling, caught only because C1 and C3 were
pre-specified as separate claims.

**PART G — verdict table** on the rules fixed in the header, then `saveRDS` to
`results/orthotopic_specificity_checks.rds`.

---

## Reading rules — fixed here, before any number is seen

**R1 — `Mcl1` in C1.** Flat, or moving far less than `Bcl2l1` → C1 is a
selective event and the guardian-dissociation sentence is available. Moving with
`Bcl2l1` at comparable magnitude → C1 is a general anti-apoptotic response, the
specificity sentence comes out, and C1 is weaker.

**R2 — `Bbc3` in C1.** Flat or falling while `Bcl2l1` rises → the
guardian–sensitiser gap widens by intervention and C1 strengthens. Rising
proportionally → the gap is unchanged and C1 is a general BCL2-family response.

**R3 — the asymmetry.** Both rulers flat in BclxL vs EV → the asymmetry holds
and the directed-edge argument is available. `ox_lvl` up with `ox_rel` flat →
ruler-specific; wording restricted to respiratory *share*, not respiration.
`ox_lvl` up → no asymmetry and the argument comes out. **In all three cases**
this must be reconciled against the standing experimental finding that Bcl-xL
overexpressing tumours have higher OXPHOS — flag in the note whether that
finding was transcriptional or functional, because a functional-versus-transcript
discordance is a different statement from an asymmetry.

> **WITHDRAWN 2026-09-09 — THE FINDING NAMED IN THE PARAGRAPH ABOVE DOES NOT EXIST.**
> "The standing experimental finding that Bcl-xL overexpressing tumours have higher
> OXPHOS" originated **here**, in this brief, and was never in the data or in the
> manuscript. The manuscript has Bcl-xL as the **permissive** partner throughout — it
> rescues the Ndi1 growth defect and makes raised respiration advantageous — and nowhere
> claims it raises respiration. Instructed to reconcile against it, script 51's note
> faithfully wrote a blocker forbidding the R3 argument; that blocker is deleted (see
> `docs/2026-09-08_orthotopic_specificity_checks.md` section 3) and R3 stands, demoted.
> Recorded here rather than quietly removed, because a brief is worth more when it carries
> its own errors. See `docs/2026-09-09_orthotopic_identity_correction.md` section 1.

**R4 — C2 under median-of-ratios.** Survives → the ceiling is biology and C2
stands as an observation. Attenuates materially → CPM deflation was carrying it
and C2 weakens or comes out. State which before quoting the number anywhere.

**R5 — PART E.** Wide `Ppargc1a` spread with flat `Bcl2l1` → threshold reading,
consistent with selection. `Bcl2l1` tracking `Ppargc1a` → dose-responsive
regulatory reading. Narrow `Ppargc1a` spread → uninformative, and say so rather
than reading the flatness.

---

## Non-negotiables

- **These are specificity and robustness checks on results already obtained,
  not new hypothesis tests.** Nothing here is pre-registered as a claim, and
  nothing here becomes a fourth claim in the manuscript. It qualifies C1's
  wording and it corrects C3's construction.
- **C3 stays UNINFORMATIVE.** PART D.3 and PART F are note-quality work. If the
  sign stabilises, the interval still spans zero and the verdict is unchanged.
  The abstract sentence stays withdrawn.
- **No MYC × OXPHOS interaction.** Not estimable here at any n; no group term
  stands in for it.
- **N3.** Transcript associations. "Primed" appears nowhere — code, comments,
  labels or prose.
- **Species = cohort.** No pooling with the 6W/12W timeline or human. Gland
  +0.351 and human −0.31/−0.46 are directions only.
- **No per-gene FDR** across the twelve. Estimates with bootstrap intervals,
  5,000 resamples, and the C2 null exact rather than sampled.
- Input object separation: GSVA log-scale, mitoPPS linear. Neither is run here.
- MitoCarta by exact filename, Sheet 4, `mt-*` stripped as script 08 does.
- R rules: no `print(n=X)` after `head()`; `dplyr::count()`; ASCII only;
  `here::here()`; deps via `00_setup_packages.R`.

---

## Deliverables

1. `scripts/51_orthotopic_specificity_checks.R`, sourced by me in Positron.
2. `results/orthotopic_specificity_checks.rds` carrying `$load_control`
   `$arm_geometry` `$c1_panel` `$asymmetry` `$norm_diagnostic` `$c2_mor`
   `$c3_partf_mor` `$dose` `$c3_ev_pgc1a` `$verdict` `$notes`.
3. A dated `docs/` note in the house pattern, every number re-read from that
   object rather than transcribed. It records the R1–R5 outcomes, and it states
   the C3 construction correction as a correction to
   `docs/2026-09-08_orthotopic_vector_series.md` — cross-referenced from that
   note, not silently edited into it.

Also carry to the branch notes, alongside the existing fat-pad corollaries: the
gate's pooled |rho| threshold has a blind spot in **both** directions. §2 of the
09-08 note showed pooled-material with per-group disagreeing; `ox_lvl` against
proliferation shows the mirror — pooled −0.323 with EV/BclxL/Pgc1a at
−0.643 / −0.679 / −0.964, three arms agreeing strongly and the pooled value
below threshold. Read the per-group breakdown before believing either verdict.
