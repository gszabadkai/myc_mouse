---
date: 2026-10-28
tags: [project/myc_mouse, nat-metab, section-3, figure-3, orthotopic, script-57, handoff]
status: >
  Read-out complete. Script 57 written and run twice (content verified). Panels are
  DRAFTS for content control only -- they are not formatted and must not go into the
  figure as drawn. Four decisions open, listed in section 5. Nothing is committed.
purpose: >
  The digest of what script 57 returned, what it changes in the manuscript, and what
  is left to do. Written for picking the work back up cold.
relates-to:
  - scripts/57_section3_panel_readout.R
  - docs/2026-10-28_section3_panel_readout.md  (the auto-generated record -- REGENERATED ON EVERY RUN, never hand-edit)
  - results/orthotopic_section3_panels.rds
  - docs/2026-10-27_section3_final_narrative_and_panel_spec.md (v2) + the author's 2026-10-28 build list
  - docs/2026-09-23_section3_narrative_and_figure3.md section 2.5 (limitations)
  - figures/natmetab/FIGURE_RULES.md "The cell-line constructs" (the colours Fig. 3 should use)
---

# Section 3 read-out: what script 57 returned

## 1. The one error to fix in the manuscript

**Fig. 3g is a 796-fold difference, not ~1,100-fold.** The 1,100 is a different panel:
it is `Pgc1a_BclxL` vs `BclxL` (301.6 / 0.273 = 1,104). Fig. 3g is `NTPgc1a` vs `NTEV`,
which is 128.6 / 0.162 = **796**. Check wherever that number already appears.

The P value is settled: **Mann-Whitney P = 0.000583**, exactly the two-sided floor at
7 vs 7. Prism reproduces it from the pasted values. The placeholder "P < 0.0001" is not
reachable at this n and comes off. **The two pasted columns are not paired** -- row 1 of
each is a different animal -- so the test must be set to unpaired.

## 2. Both open numbers came back negative, and both close an open item

**Extended Data a is not a dose-response (spec section 6, item 1).** Within
`Pgc1a_BclxL`, across its seven animals, `Bbc3` does not track `Ppargc1a`: Pearson
r **-0.165 [-0.817, +0.671]**, P = 0.72; Spearman rho -0.214. The sign matches the
across-arm direction but the interval covers nearly the whole range.

So the panel is a **three-point comparison between three backgrounds** that differ in
genotype as well as in retained dose. The Extended Data sentence must say that rather
than implying a dose-response within tumours. At n = 7 this is weak evidence of absence:
the right wording is that the data do not resolve a within-arm relation, not that there
is none. Part of the reason there is no power: `NTPgc1a`'s retained dose spans five-fold
across its animals (52, 81, 105, 129, 132, 192, 272 CPM), against `Pgc1a_BclxL`'s tight
242-364. "A partial dose" is a median over a wide spread.

**ED d is not built (spec section 7).** `Myc` in `KOPgc1a` vs `KOEV`: median
**-193 [-475, +431]** CPM, Wilcoxon P = 1. The ceiling stays what the September note
called it, the weakest of the three sequenced results.

The instinct about the outlier was right, and sharper than expected: the **mean is
+271** -- the wrong sign -- and dropping the 4,611 CPM animal moves it to **-73**. A mean
would have reversed the conclusion. This is why the spec said median or mean with an
interval and why the 80th-percentile and maximum contrasts were ruled out.

## 3. The summary numbers, each arm against its own control

| route | arm | control | Ppargc1a CPM | ox_rel vs control | mitoPPS vs control | MitoCarta share | Bbc3 log2 vs control |
|---|---|---|---|---|---|---|---|
| no transgene (zero arm) | `BclxL` | `EV` | 0.2732 | +0.000 [-0.199, +0.161] | -0.021 [-0.050, +0.087] | 7.8% vs 8.0% | not fitted |
| effector blocked | `Pgc1a_BclxL` | `BclxL` | 301.6 | +1.305 [+1.165, +1.452] | +0.570 [+0.475, +0.622] | 14.7% vs 7.8% | -0.731 [-1.120, -0.450] |
| pathway intact | `NTPgc1a` | `NTEV` | 128.6 | +0.176 [+0.039, +0.348] | +0.053 [+0.030, +0.116] | 8.5% vs 8.3% | -0.345 [-0.890, -0.010] |
| checkpoint removed | `KOPgc1a` | `KOEV` | 0.1677 | +0.271 [-0.139, +0.441] | +0.071 [-0.064, +0.105] | 8.6% vs 8.1% | +0.044 [-0.231, +1.182] |

All five values the escape-series note records reproduce from the object to rounding.
The rebuild gate passed at **exactly 0** on all four quantities, so the per-animal layer
sits on precisely script 52's ruler.

**`Bbc3` for the zero arm is blank** because script 52 fitted the p53 panel only on the
three clean within-series contrasts. Fig. 3i does not need it. It is a new contrast if
it is ever wanted.

**The `NTPgc1a` sentence holds, and the tilt reading is confirmed.** `ox_rel` is
+0.176 [+0.039, +0.348], excluding zero at transcript level, in the arm whose western
shows ETC subunits unchanged -- while the compartment does not expand (MitoCarta share
8.5% vs 8.3%). A tilt, not mass expansion. That is why the prose says "OXPHOS subunit
protein" in the 3g sentence and makes 3i a dose-scaled relation.

## 4. For Methods

- **The CRISPR arms are pooled, not clones.** p19ARF-disrupted, not null; edited fraction
  checked and p19ARF western QC passed. `data/orthotopic_series/README.md` is corrected.
  Script 52's header and script 53's caption still say "clone" and were deliberately left
  alone: they are the committed analysis behind a saved object and a published PDF.
- **No composition adjustment was applied to anything.** No axis reached |rho| >= 0.5
  pooled on `ox_rel`, so there is no adjusted estimate to report. The threshold sits in
  script 52's parameter block above every computation, saved as `$params$CUT_MAT`, and is
  called pre-registered at `docs/2026-09-08_orthotopic_vector_series.md:100`.
- **Say the gate was pooled, because the per-arm values disagree.** The strongest axis,
  endothelial, is +0.472 pooled -- under threshold -- but runs from **-0.643 (`BclxL`) to
  +0.750 (`NTEV`)**, opposite signs. The September note's warning about pooled loadings
  applies here too.
- **Use the 49-sample loadings, never the 21-sample ones.** The spec's adipose figures
  (+0.621 pooled, +0.036 / -0.393 / +0.714) are the 21-sample vector run and its three
  columns are EV / BclxL / Pgc1a_BclxL, not the three PGC-1a arms. The 49-sample adipose
  on `ox_rel` is +0.329 pooled, `Pgc1a_BclxL` +0.429, `NTPgc1a` +0.036, `KOPgc1a` +0.214.
  Note the trap: `NTPgc1a`'s +0.036 happens to equal EV's value in the 21-sample table.
- The section 2.5 limitations carry over unchanged (survivor populations; no
  transcript-level quantification, so `Bcl2l1` pools Bcl-xL with Bcl-xS; retained, never
  assigned).

## 5. What is open

1. **The panels are drafts.** They were built to control content and are not formatted.
   See section 6 for what formatting they need.
2. **The facility manifest never landed.** Arm assignment still rests on the column
   headers alone. Drop the file at `docs/orthotopic_facility_manifest.csv` (or
   `docs/orthotopic_sample_manifest.csv`, or
   `data/orthotopic_series/facility_manifest.csv`) and re-run 57; PART 3 verifies it per
   sample and reports any mismatch. The construct check already passed: no sample sits
   more than 5 robust MADs from its own arm, and the only designed-vs-observed mismatch
   is `KOPgc1a`, which is the recorded transgene loss.
3. **Extended Data Fig. 5 is still a placeholder number.** The manuscript has ED 1-4 plus
   two unnumbered figures in the human section; the numbering has to be settled.
4. **Nothing is committed.** See the handoff below.

## 6. Figure formatting and colour -- where the rules already are

**The authority is `figures/natmetab/FIGURE_RULES.md`.** `docs/STYLE.md` is prose only and
lists figure legend style under "9. Not yet covered". There is no human-section figure
layer in this repo, so nothing was set there beforehand -- the human panels are built
outside it.

**The construct palette already covers most of Fig. 3** (FIGURE_RULES.md, "The cell-line
constructs"; `construct_cols` in `figures/natmetab/_style.R`), set 2026-10-19 and extended
2026-10-23:

| Fig. 3i arm | colour |
|---|---|
| `BclxL` (zero arm) | black `#000000` |
| `Pgc1a_BclxL` (effector blocked) | dark green `#00573F` |
| `NTPgc1a` (pathway intact) | bluish green `#009E73` |
| `KOPgc1a` (checkpoint removed) | **no colour set -- see below** |

**The gap, and it is a real one.** `NTPgc1a` and `KOPgc1a` are the same construct
(PGC-1a) on different backgrounds, so under the current palette both take `#009E73` --
and FIGURE_RULES says two constructs share a hue only where they share a mechanism and
**never in one panel**. Fig. 3 distinguishes arms by **genetic background**, a dimension
the Fig. 2 palette does not encode at all. Three ways out, for the author to pick:

- give the background its own encoding (fill for construct, outline or hatch for
  background), which keeps one meaning per hue;
- treat the three routes as three categories and give the p19ARF-disrupted arm a fourth
  colour, which costs the construct-hue rule;
- drop colour from Fig. 3i entirely and let the route labels carry it, since the panel is
  already ordered by route and the axis names each arm.

The same gap applies to Extended Data a, which plots the same three arms.

Also still open from 2026-10-23: **EV-EV's colour**, and whether the no-Dox tints stay
plain or take a hatch (`outputs/natmetab/proposals/p23_nodox_tints.pdf`, pages 1 and 2).

---

# Handoff

**State.** Script 57 is written and has been run twice; the second run carries the display
fixes. Content is verified. **Nothing is committed.**

Uncommitted, all of it:

| file | state |
|---|---|
| `scripts/57_section3_panel_readout.R` | new |
| `data/orthotopic_series/README.md` | modified (pooled CRISPR correction) |
| `docs/2026-09-23_section3_narrative_and_figure3.md` | new (from Project knowledge) |
| `docs/2026-10-23_MYAZ_proteomics_mass_redox_and_bcl2_family.md` | new (from Project knowledge) |
| `docs/2026-10-27_section3_final_narrative_and_panel_spec.md` | new (the spec) |
| `docs/STYLE.md` | new (from Project knowledge) |
| `docs/2026-10-28_section3_panel_readout.md` | new, auto-generated by 57 |
| `docs/2026-10-28_section3_results_digest.md` | new (this file) |

`outputs/orthotopic/57_*` (4 CSVs, 2 PDFs) are gitignored, as with script 53.

**To paste into Prism now:** `outputs/orthotopic/57_fig3g_ppargc1a_prism.csv` -- two
columns, seven rows, unpaired Mann-Whitney.

**To pick the work back up:**

1. Read this file, then `docs/2026-10-28_section3_panel_readout.md` for the full record.
2. The panels need formatting against `figures/natmetab/FIGURE_RULES.md`, and the colour
   gap in section 6 needs a decision first.
3. Re-run 57 if the manifest lands; nothing else needs re-running.

**Do not:** re-fit anything script 52 produced; quote across the 21-sample and 49-sample
scoring runs; write a causal verb for the retained-dose relation (every arm is a survivor
population); or use the word "primed" of a transcript.
