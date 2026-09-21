---
date: 2026-09-21
tags: [project/myc_mouse, experimental-cohorts, figure-1, verification, pre-declared-rules]
status: record -- every number reconciles against results/two_timeline_verification.rds (script 54, author-run 2026-09-21 11:53 UTC) or the saved object named beside it
relates-to:
  - scripts/54_two_timeline_verification.R
  - figures/panels/PANELS.md (section "The Figure-1 verification pass")
  - docs/handoff.md
  - docs/2026-09-02_myc_oxphos_priming_gate_model.md
---

# Figure 1 verification: three checks on the two-timeline plane, and the coupling on one scale

The closing paragraph of the Figure 1 block rules out three alternative explanations for MYC
no longer killing the adult gland and leaves one hypothesis standing. Before the text was
finalised, its numbers were checked against **rules fixed before any number was retrieved**
(2026-09-21) and **not adjusted afterwards**. This note records every outcome against its rule.
The results that go against the draft come first.

---

## 1. Findings against the pre-declared rules

### 1.1 CHECK 1 FAILED: a negative result against a pre-declared threshold

**The rule, fixed before retrieval.**
- A gene is *on the diagonal* if |interaction log2FC| < 0.20 **and** its raw interaction p > 0.05.
- Check 1: if Bax and Bcl2l1 are both on the diagonal, the text stands. If either is off it with Bbc3's sign, "only the PUMA to BCL-XL balance behaved differently" is false.

**The outcome.**

| gene | interaction (SE) | against the rule | verdict |
|---|---|---|---|
| Bcl2l1 | +0.179 (0.213), raw p 0.40 | inside 0.20, p > 0.05 | **on the diagonal** |
| Bax | **-0.231** (0.150) | **beyond 0.20 by 0.031, with Bbc3's sign** (Bbc3 -0.541) | **off the diagonal** |

**"Only the PUMA to BCL-XL balance behaved differently" is WITHDRAWN** (author, 2026-09-21).
Bax is **MYC-specific on the plane by the declared rule**: its interaction is negative and
exceeds 0.20. The author accepted this against the threshold as declared.

Bax's raw interaction p (0.123) is recorded here as part of the record. It is **not for the
manuscript**, because Bax is exploratory (ruling 2, section 3).

### 1.2 Fig. 2G: no "off the line" criterion was declared, so the residual is reported without a verdict

Ruling 4 replaced the imported retention rate with a line fitted **inside** the ratio set: a
fit through the origin to the eight pro:anti ratios that are not PUMA:Bcl-xL.

| | value |
|---|---|
| slope over the other eight | **0.429** (SE 0.137); residual SD 0.221 on 7 df |
| R2 (squared Pearson, the Fig. 1G definition) | 0.43 (the fit's uncentred R2, 0.58, is not quoted) |
| PUMA:Bcl-xL, six-week effect | +0.674 |
| predicted twelve-week effect / observed | +0.289 / **-0.061** |
| **residual** | **-0.350 = -1.58 residual SDs**, the largest of the nine |
| 95% prediction interval at its six-week effect | [-0.277, +0.855]: **PUMA:Bcl-xL lies inside it** |
| the next three residuals | Pmaip1:Bcl2l1 +0.308, Bmf:Bcl2l1 -0.307, Bbc3:Mcl1 -0.260 |
| sensitivity, the seven that are not a PUMA ratio (reported, not drawn) | slope 0.456; residual -0.368 (-1.73 SD); inside [-0.258, +0.873] |

**No criterion for "off the line" was declared before the numbers were seen, so none is
applied now.**
- PUMA:Bcl-xL is the **extreme of nine**, and it lies **inside** the interval.
- The earlier legend's "deviating significantly from the expected pattern" is not what this line says.

Source: `$ratio_line`, `$ratio_resid`, `$ratio_band`.

### 1.3 Bax's 0.5495 was computed, not compared

- **The value.** Bax:Bcl-xL retention = d12 / d6 = 0.495411 / 0.901492 = **0.5495**. d6 and d12 are script 42's genotype coefficients from per-animal OLS on log2(pro) - log2(anti) (`priming_arm_teb.rds$priming`; retention is set at `42:231`). Script 54 asserts saved = recomputed to 1e-12 for all nine ratios (`$retention_provenance`).
- **Where the constant enters.** `GLOBAL_RATE = 0.55` is hard-coded at `42:61` and `44:89`. In script 42 it enters only as:
  - the +/-0.15 label bands of `machinery$vs_global`;
  - the target of the matched-pair null control (null median retention 0.551).

  It is never the value of a retention.
- **So the number is not circular.** The claim "fades at exactly the programme-wide rate" compared a ratio-level OLS retention with a rounded pathway-level slope (script 40's 0.552 on the content ruler). That is two estimators, and ruling 4 retires the comparison.
- **Flagged, not removed:** the hard-coded 0.55 in scripts 42 and 44. This is the commit-97348a8 problem: a number copied into code.

### 1.4 Check 2 passes on the amended rule

- **Ruling 1** replaced a rule keyed on the sign of y. A gene is *MYC-specific* if its interaction is negative and exceeds 0.20, whatever the sign of either arm.
- **Foxo3:** wild type +0.473 (SE 0.128), Myc+ **-0.013** (SE 0.128), interaction **-0.486** (SE 0.181). MYC-specific.
- **The sentence** is "**rose in the normal gland and did not rise under MYC**", not "fell under MYC": the Myc+ arm is flat.
- **The pair that both miss 0.05** is the *genotype* pair: 6W +0.243 (p 0.057), 12W -0.242 (p 0.059). The significance prohibition binds those, not the temporal arms (ruling 1).
- **Foxo3 is exploratory:** it was not pre-specified and was found in script 44's scan. Its p-values are in section 2 and are not quoted in the manuscript (ruling 2).

### 1.5 Check 3 passes

The six-week MYC effect across **all nine** ratios:
- n = 9; median +0.487; Q1 +0.243; Q3 +0.688.
- **IQR 0.444 >= 0.30**, so the panel stays a scatter.
- Bax:Bcl-xL is at +0.901 and PUMA:Bcl-xL at +0.674.

The two ratios the 2026-08-05 panel did not draw were **Bax:Mcl1 and Bbc3:Mcl1**. The reason
given was that "a second denominator needs a second encoding" (`fig2_priming_ratios.R:75-77`).
Ruling 3: all nine.

### 1.6 The coupling: the claim is the interaction (ruling 7)

**One scale.** Both variables are z-scored over the 24 animals. Slopes are SD of the PUMA:Bcl-xL
log ratio per SD of the OXPHOS mitoPPS score, with n = 12 per genotype (6 per age).

| fit | wild type | Myc+ | interaction |
|---|---|---|---|
| unadjusted | -0.61 (p 0.053) | +0.56 (p 0.042) | **+1.17 (SE 0.37), p 0.0046**, 20 df |
| pre-specified, epi + imm, covariates shared | +0.05 (p 0.82) | +0.83 (p 0.0001) | **+0.78 (SE 0.25), p 0.0052**, 18 df |
| epi + imm within each genotype (8 df each) | +0.04 (p 0.88) | +0.83 (p 0.0036) | -- |
| within-timepoint permutation, 5,000 draws | | | unadjusted: percentile 97.4 (p 0.026); adjusted: percentile 91.5 (p 0.085; script 43 had 91.7) |

- **Withdrawn:** "respiratory priority tracked this ratio only where MYC was present". It is false on the unadjusted fit, where the slopes are opposite and nearly equal. Ruling 7's replacement ended in a robustness clause, "and the difference did not depend on adjustment for epithelial and immune composition". That clause is now withdrawn as well (next paragraph).
- **Withdrawn:** the redox sentence. Its replacement, legend-only, is that the two redox slopes do not differ: interaction +0.11 (p 0.84) unadjusted, +0.21 (p 0.72) adjusted. As drawn, redox tracked the ratio in wild type (p 0.0003) and not detectably in Myc+ (p 0.33); after adjustment, in neither.
- **The three numbers that were never comparable.**
  - -4.80: SD of the VST ratio per *unit* of mitoPPS, unadjusted (the old Fig. 2I).
  - -2.19: the same fit on the unscaled log2 ratio (script 48 `simple_slopes`).
  - -0.01: a *different predictor* (`ox_rel`), unscaled response, adjusted.

  The "factor of five hundred" was units and predictor. On one scale the adjustment moves **both** slopes.
- **Why the wild-type slope moves.** Of the -0.65 change, -0.57 runs through the immune composite and -0.08 through the epithelial one (an exact identity). Within the wild-type animals, the ratio correlates +0.90 with the immune composite and the OXPHOS score correlates -0.63 with it.
- **Fragility, which the legend carries:**
  - the covariates are RNA surrogates from the same count matrix as the ratio;
  - they correlate -0.84 (wild type) and -0.92 (Myc+);
  - within-genotype fits have 8 residual df;
  - **no positive control was carried through the adjustment.** The rule that now requires one (handoff section 5, item 7) postdates script 48 by five days.
- **When the adjustment was specified.**
  - 2026-07-25, abce8c8 / 5503aeb: the covariates were fixed with the model's first version. No unadjusted version of the model was ever committed.
  - 2026-08-05, b4f2a93: the first per-genotype slopes on record, unadjusted.
  - 2026-09-02, 3effad1: the adjusted per-genotype slopes were first computed, as a correction of the unadjusted ones.

  **Neither per-genotype fit was pre-specified.** The panel now draws both.

**THE ROBUSTNESS CLAUSE IS WITHDRAWN: an impasse at n = 24, not a choice** (author, 2026-09-21, second round).
- Against the within-timepoint permutation null, the **unadjusted** interaction clears (percentile 97.4, p 0.026) and the **adjusted** one does not (percentile 91.5, p 0.085).
- **So the fit that clears is the contaminated one.** The fit exposed to composition clears; the fit that removes composition does not.
- This is recorded as an **impasse at n = 24, not a choice** between the two fits. Neither is preferred.
- **In Figure 1 the coupling is hypothesis-generating and asserts nothing.**
- Both permutation results are in the Fig. 2I legend, which now states the impasse; 2H+I (alt) part B says the same. Both panels assert the pattern (unadjusted p < 0.05, adjusted p >= 0.05), so a re-run that changes it stops the panel instead of leaving the sentence behind.
- Each half of the drawn panel prints its parametric and its permutation p side by side (author, third round): 0.0046 and 0.026 unadjusted, 0.0052 and 0.085 adjusted. The parametric p alone is anti-conservative at this n, and the permutation p alone hides the discrepancy. Together they are the impasse, shown on the page.

### 1.7 The arms on both rulers (ruling 5)

The 0.20 magnitude is applied as declared, with no arm-specific threshold, on the content ruler
(set-average log2FC, unweighted) and on mitoPPS.

| arm | content: WT / Myc+ | content interaction | matched-null percentile of the interaction | mitoPPS interaction |
|---|---|---|---|---|
| **OXPHOS subunits** | -0.255 / -0.405 | **-0.150** (inside 0.20 by 0.050) | **1.2** (null median -0.058; 95% -0.133 to +0.021) | **-0.044** |
| OXPHOS assembly | +0.001 / -0.175 | -0.176 | 0.6 | -0.027 |
| nucleotide metabolism | +0.000 / -0.183 | -0.183 | 2.45 | -0.004 |
| mitoribosome | -0.024 / -0.252 | **-0.228 (beyond 0.20)** | 0.0 | -0.048 |
| TCA cycle | +0.041 / -0.240 | **-0.280 (beyond)** | 0.1 | -0.040 |
| amino-acid metabolism | +0.184 / -0.057 | **-0.241 (beyond)** | 0.1 | -0.018 |
| lipid metabolism | +0.122 / -0.028 | -0.150 | 0.65 | +0.000 |
| proliferation (PROLIF_* pooled) | -0.047 / -0.149 | -0.102 | 0.0 | -- (not MitoCarta) |
| TEB vs ductal (HS) | -0.422 / -0.220 | **+0.202 (beyond, above)** | 100 | -- |

- **"OXPHOS on the diagonal" holds by the declared magnitude on both rulers.** On the content ruler, however, the respiratory chain's drop is **beyond comparably expressed genes**, at percentile 1.2 of its matched null. Both are reported (ruling 5): the content ruler sits nearer the declared boundary (margin 0.050, against 0.156 on mitoPPS).
- **On the content ruler, all seven mitochondrial arms lie outside their own null bars**, below the diagonal. The three beyond 0.20 are the three MYC raised most at six weeks: +0.61, +0.57 and +0.56 (asserted in the panel).
- **On mitoPPS, every mitochondrial arm lies within 0.048 of the diagonal.**

**WITHDRAWN: "MYC neither causes nor prevents the respiratory withdrawal"** (author, 2026-09-21, second round). **Replacement: most of the withdrawal occurs on both timelines.**
- On the content ruler the respiratory chain (OXPHOS subunits) falls -0.255 across the wild-type window and -0.405 across the Myc+ one. So **63% of the Myc+ fall (0.630) is also on the wild-type timeline.**
- **All seven mitochondrial arms exceed their matched nulls.** Each lies below the diagonal, beyond the 2.5th percentile of its 2,000 expression-matched sets (percentiles 0.0 to 2.45). The interactions run from -0.150 (lipid metabolism, with OXPHOS subunits also at -0.150) to -0.280 (TCA cycle).
- **So the MYC-specific component is compartment-wide rather than OXPHOS-specific.** The respiratory chain's own is among the smallest of the seven.
- What *is* specific to the respiratory chain is the wild-type withdrawal. It is the only mitochondrial arm that falls across the wild-type window (percentile 0.0 of its matched null); the other six run from -0.024 to +0.184.
- **Wild-type-axis percentiles, redrawn by script 54:** respiratory chain 0.0 and proliferation 1.7, against script 43's 0.00 and 1.30. The largest difference over the ten arms is 3.2 points, at nucleotide metabolism. **The text quotes script 43's.**

Source: `$arms_content`, `$arm_null`, `$arms_mitopps`, `$arm_diagonal`.

---

## 2. The four genes (the Task B table)

Source: `results/interaction_results.rds` (script 03), via script 54 `$four_genes`.
- **Values:** raw (unshrunken) MLE log2 fold changes from `~ timepoint * myc_status`, with standard errors and **raw** Wald p.
- **The adjusted p is IHW** ("Weighted BH adjusted p-values", the column's own description). It is carried for completeness; no verdict uses it.
- **The identity** c_tp = c_tn + c_int holds to <= 5.6e-17 for these four genes, and to 1.8e-15 across all 18,523.

| gene (Ensembl) | baseMean | wild type 6>12W: LFC (SE), raw p, IHW padj | Myc+ 6>12W: LFC (SE), raw p, IHW padj | interaction: LFC (SE), raw p, IHW padj | verdict | licence |
|---|---|---|---|---|---|---|
| Bbc3 (ENSMUSG00000002083) | 152 | +0.062 (0.144), 0.668, 0.881 | -0.479 (0.145), 0.00094, 0.016 | **-0.541 (0.204), 0.0081, 0.844** | MYC-specific | pre-specified |
| Foxo3 (ENSMUSG00000048756) | 3265 | +0.473 (0.128), 0.00022, 0.0078 | -0.013 (0.128), 0.922, 0.977 | **-0.486 (0.181), 0.0074, 1.000** | MYC-specific | exploratory: position only in the manuscript |
| Bax (ENSMUSG00000003873) | 481 | -0.145 (0.107), 0.173, 0.489 | -0.376 (0.105), 0.00033, 0.0077 | **-0.231 (0.150), 0.123, 1.000** | **MYC-specific: Check 1 fails** | exploratory: position only in the manuscript |
| Bcl2l1 (ENSMUSG00000007659) | 870 | -0.107 (0.151), 0.479, 0.790 | +0.073 (0.151), 0.631, 0.834 | **+0.179 (0.213), 0.401, 1.000** | on the diagonal | denominator of the pre-specified pair |

---

## 3. Rulings of record, 2026-09-21

1. **Check 2 keys on the interaction, not the sign of y.** The significance prohibition binds the genotype contrasts. The wild-type temporal arm may be marked.
2. **Raw p at the declared threshold for all four genes, with the licence stated gene by gene.** Bax and Foxo3 are exploratory: their positions are plotted, and no p-value for either is quoted in the manuscript.
3. **Check 3 runs on all nine ratios, not seven.**
4. **Fig. 2G: the line is fitted to the other eight ratios.** No imported rate, and no hard-coded 0.55.
5. **Arms:** apply 0.20, report both rulers, draw mitoPPS, and say that the content ruler sits nearer the boundary.
6. **Keep IHW and label it correctly.**
7. **The coupling claim is the interaction.** Both fits are drawn. "Only where MYC was present" and the redox sentence are withdrawn.
8. **Re-source 32, 43 and 44 after Task D's asserts.**

Afterwards:
- **`ortholog_table.rds`** stays as it is (section 4).
- **Check 1's failure is accepted.**
- **Fig. 2G is built as ruling 4 specifies.**
- **The four-gene panel is a three-way contrast.**
- **The open item in section 5 is recorded and not pursued.**

Second round, the same day:
1. **No p-value for Foxo3 or Bax anywhere a reader could read it as a test.** Neither was pre-specified. 2H+I (alt) part A keeps Foxo3's brackets and loses their p-values. **Over-extended: partly reversed in the third round, below.**
2. **The flags on `combined_df_*` and on the outputs of scripts 11 and 16 are false positives of the timestamp rule** (section 4). Render the panels. The headers of 2H and S2C credited script 03 for `combined_df_annotated_raw.rds`; 2H's temporal fold changes are repointed to `interaction_results.rds`, which script 54 read.
3. **"Priming" is cleared from the six legend blocks.** `priming_arm_teb.rds` stays as a filename and is no longer printed in captions.
4. **Commit in five groups:** script 54 and this note; the Fig. 2I rebuild, the `_panel_common.R` changes and 2H+I (alt) part B; the three `plane_` panels; the IHW relabels; the priming clearance.
5. **Two text changes:** "MYC neither causes nor prevents the respiratory withdrawal" is withdrawn (section 1.7), and so is the coupling's robustness clause (section 1.6).
6. **Fig. 2G is built next session** from script 54's object.

Third round, the same day: **a correction** (author).
- **Second-round ruling 1 over-extended ruling 2, and is partly reversed.** Ruling 2 was about the manuscript TEXT, not the panels.
  - **Pre-specification governs what may be CLAIMED, not what may be SHOWN.** An exploratory p-value, reported and labelled as exploratory, is normal practice.
  - A concealed one is worse than either alternative: the reader cannot calibrate the observation, and it reads as selective reporting.
  - The incoherence showed in the second round's own output: a sentence that used Foxo3's p-value while declining to print it.
- **Reversed:**
  - **Foxo3's interaction p** is back in the legends with its status attached, in the author's form: interaction -0.486, nominal p 0.0074; not pre-specified; does not survive correction across the 26 (Bonferroni 0.19); reported as exploratory. The ranking sentence has its anchor again.
  - **The symbol changes are undone.** Bax is an ordinary symbol inside S2C's key, and Foxo3's brackets in 2H+I (alt) are back in normal ink with their p-values.
  - **2H** carries Foxo3's interaction p and the two six-week adjusted p-values again. Those two are now computed rather than typed.
  - **D4** carries Foxo3's ring and padj again.
  - **S2D (alt)** carries Foxo3's padj on its label again.
- **Stays removed, for other reasons:**
  - **2G (alt):** Bax's dark ring, red label and padj. They marked significance on the Myc+ ARM, which is not the contrast any claim rests on. That is unrelated to pre-specification.
  - **S2D (alt):** the "reference genes" key. A class defined by a threshold was the problem there.
- **The manuscript text still carries no p-value for Foxo3 or Bax.** That part of ruling 2 stands.
- **Reversed on the same principle, though not named in the correction.** `plane_four_genes.R` (first round) withheld Bax's and Foxo3's p-values from its legend, although that legend says the rule's p half "is in the detail above". It now prints each gene's interaction p (Bax 0.123, Foxo3 0.0074), labelled exploratory. Bax's Myc+-arm p is not printed, on the 2G (alt) ground.
- **One word of the author's form is changed:** "has the largest genotype difference" became "separates the genotypes more sharply than any".
  - By magnitude, Foxo3 is fifth of the 26: Chac1 +0.805, Htra2 -0.546, Bbc3 -0.541 and Endog -0.533 are larger.
  - It is first by p, and by |z| (2.68, against Bbc3's 2.65), which is what "more sharply" states.
- **The two open questions, answered in the same message:**
  - **On-panel p-values:** print BOTH the parametric and the permutation values, side by side (section 1.6). Fig. 2I and 2H+I (alt) part B now read "p 0.0046 parametric, 0.026 permutation" (unadjusted) and "p 0.0052 parametric, 0.085 permutation" (adjusted). The label sits in a band of its own above the data.
  - **The N3 residue: leave both.** A filename is a filename, and `biogax_rescue_licence.R` quotes a claim in order to examine it, not to assert it.
- **Housekeeping:** the sixth commit stays where it is, the reversal is a commit of its own, and the branch is pushed. Fig. 2G is next session, from script 54's object, including its "priming" clearance.

---

## 4. Provenance and the freshness rule

**Script 54.**
- `results/two_timeline_verification.rds` was written 2026-09-21 11:53:13 UTC, after the script's last modification (11:47:48 UTC). The script is uncommitted.
- The first run stopped at a names-only check. `as.data.frame()` drops a tibble column's names once tibble is attached, while script 48's saved vector keeps them. The check now compares values and asserts sample order explicitly. The rerun is the object of record.
- Controls that passed, each against an existing record:
  - the plane's identity (1.8e-15);
  - two independent symbol routes, and script 42's tables (1e-12);
  - script 43's arm values (1e-9) and its mitoPPS read (1e-12);
  - script 43's interaction (6.0885, p 0.00523, to 1e-9);
  - script 48's simple slopes (1e-9);
  - the old Fig. 2I numbers (-4.798, +4.347);
  - the wild-type null percentiles (within 3.2 points);
  - script 43's permutation percentile (91.5 against 91.7).

**Objects read, all fresh by the rule** (object modification time after the last commit touching its writer):
- `interaction_results.rds` (03)
- `priming_arm_teb.rds` (42)
- `substrate_specificity_tradeoff.rds` (43)
- `collapse_module_ownership.rds` (44)
- `background_vs_myc.rds` (40)
- `gate_model_verification.rds` (48)
- `gsva_scores.rds` (15)

**Flagged by the rule: three false positives (author's ruling, 2026-09-21, second round).**

On the evidence below, each flag comes from how the rule reads timestamps. None is a stale object.

| object | writer | why the rule flags it | why it is a false positive |
|---|---|---|---|
| `combined_df_annotated_raw.rds`, `combined_df_annotated.rds` | the archived main pipeline: `archive_main_pipeline/02_deseq_interaction_model.R`, last re-saved by its `04_group_comparison.R` (the saved table carries 04's `group_*` columns) | the writers' last commit (94d346c, 2026-03-13) postdates the objects (2026-02-07) | 94d346c is a byte-identical rename into `archive_main_pipeline/`. The last content change is 2026-01-26 (7368971), before the objects |
| `interaction_gene_characterisation.rds` | script 11 | written 16 minutes before its script's commit | the script file's own modification time precedes the object: edit, run, commit |
| `cell_death_de_raw.rds` | script 16 | written 21 minutes before its script's commit | as above |

- **Rendered.** The five panels that read the two combined tables (2H, 2H+I (alt), S2C, S1E and S1F) are now rendered. No panel reads the outputs of scripts 11 and 16.
- **A real defect found beside them, and fixed.**
  - 2H and S2C read `combined_df_annotated_raw.rds`, but their headers credited script 03. 2H also took Foxo3's and Bbc3's temporal fold changes from it.
  - 2H now takes those from `interaction_results.rds` (script 03), the object script 54 read. It keys them through script 44's `$collapse_genes$ens` and no longer reads the combined table.
  - The values are identical: Foxo3 +0.473164 / -0.012519 and Bbc3 +0.061892 / -0.478699 from either object, and the same as `$four_genes`.
  - S2C still takes its symbol mapping from the combined table. S2C's header and 2H+I (alt)'s header, which did not list the table at all, now credit the archived main pipeline.

**The `ortholog_table.rds` exception (author's ruling, 2026-09-21).**
- **Why the rule flags it:** it was written 2026-02-07, before `01_load_data.R`'s commits of 2026-02-14 (9f63f29) and 2026-02-16 (41af03c).
- **Neither commit changed the biomaRt query.** 9f63f29 moved the same query behind a cache; 41af03c added Hallmark loading.
- **It is the pinned mapping** every downstream object on disk was built with.
- **Regenerating it would re-query live biomaRt.** A different Ensembl release would shift gene-set membership across scripts 15-53.
- **So it is left as it is.** Script 54 reads it through the reconciler and proves its arms against script 43 to 1e-9.

---

## 5. Open item, recorded and not pursued (author, 2026-09-21)

**The unadjusted Myc+ coupling slope is largely timepoint.**
- Adjusted for timepoint alone it falls from +0.563 to +0.264 (p 0.38); with timepoint + epi + imm it is +0.621 (p 0.043).
- The wild-type slope stays negative with timepoint alone (-0.660, p 0.080).
- The composites are not timepoint proxies: within a genotype, timepoint explains at most 4.8% of either.

This bears on whether the coupling is within-animal variation or a group difference.

Source: `$coupling_timepoint`.

---

## 6. What changed in the figure layer, and what is pending

`figures/panels/PANELS.md` has the full list. In short:
- **Fig. 2I rebuilt:** two fits, one builder.
- **Fig. 2H+I (alt) part B rebuilt** on the same builder. Rendered in the second round; part B fitted the 100 mm layout until both p-values were printed, and the panel is 112 mm since the third round.
- **Three new unslotted panels:** `plane_arms_content`, `plane_arms_mitopps`, `plane_four_genes`.
- **IHW relabels** in five legends.
- **Three helpers** added to `_panel_common.R`.

**Second round, the same day:**
- **No test for Foxo3 or Bax, anywhere.** *Reversed in the third round (section 3), except 2G (alt) and S2D (alt)'s key.* Every one of these was removed:
  - **2H+I (alt):** part A keeps Foxo3's three brackets, unlabelled, in a neutral ink outside the sig/ns pair, because an "ns" grey would be read as a test too. The legend gives Foxo3 by position, and no p-value from the mechanism-gene ranking is printed.
  - **2H:** Foxo3's interaction p is gone from the legend, and so are the typed "0.22 and 0.21". The ranking-set fact is asserted instead.
  - **2G (alt):** Bax has no ring, its name is in the plain ink, and its padj is gone. The count of significant transcripts is now computed. The typed "three" had been four with Bax (Bax, Bbc3, Bmf, Htra2); it is three without Bax.
  - **S2C:** Bax is drawn as a cross outside the key, and its padj is gone.
  - **S2D (alt):** Foxo3's on-panel "padj 0.0078" is gone, and the key's "moves (padj < 0.05)" now reads "reference genes". Bbc3 keeps its padj.
  - **Discussion D4** (`biogax_factor_roster.R`): Foxo3's red ring and its padj are gone.
- **The provenance defect** (section 4) is fixed in 2H, S2C and 2H+I (alt).
- **"Priming" is cleared** from the legend blocks of 2E, 2G (alt), S2C, S2D and S2D (alt), and from 2H+I (alt)'s source line.
  - Captions now cite "Script 42's saved object". They no longer print `priming_arm_teb.rds` or script 42's own filename, which also contains the word.
  - The sixth block, Fig. 2G's, is cleared in its rebuild.
- **Fig. 2I and 2H+I (alt) part B state the impasse** (section 1.6) and no longer say the difference persists.
- **Rendered and inspected:** 2H, 2H+I (alt), S2C, S1E, S1F, 2I, 2G (alt), 2E, S2D, S2D (alt) and D4.

**Third round, the same day (the correction, section 3):**
- **Reversed** in 2H+I (alt), 2H, S2C, S2D (alt), D4 and `plane_four_genes.R`: Foxo3's and Bax's p-values are shown and labelled exploratory, and the symbols are back to the panels' own grammar.
- **Kept, with the reason corrected:** 2G (alt)'s unmarked Bax, and S2D (alt)'s "reference genes" key.
- **Rendered and inspected:** 2H, 2H+I (alt), S2C, S2D (alt), 2G (alt), D4 and `plane_four_genes`.

**Pending:**
- **Fig. 2G's ruling-4 rebuild,** next session from script 54's object (ruling 6). The committed version still draws the retired 0.487 line, and its legend block still says "priming".
- **Task D is written (section 7) and waits on the re-source of scripts 32, 43 and 44.** Those re-runs are the first time the asserts execute.
- ~~The on-panel p in Fig. 2I and 2H+I (alt) part B.~~ Answered in the third round: both p-values, side by side.
- ~~N3 residue outside the six blocks~~ (the two file slugs in the `legends.md` headings; "respiratory capacity" quoted by `biogax_rescue_licence.R`). Answered in the third round: leave both.
- **Documents that still call IHW "BH"** (narrative v3 "Bbc3 BH 0.84"; several PANELS.md entries). They are left as written and corrected here.
- **`paper/analysis_record.qmd`** is not updated, by instruction.

---

## 7. Task D: the paragraph's numbers behind asserts

Task D applies the principle of commit 97348a8: a number copied between documents is never checked, and a number in a stopifnot is checked every run.
- **Where the numbers live.** Each writer script now DECLARES the numbers the Figure-1 closing paragraph quotes, in a `PART 0 (TEXT)` block before anything is computed, and ASSERTS them where they are computed.
- **Precision.** Each is checked at the precision the paragraph quotes it.
- **Pre-checked.** Every assert was evaluated against the saved objects before commit and holds. All three scripts parse.
- **The first run.** The asserts first execute when you re-source the scripts.

**Covered**

| the paragraph quotes | script | asserted where | value now | checked to |
|---|---|---|---|---|
| **27** (upper end of "21-27%") | 32 | PART 3, after `share_stats` | +26.54%: MYC's effect on the chaperone-free mass-marker share | the integer percent |
| **0th** (respiratory chain) and **1.3rd** (pooled proliferation) | 43 | PART A, after `wt_null` | 0.00 and 1.30, of 2,000 matched sets (seeded) | one decimal |
| **0.487** (the retention rate) | 44 | PART C, after `collapse_genes`, and again before the save | 0.48723, through the origin over the 2,648 genes MYC moves | three decimals |

**The hard-coded 0.55 in script 44** is asserted to be a reference and never a computed value:
- **The fitted rate is recomputed in closed form**, so replacing it with the constant would fail. The residual and the retention-minus-rate columns are re-derived against the fitted rate.
- **The saved `defs` holds 0.55 in exactly one field,** `global_rate_assumed`.
- **The notes no longer type the constant.** Script 44's notes printed "script 40's value 0.55" as typed text; they now format it from `GLOBAL_RATE` and call it a labelled reference.

**Not covered**
- **The lower end, 21.**
  - Script 33 writes it, not script 32: `genotype_untouched$adj_pct` = 21.12%, the mass-marker effect adjusted for prep-stress and contamination.
  - Covering it means an assert in script 33 and a fourth re-source. That is the author's call.
- **The coupling p-values**, by ruling. The permutation impasse leaves no single fit to assert against.
- **Anything withdrawn today**, by ruling: the "only PUMA" claim, the coupling's robustness clause and the shared-withdrawal phrasing.

**The repository-wide audit of 0.55**
- **Script 42** (`GLOBAL_RATE`, line 61) uses it only as a reference:
  - the +/-0.15 bands of `machinery$vs_global`;
  - the target of the matched-pair null control;
  - the saved `params$GLOBAL_RATE`.

  Its prose also types it three times (lines 644, 685 and 712; the last reads "matched pairs retain 0.55, script 40's global rate exactly", a typed copy of a computed 0.551). Its nine retentions are computed: script 54 checked each against d12/d6 to 1e-12. Script 42 is not in Task D's re-source set, so it is not edited.
- **Found: the constant standing in for a computed rate.**
  - `figures/fig05_death_arm.R` (line 58) and the assembled `figures/figure2_developmental_window.R` (line 73) read `pa$params$GLOBAL_RATE` and draw it as the rescaling line, labelled "global rate x0.55" and "x0.55".
  - `fig05` also uses `machinery$vs_global`, which is banded on the constant.
  - Both are in the older exploratory and assembled layer, not `figures/panels/`. They are not changed: they should read `defs$global_rate_fitted` or be retired, and that is the author's call.
- **The panel layer** uses no `GLOBAL_RATE`. The current Fig. 2G reads the fitted 0.487, and its rebuild drops any imported rate (ruling 4).
- **The Results text** in the analysis record says "the overall 0.55-fold transcriptomic rescaling". The record's own callout says the sentence should quote 0.487 ("the genes MYC moves") or 0.450 (whole transcriptome). If the paragraph still says 0.55, it is quoting a rounded copy of the mitochondrial slope (0.552), not script 44's rate.

**Re-source:** `scripts/32_mito_content_proxies.R`, `scripts/43_substrate_specificity_and_tradeoff.R` and `scripts/44_collapse_module_and_ownership.R`.
- **Order:** any, since none reads another's output.
- **Why:** the commit that adds the asserts postdates each object, so all three are stale by the freshness rule until re-sourced.
- **What to watch for:** each prints a `PART 0 (TEXT)` line. A stop there means the paragraph's number and the analysis have drifted apart.
- **Script 54 need not be re-run** if the three reproduce. Script 43 is seeded and script 32 is deterministic. Script 44 is seeded, and its fitted rate is closed-form.
