# Why does the Myc effect attenuate? Four competing hypotheses

**Status:** hypothesis note. **Nothing is built** — no script 42, no change to
`paper/myc_mito.qmd`. The question is blocked on the pending **MYC Western blot**, whose result
re-ranks the field (see §6). This note complements, and does not supersede,
`docs/2026-07-24_background_vs_myc_interpretation.md` (script 40): that one establishes *what*
happens, this one sets out the candidate explanations for *why*.

Prompted by the author's challenge (2026-07-25) that the background's opposition to Myc might be
the cause of the attenuation rather than an incidental fact — a genuinely different model from the
"the programme simply decays" summary, and one that had not been tested.

---

## 1. What has to be explained

Five facts, all anchored. Nothing below may contradict them.

1. **The genotype gap halves, keeping its shape.** Regressing the 12W Myc effect on the 6W effect
   across 144 MitoPathways: slope **0.55** (bootstrap 0.51-0.60), R2 0.80, rho 0.84 on the content
   ruler; slope 0.64, R2 0.79, rho 0.93 on priority. Per-tier slopes **0.53-0.83** — uniform. The
   coherence sits at the 100th percentile of overlap-preserving gene-label shuffles.
2. **The driver does not fall.** Myc's genotype gap *widens*, +1.64 to +1.80 log2; Myc:MAX +1.42
   to +1.74; Myc:repressor-sum +1.65 to +1.90 — while the **output** (Hallmark MYC targets V2)
   narrows, +1.32 to +0.95. The proximal MYC/MAX/MXD network is flat (figS8).
3. **The attenuation is global and non-selective.** Script 34: no programme — death, OXPHOS or
   MYC-core — attenuates more than expression-matched background genes (all BH p > 0.14). Issue #4:
   transcriptome-wide median absolute LFC falls 0.268 to 0.204.
4. **The background opposes Myc exactly where the claim lives.** Projection of the WT 6-12W vector
   on the Myc(6W) vector = **0.064**, the 9.8th percentile of a matched null (i.e. *less* aligned
   than an arbitrary panel). Per tier, content ruler:

   | tier | Myc@6W | WT 6-12W | Myc+ 6-12W | interaction | pathways WT-down |
   |---|---|---|---|---|---|
   | OXPHOS | +0.46 | **-0.14** | -0.27 | -0.18 | 63% |
   | Protein import | +0.52 | **-0.06** | -0.28 | -0.22 | 75% |
   | Signaling | +0.20 | +0.02 | -0.08 | -0.07 | 43% |
   | Dynamics and surveillance | +0.24 | +0.03 | -0.07 | -0.12 | 38% |
   | Central dogma | +0.47 | +0.05 | -0.19 | -0.24 | 10% |
   | Small molecule transport | +0.23 | +0.07 | -0.04 | -0.10 | 20% |
   | Metabolism | +0.40 | +0.08 | -0.11 | -0.18 | 12% |

   Note the last two columns: the attenuation is near-constant (-0.18 to -0.24) whether the
   background opposes Myc (import, OXPHOS) or moves with it (central dogma, metabolism). What
   small differences exist run in the direction the algebra forces (§5), not the direction an
   antagonism model predicts, so they are **not** counted as evidence either way.
5. **The MYC stability machinery gives nothing.** Script 41: five genes down and two up on the
   Myc+ time axis, **none Myc-specific** (every interaction padj > 0.05), no excess over an
   expression-matched background, and the directions predict a *more* stable MYC (Fbxw7 -0.35,
   Huwe1 -0.24, both down).

---

## 2. The four hypotheses

### Two constraints that discipline all of them

- **An antagonist acting equally in both genotypes cancels in the genotype gap.** If a repressor
  lowers a gene by the same log amount in Myc+ and WT, the gap is unchanged. For any background
  mechanism to *cause* attenuation, its effect must scale with MYC — i.e. be non-additive. This is
  not a refutation, it is a specification: competition and sequestration satisfy it natively,
  because their effect is proportional to the activator present. A shift in set point does not.
- **A lower baseline normally leaves *more* fold-induction headroom, not less.** "The genes are
  already low" predicts a *larger* Myc effect at 12W unless the mechanism is active blockade.

### H1 — Autonomous decay

The Myc programme runs down on its own: substrate exhaustion, feedback saturation, or a cell state
that has finished the transition Myc was driving.

*Explains:* the uniform rescaling (fact 1), the globality (fact 3).
*Fails:* it names no cause. As written it is a description of fact 1 restated as a mechanism, and
it makes no prediction that any other hypothesis does not also make.
*Standing:* the default, and the weakest kind of default — unfalsifiable in its current form.

### H2 — Background antagonism (the author's model)

The genes Myc acts on are, by 12W, already held down by opposing regulators, and that opposition
blunts what Myc can add. Three sub-forms, which behave very differently:

- **H2a selective blockade** — specific induced repressors (NRs, corepressors) block specific Myc
  target promoters, the mitochondrial/biogenic ones in particular.
  *Predicts:* a selective dent — the blocked programmes attenuate more than the rest.
  *Fails:* fact 3 (no programme is in excess) and fact 1 (per-tier slopes uniform), plus the
  negative scan in §4.
- **H2b global competition** — occupancy competition at E-boxes, or limitation of a shared
  cofactor, reducing MYC's effective transactivation everywhere.
  *Predicts:* exactly the uniform x0.55 rescaling of fact 1.
  *Standing:* **the strongest form of the author's model.** It fits the uniformity better than H1
  because it supplies a cause for it. One named candidate exists (§4, the MLX arm), and its route
  is stoichiometric — invisible to RNA-seq by construction.
- **H2c lower set-point / ceiling** — no active blockade; the genes simply sit lower and Myc's
  absolute increment is smaller.
  *Fails:* the headroom constraint above. Also fails fact 4's tier pattern.

### H3 — Reduced effective MYC activity at constant message

MYC protein level, half-life, phosphorylation state or cofactor availability falls while the
transcript does not.

*Explains:* facts 1, 2 and 3 together, and more economically than anything else — a lower
effective MYC dose *is* a uniform multiplicative damper.
*Fails:* nothing yet, because it has not been measured. Fact 5 excludes only a **transcriptional**
reprogramming of the degradation machinery, which is a weak exclusion: constant Fbxw7 mRNA with
altered T58 phosphorylation would look exactly like what we see.
*Standing:* the front-runner on parsimony, and the one the pending blot addresses directly.
Note that **H2b and H3 converge** — both say effective MYC activity falls at constant message;
they differ only in whether the cause is competitive occupancy or the protein itself.

### H4 — Compositional dilution

The 12W gland contains a different mix of cells — more differentiated, fewer Myc-responsive
proliferative ones — so a bulk genotype contrast is diluted without any per-cell change in Myc's
effect.

*Explains:* facts 1 and 3 (dilution is multiplicative and therefore uniform and global — the same
signature as H2b/H3), and it now has direct marker support (§4).
*Fails:* it does not explain fact 2's *widening* driver gap without an auxiliary assumption, since
dilution should shrink the Myc message gap too. That is a real tension and worth testing.
*Standing:* under-weighted so far. Issue #6 flagged the E2F/proliferation co-decline as a hint of
compositional dilution and deferred deconvolution; the myoepithelial signature in §4 is the first
direct evidence, and it moves H4 up.

---

## 3. Hypothesis by observation

| observation | H1 decay | H2a selective | H2b competition | H3 effective MYC | H4 dilution |
|---|---|---|---|---|---|
| 1. uniform x0.55 rescaling | explains | **fails** | explains | explains | explains |
| 2. driver gap widens | silent | silent | explains | explains | **tension** |
| 3. global, no programme in excess | explains | **fails** | explains | explains | explains |
| 4. background opposes in import/OXPHOS | silent | needs it | needs it | silent | explains |
| 5. stability machinery flat | silent | silent | silent | **not excluded** | silent |
| names a testable cause | **no** | yes | yes | yes | yes |

---

## 4. What the read-only scans of 2026-07-25 already say

All values are DESeq2 **raw** log2FC with genome-wide BH padj, from
`results/interaction_results.rds`, filtered to baseMean >= 20. These were read-only explorations
(see §8) — they are not numbered pipeline scripts and carry no provenance beyond this note.

### 4.1 No E-box competitor moves — H2b has no transcriptional candidate

Scanned: Usf1, Usf2, Bhlhe40, Bhlhe41, Tfe3, Tfeb, Mitf, Max, Mnt, Mxd1, Mxd3, Mxd4, Mxi1, Mga,
Mlx, Mlxip, Mlxipl, Srebf1, Srebf2, Arnt, Hif1a.

**None is significant in the WT background and none has a significant interaction.** The closest
to a movement are Srebf1 (+0.61, padj 0.050) and Usf2 (+0.21, p 0.051), neither Myc-specific.
This reinforces figS8: the classical antagonist arm is unmoved, and the Myc:repressor ratio widens
rather than narrows.

### 4.2 The 111-gene repressor scan — H2a has no candidate either

Nuclear receptors, corepressors (NCoR/SMRT, Sin3, CoREST), Polycomb, repressive E2Fs and the DREAM
complex, circadian repressors, cell-cycle inhibitors and mammary differentiation TFs — 111 genes
resolved. **Six move in the WT background at FDR<0.05. Zero have a significant interaction.**

| gene | WT 6-12W | padj | Myc+ 6-12W | interaction padj | reading |
|---|---|---|---|---|---|
| Esrrb | +1.58 | 0.040 | +2.27 | 1.00 | an ESRR **activator** — wrong direction for repression |
| Id4 | +1.53 | 0.0013 | +1.27 | 1.00 | see 4.3 |
| Hes1 | +0.59 | 0.012 | +0.79 | 1.00 | bHLH repressor, but not a MYC competitor |
| Foxo3 | +0.47 | 0.0078 | -0.01 | 1.00 | growth restraint; shape would give an interaction if powered |
| E2f7 | -0.83 | 0.0028 | -0.85 | 1.00 | a **repressive** E2F going down — *less* repression |
| Ppard | -0.36 | 0.013 | -0.46 | 1.00 | — |

Two of the six point the wrong way for the model. None is Myc-specific.

### 4.3 The Id4/Hes1 lead is most likely composition — and it is evidence for H4

`Id4` was the single most attractive hit (an HLH sequestrator: a mechanism class whose effect
scales with the activator present, so it needs no interaction in its own expression). Two things
demote it:

- **Id1, Id2 and Id3 are flat** (-0.40, +0.23, -0.28, all ns). A regulatory shift in ID
  stoichiometry would not usually pick out one family member.
- **Id4 rises inside a myoepithelial/basal maturation signature** in the same background:

  | marker | WT 6-12W | padj |
  |---|---|---|
  | Oxtr | +1.72 | **0.00027** |
  | Trp63 | +1.17 | 0.092 |
  | Myh11 | +0.85 | 0.024 |
  | Krt5 | +0.80 | 0.31 |
  | Acta2 | +0.67 | 0.15 |

  Id4 is itself a basal/myoepithelial marker in mammary epithelium. The parsimonious reading is a
  **shift in cellular make-up**, not the induction of a MYC antagonist.

Mechanistically the lead was weak regardless: **ID proteins sequester E-proteins** (class A bHLH),
whereas MYC:MAX is bHLH-LZ and is not an ID target.

The positive consequence is the important one. **Oxtr at padj 2.7e-4, with Myh11 and Trp63
behind it, is the first direct marker evidence that the WT background move has a differentiation
character** — which is what H4 needs. (Caveat: bulk markers, n=6, on the batch-confounded axis.)

### 4.4 The one named H2b candidate: the MLX arm

The author's observation from figS8 panel B — `Mlxip` (MondoA) and `Mlxipl` (ChREBP) fall within
Myc+ across the window. This is the only movement anywhere in the MYC/MAX/MXD/MLX network.

| gene | baseMean | Myc@6W | Myc@12W | WT 6-12W | Myc+ 6-12W | interaction |
|---|---|---|---|---|---|---|
| Mlxip | 3187 | +0.07 (0.58) | -0.27 (0.061) | +0.15 | -0.19 (0.12) | -0.34 (0.83) |
| Mlxipl | 45 | +0.38 (0.80) | -0.26 (1.00) | -0.39 | -1.03 (0.037) | -0.64 (1.00) |
| Txnip | 3982 | -0.20 (0.62) | -0.18 (0.91) | -0.20 | -0.19 (0.75) | +0.02 (1.00) |

**Three problems.**

1. **The sign is wrong for competition.** MondoA and ChREBP are E-box competitors. Fewer of them
   means *less* competition, so MYC should hold up better at 12W, not worse.
2. **The interaction is null** (padj 0.83 and 1.00). figS8's brackets are *unadjusted*
   simple-effect p-values; on the genome-wide FDR that governs a claim, Mlxip's genotype effect at
   12W is padj 0.061 and its Myc+ temporal fall is padj 0.12. "Down in Myc+" is a within-genotype
   observation on the **batch-confounded** axis.
3. **Txnip is flat** (-0.20 at 6W, -0.18 at 12W, both ns). TXNIP is MondoA's most sensitive
   output; if MondoA activity were genuinely falling, TXNIP should register it. Arrdc4, the TXNIP
   paralogue, does fall within Myc+ (-0.72, padj 0.020) — but its interaction is null too.

`Mlxipl` is not usable on its own: baseMean 45, the WT falls as well, and mammary is not a ChREBP
tissue.

**The rescue that does work, and is worth recording.** MLX is a **shared** dimerisation partner:
it binds MondoA and ChREBP, and it also binds MXD1, MXD4 and MNT. Losing MondoA/ChREBP frees MLX
to form **MXD:MLX repressive complexes** at E-boxes. That is H2b — with **no repressor transcript
changing**, which is precisely why the flat MXD/MNT/MGA scan cannot refute it. It is a
**protein-stoichiometry hypothesis, so RNA-seq is structurally the wrong instrument**; only
co-immunoprecipitation or occupancy can read it (T3b, §6).

The alternative reading — MondoA/ChREBP as **co-activators** whose loss removes a second driver
from shared metabolic targets — does not survive. Those targets straddle the global 55% retention
rather than attenuating in excess: Ldha 85%, Elovl6 78%, Fasn 73%, Slc2a1 65%, Acly 51%,
Acaca 50%.

Finally, `Mlxip` is worth flagging as **the single most interesting gene in figS8 panel B** — the
only one whose shape (WT +0.15 against Myc+ -0.19) would yield an interaction if it were powered.

### 4.5 Issue #5 did not test this

Script 30 (`30_attenuation_moderation.R`) asked whether a Myc-independent axis absorbs the
attenuation and found that none does. But every axis it entered was **activating**: luminal
differentiation (`dev_luminal`), the ESRRA/GABPA/NRF1 biogenesis TFs (`tf_biogenesis`),
proliferation, and the mtDNA reprioritisation. **A repressive axis has never been entered.** The
author's model is not a re-run of Issue #5, and saying so would be wrong.

---

## 5. Two bars any future test must clear

**The natural test is algebraically forced.** "Genes where the background falls are the genes
where Myc is blunted" is `cor(WT-time, delta-Myc-effect)`. In script 40's artifact ledger the
observed value equals the value the marginals force, to the decimal: -0.185 content, -0.327
priority, because `delta = tpos - tneg` and `tneg` enters with a negative sign. Three further
statistics are in the same ledger. **None may be cited**, and any predictor built from the same
group means inherits the problem. The only escape is an independent-sample split (T1).

**The predictor lives on the confounded axis.** Batch = timepoint. The attenuation (an
interaction) is batch-clean; the WT background move is not. A background-antagonism model
therefore links a clean outcome to a confounded predictor, and no amount of analysis fixes that —
only a design with dissociation-batch metadata or a second cohort would.

---

## 6. What would discriminate

| test | what it is | what its outcome means |
|---|---|---|
| **MYC blot** (pending) | protein at 6W vs 12W in Myc+ | protein **down** at flat message: H3 wins outright and most of this note is moot. Protein **flat**: H3 is demoted to its cofactor form, the field stays open, T1 becomes the priority |
| **T1** de-confounded split | reuse script 40 PART C: split *both* WT groups so no mouse is shared between predictor and outcome, then ask whether the background's direction predicts attenuation | the **gate** for H2. If it fails, background antagonism has no transcriptomic support and T2/T3 are moot |
| **T2** TF-regulon excess scan | across the library TF lanes and the ChEA shortlist, whose regulon attenuates *more* than expression-matched background (script 34's null machinery) | a **named** repressor would make H2a publishable; silence bounds it |
| **T3** repressive-axis moderation | script 30's model with a repressive axis — the untested arm of Issue #5 | direct but underpowered at n=24, and endogenously biased |
| **T3b** MLX co-IP | MLX:MondoA against MLX:MXD1/MXD4/MNT at 6W vs 12W | the **only** test of the MLX-release route (§4.4). A bench assay, not an analysis |
| **T4** MYC CUT&RUN/ChIP | occupancy at the same promoters, 6W vs 12W | occupancy **unchanged** with output down: the block is at co-activator/elongation. Occupancy **down** at constant protein: competition, i.e. H2b |
| **T5** deconvolution / single cell | the H4 test deferred at Issue #6 | separates dilution from a per-cell change |

Priority order once the blot lands: blot, then T1 (cheap, decisive for H2), then T5 and T2.

---

## 7. For a wider audience

Between six and twelve weeks the Myc-driven programme in these mammary glands weakens by about
half, even though the oncogene itself is not switched down — if anything the difference between
Myc-carrying and control tissue grows, while the transcriptional response it produces shrinks. Two
readings of this are on the table, and they are not the same story. The first is that the
programme simply runs down: Myc keeps pushing, the tissue has finished responding, and the
weakening needs no external cause. The second, which the data cannot currently rule out, is that
the tissue itself is doing the blunting — the normal gland is maturing over exactly this window,
and the genes Myc most wants to drive are the ones it is turning down, so the oncogene arrives at
a substrate that is progressively harder to move. The two readings look identical in the RNA,
because both predict what we see: the same genes, in the same order, at half the amplitude. What
separates them is whether MYC is still doing the same thing to each cell. That is not a question
transcriptome data can answer, and it is a good example of where a measurement, not more analysis,
is the way forward: a MYC protein blot across the window, and if the protein is unchanged, a map
of where MYC actually sits on the genome at each age.

---

## 8. Appendix

### Scan rosters (2026-07-25)

**E-box binders and MYC network:** Usf1, Usf2, Bhlhe40, Bhlhe41, Tfe3, Tfeb, Mitf, Mnt, Mxd1,
Mxd3, Mxd4, Mxi1, Mga, Max, Mlx, Mlxip, Mlxipl, Srebf1, Srebf2, Arnt, Hif1a.

**HLH sequestrators and E-proteins:** Id1, Id2, Id3, Id4, Hes1, Hey1, Hey2, Heyl, Tcf3, Tcf4,
Tcf12.

**Repressors, NRs, corepressors, chromatin, DREAM (111 resolved of 117 queried):** the three
rosters above plus Arntl, Clock, Npas2, Cry1, Cry2, Per1-3, Nr1d1, Nr1d2, Rora/b/c, Ppara, Ppard,
Pparg, Ppargc1a, Ppargc1b, Esrra, Esrrb, Esrrg, Nr3c1, Pgr, Esr1, Ar, Vdr, Thra, Thrb, Rxra, Rxrb,
Nr2f1, Nr2f2, Nr2f6, Nr0b2, Nr1h2, Nr1h3, Nr4a1, Nr5a2, Ncor1, Ncor2, Sin3a, Sin3b, Rcor1, Ezh2,
Suz12, Eed, Rnf2, Bmi1, Cbx7, Hdac2, Hdac4, Hdac7, Sirt3, Sirt6, Zbtb17, Kdm5b, Setdb1, E2f4-E2f8,
Rb1, Rbl1, Rbl2, Lin9, Lin37, Lin52, Lin54, Tfdp1, Cdkn1a, Cdkn1b, Cdkn2a, Cdkn2b, Trp53, Elf5,
Gata3, Foxa1, Stat5a, Stat5b, Prlr, Sox9, Cebpb, Klf4, Klf9, Klf15, Foxo1, Foxo3, Notch1, Runx1,
Zeb1, Snai2, Twist1, Nrf1, Gabpa, Tfam. *Absent from the count matrix:* Arntl, Nr0b2, Lin54.

**Cell-identity markers:** Krt5, Krt14, Krt17, Trp63, Acta2, Myh11, Oxtr, Krt8, Krt18, Elf5, Prlr,
Csn2, Esr1, Areg, Epcam, Ptprc, Pecam1, Adipoq.

**MLX arm and targets:** Mlx, Mlxip, Mlxipl, Txnip, Arrdc4, Slc2a1, Ldha, Pklr, Fasn, Acaca, Scd1,
Elovl6, Acly, Gck, Khk, Aldob.

### Method and provenance

Selection rule: baseMean >= 20; raw (unshrunken) log2FC; genome-wide BH padj as reported by
DESeq2 for each of the five contrasts in `results/interaction_results.rds`; symbol-to-Ensembl via
`results/combined_df_annotated.rds`. Group means are median-of-ratios normalised counts from
`results/dds_int_run.rds`.

These scans were run **read-only in a scratchpad**, not as numbered pipeline scripts. They have no
provenance beyond this note and **must be re-derived inside a numbered script before any of these
numbers is cited** in the manuscript.

Other sources: `results/background_vs_myc.rds` (script 40 — facts 1, 4, the tier table, the
artifact ledger), `figures/figS8_myc_network_levels.R` (fact 2),
`outputs/myc_stability_panel/README.md` (script 41 — fact 5),
`scripts/34_death_priming_reassessment.R` (fact 3),
`scripts/30_attenuation_moderation.R` (section 4.5).
