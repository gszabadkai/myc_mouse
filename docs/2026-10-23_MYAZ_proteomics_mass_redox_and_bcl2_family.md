---
date: 2026-10-23
status: Analysis of an existing dataset. No new experiment. Claims marked as
  licensed or not licensed at the end of each section.
subject: MYAZ quantitative proteomics, PGC-1alpha / NRF1 / Ndi1 against Control
source-file: 250131_MYAZ_mouse_and_inserts_mods_candidates_clerared_filters.tsv
purpose: >
  Three questions the manuscript currently answers by inference rather than by
  measurement: (i) whether Ndi1 raises mitochondrial mass, (ii) whether the
  constructs move matrix antioxidant capacity in the same direction, and
  (iii) what the BCL-2 family does at protein level in a non-lethal condition.
relates-to:
  - Fig. 2 main text, the sentence "Ndi1 and LbNOX are single enzymes that raise
    NADH oxidation without mitochondrial biogenesis"
  - Fig. 3f,g (p19ARF and MYC under PGC-1alpha, NRF1, Ndi1) - flagged, see 5
  - Abstract clause "independently of ATP production, mitochondrial mass, and
    membrane potential"
  - 2026-10-22_six_gene_composite_provenance.md (configuration framing)
---

# MYAZ proteomics: mass, matrix redox capacity, and the BCL-2 family

## 1. The dataset

Spectronaut-style differential abundance table, protein-group level, candidates
cleared filters. Three comparisons against Control in MYAZ cells:

| comparison | proteins |
|---|---|
| Pgc1a / Control | 6,991 |
| Nrf1 / Control | 6,991 |
| Ndi1 / Control | 6,983 |

16 ratios per protein in every arm. No LbNOX arm. No dox arms, so these are the
constructs on the parental MYAZ background.

**Construct validation.** The yeast NDI1 protein (UniProt P32340, "Rotenone-
insensitive NADH-ubiquinone oxidoreductase") is **+3.89 log2** in the Ndi1 arm
and **-0.007** in the Pgc1a arm. NRF1 protein is **+1.86** in its own arm.
PPARGC1A itself is not quantified above noise in any arm (-0.47, Q 0.30 in the
Pgc1a arm), but its downstream signature is unambiguous: ESRRA **+1.21** and
TFAM **+0.33** in the Pgc1a arm only.

## 2. Method note: the drift control

The Ndi1 arm carries a small global downward shift in mitochondrial proteins.
Raw log2 ratios are therefore not comparable across arms without a reference
distribution. Every claim below is stated as a **percentile within that arm's
own mitochondrion-annotated distribution** (GO cellular component, ~1,265
proteins per arm), not against zero.

| arm | mito median | non-mito median |
|---|---|---|
| Pgc1a | **+0.238** | -0.033 |
| Nrf1 | +0.026 | +0.005 |
| Ndi1 | **-0.066** | -0.009 |

## 3. Mitochondrial mass: measured, not inferred

Median log2 ratio by complex:

| | CI (42) | CII (4) | CIII (8) | CIV (17) | CV (16) | mitoribosome (77) |
|---|---|---|---|---|---|---|
| Pgc1a | +0.68 | +0.93 | +0.60 | +0.57 | +0.73 | +0.43 |
| Nrf1 | +0.12 | +0.10 | +0.02 | -0.06 | -0.02 | +0.30 |
| Ndi1 | -0.27 | -0.07 | -0.14 | -0.26 | -0.14 | -0.06 |

Supporting markers in the Pgc1a arm: CS +0.84, VDAC1 +0.42, TOMM20 -0.04.
In the Ndi1 arm: CS -0.13, TOMM20 -0.15, VDAC1 +0.01.

**Licensed.** PGC-1alpha produces coordinated mitochondrial biogenesis in MYAZ
and Ndi1 does not. The complex I fall in the Ndi1 arm is what a complex I
bypass would predict.

**Not licensed.** Nothing about the LbNOX lines, which are absent here. Nothing
about mitochondrial mass in iMMECs.

**Flag.** NRF1 does not produce global biogenesis in this dataset. It raises
complex I subunits (NDUFB8 +0.47, NDUFA9 +0.34) and the mitoribosome, and
leaves complexes III to V flat or slightly lower. This differs from the
western-blot record, which has NRF1 raising ETC subunits generally. Either the
timepoints differ or the blot panel was complex I weighted. Worth resolving,
because a complex I selective NRF1 would be a **better** arm for the NADH
oxidation argument, not a worse one.

## 4. Matrix redox capacity: the two arms move in opposite directions

Percentile within each arm's mitochondrial distribution.

| protein | Pgc1a | Nrf1 | Ndi1 |
|---|---|---|---|
| SOD2 | **98th** (+1.40) | 10th (-0.26) | **11th** (-0.30) |
| PRDX3 | 90th (+0.82) | 61st | 39th (-0.11) |
| PRDX5 | high (+0.83) | 20th (-0.18) | 36th |
| TXN2 | 82nd (+0.61) | 84th (+0.26) | **3rd** (-0.45) |
| IDH2 | 69th (+0.43) | 19th (-0.17) | 28th (-0.16) |
| AIFM1 | 88th (+0.73) | 26th (-0.12) | 74th (+0.03) |
| NNT | 48th (-0.10) | 48th (+0.15) | 20th (-0.00) |

**This is the most useful result in the dataset for the manuscript.** PGC-1alpha
and Ndi1 both kill on a MYC background whilst moving matrix antioxidant
capacity in opposite directions, by the same logic that membrane potential is
already excluded in Fig. 2. NRF1, which also kills, sits with Ndi1 on SOD2 and
IDH2.

**Licensed.** Antioxidant capacity is not the shared variable. A reviewer
raising ROS as the alternative to NADH oxidation has to account for three
constructs that kill whilst the matrix antioxidant system goes up steeply in
one and down in the other two.

**Not licensed.** That ROS is unchanged, or that ROS plays no part. Protein
abundance is capacity, not flux. No ROS measurement is in this dataset.

**Also licensed, separately.** NNT is flat in every arm, so the route from
matrix NADH to NADPH through NNT is not transcriptionally or
translationally adjusted here. IDH2 falling in the Ndi1 and NRF1 arms is a
separate observation and is not large.

**AIFM1.** Proposed earlier as a candidate NADH sensor. Its abundance tracks
biogenesis (88th percentile under PGC-1alpha) and is unmoved by Ndi1. That does
not exclude AIFM1, whose NADH-dependent dimerisation is not an abundance
property, but abundance cannot be used as support.

## 5. BCL-2 family at protein level

Not detected at all: BBC3, BIK, BMF, PMAIP1, HRK, BCL2A1, BCL2L2. This is
expected for low-abundance BH3-only proteins in DIA, and it means **the
proteomics cannot speak to PUMA.**

Detected, as percentile within each arm's mitochondrial distribution:

| protein | peptides | Pgc1a | Nrf1 | Ndi1 |
|---|---|---|---|---|
| BCL2L1 (BCL-xL) | 7 | 15th | 93rd | **9th** (-0.32, Q 0.017) |
| BCL2 | **2** | **1st** (-0.45) | **3rd** (-0.44) | **2nd** (-0.54) |
| MCL1 | 12 | 51st | 66th | 80th |
| BAK1 | 9 | 56th (+0.29, Q 0.001) | 46th | 77th |
| BAD | 4 | 42nd | 95th (+0.45, Q 0.040) | 11th |
| BAX | 7 | 4th | 80th | 8th |
| BID | 6 | 1st | 91st | 9th |
| BCL2L11 (BIM) | 3 | 15th | 47th | 5th |
| BOK | 5 | 13th | 24th | 14th |

Three readings, in descending order of how well they are supported.

**(a) BCL2 protein falls in all three arms**, to the 1st to 3rd percentile of
the mitochondrial distribution each time, and it is the only change shared by
all three. **Caveat that must travel with it: two peptides.** This is the
weakest identification in the table and needs a blot before it is used.

**(b) Ndi1 lowers BCL-xL.** 9th percentile, Q 0.017, 7 peptides, well clear of
the arm's -0.07 drift. Attractive because it would let Ndi1 and the mouse gland
move the same ratio from opposite ends, the gland through the PUMA numerator
and Ndi1 through the BCL-xL denominator. **But see (c) before using it.**

**(c) The honest reading: Ndi1 lowers most of the detected family.** BCL-xL
9th, BAX 8th, BID 9th, BIM 5th, BAD 11th, BOK 14th, with MCL1 (80th) and BAK1
(77th) the only two holding position. Activators fall alongside guardians, so
this is a general downshift rather than a clean activator-to-guardian ratio
change. Computing every activator:guardian protein ratio confirms it: only
BAK1:BCL2L1 moves positive in two arms (+0.36 Ndi1, +0.37 Pgc1a) and it runs
the other way in NRF1 (-0.40). **No ratio moves coherently across the three
arms.**

**Licensed.** That the apoptotic protein configuration is altered by all three
constructs, and that BCL2 falls in each. That Ndi1 does not act by raising the
BCL-2 family.

**Not licensed.** A statement that the constructs shift the activator-to-
guardian balance in a common direction at protein level. They do not, in this
dataset.

**Note for Fig. 2c.** The text reads "PGC-1alpha induced PUMA without inducing
BIM or BAX". That is consistent here (BIM 15th, BAX 4th). The proteomics adds
BAK1 rising (Q 0.001) and BCL2 falling, which the blot panel did not carry.
PUMA is therefore not the only member PGC-1alpha moves.

## 6. A flag against Fig. 3f,g

Fig. 3g is currently a "get panel" item, from MYAZ qPCR, stating that MYC and
p19ARF fall under PGC-1alpha and NRF1 and are unchanged under Ndi1.

**MYC protein is not quantified in this dataset at all.** CDKN2A is, on 6
peptides, and is **unchanged in every arm**: Ndi1 +0.017, NRF1 +0.175,
Pgc1a +0.109, none significant.

Two reasons this is a flag rather than a contradiction. CDKN2A encodes p16INK4a
and p19ARF from alternative reading frames sharing exons 2 and 3, so a protein
group assembled from shared peptides does not report p19ARF specifically. And
the Fig. 3f,g observation is in the **surviving fraction** after the apoptotic
wave, which this dataset may precede.

**Action.** Before Fig. 3g is built, confirm which timepoint and which fraction
the qPCR came from, and whether the p19ARF antibody or primer set distinguishes
the two frames.

## 7. Open items

1. **BCL-xL by western in the Ndi1 and LbNOX lines.** Decides whether (b) above
   survives (c). Lysates are likely banked.
2. **BCL2 by western**, given the two-peptide identification.
3. **Resolve the NRF1 biogenesis discrepancy** against the blots.
4. **Confirm the Fig. 3f,g timepoint** before the panel is built.
5. The percentile tables above were computed from the supplied TSV only. If the
   raw Spectronaut report is available, the same analysis should be rerun with
   the per-replicate values so the 16 ratios can be inspected directly.
