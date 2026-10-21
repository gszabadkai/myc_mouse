---
name: mtpct-is-three-way-confounded
description: "mt% (mtDNA-encoded read fraction) is a three-way-confounded, death-AGNOSTIC axis (real content x proliferation denominator x prep-leak); the mitonuclear imbalance IS mt%, so imbalance/mtDNA-arm couplings to priming are the mt% confound, not biology"
metadata:
  node_type: memory
  type: project
  originSessionId: 71ffe8cc-b22d-412b-9e40-f22881e23306
  modified: 2026-07-20T09:38:18.192Z
---

**What mt% is, and why it invalidates the imbalance-as-death reading** (2026-07-20, read-only from
`results/mtdna_axis_and_coupling_null.rds` + `results/mitopps_scores.rds` + `results/dds_int_run.rds`,
Spearman n=24). mt% = the fraction of reads from the 13 mtDNA-ENCODED genes (mt-Nd1..mt-Cytb;
Ensembl ENSMUSG00000064341 etc., NOT symbols in the count matrix). Range **3.4-40%, median 10.5%**,
within-GROUP variance dominant (not a clean group effect), outlier MYCF53_6d ~40%.

**It is THREE things at once, unresolvable:**
- (a) real mtDNA content / copy number;
- (b) a proliferation DENOMINATOR effect -- fast-cycling cells make more nuclear mRNA, diluting the
  mt fraction: `cor(mt%, proliferation) = -0.61`;
- (c) prep-quality / mitochondrial leak from the warm enzymatic digest: `cor(mt%, IEG) = +0.74`.
No way to separate them: enzymatic dissociation, NO QC metadata (viability/RIN), batch=timepoint
(see [[batch-equals-timepoint]], [[timepoint-cohort-confound]]).

**The imbalance and the mtDNA-arm ARE mt%:** `cor(imbalance, mt%) = -0.92`; `cor(mtDNA-arm mitoPPS,
mt%) = +0.95`. So the couplings that looked strong are the mt% confound:
- imbalance<->priming **+0.70 -> +0.34** on mt%-adjust;
- mtDNA-arm<->priming **-0.71 -> -0.20** on mt%-adjust.
The ONE mt%-ROBUST lead is **NUCLEAR OXPHOS<->priming +0.46 -> +0.48** (`cor(nuclear OXPHOS, mt%)
= -0.71`) -> **BH3 profiling** (mid-pack, script 38; see [[death-priming-change-untested]]).

**The SIGN kills the prep-stress->death reading** (the author's point, 2026-07-19/20). If mt%
variation were driven by prep-stress->DEATH, mt% would correlate POSITIVELY with death priming. It
does not: mt%<->priming is NEGATIVE, and `priming<->IEG = -0.45`. Everything collapses to ONE
proliferation<->stress axis: **growth pole = low mt% / low IEG / HIGH priming; stress pole = high mt%
/ high IEG / LOW priming.** So **mt% is DEATH-AGNOSTIC** -- prep-stress raises mt% compositionally
but does NOT raise priming.

**How to apply:** any mtDNA/imbalance-based death or coupling claim is the mt% confound -- do not
cite it. mtDNA content questions -> **mtDNA qPCR** (mt-ND2 / nuclear Pecam1, the Lesner method;
[[lesner-mtdna-external-validation]]). The one transcriptomic death lead is nuclear-OXPHOS<->priming
-> BH3 profiling. Full reasoning: `docs/2026-07-19_oxphos_axis_biology_and_mtdna_priming.md`. See
[[mitonuclear-imbalance-untested]], [[sample-level-couplings-need-a-null]].
