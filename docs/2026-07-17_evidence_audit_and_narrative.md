# Evidence audit, 2026-07-16/17 — what survives, what doesn't, and what the paper can say

Written 2026-07-17 on `paper-figures`. **Read this first, then the walkthrough**
(`docs/2026-07-13_BlockA_revision_walkthrough_and_intro_alignment.md`), which now carries a status
banner pointing back here.

**How this started and where it went.** The question was narrow: *"the paper says Myc drives a robust
increase in mitochondrial biogenesis — is that mitochondrial content per cell? Is there any
transcriptomic parameter that could assess it?"* Answering it produced one solid new result and, by
following the author's successive challenges, dismantled the bridge that held the manuscript's arc
together. Every number below is reproducible from a committed script (32/33, commits
`50ad21c`..`148dbe0`); nothing here is quoted from conversation only.

---

## 1. The ledger

### STANDS — genotype axis: powered, litter-controlled, null-independent, adjustment-robust

| claim | evidence | where |
|---|---|---|
| **Myc raises mitochondrial content +21–27%** | mass markers +26.5% → **+21.1% (p=0.0006)** after adjusting for prep-stress + contamination; nuclear MitoCarta +18.9% → +14.9%; mitoribosome +35.3% → +28.4%; import +33%. Mass panel **coherent: 13/14 genes up, 11 at p<0.05**. All 17 panels leave-one-out stable. | script 32 |
| Myc raises **absolute nuclear OXPHOS** | d=1.27 (script 29 Part B, VST z-composite) — now corroborated by a **different method**, d=1.15 (script 32, raw-count share) | 29, 32 |
| **Issue #4's attenuation is REAL — and now protected** | 12W is **not** noisier (median within-group SD 0.2935 vs 0.3515, ratio 0.83); Myc-contrast SEs ~identical (0.251 vs 0.259); so the DE collapse **1967 → 135** tracks median \|LFC\| falling **0.268 → 0.204** at constant SE = genuine magnitude attenuation, not power. The IEG/batch axis is genotype-independent, so it cancels inside each timepoint's contrast. | script 33 PART E2 |
| Issues #1, #2, #3, #5, #6 | all genotype/interaction contrasts — untouched by everything below | 26–31 |
| The **death phenotype** (Myc kills at 6W, not 12W) | external IHC / inducible model — **never depended on this RNA-seq** | — |
| mt% is **not contamination** | mixing would need a contaminant mt-share of **126–737%** (impossible at 3–12% of cells); and nuclear MEC genes track contamination *negatively* while mt% tracks it *positively* — **mixing cannot give opposite signs** | script 33 PART A |
| IEGs are **not** the TEB/proliferative programme | IEG ~ proliferation **−0.76** within 6W; ~TEB sets −0.57; and IEGs are **lower** at 6W (0.72/0.60) than 12W (1.03/1.02) — the TEB story predicts the opposite on direction **and** timing | script 33 D3 |

### UNSUPPORTED — not refuted; untested, or fails when tested (n=6–12)

| claim | what happened |
|---|---|
| **The mitonuclear imbalance** | **Has no supported contrast.** `imbalance ~ timepoint*genotype`: time **p=0.313**, genotype **p=0.523**, interaction **p=0.916**, Myc-at-6W **p=0.249**. The pattern the corpus quotes (0.130 / 0.499 / −0.456 / −0.173, *"peaks at 6W_pos, inverts by 12W"*) is **four group means reported descriptively and never subjected to a contrast**. Adjusting for the IEG axis takes the time beta −0.586 → −0.011. |
| **The death couplings** (script 25 Part B) | **Fail an empirical null.** Against all 884 library sets at 6W: imbalance +0.67 (null median +0.62, perm **p=0.49**); **bio_comp +0.83 (null median +0.79, p=0.50)**; MASC_comp +0.78 (p=0.10, marginal). The only axis that beats the null is the **prep-stress signature itself** (p=0.034). |
| **⇒ the mito→death bridge** | **Not established by these data.** Both halves fail *independently*: the imbalance has no statistics of its own, and its coupling to death priming is what any arbitrary gene set gives. |

**Read the `bio_comp` row twice.** r=0.83 looks compelling until you see the null median is **0.79**.
At n=12 with one dominant axis, 73–81% of *all* gene sets couple at |rho|>0.5.

**The structural cause.** Within 6W, **PC1 = 43%** of the variance and *every* measure loads on it:
IEG −0.71, proliferation +0.62, contamination −0.62, MYC activity +0.59, mt% −0.51. These are not
independent signals. That is why the null is flat — it is not a quirk of the death sets.

**The ceiling, stated honestly:** at n=6/group, *unsupported* is not *refuted*. The imbalance may be
real and underpowered. But it cannot be **reported as a finding** on this evidence, and "peaks at
6W_pos, inverts by 12W" reads as a demonstrated effect when it is four averages with p=0.31 between them.

### OPEN — not adjudicable with what is on disk

- **The IEG axis**: genotype-**independent** (p=0.49) but time-**associated** (p=0.0003). MECs are
  purified by **enzymatic dissociation only (no FACS)** — warm digest is exactly what induces the
  canonical dissociation signature (van den Brink 2017) — but with **no batch/viability/RIN metadata**
  the IEG panel is an *inferred* covariate. **Technical vs biological: unresolved.** It doesn't matter
  for the null (specificity fails either way) and it can't touch genotype claims.
- **Don't conflate two PC1s.** The *global* QC PCA (script 02, all 24, 35%) is **not** timepoint
  (p=0.81), **not** genotype (p=0.96), **not** IEG (rho −0.18). The 43% IEG-loaded PC1 is *within 6W*.

### RETRACTED this week — five, all caught by the author, none by my own validation

1. **"Commissioned but unbuilt"** — set-composition artifact. `MITOCARTA_TRANSCRIPTION` +43% is carried
   by **Mrpl12 +116%** (a *mitoribosomal* protein); `NUCLEOID` d=2.32 by **Atad3a +113%**. The actual
   core barely moves (Tfam +11% ns, Polg +3% ns, Twnk +1% ns).
2. **"Depth-confounded / TIME uninterpretable"** — a share is a **proportion**; depth cancels by
   construction, exactly as DESeq2's median-of-ratios does. Within cohort, depth predicts nothing.
3. **"Myc doesn't scale mtDNA = the imbalance in absolute share"** — compared a **significant p to a
   non-significant p** and called it a difference. The difference is p=0.57 raw / p=0.82 adjusted.
4. **The "direction test"** — **circular**. Arguing "nuclear mito genes should rise with mt% if it were
   mitochondrial; they fall; therefore it isn't" assumes the imbalance's *absence* to prove the
   measurement broken. The imbalance **predicts** that anticorrelation.
5. *(near-miss)* accepting r=0.75, p=0.005 as a mechanism **without a null**.

**Common root:** in each case a test whose premise was never checked. Recorded as a standing check in
memory (`dont-assert-a-mechanism-without-checking-it-can-work`).

---

## 2. Why this was invisible for so long

- **The imbalance and its death coupling were never tested.** Script 24 reports `mtnuc_index` as four
  group means; script 25 Part B reports couplings as raw correlations. Both were then quoted
  downstream as established. Neither had a contrast or a null.
- **The project already owned the right tool and didn't apply it.** `scripts/21_ap6_permutation_null.R`
  builds exactly this null for fGSEA — *"mito genes are highly expressed, so enrichment could be an
  artifact; test against a matched null"*. Nobody applied that logic to the **per-sample** couplings.
- **A load-bearing project fact was undocumented.** CLAUDE.md described "bulk RNA-seq of MMTV-Myc mice"
  and never mentioned MEC purification, so a stromal signal was read as tissue composition for several
  turns. Now recorded.

---

## 3. Narrative consequences — the section that matters

**The paper keeps both halves and loses the bridge.**

- **KEEP — the mitochondrial half, now stronger than before.** Myc builds mitochondria: content
  **+21–27%**, coherent across structural/import/ribosomal/OXPHOS arms, robust to composition
  adjustment. This is *new*, and it converts the vulnerable "robust increase in mitochondrial
  biogenesis" from an **enrichment-score** claim into a **quantified content** claim — with the rider
  that it is a share of transcriptome and therefore a **lower bound** (Myc amplifies total RNA/cell).
- **KEEP — the attenuation.** Real, magnitude-based, and now demonstrably not a power or batch artifact.
- **KEEP — the death phenotype.** External data. Untouched.
- **LOSE — the bridge between them.** *"At 6W the mitonuclear-imbalanced substrate is coupled to
  pro-apoptotic priming, so Myc's induction lands as killing; by 12W the imbalance resolves, the
  coupling collapses, killing fades"* → **not established**. And with it: *"the OXPHOS reprioritisation
  and the apoptosis desensitisation are two faces of one developmental mitochondrial maturation — the
  causal spine of the Introduction's closing sentence."*

**The honest arc:**

> A constitutively expressed Myc amplifies the normal pubertal proliferative/biosynthetic programme
> rather than creating a new one (#1–2), with a broad footprint whose phenotype-coupled core is
> respiratory/biosynthetic rather than biogenesis-generic (#3). It **increases mitochondrial content**
> (+21–27%) and raises absolute nuclear-OXPHOS transcript while **reallocating** the compartment (#4,
> 32). That effect **attenuates** 6W→12W as a loss of magnitude — ~2/3 oncogenic fade, ~1/3 wild-type
> convergence, sparing the buffered MYC core (#4–6). Independently, the gland becomes **resistant to
> Myc-induced apoptosis** by 12W (external phenotype). **Whether the mitochondrial state is *why* is an
> open question these bulk data cannot answer.**

**Sentences that must change** (see the walkthrough's §1B, §2, §3/4, all now marked):
- §1B "the spine, in one line" — asserts the bridge verbatim.
- §2 through-line — its final third *is* the bridge.
- §3/4 — the heading *"The apoptosis / permissive-window clause — SUPPORTED"* must become **supported as
  a phenotype, unsupported as a mechanism**. The intro clause itself survives; what dies is the claim
  that this RNA-seq explains it.

---

## 4. What would settle each open item — one experiment each

| open item | the experiment | also settles |
|---|---|---|
| Is there a mitonuclear imbalance at all? | **mtDNA qPCR** (copy number — the thing the transcriptome cannot see) | the "commissioned but unbuilt" question; whether mt% means anything here |
| Is mitochondrial content really up? | **TOMM20 / VDAC / CS blot** (or EM) | **tension A** (protein-level OXPHOS), and the lower-bound rider |
| Is the 6W substrate death-permissive? | **BH3 profiling** or caspase-by-state | the death spine's mechanism |
| Is the IEG axis technical? | a recorded **dissociation-batch / viability** variable | whether the time axis is usable for mt claims |

**Two bench experiments — mtDNA qPCR and the blot — close four open items between them**, including
tension A, which the walkthrough already routed to exactly that data.

---

## 5. Decisions for you

1. **The mitonuclear imbalance: re-test properly, or drop from the story?**
   *Recommend:* drop it as a **claim**, keep it as an observation with the ceiling stated, and let mtDNA
   qPCR decide. It cannot carry a figure at p=0.31.
2. **The death spine's Part B couplings: do they appear in a figure?**
   *Recommend:* **no**, unless they beat their null. The phenotype and the substrate description can be
   shown without the couplings.
3. **The intro's closing sentence.**
   *Recommend:* keep the apoptosis clause (external phenotype), delete the mechanistic bridge, and state
   the mito→death link as the open question that motivates the next experiment. This is a *better* paper
   position than an unsupported mechanism a reviewer would find.
4. **Does script 32's content result get a main figure?**
   *Recommend:* **yes.** It is the strongest genuinely new thing here, it answers the reviewer question
   pre-emptively, and it survives every adjustment.
5. **Whether to re-examine scripts 24/25's other claims** against the same null.
   *Recommend:* yes, but after the narrative is fixed — the null is cheap to apply (`null_for_axis()` in
   script 33 PART C) and any per-sample coupling in the corpus is exposed to it.

---

*Numbers: `results/mito_content_proxies.rds` (32), `results/mtdna_axis_and_coupling_null.rds` (33).
Verdicts print on sourcing. Sandboxes in both scripts walk every claim. Scripts 00–31 unmodified; no
re-fits; figure scripts remain 34+ and gated on the narrative.*
