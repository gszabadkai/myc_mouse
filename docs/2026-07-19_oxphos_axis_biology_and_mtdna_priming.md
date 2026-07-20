# What OXPHOS being the main axis means -- biology, loading, Myc, and the mtDNA/priming question

**Date:** 2026-07-19. **Purpose:** a biological reading of scripts 37-38 for the writeup. Covers (1)
what "OXPHOS" actually is in this dataset, (2) what it means that it is the dominant axis, and what
"loading" and the "Myc effect" mean, (3) how OXPHOS differs from mito-biogenesis in gene content and
why the two behave differently, and (4) a tested verdict on whether the **mtDNA <-> priming**
correlation is citable. Companion to `2026-07-18_what_survives_the_global_shift.md` and the five-
question synthesis. Numbers are from the committed scripts (37/38) plus a targeted mt% check run
2026-07-19.

---

## 1. What "OXPHOS" actually is in this dataset

"OXPHOS" is not one thing measured one way -- it appears in two spaces, and the paper should be
explicit about which:

- **As a gene programme:** the **nuclear-encoded** respiratory machinery -- the ~188 genes of the
  five electron-transport-chain complexes (`Ndufa/b*`, `Sdha-d`, `Uqcr*`, `Cox4-8`, `Atp5*`), their
  assembly factors, electron carriers, and cristae proteins. The **13 mtDNA-encoded** ETC subunits
  (`mt-Nd*`, `mt-Co*`, `mt-Atp*`, `mt-Cytb`) are deliberately held in a **separate** pathway
  ("mtDNA-encoded OXPHOS subunits") throughout -- this separation is load-bearing for section 4.
- **As a mitoPPS score (relative):** the **share of the mitochondrial transcriptional budget** going
  to respiration -- a within-compartment ratio, so it reads as *reallocation*, not absolute amount.
- **As a GSVA/z-score composite (the axis work, scripts 35-37):** a per-sample activity score whose
  covariation with other programmes is what "the dominant axis" is built from.

So when we say "OXPHOS is the main axis," we mean the **nuclear respiratory programme's per-sample
activity score is the single strongest-loading composite** in the 884-programme matrix.

---

## 2. OXPHOS is the main axis -- what "loading" and the "Myc effect" mean

**The finding (script 37).** Across the 24 mice the 884 pathway scores do not vary in 884
independent directions; they move along essentially **one dominant axis** (PC1 = 76% of variance,
effective dimensionality 1.7 of a possible 23). A pathway's **loading** is how tightly it moves with
that one axis. OXPHOS loads **0.968** -- it does not merely participate in the dominant axis, it
**essentially IS the axis**.

**What "loading" means, biologically.** The dominant axis is the coordinated state that rises and
falls as a whole with the things that vary between these samples -- Myc activity, proliferation,
epithelial purity, preparation stress. A high loading means OXPHOS moves in lockstep with that whole
coordinated state. It is a *description of what the axis is made of*, not a claim that OXPHOS drives
anything -- with n=24 and everything collinear, "drives" is not separable from "co-varies with."
Honest phrasing: **respiration is a defining component of the single broad transcriptional state that
distinguishes these samples.**

**Important honesty point (script 37).** OXPHOS is not *uniquely* central. It sits in a tightly
co-loading **Myc-metabolic-proliferative cluster** -- TCA (0.97), nucleotide (0.95), MYC targets
(0.93), biogenesis (0.91), proliferation (0.86) all load nearly as high. The dominant programme is
that **whole co-regulated core**, and OXPHOS is a leading member of it. This is exactly *why* the
older "OXPHOS is central to the phenotype" coupling claim is untestable: nothing in the cluster is
statistically separable from the rest at this n.

**Quantified as a ranking (script 37 PART C3 / figure `D_loading_enrichment_mito.pdf`, added
2026-07-20).** Ranking all 884 sets by loading: mito-defined sets are **44% of the library (385/884)
but 100% of the top 100** -- a real, striking enrichment (Wilcoxon |loading| mito vs rest underflows
double precision, p < 1e-300; median |loading| 0.97 mito vs 0.62 non-mito). But it is **mito-LED,
not mito-SPECIFIC**: non-mito growth programmes reach 0.95-0.98 (MYC-target/METABRIC breast-cancer
biclusters 0.98, pentose-phosphate 0.97, nucleotide 0.96) vs the redox control 0.26; and the very top
is **construction-inflated** -- 94 of the top 100 are the Gray `_MITO` / `_LE_MITO` TF lanes, which are
MitoCarta subsets *by build* (`TFT_MYC_..._LE_MITO` loads 0.99 at rank 6 while the same TF's non-mito
lane loads 0.01 at rank 881). So the ranking is **not** an argument for mitochondria-primary
regulation -- it describes a coordinated Myc anabolic-proliferative-mitochondrial growth state, and
loading is covariation, not primacy (which needs perturbation, not n=24 bulk).

**What the "Myc effect" means -- two distinct statements.**
- **Level (robust):** Myc raises OXPHOS (the content effect, +21-27% mitochondrial content; absolute
  nuclear OXPHOS up, d~1.2), and OXPHOS is a top member of the coordinated Myc-driven state. This is
  a genotype contrast and it is safe.
- **Change in loading (exploratory only):** within a genotype, does Myc make OXPHOS track the axis
  *more tightly*? The `S ~ G*genotype*timepoint` model hints yes (relative +71%), but the slopes are
  near-collinear at n=6/group and mito-biogenesis shows no such change -- so this is a **lead, not a
  result**. Do not headline it.

---

## 3. OXPHOS vs mito-biogenesis: different genes, different behaviour

**They are largely different gene programmes.** Union membership: **OXPHOS ~188 genes**, **biogenesis
~249 genes**, **overlap only 17**. They are not two labels on one list.

| | OXPHOS | mito-biogenesis |
|---|---|---|
| what the genes do | **run** respiration: ETC complex subunits, assembly factors, electron carriers, cristae | **build and maintain** the organelle: mtDNA replication/repair/nucleoid, mt-transcription, mt-translation, mitoribosome, protein import |
| examples | `Atp5a1/b/c1`, `Ndufa*`, `Cox4-8`, `Uqcr*`, `Sdha-d`, `Acad9` | `Atad3a`, `Alkbh1`, `Apex1`, `Aars2/Cars2` (mt-tRNA synthetases), `Aurkaip1`, `Coa3`, mitoribosomal proteins |
| reads as | the **functional end-product** (respiratory capacity) | the **supply/construction** line (organelle mass) |

**Why they behave differently in our data.** Both rise wholesale with Myc (the content effect) and
both load high on the dominant axis. The difference is in what survives once the shared axis / Myc
dose is removed:
- **OXPHOS retains co-variation with outcomes** -- its couplings sit *at* the correlation ceiling,
  and in the clean mitoPPS space it couples to death priming *above* the null (section 4). It behaves
  like a functional readout wired to the cell's energetic/redox/apoptotic state.
- **Biogenesis is a Myc-dose bystander** (scripts 35/36) -- once Myc/the global axis is removed its
  partial couplings fall *below* an arbitrary programme, and in mitoPPS space biogenesis->priming
  does **not** clear the null (+0.37 vs a 0.48 design-adjusted null). It behaves like a passenger that
  scales the organelle's **mass** without dictating what the organelle **does**.

**The one-line biology:** Myc expands the mitochondrial compartment (biogenesis builds mass, a
passenger), but it is **respiratory capacity (OXPHOS)** that co-varies with the proliferative and
apoptotic state. Mass and function are induced together but are not the same signal -- which is the
paper's "two mitochondria" distinction (biosynthetic organelle vs functional/apoptotic gateway).

---

## 4. Is the mtDNA <-> priming correlation citable? **No -- it is the mt% confound.**

This was tested directly (2026-07-19), because it is exactly the case the mtDNA-arm border warned
about (`2026-07-18_what_survives_the_global_shift.md`, point 3: mitoPPS cancels a global shift but
**not** variation in the mitochondrial fraction itself).

**The raw coupling looks strong and clean:**

| coupling (mitoPPS space) | raw Spearman | design-adj | percentile | vs null 0.42/0.48 |
|---|---|---|---|---|
| **mtDNA-encoded OXPHOS <-> priming** | **-0.71** | **-0.74** | 95th | well above |

At face value this is far stronger than the nuclear OXPHOS<->priming lead (+0.45) from script 38 and
clears every null. **But it is an artifact of the mitochondrial fraction:**

- `cor(mtDNA-arm, mt%) = 0.952` -- **the mtDNA-encoded mitoPPS score essentially IS mt%.** (This is
  expected: mt% is the read-share of the 13 mtDNA genes, and the mtDNA-arm pathway is those genes.)
- `cor(priming, mt%) = -0.60` -- priming also tracks mt%.
- **Adjust for mt% and the coupling collapses: -0.71 -> -0.20** (partial), below the null.
- The **imbalance** (nuclear OXPHOS - mtDNA) <-> priming behaves identically (raw -0.72; `cor(imbalance,
  mt%) = -0.94`) -- same artifact, since the imbalance is also ~mt%.

**mt% is the one axis we cannot resolve** as biological vs technical: MECs are enzymatically
dissociated (not FACS-sorted), there is no dissociation-batch / viability / RIN metadata, and batch =
timepoint. So a correlation that lives almost entirely on the mt% axis is **not citable as a
biological coupling** -- it is confounded with an unresolved technical variable.

**The contrast that makes this clean -- the nuclear coupling IS mt%-robust:**

| | cor with mt% | <-> priming raw | <-> priming adj. for mt% |
|---|---|---|---|
| **nuclear OXPHOS** | -0.71 | +0.45 | **+0.49 (survives)** |
| **mtDNA-encoded OXPHOS** | **+0.95** | -0.71 | **-0.20 (collapses)** |

The nuclear OXPHOS<->priming lead from script 38 does **not** depend on mt% (it is essentially
unchanged, +0.45 -> +0.49, when mt% is removed). The mtDNA coupling is ~72% mt%. So the two OXPHOS
arms give **opposite-signed** couplings to priming for opposite reasons -- and only the nuclear one is
a real (if weak, mid-pack) lead.

**Verdict.** **Do not cite an mtDNA (or mitonuclear-imbalance) correlation with priming.** It is the
mtDNA-fraction confound in a new guise, exactly as the border predicted; it is the reason the paper
should test the mitonuclear/mtDNA questions by **mtDNA qPCR**, which measures the physical ratio
directly and is immune to this. The citable-but-weak transcriptomic statement is the **nuclear**
OXPHOS<->priming coupling (mt%-robust, above the mitoPPS null, but only mid-pack at n=24) -> **BH3
profiling** is its bench test.

---

## 5. Bottom line for the narrative

- **Cite:** Myc raises mitochondrial content and nuclear OXPHOS (contrasts); OXPHOS (respiration) is a
  defining member of the single dominant transcriptional state, and Myc raises it (loading level);
  OXPHOS and biogenesis are distinct gene programmes (188 vs 249 genes, 17 shared) that are induced
  together but behave differently -- function co-varies with phenotype, mass does not.
- **Frame as a weak lead (not a result):** nuclear OXPHOS <-> death-priming in mitoPPS space
  (mt%-robust, above null, mid-pack) -> BH3 profiling.
- **Do not cite:** the mtDNA <-> priming / mitonuclear-imbalance <-> priming correlation (it is the
  mt% confound; -0.71 -> -0.20 on adjustment) -> mtDNA qPCR.
- **Do not cite:** any within-group loading-*change* by Myc as a result (exploratory, near-collinear).
- **Untestable here (a design answer):** MYC- vs PGC1a/ESRRA-biogenesis toward death -> ESRRA/PGC1a
  perturbation.

## 6. External validation (added 2026-07-20): Lesner et al.

The MYC nuclear-up / mtDNA-down **discordance** described above is independently established in
**TCGA (flash-frozen, no dissociation) + mouse HCC + proteomics + mtDNA qPCR** by Lesner et al.
(bioRxiv 2026.07.13.738248) -- see `docs/2026-07-20_lesner_mtDNA_external_validation.md`. This (a)
**validates the phenomenon** and kills the "our enzymatic-dissociation artifact" worry for the
imbalance; (b) **reframes the death arm** as protective **mitophagy / turnover** (MYC->NRF1->DRP1/
FUNDC1), apoptosis engaged only by forced fusion + BH3 mimetics -- consistent with our finding that
the mtDNA/imbalance<->priming coupling is the mt% confound and only **nuclear** OXPHOS<->priming is a
(weak) lead; and (c) **endorses our experimental plan** -- mtDNA qPCR (mt-ND2 / nuclear), a non-mito
MYC control, BH3 profiling, and NRF1/DRP1/FUNDC1 perturbation (the biogenesis-axis perturbation script
38 said the MYC-vs-PGC1a death question requires). Honest nuance: Lesner's is a MYC *contrast*; our
nuclear-up half reproduces as a clean genotype effect but our mtDNA-down half is genotype-independent
(p~0.49), so Lesner validates the **question** and pinpoints our underpowered arm, it is not a claim
we demonstrated the contrast on the mtDNA side.
