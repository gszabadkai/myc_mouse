---
title: myc_mouse finalisation plan
last_updated: 2026-07-02
tags: [project/myc_mouse, planning, finalisation, hypotheses, action-points, timeline]
status: the plan. Review, then hand to a Claude Code Plan-Mode session to produce the build spec.
purpose: >
  The single working plan for finalising the myc_mouse arm of the manuscript:
  reframed hypotheses, principles, action points, build structure, figure
  architecture, and a 5-day draft timeline. Consolidates and replaces the
  2026-06-25 and 2026-07-02 (v2) drafts, which can be deleted.
relates_to:
  - "docs/branch_manifest.md (authoritative two-pipeline / script inventory)"
  - "mammary_geneset_library: 2026-06-19_library_to_myc_mouse_handoff.md"
  - "DT-pipeline note chain: 00. final (developmental trajectory) + 00.1/00.3/00.4/00.5/00.6"
  - "Menegollo, Bentham et al., Cancer Res 2024 (CAN-23-3172) - the analytical companion paper"
note: >
  Days are Day 1-5 (relative); fill calendar dates when you start. Budget ~4.5
  productive hours + ~0.5 h buffer for the Positron debug/re-run loop. Days 1-5 are
  the Block A build plus the Figure draft; a clean Block B rewrite for submission is
  a later ~30 h pass.
---

# myc_mouse finalisation plan

## 0. Scope and the one big realisation

Goal for the 5 days: a **draft** of this paper section with **draft Figures 1 and 2
plus ~3 supplementaries**, not a submission-ready version. Submission polish is ~30 h
of work; the draft is ~20-25 h, which fits, but only if the gates run first and the
schedule is protected with hard cut-lines.

**The gene set library does NOT need rebuilding.** Everything decided in the library
replan - method tags, Category 7 algebra, TF lanes, the MB bicluster fork-sets from
the 2024 paper - is already in the built v1.0 (all 9 category GMTs exported mouse +
human; `metabric_sets.rds` present; `gray_chea_mito_tf_shortlist.csv` generated). The
only thing the new analysis needs that might not be in the library is the AP6
comparator panel, and that is already loaded inside myc_mouse (`gene_sets_list.rds` =
89 sets incl. 50 MSigDB Hallmark). So: library = done; this week is myc_mouse only.
The library is consumed as a v1.0 snapshot (section 5), not rebuilt.

**Branch reconciliation is not a git merge, and it now has three parts.**
`new-analysis` is already the consolidated trunk (scripts 00-12); main's scripts are
archived byte-identical in `scripts/archive_main_pipeline/`. Reconciliation means:
(a) the cell-death question - resolved: run **both** branches, neither canonical, on
raw LFCs (AP-CD); (b) whether cell death stays a parallel branch or integrates with
the mitoPPS axis - it stays parallel and reports convergence; (c) the structural
layer - the Block A / Block B split, a git strategy, and the Step 1 / Step 2 pipeline
separation (section 5).

---

## 1. The reframed scientific question

Drop the **"hijack"** framing - it smuggles in an unproven mechanism (active
co-option). The clean backbone is **two trajectories**: a Myc- and a Myc+ mammary
trajectory, the latter terminating in tumour. Descriptive first: how do they differ,
and which difference is causal for transformation.

**Myc is on from embryonic development** (MMTV), so the 6W-12W window samples a
*slice* of an already-running process - an important slice (early tumourigenesis) but
far from the whole picture. Permission and selection are not mutually exclusive: a
developing tissue plausibly runs both. Development generates the *variation* (the menu
of states present); Myc acts on it through two coupled forces - instructively
expanding competent states and selectively culling vulnerable ones. Survivors are
both developmentally competent *and* stress-resistant. The open question is the
**relative weight** of expansion vs culling, and whether either runs through the
mitochondrion.

**Myc dose is stable across 6W-12W.** Our own data (Myc mRNA and blot, 6W vs 12W) and
the literature agree: the MMTV-Myc transgene is not induced at puberty and does not
decline in the adult. It is stably expressed at both stages. This makes the
decomposition clean - across the window one input is held fixed (the oncogenic
driver) and one input moves (the developmental/cellular context). So the interaction
term (the part of the Myc+ timecourse not explained by normal development, and the
"fading selectivity" seen in mitoPPS) cannot be attributed to a changing Myc dose;
whatever changes in what Myc *does* over the window must reflect the changing
substrate it acts on. That is exactly the "developmental state shapes the oncogenic
output" claim, on firmer ground. It reinforces principle 3.

---

## 2. Principles (carry into the analysis and the prose)

1. **Drop "hijack."** Two-trajectory, descriptive-first framing instead.
2. **Permission and selection are coupled, not exclusive.** Report relative weight,
   not a winner.
3. **The window is a slice.** Myc is on from embryogenesis; 6W-12W cannot witness the
   divergence event, only an already-diverged state. Reinforced by the stable Myc
   dose: constant driver, moving substrate.
4. **Decisive vs soft claims about "what precedes the window."** Decisive ("the
   competent state was set up before the window") is NOT answerable - no pre-pubertal
   samples. Soft ("already diverged at 6W"; "the 6W gland sits at developmental
   position X on the atlas") IS answerable. Single-cell signatures give trajectory
   *context*, not trajectory *history*.
5. **Readout vs decision.** Mitochondria are both sensor and effector of state. A
   clean mito signature can be a faithful *readout* of a decision made elsewhere.
   "Decided at the mitochondrion" is defensible only if the dominant force is
   selective culling running through mitochondrial apoptotic priming (BCL2-family
   setpoint, MOMP). Otherwise the honest claim is "mitochondria integrate and report
   the developmental-oncogenic state" - a readout claim, still real. Do not assert
   decisional/causal primacy on the strength of a metabolic readout.
6. **Justify the mito focus empirically, by decomposition** - not one whole-MitoCarta
   NES (see AP6).
7. **Coordination, not magnitude, is likely the true claim.** The signal is
   coordinated-but-weak-per-gene (gene-level interaction barely passes; fGSEA rescues
   a pathway signal). Lead with the coordination test; run magnitude too. The related
   abundance-vs-reprioritisation question (AP-abund) is time-dependent - do not
   pre-commit the prose to either word.
8. **Entropy is rejected here.** Shannon entropy conflates magnitude/direction/
   composition, is direction-blind, and does not measure coordination. The
   interpretable version of any heterogeneity claim is cross-sample variance/CV of a
   signature score. Reach for a diversity measure only if a single-cell/deconvolution
   arm ever materialises.
9. **Bulk cannot separate composition from per-cell regulation.** A signature score
   conflates "more competent cells" with "same cells upregulating the program."
   Resolving needs deconvolution (CIBERSORTx/MuSiC/Bisque) or orthogonal data
   (flow/IHC) - flagged, out of scope.
10. **mitoPPS is pairwise-ratio-based and needs linear-scale normalised counts.**
    Verified against the Monzel source: compute all pairwise pathway ratios per
    sample, normalise each by its global average, aggregate per pathway. Cancels total
    mito content and intrinsic pathway-scale differences. mtDNA-encoded subunits are
    isolated in a separate synthetic pathway. GSVA, by contrast, takes log-scale
    input - do not mix the two.
11. **Decisions live in dated notes; planning and implementation are separated.**
    Newer notes supersede older. Plan first (this document), then a Claude Code
    Plan-Mode session produces the build spec, then scripts.
12. **Distinguish verified from unverified claims, and present interpretations as
    hypotheses.** Do not import literature assumptions that have not been checked
    against this dataset; flag gaps honestly in the prose.

---

## 3. The justification architecture

Layer 0: Menegollo 2024 is the human prevalence/prognosis axis (mito programs as a
prevalent, prognostic, cell-of-origin-linked switch axis in human breast cancer);
myc_mouse is the experiment that closes the seeding-circularity gap from outside. The
logic is reciprocal validation - the human data motivate the axis, the mouse data
test it experimentally. This is AP0.

Connective point for the intro (from the ER/Myc/PGC-1a comparison): Menegollo's
mito-centred clusters separate along ER and Myc signalling, and ER-driven and
Myc-driven mitochondrial biogenesis are mechanistically distinct - ER works largely
through PGC-1a/NRF-1, Myc binds E-boxes directly and couples biogenesis to a
glycolytic/glutaminolytic anabolic program. Intro/discussion framing, not a new
analysis.

The empirical mito-focus justification (AP6) rests on two complementary tests:
- **Direction-agnostic magnitude:** an expression/dispersion-matched permutation null
  on mean |LFC| (bin all genes by mean expression x dispersion, sample matched random
  sets, locate MitoCarta in the null). Both up and down genes contribute, so a switch
  registers as large. Defeats the "mito genes are just high-expression housekeeping"
  confound.
- **Bidirectional omnibus (optional):** ROAST "mixed" mode and/or CAMERA (limma) -
  purpose-built for either-direction change and inter-gene correlation. The principled
  alternative to entropy.

---

## 4. Action points (the menu)

Status key: **[done]** computed and saved; **[reframe]** read/re-express existing
outputs; **[new]** new build (mostly library-dependent GSVA).

**AP0 - Intro framing (prose, not analysis).** Lead the justification with Menegollo
2024; position myc_mouse as the experimental test that closes the seeding-circularity
gap from outside. Turns "we studied mitochondria" from assertion into a
cross-validated rationale. **[reframe / prose]**

**AP1 - WT developmental baseline.** GSVA of the dev cell-state sets (`MG_*`) on Myc-
samples 6W vs 12W; mitoPPS on Myc- 6W vs 12W. Establishes the moving developmental
backdrop without which there is no "developmental change" for Myc to depend on.
Supports permission if a lineage/state shifts and mito priorities move. **[new GSVA;
mitoPPS exists]**

**AP2 - Which developmental state Myc+ resembles (core composition test).** GSVA of
`MG_*` on all samples (one run, log-scale); genotype effect (Myc+ vs Myc-) on the
scores, timepoint-controlled. The central permission test. Supports permission if Myc+
is enriched for a competent lineage (luminal-progenitor-like is the cell-of-origin
prior for MYC/basal-like disease; possibly TEB/proliferative). Points to selection if
instead Myc+ shows loss of a population carrying a stress/death signature. **[new
GSVA]**

**AP3 - Does the Myc effect track the developmental axis?** Interaction on the dev-set
GSVA scores and mitoPPS; correlate the Myc-vs-WT difference vector against the WT
developmental-shift vector (AP1). Supports permission if Myc amplifies the direction
development is already moving rather than imposing an orthogonal one. Weakest of the
three - gene-level interaction is weak, so AP1+AP2 carry the argument; AP3 is
corroboration. **[reframe + new vector corr]**

**AP4 - Is Myc endogenously part of the competent/pubertal program?** Myc dose is
constant (section 1), so there is no transgene-level dose/state coupling to test. Two
membership/co-variation tests: (i) is Myc a member of the competent-state dev
signatures (TEB/proliferative/pubertal `MG_*`); (ii) in the WT trajectory, does
endogenous MYC-pathway activity (Felsher GSVA) co-vary with the competent-state score.
If yes, the transgene amplifies a program competent cells natively run, with dose held
fixed - the cleaner mechanistic statement. **[new, light - membership lookup + one
GSVA corr]**

**AP5 - Divergence timing (GATE 1).** (i) Genotype difference at 6W alone - already
divergent at the earliest sampled point? (ii) Interaction direction across 6W->12W -
grows, flat, or shrinks? (iii) Project all samples onto stage-annotated dev sets to
place 6W/12W on the trajectory. Gating: if divergence *grows* across the window,
"developmental change licenses Myc" weakens and the title must move; if
already-large-and-flat at 6W, the soft "already diverged" claim holds. Script 11
Part 2 already computes interaction direction (negative = Myc effect weakens at 12W);
the 357-gene p<0.05 set and its directional bias exist (`interaction_sig_genes_p05.csv`);
the 6W contrast is in `interaction_results.rds`. **[reframe - reading existing outputs
+ the projection from AP2]**

**GATE 2 - Apoptosis selection-vs-remodelling readout.** Extract the PRO (25) and ANTI
(9) apoptosis group-shift t-tests from `interaction_gene_characterisation.rds` (script
11 Part 5B saves these as `interaction_by_geneset$apoptosis_pro / apoptosis_anti`);
confirm/compute both one-sample t-tests; interpret. Decides whether the selection force
is mito-apoptotic and module-wide, i.e. whether the mitochondrion is a decision point
(principle 5) or only a readout. Current evidence (apoptosis padj ~0.91 at RNA;
protein-level decline with RNA discordance) tilts toward **readout (metabolic)**, with
the decisional claim resting on protein + inference. **First, confirm the direction and
PRO/ANTI analyses used raw / Wald-based values, not shrunken interaction LFCs** - if
shrunken, recompute on raw before this gate decides anything (same issue as AP-CD;
tension T3). **[reframe - readRDS + interpret; code-read check first]**

**AP6 - Justify the mito focus empirically, by decomposition.** (i) directional fGSEA
on sub-pathways + MB fork-sets + dev/TF sets through the interaction ranking; (ii)
expression/dispersion-matched permutation null on mean |LFC| for the compartment-level
magnitude claim; (iii) the same against a comparator panel (proliferation, metabolism,
secretory, immune - already in `gene_sets_list.rds` as Hallmark) to show *preferential*
not just absolute alteration; (iv) optional ROAST-mixed/CAMERA omnibus. Report whether
the signal is magnitude or coordination; frame as licensing *focus*, not causal
primacy. Note: proliferation will likely also score high in a MYC model - expected; use
it to make the sharper point that the mito signal is comparable to/above proliferation
AND that the developmental/stage interaction is mito-specific where proliferation is
genotype-wide. **[new - cheap, off existing DESeq2 results]**

**AP7 - MB-fork validation (CENTREPIECE).** GSVA-project the mouse tumours onto the
MB1/MB2/MB12 fork signatures from the library; test whether Myc+ progression shifts
samples toward **MB2_UF** (the human LP-state, Myc/miRNA-signalling, biogenesis-high
fork). The direct experimental validation of the analytical paper and the strongest
answer to "why mitochondria" - mechanistic and cross-validated, not descriptive.
Falsifiable: Myc+ should land in MB2_UF specifically. Make this a named figure, not a
supplementary check. **[new - needs library snapshot + GSVA; sets confirmed present]**

**AP8 - Heterogeneity / survivor spread (optional).** Cross-sample variance/CV of the
mito-fork score, Myc+ vs Myc-. If selection spreads survivor states, CV is higher in
Myc+. The interpretable place a "diversity" idea belongs - done as variance, not
entropy. **[new, optional, last]**

**AP-CD - Cell death: both branches, neither canonical, on raw LFCs.** The two
analyses are complementary and both go into the final analysis. Report where they
converge and diverge.
- Branch 1 (main `09`, ported forward from `archive_main_pipeline/09_cell_death_pathway_summary.R`):
  binomial directional test on the curated `cell_death_genes_consolidated` list,
  metric ΔLFC = `myc_6W_log2FC - myc_12W_log2FC`, on **raw** LFCs. Directional-
  hypothesis lens. The raw combined dataframe already exists
  (`results/combined_df_annotated_raw.rds`); point the branch at it rather than the
  shrunken `combined_df_annotated.rds`, and regenerate the hypothesis outputs.
- Branch 2 (new `12`): fGSEA across the Tang 2024 15 RCD modalities (`data/cell-death/`)
  on the Wald z, per contrast including the interaction. Broad-enrichment lens.
- Shrinkage scope: fGSEA (branch 2) ranks on the Wald stat from the unshrunken fit, so
  it is unaffected. The correction is specific to branch 1's ΔLFC test and to any
  averaged-LFC visual. **[reframe branch 1 to raw; branch 2 done]**

**AP-mtDNA - mtDNA-encoded vs nuclear OXPHOS dissociation.** Feature the dissociation
as a named mitoPPS result: mtDNA-encoded subunits are prioritised at 12W (both
genotypes) while nuclear-encoded OXPHOS complexes (I, III, IV) are de-prioritised, and
mtDNA-encoded genes are not induced by Myc; Complex I shows the strongest effect. The
synthetic "mtDNA-encoded OXPHOS subunits" pathway already isolates this in the mitoPPS
build, so this is a plotting/interpretation task: show the mtDNA-vs-nuclear split, with
per-complex correlations in a supplementary. **[reframe - off existing mitoPPS output]**

**AP-abund - Abundance vs reprioritisation.** State the question explicitly and answer
it with the fGSEA-vs-mitoPPS contrast: fGSEA (rank-relative) reads more like overall
abundance shifts; mitoPPS (within-compartment, absolute) reads selective
reprioritisation among mito pathways. The honest answer is time-dependent - at 6W Myc
drives selective reprioritisation; by 12W the selectivity blurs toward uniform
biogenesis (an abundance-like change). Do not collapse to a single word; the 6W/12W
split is the finding. **[reframe - off existing fGSEA + mitoPPS]**

**AP-retention - Developmental retention index.** For pathways enriched in normal
development, compute the retention index NES_pos / NES_neg (proportional) and
NES_pos - NES_neg (absolute deviation) on the temporal contrasts; use the difference as
the biologically meaningful one. Carry the Felsher interpretation guardrail into the
prose: a stronger negative timecourse NES in Myc+ is **not** loss of Myc activity,
because the cross-sectional Myc enrichment is stable (slightly increasing) at both ages
- the correct reading is "developmental decline plus a stable Myc offset". Supporting
plots: slope of pos vs dev timecourse LFC (Felsher genes), cross-sectional stability
(myc_6W vs myc_12W LFC), interaction fGSEA. **[reframe - defines a metric on existing
fGSEA]**

---

## 5. Build structure: Block A / Block B, git, and the Step 1 / Step 2 split

**Two phases.**
- **Block A - broad exploratory.** Run the gates and the full AP menu across all
  relevant sets, methods, and contrasts. Both cell-death branches, the retention index,
  the abundance-vs-reprioritisation contrast, the mtDNA split, and the optional
  composite contrast all get computed here. Then review the outputs against the
  hypotheses and decide what the paper says.
- **Block B - focused for publication.** Re-code a narrower script set producing only
  the analyses and figures the narrative uses. The submission-quality pass, later - not
  this week.

**Git strategy.**
- `new-analysis` - stable trunk for the consolidated pipeline.
- `analysis-exploratory` - branched off `new-analysis` for Block A. All broad runs and
  the raw re-runs land here. Commit per verified analysis.
- `paper-figures` - branched off `analysis-exploratory` at the reviewed commit for
  Block B, after the Block A review. Only narrative-selected code moves forward.
- Merge `paper-figures` back to `new-analysis` when the draft figures are stable; keep
  `analysis-exploratory` as the exploratory record.
- Do not check out `main` in the working directory - use the `../myc_mouse_main`
  worktree or `git show main:<path>`.

**Library snapshot.** Snapshot the library into `data/genesets_from_library/`, pinned
to the `v1.0` tag/commit hash, with a provenance README. This depends on the `v1.0`
tag existing first - do the library close-out (merge `dev`->`main`, tag `v1.0`, push
`--follow-tags`, keep `dev`) before the snapshot. The library output is
**category-stratified**, and the structure is preserved on snapshot: copy the nine
mouse category GMTs into `data/genesets_from_library/by_category/` (01_mitocarta,
02_myc_signatures, 03_mammary_development, 04_metabolism, 05_proliferation,
06_tf_targets, 07_biogenesis_discrimination, 08_apoptosis,
09_biogenesis_apoptosis_intersections), the combined master
(`mammary_mito_myc_metab_v1_mouse.gmt`) at the top level, and `provenance_table.csv`
(which carries the per-set method tag). Do not copy the `human/` tree - this is a
mouse dataset and a wrong-species GMT is an easy mistake to make. `metabric_sets.rds`
(from the library `results/`) and `gray_chea_mito_tf_shortlist.csv` (from the library
`outputs/qc/`) come along too.

**Category GMTs are used per-category, not as one master NES.** Consuming the master
GMT in a single fGSEA run is exploratory-only (Block A omnibus scan). The publication
claims run fGSEA per category, because that is what makes AP6's "preferential, not
just absolute" contrast interpretable and keeps the multiple-testing scope legible.
Each category maps to the action points it serves:

- `01_mitocarta` -> AP6 (mito-focus decomposition), AP-mtDNA (mtDNA-vs-nuclear split),
  and the mitoPPS pathway partition.
- `02_myc_signatures` -> AP4 (Felsher membership/co-variation), AP-retention (the
  retention index and Felsher guardrail).
- `03_mammary_development` -> AP1, AP2, AP3, AP5 (the `MG_*` developmental sets, via
  GSVA and the dev-set projection).
- `04_metabolism`, `05_proliferation` -> AP6 comparator panel (show the mito signal is
  comparable-to/above proliferation, and that the stage interaction is mito-specific
  where proliferation is genotype-wide).
- `06_tf_targets` -> AP6 decomposition (dev/TF lanes through the interaction ranking).
- `07_biogenesis_discrimination` -> the Category 7 MYC_SPECIFIC / CORE / DEVELOPMENTAL
  partition on Figure 2.
- `08_apoptosis` -> AP-CD branch 2 (alongside the Tang `data/cell-death/` sets) and
  Gate 2.
- `09_biogenesis_apoptosis_intersections` -> the biogenesis-vs-apoptosis intersection
  claim (Supp/Figure 2 support).

Open design point for the build spec: confirm per-category vs master fGSEA scope and
the multiple-testing correction that follows (correct within-category, or pool across
categories). This is a scientific choice - the build spec should propose it explicitly
and surface it for approval, not assume it. The AP6 comparator panel already loaded
in-repo (`gene_sets_list.rds`, 89 sets incl. 50 Hallmark) remains the Hallmark source;
the `04_metabolism`/`05_proliferation` category GMTs are the mammary-specific
comparators - decide whether to use one, the other, or both.

**Set-design rationale is in `docs/library_reference/`.** The per-set provenance and
method tags are machine-readable in `provenance_table.csv`, but the *why* behind each
category - what "biogenesis discrimination" operationalises, the Category 7 algebra,
which developmental atlases feed the `MG_*` sets, why the pass-1 TF roster is what it
is - lives in the library reference docs snapshotted into `docs/library_reference/`
(read-only; see that folder's README). Consult them when a set's meaning or method
matters for a figure claim or a method-section sentence. They are not a build spec:
sets are consumed from the snapshot, never rebuilt from these docs.

**Step 1 / Step 2 pipeline split.**
- Step 1 = data load + DESeq2 (basic contrasts + interaction) + QC (basic and
  interaction). Produces the DESeq2 objects and contrast results.
- Step 2 = pathway analyses (fGSEA, mitoPPS, GSVA, cell death, characterisation), built
  on Step 1 outputs + the library GMTs.
- Keep both the interaction and the group-based analyses. Route gene-set loading
  through the snapshot; do not rebuild sets in `01_load_data.R`. The split lets Step 2
  be swapped without re-running DESeq.

**Optional composite contrast.** 12W_pos vs 6W_neg can be added in Block A, but label
it correctly as a composite state-to-state contrast (developmental change + Myc effect
at 6W + interaction), not a clean effect, and interpret it only alongside the clean
contrasts. Probably adds little; exploratory-only and a cut candidate.

---

## 6. Execution logic: gates first (inside Block A)

The APs are a dependency tree, not a flat queue. Run the two cheap gates first; their
results decide which later APs to build and in what form.
- Gates draw only on existing myc_mouse outputs (interaction model, contrasts,
  `interaction_gene_characterisation.rds`) - they need none of the new library GSVA
  infrastructure, so they run essentially first. Confirm the Gate 2 inputs are raw /
  Wald-based before trusting them (T3).
- AP7 is load-bearing and cheap once GSVA exists - run it early-and-high regardless of
  the gates.
- AP1-4 and AP6 are shaped by what the gates return. If Gate 2 (PRO/ANTI) is null, AP4
  and the apoptosis framing shrink and you lean on the metabolic-readout
  interpretation; if Gate 1 shows growing divergence, AP1 becomes less central and the
  title moves off "licenses." Building AP1-4 before the gates risks building the wrong
  version.
- AP8 only if the selection arm survives.

Natural sequence: gates (existing data) -> AP7 + permission/justification build
(library-dependent) -> AP8 (optional). All inside Block A; the figure-selected subset
is re-coded in Block B.

---

## 7. Figure architecture (locks after Day 1 gates)

Provisional - if Gate 1 shows divergence growing *within* the window rather than
pre-diverged, Figure 1's framing changes.
- **Figure 1 - "Mitochondrial programs are a principal axis engaged during Myc
  tumourigenesis":** MB-fork projection (AP7), preferential-alteration panel (AP6), WT
  developmental mito shift (AP1).
- **Figure 2 - "Divergent Myc+/Myc- trajectories and the developmental-vs-oncogenic
  mito partition":** lineage-composition GSVA (AP2), Category 7 discrimination
  (MYC_SPECIFIC / CORE / DEVELOPMENTAL), divergence timing (AP5).
- **Supp 1:** apoptosis PRO/ANTI detail and the selection-vs-remodelling readout
  (Gate 2, script 11), plus the two-branch cell-death convergence (main 09 raw + new 12).
- **Supp 2:** matched-null permutation, comparator panel, QC, provenance.
- **Supp 3 (small):** the mtDNA-vs-nuclear OXPHOS dissociation and per-complex mitoPPS
  split (AP-mtDNA), plus the abundance-vs-reprioritisation fGSEA/mitoPPS panel
  (AP-abund). A supporting point, not a figure claim.

---

## 8. The 5-day timeline (de-risked)

**Day 1 - Decision layer + setup (existing-data, lowest risk).**
- Library close-out (in `mammary_geneset_library`): merge dev->main, tag v1.0, push
  --follow-tags, resync dev (0.5h).
- Create `analysis-exploratory` off `new-analysis`; snapshot library GMTs + provenance
  into `data/genesets_from_library/` pinned to v1.0, with a provenance README (0.75h).
- Code-read check: confirm the script-11 direction and PRO/ANTI analyses use raw /
  Wald-based values, not shrunken interaction LFCs (feeds Gate 2). No execution (0.25h).
- GATE 1 (AP5 divergence timing): read script 11 Part 2 outputs + the 6W contrast; add
  the dev-set projection placeholder (1.0h).
- GATE 2 (PRO/ANTI apoptosis): readRDS `interaction_gene_characterisation.rds`, extract
  apoptosis_pro/anti, confirm both t-tests, interpret (0.5h).
- Dated decision note: what the gates returned, which framing/APs survive, figure lock
  (0.5h).
*End: permission-vs-divergence fork resolved; selection arm confirmed or demoted;
figure plan locked; Block A branch + snapshot in place.*

**Day 2 - Cell death (both branches, raw) + GSVA pipeline.**
- Port main 09 forward onto the interaction model and point it at
  `combined_df_annotated_raw.rds`; regenerate the binomial hypothesis outputs on raw;
  run new 12 as-is; assemble the convergence read (1.0h - less than a from-scratch
  re-run because the raw combined df already exists).
- Verify the existing `outputs/heatmaps_int/*_raw.pdf` pathway heatmaps are the figure
  source (the raw variants are already generated) (0.25h).
- Build the log-scale GSVA scoring script; validate in Positron. Guard the traps:
  log-scale (not mitoPPS linear); all samples one run (2.5h).
*End: two-branch cell-death convergence recorded on raw; validated set x sample GSVA
matrix.*
*Wildcard: if cell death drags, time-box it - it is now a mechanical "run both on raw"
task, not a judgement call, and it is Supp-level, not Figure 1/2.*

**Day 3 - Harvest the GSVA-dependent results (analytical core).**
- AP7: project tumours onto MB1/MB2/MB12; test shift toward MB2_UF (1.5h).
- AP1 + AP2: WT developmental shift; lineage-composition genotype effect (1.5h).
- AP6 part 1: decomposed directional fGSEA on sub-pathways, MB fork-sets, dev/TF sets
  through the interaction ranking (1.5h).
*End: validation, composition, decomposed enrichment - computed.*

**Day 4 - Finish justification + Figure 1.**
- AP6 part 2: expression-matched permutation null + comparator panel (2.0h).
- Assemble draft Figure 1 (2.0h). Begin Supp panels (0.5h).
*End: focus justification complete; Figure 1 in draft.*

**Day 5 - Figure 2 + supplementaries + draft text.**
- Assemble draft Figure 2 (2.0h). Supp 1 + Supp 2 + Supp 3 to draft, folding in
  AP-mtDNA, AP-abund, AP-retention (1.5h).
- Draft section text: AP0 intro framing, results prose with gaps flagged, figure
  legends (1.0h).
*End: Figures 1-2 + 3 supplementaries in draft; section text drafted with gaps marked.*

---

## 9. Cut-lines (drop first if behind; add back in reverse if ahead)

1. AP3 and AP4 fine detail.
2. AP8 (cross-sample CV).
3. ROAST/CAMERA omnibus (the matched-null already carries the bidirectional point).
4. Comparator-panel breadth (keep proliferation only).
5. Supp 3 breadth (keep the mtDNA split; drop the per-complex correlation detail first).
6. All figure polish.

The cell-death branch-1 raw re-run is not a cut-line - it is a correction that feeds
Gate 2, so it runs regardless. If time collapses, drop Figure 2 to a single
discrimination panel and let it slip, rather than compressing the gates (AP5/Gate 2) or
the validation (AP7) - those are what the argument stands on.

---

## 10. Housekeeping flags

Resolved (recorded for provenance; no action):
- **Worktree phantom:** the old `../myc_mouse_new` reference never existed on disk;
  CLAUDE.md has been corrected. `../myc_mouse_main` holds `main`.
- **MitoCarta file:** `08_mitoPPS_analysis.R` loads `Mouse.MitoCarta3.0.xls` (dotted,
  new-only) from Sheet 4. `Mouse_MitoCarta3_0.xls` (underscore) is the older shared
  variant, different MD5. Load by exact name; never glob `Mouse*MitoCarta*`.
- **Count matrix:** `01_load_data.R` loads `FULL.DAT.csv` (sample-code headers).
  `FULL.DAT_copy.csv` is the same numeric data with group-label headers and different
  column order - not a byte-duplicate.

Open (fix before the writeup):
- **`NES_paradox_explanation.md` mislabels the model as "liver"** - it is the mammary
  6W/12W Myc+/Myc- data. Correct the file.
- **mitoPPS outlier:** one 12W Myc- sample is an outlier on the mitoPathway (abundance)
  score; a 6W-pos sample is a PCA outlier - not the same sample. Identify the driver,
  decide on exclusion, record the call.
- **Stray `output/` (singular) directory** holding only `gprofiler/` sits beside the
  real `outputs/`. Tidy or gitignore to avoid confusion.
- **`docs/branch_manifest.md` is dated 2026-04-22** and still frames the two-pipeline
  consolidation as ongoing; this plan's section 5 and CLAUDE.md's "current phase"
  supersede that framing. Refresh the manifest when convenient.

---

## 11. Current-state ledger

Already computed and on disk (reframe/read, do not rebuild): the DESeq2 interaction
model and seven contrasts (`interaction_results.rds`, `group_results.rds`); temporal
and cross-sectional fGSEA (`fgsea_results.rds`, `fgsea_xs_results.rds`); mitoPPS scores
(`mitopps_scores.rds`); fGSEA-vs-mitoPPS comparison (script 09/10); interaction gene
characterisation incl. PRO/ANTI apoptosis (`interaction_gene_characterisation.rds`);
new-12 cell-death fGSEA (`cell_death_fgsea.rds`); the binomial-branch hypothesis
outputs (`cell_death_hypothesis_results.csv` + `results/figures/`); the raw combined df
(`combined_df_annotated_raw.rds`); raw + shrunk per-pathway heatmaps
(`outputs/heatmaps_int/`); the AP6 comparator panel (`gene_sets_list.rds`); the
biomaRt ortholog table (`ortholog_table.rds`).

To build: the log-scale GSVA scoring script and all GSVA-dependent APs (AP1, AP2, AP7,
parts of AP3/AP4); the AP6 permutation null; the raw wiring of cell-death branch 1; the
figure assembly. Nearly everything else is reframe/read.

---

## 12. Open decisions

- **Cell-death canonical:** resolved - use both branches, run branch 1 on raw. Not a
  choice of method (AP-CD).
- **Myc dose ramp through puberty:** answered - no, stable at both stages (section 1).
- **Title verbs:** decide after the gates. The stable Myc dose makes "licenses" /
  "sets the stage for" more defensible (constant driver, moving substrate), but keep it
  gated on Gate 1.
- **How divergent the two cell-death analyses are:** judged on the raw re-run, so it
  cannot be settled until Day 2.
- **Where the `01` gene-set tidy belongs:** routing the legacy `mitocarta_pathways.csv`
  / `.gmx` loads through the library snapshot is a real code change touching fGSEA
  inputs. Decide whether it goes in Block A now or is deferred.
- **fGSEA scope - per-category vs master (section 5):** per-category is the default for
  the publication claims, master for exploratory omnibus only. Confirm this, and the
  multiple-testing correction that follows (within-category vs pooled). Also decide
  which comparators feed AP6 - the in-repo Hallmark panel, the mammary-specific
  `04_metabolism`/`05_proliferation` category GMTs, or both.

---

## 13. Tensions to reconcile

None are hard data contradictions; they are conflicts of emphasis or sequencing.

**T1 - cell-death framing.** Earlier framing asked which cell-death analysis is
canonical; that premise is rejected (run both). Recorded so it does not creep back into
the writeup.

**T2 - cell death central vs supplementary.** The narrative arc wants apoptotic
sensitivity central to the ER-vs-Myc story, but current RNA evidence (padj ~0.91) points
to readout, not decision (principle 5). Decide cell death's figure weight *after* the
raw re-run, not on the shrunken outputs. If the RNA signal stays weak, keep the
apoptotic-decision claim resting on protein + inference and hold cell death at Supp 1.

**T3 - Gate 2 inputs and shrinkage.** Gate 2 reads
`interaction_gene_characterisation.rds`. If its direction/PRO-ANTI analyses used shrunken
interaction LFCs, recompute on raw before the gate decides anything. First-step
code-read check on Day 1.

**T4 - abundance vs reprioritisation.** "Abundance" and "selective reprioritisation
blurring to uniform" are reconciled by making the answer time-dependent (AP-abund); do
not over-commit the prose.

**T5 - full narrative arc vs draft scope.** The full arc runs Timeline -> ER/Myc
deconvolution -> PGC-1a as an ER proxy (experimental) -> human/TNBC. This draft
(Fig 1-2 + supp) covers the transcriptomic timeline and the ER/Myc deconvolution. The
PGC-1a-proxy experimental arm and the human/TNBC extension are downstream (other paper
sections or the Hannon human analysis), not this draft. Named so they are not lost and
the draft is not over-scoped.
