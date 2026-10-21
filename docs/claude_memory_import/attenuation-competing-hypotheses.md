---
name: attenuation-competing-hypotheses
description: "Four candidate mechanisms for the 6W->12W Myc attenuation, which scans are already done and negative, and that the whole question is blocked on the pending MYC blot"
metadata: 
  node_type: memory
  type: project
  originSessionId: 71ffe8cc-b22d-412b-9e40-f22881e23306
  modified: 2026-07-25T15:36:31.130Z
---

**RESOLVED 2026-07-25 by the MYC blot: protein falls ~50% at constant transcript, so H3 wins outright and the MMTV attenuation is DOSE — see [[dose-vs-competence]].** H2/H4 survive only as answers to the *competence* question (why the same MYC kills less), not the *dose* question. This memory is kept for the four-way framing and, above all, for the negative scans below — they remain the authority on what was tested.

Why the Myc effect halves 6W->12W was open, and "the programme decays" ([[background-vs-myc-is-additive]]) is a description, not a mechanism. Four hypotheses, written up 2026-07-25 in `docs/2026-07-25_attenuation_competing_hypotheses.md`:

- **H1 autonomous decay** — unfalsifiable as stated; the default.
- **H2 background antagonism** (the author's model — the maturing gland holds Myc's targets down): **H2a selective** blockade FAILS (attenuation is global, no programme in excess, per-tier slopes uniform 0.53-0.83); **H2b global competition** is the strong form and fits the uniform x0.55 rescaling; **H2c lower set-point** fails on headroom.
- **H3 reduced effective MYC activity at constant message** — front-runner on parsimony; script 41 excluded only a *transcriptional* change in the degradation machinery.
- **H4 compositional dilution** — under-weighted; now has direct marker support.

**Scans ALREADY DONE and negative (do not redo, cite this):** every E-box binder is flat (Usf1/2, Bhlhe40/41, Tfe3/Tfeb/Mitf, the whole MAX/MXD/MNT/MGA/MLX network) — no interaction anywhere; a 111-gene NR/corepressor/Polycomb/DREAM/repressive-E2F scan gives 6 background movers at FDR<0.05, **zero** with a significant interaction, and two point the wrong way (Esrrb an activator UP, E2f7 a repressive E2F DOWN).

**Two corrections worth keeping.** (1) `Id4` +1.53 (padj 0.0013) looked like the one HLH-sequestrator lead but is **composition**: it rises inside a myoepithelial signature (Oxtr +1.72, padj 2.7e-4; Myh11 +0.85; Trp63 +1.17) while Id1/2/3 are flat — and ID proteins sequester E-proteins, not bHLH-LZ MYC:MAX. That Oxtr signature is the first direct evidence the background move is a differentiation shift, which is what H4 needs. (2) The author's `Mlxip`/`Mlxipl` observation is the ONE named H2b candidate, but the sign is wrong for competition, the interaction is null (padj 0.83/1.00 — figS8's brackets are unadjusted), and Txnip is flat. It survives only via **MLX release** (MLX is shared with MXD1/MXD4/MNT), which is a protein-stoichiometry route that RNA-seq cannot test.

**Issue #5 did NOT test this** — script 30's axes were all *activating* (luminal, ESRRA/NRF1/GABPA, proliferation). A repressive axis has never been entered.

**Why:** the author asked whether background repression, not autonomous decay, explains the attenuation; the negative scans bound the answer and would otherwise be lost with the session.

**How to apply:** do not build script 42 until the MYC blot lands — protein down = H3 wins and the rest is moot; protein flat = T1 (de-confounded split, reusing script 40 PART C) is the priority. Never test H2 with `cor(WT-time, delta-Myc-effect)`: it is in the artifact ledger, algebraically forced ([[background-vs-myc-is-additive]]). And check any new mechanism against these scans before asserting it ([[dont-assert-a-mechanism-without-checking-it-can-work]]).
