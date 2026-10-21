---
name: hypothesis-generating-not-failed-tests
description: "The author's reframe: a number that fails its null can still be a good hypothesis-generator -- present the competing readings side by side with the experiment that adjudicates each, rather than reporting the test as dead"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 71ffe8cc-b22d-412b-9e40-f22881e23306
  modified: 2026-07-27T14:55:28.918Z
---

The author, 2026-07-27, on the OXPHOS <-> PUMA-priming interaction after I had reported it as a
bounded near-miss: *"Understand the stats. But we are generating hypothesis. Can those numbers
(ranking exercise as you call it) actually generate the hypothesis we prove in the in vitro exps? It
might indicate that indeed OXPHOS is critical ... it might indicate that it is the TEB-Ductal axis
... and it might indicate that controlling OXPHOS is a trade-off for proliferation and PGC1 only
increases proliferation on a Myc background if the cell death is blocked? These kind of panels I'm
looking for on the figure to set the stage for the in vitro data."*

**This is the right frame for this corpus and it changes panel DESIGN, not statistics.** Every panel
I had built asked *"does this survive its null?"*, so each led with a hedge. A panel that asks *"what
does this nominate for testing, and what would settle it?"* uses the same numbers and reads as a
result. **The bioinformatics figure's job is to set up the perturbation experiments, not to stand
alone.**

**How to apply.**
- When a scan or coupling fails significance, do not stop at "negative". Enumerate the **competing
  readings the data licenses**, and give each the **experiment that adjudicates it** and whether it is
  done or open. Mark them as rankings on the panel's face, not only in the legend.
- **Report the reading that looks like a nuisance control.** The in-vivo proliferation negative had
  been filed as a control; in fact it is the *negative limb of the sign reversal* the cells complete
  (mouse: death intact, respiration costs growth; MYAZ + Bcl-xL: death blocked, PGC-1a buys growth).
  **Same intervention, opposite sign, switch = apoptotic competence.** That pairing is the hand-off.
- **Keep the rival hypothesis in the figure**, even when it is not separable. TEB-vs-respiration is
  not resolvable at n=24; showing both is more honest and more useful than picking one.
- Do not let this soften the limits: every row still carries its control and its percentile, and
  causation still belongs to the perturbations ([[prespecify-from-the-perturbations]],
  [[dose-vs-competence]]).

Full table with evidence, controls, limits and experiments: **section 0.10** of
`docs/2026-07-26_introduction_alignment_and_the_question.md`; rendered as **Figure 2F**
(`figures/figure2_developmental_window.R`). Related: [[figure-phase-architecture]],
[[sample-level-couplings-need-a-null]], [[correlation-ceiling]].
