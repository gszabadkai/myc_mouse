---
name: use-the-defined-scheme-dont-invent-one
description: "Figure panels must reuse the project's declared palette, labels and layer idiom; a per-panel variation reads as an error, not a choice"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 71ffe8cc-b22d-412b-9e40-f22881e23306
  modified: 2026-07-30T15:58:17.904Z
---

**2026-07-30.** Twice in one session the author rejected panel encodings that were *close to* the
house style but not it. On figS1C: *"labels and colours do not fit with the generic scheme. use as
it was defined."* The offences were small — box fills at `alpha = 0.5` instead of the established
0.28, points with a grey-25 stroke instead of solid group-coloured points, and rotated x-axis group
labels where the project's panels blank the x axis and carry a group key instead. Earlier the same
day: one sample palette fixed project-wide (Okabe-Ito, hue = genotype, lightness = age), one
diverging ramp (`ms_diverging`), one contrast vocabulary (`myc_6W`/`myc_12W`/`6>12W_wt`/`6>12W_myc`),
one spelling rule (mouse **Myc**, not MYC, except genuinely human objects).

**Why.** Across a figure set the reader learns the encoding once. A panel that varies it spends the
reader's attention on decoding rather than on the data, and a near-miss is worse than an obvious
difference because it reads as a mistake. The author is composing these into one manuscript figure,
so consistency is a property of the deliverable, not a nicety.

**How to apply.** Before writing a new panel, read `figures/panels/_panel_common.R` and
`figures/theme_myc.R` for the declared palettes and helpers, and grep `figures/*.R` for an existing
panel of the same *kind* (per-group distribution, contrast lollipop, stacked composition) and copy
its layer idiom verbatim — `fig01_mito_content.R:118-170` is the reference for a per-group boxplot
with points. If something genuinely needs a new encoding, declare it in `_panel_common.R` with the
reason, never inline in one panel. `figures/panels/PANELS.md` records the fixed schemes. Related:
[[panels-follow-the-narrative]], [[figure-phase-architecture]].

**`priority` -> `mitoPPS` (author, 2026-08-09, "for clarity"), project-wide in `figures/panels/`.**
The quantity is called **mitoPPS** wherever a reader meets it -- drawn labels and legend blocks,
including "the mitoPPS ruler" where the text said "the priority ruler". Three exceptions, each for a
reason: **data** (`ruler == "priority"`, `$comparator_priority`, `comp$priority` are column values
and object fields in the saved rds -- renaming breaks lookups silently); **the phenomenon**
(reprioritisation, a pathway being demoted -- that is what happens, not what is measured); and **the
definition**, which turns circular once the score is called mitoPPS and is reworded instead ("a
relative share of the compartment's resources"). `figures/*.R`, the superseded exploratory layer,
is NOT swept -- it bakes titles into the plots and uses the word as running prose.

**Significance ink (`sig_cols`) now covers TEXT, not just brackets**: a bracket's p-value label, and
a gene NAME on a per-gene scatter whose padj clears the threshold. The rule behind it still binds --
where the geometry carries the effect size, significance never becomes a mark's fill.
