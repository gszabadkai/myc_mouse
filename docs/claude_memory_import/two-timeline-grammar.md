---
name: two-timeline-grammar
description: "Figure 2's alternatives all share one plane -- WT 6>12W on x, Myc+ 6>12W on y, identity line = development alone -- and the vertical drop from that line IS the interaction, exactly"
metadata: 
  node_type: memory
  type: project
  originSessionId: 71ffe8cc-b22d-412b-9e40-f22881e23306
  modified: 2026-08-08T22:51:53.416Z
---

The grammar the author chose (2026-08-09) for the Figure 2 alternatives, and the arithmetic that
makes it more than a nice picture.

**THE PLANE.** `x = 6>12W_wt`, `y = 6>12W_myc`, dashed identity line = *development alone*.
Declared once as **`two_timeline_base()`** in `figures/panels/_panel_common.R`; four panels use it
(`fig2_two_timelines_mito_alt.R`, `fig2_death_two_timelines_alt.R`, `fig2_puma_chain_alt.R`,
`figS2_p53_arm_alt.R`).

**WHY IT IS EXACT.** Script 40's ruler satisfies `c_tp = c_tn + c_int` to floating point (1.2e-16
on content, **bit-identical** on priority; asserted in the panels). So **a vertical drop from the
diagonal IS the interaction** — not an approximation of it. That is what turns "OXPHOS repression
is more than the overall MYC effect decline" into geometry instead of a sentence. `coord_equal()`
is load-bearing: both axes are the same quantity, so equal scaling is honest *and* it pins the
diagonal at exactly 45°, which is what lets its label lie along it without the panel-aspect
arithmetic Fig. 2G needed.

**THE NUMBERS THE PLANE SHOWS.**
- **OXPHOS is additive**: subunits **−0.255** in development, **−0.150** more from Myc, **−0.405**
  under Myc. Assembly factors of the same complexes: **+0.001** in development, then only the
  global term. The compartment as a whole *rises* in development (median **+0.041**, 72 % above
  zero) and falls under Myc (median −0.135); **93 %** of the 143 sit below the diagonal, which is
  the global fade, not a mitochondrial result.
- **The chain is an interaction, not a co-fall** (group means of z, per animal):

  | | OXPHOS | Foxo3 | Bbc3 | Bcl2l1 | PUMA:Bcl-xL |
  |---|---|---|---|---|---|
  | 6W_wt | +0.07 | −0.90 | −0.12 | +0.74 | −0.55 |
  | 12W_wt | −0.78 | **+0.89** | +0.18 | +0.39 | −0.16 |
  | 6W_myc | +0.95 | +0.01 | +0.84 | −0.71 | +0.94 |
  | 12W_myc | −0.25 | −0.00 | **−0.90** | −0.42 | −0.23 |

  Read **down the timelines**: WT — OXPHOS falls, **Foxo3 rises**, PUMA untouched. Myc — OXPHOS
  falls from higher, **Foxo3 fails to rise**, **PUMA collapses**. `Bcl2l1` follows neither, so the
  ratio moves with its numerator. I first wrote this panel claiming "the 12W_myc block falls
  together"; **the data contradicted that and the panel is now asserted against inversion.**
- **Bmf is the panel's own control**: +1.02 WT / +0.84 Myc, significant in both, sitting **on** the
  diagonal — the one death gene the window itself moves.

**WHAT THE 2026-08-09 REVISION CHANGED.** The plane survives in **three** panels (2F, 2G, S2D
alternatives); **2H+I left it**. The author's verdict on the per-animal heatmap was *"the heatmap
does not work, the individual data are too much variable"* -- right about the data, wrong as a
display: at n = 6 the reader had to average six noisy tiles by eye in four places. It is now
distributions in Fig. 1E's idiom (box + every animal), in the order the chain runs, with the two
members of the ratio beside them as a 2 x 4 median heat map. The numbers it draws, one OLS
instrument throughout (`lm(value ~ timepoint)` within genotype; a z-score is linear so these are the
p-values of the untransformed quantity):

| | wild type 6>12W | Myc+ 6>12W | difference |
|---|---|---|---|
| OXPHOS mitoPPS | -0.11, p 0.108 | **-0.15, p 0.027** | p 0.615 |
| `Foxo3` | **+0.46, p 0.0044** | -0.00, p 0.982 | **p 0.013** |
| PUMA:Bcl-xL | +0.17, p 0.517 | **-0.50, p 0.019** | **p 0.042** |

Read across: the respiratory arm falls on BOTH timelines, `Foxo3` rises ONLY in the wild type, the
priming ratio collapses ONLY under Myc. The **difference** column is the only one `batch = timepoint`
leaves clean, and it is drawn as a second-tier bracket from 1.5 to 3.5 -- midpoint of one lower
bracket to midpoint of the other, which is what a difference of differences is.

**How to apply:** any new Figure-2 panel starts from `two_timeline_base()`, never a re-typed
diagonal. Read a point's *vertical distance from the line*, not its position along the axes — both
axes are batch-confounded ([[batch-equals-timepoint]]) and only the interaction is batch-clean.

Related: [[foxo3-puma-pair-and-the-nonlinear-dedosing]], [[background-vs-myc-is-additive]],
[[main-figure-panels-must-carry-the-argument]], [[figure-phase-architecture]].
