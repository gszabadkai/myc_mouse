---
name: main-figure-panels-must-carry-the-argument
description: "The author's 2026-08-09 criticism: the built Figure 2 panels were 'overly correct' and read as supplementary -- a main panel must deliver the message alone, because people do not read the text"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 71ffe8cc-b22d-412b-9e40-f22881e23306
  modified: 2026-08-08T22:52:11.429Z
---

After all fourteen commissioned panels were built and signed off one at a time, the author read the
assembled Figure 2 and said the panels were **hard to understand** — *"overly correct, but they look
more like supplementary panels"* — and set the standard: **"the figures have to talk by themselves,
people do not often read the text."**

**Why:** building panel-by-panel against one sentence each optimises every panel for *its own
sentence* and nothing for the argument. Each was defensible; none carried the section. The section's
message is a chain (development strips OXPHOS → Myc strips it again → the PUMA:Bcl-xL trigger is
lost → the adult gland cannot die) and no single-sentence panel showed a chain.

**How to apply:**

- **Ask what the section claims before asking what the sentence says.** A main-figure panel earns
  its place by advancing the argument; a panel that only survives its own caveats belongs in the
  supplement.
- **Prefer one grammar across a figure** so the reader learns the axes once ([[two-timeline-grammar]]).
  Variety per panel is a cost, not a feature.
- **Make panels talk WITHOUT breaking the publication rule.** `theme_panel()` still blanks
  title/subtitle/caption. What carries meaning is **naming drawn elements** — `development alone`
  along a diagonal, `below the line = lost under Myc as well` in a corner, `global rescaling
  (slope 0.49)` on an imposed line. Naming a drawn element is what a key does; it is not prose.
- **Legibility decisions are forced by the data, not taste, and belong on the record.** From this
  round: clip a per-animal fill at 2 SD when one outlier would leave every tile pale; print the
  value in a group-mean cell, because a mean of six z-scores is necessarily paler than its animals
  and would otherwise be the weakest mark exactly where the claim lives; rotate that number when
  the tile is 2.5 mm wide.
- **Alternatives are additive.** The author asked to KEEP every existing panel and build
  alternatives beside them, choosing per slot. `panels_to_pdf.R` writes a separate
  `Figure2_alternatives.pdf` so the committed circulation documents do not change.

**Six things the 2026-08-09 revision settled, all reusable:**

- **A key is laid out like the axis it repeats.** Fill = the x axis -> horizontal bar; fill = the y
  axis -> vertical bar. Put it inside, in a region the script ASSERTS is empty.
- **`theme_classic` sets `axis.text` to `rel(0.8)` of `base_size`**, so a legend sized from
  `base_size` prints LARGER than the axis it repeats. Set key text explicitly (4.8 at base 6).
- **Zero on `ms_diverging` is `#FAFAFA`.** Points carrying that fill must be outlined circles
  (`shape = 21`) or everything near the median vanishes into the page.
- **Written-out column labels are expensive**: "12W_myc" rotated under a 6 mm column costs ~10 mm of
  height, and patchwork then pads the neighbouring panel's blank axis to match. A strip of the
  sample colours costs 2 mm and keys off the legend already on the figure.
- **`coord_fixed()` on a small heat map** or its cells stretch to whatever height the panel beside
  them needs.
- **Never draw one half of an interaction without its statistic.** See
  [[foxo3-puma-pair-and-the-nonlinear-dedosing]] point 6.

Related: [[panels-follow-the-narrative]], [[use-the-defined-scheme-dont-invent-one]],
[[hypothesis-generating-not-failed-tests]], [[two-timeline-grammar]].
