---
date: 2026-10-17
tags: [project/myc_mouse, nat-metab, figure-2, colour, proposal]
status: ADOPTED (author, 2026-10-19) -- Option A, with bars at 100%; in _style.R (construct_cols) and FIGURE_RULES.md. Extended 2026-10-23: PGC-1α + BCL-xL #00573F, and the no-Dox tints (construct_nodox); plain or hatched open
relates-to:
  - figures/natmetab/proposals/p2_construct_colours.R   (draws the page)
  - outputs/natmetab/proposals/p2_construct_colours.pdf (three pages: in use, Option A, Option B)
  - figures/natmetab/proposals/p2_bar_fills.R           (draws the bar-fill page, 2026-10-19)
  - outputs/natmetab/proposals/p2_bar_fills.pdf         (bars at 28%, 60%, 100%; 100% under three CVDs)
  - figures/natmetab/proposals/p3_pgc1a_bclxl_colours.R (the PGC-1α + BCL-xL candidates, 2026-10-23)
  - figures/natmetab/proposals/p23_nodox_tints.R        (the recorded palette with the no-Dox tints, plain and hatched)
  - figures/natmetab/_style.R, figures/natmetab/FIGURE_RULES.md (the colours in use)
---

# Fig. 2: colours for the cell-line constructs (adopted 2026-10-19)

## Decision (author, 2026-10-23): PGC-1α + BCL-xL, and the no-Dox tints

**PGC-1α + BCL-xL (Figs 3-4) is `#00573F`**, the biogenesis green darkened to 55%: candidate (b)
of five on `p3_pgc1a_bclxl_colours.pdf` ("yes, b is good"). Under deuteranopia, protanopia and
tritanopia it stays at least CIELAB dE 25 from PGC-1α `#009E73` and 31 from BCL-xL black, the most
even balance of the candidates. The others: `#005F45` (x 0.60, 22 from PGC-1α), `#004734` (x 0.45,
25 from BCL-xL), Paul Tol's green `#117733` (17 from PGC-1α under tritanopia), and the PGC-1α green
with a black hatch.

**No Dox.** The author uses `#B8E4D8`, which is the PGC-1α green at 28% over white (Option A's tint
column below), for PGC-1α's no-Dox bar on some bar graphs. The same rule gives Ndi1 / LbNOX matrix
`#F1D9E6` and LbNOX cytosol `#FBF7CA`; all three are `construct_nodox` in `_style.R`. No tint is
set for EV, BCL-xL or PGC-1α + BCL-xL. Prism: the no-Dox bar takes the tint, with the same black
0.9-pt border, white symbols and upward error bars.

What the tints cost, CIELAB dE under deuteranopia (computed 2026-10-23; `p23_nodox_tints.pdf` has
the table under all four visions):
- PGC-1α no Dox against EV `#CCCCCC`: **4** (17 in normal vision). Ndi1 no Dox against EV: **8**.
- PGC-1α no Dox against Ndi1 no Dox: **4**.
- LbNOX cytosol no Dox against the white page: 7 under tritanopia, so its black border carries it.
- No tint of the green escapes EV under deuteranopia (28-60%: dE 1-10), because EV's lightness lies
  between the tint and the green.

So a tint is safe only where it does not share a graph with EV or another tint. **Open: plain, or
hatched in the full colour** (0.5 pt), which the pattern separates where lightness cannot (page 2
of `p23_nodox_tints.pdf`). Also open: where Dox switches MYC (iMMEC rtTA-MYC), the tint marks MYC
off, an exception to "MYC state takes no colour" below; the axis still names the state.

## Decision (author, 2026-10-19)

**Option A, with bars filled at 100%.** Making the MYAZ PGC-1α bar charts, the author found the 28%
fills of the proposal "almost invisible, the fill colour is too light": the 28% rule is Fig. 1's
box rule, where the dots carry the hue; in a bar the fill is the only mark. Shown the bars at 28%,
60% and 100% (`p2_bar_fills.pdf`): "the 100% bar fill, border and symbols suggestion looks OK".

**Prism settings for a bar chart:**

| construct | bar fill |
|---|---|
| EV | `#CCCCCC` |
| BCL-xL | `#000000` |
| PGC-1α, NRF1 | `#009E73` |
| Ndi1, LbNOX matrix | `#CC79A7` |
| LbNOX cytosol | `#F0E442` |

- bar border black, 0.9 pt (the Prism format's every line);
- symbols: circles, white fill, black border (they read on every fill, black and yellow included);
- error bars upward only, so none hides inside a dark bar;
- MYC state on the axis, no colour.

The 60% fills, for the record: EV `#E0E0E0`, BCL-xL `#666666`, PGC-1α `#66C5AB`, Ndi1 `#E0AFCA`,
LbNOX cytosol `#F6EF8E` -- readable, but EV stays faint and the hues separate less under CVD.
Box plots, if Fig. 2 has any, keep Fig. 1's rule (28% fill, coloured dots).

In the code: `construct_cols` in `figures/natmetab/_style.R` (ASCII keys; labels are set in Prism
and Illustrator), with a guard that no construct takes a mouse group's colour; the rules in
`figures/natmetab/FIGURE_RULES.md`, "The cell-line constructs".

---

The proposal as written on 2026-10-17 follows.

**The author's question (2026-10-17):** "can you recall what is the current colour scales, what
colours are being used? I'm working on Fig.2 where groups are different from the mouse so I need to
use different colours. these should adhere to colourblind schemes but also Nature prefereed colours?"

## What Nature asks

The research figure guide
([building and exporting figure panels](https://research-figure-guide.nature.com/figures/building-and-exporting-figure-panels/)):
avoid red/green combinations and rainbow scales; keys in the figure, not colour descriptions in the
caption; no coloured text; RGB, not CMYK. Its one **example colour-blind palette is Okabe-Ito**:
black `#000000`, orange `#E69F00`, sky blue `#56B4E9`, bluish green `#009E73`, yellow `#F0E442`,
blue `#0072B2`, vermillion `#D55E00`, reddish purple `#CC79A7`.

## What is in use (Fig. 1, ED 1-3)

- **The four mouse groups** (hue = genotype, lightness = age): WT 6W `#0072B2`, WT 12W `#56B4E9`,
  MYC+ 6W `#D55E00`, MYC+ 12W `#E69F00`. Fills, the hue at 28% over white: `#B8D8E9`, `#D0EAF9`,
  `#F3D2B8`, `#F8E4B8`. Box lines `#595959`. Contrasts take the same blue and orange.
- **Heatmaps:** diverging espresso `#4A3525` - white `#FAFAFA` - mint `#2A8A6D`; sequential
  `#F2F4F6` - `#8FA3B0` - `#243642`.
- **ED 2b-c**, mitochondrial class: `#8E2C74`, `#CC79A7`, grey `#D9D9D9`.
- **ED 2f-g**, the seven MitoPathway tiers: all of Okabe-Ito, Metabolism grey `#B8B8B8`.

So Fig. 1 uses four of Nature's eight: **black, bluish green, yellow and reddish purple are free.**

## Fig. 2's groups (from the Doc, "Nat metab version")

Constructs: EV, BCL-xL, PGC-1α, NRF1, Ndi1, LbNOX (matrix), LbNOX (cytosol); each with MYC off/on
(iMMEC Tre-MYC) or basal/raised (MYAZ). A panel holds EV and at most two constructs (Fig. 2k: EV,
LbNOX matrix, LbNOX cytosol).

## Recommended: Option A, the half of Nature's palette Fig. 1 does not use, one hue per mechanism

| construct | colour | fill (28% over white; superseded for bars on 2026-10-19 -- bars take the colour) |
|---|---|---|
| EV | light grey `#CCCCCC` | `#F1F1F1` |
| BCL-xL (anti-apoptotic) | black `#000000` | `#B8B8B8` |
| PGC-1α, NRF1 (mitochondrial biogenesis) | bluish green `#009E73` | `#B8E4D8` |
| Ndi1, LbNOX matrix (matrix NADH oxidation) | reddish purple `#CC79A7` | `#F1D9E6` |
| LbNOX cytosol (cytosolic NADH oxidation) | yellow `#F0E442` | `#FBF7CA` |

- **MYC state takes no colour**: it is on the axis (- / +), and colour that encodes nothing is not
  used (FIGURE_RULES.md).
- Constructs share a hue only where they share a mechanism, and never in the same panel; the axis
  names the construct.
- Within every panel the colours stay apart under deuteranopia, protanopia and tritanopia (Fig. 2k's
  EV / matrix / cytosol: closest pair dE >= 19, CIELAB).
- **EV is light grey, not mid-grey**: `#999999` has reddish purple's lightness and merges with it
  under deuteranopia (dE 7).
- **Yellow needs a dark outline** on dots (`grey20`, as on the page) and cannot carry a line on white.

## Rejected: Option B, Paul Tol's "muted", one hue per construct

EV `#CCCCCC`, BCL-xL indigo `#332288`, PGC-1α teal `#44AA99`, NRF1 green `#117733`, Ndi1 wine
`#882255`, LbNOX matrix purple `#AA4499`, LbNOX cytosol sand `#DDCC77`. Under deuteranopia and
protanopia the indigo and purple turn navy -- Fig. 1's WT blue -- and the teal turns grey, next to EV.

## Open

1. ~~**Figs 3-4** reuse EV, BCL-xL and PGC-1α and add **PGC-1α + BCL-xL** (and EV-EV): the palette
   should be set paper-wide; the combination needs its own colour (a darker green is one candidate).~~
   PGC-1α + BCL-xL done 2026-10-23 (`#00573F`, the decision at the top); EV-EV still open.
2. Bluish green sits near Fig. 1g's mint heatmap end, and reddish purple is ED 2b-c's mitochondrial
   class: other figures and other encodings, so acceptable, but noted.
3. ~~When the author chooses, add the construct palette to `_style.R` and `FIGURE_RULES.md`.~~
   Done 2026-10-19 (see the decision above).

Sources: [Nature research figure guide](https://research-figure-guide.nature.com/figures/building-and-exporting-figure-panels/);
[Nature blog, colour-blind readers](https://blogs.nature.com/blog/post_4-4/);
[Elsevier, in the eye of the beholder](https://www.elsevier.com/connect/in-the-eye-of-the-beholder).
