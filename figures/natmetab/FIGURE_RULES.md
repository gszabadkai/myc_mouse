# Figure rules: GS lab, Nature Metabolism panels

A running record of the rules for figure panels, started 2026-10-09 and added to as decisions are
made. `_style.R` implements them. Rules marked **[checked]** are enforced by `save_panel()`, which
refuses a panel that breaks one, and `_style_selftest.R` proves each check works. This file is the
seed for the GS-lab figure skill planned for the end of the figure sessions.

## Workflow

- **One R script per panel** in `figures/natmetab/`. Each draws one panel at its final printed size.
- **Assembly happens in Illustrator**, together with blots, images and Prism graphs. Place every
  panel at 100%, because scaling a panel changes the size of its text. Linked placement lets a
  re-rendered panel update in the figure.
- **Panel letters are added in Illustrator**: 8 pt, bold, lowercase (a, b, c). R panels carry no
  letters.
- **Output:** `outputs/natmetab/<Fig>/<Fig><letter>_<name>.pdf`, e.g.
  `outputs/natmetab/Fig1/Fig1b_myc_transcript.pdf`. Figures run `Fig1` to `Fig4` and `EDFig1` to
  `EDFig10`. Two R pieces of one panel share the letter and differ by name.
- **Content first:** each panel is discussed before it is built, including any sentence the data
  do not support.

## Page (Nature)

- Width 89 mm (single column) or 183 mm (double). Height at most 170 mm, which leaves room for the
  legend. **[checked]**
- Extended Data: up to 10 figures, one page each.
- Final submission: main figures as PDF; Extended Data as TIFF or JPEG, at most 300 dpi. PNG is not
  used. For the initial submission a PDF of the figures is enough, but rules that cost little are
  followed anyway.

## Type

- **Helvetica throughout:** R panels, Prism graphs and Illustrator labels. Nature accepts Helvetica
  or Arial. The manuscript's Arial does not need matching, because Nature typesets the text.
- **Sizes:** 6 pt for tick labels, keys, on-panel labels and P values; 7 pt for axis titles and
  strips. Nothing below 5 pt **[checked]**, superscripts included **[checked]**.
- **Black text only.** Colour belongs to marks and keys, never to text. **[checked]**
- Lettering is lower case with the first letter capitalised (Nature), e.g. "Myc mRNA (normalized
  counts)".
- Units go in parentheses. Thousands take commas (1,000). Zero is printed as "0".
- Gene symbols are italic (mouse: *Myc*); proteins are in capitals (MYC).
- **Export:** base `pdf()` with Helvetica, not embedded; Illustrator embeds the font on save.
  - Never `cairo_pdf`. It embeds the font but places every letter separately, and at 5-7 pt Apple's
    viewers collapse the spaces between words (commit a58fa1c; re-checked 2026-10-09).
  - Not "Arial through renamed Helvetica metrics" either: Apple's viewers then draw the text in
    Times.

## Axes

- **The y axis starts at 0 whenever an x axis is drawn.** **[checked]** So a magnitude never goes
  on a log axis. **[checked]**
- Signed quantities (log2 fold changes) run from negative to positive, with 0 inside the axis and
  marked by a line. **[checked:** the axis must reach 0**]**
- **Both ends of every continuous axis sit on a major or minor tick**, with no padding.
  **[checked]** `scale_y_nm()` and `scale_x_nm()` do this.
- An axis may end on a minor tick only where minor ticks are **drawn**. ggplot computes minor
  breaks it does not draw, so ending on one would look tickless. **[checked]**
  `scale_y_nm(minor = TRUE)` draws them halfway between labelled ticks, which trims the empty
  band that rounding up to a major tick leaves.

## Statistics on panels

- **P values come from the analysis of record** (DESeq2, or the model of record) and are passed in.
  A plot never computes them: no `ggpubr::stat_compare_means()`, no ggstatsplot. At n = 6 against 6,
  a Wilcoxon test cannot go below P = 0.0022, while DESeq2 gives 2.7 x 10^-17 for *Myc*.
- P values are exact, to two significant figures, with an italic capital *P* (the journal's house
  style). Below 0.001 they are written as x 10^k, with the exponent at full size.
- DESeq2's adjusted P values are IHW-adjusted (weighted BH). Legends say "IHW-adjusted P, DESeq2
  Wald test".
- Brackets are drawn with `ggpubr::geom_bracket()` through `brackets()` in the style file. That is
  the only use of ggpubr.
- **Intervals** on effects are 95% (log2FC +/- 1.96 SE, Wald), at every age drawn. They are
  unadjusted, so where an interval excludes 0 but no adjusted P passes, the legend says so (ED 1a:
  *Mxi1*, *Mlx*, *Mlxip*).

## Showing data

- Every animal is shown (n = 6 per group) over a box. No bars with error bars, and no violins at
  n = 6.
- The chart form is checked against data-to-viz (https://www.data-to-viz.com/) and its caveats.
  Each script's header notes the caveats checked.
- **No single central value on a group with values in both directions.** A median or mean of a
  bimodal group sits inside one mode or between them and moves with the balance of membership, not
  with the shape. On ED 1c the Apoptosis median went +0.74 -> -0.94 and the means of three more rows
  changed sign, while the shape barely moved (2026-10-10). Draw the points large enough to carry
  the result instead.
- **A display window that hides points** (ED 1b, +/-3 log2) is allowed only if every number is
  computed on all points and the legend gives the count left out.

## Colour

| group | label | points, lines | box or bar fill |
|---|---|---|---|
| WT, 6 weeks | WT / 6W | `#0072B2` | `#B8D8E9` |
| WT, 12 weeks | WT / 12W | `#56B4E9` | `#D0EAF9` |
| MYC+, 6 weeks | MYC+ / 6W | `#D55E00` | `#F3D2B8` |
| MYC+, 12 weeks | MYC+ / 12W | `#E69F00` | `#F8E4B8` |

- **Encoding:** hue is genotype and lightness is age (Okabe-Ito).
- **Fills** are the hue at 28% over white, as solid colours, so Prism uses identical values.
- **Neutrals:** box outlines and medians are `#595959`; axes are black.
- **Genotype alone:** WT `#0072B2`, MYC+ `#D55E00`.
- **Genotype effects** take the MYC+ hue of their age: "MYC+ vs WT, 6W" `#D55E00`, "MYC+ vs WT,
  12W" `#E69F00` (`contrast_cols` in the style file).
- **Diverging fills:** `#4A3525` (espresso) to `#FAFAFA` (white, pinned at zero) to `#2A8A6D`
  (mint). Not blue-red, which would read as genotype.
- **Unsigned levels:** `#F2F4F6` to `#8FA3B0` to `#243642`.
- **Consistency:** the same group has the same colour in every panel, R or Prism. Colour that
  encodes nothing is not used.

## Labels

- **Groups:** two lines, "WT / 6W" and "MYC+ / 12W". MYC+ matches the manuscript text.
- **Contrasts (effects):** what is compared, then where. Genotype effects: "MYC+ vs WT, 6W" and
  "MYC+ vs WT, 12W" (in use from ED 1a-c). Change across the window, still a draft: "12W vs 6W,
  MYC+".
- Gene symbols on axes are italic; set and programme names are not. Row labels start with a
  capital (Nature lettering).

## Built (2026-10-10)

| slot | script | content |
|---|---|---|
| Fig. 1b | `fig1b_myc_transcript.R` | *Myc* mRNA per animal, beside the MYC western |
| ED 1a | `edfig1a_myc_network.R` | MYC/MAX/MXD network: MYC+ vs WT at 6W and 12W, 95% intervals |
| ED 1b | `edfig1b_myc_effect_rescaled.R` | MYC effect 12W vs 6W, 2,648 MYC-responsive genes, slope 0.49 |
| ED 1c | `edfig1c_nes_by_programme.R` | NES of every gene set at 6W and 12W, by programme |

## Pending

- The development-contrast labels and colours (first needed at Fig. 1e/f).
