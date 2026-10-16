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
- **Show options, don't describe them** (author, 2026-10-12: "difficult to make decisions just by
  reading them"). Proposals are real scripts drawn at their printed size, kept in
  `figures/natmetab/proposals/` (outside the dry-run set) and saved to
  `outputs/natmetab/proposals/`; a page mock-up places them with the figure's other panels on the
  page (`proposals/page_mock_247.R` since 2026-10-14; `page_mock.R` is the 170-mm one). A proposal
  moves into `figures/natmetab/` under its final letter once chosen.
- **Two switches in `_style.R` for mock-ups**, both off by default:
  - `options(natmetab.txt = 5, natmetab.txt_title = 6)`, set before a figure's panels are sourced,
    draws them one point smaller. The 5-pt floor does not move.
  - `options(natmetab.capture = TRUE)` makes `save_panel()` keep the checked, drawable panel in
    `.nm_captured`, so a page mock-up can take each panel exactly as its script saves it.

## Page (Nature)

- Width 89 mm (single column) or 183 mm (double). Height at most 170 mm, which leaves room for the
  legend. **[checked]**
- **For the initial submission the page is 183 x 247 mm** (Extended Data 180 x 247), the author's
  call on 2026-10-14. The 170-mm height applies at acceptance, so rows stay modular: a figure must
  be re-flowable onto 170 mm by moving rows, not by redrawing panels. `save_panel()` still checks
  each PANEL against 170 mm, which no single panel approaches.
- **Where font size is the limit, use the minimum, 5 pt** (author, 2026-10-14), and nothing below
  it: the floor stays 5 pt **[checked]**.
- Extended Data: up to 10 figures, one page each, **at most 180 mm wide** by 170 mm tall (Nature
  research figure guide, Extended Data formatting, checked 2026-10-12). **[checked:** an `EDFig`
  panel wider than 180 mm is refused**]** A row of panels must also leave room for the gaps
  between them.
- Final submission: main figures as PDF; Extended Data as TIFF or JPEG, at most 300 dpi. PNG is not
  used. For the initial submission a PDF of the figures is enough, but rules that cost little are
  followed anyway.

## Type

- **Helvetica throughout:** R panels, Prism graphs and Illustrator labels. Nature accepts Helvetica
  or Arial. The manuscript's Arial does not need matching, because Nature typesets the text.
- **Sizes:** 6 pt for tick labels, keys, on-panel labels and P values; 7 pt for axis titles and
  strips. Nothing below 5 pt **[checked]**, except a reduced exponent (next point). A panel in the
  Prism format is a point smaller (5 and 6 pt; see below), and so are Fig. 1d to 1i (author,
  2026-10-15: 5 pt for tick and row labels, keys, P values and tile numbers; 6 pt for strip and
  axis titles, group names and a colour bar's title).
- **A P value's exponent prints at the label's size, unless the author asks for it reduced**
  (Fig. 1d, 2026-10-15: "superscript size should be reduced"; Fig. 1i's *R*2 the same).
  `brackets(exponent = "reduced")`
  sets it in plotmath's script style, 70% of the label, which prints at 4 pt inside a 5-pt label:
  the one text allowed under 5 pt. A superscript in no named style is refused **[checked]**,
  because plotmath prints it at 70% with a gap after its minus.
- **Base `pdf()` rounds every font size to a whole point** (read from the PDF, 2026-10-15): 5.4 pt
  prints at 5, 5.5 at 6, and a 70% exponent in a 5-pt label at 4. Use whole-point sizes.
- **Black text only.** Colour belongs to marks and keys, never to text. **[checked]**
- Lettering is lower case with the first letter capitalised (Nature), e.g. "Myc mRNA (x1,000,
  normalised counts)".
- Units go in parentheses. Thousands take commas (1,000). Zero is printed as "0".
- Gene symbols are italic (mouse: *Myc*); proteins are in capitals (MYC).
- **Export:** base `pdf()` with Helvetica, not embedded; Illustrator embeds the font on save.
  - Base `pdf()` writes the page size in whole points, rounded down: a 56-mm page is 158 pt, which
    Illustrator shows as 55.74 mm (checked 2026-10-15). The panel is laid out on the full 56 mm, so
    what is cut is a sliver of empty margin, top and right.
  - Never `cairo_pdf`. It embeds the font but places every letter separately, and at 5-7 pt Apple's
    viewers collapse the spaces between words (commit a58fa1c; re-checked 2026-10-09).
  - Not "Arial through renamed Helvetica metrics" either: Apple's viewers then draw the text in
    Times.
- **Base `pdf()` draws the ASCII hyphen as a minus sign.** So labels avoid hyphens ("Folate and 1C
  metabolism", not "1-C"), and a hyphen is typed only where a minus is meant ("OXPHOS - proliferation"
  prints as OXPHOS − proliferation, the retired ED 2f).
- **For a true hyphen, write the Latin-1 soft hyphen `"­"`**: base `pdf()` prints it as a
  hyphen ("BCL­XL" prints BCL-XL; checked in Quick Look, 2026-10-12). Use it wherever a label
  needs one (protein names, "Pro-apoptotic").
- **Preview with Apple's engine** (Quick Look, `qlmanage -t`), not poppler's `pdftoppm`, which
  draws Helvetica-Oblique upright, so an italic *P* looks lost when it is not (2026-10-11).

## Axes

- **The y axis starts at 0 whenever an x axis is drawn** **[checked]**, except under a box plot (next
  point). So a magnitude never goes on a log axis. **[checked]**
- **A box plot's y axis need not start at 0** (author, 2026-10-15), so for resolution it runs from the
  tick closest below the lowest whisker to the first tick above the highest mark, a bracket's *P*
  label included. Whiskers run minimum to maximum, so the lowest and highest animals set the ends.
  `box_axis()` works the ends out in printed millimetres on a graph of fixed height, `label_mm()`
  measuring the label, and a point on the lower end is drawn whole (`coord_cartesian(clip = "off")`).
  **[checked:** a box plot is exempt from the from-0 check**]** First used on Fig. 1d's module shares.
- Signed quantities (log2 fold changes) run from negative to positive, with 0 inside the axis and
  marked by a line. **[checked:** the axis must reach 0**]**
- **Both ends of every continuous axis sit on a major or minor tick**, with no padding.
  **[checked]** `scale_y_nm()` and `scale_x_nm()` do this.
- An axis may end on a minor tick only where minor ticks are **drawn**. ggplot computes minor
  breaks it does not draw, so ending on one would look tickless. **[checked]**
  `scale_y_nm(minor = TRUE)` draws them halfway between labelled ticks, which trims the empty
  band that rounding up to a major tick leaves.
- A **schematic draws no axes**, so its scales need not end on ticks: the tick check applies to
  drawn axes only (ED 2a, whose canvas is in millimetres).
- Breaks are rounded to 10 decimals inside the scale helpers: an unrounded 0.4000000000000001
  fell outside a limit of 0.4 and was dropped (Fig. 1f, 2026-10-11; a self-test case covers it).
- **Pad the data range by 3-4% before the axis rounds to its end tick**, so a point at the
  extreme is not cut by the axis end (a box plot instead starts on the tick closest below its
  lowest whisker; see above). Choose the pad so that, on facets side by side, the last
  label of one facet cannot run into the first of the next (ED 2f, the old 2g, ends on a minor tick for this).

## Statistics on panels

- **P values come from the analysis of record** (DESeq2, or the model of record) and are passed in.
  A plot never computes them: no `ggpubr::stat_compare_means()`, no ggstatsplot. At n = 6 against 6,
  a Wilcoxon test cannot go below P = 0.0022, while DESeq2 gives 2.7 x 10^-17 for *Myc*.
- P values are exact, to two significant figures, with an italic capital *P* (the journal's house
  style). Below 0.001 they are written as x 10^k, with the exponent at full size unless the panel
  reduces it (Type, above).
- DESeq2's adjusted P values are IHW-adjusted (weighted BH). Legends say "IHW-adjusted P, DESeq2
  Wald test".
- Brackets are drawn with `ggpubr::geom_bracket()` through `brackets()` in the style file. That is
  the only use of ggpubr.
- **Bracket spacing is set in millimetres of the printed panel**, not as a fraction of the data,
  so a label clears the bracket above it whatever the scale (about 3.4 mm per level at 6 pt).
- **A *P* below 0.001 prints with an exponent and is about 16 mm wide at 6 pt.** Stagger it with
  its neighbour rather than shrinking the type, and let it overhang the y axis
  (`coord_cartesian(clip = "off")`) rather than be clipped (Fig. 1h proposal, *Bmf*, 2026-10-13).
- **Intervals** on effects are 95% (log2FC +/- 1.96 SE, Wald), at every age drawn. They are
  unadjusted, so where an interval excludes 0 but no adjusted P passes, the legend says so (ED 1a:
  *Mxi1*, *Mlx*, *Mlxip*).

## Showing data

- Every animal is shown (n = 6 per group) over a box. No bars with error bars, and no violins at
  n = 6.
- **Boxes: median and interquartile range; whiskers minimum to maximum** **[checked]** (author,
  2026-10-14). Prism draws every point only over min-max whiskers, so R matches it and one legend
  sentence covers R and Prism panels: "boxes show the median and interquartile range, whiskers the
  minimum and maximum, and dots individual animals". In R: `geom_boxplot(coef = Inf)`.
- The chart form is checked against data-to-viz (https://www.data-to-viz.com/) and its caveats.
  Each script's header notes the caveats checked.
- **No single central value on a group with values in both directions.** A median or mean of a
  bimodal group sits inside one mode or between them and moves with the balance of membership, not
  with the shape. On ED 1c the Apoptosis median went +0.74 -> -0.94 and the means of three more rows
  changed sign, while the shape barely moved (2026-10-10). Draw the points large enough to carry
  the result instead.
- **A display window that hides points** (ED 1b, +/-3 log2) is allowed only if every number is
  computed on all points and the legend gives the count left out.
- **A panel may be stacked parts in one PDF** (`stack_parts()`): each part is checked on its own,
  and `collect_guides = TRUE` gives one shared key (Fig. 1d). Parts side by side use `ncol = 2`
  (Fig. 1f, ED 2d). Patchwork lays out a shared key with its own theme (11 pt), so the helper gives
  it the panel theme and no outer margin. Stacked parts that share an x axis get the same
  `scale_x_nm()` call, so their zero lines align (the retired ED 2f).
- **Name points on a scatter only where the names fit off the cloud.** A `ggrepel` layer is
  given every point, the unnamed ones with an empty label, so names keep off all points, not just
  the named ones. Where names do not fit, the points are marked and the legend names them (ED 2e:
  eleven names at 62 mm crossed the cloud).
- **Significance glyphs mean one thing per figure** (author, 2026-10-11). Filled = *P* < 0.05 by
  the panel's test of record, open = not, and the key and legend name the test: Fig. 1e, the
  adjusted *P*; Fig. 1f, matched sets for expression and script 08's BH for priority (since
  2026-10-14); ED 3b, the 6W effect on each ratio (script 42, unadjusted), where only the 6W
  points are tested, the key draws the symbols in the 6W colour and the 12W points stay filled
  (2026-10-16). A panel without a test of record keeps its dots filled and draws no open ones.
- **Size marks the leading contrast** (author, 2026-10-11). Of two contrasts drawn per row, the one
  the text leads with is larger (1.6) and the other smaller (about 1.1), and the key shows both
  sizes: Fig. 1e, 6W over 12W; Fig. 1f, WT over MYC+ (WT is the gland's own programme).
  A panel that compares the two directly keeps them the same size (ED 1a, ED 1c: 6W against 12W).
- **A number the text quotes is asserted where it is computed** (the Task D principle), e.g. Fig.
  1d's "20-33%".

## The Prism format (an R panel among Prism graphs)

Fig. 1b sits in the bench row, beside the western and the Prism graphs, so it takes their format
rather than the R panels' (author, 2026-10-15: "Prism is less flexible, so will need to adapt to
that"). It is the first such panel; another R panel placed among Prism graphs would take the same.
- **The graph is 20 x 20 mm.** This is the data area, which Prism calls the graph size, fixed
  through `stack_parts()` with grid units. The y axis ends at the tick just above the top point
  (Fig. 1b: 12.5 thousand against a top animal of 12,180, a minor tick).
- **Large counts read in thousands**: tick labels 0, 5, 10 and "x1,000" in the title (author,
  2026-10-15), so the labels stay narrow beside a 20-mm axis.
- **No x axis and no group labels.** Illustrator sets the labels, as on the Prism graphs.
- **Axis title 6 pt; tick labels and P values 5 pt**: `theme_nm(txt = 5, txt_title = 6, ...)` and
  `brackets(txt = 5, ...)`.
- **Every line 0.9 pt**: axis, ticks, box outlines, whiskers, medians and brackets (author,
  2026-10-15: 1 pt read heavier than the Prism graphs). `lw_pt()` converts printed points to
  ggplot's linewidth: linewidth 1 prints 2.13 pt and `NM_LINE` prints 0.53 pt, read from the PDF.
  ggplot draws a box's median at twice its linewidth unless `median.linewidth` is given.
- **Brackets are straight lines, without the downward tips.** `brackets(tips = FALSE)` draws the
  bar alone, because ggpubr would keep the tips as zero-length segments, which are stray points in
  Illustrator.
- **Each bracket sits just above the animals it compares**, as Prism places them, not stacked
  above the whole graph (author, 2026-10-15). Fig. 1b: the 12W genotype bracket inside the graph
  over the 12W groups, the 6W one just above the axis end, and the MYC+ window above that.
  Brackets above the axis end are kept by `scales::oob_keep` on the y scale, drawn by
  `coord_cartesian(clip = "off")`, and held by a top margin. A *P* the author leaves off the
  panel goes in the legend (Fig. 1b: the WT window, *P* = 0.48).
- **The axis title is close to the axis, with its lines close together** (author, 2026-10-15):
  2.4 mm baseline to baseline (6 pt at 1.15), 0.6 mm from the tick labels. Plotmath's `atop()`
  spaces two lines about twice that far apart, so the title is set as a narrow column of its own,
  a part of the stack beside the graph, and centred on it as an axis title is.
- **The saved page is the title, the graph with its axis, and the brackets** (`natural_size()`), so
  there is no margin to trim. **[checked:** `save_panel()` refuses any fixed-size panel saved
  smaller than it needs, whose edges would be cut (ED 2d and 2e too)**]**
- At 0.9 pt, a box whose interquartile range is under about 1 mm on the axis prints as a nearly
  solid bar (Fig. 1b, both WT groups). Prism draws it the same way.

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
- **Changes across the window** (12W vs 6W) take the genotype colour: "12W vs 6W, WT" `#0072B2`,
  "12W vs 6W, MYC+" `#D55E00` (`geno_cols` with `dev_labels`; Fig. 1f, ED 2f).
- **Diverging fills:** `#4A3525` (espresso) to `#FAFAFA` (white, pinned at zero) to `#2A8A6D`
  (mint). Not blue-red, which would read as genotype.
- **Unsigned levels:** `#F2F4F6` to `#8FA3B0` to `#243642`.
- **Mitochondrial class of a gene set** (ED 2b-c): MitoCarta `#8E2C74`, curated (with mitochondria)
  `#CC79A7`, not mitochondrial `#D9D9D9` (`mito_class_cols`). One hue for the two mitochondrial
  classes (reddish purple, darker = defined on MitoCarta) against neutral grey, so mitochondrial
  against not reads first (author, 2026-10-14; the slate shades of 2026-10-12 did not separate
  them). The old red/orange/blue is out because blue and orange mean genotype; green is the
  diverging scale's end.
- **A schematic of the design teaches the contrast colours** (ED 2a): its arrows take the colour
  each contrast has on every effect panel; its text stays black.
- **Consistency:** the same group has the same colour in every panel, R or Prism. Colour that
  encodes nothing is not used.

## Labels

- **Groups:** two lines, "WT / 6W" and "MYC+ / 12W". MYC+ matches the manuscript text.
- **Contrasts (effects):** what is compared, then where. Genotype effects: "MYC+ vs WT, 6W" and
  "MYC+ vs WT, 12W" (from ED 1a-c). Change across the window: "12W vs 6W, WT" and "12W vs 6W,
  MYC+" (from Fig. 1f). Where the key names the contrast, the axis title names only the ruler and
  its unit ("Expression (log2 fold change)", "Priority (mitoPPS)").
- Gene symbols on axes are italic; set and programme names are not. Row labels start with a
  capital (Nature lettering).
- **"z-score", hyphenated, everywhere** (author, 2026-10-15), with the soft hyphen so it does not
  print as a minus (Type, above). An italic *z* for a test statistic stays as it is (Fig. 1h, ED 3a, 3c).

## Built

The panels of record, promoted 2026-10-14 (author: "promote them to final panels"). Each writes
`outputs/natmetab/<Fig>/<Fig><letter>_<name>.pdf` at the size it is placed. Sizes are for the
183 x 247 mm submission page (ED 180 x 247).

| slot | script | content; size (mm) |
|---|---|---|
| Fig. 1b | `fig1b_myc_transcript.R` | *Myc* mRNA per animal (x1,000), beside the MYC western, in the Prism format (graph 20 x 20 mm, no x axis, 0.9-pt lines, brackets without tips over the animals they compare, the WT window's *P* in the legend; since 2026-10-15); 30 x 28.5 |
| Fig. 1d | `fig1d_mito_content.R` | module shares (4 modules, every animal) over the MYC effect on 143 MitoPathways at 6W and 12W; one PDF, shared key; tick labels, key and *P* 5 pt, titles 6 pt, exponent reduced, y title "Fraction of nuclear-encoded transcriptome (%)"; upper graphs 28 mm tall, each axis from the tick below its lowest animal to the tick above its *P* label (since 2026-10-15); 100 x 60 |
| Fig. 1e | `fig1e_myc_reallocation.R` | mitoPPS MYC effect at 6W and 12W (12W smaller), 11 rule-chosen rows in three groups, filled = adjusted *P* < 0.05; row labels and key 5 pt, group names and axis title 6 pt; 78 x 55 (since 2026-10-15) |
| Fig. 1f | `fig1f_window_reallocation.R` | 12W vs 6W in WT and MYC+ (MYC+ smaller) on two rulers, with the old ED 2f's tests: filled = *P* < 0.05 (2,000 expression-matched sets for expression, script 54 and script 56; script 08's BH for priority, every point open); OXPHOS minus proliferation bracketed (WT, *P* < 0.0005, script 43); no median ticks; 5/6 pt (since 2026-10-15); 100 x 56 |
| Fig. 1g | `fig1g_bcl2_heatmap.R` | 16 BCL-2 family genes by role, group means of z, espresso-white-mint +/-2; 5/6 pt, the key titled "z-score", unframed, at the top beside the two short role names (since 2026-10-15); 66 x 56 |
| Fig. 1h | `fig1h_bcl2_fade.R` | the BCL-2 family (16 genes) by departure from MYC's fade (script 55), lines at *P* = 0.05 and 0.01, fill = 6W effect +/-0.5; the old ED 3b's top half, moved here on 2026-10-16 (author); 5/6 pt, the x title a strip of its own (two lines 2.4 mm apart), the bar unframed; 52 x 56 |
| Fig. 1i | `fig1i_coupling.R` | OXPHOS priority against *Bbc3:Bcl2l1*, adjusted for immune content; WT dashed and thinner; MYC+ *R*2 on the panel, *P* values in the legend; 5/6 pt, the square data area as large as the 56-mm height allows (since 2026-10-15); 56 x 56 |
| ED 1a | `edfig1a_myc_network.R` | MYC/MAX/MXD network: MYC+ vs WT at 6W and 12W, 95% intervals; 80 x 64 |
| ED 1b | `edfig1b_myc_effect_rescaled.R` | MYC effect 12W vs 6W, 2,648 MYC-responsive genes, slope 0.49; 58 x 58 |
| ED 1c | `edfig1c_nes_by_programme.R` | NES of every gene set at 6W and 12W, by programme; 112 x 80 |
| ED 2a | `edfig2a_design.R` | the 2 x 2 design, n = 6 per group, arrows in the contrast colours; 46 x 40 |
| ED 2b | `edfig2b_geneset_library.R` | the 986-set library by category, filled by mitochondrial class (purple against grey); 58 x 40 |
| ED 2c | `edfig2c_pathway_axis.R` | pathway-level PCA, PC1 broken for the one far animal, equal units, beside the ranked PC1 loadings with six sets named; 176 x 46 |
| ED 2d | `edfig2d_window_independent.R` | WT window against the MYC effect (r -0.02); MYC+ window against WT window (r 0.69); 44-mm data areas; 114 x 56 |
| ED 2e | `edfig2e_myc_reallocation_ages.R` | mitoPPS MYC effect 12W vs 6W, 143 pathways, fit slope 0.64, Fig. 1e rows in black; 44-mm data area; 60 x 56 |
| ED 2f | `edfig2f_window_two_rulers.R` | window, expression against priority, 143 pathways, WT and MYC+ facets (ED 2g until 2026-10-14); 84 x 64 |
| ED 3a | `edfig3a_fade_downstream.R` | the downstream genes by departure from MYC's fade, the same lines, fill = 6W effect +/-0.5, *Apaf1* left out; the old ED 3b's bottom half (2026-10-16; the old script, renamed); 86 x 44 |
| ED 3b | `edfig3b_ratio_expected.R` | each BCL-2 ratio's 6W effect x 0.43 (the common shrink of the other eight) against its 12W effect; rows by the 6W effect, filled = 6W effect *P* < 0.05 (script 42); the old Fig. 1h, moved and changed on 2026-10-16 (author); 86 x 45 |
| ED 3c | `edfig3c_puma_regulators.R` | the PUMA regulators by departure, the same lines, four beyond *P* = 0.01 labelled with their empirical *P*; 86 x 92 |

**Retired 2026-10-14** (`figures/natmetab/retired/`, outside the dry-run set; their old PDFs are
in `outputs/natmetab/_retired/`): the 2026-10-11 Fig. 1f (`fig1f_window_reallocation.R`), the
narrow ED 2c (`edfig2c_pathway_axis.R`), and the old ED 2f (`edfig2f_window_matched_controls.R`,
folded into Fig. 1f).

**Retired 2026-10-16** (author: "they do not add important details to the current text version"):
ED 3a, OXPHOS priority and the ratios per animal (`edfig3a_per_animal.R`), and ED 3d, the coupling
unadjusted and adjusted (`edfig3d_coupling_both_fits.R`). Both are in `retired/`, and their PDFs,
with the old Fig. 1h's and the old ED 3b's, in `outputs/natmetab/_retired/`.

**The pages** (`page_mock_247.R`, a check, not an Illustrator input; writes
`outputs/natmetab/Fig1_ED2_ED3_247.pdf`):
- **Fig. 1**, 237 of 247 mm: bench row a-c 56 (placeholder); d 100 x 60 + e 78 x 55 (60); f 100 x 56;
  g 66 + h 52 + i 56, all 56 tall, one row of 180 (author, 2026-10-15); 3-mm gaps.
- **ED 2**, about 215 mm: a 46 + b 58 (40); c 176 x 46; d 114 + e 60 (56); f 84 x 64.
- **ED 3**, 92 mm (2026-10-16): a 86 x 44 over b 86 x 45, beside c 86 x 92; 3-mm gaps.

**Proposals** (`figures/natmetab/proposals/`): all superseded, kept for the record --
`p1g_bcl2_balance.R`, `p1g_alt_ratio_rows.R`, `p1h_puma_members.R`, `p1h_alt_bcl2_plane.R`,
`p_fade_ranking.R`, and the sheets `page_mock.R` (170 mm) and `sheet_1g_compare.R`.

## Pending

- **Fig. 1 page budget: settled for submission by the 247-mm page** (2026-10-14). Fig. 1 holds the
  bench row a-c and d-i in three rows of two (`proposals/page_mock_247.R`; sizes in the proposals
  table are as placed there). At acceptance (170 mm) whole rows move out; the A / A2 / B layouts of
  2026-10-12 are the record of what that costs.
- **Text size: settled 2026-10-14** -- the minimum, 5 pt, where font size is the limit, and nothing
  below it (the author had asked on 2026-10-12 about going below 5 pt for the submission; not done,
  because Nature requires 5 pt at acceptance and the mock-ups show that area, not type, is the
  limit).
- **The second-round proposals await the author** (handoff 2a.31): the merged 1f, the purple class
  colours and ED 2c's labels. ED 3a and ED 3d were retired on 2026-10-16.
