---
date: 2026-11-07
tags: [project/myc_mouse, nat-metab, figures, comparison, read-only]
status: >
  Record of the old-vs-current comparison the 2026-11-07 handoff scheduled. Written by a
  read-only subagent against the Drive PDFs (local mount
  /Users/gs/G/data/MK_myc_2022/paper figures/NM_Figs_to_submit/) and the Doc; reproduced
  verbatim below. Nothing in either repo, on Drive or in the Doc was changed. Three items
  were re-checked by the main session against the files: the Fig. 3i Bcl-xL value
  (-0.0207 in outputs/natmetab/Fig3/Fig3i_prism_mitopps_contrast.csv, drawn positive),
  Fig. 1b's exponents (the PDF text reads "2.7 x 1017", "2.1 x 1020"), and the zero-CPM
  NTEV tumour in outputs/orthotopic/57_fig3g_ppargc1a_prism.csv. The rest is the
  subagent's reading, not re-verified.
---

# Current figures vs the myc_mouse repo (2026-11-07, read-only)

**Paths.** repo = `/Users/gs/code/myc_mouse` (scripts `figures/natmetab/`, outputs `outputs/natmetab/`); HX = `/Users/gs/code/myc_human_exploratory`; scratch = `/private/tmp/claude-501/-Users-gs-code-myc-mouse/36bc2107-e228-4702-9955-ed6994534be7/scratchpad/compare/` (PDF copies, renders, Doc text) and `.../scratchpad/compare_tools/` (my readers).
**How checked.** The Drive API delivered only NM_Fig1 (Figs 2-4 are 8.5-10.3 MB, over the connector's limit), so all ten PDFs were copied from the local Drive mount (`/Users/gs/G/data/MK_myc_2022/paper figures/NM_Figs_to_submit/`): sizes and times equal the Drive metadata, NM_Fig1's md5 equals the API download. Rendered with `qlmanage`. Repo panels located on each page by matching their `pdftotext -bbox` words at one common scale + offset (`compare_tools/place_match.py`, `axis_scales.py`); error bars and symbols read from vector geometry (HX `docs/legends/tools/pdf_error_bars.py`, `compare_tools/pdf_marks.py`); Prism values read from Drive's `NM_all_figs_corrected.prism` (11-06 17:36); Doc read fresh (tab "Nat metab version" = lines 1-197).
**Since the HX inventory (its §9 read the 18:01 PDFs):** Fig. 2, ED 4, ED 6 are byte-identical to those copies. Fig. 3 (23:16): text identical; all 43 changed paths lie in 3i's left graph (the F13 fix, §4). Fig. 4 (00:27): text identical and the renders match by eye except 4c's lung micrographs (re-embedded at ~1,050-1,240 ppi instead of 300-ppi JPEG; short scale bars now visible).

## 1. Current structure
**Fig. 1** (10-16 16:46 BST)
- a — MMTV-Myc schematic, IHC sections, % Ki67+ / % CC3+ boxes — BioRender + images + Prism `1a Ki67`, `1a CC3`
- b — nuclear MYC intensity, *Myc* mRNA per mouse, immunoblot (MYC, 5 OXPHOS subunits, VCL) — Prism `1b Myc WB Quant`, Prism `1b Myc mRNA` (repo values; **not** `Fig1b_myc_transcript.pdf`, §2), blot
- c — MYC-ER CC3 % and BrdU %, section images — Prism `1c MycER CC3/BrdU` + images
- d — mito module shares + per-MitoPathway MYC effect — `fig1d_mito_content.R` → `Fig1/Fig1d_mito_content.pdf`
- e — mitoPPS MYC effect, 11 rows — `fig1e_myc_reallocation.R` → `Fig1/Fig1e_myc_reallocation.pdf`
- f — 12W vs 6W expression ruler + window scatter — `fig1f_window_reallocation.R` → `Fig1/Fig1f_window_reallocation.pdf`
- g — BCL-2 family heatmap — `fig1g_bcl2_heatmap.R` → `Fig1/Fig1g_bcl2_heatmap.pdf`
- h — BCL-2 departure from MYC's fade — `fig1h_bcl2_fade.R` → `Fig1/Fig1h_bcl2_fade.pdf`
- i — *Bbc3:Bcl2l1* vs OXPHOS priority — `fig1i_coupling.R` → `Fig1/Fig1i_coupling.pdf`

**Fig. 2** (11-06 18:02; HX inventory §2/§9) — no repo panel
- a — immunoblot, iMMEC + Tre-Myc, EV/Bcl-xL, Dox 0-48 h — blot
- b — % apoptotic vs Dox time, EV vs Bcl-xL — Prism `2b`
- c — immunoblot, MYAZ EV/Pgc-1α — blot
- d — % apoptotic MYAZ, EV vs Pgc-1α — Prism `2d`
- e — % confluency over 64 h — Prism `2e`
- f — immunoblot, iMMEC EV/Pgc-1α ± Dox — blot
- g — % apoptotic, EV/Pgc-1α × Dox — Prism `2g`
- h — % apoptotic, EV/Ndi1 × Dox — Prism `2h`
- i — % apoptotic, MYAZ + Tre-Myc, EV/Ndi1/mt-/c-LbNOX × Dox — Prism `2i`
- j — apoptosis vs ABT-737, iMMEC vs MYAZ — Prism `2j` (P from a separate mixed-effects file)

**Fig. 3** (11-06 23:16)
- a — % apoptotic iMMEC EV/Bcl-xL, GFP/48 h/cont — Prism `3a`
- b — immunoblot of those cells — blot
- c — Myc relative to GFP — Prism `3c`
- d — immunoblot, iMMEC EV/Pgc-1α ± Dox — blot
- e — % apoptotic after adaptation — Prism `3e`
- f — *p19arf* and *Cdkn1a* qPCR (two graphs) — Prism `3f`
- g — tumour volume — Prism; immunoblot — blot; *Ppargc1a* CPM (log axis) — Prism `3g ... RNAseq counts`, data = `outputs/orthotopic/57_fig3g_ppargc1a_prism.csv` (script 57; NTEV vs NTPgc1a; points match)
- h — tumour volume — Prism; immunoblot — blot; PUMA protein — Prism; *Cdkn2a* CPM, 4 arms — Prism `3h ... Cdkn2a cpm`, data = `cdkn2a_prism_table.R` → `EDFig5/EDFig5_prism_cdkn2a_per_animal.csv` (means 218/236/235/183 match)
- i — left: OXPHOS mitoPPS vs matched control, bootstrap bars — Prism `3i ... mitoPPS` from `fig3i_prism_tables.R` → `Fig3/Fig3i_prism_mitopps_contrast.csv`; right: *Ppargc1a* CPM per tumour — Prism `3i ... retained` = `Fig3/Fig3i_prism_retained_dose.csv` (identical values)
- j — tumour volume NT/p19arf-KO × EV/Pgc-1α — Prism `3j`

**Fig. 4** (11-07 00:27)
- a — lung colonisation, EV+EV/Bcl-xL+EV/Bcl-xL+Pgc-1α — Prism `4a`
- b — orthotopic volume, EV/Ndi1 × EV/Bcl-xL — Prism `4b`
- c — spontaneous metastasis (% area) + representative lungs — Prism `4c` + images
- d — lung colonisation over time, EV vs Ndi1 — Prism `4d`
- e — relative burden, EV/mt-/c-LbNOX — Prism `4e`
- f — BCL-2 family standardised β heatmap, TCGA/SCAN-B — HX E50
- g — *BBC3:BCL2L1* vs OXPHOS priority — HX E52
- h — co-amplification odds ratios — HX E54
- i — pCR forest — HX E57
- j — survival forest — HX E57

**ED 1** (10-16 10:40)
- a MYC/MAX/MXD network LFC ± 95% CI — `edfig1a_myc_network.R` → `EDFig1/EDFig1a_myc_network.pdf`
- b MYC effect 12W vs 6W, slope 0.49 — `edfig1b_myc_effect_rescaled.R` → `EDFig1/EDFig1b_myc_effect_rescaled.pdf`
- c NES by programme — `edfig1c_nes_by_programme.R` → `EDFig1/EDFig1c_nes_by_programme.pdf` (placed at 96.8%, §2)
- d MYC-ER transgene schematic — BioRender

**ED 2** (10-16 17:01 BST)
- a design — `edfig2a_design.R` → `EDFig2/EDFig2a_design.pdf`
- b gene-set library — `edfig2b_geneset_library.R` → `EDFig2b_geneset_library.pdf`
- c pathway PCA + PC1 loadings — `edfig2c_pathway_axis.R` → `EDFig2c_pathway_axis.pdf`
- d MYC effect 12W vs 6W (mitoPPS), slope 0.64 — `edfig2e_myc_reallocation_ages.R` → `EDFig2e_myc_reallocation_ages.pdf`
- e window vs MYC, two scatters (r −0.02, 0.69) — `edfig2d_window_independent.R` → `EDFig2d_window_independent.pdf`
- f 143 MitoPathways ranked by the 6W MYC effect — `edfig2f_reallocation_ranked.R` → `EDFig2f_reallocation_ranked.pdf` (105%, §2)
- g ranked by the WT window — `edfig2g_reallocation_window_ranked.R` → `EDFig2g_reallocation_window_ranked.pdf` (105%, §2)

**ED 3** (10-16 10:42)
- a downstream genes' departure — `edfig3a_fade_downstream.R` → `EDFig3/EDFig3a_fade_downstream.pdf`
- b BCL-2 ratios, 6W × 0.43 vs 12W — `edfig3b_ratio_expected.R` → `EDFig3b_ratio_expected.pdf`
- c PUMA regulators' departure — `edfig3c_puma_regulators.R` → `EDFig3c_puma_regulators.pdf`

**ED 4** (11-06 18:02; HX §3/§9) — no repo panel
- a — *Cox4i1* qPCR, EV/Pgc-1α/Ndi1 — Prism `ED4a`
- b — *Ndufs3* qPCR — Prism `ED4b`
- c — OCR basal/maximal/ATP-linked, EV vs Pgc-1α — Prism `ED4c`
- d — *Bbc3* qPCR — Prism `ED4d`
- e — basal OCR, OCR ± rotenone, ATP/ADP (EV vs Ndi1) — Prism `ED4e` (three sheets)
- f — ΔΨm (TMRM) — Prism `ED4f`
- g — mitochondrial NAD(P)H/NAD(P)+ — Prism `ED4g`

**ED 5** (11-06 01:58; not in the HX inventory; no legend in the Doc)
- a spontaneous metastasis % area, EV+EV/Bcl-xL+EV/Bcl-xL+Pgc-1α — Prism `ED5a`
- b representative lungs — images
- c circulating tumour cells (fraction) — Prism `ED5c`
- d two TF-activity bar charts — **being redone (TF activity)**; current rasters from untracked `scripts/external/Pgc bclxl analysis.R`. A parallel session wrote untracked `scripts/58_orthotopic_tf_activity.R`, `figures/natmetab/proposals/p5d_tf_colours.R` and `outputs/natmetab/proposals/p5d/` at 22:23-22:25 today (not assessed)
- e immunoblots Twist1, MMP9, MMP2, β-actin — blot

**ED 6** (11-06 18:01; HX §6) — no repo panel
- a — correlations among OXPHOS, MYC and proliferation scores — HX E53
- b — *BCL2L1* and *BBC3* association with OXPHOS, and their difference — HX E53
- c — unadjusted Spearman, 15 genes + signature × 7 measures — HX E55
- d — 15 genes over the MYC × OXPHOS groups — HX E56

## 2. Fig. 1 and ED 1-3 against the repo
- **Present, current, at 100%:** Fig. 1d-i, ED 1a-b, ED 2a-e, ED 3a-c — every unrotated word of each current PDF is found at one offset (e.g. Fig. 1g 98/98, ED 1a 36/36, ED 2b 36/36); the only misses are rotated axis titles and plotmath pieces, which pdftotext splits differently. Every repo PDF is newer than its script, and both re-saved pages (Fig. 1 16:46, ED 2 17:01) post-date the last render (10:12).
- **Swaps done:** ED 2d = `EDFig2e_myc_reallocation_ages.pdf` (left), 2e = `EDFig2d_window_independent.pdf` (right); 2f = `EDFig2f_reallocation_ranked.pdf` (upper), 2g = `EDFig2g_reallocation_window_ranked.pdf` (lower). (`_final_pdfs/NM_Extended_Data_Fig2.pdf` still has both pairs the old way.)
- **Mismatch — Fig. 1b's *Myc* mRNA graph is not the repo panel.** `Fig1b_myc_transcript.pdf` (the Built table's 1b) is not on the page (5/22 words, only shared ticks). The placed graph is the Prism sheet `1b Myc mRNA`: same mice (group means 2,542/8,141/2,142/7,577 = repo; 24 symbols measured, within ~30 counts), but y axis 0-15,000 "(normalised counts)" not "×1,000", and P values typed by hand (exponents wrong, §4).
- **Mismatch — ED 2f and 2g placed at 105.1%** (x and y scale fitted from 11 words each): 179.3 × 47.0 mm instead of 170.7 × 44.8; 7-pt type prints at 7.4. The relinked 171-mm PDFs kept the old 180-mm frame (the same in `_final_pdfs`).
- **Mismatch — ED 1c placed at 96.8%** (24 words): 108.3 × 77.2 mm instead of 111.8 × 79.7; 6-pt type prints at 5.8.
- No retired, proposal or old-ED-5 PDF is placed anywhere: `_retired/`, `EDFig5/`, `Fig3/` and the proposals match only shared labels (≤12 words; the earlier `proposals/Fig1g_heatmap.pdf` 32 of 101 gene names), except two proposals of the same size and text as the panels they were promoted to (`proposals/EDFig2c_pathway_axis_wide.pdf`, `proposals/EDFig3c_puma_regulators.pdf`).
- **Repo records stale against the page:** FIGURE_RULES "Built" still lists Fig. 1b as the R panel and ED 2d/2e with the pre-swap scripts, and its ED 2 page line reads "d 114 + e 60 ... f 180 x 45; g 180 x 45"; `page_mock_247.R` draws the pre-swap ED 2. LEGENDS.md and handoff 2a.37 carry the swap.

## 3. Repo panels not in the current structure
| file (scripts in `figures/natmetab/`, outputs in `outputs/natmetab/`) | why unused |
|---|---|
| `fig1b_myc_transcript.R` → `Fig1/Fig1b_myc_transcript.pdf` | 1b's mRNA graph is now the Prism sheet `1b Myc mRNA` |
| `edfig5a_bbc3_routes.R` → `EDFig5/EDFig5a_bbc3_routes.pdf`, `EDFig5a_p53_specificity.pdf` | old ED 5 (Bbc3 across routes), superseded (author, handoff 2026-11-07) |
| `edfig5b_ox_rel_routes.R` → `EDFig5/EDFig5b_ox_rel_routes.pdf` | old ED 5b (ox_rel ruler), superseded |
| `fig3i_prism_tables.R` → `EDFig5/EDFig5b_prism_ox_rel_contrast.csv`, `EDFig5b_prism_ox_rel_per_animal.csv`, `EDFig5_prism_mito_share_per_animal.csv`, `Fig3/Fig3i_prism_mitopps_per_animal.csv` | tables for the superseded ED 5 / not drawn (3i uses the contrast table). Its other two CSVs feed 3i |
| `cdkn2a_prism_table.R` → `EDFig5/EDFig5_prism_cdkn2a_normcount_per_animal.csv`, `EDFig5_cdkn2a_per_animal_long.csv` | not drawn; its CPM CSV **is** used, by Fig. 3h, despite the EDFig5 name |
| `fig3i_route_comparison.R` → `Fig3/Fig3i_retained_dose.pdf`, `Fig3i_mitopps.pdf` | R record of 3i; the placed 3i is Prism |
| `Fig3/Fig3i_route_comparison.pdf` (10-28 16:43) | orphan from before the script split into two panels (it now writes only the two above) |
| `page_mock_247.R` → `Fig1_ED2_ED3_247.pdf` | a page check, never placed (pre-swap ED 2) |
| `_final_pdfs/` | first assemblies (10-16 10:40-10:42): ED 1 and ED 3 byte-identical to Drive; Fig. 1 and ED 2 superseded by the Drive saves |
| `_demo/` (6), `_retired/` (9), `proposals/` (incl. `NM_Fig3.pdf` + `NM_Fig2_3.prism` of 10-28, an older assembled Fig. 3) | not current by design |
| `outputs/orthotopic/57_fig3i_route_comparison.{csv,pdf}`, `57_eda_bbc3_dose.{csv,pdf}`, `57_per_animal_long.csv` | script 57 drafts/records; only `57_fig3g_ppargc1a_prism.csv` is used (3g) |

## 4. Open issues touching repo outputs
**F13 (Fig. 3i left, NM_Fig3.pdf 23:16).** Axis from its own ticks: 0.0 at 100.14 pt, 97.65 pt per unit. Symbol centres and bar ends:

| arm | drawn now | repo (`Fig3i_prism_mitopps_contrast.csv`) | verdict |
|---|---|---|---|
| NT + Pgc-1α | 0.053 [0.030, 0.115] | 0.053 [0.030, 0.116] | offsets fixed |
| Bcl-xL | **+0.020** [−0.010, +0.127] | **−0.021** [−0.050, +0.087] | **value's sign lost** |
| Pgc-1α + Bcl-xL | 0.569 [0.473, 0.622] | 0.570 [0.475, 0.622] | offsets fixed |
| p19ARF-KO + Pgc-1α | 0.071 [−0.050, 0.105] | 0.071 [−0.064, 0.105] | fixed; low end clipped at the axis foot (−0.05) |

- The exchanged ±error is fixed on all four arms. **New:** Bcl-xL is drawn 0.041 too high. Drive's Prism sheet `3i Orthotopic OXPHOS mitoPPS` reads `Bcl-xL,2,0.0207,0.0295,0.1076` (positive value; offsets still exchanged on every row), so the 23:16 PDF came from a later Prism state that is not on Drive, and the sign slip survived. HX F13's own table has the same slip (its "table's interval" for Bcl-xL, [−0.0088, +0.1283], is +0.0207 ± offsets; script 57's readout gives −0.021 [−0.050, +0.087]). Bcl-xL spans zero either way. HX's "ED 5b needs the same check" is moot: the current ED 5 has no ox_rel panel.
- **Fig. 1b:** the mRNA P values print "2.7 x 10<sup>17</sup>" and "2.1 x 10<sup>20</sup>" (text layer and render); the repo's are 2.7 × 10<sup>−17</sup> and 2.1 × 10<sup>−20</sup>. P = 0.76 matches.
- **Fig. 3g right:** one NTEV tumour has *Ppargc1a* = 0 CPM (`57_fig3g_ppargc1a_prism.csv`) and cannot sit on the log axis: 6 EV symbols drawn for "n = 7"; the median line (0.16) is the median of all 7.
- **Fig. 3h *Cdkn2a*:** Dunnett against NT + EV puts brackets across series (NT+EV vs Bcl-xL+EV, P = 0.46; vs Bcl-xL+Pgc-1α, P = 0.044), which `cdkn2a_prism_table.R` (header, lines 24-42) says are descriptive only, brackets to stay within a series, CPM read beside the normalised count. The Doc legend says "Dunnett's test against EV + EV". The text "*Cdkn2a* also fell in the Pgc-1α expressing tumours (Fig. 3h)" holds only for Pgc-1α + Bcl-xL (183 CPM); NT + Pgc-1α is 236 against 218 (P = 0.42).
- **Fig. 3i/3j against the repo's recorded spec** (FIGURE_RULES "cell-line constructs", `fig3i_route_comparison.R`): the p19ARF-KO arm is orange, CMYK 7/33/100/0, the exact value of MYC+ 12W on Fig. 1a, not `#AA4499`; the arms run NT-Pgc-1α, Bcl-xL, Pgc-1α+Bcl-xL, KO, not the recorded escalation; every symbol is filled (the R record opens those whose interval spans 0); the Doc's 3i legend lacks the "asymmetric" sentence `fig3i_prism_tables.R` requires.
- **Doc legends against LEGENDS.md / the scripts:** ED 1b drops "22 genes beyond ±3 not shown, included in the fit" (FIGURE_RULES makes that count mandatory); ED 3a drops "*Apaf1* (z = +2.14, P = 0.018) not shown" and the clamped *Htra2*/*Birc5* fills; ED 3b calls 0.43 "the slope of the overall change of all ratios" (the script fits the eight ratios other than *Bbc3:Bcl2l1*, `edfig3b_ratio_expected.R` lines 25, 47); Fig. 1i names only immune marker scores (the model also carries epithelial ones, `fig1i_coupling.R` lines 16, 106-108); Fig. 1f's null loses "from the same expression ventiles"; ED 2c loses the broken PC1 axis.
- **Known, still in the text (handoff 2a.25):** "OXPHOS ... priority score was significantly reduced between 6 and 12 weeks" and "fatty acid oxidation ... gained priority": neither window change passes BH (2a.25: OXPHOS adjusted P 0.46 WT, 0.17 MYC+; ED 2g: no WT change passes at all).

## 5. Doc citations against the assembled letters
- Scripted over every "Fig. Nx" / "ED Fig. Nx" in lines 1-197 (`compare_tools/cite_check.py`; 76 citations, ranges and lists expanded).
- **Only one points to a missing letter: "ED Fig. 5a-f" (line 108); ED 5 has a-e. Known; the author is correcting it.**
- All Doc figure-title links point to the current Drive files. ED 5 has a title link but no legend.
- Fig. 3d is the one panel never cited in the prose.
- Not asked, noticed (the letter exists but shows something else): "iMMECs carrying Myc and Pgc-1α (*see Fig. 2e*)" (2e is MYAZ confluency; iMMEC + Pgc-1α is 2f-g); "we expressed Ndi1 in the orthotopic model (*see Fig. 2d,e*)" (2d-e are Pgc-1α); "Pgc-1α ... raised ... ΔΨm (Fig. 2c, ED Fig. 4a-c)" (no Pgc-1α ΔΨm panel; ED 4f is Ndi1); "(Fig. 1f, ED Fig. 2f, g)" for the 6-to-12-week decline (2f is MYC's 6W ranking, not the window).
