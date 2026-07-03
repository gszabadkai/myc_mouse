# Mammary developmental sets (`MG_*`) catalog and consensus working doc

**Purpose.** This file inventories every `MG_*` gene set emitted by the readers in
`R/06_build_mammary_dev_sets.R` (files 1-17; Chung Suppl2 deferred to script 07),
records what each set *biologically refers to*, and lays out the evidence for
deciding which `MG_<STATE>_CONSENSUS` sets are worth building.

**How to use it.** Read, then edit in place. The decisions that matter:
1. Which canonical cell states get a `_CONSENSUS` set at all.
2. For each, which sources/sets are allowed to vote (the `consensus_state_map`).
3. The combination rule: strict cross-source **intersection** vs **majority vote**
   (>= half of voting sources) vs union. (Empirical sizes for both are below.)
4. Whether to confront these groupings with the recent field reviews and adjust
   the state vocabulary itself.

Nothing here is committed to code beyond the per-source readers. The
`consensus_state_map` in the script currently holds only the Henry + Garcia-Sola
(Suppl1) mappings; the cross-source picture below is **provisional, for your
review**.

Status: 139 sets total, all mouse symbols (MGI). `method` is `fgsea` (size-capped
DE marker lists, truncated to top-N by logFC, padj < 0.05) or `gsva` (published
membership lists kept whole, size-exempt).

---

## 1. Canonical cell-state vocabulary (current)

| State | Meaning |
|---|---|
| `MASC` | Mammary stem cell (note: fetal MaSC vs adult bipotent stem differ markedly) |
| `BASAL` | Mature basal / myoepithelial |
| `BASAL_PROG` | Basal progenitor |
| `LP` | Luminal progenitor (field-ambiguous: see notes) |
| `HR` | Luminal hormone-sensing, mature (ER+/PR+) |
| `ALV` | Luminal alveolar / secretory |

**Tensions to resolve against the literature** (these directly affect what a
consensus means):
- **MASC**: fetal mammary stem (Giraddi fMaSC, embryonic) is transcriptionally
  distinct from adult bipotent/basal stem (Henry CD1d+, Saeki Stem pseudotime).
  Pooling them may be wrong.
- **LP**: "luminal progenitor" is used inconsistently. Henry CD61+ luminal =
  classic LP; Saeki `LAPRO` = luminal *alveolar* progenitor; Garcia-Sola
  `Luminal_Progenitor` = global progenitor; Gray `AP` = alveolar progenitor.
  These are not the same population.
- **HR vs ALV progenitor split**: several sources separate progenitor from mature
  within both HR and ALV (Saeki LHPRO/LHOR, LAPRO/LALV; Garcia-Sola HormSensProg/
  Int/Dif). Decide whether consensus collapses progenitor+mature or keeps them apart.

---

## 2. Empirical cross-source overlap (PROVISIONAL voter mapping)

Computed by: union the same-state sets *within* each source, then intersect
*across* sources. This is the single most important table for deciding which
consensus sets are real. `int(all)` = strict intersection across all voting
sources; `maj(>=half)` = genes appearing in at least half the voting sources.

| State | Voting sources | `int(all)` | `maj(>=half)` | Verdict |
|---|---|---|---|---|
| **BASAL** | GRAY, PAL, SAEKI, GARCIASOLA, HENRY, GIRADDI (6) | **17** | 195 | Strong - strict intersection viable |
| **HR** | GRAY, SAEKI, GARCIASOLA, HENRY (4) | **68** | 264 | Strong |
| **ALV** | SAEKI, GARCIASOLA, HENRY, GRAY (4) | **37** | 247 | Strong |
| **LP** | PAL, SAEKI, GARCIASOLA, HENRY (4) | **0** | 43 | **Fails strict** - definitions diverge |
| **MASC** | HENRY, SAEKI, GARCIASOLA (3) | **0** | 31 | **Fails strict** - fetal vs adult stem |

**BASAL int(all) = 17 genes** (the canonical myoepithelial core):
`Acta2, Actg2, Aebp1, Cnn1, Cpne8, Krt14, Krt17, Krt5, Lhfp, Myh11, Myl9, Mylk,
Nexn, Palld, Scn7a, Tagln, Tpm2`

**ALV int(all) = 37 genes** (secretory/alveolar core):
`Aldh1a3, Aldoc, Atp6v1b1, Barx2, Btn1a1, Car2, Cck, Cd14, Cd81, Ceacam1, Cmas,
Crabp2, Csn3, Ctsh, Cytip, Dbi, Ehf, Elf5, Fcgbp, Kcnn4, Kit, Lipa, Lrg1, Ltf,
Lurap1l, Mfge8, Mgst1, Muc15, Ogfrl1, Phlda1, Plb1, Plet1, Plscr2, Tnfaip2, Trf,
Wfdc18, Wnt5b`

**HR int(all) = 68 genes** (not listed here for length; includes the expected
Esr1/Prlr/Foxa1/Cited1 hormone-sensing programme).

Pairwise intersections (to see which sources drive or break consensus):
- BASAL: SAEKI x GARCIASOLA 188, GRAY x SAEKI 170, GRAY x PAL 168; GIRADDI is the
  weakest partner (27-59). Dropping GIRADDI would lift int(all) above 17.
- HR: all pairs 108-138 - uniformly strong.
- ALV: GARCIASOLA x HENRY 150, SAEKI x HENRY 129; GRAY (AP, 120 genes) is weakest
  (46-66) - GRAY `AP` is alveolar *progenitor*, not mature alveolar.
- LP: PAL x SAEKI 1, PAL x HENRY 0, GARCIASOLA x HENRY 3 - essentially disjoint;
  only SAEKI x GARCIASOLA 38 overlaps. The PAL LP set is tiny (15).
- MASC: HENRY x SAEKI 1, HENRY x GARCIASOLA 4, SAEKI x GARCIASOLA 26 - disjoint.

**Implication for the combination rule.** Strict intersection is defensible for
BASAL/HR/ALV. For LP and MASC it returns nothing - they need either a majority
rule, a narrower (within-source-type) definition, or no consensus set. This is
where the field reviews should adjudicate.

---

## 3. Full per-source catalog

Columns: set | method | stage | size | biological referent.

### GRAY et al. 2023 (15 sets) - FACS major types + DE contrasts
| Set | method | stage | size | Refers to |
|---|---|---|---|---|
| `MG_AP_GRAY` | gsva | adult | 120 | Alveolar progenitor major type |
| `MG_HS_GRAY` | gsva | adult | 556 | Hormone-sensing major type |
| `MG_BASAL_GRAY` | gsva | adult | 603 | Basal major type |
| `MG_HEVSLE_AP_GRAY_UP/_DN` | fgsea/gsva | adult | 150/150 | High- vs low-estrogen, AP lineage |
| `MG_HEVSLE_HS_GRAY_UP/_DN` | fgsea/gsva | adult | 150/150 | High- vs low-estrogen, HS lineage |
| `MG_HEVSLE_BA_GRAY_UP/_DN` | fgsea/gsva | adult | 150/150 | High- vs low-estrogen, basal lineage |
| `MG_TEB_VS_DUCTAL_AP_GRAY_UP/_DN` | fgsea/gsva | puberty | 150/150 | TEB vs ductal, AP lineage |
| `MG_TEB_VS_DUCTAL_HS_GRAY_UP/_DN` | fgsea/gsva | puberty | 150/150 | TEB vs ductal, HS lineage |
| `MG_TEB_VS_DUCTAL_BA_GRAY_UP/_DN` | fgsea/gsva | puberty | 150/150 | TEB vs ductal, basal lineage |

Consensus candidates: `MG_BASAL_GRAY`->BASAL, `MG_HS_GRAY`->HR, `MG_AP_GRAY`->ALV
(progenitor; weakest ALV partner - decide). HEVSLE/TEB contrasts are directional
DE, not state signatures - exclude from consensus.

### PAL et al. 2017 (15 sets) - timeline + FACS clusters
| Set | method | stage | size | Refers to |
|---|---|---|---|---|
| `MG_BASAL_PAL2017` | fgsea/gsva | adult | 200 | Basal (top-200 markers) |
| `MG_LUMINAL_PAL2017` | fgsea/gsva | adult | 200 | Luminal (top-200) |
| `MG_LP_PAL2017_CIII` | fgsea/gsva | pubertal | 15 | Luminal progenitor (CD55+) |
| `MG_LUMINAL_PAL2017_CI/CII` | fgsea/gsva | pubertal/adult | 31/23 | Luminal1 / Luminal2 |
| `MG_BASAL_PAL2017_CIV/CVI/CVII` | fgsea/gsva | mixed | 26/43/30 | Basal subclusters |
| `MG_MIXED_PAL2017_CV` | fgsea/gsva | pubertal | 57 | Mixed-lineage |
| `MG_CLUSTER_I..VI_PAL2017` | fgsea/gsva | embryonic..adult | 8-41 | Developmental timeline clusters |
Consensus candidates: `MG_BASAL_PAL2017`->BASAL, `MG_LP_PAL2017_CIII`->LP (tiny,
breaks LP). Timeline/mixed clusters: descriptive, exclude.

### SAEKI et al. 2021 (10 sets) - integrated markers + scGSVA pseudotime
| Set | method | stage | size | Refers to |
|---|---|---|---|---|
| `MG_BASALPRO_SAEKI` | fgsea/gsva | integrated | 230 | MaSC / basal progenitor |
| `MG_BASAL_SAEKI` | fgsea/gsva | integrated | 303 | Basal |
| `MG_LAPRO_SAEKI` | fgsea/gsva | integrated | 232 | Luminal alveolar progenitor |
| `MG_LALV_SAEKI` | fgsea/gsva | integrated | 212 | Luminal alveolar (mature) |
| `MG_LHPRO_SAEKI` | fgsea/gsva | integrated | 210 | Luminal hormone-sensing progenitor |
| `MG_LHOR_SAEKI` | fgsea/gsva | integrated | 215 | Luminal hormone-sensing (mature) |
| `MG_STEM_SAEKI_GSVA` | gsva | integrated | 160 | Stem (S5 pseudotime top-160) |
| `MG_BASAL_SAEKI_GSVA` | gsva | integrated | 240 | Basal (S4 pseudotime) |
| `MG_ALV_SAEKI_GSVA` | gsva | integrated | 500 | Alveolar (S2 pseudotime) |
| `MG_HOR_SAEKI_GSVA` | gsva | integrated | 200 | Hormone-sensing (S1 pseudotime) |
Consensus candidates: BASAL (`MG_BASAL_SAEKI`), HR (`MG_LHOR_SAEKI`), ALV
(`MG_LALV_SAEKI`), LP (`MG_LAPRO_SAEKI`?), MASC (`MG_STEM_SAEKI_GSVA`,
`MG_BASALPRO_SAEKI`). Note Saeki uniquely supplies progenitor/mature splits and
two parallel definitions (marker vs pseudotime) per state.

### HENRY et al. 2021 (36 sets) - FACS populations + mEC/mNEC clusters
FACS lineage signatures (kept whole):
**as for method use both fgsea/gsva for ALL**
| Set | size | Refers to | -> state |
|---|---|---|---|
| `MG_MASC_CD1DPOS_HENRY` | 63 | CD1d+ MaSC | MASC |
| `MG_BASAL_CD61POS_HENRY` | 46 | Basal CD61+ | BASAL_PROG |
| `MG_BASAL_CD61NEG_HENRY` | 63 | Basal CD61- | BASAL |
| `MG_LUMINAL_CD61POS_HENRY` | 63 | Luminal CD61+ | LP |
| `MG_LUMINAL_CD133POS_HENRY` | 40 | Luminal CD133+ | HR |
| `MG_LUMINAL_CD133NEG_HENRY` | 48 | Luminal CD133- | ALV |
mEC clusters (epithelial; states per user's UMAP):
- BASAL: `MG_MEC_C1_HENRY` (200), `MG_MEC_C8_HENRY` (199)
- ALV: `MG_MEC_C2_HENRY` (112), `C4` (116), `C5` (144), `C10` (199)
- HR: `MG_MEC_C3_HENRY` (200), `C6` (181), `C7` (70), `C9` (195)
- unassigned/descriptive: `MG_MEC_C11..C15_HENRY` (32-200)
mNEC clusters (`MG_MNEC_C1..C15_HENRY`, 54-200): non-epithelial / contamination
controls - NO consensus, descriptive only.

### GIRADDI et al. 2018 (12 sets) - NMF signatures (method: only gsva, large)
| Set | stage | size | Refers to |
|---|---|---|---|
| `MG_FMASC_GIRADDI` | embryonic | 857 | Fetal MaSC NMF |
| `MG_FMASC_5A/5B_GIRADDI` | embryonic | 1253/683 | fMaSC subclusters |
| `MG_MMPR_GIRADDI` | postnatal | 234 | Multipotent mammary progenitor |
| `MG_MMPR_6A/6B_GIRADDI` | postnatal | 1022/1049 | MMPr subclusters |
| `MG_BASAL_GIRADDI` | adult | 187 | Adult basal NMF |
| `MG_LUMINAL_GIRADDI` | adult | 303 | Adult luminal NMF |
| `MG_BASAL_BALANCER_GIRADDI` | adult | 937 | Basal balancer |
| `MG_LUMINAL_BALANCER_GIRADDI` | adult | 477 | Luminal balancer |
| `MG_MATRIX_GIRADDI` | developmental | 259 | Matrix/stromal NMF |
| `MG_IMMUNE_GIRADDI` | developmental | 282 | Immune NMF |
Consensus candidates: `MG_BASAL_GIRADDI`->BASAL (weakest BASAL partner). fMaSC ->
MASC *only if* fetal stem is pooled with adult (contested). Balancer/Matrix/Immune
/MMPr: descriptive, exclude.

### SCHEELE et al. 2017 (9 sets) - StemID pubertal clusters (EXCLUDED from consensus)
`MG_CLUSTER_1..9_SCHEELE` (26-182). Inferred identities carried in `category` for
readability only (C1=HR, C2/C3=ALV, C6/C9=BASAL, C4=stromal, C5=luminal,
C7=proliferating, C8=mixed). TEB/puberty-specific context - user decision: do NOT
let them vote against adult/integrated consensus.
**method: only gsva**

### GARCIA-SOLA et al. 2021 (31 sets) - global types + trajectories
**method: both fgsea/gsva for ALL**
Named cell-type markers (Suppl1): 
| Set | size | -> state |
|---|---|---|
| `MG_BASAL_BASAL_GARCIASOLA` | 200 | BASAL |
| `MG_BASAL_MYOEPITHELIAL_GARCIASOLA` | 200 | BASAL |
| `MG_LUMINAL_PROGENITOR_GARCIASOLA` | 200 | LP |
| `MG_LUMINAL_ALVPROG_GARCIASOLA` | 200 | ALV |
| `MG_LUMINAL_ALVSEC_GARCIASOLA` | 200 | ALV |
| `MG_LUMINAL_HORMSENSPROG/INT/DIF_GARCIASOLA` | 200 each | HR |
| `MG_STEM_GARCIASOLA` | 200 | (MASC? top markers stromal Col1a1/Dcn/Fabp4 - EXCLUDED) |
| `MG_INMUNO_GARCIASOLA` | 200 | immune - EXCLUDED |
Trajectory subclusters (descriptive, EXCLUDED from consensus - would double-vote):
- `MG_HS_C0..C9_GARCIASOLA` (186-200): HS pseudotime trajectory
- `MG_ALV_C0..C10_GARCIASOLA` (200): alveolar pseudotime trajectory

### CHUNG et al. 2019 (11 sets) - ATAC open/closed (EXCLUDED from consensus)
**method: only GSVA for ALL**
`MG_{FETAL,BASAL,LP,ML}_{OPEN,CLOSED}_CHUNG_ATAC`,
`MG_CLUSTER6_OPEN_CHUNG_ATAC`, `MG_FETAL_{BA,LP}_LIKE_OPEN_CHUNG_ATAC` (19-300).
Chromatin accessibility, not expression - a different modality. Descriptive only;
would dilute the scRNA-marker intersections. Decide separately whether the `_open`
sets (Basal_open, LP_open, ML_open) are worth a *separate* ATAC-consensus.

---

## 4. Sets excluded from consensus (summary)

- All directional DE contrasts (`*_UP`/`*_DN`, HEVSLE, TEB_VS_DUCTAL).
- All Scheele clusters (puberty-specific).
- Garcia-Sola HS/Alv trajectory subclusters (intra-lineage pseudotime; double-vote).
- Garcia-Sola Stem (stromal markers) + Inmuno (immune).
- All Chung ATAC sets (different modality).
- Henry mNEC (non-epithelial) and mEC C11-C15 (unassigned).
- Giraddi Balancer / Matrix / Immune / MMPr subclusters; fMaSC pending MASC decision.
- Pal timeline clusters I-VI and Mixed.

---

## 5. Open decisions for you

1. **Build which consensus sets?** Evidence supports BASAL, HR, ALV cleanly.
   LP and MASC fail strict intersection - drop, or redefine, or use majority vote?
	1. *Keep only BASAL, HR and ALV*
2. **Combination rule:** strict intersection (robust, small) vs majority `>=half`
   (larger, looser) vs per-source-type. Pick per state or globally.
	1. *use `maj(>=half)` defined gene set for all 3 consensus sets*
3. **Progenitor vs mature:** collapse (HR = HORMSENS{Prog,Int,Dif} pooled) or keep
   separate consensus sets (`MG_HR_PROG_CONSENSUS` vs `MG_HR_MATURE_CONSENSUS`)?
	1. *do not generate consensus for progenitor states, keep them individual*
4. **fetal MaSC:** include Giraddi fMaSC in a MASC consensus, or keep fetal stem
   as its own (single-source) descriptive set?
	1. *NA* - *no MASC consensus is generated*
5. **Minimum size / minimum #sources** for emitting a consensus set (plan said
   size >= 15; how many sources must agree?).
	1. *NA - all sets are >15*
6. **ATAC consensus:** any value in a separate Chung-only open-chromatin consensus?
	1. No
7. Confront 1-4 with the recent field reviews and adjust the state vocabulary if
   the canonical states themselves need revising.
	1. Voacbulary to use (only for the consensus, do not change the individual set names): **BASAL -> BMYO (basal myoepithelial); ALV -> LASP (luminal adaptive secretory precursor); HR -> LHS (luminal hormone sensitive)** as per [Gray-et-al-2025-review](https://doi.org/10.1016/j.devcel.2025.06.032)
