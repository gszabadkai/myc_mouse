# Human datasets: sources, versions and pins (both human repos)

2026-11-11. Every human dataset that enters `myc_human_validation` or
`myc_human_exploratory`, with where it came from, which release, how it was retrieved, the
checksum where one was pinned, the sample numbers, and the scripts and figures that use it.
For the Methods and the Data availability statement.

**How this was compiled.** Read-only, from:
- each repo's `data/*/README.md`;
- the input paths that each repo's scripts reference;
- the analysis-set notes in `myc_human_exploratory/docs/`;
- the `!Series_*` header lines of the local GEO series matrices.

The repos were read at:
- `myc_human_validation` `d3ac60e`, frozen;
- `myc_human_exploratory` `2c49de2` on `main`, plus the unmerged branch `bh3-mimetic-oxphos`
  for GDSC.

No file in either repo was changed. Checksums are as each README records them. They were not
recomputed for this note.

Long checksums are shortened here to their first and last characters. The full values are in
the README that each section names.

Citations: PMIDs and DOIs were checked in PubMed on 2026-11-11. For the GEO series, the papers
are the ones each series' own `!Series_pubmed_id` lines name.

---

## 1. What the paper uses

| Panel | Dataset | n drawn | Section |
|---|---|---|---|
| Fig. 4f, 4g; ED 6a–d | TCGA-BRCA RNA-seq (GDC STAR counts, via UCSC Xena) | 938 | 2.1 |
| Fig. 4f, 4g; ED 6a–d | SCAN-B, GEO **GSE202203** (raw counts) | 3,143 | 2.4 |
| Fig. 4f, 4g (TCGA covariates) | PAM50 (cBioPortal, PanCanAtlas 2018); ABSOLUTE purity; leukocyte fraction (PanImmune) | 938 | 2.2 |
| Fig. 4h | TCGA-BRCA GISTIC2, **ISAR-corrected** (PanCanAtlas); aneuploidy score and fraction genome altered (cBioPortal) | 1,043 | 2.3, 2.2 |
| Fig. 4i; ED 6a | Neoadjuvant: **GSE194040** (I-SPY2-990), **GSE164458** (BrighTNess), **GSE25066** (Hatzis) | 988 / 482 / 465 | 2.5 |
| Fig. 4j | **METABRIC** (author's copy; source Synapse `syn1757063`) | 1,981; strata 1,505 / 1,512 / 435 | 2.6 |
| Scores in all of the above | Human MitoCarta 3.0; MSigDB 2026.1.Hs Hallmark (proliferation); MYC compendium GMX (`MYC_UP.V1_UP`, the reference MYC score); Felsher integrative signature; CollecTRI human | – | 4 |

pCR events in Fig. 4i: 319 / 236 / 93. Breast-cancer deaths in Fig. 4j: ER-positive (IHC)
1,505 / 431; ER-positive (`ER.Expr`) 1,512 / 432; ER-negative (IHC) 435 / 183; all patients
1,981 / 623.

The other datasets below (DepMap, PRISM, GDSC, MC3, RPPA, the Firehose GISTIC run, the TCGA
CDR and the Menegollo forkscale objects) feed validation gates or exploratory analyses that
are not drawn in the Letter.

---

## 2. Patient cohorts

### 2.1 TCGA-BRCA expression

| | |
|---|---|
| What | GDC harmonised **STAR – Counts** (GENCODE v36), `unstranded` counts |
| Retrieved from | UCSC Xena GDC hub: `https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-BRCA.star_counts.tsv.gz` |
| Retrieved | 2026-08-28 |
| File | `myc_human_validation/data/raw/tcga_expression_xena/TCGA-BRCA.star_counts.tsv.gz`, 138,503,699 bytes |
| SHA-256 | `058d79121460c535b73312247a55d3108d18ea6ddd0ccc1070b5dd93ceeedaa4` |
| Shape | 60,660 genes x 1,226 samples; values `log2(count + 1)` |
| Record | `myc_human_validation/data/tcga_expression_xena/README.md` |
| Read by | validation script 01. Its cached matrix, `data/raw/tcga_brca_expression/tcga_brca_star_counts.rds`, is also read directly by exploratory E37. |

**Why Xena and not the GDC API.**
- Xena re-hosts the same GDC STAR-Counts quantification.
- The GDC per-file route was abandoned: 1,231 files, 5.2 GB, about 7 h.
- 42 of those GDC files remain on disk, under `data/raw/tcga_brca_expression/` with
  `gdc_manifest.txt`.
- Four of them match their inverted Xena columns exactly, so Xena equals GDC `unstranded`
  under `log2(x+1)`.
- One of them also supplies the GENCODE v36 gene annotation, because Xena's probeMap is
  GENCODE v22.

**Sample handling, by script 01:**
- Counts are recovered as `round(2^x - 1)`.
- **Five `-01A` columns are averages of replicate aliquots.** Xena averaged them; they carry
  half-integer counts. They are dropped in favour of the patient's clean `-01B`/`-01C`.
- Primary tumour (`-01`) only.

**The cohort at each stage:**

| Stage | n |
|---|---|
| Primary tumours | **1,095** |
| Complete on PAM50, purity and leukocyte fraction: the Fig. 4 and ED 6 set | **938** |
| Solid-tissue normals, restored by exploratory E37 for E36 (not drawn in the paper) | **113** |

**Cite** the TCGA Network 2012 and, for the harmonised data, the GDC.
- Cancer Genome Atlas Network, *Nature* 2012;490:61. PMID 23000897,
  doi:10.1038/nature11412

### 2.2 TCGA-BRCA clinical and genomic annotation

All files are open access with no token.

**From the GDC API, committed in `myc_human_validation/data/tcga_clinical/`:**

| File | Query | n | SHA-256 |
|---|---|---|---|
| `gdc_brca_cases_2026-08-28.tsv` | `api.gdc.cancer.gov/cases`, `project_id = TCGA-BRCA` | 1,098 cases | `15aa4689…c3d9d9` |
| `gdc_brca_rnaseq_aliquots_2026-08-28.tsv` | `api.gdc.cancer.gov/files`, STAR – Counts aliquots; gives the sequencing plate | 1,231 | `703ac03b…128ad` |

**From cBioPortal, 2026-08-28, committed as `tcga_brca_clinical_snapshot_2026-08-28.tsv`**
(SHA-256 `40fd7403…ec181`). The query was
`www.cbioportal.org/api/studies/<study>/clinical-data`, one call per attribute.

| Attributes | Study | Level |
|---|---|---|
| `PAM50_SUBTYPE` | `brca_tcga_pan_can_atlas_2018` | patient |
| `ANEUPLOIDY_SCORE`, `FRACTION_GENOME_ALTERED` | `brca_tcga_pan_can_atlas_2018` | sample |
| `ER_STATUS_BY_IHC`, `PR_STATUS_BY_IHC`, `IHC_HER2`, `HER2_FISH_STATUS` | `brca_tcga` (legacy study) | patient |

The legacy study is used for receptor status because the PanCanAtlas study carries no receptor
IHC.

**From the GDC PanCanAtlas publication pages**, `api.gdc.cancer.gov/data/<UUID>`, 2026-08-28.
The files are gitignored, so the README is the record:
`myc_human_validation/data/tcga_pancanatlas/README.md`.

| File | GDC UUID | SHA-256 | Used for | Cite |
|---|---|---|---|---|
| `TCGA_mastercalls.abs_tables_JSedit.fixed.txt` (ABSOLUTE) | `4f277128-f793-4354-a13d-30cc7fe9f6b5` | `f430a975…63cd` | **tumour purity, Fig. 4 covariate** | Carter 2012, *Nat Biotechnol* 30:413, PMID 22544022, doi:10.1038/nbt.2203 |
| `TCGA_all_leuk_estimate.masked.20170107.tsv` | `6f75c9d7-5134-4ed1-b8f3-72856c98a4e8` | `5a8268ca…82d9` | **leukocyte fraction, Fig. 4 covariate** | Thorsson 2018, *Immunity* 48:812, PMID 29628290, doi:10.1016/j.immuni.2018.03.023 |
| `PANCAN_ArmCallsAndAneuploidyScore_092817.txt` | `4c35f34f-b0f3-4891-8794-4840dd748aad` | `c8a7f3dd…e23c` | per-arm calls (G2 conditioning) | Taylor 2018, *Cancer Cell* 33:676, PMID 29622463, doi:10.1016/j.ccell.2018.03.007 |
| `mc3.v0.2.8.PUBLIC.maf.gz` | `1c8cfe5f-e52d-41ba-94da-f15ea1337efc` | `cb0caa0e…6b81` | TP53/PIK3CA/PTEN calls (validation covariates) | Ellrott 2018, *Cell Syst* 6:271, PMID 29596782, doi:10.1016/j.cels.2018.03.002 |
| `TCGA-RPPA-pancan-clean.txt` | `fcbb373e-28d4-4818-92f3-601ede3da5e1` | `06246573…00de` | validation Block E (protein); no PUMA on the panel | – |

PAM50 is the PanCanAtlas call: Hoadley 2018, *Cell* 173:291, PMID 29625048,
doi:10.1016/j.cell.2018.03.022.

**TCGA Clinical Data Resource** (`data/tcga_cdr/TCGA-CDR-SupplementalTableS1.xlsx`):
- GDC UUID `1b5f413e-a8d1-4d10-92eb-7c4ae739ed81`, SHA-256 `ea594c0f…5a5a`, snapshot
  2026-08-29.
- Liu 2018, *Cell* 173:400, PMID 29625055, doi:10.1016/j.cell.2018.02.052.
- **Snapshotted but read by no script in either repo.** The plan forbade TCGA survival
  analysis.

### 2.3 TCGA-BRCA copy number (GISTIC2 thresholded calls)

The files are gitignored. The record is `myc_human_validation/data/gistic_tcga_brca/README.md`,
retrieved 2026-08-28.

| Run | Source | SHA-256 | BRCA n | Role |
|---|---|---|---|---|
| **A: ISAR-corrected, PanCanAtlas** | `api.gdc.cancer.gov/data/a9dae2ab-9462-4f9a-9730-5bad520ab2d7` (`ISAR_GISTIC.all_thresholded.by_genes.txt.gz`; 24,203 genes x 9,991 pan-cancer samples) | gz `bbd82167…0331`; raw `dde3270d…02a3` | **1,043** | **primary; Fig. 4h** |
| B: Broad Firehose BRCA-TP, Level 4 | `gdac.broadinstitute.org/runs/analyses__2016_01_28/…/gdac.broadinstitute.org_BRCA-TP.CopyNumber_Gistic2.Level_4.2016012800.0.0.tar.gz` | `2ee43313…e597` | 1,080 | pre-specified sensitivity |

How the BRCA set was defined:
- membership by the GDC case list (section 2.2), not by either file;
- sample type `-01` only;
- one column per patient.

Run A's 1,043 patients are a strict subset of run B's 1,080.

**Cite** Taylor 2018 for ISAR GISTIC (as in section 2.2).

**Where Fig. 4h's numbers come from.** They are the saved G2 result,
`g2_cnv_cooccurrence.rds`:
- validation script 05, 2026-08-28;
- copied into `myc_human_exploratory/data/from_validation_g2/`, md5
  `33da579ae1777c6bd246f5e8518d967d`;
- drawn by exploratory E54.

### 2.4 SCAN-B: GEO GSE202203

| | |
|---|---|
| Accession | **GSE202203**, n = 3,207. **Not GSE96058**: that deposit is FPKM, already logged, and cannot carry mitoPPS. |
| Platforms | GPL11154 (HiSeq 2000; 2,913 samples) + GPL18573 (NextSeq 500; 294 samples) |
| Files | `GSE202203_RawCounts_gene_3207.tsv.gz`, plus the two series matrices |
| URLs | `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE202nnn/GSE202203/suppl/…` and `…/matrix/…` |
| Retrieved | 2026-08-31 |
| SHA-256 | counts `db38c17f…3c5709`; GPL11154 matrix `59613566…1350`; GPL18573 matrix `762b437c…2079` |
| As deposited | GRCh38; HISAT2 2.1.0; StringTie 1.3.3b `-e`; protein-coding only; symbols from UCSC knownGenes of 2014-09-22 |
| Record | `myc_human_validation/data/scanb/README.md` |
| Read by | validation scripts 16 and 17. The exploratory repo reads the copies in `data/from_validation/` (`scanb_vst.rds`, `scanb_linear.rds`, `scanb_pheno.rds`). |

**Analysis set: n = 3,143**, which is 3,207 minus the 64 tumours with no PAM50 call.

**Handling:**
- The counts are StringTie estimates and are **rounded** to integers.
- **19 of 89 OXPHOS-subunit symbols are 2014 names.** They are harmonised through MitoCarta's
  `Synonyms` column, which takes coverage from 0.775 to 0.989.
- The deposit's survival, treatment and ESR1/ESR2 fields are **dropped on load**
  (`SCANB_FORBIDDEN`).
- SCAN-B has no purity or leukocyte estimate, so its models carry subtype only.

**Cite** Dalal 2022, *Sci Rep*, PMID 35304506, doi:10.1038/s41598-022-08210-3 (the series'
paper). SCAN-B is ClinicalTrials.gov NCT02306096.

### 2.5 Neoadjuvant chemotherapy cohorts (pCR)

The cohorts were selected by a documented GEO search (447 series screened, 34 passing):
`myc_human_validation/docs/2026-08-28_D5_cohort_selection.md`.

The files are gitignored. The record is `myc_human_validation/data/neoadjuvant/README.md`.
All were downloaded on 2026-08-30 from `https://ftp.ncbi.nlm.nih.gov/geo/…`.

They are read by:
- validation script 12, which fetches and harmonises them, saving
  `results/neoadjuvant_cohorts.rds`;
- validation script 13, the outcome models, saving `results/h4_outcome_models.rds`;
- exploratory E39–E45 and E57, which read those two validation `results/` objects in place.

| Cohort | Accession | Platform; data as deposited | Files (SHA-256) | n in Fig. 4i (pCR) | Series papers (`!Series_pubmed_id`) |
|---|---|---|---|---|---|
| **I-SPY2-990** | **GSE194040** | Agilent arrays GPL20078 + GPL30493; gene-level log2, **ComBat-adjusted as deposited** | `GSE194040_ISPY2ResID_AgilentGeneExp_990_FrshFrzn_meanCol_geneLevel_n988.txt.gz` (`9d8d1d35…72c4`); two series matrices (`8a7e5d56…a702`, `20aa7bd2…cb19`) | **988 (319)** | Wolf 2022, *Cancer Cell*, PMID 35623341, doi:10.1016/j.ccell.2022.05.005; also PMID 40905691 |
| **BrighTNess** (TNBC, phase III) | **GSE164458** | RNA-seq GPL21290; deposited log2 | `GSE164458_BrighTNess_RNAseq_log2_Processed_ASTOR.txt.gz` (`7209cd0e…d9e8`); series matrix (`8e80fd0c…7d24`) | **482 (236)** | Loibl 2018, *Lancet Oncol* 19:497, PMID 29501363, doi:10.1016/S1470-2045(18)30111-6; Metzger-Filho 2021, *NPJ Breast Cancer* 7:142, PMID 34764307, doi:10.1038/s41523-021-00349-y |
| **Hatzis** (taxane–anthracycline; pooled, includes I-SPY1) | **GSE25066** | Affymetrix HG-U133A, GPL96; normalised table inside the series matrix | `GSE25066_series_matrix.txt.gz` (`beef319f…60fb`); `GPL96.annot.gz` (`88e0b223…51c4`) for the probe collapse | **465 (93)** | Hatzis 2011, *JAMA* 305:1873, PMID 21558518, doi:10.1001/jama.2011.593; also PMIDs 24337596, 36293478 |

**Notes:**
- **GSE25066's 465:** 508 → 470 complete cases → 465 after dropping the 5 HER2-positive
  patients. The same 465 carry 103 DRFS events (E39/E40, not drawn).
- **No cohort supports mitoPPS.** All three are deposited on a log scale. GSVA only.
- **BBC3 is absent from GSE25066.** Its only probe, `211692_s_at`, is multi-mapping and is
  dropped by the pre-fixed collapse rule.
- **GSE25066's four published predictors are forbidden as covariates.** They are derived from
  the same expression matrix.

### 2.6 METABRIC

| | |
|---|---|
| File | `myc_human_exploratory/data/menegollo_biclusters/METABRIC_DATA.RData`, 331,447,619 bytes, md5 `1d1069ef983b6609e5c9037d1f584a2d`, mtime 2017-07-21. **Untracked**; ignored by name in `.gitignore`. |
| Source | **Synapse `syn1757063`, collection `syn1688369`**, as recorded on 2026-10-03 on the author's authority, per the Methods of Menegollo et al. 2024 |
| Copied in | 2026-10-02, by the author |
| Contents | `METABRIC.data`: symbol-keyed, 23,043 genes, log2, per-gene median-centred. Clinical tables: `Complete_METABRIC_Clinical_Features_Data`, `Complete_METABRIC_Clinical_Survival_Data_from_METABRIC` (623 cause-specific events), `Clinical_Overall_Survival_Data_from_METABRIC` (888 all-cause events) |
| n | **1,981** |
| Fig. 4j strata | ER-positive (IHC) **1,505**; ER-positive (`ER.Expr`) **1,512**; ER-negative (IHC) **435**; endpoint breast-cancer-specific death |
| Companion file | `METABRIC_starting_data_final_corr_groups.Rdata`: 735,646,872 bytes, md5 `6b64b2b396ec4d6af6ad506501db5fa0`, dated 2024-02-26. It is the Menegollo compile output (`all.clinical.df`: sample IDs plus `MB{1,2,3}.forkscale`). It comes from the author's Drive folder `…/Cell Resource Paper 2020/Fig1_gene_expression/METABRIC/Data compilation/`. **A different file of the same name exists in `…/Data visualisation/`** (md5 `faba17e5…`); it is not this one. |
| Record | `myc_human_exploratory/data/menegollo_biclusters/README.md` (addendum of 2026-10-02/03) and `docs/2026-10-02_E39_respiratory_axis_decomposition_declaration.md` section 4.5 |
| Read by | exploratory E41 and E43 to E45 (the METABRIC ladder and the Fig. 4j values), and E57 (drawing) |

**Handling:**
- Array data, so **mitoPPS is unavailable**. GSVA on the matrix as supplied.
- OXPHOS subunits are recovered for 69 of 89 genes; the proliferation set for 293 of 318.
- **BBC3 is absent** from the matrix.

**Cite:**
- Curtis 2012, *Nature* 486:346, PMID 22522925, doi:10.1038/nature10983;
- Pereira 2016, *Nat Commun* 7:11479, PMID 27161491, doi:10.1038/ncomms11479.

### 2.7 Menegollo / Bentham bicluster snapshot (TCGA and METABRIC forkscale)

| | |
|---|---|
| Source | `github.com/gszabadkai/Menegollo_Bentham` (private), commit `8fdbb3437ae5537055d5d5429411bdb3b333c04a` (2024-06-18) |
| Snapshot | 2026-08-29, by `gh api` blob fetch; each file re-hashed and matched to its upstream blob SHA |
| Files | `TCGA_MB{1,2,3}_RNASeq_data_nonorm.RData` (1,037 TCGA primary tumours each); `TCGA.all.biclusters.RNAseq.Rdata` (849, cross-check only); `METABRIC_sort_data.RData`; `METABRIC_PC1_GSEA.RData` (n = 1,981) |
| Where | identical copies in `data/menegollo_biclusters/` of both repos (tracked) |
| Read by | validation script 15 (F3 / Block D); exploratory E32 and E33 |
| In the paper | not drawn |

**Cite** Menegollo 2024, *Cancer Res*, PMID 38924467, doi:10.1158/0008-5472.CAN-23-3172.

---

## 3. Cell-line resources (not drawn in the Letter)

### 3.1 DepMap 26Q1 and PRISM Repurposing 24Q2

The files are held in `myc_human_validation/data/raw/depmap/` (gitignored).
`myc_human_exploratory/data/raw/depmap` is a **symlink** to it.

The records are `myc_human_validation/data/depmap/README.md` and
`myc_human_exploratory/data/raw/depmap_README.md`.

| File | Release | Route | Bytes | SHA-256 |
|---|---|---|---|---|
| `Model.csv` | 26Q1 | browser, depmap.org/portal (behind Cloudflare; `curl` returns an HTML page) | 697,455 | `ea4e0b2a…83f8` |
| `OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv` (unstranded) | 26Q1 | browser | 305,007,605 | `0377be80…021a` |
| `CRISPRGeneEffect.csv` | 26Q1 | browser | 440,646,050 | `e610a4ce…89e` |
| `Repurposing_Public_24Q2_Extended_Primary_Data_Matrix.csv` | PRISM 24Q2 | figshare file 46630984 (DOI 10.25452/figshare.plus.25917643) | 72,456,953 | `3b6554cf…bc3d` |
| `Repurposing_Public_24Q2_Extended_Primary_Compound_List.csv` | PRISM 24Q2 | figshare file 46630981 | 719,567 | `7e78f590…a5a463` |

**Notes:**
- `CRISPRGeneEffect_24Q4_figshare_SUPERSEDED.csv` is kept on disk and must not be used. A byte
  guard stops any script that is handed it.
- Breast lines: 51 with both expression and CRISPR. Pan-cancer: 1,140.
- **Read by:** validation script 14 (Block G); exploratory E17 to E19.

### 3.2 GDSC release 8.5 and Cell Model Passports

The files are gitignored, in `myc_human_exploratory/data/raw/gdsc/`, downloaded 2026-09-28.

The README is on the **unmerged branch** `bh3-mimetic-oxphos`:
`git show origin/bh3-mimetic-oxphos:data/raw/gdsc_README.md`.

| File | Host | SHA-256 |
|---|---|---|
| `GDSC2_fitted_dose_response_27Oct23.xlsx` | `cmp.cog.sanger.ac.uk/download/` | `0e99120a…c899` |
| `GDSC1_fitted_dose_response_27Oct23.xlsx` | same | `6dfad242…b6a5` |
| `screened_compounds_rel_8.5.csv` | same | `d3cbb8b5…e782` |
| `model_list_20260921.csv` | `cog.sanger.ac.uk/cmp/download/` | `419f65d1…7006` |
| `cnv_summary_20260804.zip` / `.csv` | same | `17c816f5…7512277` / `9cefb776…cdd` |
| `mutations_summary_20260724.zip` / `.csv` | same | `f8f97a73…dce` / `998344a9…6ef` |

**Notes:**
- The old `cancerrxgene.org` bulk-download URL now returns HTTP 410.
- **Read by:** exploratory E31 (BH3-mimetic analysis), **on that branch only**.

---

## 4. Reference catalogues and gene sets

These are not datasets, but they decide what every score measures, and each is pinned.

| Resource | Version / pin | Where | Used for | Cite |
|---|---|---|---|---|
| **Human MitoCarta 3.0** | `Human.MitoCarta3.0.xls`, MD5 `3c0bd24e362238216e142bc708e41286`, 10,179,584 bytes; a byte copy of the Broad download (snapshot 2026-08-28) | `data/mitocarta_human/` in both repos | OXPHOS subunits (89 scored after the mtDNA strip), all mito arms, mitoPPS, symbol harmonisation | Rath 2021, *Nucleic Acids Res* 49:D1541, PMID 33174596, doi:10.1093/nar/gkaa1011 |
| **CollecTRI, human** | OmniPath web service, `organisms=9606`, 2026-08-28; SHA-256 (uncompressed) `d72660703f7ccb8068b75994a1e74986b451d0c098e36d62efd7a88e631d287d`; 64,723 interactions | `data/collectri_human/` in both repos | MYC activity estimator `M_b` (ULM); FOXO3 activity | Müller-Dott 2023, *Nucleic Acids Res*, PMID 37843125, doi:10.1093/nar/gkad841 |
| **MSigDB** via `msigdbr` | **2026.1.Hs**, read at runtime | – | Hallmark E2F + G2M → `PROLIF_STD` (327) and `PROLIF_DISJOINT` (318, without the Felsher genes): the **proliferation score**; `HALLMARK_MYC_TARGETS_V1`; Reactome comparators | – |
| **MYC signature compendium** | `myc_signature_genesets.gmx`, 16 signatures, from `mammary_geneset_library` tag `v1.0` (`cbd8f16d`), path `data/raw/from_myc_mouse/`; human-native | `myc_human_exploratory/data/genesets_myc_human/` | **`MYC_UP.V1_UP`, the reference MYC activity score in Fig. 4f and 4g**; the ED 6a MYC panel | – |
| **Felsher integrative MYC signature** | `felsher_integrative_signature.csv`, 67 genes (61 after the MitoCarta strip), library `v1.0`, blob `193a18f5…` | `data/genesets_from_library_human/` in both repos | MYC estimator `M_a` (GSVA) | – |
| **Curated metabolic list** (`GS_metabolic`, human sheet) | library `v1.0`, blob `21f79d9a…`; SHA-256 `891427de…0bb8` | `data/genesets_metabolic_human/` in both repos | second opinion only, never scored | – |
| **Cell-death sets** | (A) `cell_death_genes_consolidated.csv`, human columns, from `myc_mouse` commit `6a9c7dd5…`; (B) Tang 2024, 15 modalities, library `v1.0` | `myc_human_exploratory/data/genesets_celldeath_human/` | exploratory death axes | Tang 2024, doi:10.1016/j.csbj.2024.08.012 |
| **BH3 binding code** | hand-written table, 7 rows | `myc_human_exploratory/data/annotations/bh3_binding_code.csv` | Fig. 4g binding annotation (E49) | Chen 2005, PMID 15694340; Ku 2011, PMID 21060336; Wang 2023, PMID 37560128 |

---

## 5. Mouse material inside the human repos

Not human data. Listed so it is not mistaken for human data.

- **`myc_human_validation/data/from_myc_mouse/`**: result tables only, read beside human
  results and never as inputs.
  - Interaction genes, mitoPPS ANOVA, the transcribed OXPHOS–PUMA estimates.
  - Taken from `myc_mouse` `paper-final` `6a9c7dd` on 2026-08-27.
- **`myc_human_exploratory/data/from_myc_mouse/mitopps_scores.rds`**: this repo's mitoPPS
  object, for the rank comparison (E24 to E30).
  - md5 `b8a125af4bc0d5909ad03f6e126e2890`, copied 2026-09-10.

---

## 6. How the exploratory repo reaches the data

**Copies:**
- **`data/from_validation/`**: derived objects copied from validation `d3ac60e` on
  2026-08-31.
  - TCGA and SCAN-B VST and linear matrices, mito and MYC scores, covariates, SCAN-B
    phenotype and symbol map.
  - Gitignored, about 560 MB; regenerable by re-copying at that SHA.
- **`data/from_validation_g2/`**: the G2 object for Fig. 4h, tracked.

**Through a symlink:**
- **`data/raw/depmap`**, to the validation repo's `data/raw/depmap`.

**Direct reads of the validation repo, at run time.** These contradict the "nothing is read
across repos at runtime" line in `data/from_validation/README.md`.
- From `results/`:
  - `state_definition.rds` (E32 to E37);
  - `forkscale_models.rds` (E32);
  - `neoadjuvant_cohorts.rds` and `h4_outcome_models.rds` (E39, E40, E42, E43);
  - `g1_overlap_audit.rds` (E41).
- From `data/raw/`:
  - `tcga_brca_expression/tcga_brca_star_counts.rds` (E37).

**These files are in no commit of the validation repo,** because it ignores `results/` and
`data/raw/`. So the Fig. 4i inputs depend on that working tree.

---

## 7. Open points for Methods and Data availability

1. **The METABRIC source rests on an attribution, not a check.**
   - Nothing inside `METABRIC_DATA.RData` names Synapse.
   - The md5 is the only fingerprint.
   - The array platform and normalisation are not recorded in either repo.
   - Before submission, confirm `syn1757063` against the Menegollo Methods and the
     Synapse/METABRIC terms of access.
2. **Raw downloads live only on this machine.**
   - Validation `data/raw/` is about 2.8 GB; the METABRIC pair is 1.07 GB; GDSC is about
     129 MB.
   - Every public file can be re-downloaded and checked: URL plus SHA-256.
   - The two METABRIC files and the validation `results/` objects cannot. They need a backup.
3. **GDSC's README is only on an unmerged branch.** Merge or copy it if E31 stays in the
   record.
4. **Accessions for the Data availability statement:**
   - TCGA-BRCA: GDC project `TCGA-BRCA`, PanCanAtlas files by UUID;
   - SCAN-B: GSE202203;
   - I-SPY2: GSE194040;
   - BrighTNess: GSE164458;
   - Hatzis: GSE25066;
   - METABRIC: Synapse `syn1757063`;
   - DepMap 26Q1; PRISM Repurposing 24Q2 (DOI 10.25452/figshare.plus.25917643);
   - GDSC release 8.5.
5. **Processing to state in Methods, because it comes from the deposits:**
   - GSE194040 is ComBat-adjusted as deposited.
   - SCAN-B counts are non-integer StringTie estimates, rounded.
   - TCGA reached the analysis via Xena; five averaged aliquot columns were dropped.
   - BBC3 is absent from GSE25066 and from METABRIC.
