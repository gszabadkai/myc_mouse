# Gene set library snapshot

A frozen copy of the `mammary_geneset_library` deliverable, used to feed the
myc_mouse fGSEA / GSVA / mitoPPS analyses. Consume as-is. **Do not rebuild these sets
in myc_mouse**, and do not edit them here - if the library changes, re-snapshot from a
new tag rather than patching files in place.

## Provenance

- Source repo: `mammary_geneset_library`
- Version: tag `v1.0`
- Commit (dereferenced): `cbd8f16d2b0f95c5d4e86bed6aa112e42538a34b`
- Snapshot date: 2026-07-02
- Snapshotted by: manual copy from the library `outputs/` and `results/` at tag `v1.0`

The tag is the point of truth: it pins this analysis to one reproducible library
build rather than "whatever was in the folder that day". To reproduce, check out
`v1.0` in the library repo and re-export.

## Contents

```
data/genesets_from_library/
  README.md                              # this file
  mammary_mito_myc_metab_v1_mouse.gmt    # combined master GMT (mouse) - omnibus/exploratory only
  provenance_table.csv                   # per-set provenance + METHOD TAG (fgsea / gsva / both)
  metabric_sets.rds                      # METABRIC bicluster fork-sets (from library results/)
  gray_chea_mito_tf_shortlist.csv        # Gray 2023 + ChEA mito-TF shortlist
  by_category/                           # the nine mouse category GMTs (the primary inputs)
    01_mitocarta_mouse.gmt
    02_myc_signatures_mouse.gmt
    03_mammary_development_mouse.gmt
    04_metabolism_mouse.gmt
    05_proliferation_mouse.gmt
    06_tf_targets_mouse.gmt
    07_biogenesis_discrimination_mouse.gmt
    08_apoptosis_mouse.gmt
    09_biogenesis_apoptosis_intersections_mouse.gmt
```

The human GMT tree from the library is **deliberately not snapshotted** - this is a
mouse dataset, and keeping human GMTs out of reach removes an easy wrong-species
mistake. If a human arm is ever needed (e.g. the Hannon analysis), snapshot it into a
separate, clearly named directory at that time.

## How to use it

- **Method tag governs routing.** Every set in `provenance_table.csv` carries a tag:
  `fgsea`, `gsva`, or `both`. Route each set to the method its tag names. The GSVA
  large-set size-filter exemption applies as recorded in the library decisions.
- **Category GMTs are the primary inputs, used per-category** - not pooled into one
  master NES. Running fGSEA per category is what keeps the multiple-testing scope
  legible and makes the "preferential, not just absolute" mito-focus argument (AP6)
  interpretable. The master GMT is for exploratory omnibus scans only (Block A).
- **Category -> action-point map** (see the finalisation plan section 5 for detail):
  - `01_mitocarta` -> AP6, AP-mtDNA, mitoPPS partition
  - `02_myc_signatures` -> AP4, AP-retention (Felsher)
  - `03_mammary_development` -> AP1, AP2, AP3, AP5 (`MG_*` via GSVA + projection)
  - `04_metabolism`, `05_proliferation` -> AP6 comparator panel
  - `06_tf_targets` -> AP6 decomposition (dev/TF lanes)
  - `07_biogenesis_discrimination` -> Figure 2 Category-7 partition
  - `08_apoptosis` -> AP-CD branch 2 + Gate 2 (with the Tang `data/cell-death/` sets)
  - `09_biogenesis_apoptosis_intersections` -> biogenesis-vs-apoptosis intersection claim
- The in-repo `gene_sets_list.rds` (89 sets incl. 50 MSigDB Hallmark) remains the
  Hallmark comparator source; the `04`/`05` category GMTs are the mammary-specific
  comparators.

## Not included / deferred in the library v1.0

Recorded so their absence is not mistaken for an error: the human Phase C arm
(Reed/Kumar/Pal-EMBO + Hannon), pass-2 TFs (PGR, STAT5), the Pal 2021 TEB
re-analysis, GRN confirmation, and deconvolution were deferred in the library and are
not in this snapshot. This is why the library keeps its `dev` branch.
