# Library reference docs (read-only)

Design-rationale documents copied from `mammary_geneset_library` **v1.0**
(tag `v1.0`, commit `cbd8f16d2b0f95c5d4e86bed6aa112e42538a34b`), for interpreting the
gene sets in `data/genesets_from_library/`.

**These are reference, not a build spec.** They explain *what each set means and why
it was built that way*. Gene sets are consumed from the snapshot in
`data/genesets_from_library/` and are never rebuilt from these documents. Do not run
or port the library's build scripts.

For the machine-readable per-set source of truth (provenance + `fgsea`/`gsva`/`both`
method tag), use `data/genesets_from_library/provenance_table.csv`. The docs here are
for the design logic the CSV cannot carry.

## What to consult, and when

**Set-design decisions (start here for any "why is this set like this" question):**
- `2026-06-15_phase_c_replan_and_library_decisions.md` — the primary decisions record:
  the category roster, the fGSEA vs GSVA split rationale, the rescoped Category 7
  algebra (oncogenic vs developmental biogenesis: MYC_SPECIFIC / CORE / DEVELOPMENTAL),
  the TF lanes, the naming grammar. The single most useful doc.
- `2026-06-15_phase_c_dataset_review.md` — which developmental atlases were reviewed
  and which were chosen as sources, with the reasoning.

**Category catalogs (what is actually in each set):**
- `mammary_dev_sets_catalog_w_decisions_on_consensus_and_method.md` — the developmental
  (`MG_*`) sets, with the consensus/union decisions and per-set method choices. Consult
  before putting a developmental set on a figure (AP1/AP2/AP3/AP5).
- `mammary_dev_sets_catalog.md` — the plainer catalog of the same sets.
- `Gray_et_al_developmental_TFS_selection.md` — how the developmental TF roster was
  selected from Gray 2023; the basis for the pass-1 TF lane (E2f1, Esr1, Esrra, Gabpa,
  Myc, Nrf1) and the deferred pass-2 (PGR, STAT5). Consult for AP6's TF decomposition.
- `dropped_small_sets.md` — sets removed at QC (size filter) and why; check here if a
  set you expected is missing.

**Provenance / literature (for the methods section and set attribution):**
- `gene set resources table with developmental stages and supplementary information.md`
  — the literature review and per-source provenance behind the sets.
- `gene_set_history.md` — provenance of the gene-set assets inherited from myc_mouse.
- `mammary_dev_signalling.md` — background on the signalling pathways of mammary
  development (estrogen / progesterone / prolactin axes). Context for interpreting the
  `MG_*` developmental sets and the ER-vs-Myc deconvolution, and useful for the
  intro/methods prose. This is biological *background*, not a set catalog — it explains
  the program the sets track, not the contents of any GMT.
- `README.md` (the library's own) — high-level overview of the library deliverable.

**Reference-only, do NOT act on:**
- `PLAN.md` — the library's build recipe (scripts 00-14). Included only so the
  construction logic is inspectable. It is **not** a build spec for myc_mouse; do not
  reconstruct sets from it. If the library ever needs changing, that happens in the
  library repo and produces a new tag, after which we re-snapshot here.

## Not copied

Deliberately left out of this reference folder: `PROMPT_DAY1.md` (a library-build
opening prompt, irrelevant here), `CHANGELOG.md` (marginal), the library's
`docs/reference/myc_mouse/` (read-only copies of our own scripts — circular), and the
human GMT tree (this is a mouse dataset; keeping human sets out of reach avoids a
wrong-species mistake).

## Deferred in library v1.0

Recorded so their absence is not mistaken for an error: the human Phase C arm
(Reed/Kumar/Pal-EMBO + Hannon), the pass-2 TF lane (PGR, STAT5 — hormone-signalling
developmental TFs, a scoped-but-unbuilt second batch), the Pal 2021 TEB re-analysis,
GRN confirmation, and deconvolution. Each is a known, already-scoped addition if a
later decision needs it; none is a bug.
