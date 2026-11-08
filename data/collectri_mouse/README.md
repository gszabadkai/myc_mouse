# CollecTRI mouse TF regulons, pinned

Transcription-factor regulons for the orthotopic TF-activity analysis behind Extended Data
Fig. 5d (`scripts/58_orthotopic_tf_activity.R`). Pinned on the model of
`../myc_human_validation/data/collectri_human/README.md`.

**Do not edit in place.** Re-snapshot and update this file if a newer retrieval is needed.

## Files

| file | what |
|---|---|
| `collectri_mouse_omnipath.tsv.gz` | the OmniPath response, gzipped with `gzip -n -9` and otherwise unmodified |
| `uniprot_accessions_2026-11-07.tsv` | UniProt's gene names for the 1,077 network nodes OmniPath names only by UniProt accession |

## Provenance

- Upstream: the OmniPath web service, https://omnipathdb.org. CollecTRI: Müller-Dott et al.,
  *Nucleic Acids Res* 2023. The mouse regulons are OmniPath's own translation of the human
  resource (`organisms=10090`).
- Query URL (verbatim):
  `https://omnipathdb.org/interactions?resources=CollecTRI&genesymbols=yes&fields=sources,references,curation_effort&format=tsv&organisms=10090`
- Retrieved 2026-11-07T22:12:10Z by `curl` (HTTP 200). A first retrieval the same day at
  01:07 (local) returned a byte-identical file.
- SHA-256 of the uncompressed TSV:
  `76a64b108adc9787a7a0501cfd1c3365382629397f5aa6e230a6a06bb59afee6`
  (12,127,342 bytes; a header and 39,849 interactions; no newline after the last row).
  `gzcat collectri_mouse_omnipath.tsv.gz | shasum -a 256` reproduces it. MD5
  `d138f1f2689d5dee6416ba302a2f847f`, which script 58 checks with base R's `tools::md5sum()`.
- Columns: `source`, `target` (UniProt), `source_genesymbol`, `target_genesymbol`,
  `is_directed`, `is_stimulation`, `is_inhibition` (written `True`/`False`),
  `consensus_direction`, `consensus_stimulation`, `consensus_inhibition`, `sources`,
  `references`, `curation_effort`.

## Why a snapshot and not a live decoupleR call

`decoupleR::get_collectri()` and `OmnipathR::collectri()` fail in this environment
(OmnipathR 4.0.0, Bioconductor 3.23) for every organism, human included:
"Join columns in `y` must be present in the data. Problem with `ncbi_tax_id`", raised while
OmnipathR builds its organism table from Ensembl's species page. Its fallback,
`static_table()`, then fails in `unnest_evidences()`. The human repo met the same failure on
2026-08-28 (`../myc_human_validation/data/collectri_human/README.md`). The web service answers
the query above directly. A dated file with a checksum is in any case reproducible, which a
live call is not.

## How the analysis reads it

`scripts/58_orthotopic_tf_activity.R` applies `decoupleR::get_collectri()`'s own
post-processing (decoupleR 2.17.0, `split_complexes = FALSE`), the call the original analysis
made: `True`/`False` to 1/0; `COMPLEX` sources mapped to AP1/NFKB (the mouse table has none);
one edge per source-target pair; mode of regulation +1 for stimulation and -1 otherwise.

## Accession-named nodes

OmniPath names 1,077 mouse nodes by UniProt accession because it found no gene symbol for them
(226 as a source, 1,003 as a target). An edge whose target carries such a name can never match
a gene in the expression table, so without resolution 5,870 of the 39,521 edges are dropped
silently, and 6,386 edges have an unreadable TF name (A0A979HLR9 for Sp1, A3KN90 for Nr3c2).

`uniprot_accessions_2026-11-07.tsv` holds UniProt's answer for each, retrieved 2026-11-07 from
`https://rest.uniprot.org/uniprotkb/stream` (query `accession:<id>`, fields `accession,
gene_primary, gene_names, organism_id, reviewed`). The batch search did not return 8 of them:
each is an inactive accession that UniProt merged into a Swiss-Prot entry, and its row carries
that entry's gene name, with the merge stated in `status`. SHA-256
`e7fe00834e0a44cf65796ee0e00c6067a686e63ba7a21e6ce18efe391d4d63e4`, MD5
`54e5e309d44ec6ddca41ca35b95210d2`. Script 58 resolves each name
by these rules, in order, and reports the counts:

1. `org.Mm.eg.db` (the annotation the expression table uses) gives one symbol: use it.
2. It gives several and UniProt's primary name is among them: use that.
3. It gives none and UniProt's primary name is a current `org.Mm.eg.db` symbol: use that.
4. Otherwise the accession stays, and cannot match a gene.

Checked on 2026-11-07: UniProt names 1,074 of the 1,077 and `org.Mm.eg.db` 857; where both
answer (854) they agree on 840. The 14 disagreements are pseudogene or paralog names on the
UniProt side (Atg4a-ps for Atg4a, Gapdhrt2 for Gapdh) or recent renamings (Ca12 for Car12),
which is why `org.Mm.eg.db` comes first.
