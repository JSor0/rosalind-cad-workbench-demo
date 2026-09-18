# Provenance ledger for the review repository

This curated ledger derives from the integrated 17 September 2026 project ledger. It describes actual operation classes and links only included files. Original work remains preserved separately.

| Evidence | Actual origin | Included record |
|---|---|---|
| X5S1 hypothesis and fixed overlay | Jan's earlier investigation; historical US-align comparison, not a blind new discovery | [Coordinate extraction](COORDINATE_EXTRACTION.json), [exact queries](../analysis/foldseek/queries/) |
| Forward and reciprocal search | Two completed native Tamarind Foldseek jobs, 16 September 2026 | [Completion](../analysis/foldseek/COMPLETION.json), [result summary](../analysis/foldseek/RESULT_SUMMARY.json), full tables/settings/logs |
| Family alignment | Completed native Tamarind FoldMason job, 17 September 2026 | [Completion](../analysis/foldmason/provenance/cloud_completion.json), [native outputs](../analysis/foldmason/native/), [input manifest](../analysis/foldmason/provenance/INPUT_MANIFEST.json) |
| Sequence Viewer | Inspected the 22-row AA alignment, set mouse reference, queried H313 | [Rows](../analysis/foldmason/provenance/viewer_rows_summary.json), [reference](../analysis/foldmason/provenance/viewer_reference.json), [column query](../analysis/foldmason/provenance/viewer_H313.json) |
| Molecular Structure Viewer | Genuine interactive loading/styling; later specialist completed 12 fits and six native exports | Separate from the fixed custom-rendered main figures; [broad fit receipts](VIEWER_ALIGNMENT_RESULTS.json), [C3 fit receipts](VIEWER_C3_ALIGNMENT_RESULTS.json) and [pipeline](WORKBENCH_PIPELINE.md) |
| Public retrieval | PDB, AlphaFold, UniProt, InterPro, NCBI, GEO | Included input manifests, annotation responses and source URLs |
| Sequence comparison | Local MMseqs2 18-8cc5c | [Run receipt](../analysis/sequence/RUN_RECEIPT.json), exact FASTAs/full output |
| Expression | Local calculation from deposited processed Kallisto matrices | [Cohort](../analysis/expression/COHORT_MANIFEST.tsv), [saved validation](../analysis/expression/validation_and_provenance.json) |
| Genomic locus | Public deposited records plus local translation checks | [Audit](../analysis/locus/GENOME_RECORD_AUDIT_20260915.json) and its versioned inputs |
| Main figures/GIF | R composition and custom Mol* rendering of fixed real coordinates | [Code](../code/README.md), [animation receipt](ANIMATION_RENDER.json) |

## Immutable historical source

Archived viewer: [cad-c3-folddisco-overlay-molstar-v2.3.2.html](https://github.com/JSor0/cad-c3-molstar-viewer/blob/main/cad-c3-folddisco-overlay-molstar-v2.3.2.html). SHA256: `4bdd61fc74fe14b22de5632354190ec413601f6d299236e72531732b9e9eea31`. It is not bundled or modified. Access depends on the archive's own permissions. Its bundled application and unrelated dimer–DNA contexts are excluded.

Exact mouse query SHA256: `5debaa03c463e357fc520a03988e0db9e86adceb52687c62ceca0133f7471633`. Exact full archived X5S1 query SHA256: `0ab38c3bf658a063f45a28464ba2de7604694f9026ae624ba8b0a7d77cfc49d1`. No refit, renumbering or model replacement occurred in packaging.

## Portability and integrity

[SOURCE_FILES.tsv](SOURCE_FILES.tsv) records original workspace-relative origins, source hashes, packaged hashes and adaptations. Metadata has explicit account-namespace/local-path sanitization; raw scientific tables, FASTAs, PDBs, compressed matrices and final media remain byte-identical unless individually recorded. `archived-external-location:` denotes an omitted historical execution location, not a promised included file. Source-hash fields inside historical receipts refer to their original unsanitized bytes; use [VERSION_MANIFEST.tsv](VERSION_MANIFEST.tsv) to check packaged bytes.

Original historical files still missing: original Foldseek JSONs, standalone US-align alignment, exact historical proteome snapshot and original scripts. Fresh/current and historical results remain distinct. Failed operation payloads, superseded outputs, private session material and caches are excluded; omissions do not relabel failures as successes.
