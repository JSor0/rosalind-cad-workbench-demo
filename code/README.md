# Used code and reproducibility boundaries

The four Python analysis scripts were used in the completed project; their original hashes match saved execution/audit receipts. This review copy adapts only input/output paths and the explicitly described local-binary check. See [source hashes and changes](../provenance/SOURCE_FILES.tsv). Python uses the standard library. The native Foldseek/FoldMason jobs are represented by their exact inputs, returned data, settings and logs; these scripts do not resubmit jobs.

## Offline inspection

From the repository root:

```sh
python3 code/summarize_mmseqs_direct_20260916.py --validate-only
```

This validates the saved table/FASTA content and produces no search or changed results. Original binary identity is recorded, not reverified; executable bytes are not included. Preparation/receipt hash linkage refers to the original recorded source, while the repository manifest checks sanitized copies.

## Reproduce supporting calculations when needed

These commands are documented for a reviewer and were not used to replace the canonical results during handoff. Use new, absent output directories inside the repository:

```sh
python3 code/summarize_native_foldseek_20260916.py --output rebuild/foldseek
python3 code/reanalyze_geo_10k_expression.py --output rebuild/expression
python3 code/audit_genome_records_20260915.py --output rebuild/locus
```

They validate/summarize preserved inputs; no network or new structural/sequence search is initiated. JSON reports will have current script hashes, output locations and dates. Historical values remain in `analysis/`. MMseqs2's original search command is preserved in the sequence receipt; re-execution would be a separate analysis, not required for review.

## Compose figures from frozen real render assets

Existing final media under `figures/` and `animation/` are authoritative. The four final R composition scripts are unchanged used sources. From the repository root, using installed R with `grid`, `png`, `jsonlite`, Helvetica and macOS Quartz:

```sh
mkdir -p rebuild/figures
Rscript code/figures/render_workbench_discovery.R code/figures rebuild/figures
Rscript code/figures/render_hero.R code/figures rebuild/figures
Rscript code/figures/render_expression_context.R code/figures rebuild/figures
Rscript code/figures/render_workflow.R rebuild/figures
```

No style, colors, arrow sizes, molecular geometry or scientific values were edited. Molecular PNG inputs and approved icons are included unchanged. These scripts compose existing real Mol* renders; they do not regenerate molecular coordinates or rerun alignments. Cross-platform graphics devices/fonts may alter text metrics; no claim of identical Linux rendering is made. The complete archived Mol* application and animation renderer are deliberately not bundled.

## Handoff verification

All four R compositions were rerun into an ignored scratch directory on macOS. Every output PNG was pixel-identical to the corresponding frozen final figure. The original final PNG/GIF bytes are retained. The portable MMseqs2 offline audit also passed all 75 saved rows; no new search was executed.
