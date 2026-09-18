# Methods and limitations

## Three completed native computations

Two Tamarind Foldseek jobs completed on 16 September 2026. Forward input was experimental 1V0D chain A132–328 (197 residues); reciprocal input was the full 423-residue archived X5S1 predicted monomer. Exact inputs are in [queries](../analysis/foldseek/queries/). Settings/logs and complete returned tables are included for each direction.

Foldseek 10.941cd33; sensitivity 9.5; monomer search; alignment type 2 (3Di+AA); E cutoff 10; coverage cutoff 0; maxSeqs 1000; backtraces retained. The 1,000 value is a prefilter candidate limit, not a guarantee of 1,000 returned hits or exhaustive coverage. Both prefilters reached that limit. Logged database counts were 66,725,340 for `alphafold-uniprot50` and 324,204 for `pdb`; snapshot dates were undisclosed. Provider order determines reported rank; X5S1's E-value rank is 285, while provider rank is 290. Reported identity is retained as supplied by the engine.

The native FoldMason job `cad_foldmason_22_20260917` completed on 17 September 2026. It used 22 existing predicted structures, 14 species, including six archived Hypsibius models and 16 retrieved AFDB v6 models. The mouse anchor here is predicted full-length O54788, not experimental 1V0D. Native command: `easy-msa inputs out/result.fasta tmp --report-mode 1`; engine dd3c23500ce0893bc1b5f9b4a3e17c63796a5d0b; 64 threads; no refinement iterations; B-factor masking threshold 0. The AA alignment has 22 rows and 746 columns. A guide tree organizes alignment construction; it is not a statistically supported phylogeny. Native AA/3Di/guide files were byte-identical to the separately recorded local fallback, without making the versions, logs or runtimes identical. Native HTML remains in the original archive; this compact repository includes its alignment data, settings and log without bundling its report application.

## Local sequence control

MMseqs2 18-8cc5c searched eight full-length queries against 20,832 H. exemplaris UP000192578 proteins plus human/mouse DFFB (20,834 total) from UniProt 2026_03. Exact FASTAs, full 75-row output, preparation and run receipt are included. Core settings: `easy-search --format-mode 4 --alignment-mode 3 -s 7.5 -e 10 --max-seqs 1000 -v 1`. The actual run additionally specified four threads, retained backtraces and a documented 24-field output; these are explicit fresh operational choices. No iterative profile search was executed. A missing pair means not reported under this search, not non-homology or a measured E-value above 10.

## Coordinates and correspondence

Main molecular graphics retain the historical whole-domain superposition and residue numbering. They were rendered from real coordinates using custom/local Mol* 5.10.1. US-align was not rerun. Native specialist Mol* 5.11 fits are separate operations and do not replace the fixed mapping. The barcode uses the forward Foldseek alignment: 211 columns, 32 identical, 139 differing and 40 gap columns. Saved selected functional pairs use the historical mapping, checked against the separate reciprocal alignment. These different alignments must not be conflated.

Predicted local geometry is uncertain. C229/C238/H242/C307 map to C296/C311/H316/C392; D262/H263/H308 to D337/H338/H393; N299 to K384. These are selected noncontiguous positions, not contiguous motifs. The mapped nonconservative substitution and uncertain local confidence preclude claims of an intact active site or metal coordination. H313 mapping is sensitive to alignment method. Architecture is source-defined: X5S1's historical C3-like span is 204–403; no canonical CIDE-N/C1 is reported, while C2-like elements remain possible.

## Genomic and expression checks

Deposited genomic records establish BV898_03327/OQV22896.1 and an exact 423-aa CDS translation on MTYJ01000015.1 in GCA_002082055.1. The alternative boundary is an alternative **donor**. Individual historical bHd16413.1/.2 protein assignments across expression/genome reference builds remain unresolved.

Only the pooled 10,000-adult cohort is used: GSM2472501–2503 active, GSM2472504–2506 tun. Sum bHd16413.1/.2 TPM per library from deposited Kallisto tables. Saved active values are 1.565613, 0.997161, 1.597606; tun values 7.912150, 11.607570, 10.530630. Recorded means are 1.3867933333333333 and 10.016783333333333; ratio 7.222982035294852, displayed as 1.387, 10.017, 7.22×. No raw-read reanalysis, new significance test, native NGS run or transferred published FDR. Single-adult and 30-adult cohorts remain separate.

## Historical versus fresh

The X5S1 hypothesis was known before these guided reproductions. Recovered S1/S4 PDFs retain historical errors and are not redistributed here. Original Foldseek JSONs, standalone historical US-align alignment, historical proteome snapshot and original historical scripts remain missing. A current-database reproduction is not exact historical replay. The older locus JSON reflects what had been recovered at its original execution date; S1/S4 subsequently became available.
