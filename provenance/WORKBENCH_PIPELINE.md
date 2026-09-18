# Actual Workbench pipeline

1. **Inputs:** experimental mouse CAD C3 (1V0D A132–328); existing archived predicted X5S1 monomer. No new predictions.
2. **Rosalind-connected Tamarind Foldseek:** whole-domain search against `alphafold-uniprot50`; reverse full-X5S1 search against provider `pdb`. Both native jobs completed.
3. **Rosalind-connected Tamarind FoldMason:** 22 existing predicted structures, 14 species. Native AA/3Di alignments, guide tree, report, settings and log completed. The compact handoff carries data outputs, settings and log; the report application remains in the source archive.
4. **Molecular Structure Viewer:** genuine loading, styling, inspection and later specialist Mol* 5.11 fits/exports. These differ from the historical fixed US-align comparison used for the main custom-rendered illustrations.
5. **Sequence & Alignment Viewer:** successful 22-row inspection, mouse reference selection and H313 column query. It did not generate FoldMason or phylogenetic results.
6. **Life Sciences Databases:** later successful AlphaFold/UniProt/NCBI skill-helper retrievals. Earlier failed attempts and ordinary public requests retain separate attribution. The low-identity annotation audit used ordinary public requests.
7. **Supporting local computation:** MMseqs2, genomic CDS translation/identifier checks and aggregation of deposited Kallisto TPMs. R/Python composed the figures; custom Mol* produced the main molecular renders/GIF.

**Not completed or claimed:** native NGS execution, raw RNA-seq reprocessing, new AlphaFold predictions, new US-align run, experimental validation. The specialist tasks did not invoke the app-only Rosalind launcher or switch to a GPT-Rosalind model. Selecting a plugin or reading its instructions alone is not counted as execution.

This attribution follows the saved operation records and integrated ledger. The repository is for review of the existing analysis and provenance, not a new run.
