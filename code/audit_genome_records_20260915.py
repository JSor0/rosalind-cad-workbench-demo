#!/usr/bin/env python3
"""Validate two deposited CDS products against a bounded genomic slice. Stdlib only."""
import argparse, datetime, difflib, hashlib, itertools, json, re
from pathlib import Path
import xml.etree.ElementTree as ET
ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / 'analysis/locus/inputs'
ap=argparse.ArgumentParser(description=__doc__)
ap.add_argument('--output',required=True,type=Path)
OUT=ap.parse_args().output.resolve()
OUT.mkdir(parents=True,exist_ok=False)
DERIVED=OUT/'derived'
DERIVED.mkdir(parents=True, exist_ok=True)
START, END, SCAFFOLD = 723181, 728004, 'MTYJ01000015.1'
files = ['OQV22896.1_OQV22897.1.ncbi-protein.xml', 'A0A1W0X5S1.uniprot.json', 'A0A1W0X5M8.uniprot.json', 'MTYJ01000015.1_723181-728004.ncbi-nucleotide.xml', 'GCA_002082055.1.ncbi-assembly-summary.json', 'MTYJ01000015.1.ncbi-nucleotide-summary.json']
def sha(b): return hashlib.sha256(b).hexdigest()
def put(path, content):
    data = content.encode() if isinstance(content,str) else content
    if path.exists():
        if path.read_bytes() == data: return
        raise FileExistsError('Refusing to replace different output: ' + str(path))
    path.write_bytes(data)
def quals(feature):
    return {q.findtext('GBQualifier_name'):q.findtext('GBQualifier_value') for q in feature.findall('GBFeature_quals/GBQualifier')}
def fasta(header, seq):
    return '>'+header+'\n'+'\n'.join(seq[i:i+70] for i in range(0,len(seq),70))+'\n'
nt_record = ET.parse(RAW / files[3]).getroot().find('GBSeq')
nt = nt_record.findtext('GBSeq_sequence').upper()
assert len(nt)==END-START+1
source=next(quals(f) for f in nt_record.findall('GBSeq_feature-table/GBFeature') if f.findtext('GBFeature_key')=='source')
gene=next(quals(f) for f in nt_record.findall('GBSeq_feature-table/GBFeature') if f.findtext('GBFeature_key')=='gene' and quals(f).get('locus_tag')=='BV898_03327')
codons=[''.join(c) for c in itertools.product('TCAG', repeat=3)]
aa='FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
table=dict(zip(codons,aa))
def translate(s):
    assert len(s)%3==0 and set(s)<=set('ACGT')
    return ''.join(table[s[i:i+3]] for i in range(0,len(s),3))
pairs={'OQV22896.1':'A0A1W0X5S1','OQV22897.1':'A0A1W0X5M8'}
results={}; cds_fastas=[]; protein_fastas=[]; exon_rows=[]
for record in ET.parse(RAW / files[0]).getroot().findall('GBSeq'):
    pid=record.findtext('GBSeq_accession-version'); accession=pairs[pid]
    protein=record.findtext('GBSeq_sequence').upper()
    feature=next(f for f in record.findall('GBSeq_feature-table/GBFeature') if f.findtext('GBFeature_key')=='CDS')
    q=quals(feature); coded_by=q['coded_by']
    assert 'complement' not in coded_by
    exons=[tuple(map(int,m)) for m in re.findall(re.escape(SCAFFOLD)+r':(\d+)\.\.(\d+)',coded_by)]
    assert len(exons)==7
    cds=''.join(nt[a-START:b-START+1] for a,b in exons)
    trans=translate(cds)
    assert trans.endswith('*') and '*' not in trans[:-1] and trans[:-1]==protein
    up=json.loads((RAW/(accession+'.uniprot.json')).read_text())
    assert up['primaryAccession']==accession and up['sequence']['value']==protein
    nf=next(f for f in nt_record.findall('GBSeq_feature-table/GBFeature') if f.findtext('GBFeature_key')=='CDS' and quals(f).get('protein_id')==pid)
    assert quals(nf)['translation']==protein
    local_exons=[tuple(map(int,m)) for m in re.findall(r'(\d+)\.\.(\d+)',nf.findtext('GBFeature_location'))]
    assert [(a+START-1,b+START-1) for a,b in local_exons]==exons
    assert q['locus_tag']==gene['locus_tag']=='BV898_03327'
    results[accession]={'protein_id':pid,'locus_tag':q['locus_tag'],'protein_length':len(protein),'cds_length_including_stop':len(cds),'start_codon':cds[:3],'stop_codon':cds[-3:],'strand':'+','cds_exons_1_based_inclusive':exons,'cds_ambiguous_bases':len(re.sub('[ACGT]','',cds)),'translation_matches_ncbi_protein':True,'translation_matches_uniprot':True,'translation_matches_nucleotide_annotation':True,'cds_sha256':sha(cds.encode()),'protein_sequence_sha256':sha(protein.encode()),'uniprot_entry_audit':up['entryAudit'],'protein_record_update_date':record.findtext('GBSeq_update-date'),'annotation_method':q['note']}
    cds_fastas.append(fasta(pid+'|'+accession+'|BV898_03327|CDS|'+SCAFFOLD+'|strand=+|includes_stop',cds))
    protein_fastas.append(fasta(accession+'|'+pid+'|BV898_03327|translated_from_versioned_genomic_CDS',protein))
    for i,(a,b) in enumerate(exons,1): exon_rows.append([accession,pid,'BV898_03327',SCAFFOLD,'+',i,a,b,b-a+1])
first=json.loads((RAW/'A0A1W0X5S1.uniprot.json').read_text())['sequence']['value']
second=json.loads((RAW/'A0A1W0X5M8.uniprot.json').read_text())['sequence']['value']
changes=[{'operation':tag,'X5S1_start_1_based':i+1,'X5S1_end_1_based':j,'X5S1_residues':first[i:j],'X5M8_start_1_based':k+1,'X5M8_end_1_based':l,'X5M8_residues':second[k:l]} for tag,i,j,k,l in difflib.SequenceMatcher(None,first,second,autojunk=False).get_opcodes() if tag!='equal']
assert len(changes)==1 and changes[0]['X5S1_residues']=='E' and changes[0]['X5M8_residues']=='GPITCK'
asm=json.loads((RAW/files[4]).read_text())['result']['1063911']
summ=json.loads((RAW/files[5]).read_text())['result']['1174754726']
assert asm['assemblyaccession']=='GCA_002082055.1' and asm['wgs']=='MTYJ01'
assert summ['accessionversion']==SCAFFOLD
boundary={'changed_coding_exon_number':4,'X5S1_interval_1_based_inclusive':[725933,726089],'X5M8_interval_1_based_inclusive':[725933,726104],'added_genomic_interval_1_based_inclusive':[726090,726104],'added_nucleotides':nt[726090-START:726104-START+1],'added_nt_length':15,'interpretation':'alternative 5-prime splice donor boundary in the deposited plus-strand model; not an alternative acceptor','X5S1_donor_dinucleotide':nt[726090-START:726091-START+1],'X5M8_donor_dinucleotide':nt[726105-START:726106-START+1],'shared_acceptor_dinucleotide':nt[726200-START:726201-START+1],'coding_nt_count_through_X5S1_exon4':sum(b-a+1 for a,b in results['A0A1W0X5S1']['cds_exons_1_based_inclusive'][:4]),'protein_changes':changes,'net_protein_length_difference':5}
assert boundary['X5S1_donor_dinucleotide']==boundary['X5M8_donor_dinucleotide']=='GT' and boundary['shared_acceptor_dinucleotide']=='AG'
put(DERIVED/'BV898_03327.deposited_CDS.fasta',''.join(cds_fastas))
put(DERIVED/'BV898_03327.translated_proteins.fasta',''.join(protein_fastas))
put(DERIVED/'MTYJ01000015.1_723181-728004.fasta',fasta(SCAFFOLD+':723181-728004|1_based_inclusive|strand=+',nt))
put(DERIVED/'BV898_03327.cds_exons.tsv','uniprot\tprotein_id\tlocus\tscaffold\tstrand\texon\tstart_1based\tend_1based\tlength\n'+'\n'.join('\t'.join(map(str,row)) for row in exon_rows)+'\n')
report={'operation_class':'same-input recomputation from newly retrieved versioned public sequence records; not a replay of missing historical scripts','date_utc':datetime.datetime.now(datetime.timezone.utc).date().isoformat(),'input_record_hashes':{f:sha((RAW/f).read_bytes()) for f in files},'script_sha256':sha(Path(__file__).read_bytes()),'source_links':{'X5S1':'https://www.ncbi.nlm.nih.gov/protein/OQV22896.1','X5M8':'https://www.ncbi.nlm.nih.gov/protein/OQV22897.1','scaffold':'https://www.ncbi.nlm.nih.gov/nuccore/MTYJ01000015.1','assembly':'https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_002082055.1/','uniprot_X5S1':'https://www.uniprot.org/uniprotkb/A0A1W0X5S1/entry','uniprot_X5M8':'https://www.uniprot.org/uniprotkb/A0A1W0X5M8/entry'},'organism':source['organism'],'strain':source['strain'],'taxon':source['db_xref'],'source_record_comment':nt_record.findtext('GBSeq_comment'),'assembly':{'accession':asm['assemblyaccession'],'assembly_name_as_deposited':asm['assemblyname'],'wgs_prefix':asm['wgs'],'bioproject':'PRJNA360553','biosample':asm['biosampleaccn']},'scaffold':{'accession':SCAFFOLD,'length':summ['slen'],'slice_start_1_based':START,'slice_end_1_based':END,'slice_length':len(nt),'slice_ambiguous_bases':len(re.sub('[ACGT]','',nt)),'bases_from_cds_end_to_scaffold_end':summ['slen']-727004},'gene_note_crosswalk':{'locus_tag':gene['locus_tag'],'note':gene['note'],'status':'directly supported by retrieved nucleotide gene annotation','individual_bHd16413_1_2_assignments':'not present in retrieved annotations; not sequence-verified across nHd.3.0 expression and nHd_3.1 genome builds'},'models':results,'alternative_boundary':boundary,'tool_record':{'Life_Sciences_Databases_plugin':'documented NCBI-Entrez and UniProt scripts attempted with /usr/bin/python3, both failed','exact_error':{'ok':False,'error':{'code':'missing_dependency','message':"`requests` is required: No module named 'requests'"},'warnings':[]},'successful_retrieval':'Python standard-library urllib to NCBI E-Utilities and UniProt REST; separately credited local retrieval','calculation':'this script parses coded_by intervals, joins genomic bases, uses standard genetic code 1, checks deposited translations and compares sequences','new_installations':False},'limitations':['Original S1/S3/S4 reports and historical genome-audit scripts remain unavailable.','The public records retain versioned sequence identifiers, but current annotation metadata were retrieved now.','No raw RNA reads were aligned; annotated splicing and translation are not experimental validation.','The 4,824-bp slice verifies no ambiguous bases within that slice and both CDSs only; no 50-kb flank gap claim is reproduced.','mRNA and gene annotations have uncertain terminal bounds (< and >); complete coding models do not establish transcript ends.','bHd16413 is directly annotated at the gene, but individual .1/.2 transcript identities across expression/genome reference builds remain unverified.']}
json_path=OUT/'GENOME_RECORD_AUDIT_20260915.json'
put(json_path,json.dumps(report,indent=2)+'\n')
md='''# CAD locus audit: deposited sequence records

Public records retrieved 15 September 2026. This is a focused recomputation from versioned public records and current annotation metadata. Original historical audit scripts were not recovered.

## Verified findings

- The two deposited proteins are **A0A1W0X5S1 / OQV22896.1, 423 aa**, and **A0A1W0X5M8 / OQV22897.1, 428 aa**. They share locus **BV898_03327**. Both CDSs translate exactly to their respective NCBI protein and UniProt sequence, with **ATG** starts and **TAG** stops. CDS lengths including stop are **1,272 nt** and **1,287 nt**. [NCBI X5S1](https://www.ncbi.nlm.nih.gov/protein/OQV22896.1), [NCBI X5M8](https://www.ncbi.nlm.nih.gov/protein/OQV22897.1).
- Both coding models occupy **724181–727004 on the plus strand of MTYJ01000015.1**, using seven coding exons. Coordinates in this report are **1-based, inclusive**. The fetched genomic slice contains no ambiguous bases, and both extracted CDSs contain no ambiguous bases. [Nucleotide record](https://www.ncbi.nlm.nih.gov/nuccore/MTYJ01000015.1).
- The nucleotide annotation explicitly gives the gene **BV898_03327** the note **bHd16413**. The individual transcript labels **bHd16413.1 / bHd16413.2** are absent from these retrieved annotations. Their exact assignment across the expression reference **nHd.3.0** and genomic reference **nHd_3.1** is not sequence-verified here.
- Source records name **Hypsibius exemplaris, strain Z151**. The nucleotide comment states that the record was renamed from **H. dujardini strain Z151** in March 2022 following strain reclassification. These labels describe the same deposited lineage. [Nucleotide record](https://www.ncbi.nlm.nih.gov/nuccore/MTYJ01000015.1).
- NCBI assembly metadata identifies **GCA_002082055.1**, deposited assembly name **nHd_3.1**, WGS prefix **MTYJ01**, BioProject **PRJNA360553**, and BioSample **SAMN06212344**. Scaffold metadata gives **832,999 bp**, placing the coding end **105,995 bp** from the scaffold end. [Assembly record](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_002082055.1/).

## Correction to the updated brief

The retrieved annotations show an **alternative donor boundary**, not an alternative acceptor. On this plus-strand locus, exon 4 ends at **726089** in X5S1 and **726104** in X5M8; the next exon starts at **726202** in both. The longer model retains **726090–726104**, the 15-nt sequence **GTCCGATTACTTGCA**. Both alternative introns start **GT** and share the downstream **AG** acceptor. This is an annotation-based splice interpretation, not a raw-read junction validation. [NCBI X5S1](https://www.ncbi.nlm.nih.gov/protein/OQV22896.1), [NCBI X5M8](https://www.ncbi.nlm.nih.gov/protein/OQV22897.1).

The protein difference is **X5S1 E199 → X5M8 GPITCK199–204**, a net five-amino-acid length increase coupled to a changed junction residue. Calling it only a five-amino-acid insertion obscures that change. The coding sequence through the shorter fourth exon is 595 nt, ending after the first base of a codon; changing the donor therefore changes the split codon while preserving the downstream frame.

## Coding exons

| Exon | X5S1 coordinates | X5M8 coordinates |
|---|---|---|
'''
for i,(a,b) in enumerate(results['A0A1W0X5S1']['cds_exons_1_based_inclusive']):
    c,d=results['A0A1W0X5M8']['cds_exons_1_based_inclusive'][i]
    md+=f'| {i+1} | {a}–{b} | {c}–{d} |\n'
md+='''
## Actual operations and limitations

The documented Life Sciences Databases NCBI-Entrez and UniProt scripts were attempted with `/usr/bin/python3`; both returned `missing_dependency: requests is required: No module named requests`. No package was installed. Successful retrieval used Python's standard-library `urllib`, and the locus/CDS calculations used this saved local script. This is not claimed as successful plugin execution.

Only a **4,824-bp genomic slice** and small protein/assembly metadata records were retrieved. The historical claim of no assembly gap within 50 kb was not tested. CDS translation confirms annotation consistency; it does not prove expression, transcript ends, protein production, or enzyme activity. The gene/mRNA annotations have uncertain terminal bounds, whereas the CDS intervals have explicit starts and stops. No raw reads were analyzed.

## Reproducibility files

- Raw records and full URLs, retrieval times and SHA-256 hashes: `analysis/locus/inputs/retrieval_ledger_20260915.json`.
- Machine-readable findings and input hashes: `outputs/GENOME_RECORD_AUDIT_20260915.json`.
- CDS/protein FASTAs, bounded genomic FASTA, and exon TSV: `analysis/locus/derived/`.
- Exact local calculation: `/usr/bin/python3 scripts/audit_genome_records_20260915.py`.

Every original report, prior output, and retrieved source record was preserved.
'''
put(OUT/'GENOME_RECORD_AUDIT_20260915.md',md)
print(json.dumps({'outputs':[str(json_path),str(OUT/'GENOME_RECORD_AUDIT_20260915.md')],'models':results,'alternative_boundary':boundary,'gene_note':gene},indent=2))
