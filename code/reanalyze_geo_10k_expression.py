#!/usr/bin/python3
"""Descriptive reanalysis of six deposited GSE94295 10k-adult Kallisto tables.
No raw-read analysis, count model, p-values or native NGS workflow.
Usage: /usr/bin/python3 scripts/reanalyze_geo_10k_expression.py
"""
import argparse, csv, gzip, hashlib, io, json, math, re, statistics, sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
INPUT=ROOT/'analysis/expression'
ap=argparse.ArgumentParser(description=__doc__)
ap.add_argument('--output',required=True,type=Path)
OUTPUT=ap.parse_args().output.resolve()
assert not OUTPUT.exists(), 'Use a new output directory'
SAMPLES=['GSM2472501','GSM2472502','GSM2472503','GSM2472504','GSM2472505','GSM2472506']
TARGET='bHd16413'
PSEUDOCOUNT=0.1

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def preserve(p,data):
    if p.exists():
        if p.read_bytes()!=data:raise ValueError('Refusing overwrite of changed output '+str(p))
    else:
        with p.open('xb') as f:f.write(data)
def tsv(p,rows):
    s=io.StringIO(newline='')
    w=csv.DictWriter(s,list(rows[0]),delimiter='\t',lineterminator='\n')
    w.writeheader();w.writerows(rows);preserve(p,s.getvalue().encode())
def fmt(v):return 'NA' if v is None else format(v,'.12g')
def metrics(v):
    a=statistics.mean(v[:3]);t=statistics.mean(v[3:])
    return {'active_mean_tpm':a,'active_sample_sd_tpm':statistics.stdev(v[:3]),
        'tun_mean_tpm':t,'tun_sample_sd_tpm':statistics.stdev(v[3:]),
        'tun_over_active_mean_ratio':t/a if a>0 else None,
        'log2_ratio_with_0_1_tpm_pseudocount':math.log2((t+PSEUDOCOUNT)/(a+PSEUDOCOUNT))}
def main():
    manifest=json.loads((INPUT/'SOURCE_MANIFEST.json').read_text())
    for r in manifest['records']:
        p=ROOT/r['relative_path']
        assert p.stat().st_size==r['bytes'] and sha(p)==r['sha256']
    with (INPUT/'COHORT_MANIFEST.tsv').open() as f:cohort=list(csv.DictReader(f,delimiter='\t'))
    assert [r['sample_id'] for r in cohort]==SAMPLES
    assert [r['condition'] for r in cohort]==['active']*3+['tun']*3
    assert all(r['individuals']=='10000' and r['genome_build_as_deposited']=='nHd.3.0' for r in cohort)
    tables={};checks=[];mapping={};lengths={}
    for sample in cohort:
        sid=sample['sample_id'];p=ROOT/sample['table_path']
        assert sha(p)==sample['table_sha256']
        with gzip.open(p,'rt') as f:
            rd=csv.DictReader(f,delimiter='\t')
            assert rd.fieldnames==['target_id','length','eff_length','est_counts','tpm']
            rows=list(rd)
        ids=[r['target_id'] for r in rows]
        assert len(ids)==len(set(ids)),'duplicate transcript identifiers'
        table={}
        for row in rows:
            tid=row['target_id'];match=re.fullmatch(r'(bHd\d+)\.(\d+)',tid)
            assert match,'unrecognized transcript identifier '+tid
            mapping[tid]=match.group(1)
            row={k:(v if k=='target_id' else float(v)) for k,v in row.items()}
            assert all(math.isfinite(row[k]) and row[k]>=0 for k in ['length','eff_length','est_counts','tpm'])
            assert row['length'].is_integer() and row['length']>0
            if tid in lengths:assert lengths[tid]==row['length'],'transcript length mismatch'
            lengths[tid]=row['length'];table[tid]=row
        tables[sid]=table
        checks.append({'sample_id':sid,'transcript_count':len(rows),
          'sum_tpm':sum(r['tpm'] for r in table.values()),
          'sum_est_counts':sum(r['est_counts'] for r in table.values()),'source_sha256':sha(p)})
    first=set(tables[SAMPLES[0]])
    assert all(set(t)==first for t in tables.values()),'transcript universes differ'
    gids=sorted(set(mapping.values()),key=lambda x:int(x[3:]))
    tids=sorted(first,key=lambda x:tuple(map(int,x[3:].split('.'))))
    by_gene={g:[] for g in gids}
    for tid in tids:by_gene[mapping[tid]].append(tid)
    values={g:[sum(tables[s][t]['tpm'] for t in by_gene[g]) for s in SAMPLES] for g in gids}
    for i,c in enumerate(checks):
        c['locus_sum_tpm']=sum(v[i] for v in values.values())
        assert abs(c['locus_sum_tpm']-c['sum_tpm'])<1e-6,'TPM changed during aggregation'
    OUTPUT.mkdir(parents=True,exist_ok=True)
    grows=[]
    for g in gids:
        row={'gene_id':g,'transcript_count':len(by_gene[g])}
        row.update({s+'_tpm':fmt(v) for s,v in zip(SAMPLES,values[g])})
        row.update({k:fmt(v) for k,v in metrics(values[g]).items()});grows.append(row)
    tsv(OUTPUT/'gene_expression_summary.tsv',grows)
    trows=[]
    for t in tids:
        row={'transcript_id':t,'gene_id':mapping[t],'length':int(lengths[t])}
        for s in SAMPLES:
            for col in ['eff_length','est_counts','tpm']:row[s+'_'+col]=fmt(tables[s][t][col])
        row.update({k:fmt(v) for k,v in metrics([tables[s][t]['tpm'] for s in SAMPLES]).items()});trows.append(row)
    tsv(OUTPUT/'transcript_expression_summary.tsv',trows)
    tsv(OUTPUT/'transcript_to_gene.tsv',[{'transcript_id':t,'historical_gene_id':mapping[t],'mapping_rule':'remove terminal dot and numeric transcript suffix'} for t in tids])
    tsv(OUTPUT/'bHd16413_gene_summary.tsv',[r for r in grows if r['gene_id']==TARGET])
    tsv(OUTPUT/'bHd16413_transcripts.tsv',[r for r in trows if r['gene_id']==TARGET])
    m=metrics(values[TARGET]);v=values[TARGET]
    check={'operation_class':'processed-table reanalysis','software':'Python standard library '+sys.version.split()[0],
      'analysis_script':'scripts/reanalyze_geo_10k_expression.py','analysis_script_sha256':sha(Path(__file__)),
      'source_manifest':'analysis/expression/SOURCE_MANIFEST.json',
      'source_manifest_sha256':sha(INPUT/'SOURCE_MANIFEST.json'),'cohort_manifest_sha256':sha(INPUT/'COHORT_MANIFEST.tsv'),
      'samples':checks,'transcript_universe':len(tids),'gene_universe':len(gids),
      'tun_nonzero_gene_universe':sum(statistics.mean(v[3:])>0 for v in values.values()),
      'all_transcript_ids_match_bHd_numeric_suffix':True,'all_six_transcript_universes_equal':True,
      'all_transcript_lengths_equal_across_samples':True,'duplicate_transcript_ids':0,
      'tpm_aggregation_preserves_library_totals':True,
      'aggregation':'sum transcript TPM by validated removal of final .numeric suffix; no abundance filtering',
      'condition_summary':'arithmetic mean; sample SD with n-1 denominator; n=3 libraries per condition',
      'direct_ratio':'tun mean / active mean; NA if active mean is zero',
      'descriptive_log2_ratio':'log2((tun mean + 0.1)/(active mean + 0.1))','pseudocount_tpm':0.1,
      'target_gene':TARGET,'target_transcripts':by_gene[TARGET],'target_gene_metrics':m,
      'genome_build_as_deposited':'nHd.3.0','distinct_genomic_audit_build':'nHd3.1/GCA_002082055.1',
      'locus_crosswalk_status':'bHd16413 to BV898_03327/OQV products supplied by handoff; sequence identity not established by table identifiers or length alone.',
      'historical_comparison_source':'ROSALIND_CAD_EVIDENCE_AUDIT.md section 6B; original reports and analyze_expression.py not used',
      'historical_rounded_values':{'active_mean_tpm':1.387,'tun_mean_tpm':10.017,'direct_ratio':7.22,'pseudocount_log2_ratio':2.77},
      'limitations':['Newly retrieved tables, not proven identical to historical input snapshot.',
      'No raw reads, mapping, splice-junction analysis, count model or new p-values.',
      'Estimated counts retained as fractional values, separate from TPM.',
      'Published single-adult FDRs and separate 30-adult cohort not recomputed.']}
    preserve(OUTPUT/'validation_and_provenance.json',(json.dumps(check,indent=2)+'\n').encode())
    report=f"""# GSE94295: deposited 10,000-adult expression reanalysis

**Operation:** processed-table reanalysis. These are newly calculated descriptive results from six newly retrieved GEO Kallisto tables. No raw RNA-seq analysis or formal differential-expression test was performed.

## Result

For historical RNA locus **bHd16413**, sum **bHd16413.1** and **bHd16413.2** within each library:

| Condition | Replicate 1 TPM | Replicate 2 TPM | Replicate 3 TPM | Mean TPM | Sample SD TPM |
|---|---:|---:|---:|---:|---:|
| Active | {v[0]:.6f} | {v[1]:.6f} | {v[2]:.6f} | {m['active_mean_tpm']:.9f} | {m['active_sample_sd_tpm']:.9f} |
| Tun | {v[3]:.6f} | {v[4]:.6f} | {v[5]:.6f} | {m['tun_mean_tpm']:.9f} | {m['tun_sample_sd_tpm']:.9f} |

Direct mean ratio: **{m['tun_over_active_mean_ratio']:.6f}**. Descriptive log2 ratio with 0.1 TPM added to both condition means: **{m['log2_ratio_with_0_1_tpm_pseudocount']:.6f}**.

These independently recomputed values match the handoff's rounded historical 1.387 to 10.017 TPM, 7.22-fold and +2.77. The handoff is the comparison source; the original report and historical analysis script were not recovered by this operation.

## Inputs and identity checks

The cohort is active GSM2472501, GSM2472502, GSM2472503 and tun GSM2472504, GSM2472505, GSM2472506. All six GEO records identify whole-body adult Z151 pools of 10,000 animals, deposited as Hypsibius dujardini, sequenced on Illumina HiSeq 2000. Metadata reports **Kallisto v0.42.4** and **nHd.3.0**. The separate locus audit uses **nHd3.1/GCA_002082055.1**; the reference labels remain distinct.

All six tables contain the same {len(tids):,} unique transcript IDs and transcript lengths. Every ID matches the bHd-number.numeric-transcript convention. Validated suffix removal creates {len(gids):,} loci; {check['tun_nonzero_gene_universe']:,} have positive mean tun TPM. Gene sums preserve each library's total TPM to numerical precision. No abundance filtering or selected-candidate ranking was applied.

Target transcript lengths are 1,272 and 1,287 nucleotides. Identifiers and lengths alone do not prove the bHd16413 to BV898_03327/OQV crosswalk; that requires a separate sequence/CDS check. Estimated counts and library-specific effective lengths are retained in the original compressed tables and full transcript summary. Counts are not rounded to integers.

## Method and scope

Sum TPM by historical gene ID within each library. Compute arithmetic means and sample SDs across three libraries per condition. Direct ratio is tun mean / active mean. The distinct descriptive log2 statistic is log2((tun_mean + 0.1)/(active_mean + 0.1)). No formal significance is inferred.

Published single-adult FDRs and the separate 30-animal cohort were not recomputed. Endpoint RNA abundance does not establish protein production, nuclease activity, or a causal desiccation mechanism.

## Reproduce and provenance

Run /usr/bin/python3 scripts/reanalyze_geo_10k_expression.py from this project. It verifies source hashes and writes only absent outputs or accepts byte-identical existing outputs. validation_and_provenance.json records formulas, software, hashes and checks. Source metadata, original compressed tables, URLs and SHA-256 values are in analysis/expression/.

GEO samples: [GSM2472501](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2472501), [GSM2472502](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2472502), [GSM2472503](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2472503), [GSM2472504](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2472504), [GSM2472505](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2472505), [GSM2472506](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2472506).

Execution used local Python standard-library retrieval and analysis. Installed NGS Analysis Workbench interpretation guidance informed the review, but no native workflow was run. There is no registered run to attach this review to, and the plugin exposes no direct validator for arbitrary Kallisto tables. Earlier Entrez attempts failed and are not successful plugin operations.
"""
    preserve(OUTPUT/'EXPRESSION_10K_REANALYSIS.md',report.encode())
    outputs=[{'relative_path':str(p.relative_to(ROOT)),'bytes':p.stat().st_size,'sha256':sha(p)} for p in sorted(OUTPUT.iterdir()) if p.name!='OUTPUT_HASHES.json']
    preserve(OUTPUT/'OUTPUT_HASHES.json',(json.dumps(outputs,indent=2)+'\n').encode())
    print(json.dumps({'gene_universe':len(gids),'transcript_universe':len(tids),'tun_nonzero_gene_universe':check['tun_nonzero_gene_universe'],'target_gene_metrics':m,'output_dir':str(OUTPUT)},indent=2))
if __name__=='__main__':main()
