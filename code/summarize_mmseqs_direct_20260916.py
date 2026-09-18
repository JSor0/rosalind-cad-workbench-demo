#!/usr/bin/python3
"""Offline, create-only audit of an existing MMseqs search; never launches a search."""
import argparse
import collections
import csv
import datetime
import hashlib
import io
import json
import math
from pathlib import Path

RUN = 'mmseqs_reproduction_20260916_uniprot_current'
FIELDS = 'query,target,evalue,bits,pident,nident,alnlen,mismatch,gapopen,qstart,qend,qlen,tstart,tend,tlen,qcov,tcov,cigar,qseq,tseq,qaln,taln,qheader,theader'.split(',')
PANEL = ['O76075', 'O54788', 'A0A1W0WC24', 'A0A9X6RM91', 'A0A9X6RN34', 'A0A9X6NI93', 'A0A1W0X5S1', 'A0A1W0X6E7']
LABELS = dict(zip(PANEL, ['human DFFB', 'mouse DFFB', 'WC24', 'RM91', 'RN34', 'NI93', 'X5S1', 'X6E7']))
LABELS['A0A1W0X5M8'] = 'X5M8 alternative product'
MISSING = 'not reported under this search at E<=10'
INT_FIELDS = 'bits nident alnlen mismatch gapopen qstart qend qlen tstart tend tlen'.split()
NUMERIC_FIELDS = 'evalue pident qcov tcov'.split()
METRICS = FIELDS[:18]

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def fasta(path):
    records, headers = {}, {}
    current = None
    for line in path.read_text().splitlines():
        if line.startswith('>'):
            current = line[1:].split()[0]
            assert current not in records, ('duplicate FASTA identifier', current)
            records[current], headers[current] = '', line[1:]
        elif line.strip():
            assert current is not None
            records[current] += line.strip()
    return records, headers

def metrics(row):
    result = {key: row[key] for key in METRICS}
    for key in INT_FIELDS:
        result[key] = int(result[key])
    for key in NUMERIC_FIELDS:
        result[key] = float(result[key])
    result['reported_evalue_text'] = row['evalue']
    result['reported_pident_text'] = row['pident']
    result['identity_from_counts_percent'] = 100 * result['nident'] / result['alnlen']
    return result

def audit(project):
    out = project / 'analysis/sequence'
    prep_path, receipt_path, raw = [out / x for x in ['PREPARATION.json', 'RUN_RECEIPT.json', 'direct_e10.tsv']]
    prep, receipt = json.loads(prep_path.read_text()), json.loads(receipt_path.read_text())
    assert receipt['status'] == 'completed' and receipt['returncode'] == 0
    assert receipt['preparation_sha256'] == prep['original_source_sha256']  # Original linked receipt, before explicit path sanitization
    assert sha(raw) == receipt['output']['sha256']
    assert raw.stat().st_size == receipt['output']['bytes']
    assert receipt['argv'] == prep['command']['argv']
    query_path, target_path = [project/receipt['argv'][i] for i in [2, 3]]
    for path in [query_path, target_path]:
        assert sha(path) == prep['prepared_files'][str(path.relative_to(project))]['sha256']
    assert len(prep['mmseqs_executable_sha256']) == 64  # Original binary identity is recorded; binary is not bundled or reverified here
    queries, qheaders = fasta(query_path)
    targets, theaders = fasta(target_path)
    assert list(queries) == PANEL and len(targets) == 20834
    assert len(set(targets)) == 20834 and sum(map(len, targets.values())) == 8999988
    assert all(queries[q] == targets[q] for q in PANEL)
    lines = raw.read_text().splitlines()
    header = lines[0].removeprefix('#').split('\t')
    assert header == FIELDS, header
    assert all(len(line.split('\t')) == len(FIELDS) for line in lines[1:])
    rows = [dict(zip(FIELDS, line.split('\t'))) for line in lines[1:]]
    assert len(rows) == 75 and set(r['query'] for r in rows) == set(PANEL)
    pairs, max_identity_difference, max_coverage_difference = {}, 0.0, 0.0
    checks = collections.Counter()
    for row_number, row in enumerate(rows, 2):
        q, t = row['query'], row['target']
        assert q in queries and t in targets, row_number
        assert (q, t) not in pairs, ('duplicate directed pair', row_number)
        pairs[q, t] = row
        for side, acc, seqs, headers in [('q', q, queries, qheaders), ('t', t, targets, theaders)]:
            start, end, length = [int(row[side + suffix]) for suffix in ['start', 'end', 'len']]
            assert 1 <= start <= end <= length == len(seqs[acc]), row_number
            assert row[side+'seq'] == seqs[acc], row_number
            assert row[side+'header'] == headers[acc], row_number
            assert row[side+'aln'].replace('-', '') == seqs[acc][start-1:end], row_number
            cov_error = abs((end-start+1)/length - float(row[side+'cov']))
            assert cov_error <= 0.0005000001, (row_number, 'coverage', cov_error)
            max_coverage_difference = max(max_coverage_difference, cov_error)
        a, b = row['qaln'], row['taln']
        assert len(a) == len(b) == int(row['alnlen']), row_number
        assert all(x != '-' or y != '-' for x, y in zip(a, b)), row_number
        identity = sum(x == y and x != '-' for x, y in zip(a, b))
        mismatch = sum(x != y and x != '-' and y != '-' for x, y in zip(a, b))
        gapopens = sum(x == '-' and (i == 0 or a[i-1] != '-') for i, x in enumerate(a)) + sum(x == '-' and (i == 0 or b[i-1] != '-') for i, x in enumerate(b))
        assert [identity, mismatch, gapopens] == [int(row[k]) for k in ['nident', 'mismatch', 'gapopen']], row_number
        # This result uses floor(identity/alnlen*1000)/10 percent, not rounding to 3 decimals.
        assert abs(math.floor(identity/len(a)*1000)/10 - float(row['pident'])) < 1e-9, row_number
        max_identity_difference = max(max_identity_difference, 100*identity/len(a)-float(row['pident']))
        ops = ['D' if x == '-' else 'I' if y == '-' else 'M' for x, y in zip(a, b)]
        grouped = []
        for op in ops:
            if grouped and grouped[-1][1] == op:
                grouped[-1][0] += 1
            else:
                grouped.append([1, op])
        cigar = ''.join(str(count)+op for count, op in grouped)
        assert cigar == row['cigar'], row_number
        assert math.isfinite(float(row['evalue'])) and 0 <= float(row['evalue']) <= 10, row_number
        for check in ['full_sequences_and_headers_match_fasta', 'sequence_lengths_and_coordinate_bounds', 'ungapped_alignment_matches_source_spans', 'alignment_lengths_no_double_gaps', 'identities_mismatches_gap_opens', 'reported_identity_precision', 'coverage_rounding', 'exact_cigar_reconstruction', 'evalue_within_reporting_cutoff']:
            checks[check] += 1
    selfhits = []
    for acc in PANEL:
        row = pairs[acc, acc]
        length = len(queries[acc])
        assert row['qaln'] == row['taln'] == queries[acc]
        assert [int(row[k]) for k in ['qstart', 'tstart']] == [1, 1]
        assert [int(row[k]) for k in ['qend', 'tend', 'nident', 'alnlen']] == [length]*4
        assert float(row['pident']) == 100 and float(row['qcov']) == float(row['tcov']) == 1
        selfhits.append(metrics(row))
    pair_columns = ['query', 'target', 'query_label', 'target_label', 'status'] + FIELDS[2:18]
    buf = io.StringIO(newline='')
    writer = csv.DictWriter(buf, fieldnames=pair_columns, delimiter='\t', lineterminator='\n')
    writer.writeheader()
    for q in PANEL:
        for t in PANEL:
            row = pairs.get((q, t))
            entry = {'query': q, 'target': t, 'query_label': LABELS[q], 'target_label': LABELS[t], 'status': 'reported' if row else MISSING}
            if row:
                entry.update({k: row[k] for k in FIELDS[2:18]})
            writer.writerow(entry)
    pair_text = buf.getvalue()
    selected = [('O76075','O54788'), ('O54788','O76075'), ('O76075','A0A1W0WC24'), ('O76075','A0A9X6RM91'), ('O76075','A0A9X6RN34'), ('A0A9X6RN34','O76075'), ('A0A9X6RN34','A0A9X6NI93'), ('A0A9X6NI93','A0A9X6RN34'), ('A0A1W0WC24','A0A9X6NI93'), ('A0A1W0X5S1','A0A1W0X6E7'), ('A0A1W0X6E7','A0A1W0X5S1'), ('A0A1W0X5S1','A0A1W0X5M8'), ('A0A1W0X5S1','A0A9X6RM91'), ('A0A1W0X5S1','A0A9X6NI93')]
    absent_reference_pairs = [(q,t) for q,t in [('A0A1W0X5S1','O76075'), ('O76075','A0A1W0X5S1'), ('A0A1W0X5S1','O54788'), ('O54788','A0A1W0X5S1')] if (q,t) not in pairs]
    assert len(absent_reference_pairs) == 4
    n_panel = sum(t in PANEL for q,t in pairs)
    limits = [
        'The search is a fresh, guided calculation using archived queries and a current frozen full proteome. It is not blind rediscovery.',
        'Local MMseqs2 execution is not native Rosalind/Life Sciences plugin execution. A sequence-viewer display is a separate operation.',
        'MMseqs2 release 18-8cc5c differs from historical build cb12a2d75a9808ee61721029d064a7bd80af6fec; unspecified defaults can differ.',
        'Current UniProt release 2026_03 has the same 20,832-protein count as reported historically. Exact historical target identity remains unknown because the original target hash/snapshot was not recovered.',
        'All eight archived query sequences exactly match their current target records; this does not establish identity of the other 20,826 targets across releases.',
        'S1 section 14.6 (PDF page 26) clips the output-field list after qle. The complete output fields, --threads 4 and -a 1 are explicit fresh operational choices.',
        'Missing pairs are not reported under this search at E<=10. Missingness does not imply E>10 or non-homology, because no alignment/E-value was returned.',
        'No iterative profile search was executed here. The historical profile results remain source-reported, not reproduced by this direct search.',
        'Percent identity is local alignment identity. It is not forced global, multiple-alignment or structural-pair identity.',
        'Reported E=0.000E+00 for the RN34 self hit is preserved as engine output; it is not a claim of an exact zero statistical probability.',
        'These sequence results do not establish nuclease activity, domain-loss history, ICAD independence or biochemical function.'
    ]
    summary = {
        'schema_version': 1, 'generated_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'analysis_class': 'independent offline validation and summary of the completed fresh proteome-calibrated direct sequence search',
        'parser': {'path': str(Path(__file__).resolve()), 'sha256': sha(Path(__file__))},
        'input_files': {str(path): {'sha256': sha(path), 'bytes': path.stat().st_size} for path in [raw, prep_path, receipt_path, query_path, target_path]},
        'run_started_utc': receipt['started_utc'], 'run_completed_utc': receipt['completed_utc'], 'run_elapsed_seconds': receipt['elapsed_seconds'], 'returncode': receipt['returncode'],
        'actual_execution': receipt['actual_execution'], 'native_plugin_used_for_search': False,
        'mmseqs_release': prep['mmseqs_release'], 'mmseqs_version': prep['mmseqs_reported_version'], 'historical_build': prep['historical_mmseqs_build'],
        'uniprot_release': prep['uniprot_release'], 'uniprot_release_date': prep['uniprot_release_date'],
        'historical_target_exact_identity_verified': False, 'all_eight_archived_queries_match_current_records': True,
        'iterative_profiles_executed': False, 'query_order': PANEL, 'query_labels': {q:LABELS[q] for q in PANEL},
        'query_count': len(queries), 'target_count': len(targets), 'target_total_residues': sum(map(len, targets.values())),
        'raw_header': header, 'header_has_leading_hash': lines[0].startswith('#'), 'data_row_count': len(rows), 'unique_directed_pairs': len(pairs), 'distinct_target_count': len({t for q,t in pairs}),
        'rows_per_query': {q:sum(r['query']==q for r in rows) for q in PANEL},
        'directed_panel_pairs': {'total': 64, 'reported': n_panel, 'not_reported': 64-n_panel, 'off_panel_rows': len(rows)-n_panel, 'missing_status_text': MISSING},
        'validation': {'status': 'passed', 'all_rows_valid': True, 'checks_passed_by_row': dict(checks), 'self_hits_full_length_100_percent_identity': len(selfhits), 'maximum_coverage_rounding_difference': max_coverage_difference, 'maximum_reported_identity_truncation_percentage_points': max_identity_difference, 'reported_identity_convention_verified': 'floor(nident/alnlen*1000)/10 percent for all 75 rows', 'cigar_convention': 'M = paired residues; I = query residue with target gap; D = target residue with query gap'},
        'self_hits': selfhits, 'key_reported_pairs': [metrics(pairs[key]) for key in selected],
        'x5s1_human_mouse_missing_directed_pairs': [{'query':q, 'target':t, 'status':MISSING, 'evalue':None} for q,t in absent_reference_pairs],
        'limitations': limits,
        'generated_pair_table': {'path':str(out/'PANEL_PAIRS.tsv'), 'sha256':hashlib.sha256(pair_text.encode()).hexdigest(), 'rows':64}
    }
    md = [
        '# Fresh proteome-calibrated sequence results', '',
        'The completed local MMseqs2 search returned **75 valid alignments** from eight full-length archived queries against **20,834 target proteins**. All eight self hits span the full sequence at 100% identity. Of the 64 directed panel comparisons, 35 were reported and 29 were not reported under this search at E<=10. These statistics come from the full calibrated target database, not a search of an eight-protein target panel.', '',
        '## Main findings', '',
        '- X5S1 and human/mouse DFFB have no reported match in any of the four directed comparisons under this search at E<=10. This is a bounded search result, not proof of non-homology.',
        '- The expected controls work: human and mouse strongly match each other, X5S1 matches X6E7 strongly, and X5S1 matches its X5M8 alternative product across both full products.',
        '- The panel is not uniformly sequence-dark. Human to RN34 is E=6.236e-6, RN34 to human is E=1.132e-4, and RN34/NI93 match in both directions at approximately 1.65e-13. Human also strongly matches WC24 and RM91.',
        '- X5S1 to RM91 (E=0.8462, 44 columns) and X5S1 to NI93 (E=1.950, 134 columns) are reported but weak. They are not significant sequence bridges.',
        '- No iterative profile search was run; the historical three-iteration X5S1 and NI93 results are not reproduced by this direct-search calculation.', '',
        '## Selected controls and cross-group comparisons', '',
        '| Query → target | Reported E-value | Identity (%) | Alignment columns | Query span / full length | Target span / full length |',
        '|---|---:|---:|---:|---|---|'
    ]
    for key in selected:
        row = pairs[key]
        md.append('| '+LABELS[key[0]]+' → '+LABELS[key[1]]+' | '+row['evalue']+' | '+row['pident']+' | '+row['alnlen']+' | '+row['qstart']+'–'+row['qend']+' / '+row['qlen']+' | '+row['tstart']+'–'+row['tend']+' / '+row['tlen']+' |')
    md += ['', 'Reported identity uses the engine precision: all 75 rows equal floor(nident/alnlen × 1000)/10 percent. Aligned identity counts were independently checked; the raw output is unchanged.', '', '## All directed panel pairs', '', 'Rows are queries; columns are targets. NR means **not reported under this search at E<=10**. It does not mean an inferred E-value above 10. E-values are direction-specific; missingness can also be directional.', '', '| Query / target | Human | Mouse | WC24 | RM91 | RN34 | NI93 | X5S1 | X6E7 |', '|---|---:|---:|---:|---:|---:|---:|---:|---:|']
    for q in PANEL:
        md.append('| '+LABELS[q]+' | '+' | '.join(pairs[q,t]['evalue'] if (q,t) in pairs else 'NR' for t in PANEL)+' |')
    md += ['', '## What was executed and validated', '',
        '- Actual execution: **local MMseqs2**, release '+prep['mmseqs_release']+'; no native plugin executed this search. A native Sequence & Alignment Viewer operation, if performed, is recorded separately.',
        '- UniProt '+prep['uniprot_release']+' ('+prep['uniprot_release_date']+'): UP000192578, 20,832 unique accessions, plus human O76075 and mouse O54788. Distinct accession records were retained without sequence-based deduplication.',
        '- The eight archived query sequences are all exact amino-acid matches to their current records. Search elapsed time was '+format(receipt['elapsed_seconds'],'.6f')+' seconds; exit status 0.',
        '- All 75 alignments independently pass full query/target/header equality, sequence lengths, 1-based coordinate bounds, ungapped alignment-to-sequence mapping, alignment length, identities, mismatches, gap openings, exact CIGAR reconstruction and coverage precision checks. There are 42 distinct targets and no duplicated directed pair rows.',
        '- Reporting threshold: E<=10; sensitivity 7.5; alignment mode 3; max-seqs 1000; format mode 4; verbosity 1. Four threads bound local resources; -a 1 retains backtraces for alignment output.', '',
        '## Historical method versus fresh implementation', '',
        'Recovered S1 §14.6 (PDF page 26) documents the direct search and the separate iterative profiles. Its MMseqs2 build, cb12a2d75a9808ee61721029d064a7bd80af6fec, differs from the installed 18-8cc5c build '+prep['mmseqs_reported_version']+'. The original target snapshot/hash and complete command ledger were not recovered. Identical target counts and matching rounded benchmark results do not establish identical full databases.', '',
        'S1 clips its output-field list after qle. The complete 24-field list in this run was openly specified from installed help; alignment strings, full sequences and headers support residue-level inspection. The explicit four-thread limit and backtrace flag are also fresh operational choices. This is a guided current-database reproduction implementation, not an exact replay of the historical environment.', '',
        'The current X5S1→X6E7, X5S1→RM91, X5S1→NI93 and human→RN34 values agree with the rounded historical benchmarks summarized in the scientific brief. Original raw direction-specific search output remains unavailable for a full byte-level result comparison.', '',
        'Historical claims about iterative profile negatives remain source-reported. These direct results do not demonstrate nuclease activity, ICAD independence or an evolutionary loss event.', '',
        '## Files and reproducibility', '',
        '- [All 64 directed panel pairs]('+str(out/'PANEL_PAIRS.tsv')+')',
        '- [Machine-readable validation and summary]('+str(out/'RESULT_SUMMARY.json')+')',
        '- [Unmodified 75-row MMseqs output]('+str(raw)+')',
        '- [Run receipt with exact command]('+str(receipt_path)+')',
        '- [Frozen input preparation and hashes]('+str(prep_path)+')',
        '- [Offline create-only parser]('+str(Path(__file__).resolve())+')',
        'Recovered S1 §14.6 / PDF page 26 supplied the historical method; copyrighted report is retained outside this review repository.', '',
        'Raw output SHA-256: `'+sha(raw)+'`.', '',
        'Revalidate without writing or rerunning the search:', '', '```bash',
        '/usr/bin/python3 '+str(Path(__file__).resolve())+' --validate-only', '```', ''
    ]
    return out, summary, pair_text, '\n'.join(md)

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--project', type=Path, default=Path(__file__).resolve().parents[1])
    ap.add_argument('--validate-only', action='store_true')
    args = ap.parse_args()
    out, summary, pairs, report = audit(args.project.resolve())
    files = {out/'PANEL_PAIRS.tsv': pairs, out/'SEQUENCE_RESULTS.md': report}
    summary['generated_report'] = {'path':str(out/'SEQUENCE_RESULTS.md'), 'sha256':hashlib.sha256(report.encode()).hexdigest()}
    files[out/'RESULT_SUMMARY.json'] = json.dumps(summary, indent=2)+'\n'
    if not args.validate_only:
        assert not any(path.exists() for path in files), 'Refusing to overwrite an existing output; use --validate-only.'
        for path, content in files.items():
            with path.open('x') as f:
                f.write(content)
        print(json.dumps({'created': {str(path):sha(path) for path in files}}, indent=2))
    print(json.dumps({'validation':'passed', 'rows':summary['data_row_count'], 'self_hits':8, 'panel_pairs':summary['directed_panel_pairs']}, indent=2))

if __name__ == '__main__':
    main()
