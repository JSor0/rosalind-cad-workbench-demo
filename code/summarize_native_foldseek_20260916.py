#!/usr/bin/env python3
"""Summarize preserved native Foldseek tables; no search/alignment is performed."""
import csv, hashlib, json, pathlib, re, argparse, shutil
ROOT = pathlib.Path(__file__).resolve().parents[1]
ap=argparse.ArgumentParser(description=__doc__)
ap.add_argument('--output',required=True,type=pathlib.Path)
args=ap.parse_args()
RUN=args.output.resolve()
RUN.mkdir(parents=True,exist_ok=False)
for sub in ['forward_afdb50','reciprocal_pdb']:
    (RUN/sub).mkdir()
    for name in ['aln.csv','output.log']:
        shutil.copyfile(ROOT/'analysis/foldseek'/sub/name,RUN/sub/name)

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def summarize(folder, needle):
    path = folder / "aln.csv"
    with path.open() as f:
        reader = csv.DictReader(f)
        rows = list(reader)
        fields = reader.fieldnames
    result = []
    for rank, r in enumerate(rows, 1):
        e = float(r["evalue"])
        e_rank = 1 + sum(float(other["evalue"]) < e for other in rows)
        q, t = r["qaln"], r["taln"]
        assert len(q) == len(t) == int(r["alnlen"]), (rank, "alignment length")
        assert q.replace("-", "") == r["qseq"][int(r["qstart"])-1:int(r["qend"])], (rank, "query span")
        assert t.replace("-", "") == r["tseq"][int(r["tstart"])-1:int(r["tend"])], (rank, "target span")
        paired = sum(a != "-" and b != "-" for a,b in zip(q,t))
        identities = sum(a == b and a != "-" for a,b in zip(q,t))
        assert not any(a == b == "-" for a,b in zip(q,t))
        result.append(dict(provider_rank=rank,evalue_rank=e_rank,paired_columns=paired,
                           identical_columns=identities,**r))
    with (folder/"returned_hits_ranked.tsv").open("x") as f:
        writer=csv.DictWriter(f,fieldnames=["provider_rank","evalue_rank","paired_columns","identical_columns"]+fields,delimiter="\t")
        writer.writeheader(); writer.writerows(result)
    log=(folder/"output.log").read_text()
    versions=sorted(set(re.findall(r"MMseqs Version:\s*(\S+)",log)))
    sizes=sorted(set(map(int,re.findall(r"Target database size:\s*(\d+)",log))))
    significant = {str(t):sum(float(r["evalue"]) <= t for r in rows) for t in (0.001,1,10)}
    selected=[]
    for row in result:
        if needle.lower() in row["target"].lower():
            selected.append({k:v for k,v in row.items() if k not in ("qseq","tseq","qaln","taln")})
    return dict(csv=str(path),csv_sha256=digest(path),
                rows=len(rows),unique_targets=len({r["target"] for r in rows}),
                query_lengths=sorted({int(r["qlen"]) for r in rows}),
                count_by_evalue_threshold=significant,
                selected_hits=selected,
                foldseek_log_version=versions,target_database_counts=sizes,
                database_snapshot_date=None,
                database_date_status="Not disclosed in native schema/settings/log; execution date is not snapshot date.",
                rank_definition="provider_rank preserves native CSV order; evalue_rank is 1 + count of rows with strictly lower E-value (ties share rank).",
                warnings=[line for line in log.splitlines() if re.search(r"warning|error|overflow",line,re.I)],
                native_command=[line for line in log.splitlines() if line.startswith("search ")],
                validation="Every native row passed alignment length, sequence-span and no-double-gap checks.",
                query_limit=1000,query_limit_semantics="maxSeqs candidate/prefilter cap; not a promise of 1000 qualifying returned rows.")

summary={"operation_class":"updated-database search","scientific_status":"fresh reproduction; not exact historical replay",
         "computation_attribution":"Foldseek executed by native Tamarind Bio plugin; this script only validates/summarizes returned files.",
         "script_sha256":digest(pathlib.Path(__file__)),
         "forward":summarize(RUN/"forward_afdb50","A0A1W0X5S1"),
         "reciprocal":summarize(RUN/"reciprocal_pdb","1v0d")}
with (RUN/"RESULT_SUMMARY.json").open("x") as f: json.dump(summary,f,indent=2)
print(json.dumps(summary,indent=2))
