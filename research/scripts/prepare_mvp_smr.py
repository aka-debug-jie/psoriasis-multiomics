"""Convert the complete original MVP GWAS to SMR format without target filtering."""
from pathlib import Path
from collections import Counter
from datetime import datetime, timezone
import csv
import gzip
import hashlib
import json
import math
import os
import re
import subprocess

ROOT = Path(__file__).resolve().parents[2]
HEADER = ['SNP','A1','A2','freq','b','se','p','n']


def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):
            h.update(block)
    return h.hexdigest()


def transform(row):
    if not re.fullmatch(r'rs[0-9]+',row['rsid']):
        return None,'invalid_rsid',False
    a1,a2=row['effect_allele'].upper(),row['other_allele'].upper()
    if a1 not in {'A','C','G','T'} or a2 not in {'A','C','G','T'} or a1==a2:
        return None,'non_biallelic_SNV',False
    try:
        odds,freq,p,n=[float(row[k]) for k in ['odds_ratio','effect_allele_frequency','p_value','n']]
    except (ValueError,TypeError):
        return None,'missing_effect_frequency_p_or_n',False
    if not all(math.isfinite(x) for x in (odds,freq,p,n)) or odds<=0 or not 0<freq<1 or not 0<=p<=1 or n<=0:
        return None,'invalid_effect_frequency_p_or_n',False
    try:
        se=float(row['standard_error'])
    except (ValueError,TypeError):
        se=float('nan')
    origin='reported_SE'
    if not math.isfinite(se) or se<=0:
        try:
            low,high=float(row['ci_lower']),float(row['ci_upper'])
        except (ValueError,TypeError):
            return None,'missing_SE_and_CI',False
        if not math.isfinite(low) or not math.isfinite(high) or not 0<low<high:
            return None,'invalid_SE_and_CI',False
        se=(math.log(high)-math.log(low))/(2*1.959963984540054)
        origin='SE_from_rounded_OR_95CI'
    clamped=p==0
    values=[freq,math.log(odds),se,1e-300 if p==0 else p,n]
    result=[row['rsid'],a1,a2]+[format(x,'.17g') for x in values]
    return result,origin,clamped


def deduplicate(sorted_path, output, conflicts):
    counts=Counter()
    with sorted_path.open() as stream, output.open('w') as out, conflicts.open('w') as audit:
        out.write('\t'.join(HEADER)+'\n')
        audit.write('SNP\tinput_rows\treason\n')
        key=first=None; size=0; conflict=False
        def flush():
            if key is None:return
            if conflict:
                counts['conflicting_ids']+=1
                counts['conflicting_rows_removed']+=size
                audit.write(f'{key}\t{size}\tconflicting_standardized_effects\n')
            else:
                counts['output_rows']+=1
                counts['identical_duplicates_collapsed']+=size-1
                out.write(first)
        for line in stream:
            identifier=line.split('\t',1)[0]
            if identifier!=key:
                flush(); key=identifier; first=line; size=1; conflict=False
            else:
                size+=1; conflict |= line!=first
        flush()
    return counts


def main():
    source=ROOT/'research/data/raw/gwas/GCST90476190.tsv.gz'
    out=ROOT/'research/data/processed/smr_gwas';out.mkdir(parents=True,exist_ok=True)
    tmp=ROOT/'research/cache/tmp/mvp_smr';tmp.mkdir(parents=True,exist_ok=True)
    unsorted=tmp/'mvp.unsorted.tsv';ordered=tmp/'mvp.sorted.tsv'
    counts=Counter();started=datetime.now(timezone.utc).isoformat()
    with gzip.open(source,'rt') as stream,unsorted.open('w') as target:
        for row in csv.DictReader(stream,delimiter='\t'):
            counts['input_rows']+=1
            result,reason,clamped=transform(row)
            counts[reason]+=1
            if result is not None:
                counts['valid_rows_before_deduplication']+=1
                counts['p_zero_clamped']+=int(clamped)
                target.write('\t'.join(result)+'\n')
            if counts['input_rows']%1000000==0:
                print('PROCESSED',counts['input_rows'],'VALID',counts['valid_rows_before_deduplication'],flush=True)
    subprocess.run(['sort','--stable','--parallel=2','-S','512M','-T',str(tmp),'-t','\t','-k1,1',str(unsorted),'-o',str(ordered)],
                   check=True,env=dict(os.environ,LC_ALL='C'))
    final=out/'mvp_full.ma';conflicts=out/'conflicting_rsids.tsv'
    unique=deduplicate(ordered,final,conflicts)
    assert counts['valid_rows_before_deduplication']==unique['output_rows']+unique['identical_duplicates_collapsed']+unique['conflicting_rows_removed']
    counts.update(unique)
    record=dict(status='complete_GWAS_converted_not_SMR_analysis',started_utc=started,finished_utc=datetime.now(timezone.utc).isoformat(),counts=dict(counts),
                source=str(source.relative_to(ROOT)),source_sha256=digest(source),script_sha256=digest(Path(__file__)),
                outputs={p.name:digest(p) for p in [final,conflicts]},
                source_build='GRCh38; no coordinate columns or cross-build position matching in MA',
                selection='All valid biallelic SNVs from complete source; no gene/target/P/MAF threshold filter',
                duplicate_rule='Identical standardized eight-column rows collapsed; entire conflicting rsID groups removed',
                se_approximation='Missing/invalid SE reconstructed from rounded OR 95% CI; count recorded',
                p_zero_rule='Reported zero P represented as 1e-300 for SMR input; count recorded')
    (out/'run.json').write_text(json.dumps(record,indent=2)+'\n')
    # Only this script's temporary conversion products are disposable.
    unsorted.unlink();ordered.unlink()
    print(json.dumps(record,indent=2))


if __name__=='__main__':main()
