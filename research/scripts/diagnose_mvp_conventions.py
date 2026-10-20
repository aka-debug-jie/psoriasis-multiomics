"""Test a source-code-motivated 2P signature, without correcting source statistics."""
from pathlib import Path
from decimal import Decimal, InvalidOperation
import csv, json, hashlib

ROOT=Path(__file__).resolve().parents[2]
SOURCE=ROOT/'research/results/mvp_precision_audit/mvp_precision_audit.tsv'
OUT=ROOT/'research/results/mvp_convention_resolution'


def number(s):
    try:
        d=Decimal(s)
        return d if d.is_finite() else None
    except InvalidOperation: return None


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    with SOURCE.open() as f: data=[r for r in csv.DictReader(f,delimiter='\t') if r['currently_used_rsid'].lower()=='true']
    assert len(data)==len({r['original_rsid'] for r in data})
    rows=[]
    for r in data:
        p=number(r['original_p_value'])
        ends=[number(r[k]) for k in ('p_ci_feasible_lower','p_ci_feasible_upper','p_rounding_lower','p_rounding_upper')]
        d={k:r[k] for k in ('original_rsid','original_p_value','status','source_record','source_file','p_ci_feasible_lower','p_ci_feasible_upper')}
        d['two_P_low']=d['two_P_high']='NA'
        if p is None or not 0<p<Decimal('.5'): d['two_P_status']='not_evaluated_P_outside_range'
        elif any(x is None for x in ends): d['two_P_status']='not_evaluated_CI_interval_unresolved'
        else:
            cl,cu,pl,pu=ends
            d['two_P_low']=str(2*pl);d['two_P_high']=str(2*pu)
            d['two_P_status']='intersects' if max(cl,2*pl)<=min(cu,2*pu) else 'does_not_intersect'
        rows.append(d)
    target=OUT/'two_P_diagnostic.tsv'
    with target.open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
    counts={s:{t:sum(r['status']==s and r['two_P_status']==t for r in rows) for t in sorted({r['two_P_status'] for r in rows})} for s in sorted({r['status'] for r in rows})}
    record=dict(status='source_convention_unresolved_after_signature_test',main_unique_rsids=len(rows),counts=counts,
                signature='If SE=abs(beta)/normal_isf(P), two-sided Wald P becomes 2P for 0<P<0.5',
                interpretation='An intersecting rounding interval does not identify the generating source-code branch or CI convention; source statistics unchanged',
                inputs={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (SOURCE,Path(__file__))},
                outputs={target.name:hashlib.sha256(target.read_bytes()).hexdigest()})
    (OUT/'run.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record,indent=2))


if __name__=='__main__':main()
