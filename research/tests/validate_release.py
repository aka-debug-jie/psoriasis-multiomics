"""Verify sensitivity calculations, table curation, and exported figure artifacts."""
from pathlib import Path
import csv,hashlib,json,math,re,struct,subprocess
import numpy as np
import pandas as pd
from scipy.stats import norm
ROOT=Path(__file__).resolve().parents[2]
RELEASE=ROOT/'research/release'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def close(a,b):np.testing.assert_allclose(a,b,rtol=1e-8,atol=1e-10,equal_nan=True)
def png_info(path):
    d=path.read_bytes();assert d[:8]==b'\x89PNG\r\n\x1a\n';offset=8;info={}
    while offset<len(d):
        n=struct.unpack_from('>I',d,offset)[0];kind=d[offset+4:offset+8];payload=d[offset+8:offset+8+n]
        if kind==b'IHDR':info['width'],info['height']=struct.unpack_from('>II',payload)
        if kind==b'pHYs':info['xppm'],info['yppm'],info['unit']=struct.unpack('>IIB',payload)
        offset+=12+n
        if kind==b'IEND':break
    return info

def main():
    report={};raw_manifest=json.loads((RELEASE/'data/manifest.json').read_text())
    for item in raw_manifest['files']:
        source=ROOT/'research'/item['source'];target=RELEASE/'data'/item['output']
        assert sha(source)==item['source_sha256']==sha(target)==item['output_sha256']
        d=pd.read_csv(target,sep='\t');assert d.shape==(item['rows'],item['columns'])
        curated=RELEASE/'tables'/item['output']
        with target.open(newline='') as a,curated.open(newline='') as b:
            x=csv.reader(a,delimiter='\t');y=csv.reader(b,delimiter='\t');hx=next(x);hy=next(y);assert len(hy)==len(set(hy))
            for row in x:assert row==next(y)
            assert next(y,None) is None
    report['tables']=dict(count=len(raw_manifest['files']),source_bytes_preserved=True,curated_values_preserved=True,curated_headers_unique=True)
    audit=ROOT/'research/results/shared_outcome_audit';frozen=pd.read_csv(audit/'single_IV_frozen.tsv',sep='\t')
    original=pd.read_csv(ROOT/'research/data/processed/substitute_mr/frozen_instruments.tsv',sep='\t');original=original[original.eligibility=='eligible']
    expected=original.sort_values(['F_statistic','cohort','assay','rsid'],ascending=[False,True,True,True]).drop_duplicates('uniprot')
    assert len(frozen)==2709 and set(zip(frozen.cohort,frozen.uniprot,frozen.rsid))==set(zip(expected.cohort,expected.uniprot,expected.rsid))
    follow=pd.read_csv(audit/'followup_comparison.tsv',sep='\t');assert len(follow)==17
    unavailable=~follow.mvp_single_calculated;assert unavailable.sum()==3
    assert follow.loc[unavailable,['single_IV_discovery','single_IV_followup']].isna().all().all()
    assert int(follow.single_IV_followup.sum())==6
    d=pd.read_csv(audit/'clinical_comparison.tsv',sep='\t');valid=d.single_p.notna();assert len(d)==374 and valid.sum()==308
    p=d.loc[valid,'single_p'].to_numpy();order=np.argsort(p);q=np.minimum.accumulate((p[order]*374/np.arange(1,len(p)+1))[::-1])[::-1].clip(max=1)
    close(q,d.loc[valid,'single_BH_all_requested_pairs'].to_numpy()[order]);assert d.loc[~valid,'single_BH_all_requested_pairs'].isna().all()
    report['single_IV']=dict(frozen=len(frozen),retained=6,evaluable_not_retained=8,not_evaluable=3,clinical_evaluable=int(valid.sum()),clinical_BH05=int((q<.05).sum()),scope='post_hoc_sensitivity_not_replacement_of_primary_method')
    precision=pd.read_csv(ROOT/'research/results/mvp_precision_audit/mvp_precision_audit.tsv',sep='\t');computable=precision.se_min.notna()
    close(precision.loc[computable,'p_ci_feasible_lower'],2*norm.sf(precision.loc[computable,'beta_abs_max']/precision.loc[computable,'se_min']))
    close(precision.loc[computable,'p_ci_feasible_upper'],2*norm.sf(precision.loc[computable,'beta_abs_min']/precision.loc[computable,'se_max']))
    # Explicit relative-tail check prevents cancellation-to-zero from hiding under absolute tolerance.
    expected_p=2*norm.sf(precision.loc[computable,'beta_abs_max']/precision.loc[computable,'se_min']);observed=precision.loc[computable,'p_ci_feasible_lower'].to_numpy()
    tail=expected_p>1e-300;np.testing.assert_allclose(observed[tail]/expected_p[tail],1,rtol=1e-10,atol=1e-10)
    summary=json.loads((ROOT/'research/results/mvp_precision_audit/summary.json').read_text());assert not summary['current_rsids_still_missing']
    report['MVP_CI']=summary['scopes']['current_mvp_table_rsid_matches']
    manifest=json.loads((RELEASE/'figures_manifest.json').read_text());assert len(manifest['figures'])==8
    for rel,expected_hash in manifest['inputs'].items():assert sha(ROOT/rel)==expected_hash,rel
    for name,record in manifest['panel_tables'].items():
        p=RELEASE/record['file'];assert sha(p)==record['sha256'];tab=pd.read_csv(p,sep='\t');assert len(tab)==record['rows'];assert list(tab)==record['columns']
    checks=[]
    for name,record in manifest['figures'].items():
        for rel,expected_hash in record['outputs'].items():assert sha(RELEASE/rel)==expected_hash,rel
        pdf=RELEASE/'figures'/(name+'.pdf');info=subprocess.check_output(['pdfinfo',str(pdf)],text=True)
        assert re.search(r'Pages:\s+1\b',info)
        width,height=map(float,re.search(r'Page size:\s+([\d.]+) x ([\d.]+)',info).groups())
        # grDevices::cairo_pdf creates integer-point MediaBoxes; check the exact
        # documented device quantization, not an arbitrary widened tolerance.
        assert width==math.floor(72*record['width_in']+1e-8) and height==math.floor(72*record['height_in']+1e-8)
        fonts=subprocess.check_output(['pdffonts',str(pdf)],text=True).splitlines()[2:]
        assert fonts and all(row.split()[-5]=='yes' for row in fonts),fonts
        png=png_info(RELEASE/'figures'/(name+'.png'))
        assert png['width']==round(record['width_in']*600) and png['height']==round(record['height_in']*600)
        assert png.get('unit')==1 and abs(png['xppm']*.0254-600)<.1
        checks.append(dict(figure=name,pages=1,nominal_width_inches=record['width_in'],nominal_height_inches=record['height_in'],pdf_width_points=width,pdf_height_points=height,dpi=600,fonts_embedded=True))
    panel=lambda name:pd.read_csv(RELEASE/'figure_data'/(name+'.tsv'),sep='\t')
    assert len(panel('F01b_forest'))==51 and len(panel('F02b_followup_sensitivity'))==17
    assert panel('F03a_cells').cell_id.nunique()==46225 and panel('F03a_cells').donor_id.nunique()==22
    composition=panel('F03c_donor_composition');assert len(composition)==154 and composition.cells.sum()==46225
    close(composition.groupby('donor').fraction.sum().to_numpy(),np.ones(22))
    for name in ['F04b_cell','F04c_donor']:
        d=panel(name);assert len(d)==102 and not d.duplicated(['gene','cell_type']).any();assert d.loc[~d.tested,'display_FC'].isna().all()
    d=panel('F05a_coloc_windows');assert len(d)==170 and not d.duplicated(['gene_id','source_label','window_kb']).any()
    d=panel('F06_clinical');assert len(d)==374 and (d.support=='Single IV unavailable').sum()==66
    assert (d.support=='Single-IV BH < 0.05').sum()==31
    d=panel('F07a_enrichment');assert len(d)==4 and (d.fdr_bh<.05).all()
    report.update(status='passed',figures=checks,scientific_inference='independent_weight_pooling_remains_a_method_risk',validation_scope='arithmetic, data preservation, export geometry, fonts and source hashes; visual review recorded separately')
    report['inputs']={str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),RELEASE/'data/manifest.json',RELEASE/'tables/manifest.json',RELEASE/'figures_manifest.json',audit/'run.json',ROOT/'research/results/mvp_precision_audit/summary.json']}
    (RELEASE/'qa/validation.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))

if __name__=='__main__':main()
