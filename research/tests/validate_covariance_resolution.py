"""Independent GLS arithmetic, raw BED LD, provenance and frozen-family checks."""
from pathlib import Path
import hashlib,json,subprocess,re,math
import numpy as np
import pandas as pd
from scipy.linalg import cho_factor,cho_solve
from scipy.stats import norm,false_discovery_control
from validate_release import png_info

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'research/results/covariance_resolution'


def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
    return h.hexdigest()


def main():
    run=json.loads((OUT/'run.json').read_text());selection=json.loads((OUT/'selection.json').read_text())
    for record in (run,selection):
        for n,h in record['inputs'].items():assert sha(ROOT/n)==h,n
        for n,h in record['outputs'].items():assert sha(OUT/n)==h,n
    plan=pd.read_csv(OUT/'frozen_analysis_plan.tsv',sep='\t');picked=plan[plan.selected]
    assert not picked.duplicated(['branch','uniprot','variant_key']).any()
    assert not picked.duplicated(['branch','uniprot','cohort']).any()
    assert picked.groupby('branch').uniprot.nunique().eq(2709).all()
    r=np.loadtxt(OUT/'LD_signed.ld');bim=pd.read_csv(OUT/'LD_variant_order.tsv',sep='\t')
    n_people=len((OUT/'LD_subset.fam').read_text().splitlines());assert n_people==503
    bed=(OUT/'LD_subset.bed').read_bytes();assert bed[:3]==b'\x6c\x1b\x01'
    stride=(n_people+3)//4
    def genotype(i):
        block=np.frombuffer(bed[3+i*stride:3+(i+1)*stride],dtype=np.uint8)
        codes=((block[:,None]>>np.array([0,2,4,6]))&3).ravel()[:n_people]
        return np.array([2.,np.nan,1.,0.])[codes]
    pairs=[]
    for _,g in picked[picked.branch=='unique_variant_givw'].groupby('uniprot'):
        if len(g)>1:pairs.append(tuple(g.LD_index.to_numpy(int)[:2]))
    pairs=pairs[:40];ld_errors=[]
    for i,j in pairs:
        x,y=genotype(i),genotype(j);keep=np.isfinite(x)&np.isfinite(y)
        exact=np.corrcoef(x[keep],y[keep])[0,1];ld_errors.append(abs(exact-r[i,j]))
    assert max(ld_errors)<5e-6  # PLINK's text matrix uses limited significant digits.
    data=pd.read_csv(OUT/'all_outcome_results.tsv',sep='\t');audit=pd.read_csv(OUT/'harmonization_and_coverage.tsv',sep='\t')
    assert len(data)==2709*2*5
    lookup=data.set_index(['branch','outcome','scenario','uniprot'])
    errors=[];checked=0
    for key,g in audit[audit.included].groupby(['branch','outcome','scenario','uniprot'],sort=False):
        assert not g.variant_key.duplicated().any() and not g.rsid.duplicated().any()
        x=g.beta_x.to_numpy();y=g.y_aligned.to_numpy();sy=g.sy_used.to_numpy()
        if len(g)==1:rr=np.eye(1)
        else:
            ix=g.LD_index.to_numpy(int);s=g.LD_sign.to_numpy();rr=r[np.ix_(ix,ix)]*np.outer(s,s)
        omega=rr*np.outer(sy,sy);factor=cho_factor(omega,lower=True)
        wx=cho_solve(factor,x);den=x@wx;beta=(wx@y)/den;se=den**-.5
        actual=lookup.loc[key]
        np.testing.assert_allclose([actual.beta,actual.se,actual.p],[beta,se,2*norm.sf(abs(beta/se))],rtol=1e-8,atol=1e-10)
        errors.append(abs(actual.beta-beta));checked+=1
    assert checked==int(data.p.notna().sum())
    comparison=pd.read_csv(OUT/'all_proteins_comparison.tsv',sep='\t')
    assert comparison.loc[comparison.beta.isna(),'external_concordant'].isna().all()
    clinical=pd.read_csv(OUT/'clinical_original_17.tsv',sep='\t')
    for _,d in clinical.groupby('branch'):
        assert len(d)==374
        p=d.p.to_numpy();valid=np.isfinite(p);padded=np.where(valid,p,1.)
        np.testing.assert_allclose(d.loc[valid,'BH_fixed_374'],false_discovery_control(padded)[valid],rtol=1e-8,atol=1e-10)
        assert d.loc[~valid,'BH_fixed_374'].isna().all()
    baseline=json.loads((ROOT/'research/records/method_resolution_20261017/first_run_statistics.json').read_text())
    repeated={n:sha(OUT/n)==h for n,h in baseline.items()};assert all(repeated.values())
    original=json.loads((ROOT/'research/records/final_review_20261017/result_baseline.json').read_text())
    assert all(sha(ROOT/n)==h for n,h in original.items())
    convention=ROOT/'research/results/mvp_convention_resolution'
    convention_run=json.loads((convention/'run.json').read_text())
    for n,h in convention_run['inputs'].items():assert sha(ROOT/n)==h,n
    for n,h in convention_run['outputs'].items():assert sha(convention/n)==h,n
    heterogeneity=OUT/'heterogeneity';hr=json.loads((heterogeneity/'run.json').read_text())
    for n,h in hr['inputs'].items():assert sha(ROOT/n)==h,n
    for n,h in hr['outputs'].items():assert sha(heterogeneity/n)==h,n
    hd=pd.read_csv(heterogeneity/'all_outcome_MRE.tsv',sep='\t')
    valid=hd.se.notna();multi=valid&hd.Q_df.gt(0);single=valid&hd.Q_df.eq(0)
    np.testing.assert_allclose(hd.loc[single,'se_MRE'],hd.loc[single,'se'],rtol=1e-8,atol=1e-10)
    np.testing.assert_allclose((hd.loc[multi,'se_MRE']/hd.loc[multi,'se'])**2,np.maximum(1,hd.loc[multi,'Q']/hd.loc[multi,'Q_df']),rtol=1e-8,atol=1e-10)
    dd=hd[(hd.branch=='unique_variant_givw')&(hd.scenario=='CI_Wald')&(hd.gene=='DDX58')].iloc[0]
    assert np.isclose(dd.p_MRE,.01083406,rtol=1e-6)
    figure=OUT/'figures';fm=json.loads((figure/'manifest.json').read_text())
    for n,h in fm['inputs'].items():assert sha(ROOT/n)==h,n
    for n,h in fm['outputs'].items():assert sha(figure/n)==h,n
    pdf=figure/'METHOD_COMPARISON.pdf';info=subprocess.check_output(['pdfinfo',str(pdf)],text=True)
    assert re.search(r'Pages:\s+1\b',info)
    width,height=map(float,re.search(r'Page size:\s+([\d.]+) x ([\d.]+)',info).groups())
    assert width==math.floor(72*fm['width_in']+1e-8) and height==math.floor(72*fm['height_in']+1e-8)
    fonts=subprocess.check_output(['pdffonts',str(pdf)],text=True).splitlines()[2:]
    assert fonts and all(row.split()[-5]=='yes' for row in fonts)
    png=png_info(figure/'METHOD_COMPARISON.png')
    assert png['width']==round(fm['width_in']*600) and png['height']==round(fm['height_in']*600)
    assert png.get('unit')==1 and abs(png['xppm']*.0254-600)<.1
    report=dict(status='passed',GLS_rows_independently_checked=checked,max_beta_absolute_error=max(errors),
      LD_pairs_from_raw_BED=len(pairs),max_LD_rounding_error=max(ld_errors),reference_individuals=n_people,
      all_frozen_proteins=2709,no_repeated_variant_or_source=True,clinical_family_per_branch=374,
      repeated_statistical_tables=repeated,original_30_tables_unchanged=True,
      missing_MVP_direction_kept_NA=True,method_figure=dict(pages=1,fonts_embedded=True,dpi=600),
      multiplicative_dispersion_checked=True,
      inference_scope='Checks implementation and reproducibility; MVP source convention and reference/common-scale assumptions remain',
      input_run_sha256=sha(OUT/'run.json'),validator_sha256=sha(Path(__file__)))
    p=ROOT/'research/records/method_resolution_20261017/validation.json';p.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))


if __name__=='__main__':main()
