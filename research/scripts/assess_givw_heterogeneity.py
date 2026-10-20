"""Multiplicative residual-dispersion sensitivity; leaves gIVW and delta results intact."""
from pathlib import Path
import json,hashlib
import numpy as np
import pandas as pd
from scipy.stats import norm
from mr_covariance import bh_fixed

ROOT=Path(__file__).resolve().parents[2];BASE=ROOT/'research/results/covariance_resolution';OUT=BASE/'heterogeneity'


def expand(d):
    d=d.copy();d['residual_scale']=np.where(d.Q_df>0,np.maximum(1,d.Q/d.Q_df),1.)
    d.loc[d.se.isna(),'residual_scale']=np.nan
    d['se_MRE']=d.se*np.sqrt(d.residual_scale);d['p_MRE']=2*norm.sf(abs(d.beta/d.se_MRE))
    d['heterogeneity_estimable']=d.Q_df.gt(0)&d.se.notna()
    return d


def main():
    OUT.mkdir(exist_ok=True)
    source=BASE/'all_outcome_results.tsv';clinical_source=BASE/'clinical_original_17.tsv'
    d=expand(pd.read_csv(source,sep='\t'));c=expand(pd.read_csv(clinical_source,sep='\t'))
    c['BH_MRE_fixed_374']=c.groupby('branch',sort=False).p_MRE.transform(lambda p:bh_fixed(p,374))
    comparison=d[d.outcome=='mvp'].copy()
    for name in ('finngen','stuart'):
        y=d[d.outcome==name][['branch','uniprot','beta','p_MRE']].rename(columns={'p_MRE':name+'_p_MRE','beta':name+'_beta'})
        comparison=comparison.merge(y,on=['branch','uniprot'],how='left',validate='many_to_one')
    def below(p,cut):return (p<cut).astype('boolean').mask(p.isna(),pd.NA)
    comparison['MRE_discovery']=below(comparison.p_MRE,.05/2709)
    comparison['MRE_external_nominal']=below(comparison.finngen_p_MRE,.05)|below(comparison.stuart_p_MRE,.05)
    comparison['MRE_followup']=comparison.MRE_discovery&comparison.MRE_external_nominal
    tables={'all_outcome_MRE.tsv':d,'protein_MRE_comparison.tsv':comparison,'clinical_MRE.tsv':c}
    for name,tab in tables.items():tab.to_csv(OUT/name,sep='\t',index=False,na_rep='NA',float_format='%.17g')
    dd=comparison[(comparison.branch=='unique_variant_givw')&(comparison.scenario=='CI_Wald')&(comparison.gene=='DDX58')].iloc[0]
    assert np.isclose(dd.p_MRE,.01083406,rtol=1e-6)
    summary=[dict(branch=b,scenario=s,discovery=int(g.MRE_discovery.sum()),followup=int(g.MRE_followup.sum())) for (b,s),g in comparison.groupby(['branch','scenario'])]
    files=[source,clinical_source,Path(__file__),Path(__file__).with_name('mr_covariance.py')]
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    record=dict(status='multiplicative_residual_dispersion_sensitivity_computed',summary=summary,
       clinical=[dict(branch=b,BH05=int((g.BH_MRE_fixed_374<.05).sum())) for b,g in c.groupby('branch')],
       method='phi=max(1,Q/(k-1)); SE_MRE=SE_first*sqrt(phi); k=1 retains Wald; independent of delta error-propagation grid',
       caveat='Only two tools for DDX58 and IL7R; dispersion is imprecise and does not identify the cause of heterogeneity',
       inputs={str(p.relative_to(ROOT)):sha(p) for p in files},outputs={n:sha(OUT/n) for n in tables})
    (OUT/'run.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record,indent=2))


if __name__=='__main__':main()
