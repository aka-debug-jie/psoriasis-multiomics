"""Sensitivity to shared outcome estimates; preserves the paper-style reconstruction."""
from pathlib import Path
from datetime import datetime,timezone
import json
import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control
from audit_mr_selection import digest
ROOT=Path(__file__).resolve().parents[2]

def main():
    out=ROOT/'research/results/shared_outcome_audit';out.mkdir(parents=True,exist_ok=True)
    source=ROOT/'research/data/processed/substitute_mr/frozen_instruments.tsv'
    registry_path=ROOT/'research/data/processed/substitute_mr/protein_registry.tsv'
    selection_path=ROOT/'research/data/processed/substitute_mr/selection.json'
    x=pd.read_csv(source,sep='\t');reg=pd.read_csv(registry_path,sep='\t')
    x=x[x.eligibility=='eligible'].copy();x['instrument_id']=x.cohort+'::'+x.uniprot
    # Fix exactly one exposure index per protein, without viewing any outcome.
    frozen=x.sort_values(['F_statistic','cohort','assay','rsid'],ascending=[False,True,True,True]).drop_duplicates('uniprot')
    frozen.to_csv(out/'single_IV_frozen.tsv',sep='\t',index=False,na_rep='NA')
    settings=json.loads(selection_path.read_text());denom=settings['fixed_testing_denominator'];assert len(frozen)==denom
    inputs=[source,registry_path,selection_path];singles={};reuse={}
    main_dir=ROOT/'research/results/mr_fenland_substitute'
    for outcome in ['mvp','finngen','stuart']:
        path=main_dir/(outcome+'_cohort_wald.tsv');inputs.append(path);d=pd.read_csv(path,sep='\t')
        cols=['uniprot','instrument_id','rsid','cohort','gene','beta_mr','se_mr','p','OR','CI_low','CI_high']
        s=frozen[['uniprot','instrument_id']].merge(d[cols],on=['uniprot','instrument_id'],how='left',validate='one_to_one')
        s['calculated']=s.p.notna();s['status']=np.where(s.calculated,'fixed_single_IV_calculated','frozen_IV_unavailable_no_substitution')
        s.to_csv(out/(outcome+'_single_IV.tsv'),sep='\t',index=False,na_rep='NA');singles[outcome]=s
        rows=[]
        for uid,g in d.groupby('uniprot'):
            repeated=g.groupby('rsid').size();rows.append(dict(uniprot=uid,gene=g.gene.iloc[0],n_cohort_estimates=len(g),n_distinct_outcome_snps=g.rsid.nunique(),
             shared_same_SNP=bool((repeated>1).any()),repeated_rsids=';'.join(repeated[repeated>1].index),
             max_times_same_outcome=int(repeated.max()),cohort_p_min=g.p.min(),cohort_p_max=g.p.max()))
        reuse[outcome]=pd.DataFrame(rows);reuse[outcome].to_csv(out/(outcome+'_reuse.tsv'),sep='\t',index=False)
    main_path=main_dir/'all_proteins_comparison.tsv';inputs.append(main_path);m=pd.read_csv(main_path,sep='\t')
    comparison=m.merge(reuse['mvp'].drop(columns=['gene']),on='uniprot',validate='one_to_one')
    for outcome,s in singles.items():
        s=s[['uniprot','beta_mr','se_mr','p','cohort','calculated']].rename(columns={k:outcome+'_single_'+k for k in ['beta_mr','se_mr','p','cohort','calculated']})
        comparison=comparison.merge(s,on='uniprot',how='left',validate='one_to_one')
    comparison['pooled_followup']=comparison.passes_frozen_bonferroni&comparison.externally_supported
    def below(p,cutoff):return (p<cutoff).astype('boolean').mask(p.isna(),pd.NA)
    comparison['single_IV_discovery']=below(comparison.mvp_single_p,.05/denom)
    comparison['single_IV_external_nominal']=below(comparison.finngen_single_p,.05)|below(comparison.stuart_single_p,.05)
    comparison['single_IV_external_concordant']=(below(comparison.finngen_single_p,.05)&(comparison.mvp_single_beta_mr*comparison.finngen_single_beta_mr>0).astype('boolean').mask(comparison.mvp_single_beta_mr.isna()|comparison.finngen_single_beta_mr.isna(),pd.NA))|(below(comparison.stuart_single_p,.05)&(comparison.mvp_single_beta_mr*comparison.stuart_single_beta_mr>0).astype('boolean').mask(comparison.mvp_single_beta_mr.isna()|comparison.stuart_single_beta_mr.isna(),pd.NA))
    comparison['single_IV_followup']=comparison.single_IV_discovery&comparison.single_IV_external_nominal
    comparison['single_IV_assessment']=np.where(~comparison.mvp_single_calculated,'not_evaluable',np.where(comparison.single_IV_followup.fillna(False),'retained','not_retained_in_observed_data'))
    comparison['inference_status']=np.where(comparison.shared_same_SNP,'shared_outcome_correlation_omitted','cross_index_LD_covariance_not_modelled')
    comparison.to_csv(out/'all_proteins_comparison.tsv',sep='\t',index=False,na_rep='NA')
    chosen=comparison[comparison.pooled_followup];chosen.to_csv(out/'followup_comparison.tsv',sep='\t',index=False,na_rep='NA')
    # Clinical comparison uses the same exposure-only frozen choice. Missing
    # choices are not replaced by whichever instrument happens to be available.
    clpath=ROOT/'research/results/clinical_mr_substitute/cohort_wald.tsv';cmpath=ROOT/'research/results/clinical_mr_substitute/all_meta_results.tsv';inputs.extend([clpath,cmpath])
    d=pd.read_csv(clpath,sep='\t');pooled=pd.read_csv(cmpath,sep='\t')
    ids=set(frozen.instrument_id);s=d[d.instrument_id.isin(ids)].copy();assert not s.duplicated(['uniprot','phenocode']).any()
    c=pooled.merge(s[['uniprot','phenocode','instrument_id','cohort','rsid','beta_mr','se_mr','p']].rename(columns={'p':'single_p','beta_mr':'single_beta','se_mr':'single_se'}),on=['uniprot','phenocode'],how='left',validate='one_to_one')
    valid=c.single_p.notna();c['single_BH_all_requested_pairs']=np.nan
    # Retain the fixed original 17x22 family, including any uncalculable pairs.
    p=c.loc[valid,'single_p'].to_numpy();order=np.argsort(p);q=np.minimum.accumulate((p[order]*len(c)/np.arange(1,len(p)+1))[::-1])[::-1].clip(max=1)
    restored=np.empty(len(p));restored[order]=q;c.loc[valid,'single_BH_all_requested_pairs']=restored
    c['single_direction_concordant']=np.where(valid,c.single_beta*c.beta>0,np.nan)
    c.to_csv(out/'clinical_comparison.tsv',sep='\t',index=False,na_rep='NA')
    record=dict(status='shared_outcome_sensitivity_computed',created_utc=datetime.now(timezone.utc).isoformat(),
     scientific_inference_status='independent_weight_pooled_results_require_covariance_caution',
     original_analysis_preserved=True,sensitivity_rule='one index per protein by maximum source exposure F; cohort/assay/rsid lexical ties; before outcome matching; no missing-IV replacement',
     purpose='Remove cross-source repeated outcome weighting for sensitivity; not a covariance-aware replacement and not proof of causality',
     denominator=denom,shared_same_SNP_proteins=int(reuse['mvp'].shared_same_SNP.sum()),
     single_IV_calculated=int(singles['mvp'].calculated.sum()),single_IV_discovery=int(comparison.single_IV_discovery.sum()),
     pooled_followup=len(chosen),pooled_followup_shared_same_SNP=int(chosen.shared_same_SNP.sum()),
     pooled_followup_retained_single_IV=int(chosen.single_IV_followup.sum()),
     single_IV_followup_total=int(comparison.single_IV_followup.sum()),
     clinical_requested=len(c),clinical_single_calculated=int(valid.sum()),clinical_pooled_BH05=int((c.BH_all_requested_pairs<.05).sum()),clinical_single_BH05=int((c.single_BH_all_requested_pairs<.05).sum()),
     limitations=['Exposure-only deterministic rule defined for post hoc sensitivity, not a preregistered or recovered author method','Retains CI-inferred MVP SE and uncertain Fenland scaling','Different SNPs may be in LD; shared identical SNPs are only the directly demonstrated dependence','Same disease GWAS is not an independent replication across exposure cohorts'],
     inputs={str(p.relative_to(ROOT)):digest(p) for p in inputs},script_sha256=digest(Path(__file__)),outputs={p.name:digest(p) for p in out.glob('*.tsv')})
    (out/'run.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps({k:v for k,v in record.items() if k not in ['inputs','outputs']},indent=2))

if __name__=='__main__':main()
