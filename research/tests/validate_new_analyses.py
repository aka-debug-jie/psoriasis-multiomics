"""Independent arithmetic and provenance checks for the executed reconstruction."""
from pathlib import Path
from functools import lru_cache
import argparse,hashlib,json
import numpy as np
import pandas as pd
from scipy.special import logsumexp,gammaln
from scipy.stats import norm
ROOT=Path(__file__).resolve().parents[2]
VERIFIED={}

@lru_cache(None)
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
    return h.hexdigest()

def close(a,b):np.testing.assert_allclose(a,b,rtol=1e-8,atol=1e-10)

def provenance(folder):
    run=json.loads((folder/'run.json').read_text())
    scripts={'functional_substitute':'run_functional_enrichment.py','drug_network_substitute':'run_drug_network.py','evidence_substitute':'integrate_reconstructed_evidence.py'}
    if folder.name in scripts:
        path=ROOT/'research/scripts'/scripts[folder.name]
        assert sha(path)==run['script_sha256'],str(path)
        VERIFIED[str(path.relative_to(ROOT))]=sha(path)
    for name,expected in run['inputs'].items():
        path=Path(name);path=path if path.is_absolute() else ROOT/path
        assert sha(path)==expected,str(path)
        VERIFIED[str(path.relative_to(ROOT))]=expected
    for name,expected in run['outputs'].items():
        assert sha(folder/name)==expected,str(folder/name)
        VERIFIED[str((folder/name).relative_to(ROOT))]=expected
    return run

def bh(p):
    p=np.asarray(p);o=np.argsort(p);v=np.minimum.accumulate((p[o]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1].clip(max=1)
    result=np.empty(len(p));result[o]=v;return result

def validate_coloc(folder,protein=False):
    provenance(folder);tab=pd.read_csv(folder/'coloc_summary.tsv',sep='\t')
    computed=tab[tab.status.eq('computed' if protein else 'computed_full_available_cis_region')];error=0.;snperror=0.;total=0
    for _,r in computed.iterrows():
        d=pd.read_csv(folder/r.variant_table,sep='\t');assert not d.SNP.duplicated().any();assert len(d)==r.nsnps;total+=len(d)
        b=d.beta.to_numpy() if protein else d.b.to_numpy();se=d.se.to_numpy() if protein else d.SE.to_numpy()
        v=se**2;vy=d.se_gwas.to_numpy()**2;by=d.beta_gwas_aligned.to_numpy();w=(.15*(r.analysis_sdY if protein else r.sdY))**2
        a=.5*(np.log(v/(v+w))+w/(v+w)*b*b/v)
        c=.5*(np.log(vy/(vy+.04))+.04/(vy+.04)*by*by/vy)
        left=np.r_[-np.inf,np.logaddexp.accumulate(c)[:-1]];right=np.r_[np.logaddexp.accumulate(c[::-1])[::-1][1:],-np.inf]
        h=np.array([0,np.log(r.effective_p1)+logsumexp(a),np.log(r.effective_p2)+logsumexp(c),
                    np.log(r.effective_p1*r.effective_p2)+logsumexp(a+np.logaddexp(left,right)),np.log(r.effective_p12)+logsumexp(a+c)])
        posterior=np.exp(h-logsumexp(h));reported=r[[f'PP.H{i}.abf' for i in range(5)]].to_numpy(float)
        close(posterior,reported);error=max(error,float(abs(posterior-reported).max()))
        snp=np.exp(a+c-logsumexp(a+c));close(snp,d['SNP.PP.H4']);snperror=max(snperror,float(abs(snp-d['SNP.PP.H4']).max()))
        assert bool(posterior[4]>.8)==bool(r.strong_coloc)
        if protein:close(d.beta*d.sdY_estimated,d.raw_beta);close(d.se*d.sdY_estimated,d.raw_se)
    return dict(regions=len(tab),computed=len(computed),variant_rows=total,H0_H4_max_abs=error,SNP_H4_max_abs=snperror,strong_regions=int(computed.strong_coloc.sum()))

def half(a,d=6):
    a=np.asarray(a,float);out=np.zeros_like(a);mask=a!=0;out[mask]=.5*10.**(np.floor(np.log10(abs(a[mask])))-d+1);return out

def rounding_overlap(a,lo,hi,d=6):assert ((a-half(a,d)<=hi)&(a+half(a,d)>=lo)).all()

def main(without_clinical=False):
    report={};base=ROOT/'research/results';mr=base/'mr_fenland_substitute';run=provenance(mr)
    for outcome in ['mvp','finngen','stuart']:
        d=pd.read_csv(mr/(outcome+'_cohort_wald.tsv'),sep='\t');assert not d.duplicated(['cohort','uniprot']).any()
        close(d.beta_mr,d.beta_y_aligned/d.beta_exposure);close(d.se_mr,d.se/abs(d.beta_exposure));close(d.p,2*norm.sf(abs(d.beta_mr/d.se_mr)))
        m=pd.read_csv(mr/(outcome+'_meta.tsv'),sep='\t').set_index('uniprot')
        for uid,g in d.groupby('uniprot'):
            weights=1/g.se_mr**2;close(m.loc[uid,'beta'],sum(weights*g.beta_mr)/sum(weights));close(m.loc[uid,'se'],np.sqrt(1/sum(weights)))
    all_m=pd.read_csv(mr/'all_proteins_comparison.tsv',sep='\t');selected=pd.read_csv(mr/'followup_proteins.tsv',sep='\t')
    criterion=(all_m.p<.05/run['frozen_denominator'])&((all_m.p_finngen<.05)|(all_m.p_stuart<.05))
    assert set(all_m.loc[criterion,'uniprot'])==set(selected.uniprot)
    report['MR']=dict(denominator=run['frozen_denominator'],computed=len(all_m),selected=len(selected),opposite_external_directions=selected.loc[~selected.replication_direction_concordant,'gene'].tolist())
    report['mRNA_coloc']=validate_coloc(base/'coloc_mrna_substitute')
    report['protein_coloc']=validate_coloc(base/'coloc_protein_substitute',True)
    smr=[]
    for folder in sorted((base/'smr').iterdir()):
        if not (folder/'run.json').exists():continue
        provenance(folder);d=pd.read_csv(folder/'all_probe_results.tsv',sep='\t')
        bx=d.b_eQTL.to_numpy();by=d.b_GWAS.to_numpy();sx=d.se_eQTL.to_numpy();sy=d.se_GWAS.to_numpy()
        ax0=abs(bx)-half(bx);ax1=abs(bx)+half(bx);ay0=np.maximum(0,abs(by)-half(by));ay1=abs(by)+half(by)
        br=np.stack([(by+i*half(by))/(bx+j*half(bx)) for i in [-1,1] for j in [-1,1]])
        rounding_overlap(d.b_SMR,br.min(0),br.max(0))
        slo=np.sqrt(((sy-half(sy))/ax1)**2+(ay0*(sx-half(sx))/ax1**2)**2)
        shi=np.sqrt(((sy+half(sy))/ax0)**2+(ay1*(sx+half(sx))/ax0**2)**2);rounding_overlap(d.se_SMR,slo,shi)
        zlo=ay0/np.sqrt((sy+half(sy))**2+(ay0*(sx+half(sx))/ax0)**2)
        zhi=ay1/np.sqrt((sy-half(sy))**2+(ay1*(sx-half(sx))/ax1)**2);rounding_overlap(d.p_SMR,2*norm.sf(zhi),2*norm.sf(zlo),7)
        close(bh(d.p_SMR),d.SMR_FDR_all_tested_probes)
        available=d.p_HEIDI.notna();assert d.loc[~available,'nominal_SMR_and_HEIDI_rule'].isna().all()
        rule=(d.p_SMR<.05)&(d.p_HEIDI>.01)
        assert (d.loc[available,'nominal_SMR_and_HEIDI_rule'].astype(bool)==rule[available]).all()
        smr.append(dict(tissue=folder.name,rows=len(d),supported=int(rule.sum()),missing_HEIDI=int((~available).sum())))
    report['SMR']=smr
    folder=base/'multi_instrument_substitute';provenance(folder)
    instruments=pd.read_csv(folder/'harmonized_instruments.tsv',sep='\t');results=pd.read_csv(folder/'method_results.tsv',sep='\t')
    assert not instruments.duplicated(['cohort','uniprot','mode','rsid']).any()
    ncheck=0
    for (cohort,uid,mode),g in instruments.groupby(['cohort','uniprot','mode']):
        assert g.assay.nunique()==1 and (g.neg_log10_p>-np.log10(5e-6)).all()
        r=results[(results.cohort==cohort)&(results.uniprot==uid)&(results['mode']==mode)]
        x=g.beta.to_numpy();y=g.beta_outcome_aligned.to_numpy();w=1/g.se.to_numpy()**2
        ivw=r[r.method.isin(['IVW','single_IV_Wald'])].iloc[0];estimate=sum(w*x*y)/sum(w*x*x)
        q=sum(w*(y-estimate*x)**2);se=np.sqrt((max(1,q/(len(x)-1)) if len(x)>1 else 1)/sum(w*x*x))
        close(estimate,ivw.b);close(se,ivw.se);ncheck+=1
    report['multi_IV']=dict(groups_independent_IVW_checks=ncheck,method_rows=len(results),computed=int(results.status.eq('computed').sum()))
    provenance(base/'steiger_substitute');report['Steiger']=pd.read_csv(base/'steiger_substitute/sensitivity_summary.tsv',sep='\t').to_dict('records')
    folder=base/'functional_substitute';provenance(folder);func={}
    for path in [*folder.glob('GO_*.tsv'),folder/'KEGG.tsv']:
        d=pd.read_csv(path,sep='\t');close(bh(d.p_value),d.fdr_bh)
        errors=[]
        for r in d.itertuples():
            N=r.background_annotated;M=r.term_background;n=r.query_annotated;k=r.overlap
            z=np.arange(max(k,n-(N-M)),min(n,M)+1)
            def lchoose(a,b):return gammaln(a+1)-gammaln(b+1)-gammaln(a-b+1)
            p=np.exp(logsumexp(lchoose(M,z)+lchoose(N-M,n-z)-lchoose(N,n)))
            close(p,r.p_value);errors.append(abs(p-r.p_value))
        func[path.name]=dict(terms=len(d),max_absolute_error=max(errors))
    report['enrichment']=func;provenance(base/'drug_network_substitute')
    provenance(base/'evidence_substitute')
    evidence=pd.read_csv(base/'evidence_substitute/followup_evidence.tsv',sep='\t')
    assert set(evidence.uniprot)==set(selected.uniprot)
    assert evidence[['UKB_PGS_psoriasis','UKB_PGS_PsA']].isna().all().all()
    assert [c for c in evidence if 'score' in c.lower()]==['author_tier_score']
    assert not evidence.loc[evidence.gene=='TGFBI','external_direction_agreement'].any()
    report['evidence']=dict(rows=len(evidence),missing_UKB_retained=True,computed_total_score=False)
    if not without_clinical:
        folder=base/'clinical_mr_substitute';run=provenance(folder);d=pd.read_csv(folder/'cohort_wald.tsv',sep='\t');m=pd.read_csv(folder/'all_meta_results.tsv',sep='\t').set_index(['phenocode','uniprot'])
        close(d.beta_mr,d.beta_aligned/d.beta_exposure);close(d.se_mr,d.sebeta/abs(d.beta_exposure))
        for key,g in d.groupby(['phenocode','uniprot']):
            w=1/g.se_mr**2;close(m.loc[key,'beta'],sum(w*g.beta_mr)/sum(w));close(m.loc[key,'se'],np.sqrt(1/sum(w)))
        report['clinical']=dict(outcomes=run['outcomes'],computed=run['computed_pairs'])
    report['status']='passed';report['clinical_included']=not without_clinical;report['script_sha256']=sha(Path(__file__))
    report['verified_file_hashes']=VERIFIED
    target=ROOT/'research/records/new_analyses_validation.json';target.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--without-clinical',action='store_true');main(p.parse_args().without_clinical)
