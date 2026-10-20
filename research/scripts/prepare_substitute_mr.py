"""Freeze instruments before outcome lookup for the authorized Fenland replacement."""
from pathlib import Path
import argparse
import json
import re
import numpy as np
import pandas as pd
from audit_mr_selection import gene_locations, digest
from gene_identity import identity_resolver

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'research/data/processed/substitute_mr'

def freeze():
    OUT.mkdir(parents=True,exist_ok=True)
    base=ROOT/'research/data/processed/exposures/eldjarn_all_cis_unselected.tsv'
    fen=ROOT/'research/data/processed/Fenland_Olink/source_indices.tsv'
    gtf=ROOT/'research/data/metadata/gencode.v19.annotation.gtf.gz'
    x=pd.read_csv(base,sep='\t').rename(columns={'se_reconstructed':'se_exposure'})
    x['source_row']=np.arange(1,len(x)+1)
    x['unit']='source_standardized_beta'
    x['gene_id']=None
    x['sdY_estimated']=1.
    x['raw_beta']=x.beta
    x['raw_se']=x.se_exposure
    f=pd.read_csv(fen,sep='\t')
    # X dosage and PAR handling are not established for sdY.est's diploid formula.
    f['exclusion']=np.select([~f.chr.between(1,22),f.assay_multi_protein,
                             ~f.source_published_threshold_pass],
                            ['sex_chromosome_scale_unresolved','multi_protein_assay','below_published_4_5e11_threshold'],default='eligible')
    f=f.sort_values(['log10p','id','rsid'],ascending=[False,True,True])
    repeated=f.loc[f.exclusion=='eligible','Uniprot_ID'].duplicated()
    f.loc[repeated[repeated].index,'exclusion']='same_protein_secondary_panel'
    f.to_csv(OUT/'fenland_selection.tsv',sep='\t',index=False)
    chosen=f[f.exclusion=='eligible'].copy()
    z=pd.DataFrame(dict(cohort=chosen.cohort,assay=chosen.id,gene=chosen.Assay,
                       uniprot=chosen.Uniprot_ID,variant=chosen.rsid,chr=chosen.chr,pos=chosen.pos,
                       EA=chosen.Allele1,OA=chosen.Allele2,eaf=chosen.Freq1,
                       beta=chosen.beta_estimated_SD,se_exposure=chosen.se_estimated_SD,
                       raw_beta=chosen.Effect,raw_se=chosen.StdErr,
                       neg_log10_p=chosen.log10p,build='GRCh37',source_row=chosen.index+1,
                       unit='per_estimated_SD_not_observed_SD',gene_id=chosen.ensembl_gene_id,
                       sdY_estimated=chosen.sdY_estimated))
    # For Fenland, the registered ENSG provides identity; names may include panel suffixes.
    loc=gene_locations(gtf)
    loc=loc[loc.gene_chr.isin(['chr'+str(i) for i in range(1,23)])].copy()
    loc['gene_id']=loc.gene_id.str.split('.').str[0]
    byid=loc.drop_duplicates('gene_id').set_index('gene_id')
    geneinfo=ROOT/'research/data/metadata/Homo_sapiens.gene_info.gz'
    resolve,_=identity_resolver(geneinfo,loc)
    z['gene']=[byid.loc[g,'gene'] if g in byid.index else name for g,name in zip(z.gene_id,z.gene)]
    x=pd.concat([x,z],ignore_index=True)
    x['rsid']=x.variant.fillna('').map(lambda s: ','.join(sorted(t for t in str(s).split(',') if re.fullmatch('rs[0-9]+',t))))
    x['F_statistic']=(x.beta/x.se_exposure)**2
    x['source_eligible']=x.rsid.ne('') & x.EA.isin(list('ACGT')) & x.OA.isin(list('ACGT')) & x.EA.ne(x.OA) & np.isfinite(x.beta) & x.se_exposure.gt(0) & x.F_statistic.gt(10) & x.chr.astype(str).isin([str(i) for i in range(1,23)])
    eligible=x[x.source_eligible].sort_values(['neg_log10_p','assay','rsid'],ascending=[False,True,True])
    frozen=eligible.drop_duplicates(['cohort','uniprot']).copy()
    # A protein's coding-gene identity must agree across sources. Do not pick the
    # most convenient gene to determine MHC inclusion.
    registry=[]
    for prot,g in frozen.groupby('uniprot'):
        names=set(g.gene.dropna());r=dict(uniprot=prot,genes=';'.join(sorted(names)),gene_id=None)
        identities=[]
        for _,row in g.iterrows():
            identities.append(row.gene_id if pd.notna(row.gene_id) and row.gene_id in byid.index else resolve(row.gene))
        r['resolved_source_identities']=';'.join(str(v) for v in identities)
        if any(v is None for v in identities):r['eligibility']='gene_location_unresolved'
        elif len(set(identities))!=1:r['eligibility']='source_gene_conflict'
        else:
            gid=identities[0];q=byid.loc[gid]
            r.update(gene_id=gid,gene_symbol=q.gene,gene_chr=q.gene_chr,gene_start=int(q.gene_start),gene_end=int(q.gene_end))
            r['eligibility']='MHC_excluded' if q.mhc_overlap else 'eligible'
        registry.append(r)
    registry=pd.DataFrame(registry)
    frozen=frozen.merge(registry[['uniprot','eligibility']],on='uniprot',validate='many_to_one')
    x.to_csv(OUT/'source_rows.tsv',sep='\t',index=False,na_rep='NA')
    frozen.to_csv(OUT/'frozen_instruments.tsv',sep='\t',index=False,na_rep='NA')
    registry.to_csv(OUT/'protein_registry.tsv',sep='\t',index=False,na_rep='NA')
    n=int((registry.eligibility=='eligible').sum())
    record=dict(status='frozen_before_new_outcome_lookup',formal_original_data=False,
                definition='Source significant cis rows, autosomal biallelic rsID, F>10; strongest exposure association per cohort/protein before any outcome lookup; no replacement if absent',
                fenland_threshold=4.5e-11,fenland_independent_proteins=len(chosen),
                fenland_exclusions=f.exclusion.value_counts().to_dict(),
                protein_eligibility=registry.eligibility.value_counts().to_dict(),
                fixed_testing_denominator=n,bonferroni_threshold=.05/n,
                steiger_filter_applied=False,scale='Fenland beta/SE divided by estimated sdY; scale uncertainty not in SE',
                sources={str(p.relative_to(ROOT)):digest(p) for p in [base,fen,gtf,geneinfo,Path(__file__).with_name('gene_identity.py')]},
                script_sha256=digest(Path(__file__)),outputs={p.name:digest(p) for p in OUT.glob('*.tsv')})
    (OUT/'selection.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(record,indent=2),flush=True)

def extract():
    frozen=pd.read_csv(OUT/'frozen_instruments.tsv',sep='\t')
    ids={s for v in frozen.rsid for s in v.split(',')}
    specs=[('mvp','research/data/processed/smr_gwas/mvp_full.ma','SNP'),
           ('finngen','research/data/raw/gwas/finngen_R12_L12_PSORIASIS.gz','rsids'),
           ('stuart','research/data/raw/gwas/34927100-GCST90019016-EFO_0000676.h.tsv.gz','hm_rsid')]
    records=[]
    for name,rel,rscol in specs:
        path=ROOT/rel;parts=[];total=0
        for c in pd.read_csv(path,sep='\t',chunksize=200000,dtype=str,na_values=['#NA']):
            total+=len(c)
            rs=c[rscol].fillna('');keep=rs.isin(ids)
            aliases=rs.str.contains('[;,]',regex=True)
            if aliases.any():keep.loc[aliases]|=rs[aliases].map(lambda v:bool(set(v.replace(';',',').split(','))&ids))
            if keep.any():parts.append(c[keep])
        d=pd.concat(parts,ignore_index=True)
        if name=='mvp':
            d=d.rename(columns={'SNP':'rsid','A1':'EA','A2':'OA','freq':'eaf','b':'beta','p':'p_reported'})
        elif name=='finngen':
            d=d.rename(columns={'rsids':'rsid','alt':'EA','ref':'OA','af_alt':'eaf','sebeta':'se','pval':'p_reported','#chrom':'chr'})
        else:
            d=d.drop(columns=['beta'])
            d=d.rename(columns={'hm_rsid':'rsid','hm_effect_allele':'EA','hm_other_allele':'OA','hm_beta':'beta','standard_error':'se','p_value':'p_reported','hm_chrom':'chr','hm_pos':'pos'})
            d['eaf']=np.nan
        d['rsid']=d.rsid.str.replace(';',',',regex=False).str.split(',');d=d.explode('rsid')
        d=d[d.rsid.isin(ids)]
        cols=['rsid','EA','OA','eaf','beta','se','p_reported']+[k for k in ['chr','pos','n'] if k in d]
        d=d[cols]
        for k in ['beta','se','eaf','p_reported']:d[k]=pd.to_numeric(d[k],errors='coerce')
        d=d.drop_duplicates()
        bad=d[d.rsid.duplicated(keep=False)]
        bad.to_csv(OUT/(name+'_conflicting_rows.tsv'),sep='\t',index=False)
        d=d[~d.rsid.isin(bad.rsid)]
        target=OUT/(name+'_outcome.tsv');d.to_csv(target,sep='\t',index=False,na_rep='NA')
        records.append(dict(outcome=name,source_rows=total,matched_rows=len(d),conflicting_rsids=bad.rsid.nunique(),source=rel,source_sha256=digest(path),output_sha256=digest(target)))
        print(name,total,len(d),flush=True)
    (OUT/'outcome_extraction.json').write_text(json.dumps(dict(requested=len(ids),selection_sha256=digest(OUT/'selection.json'),script_sha256=digest(Path(__file__)),outcomes=records),indent=2)+'\n')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--stage',choices=['freeze','extract'],required=True)
    args=p.parse_args();freeze() if args.stage=='freeze' else extract()
