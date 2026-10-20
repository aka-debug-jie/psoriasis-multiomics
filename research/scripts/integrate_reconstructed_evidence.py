"""Integrate observed evidence without filling missing score components."""
from pathlib import Path
import json
import pandas as pd
import numpy as np
from audit_mr_selection import digest,gene_locations
from gene_identity import identity_resolver
ROOT=Path(__file__).resolve().parents[2]

def main():
    out=ROOT/'research/results/evidence_substitute';out.mkdir(parents=True,exist_ok=True);inputs=[]
    def table(rel):
        p=ROOT/rel;inputs.append(p);return pd.read_csv(p,sep='\t')
    genes=table('research/data/processed/coloc_substitute_registry/gene_registry.tsv')
    mr=table('research/results/mr_fenland_substitute/all_proteins_comparison.tsv')
    follow=table('research/results/mr_fenland_substitute/followup_proteins.tsv')
    authors=table('research/references/author_main_table1.tsv')
    bulk=table('research/results/bulk_ncbi/all_gene_differential_expression.tsv')
    sc=table('research/results/sc_de/all_gene_psoriasis_vs_normal.tsv')
    pseudo=table('research/results/sc_pseudobulk/donor_level_differential_expression.tsv')
    mrna=table('research/results/coloc_mrna_substitute/coloc_summary.tsv')
    protein=table('research/results/coloc_protein_substitute/coloc_summary.tsv')
    direct=table('research/results/drug_network_substitute/candidate_direct_targets.tsv')
    component=table('research/results/drug_network_substitute/candidate_complex_family_components.tsv')
    network=table('research/results/drug_network_substitute/candidate_network_drug_associations.tsv')
    clinical_path='research/results/clinical_mr_substitute/all_meta_results.tsv'
    clinical=table(clinical_path) if (ROOT/clinical_path).exists() else None
    gtf=ROOT/'research/data/metadata/gencode.v19.annotation.gtf.gz';info=ROOT/'research/data/metadata/Homo_sapiens.gene_info.gz';inputs.extend([gtf,info])
    loc=gene_locations(gtf);loc=loc[loc.gene_chr.isin(['chr'+str(i) for i in range(1,23)])].copy();loc.gene_id=loc.gene_id.str.split('.').str[0]
    resolve,entrez=identity_resolver(info,loc)
    authors['gene_id']=[resolve({'IL23':'IL23A','IL12':'IL12A'}.get(t,t)) for t in authors.target]
    authors=authors.rename(columns={c:'author_'+c for c in authors if c!='gene_id' and not c.startswith('author_')})
    smr={p.parent.name:table(str(p.relative_to(ROOT))) for p in sorted((ROOT/'research/results/smr').glob('*/all_probe_results.tsv'))}
    rows=[];bridges=[];contrast_exports={'cell_Wilcoxon':[],'donor_pseudobulk':[]}
    for g in genes.itertuples():
        r=dict(gene_id=g.gene_id,gene=g.symbol,our_followup=bool(g.available_followup),author_reference=bool(g.author_reference))
        m=mr[mr.gene_id==g.gene_id]
        # Do not collapse multiple assay-protein identities to a synthetic effect.
        r['MR_protein_rows']=len(m)
        if len(m)==1:
            q=m.iloc[0];r.update(uniprot=q.uniprot,MR_beta=q.beta,MR_se=q.se,MR_p=q.p,MR_cohorts=q.cohorts,
                               FinnGen_beta=q.beta_finngen,FinnGen_p=q.p_finngen,Stuart_beta=q.beta_stuart,Stuart_p=q.p_stuart,
                               external_nominal_significance=bool(q.externally_supported),external_direction_agreement=bool(q.replication_direction_concordant))
        else:r.update(uniprot=None,MR_beta=None,MR_p=None)
        eid=entrez(g.symbol);b=bulk[bulk.GeneID.astype(str)==str(eid)]
        r['bulk_tested']=len(b)==1
        if len(b)==1:q=b.iloc[0];r.update(bulk_log2FC=q.logFC,bulk_FDR=q['adj.P.Val'],bulk_significant=bool(q.significant))
        for name,data,signif,fold in [('cell_Wilcoxon',sc,'significant','log2FC'),('donor_pseudobulk',pseudo,'significant','logFC')]:
            layer_ids=set(data.gene_id);layer_id=g.gene_id if g.gene_id in layer_ids else None
            reason='exact_Ensembl_ID' if layer_id else 'not_present_in_tested_gene_set'
            if layer_id is None:
                matches=set(data.loc[data.symbol.eq(g.symbol),'gene_id'])
                if len(matches)==1:layer_id=next(iter(matches));reason='unique_gene_symbol_version_bridge'
            bridges.append(dict(canonical_gene_id=g.gene_id,symbol=g.symbol,layer=name,source_gene_id=layer_id,mapping=reason))
            d=data[data.gene_id==layer_id].copy();tested=d[d[signif].notna()]
            if g.available_followup:
                d['source_gene_id']=d.gene_id;d['gene_id']=g.gene_id;d['canonical_symbol']=g.symbol;contrast_exports[name].append(d)
            r[name+'_types_tested']=len(tested)
            pos=tested[tested[signif].astype(bool)]
            r[name+'_significant_types']=';'.join(pos.cell_type)
            r[name+'_directions']=';'.join(pos.cell_type+':'+np.where(pos[fold]>0,'up','down'))
            r[name+'_any_support']=bool(len(pos)) if len(tested) else None
        for tissue,data in smr.items():
            d=data[data.gene_id==g.gene_id];r['SMR_'+tissue+'_probes']=len(d)
            if len(d)==1:
                q=d.iloc[0];r['SMR_'+tissue+'_P']=q.p_SMR;r['HEIDI_'+tissue+'_P']=q.p_HEIDI
                r['SMR_'+tissue+'_support']=q.nominal_SMR_and_HEIDI_rule
        for name,data in [('mRNA_coloc',mrna),('protein_coloc',protein)]:
            d=data[(data.gene_id==g.gene_id)&data['PP.H4.abf'].notna()]
            r[name+'_computed_regions']=len(d);r[name+'_max_H4']=d['PP.H4.abf'].max() if len(d) else None
            r[name+'_H4_gt_08_regions']=int((d['PP.H4.abf']>.8).sum()) if len(d) else None
        uid=r['uniprot']
        # Drug analyses were executed on the 17 selected proteins, not all references.
        if r['our_followup']:
            r['drug_single_protein_rows']=int((direct.uniprot==uid).sum())
            r['drug_complex_family_component_rows']=int((component.uniprot==uid).sum())
            r['drug_network_association_rows']=int((network.candidate_uniprot==uid).sum())
            if clinical is not None:
                c=clinical[(clinical.uniprot==uid)&clinical.p.notna()]
                r['clinical_outcomes_computed']=len(c);r['clinical_nominal_P05']=int((c.p<.05).sum());r['clinical_global_BH05']=int((c.BH_all_requested_pairs<.05).sum())
        r['UKB_PGS_psoriasis']=None;r['UKB_PGS_PsA']=None;r['UKB_status']='omitted_by_user_not_replaced_by_summary_GWAS'
        rows.append(r)
    full=pd.DataFrame(rows).merge(authors,on='gene_id',how='left',validate='one_to_one')
    full.to_csv(out/'all_registered_evidence.tsv',sep='\t',index=False,na_rep='NA')
    selected=full[full.our_followup];selected.to_csv(out/'followup_evidence.tsv',sep='\t',index=False,na_rep='NA')
    agreement=full[['gene','gene_id','our_followup','author_reference','MR_p','external_nominal_significance','external_direction_agreement','author_tier_score']]
    agreement.to_csv(out/'author_candidate_comparison.tsv',sep='\t',index=False,na_rep='NA')
    pd.DataFrame(bridges).to_csv(out/'expression_gene_identity.tsv',sep='\t',index=False,na_rep='NA')
    for name,parts in contrast_exports.items():pd.concat(parts,ignore_index=True).to_csv(out/(name+'_followup_contrasts.tsv'),sep='\t',index=False,na_rep='NA')
    record=dict(status='observed_evidence_integrated',genes=len(full),our_followup=len(selected),author_reference=int(full.author_reference.sum()),
                total_reconstructed_score_computed=False,author_score_role='reported reference only, not filled from our partial evidence',
                missing='NA retained; zero only for an actually queried/tested count; PGS/Cox slots remain NA',
                clinical_included=clinical is not None,external_direction_conflicts=selected.loc[selected.external_direction_agreement==False,'gene'].tolist(),
                inputs={str(p.relative_to(ROOT)):digest(p) for p in inputs},script_sha256=digest(Path(__file__)),outputs={p.name:digest(p) for p in out.glob('*.tsv')})
    (out/'run.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps({k:v for k,v in record.items() if k not in ['inputs','outputs']},indent=2))

if __name__=='__main__':main()
