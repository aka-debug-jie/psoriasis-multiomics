#!/usr/bin/env python3
"""Integrate the fixed corrected analysis, preserving original results and test families."""
from pathlib import Path
from functools import lru_cache
from collections import defaultdict
import csv,gzip,json,hashlib,argparse
import numpy as np
import pandas as pd
from mr_covariance import bh_fixed
from run_functional_enrichment import read_ontology,enrich

ROOT=Path(__file__).resolve().parents[2]
CONFIG=ROOT/'research/config/primary_analysis.json'


def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for block in iter(lambda:f.read(4*1024*1024),b''):h.update(block)
    return h.hexdigest()


def below(p,cut):return (p<cut).astype('boolean').mask(p.isna(),pd.NA)


def main(output=None):
    cfg=json.loads(CONFIG.read_text());out=ROOT/(output or cfg['output']);data=out/'data';data.mkdir(parents=True,exist_ok=True)
    inputs=[CONFIG,Path(__file__),Path(__file__).with_name('mr_covariance.py'),Path(__file__).with_name('run_functional_enrichment.py')]
    def table(rel):
        p=ROOT/rel;inputs.append(p);return pd.read_csv(p,sep='\t')
    def save(d,name):
        p=data/(name+'.tsv');p.parent.mkdir(parents=True,exist_ok=True);d.to_csv(p,sep='\t',index=False,na_rep='NA',float_format='%.17g')
    allmr=table('research/results/covariance_resolution/all_proteins_comparison.tsv')
    branch=cfg['branch'];scenario=cfg['MVP_scenario'];den=cfg['testing_denominator']
    primary=allmr[(allmr.branch==branch)&(allmr.scenario==scenario)].copy();assert len(primary)==den
    primary=primary.rename(columns={'passes_frozen_bonferroni':'legacy_independent_weight_passes_bonferroni',
        'externally_supported':'legacy_independent_weight_external_nominal','original_p':'legacy_independent_weight_p',
        'original_followup':'legacy_independent_weight_followup'})
    primary['primary_selected']=below(primary.p,.05/den)&(below(primary.finngen_p,.05)|below(primary.stuart_p,.05))
    primary['background_primary']=primary.p.notna()&(primary.finngen_p.notna()|primary.stuart_p.notna())
    primary['background_discovery']=primary.p.notna()
    registry=table('research/data/processed/substitute_mr/protein_registry.tsv')
    primary=primary.merge(registry[['uniprot','gene_id','gene_symbol']],on='uniprot',how='left',validate='one_to_one')
    selected=primary[primary.primary_selected.fillna(False)].copy().sort_values(['p','gene'])
    assert selected.gene_id.notna().all() and not selected.gene_id.duplicated().any()
    uids=set(selected.uniprot);gids=set(selected.gene_id)
    save(primary,'primary_all_proteins');save(selected,'primary_candidates')

    # Whole-study biological tests are retained without recomputing their FDR on seven genes.
    old=table('research/results/evidence_substitute/all_registered_evidence.tsv')
    keep=['gene_id']+[c for c in old if c.startswith(('bulk_','cell_Wilcoxon_','donor_pseudobulk_','SMR_','HEIDI_','mRNA_coloc_','protein_coloc_','UKB_','author_'))]
    evidence=selected.merge(old[keep],on='gene_id',how='left',validate='one_to_one')
    assert evidence.bulk_tested.notna().all()
    biological={}
    for name,path in {
        'expression_cell':'research/results/evidence_substitute/cell_Wilcoxon_followup_contrasts.tsv',
        'expression_pseudobulk':'research/results/evidence_substitute/donor_pseudobulk_followup_contrasts.tsv',
        'expression_gene_identity':'research/results/evidence_substitute/expression_gene_identity.tsv',
        'coloc_mRNA':'research/results/coloc_mrna_substitute/coloc_summary.tsv',
        'coloc_protein':'research/results/coloc_protein_substitute/coloc_summary.tsv'}.items():
        d=table(path);idcol='canonical_gene_id' if 'canonical_gene_id' in d else 'gene_id';d=d[d[idcol].isin(gids)].copy()
        d=d.drop(columns=['available_followup','author_reference'],errors='ignore');biological[name]=d;save(d,name)
    smr=[]
    for p in sorted((ROOT/'research/results/smr').glob('*/all_probe_results.tsv')):
        d=table(str(p.relative_to(ROOT)));d=d[d.gene_id.isin(gids)].copy();d['candidate_scope']='corrected_primary';smr.append(d)
    smr=pd.concat(smr,ignore_index=True);save(smr,'SMR_HEIDI')
    save(evidence[['gene_id','gene']+[c for c in evidence if c.startswith('bulk_')]],'bulk_primary')

    # Primary candidate-set clinical analysis retains the original correction as the main marker.
    clinical=table('research/results/covariance_resolution/clinical_original_17.tsv')
    clinical=clinical[(clinical.branch==branch)&clinical.uniprot.isin(uids)].copy()
    panel=table('research/data/processed/clinical_finngen/phenotype_manifest.tsv')
    assert len(clinical)==len(selected)*len(panel) and not clinical.duplicated(['uniprot','outcome']).any()
    clinical=clinical.rename(columns={'BH_fixed_374':'BH_original_family_374','outcome':'phenocode'})
    clinical['BH_candidate_family_exploratory']=bh_fixed(clinical.p,len(clinical))
    clinical['candidate_family_size']=len(clinical)
    cmre=table('research/results/covariance_resolution/heterogeneity/clinical_MRE.tsv')
    cmre=cmre[(cmre.branch==branch)&cmre.uniprot.isin(uids)][['uniprot','outcome','se_MRE','p_MRE','BH_MRE_fixed_374']].rename(columns={'outcome':'phenocode','BH_MRE_fixed_374':'MRE_BH_original_family_374'})
    clinical=clinical.merge(cmre,on=['uniprot','phenocode'],validate='one_to_one')
    clinical['MRE_BH_candidate_family_exploratory']=bh_fixed(clinical.p_MRE,len(clinical))
    clinical=clinical.merge(panel,on='phenocode',validate='many_to_one',suffixes=('','_phenotype'))
    save(clinical,'clinical_primary');save(panel,'clinical_endpoints')
    summary=clinical.groupby('uniprot').agg(clinical_estimable=('p','count'),clinical_BH374=('BH_original_family_374',lambda s:int((s<.05).sum())),
       clinical_BH_candidate_exploratory=('BH_candidate_family_exploratory',lambda s:int((s<.05).sum()))).reset_index()
    evidence=evidence.merge(summary,on='uniprot',validate='one_to_one')

    # Recompute candidate flags and links on the exact archived graph, not an assumed new API response.
    relations=table('research/results/drug_network_substitute/all_drug_target_relations.tsv')
    relations['candidate']=relations.uniprot.isin(uids)
    mapping=table('research/results/drug_network_substitute/STRING_mapping.tsv')
    network=table('research/results/drug_network_substitute/STRING_network.tsv')
    request=ROOT/'research/data/raw/functional/string_substitute/network.request.json';rawnet=request.with_name('network.json')
    inputs.extend([request,rawnet]);cached=json.loads(request.read_text());assert sha(rawnet)==cached['sha256']
    targets=set(relations.loc[relations.organism=='Homo sapiens','uniprot']);query=uids|targets
    mapping=mapping[mapping.queryItem.isin(query)].copy();ids=set(mapping.stringId)
    net=network[network.stringId_A.isin(ids)&network.stringId_B.isin(ids)&(network.score>.4)].copy()
    reverse={key:set(g.queryItem) for key,g in mapping.groupby('stringId')};links=[]
    for edge in net.itertuples():
        for a,b in ((edge.stringId_A,edge.stringId_B),(edge.stringId_B,edge.stringId_A)):
            for candidate in reverse[a]&uids:
                for target in reverse[b]&targets:
                    for drug in relations[(relations.uniprot==target)&(relations.organism=='Homo sapiens')].itertuples():
                        links.append(dict(candidate_uniprot=candidate,target_uniprot=target,molecule_chembl_id=drug.molecule_chembl_id,
                          target_chembl_id=drug.target_chembl_id,drug_target_relation=drug.relation,STRING_score=edge.score,
                          interpretation='functional_association_not_direct_drug_binding'))
    links=pd.DataFrame(links).drop_duplicates()
    direct=relations[relations.candidate&(relations.relation=='single_protein_target')]
    components=relations[relations.candidate&(relations.relation=='complex_or_family_component')]
    for name,d in [('all_drug_relations',relations),('drug_direct',direct),('drug_components',components),('drug_network_links',links),('STRING_mapping',mapping),('STRING_induced_network',net)]:save(d,'drugs/'+name)
    missing=sorted(query-set(mapping.queryItem));save(pd.DataFrame({'uniprot':missing,'status':'unmapped_in_archived_STRING_query'}),'drugs/STRING_unmapped')
    for name,d,col in [('direct',direct,'uniprot'),('components',components,'uniprot'),('network',links,'candidate_uniprot')]:
        evidence['drug_'+name+'_rows']=evidence.uniprot.map(d.groupby(col).size()).fillna(0).astype(int)

    # Reuse verified identity mapping, then rebuild annotation sets and every hypergeometric family.
    identity=table('research/results/functional_substitute/background_mapping.tsv')
    identity=primary[['uniprot','background_primary','background_discovery']].merge(identity,on='uniprot',how='left',validate='one_to_one')
    identity['GeneID']=identity.GeneID.astype('Int64').astype('string')
    save(identity,'enrichment_background_mapping')
    eligible=identity[(identity.background_status=='eligible')&identity.GeneID.notna()]
    uid_gene=eligible.set_index('uniprot').GeneID.to_dict();assert uids.issubset(uid_gene)
    query_genes={uid_gene[u] for u in uids}
    backgrounds={name:set(eligible.loc[eligible[col],'GeneID']) for name,col in [('primary','background_primary'),('discovery_sensitivity','background_discovery')]}
    raw=ROOT/'research/data/raw/functional';obo=raw/'go-basic.obo';gaf=raw/'goa_human.gaf.gz'
    inputs.extend([obo,gaf,raw/'kegg_hsa_pathway_gene.tsv',raw/'kegg_hsa_pathway_names.tsv'])
    terms,alt,ontology_header=read_ontology(obo)
    @lru_cache(None)
    def ancestors(term):
        if term not in terms:return frozenset()
        return frozenset({term}.union(*(ancestors(p) for p in terms[term]['parents'])))
    go=defaultdict(set);unresolved=set();gaf_header=[]
    with gzip.open(gaf,'rt') as f:
        for line in f:
            if line.startswith('!'):gaf_header.append(line.strip());continue
            r=line.rstrip('\n').split('\t')
            if r[0]!='UniProtKB' or r[1] not in uid_gene or r[12].split('|')[0]!='taxon:9606' or 'NOT' in r[3].split('|') or r[6]=='ND':continue
            term=alt.get(r[4],r[4])
            if term not in terms:unresolved.add(term);continue
            for parent in ancestors(term):go[parent].add(uid_gene[r[1]])
    kegg=defaultdict(set)
    for line in (raw/'kegg_hsa_pathway_gene.tsv').read_text().splitlines():
        gene,pathway=line.split('\t');kegg[pathway.removeprefix('path:')].add(gene.removeprefix('hsa:'))
    names={k.removeprefix('path:'):v for k,v in (line.split('\t',1) for line in (raw/'kegg_hsa_pathway_names.tsv').read_text().splitlines())}
    libraries={ns:({k:v for k,v in go.items() if terms[k]['namespace']==ns},{k:v['name'] for k,v in terms.items()}) for ns in ['biological_process','molecular_function','cellular_component']}
    libraries['KEGG']=(kegg,names);enrichment={}
    for scope,background in backgrounds.items():
        enrichment[scope]={}
        for library,(sets,labels) in libraries.items():
            d,s=enrich(sets,labels,background,query_genes,*cfg['enrichment_term_size'])
            s['significant_terms']=int((d.fdr_bh<.05).sum());enrichment[scope][library]=s
            save(d,'enrichment/'+scope+'/'+library)

    # Existing sensitivity estimates are aligned to the final primary list, without voting into a score.
    sensitivity=allmr[allmr.uniprot.isin(uids)].copy();sensitivity['model']='fixed_effect'
    mre=table('research/results/covariance_resolution/heterogeneity/protein_MRE_comparison.tsv');mre=mre[mre.uniprot.isin(uids)].copy()
    mre['model']='multiplicative_random_effect';mre['p']=mre.p_MRE;mre['se']=mre.se_MRE;mre['followup']=mre.MRE_followup
    keys=['uniprot','gene','branch','scenario','model','beta','se','p','followup','observed_variants','Q','Q_p']
    sensitivity=pd.concat([sensitivity[keys],mre[keys]],ignore_index=True);save(sensitivity,'primary_sensitivity')
    delta=table('research/results/covariance_resolution/delta_error_correlation_grid.tsv');save(delta[delta.uniprot.isin(uids)],'primary_delta_grid')
    used=table('research/results/covariance_resolution/harmonization_and_coverage.tsv')
    used=used[(used.branch==branch)&(used.outcome=='mvp')&(used.scenario==scenario)&used.included&used.uniprot.isin(uids)]
    chosen_pairs=set(zip(used.instrument_id,used.rsid))
    steiger=table('research/results/steiger_substitute/instrument_direction_grid.tsv')
    steiger=steiger[[pair in chosen_pairs for pair in zip(steiger.instrument_id,steiger.rsid)]].copy()
    assert set(zip(steiger.instrument_id,steiger.rsid))==chosen_pairs
    save(steiger,'steiger_selected_instruments')
    multi=table('research/results/multi_instrument_substitute/method_results.tsv');multi=multi[multi.uniprot.isin(uids)].copy()
    multi['interpretation']='source_specific_sensitivity_not_cross_source_pooled';save(multi,'multi_instrument_methods')
    for name in ('coverage','harmonized_instruments'):
        z=table('research/results/multi_instrument_substitute/'+name+'.tsv');save(z[z.uniprot.isin(uids)],'multi_instrument_'+name)
    for label,b,s,m in [('MRE',branch,scenario,'multiplicative_random_effect'),('strict_LD','ld_pruned_givw',scenario,'fixed_effect'),('P_equivalent',branch,'P_equivalent','fixed_effect')]:
        z=sensitivity[(sensitivity.branch==b)&(sensitivity.scenario==s)&(sensitivity.model==m)].set_index('uniprot')
        evidence[label+'_p']=evidence.uniprot.map(z.p);evidence[label+'_followup']=evidence.uniprot.map(z.followup)

    # A complete role-reversed check uses native FinnGen SE and Stuart, never MVP effects or P.
    estimates=table('research/results/covariance_resolution/all_outcome_results.tsv')
    forest=estimates[(estimates.branch==branch)&estimates.uniprot.isin(uids)&(((estimates.outcome=='mvp')&(estimates.scenario==scenario))|((estimates.outcome!='mvp')&(estimates.scenario=='reported_SE')))].copy()
    save(forest,'genetic_effects')
    mrestimates=table('research/results/covariance_resolution/heterogeneity/all_outcome_MRE.tsv')
    controls=[]
    for model,source,pcol,secol in [('fixed_effect',estimates,'p','se'),('multiplicative_random_effect',mrestimates,'p_MRE','se_MRE')]:
        for b in source.branch.unique():
            f=source[(source.branch==b)&(source.outcome=='finngen')].copy();s=source[(source.branch==b)&(source.outcome=='stuart')]
            z=f[['uniprot','gene','beta',pcol,secol,'observed_variants','variants','source_cohorts','coverage']].rename(columns={pcol:'FinnGen_p',secol:'FinnGen_se','beta':'FinnGen_beta','variants':'FinnGen_variants'})
            z=z.merge(s[['uniprot','beta',pcol,'variants']].rename(columns={'beta':'Stuart_beta',pcol:'Stuart_p','variants':'Stuart_variants'}),on='uniprot',how='left',validate='one_to_one')
            z['FinnGen_discovery']=below(z.FinnGen_p,.05/den)
            direction=(z.FinnGen_beta*z.Stuart_beta>0).astype('boolean').mask(z.FinnGen_beta.isna()|z.Stuart_beta.isna(),pd.NA)
            z['Stuart_nominal_concordant']=below(z.Stuart_p,.05)&direction;z['control_followup']=z.FinnGen_discovery&z.Stuart_nominal_concordant
            z['branch']=b;z['model']=model;z['in_MVP_primary_set']=z.uniprot.isin(uids);controls.append(z)
    control=pd.concat(controls,ignore_index=True);save(control,'MVP_free_control_all')
    save(control[control.FinnGen_discovery.fillna(False)],'MVP_free_control_discoveries')
    baseline=control[(control.branch==branch)&(control.model=='fixed_effect')].set_index('uniprot')
    evidence['MVP_free_FinnGen_discovery']=evidence.uniprot.map(baseline.FinnGen_discovery)
    evidence['MVP_free_Stuart_concordant']=evidence.uniprot.map(baseline.Stuart_nominal_concordant)
    evidence['MVP_source_convention']='unresolved_CI_based_primary_with_sensitivities'
    evidence['reconstructed_total_score']=np.nan;save(evidence,'integrated_evidence')
    authors=[c for c in evidence if c.startswith('author_')];save(evidence[['gene','gene_id','uniprot']+authors],'author_reported_comparison')
    author_all=old[old.author_reference][['gene','gene_id']+authors].copy()
    author_all['in_corrected_primary_set']=author_all.gene_id.isin(gids);save(author_all,'author_all_candidates')
    missingness=[];matrix=[]
    def cell(gene,domain,component,computed,planned,support,value=np.nan,mvp=False):
        status='not_evaluable' if computed==0 else ('supported' if support else ('partial_no_observed_support' if computed<planned else 'not_supported_in_tested_data'))
        matrix.append(dict(gene=gene,domain=domain,component=component,computed_units=computed,planned_units=planned,
          support=(support if computed==planned or (computed and support) else pd.NA),
          observed_support=support if computed else pd.NA,observed_value=value,MVP_statistics_used=mvp,status=status))
    for row in evidence.itertuples():
        for layer,count in [('bulk',int(bool(row.bulk_tested))),('cell_Wilcoxon',row.cell_Wilcoxon_types_tested),('pseudobulk',row.donor_pseudobulk_types_tested),('mRNA_coloc',row.mRNA_coloc_computed_regions),('protein_coloc',row.protein_coloc_computed_regions)]:
            missingness.append(dict(gene=row.gene,layer=layer,computed_units=count,status='observed' if count>0 else 'not_evaluable',MVP_statistics_used=layer in ['mRNA_coloc','protein_coloc']))
        missingness.append(dict(gene=row.gene,layer='UKB_PGS_Cox',computed_units=0,status='omitted_by_user',MVP_statistics_used=False))
        cell(row.gene,'Expression','Bulk',int(bool(row.bulk_tested)),1,row.bulk_significant,row.bulk_log2FC)
        cell(row.gene,'Expression','Cell Wilcoxon',row.cell_Wilcoxon_types_tested,6,row.cell_Wilcoxon_any_support)
        cell(row.gene,'Expression','Pseudobulk',row.donor_pseudobulk_types_tested,6,row.donor_pseudobulk_any_support)
        for tissue in sorted(smr.tissue.unique()):
            s=smr[(smr.gene_id==row.gene_id)&(smr.tissue==tissue)]
            known=s.nominal_SMR_and_HEIDI_rule.dropna();supported=bool(known.any()) if len(known) else False
            cell(row.gene,'SMR / HEIDI',tissue,len(known),max(1,len(s)),supported,mvp=True)
            missingness.append(dict(gene=row.gene,layer='SMR_'+tissue,computed_units=len(known),status='observed' if len(known) else 'not_evaluable',MVP_statistics_used=True))
        for layer,planned in [('coloc_mRNA',8),('coloc_protein',2)]:
            z=biological[layer];z=z[(z.gene_id==row.gene_id)&z['PP.H4.abf'].notna()]
            cell(row.gene,'Colocalization',layer,len(z),planned,bool((z['PP.H4.abf']>.8).any()),z['PP.H4.abf'].max() if len(z) else np.nan,True)
    save(pd.DataFrame(missingness),'evidence_coverage')
    save(pd.DataFrame(matrix),'evidence_matrix')
    catalog=[];dictionary=[]
    for p in sorted(data.rglob('*.tsv')):
        if p.name in ('table_catalog.tsv','column_catalog.tsv'):continue
        with p.open(newline='') as f:
            reader=csv.reader(f,delimiter='\t');header=next(reader);count=sum(1 for _ in reader)
        name=str(p.relative_to(data));catalog.append(dict(file=name,rows=count,columns=len(header),sha256=sha(p)))
        dictionary.extend(dict(file=name,position=i+1,column=c) for i,c in enumerate(header))
    save(pd.DataFrame(catalog),'table_catalog');save(pd.DataFrame(dictionary),'column_catalog')
    result=dict(status='corrected_primary_integrated_with_declared_source_uncertainty',primary_candidates=len(selected),primary_genes=selected.gene.tolist(),
       MVP_estimable=int(primary.p.notna().sum()),external_opportunity_background=int(primary.background_primary.sum()),
       clinical_pairs=len(clinical),clinical_estimable=int(clinical.p.notna().sum()),clinical_BH374=int((clinical.BH_original_family_374<.05).sum()),
       clinical_BH_candidate_exploratory=int((clinical.BH_candidate_family_exploratory<.05).sum()),
       clinical_MRE_BH374=int((clinical.MRE_BH_original_family_374<.05).sum()),
       drug_direct_rows=len(direct),drug_component_rows=len(components),network_drug_rows=len(links),STRING_induced_edges=len(net),STRING_unmapped=missing,
       steiger_selected_instruments=len(chosen_pairs),steiger_grid_rows=len(steiger),steiger_all_forward=bool(steiger.forward_direction.all()),
       multi_instrument_rows=len(multi),multi_instrument_computed=int((multi.status=='computed').sum()),
       enrichment=enrichment,ontology_header=ontology_header,gaf_header=gaf_header,unresolved_GO_terms=sorted(unresolved),
       MVP_free_control=[dict(branch=b,model=m,estimable=int(g.FinnGen_p.notna().sum()),discovery=int(g.FinnGen_discovery.sum()),Stuart_concordant_followup=int(g.control_followup.sum())) for (b,m),g in control.groupby(['branch','model'])],
       original_results_preserved=True,scientific_inference_fully_validated=False,MVP_source_convention_resolved=False,
       sources_are_not_independent_votes=True,total_reconstructed_score_computed=False,
       inputs={str(p.relative_to(ROOT)):sha(p) for p in sorted(set(inputs))},script_sha256=sha(Path(__file__)),
       outputs={str(p.relative_to(out)):sha(p) for p in sorted(data.rglob('*.tsv'))})
    (out/'run.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['inputs','outputs','gaf_header','unresolved_GO_terms']},ensure_ascii=False,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output');main(p.parse_args().output)
