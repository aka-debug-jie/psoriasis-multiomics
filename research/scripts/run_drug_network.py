"""Separate ChEMBL single-protein targets, target components, and STRING links."""
from pathlib import Path
import json
from datetime import datetime,timezone
import pandas as pd
import requests
from audit_mr_selection import digest
ROOT=Path(__file__).resolve().parents[2]

def main():
    raw=ROOT/'research/data/raw/functional';out=ROOT/'research/results/drug_network_substitute';out.mkdir(parents=True,exist_ok=True)
    follow_path=ROOT/'research/results/mr_fenland_substitute/followup_proteins.tsv'
    candidates=pd.read_csv(follow_path,sep='\t');uids=set(candidates.uniprot)
    acquisition=raw/'chembl_mechanisms/run.json';run=json.loads(acquisition.read_text());assert run['status']=='complete'
    relations=[];unresolved=[];component_symbols={};inputs=[follow_path,acquisition,raw/'string_version.json']
    for molecule,m in sorted(run['mechanism_manifests'].items()):
        for response in m['responses']:
            path=raw/'chembl_mechanisms'/molecule/response['file'];assert digest(path)==response['sha256'];inputs.append(path)
            for mechanism in json.loads(path.read_text())['mechanisms']:
                tid=mechanism.get('target_chembl_id')
                if tid is None:unresolved.append(dict(molecule=molecule,reason='mechanism_target_unspecified'));continue
                targetpath=raw/'chembl_targets'/tid/'target.json';inputs.append(targetpath)
                assert digest(targetpath)==run['target_manifests'][tid]['sha256']
                target=json.loads(targetpath.read_text())
                for component in target.get('target_components',[]):
                    accession=component.get('accession')
                    if not accession:continue
                    symbols={s['component_synonym'] for s in component.get('target_component_synonyms',[]) if s['syn_type']=='GENE_SYMBOL'}
                    component_symbols.setdefault(accession,set()).update(symbols)
                    direct=target['target_type']=='SINGLE PROTEIN' and component['relationship']=='SINGLE PROTEIN'
                    relations.append(dict(molecule_chembl_id=molecule,target_chembl_id=tid,target_name=target['pref_name'],target_type=target['target_type'],organism=target['organism'],
                                           uniprot=accession,component_relation=component['relationship'],action_type=mechanism.get('action_type'),mechanism=mechanism.get('mechanism_of_action'),
                                           relation='single_protein_target' if direct else 'complex_or_family_component',candidate=accession in uids))
    rel=pd.DataFrame(relations).drop_duplicates();rel.to_csv(out/'all_drug_target_relations.tsv',sep='\t',index=False)
    rel[rel.candidate & (rel.relation=='single_protein_target')].to_csv(out/'candidate_direct_targets.tsv',sep='\t',index=False)
    rel[rel.candidate & (rel.relation!='single_protein_target')].to_csv(out/'candidate_complex_family_components.tsv',sep='\t',index=False)
    pd.DataFrame(unresolved).to_csv(out/'unresolved_mechanisms.tsv',sep='\t',index=False)
    targets=set(rel.loc[rel.organism=='Homo sapiens','uniprot']);query=sorted(uids|targets)
    version=json.loads((raw/'string_version.json').read_text())[0];base=version['stable_address']
    cache=raw/'string_substitute';cache.mkdir(exist_ok=True);requests_log=[]
    def get(name,endpoint,params):
        path=cache/(name+'.json');meta=cache/(name+'.request.json');url=base+'/api/json/'+endpoint
        request=dict(url=url,params=params)
        if path.exists() and meta.exists():
            old=json.loads(meta.read_text())
            if old['request']==request and old['sha256']==digest(path):requests_log.append(old);inputs.extend([path,meta]);return json.loads(path.read_text())
        response=requests.get(url,params=params,timeout=90);response.raise_for_status();payload=response.json()
        if not isinstance(payload,list):raise ValueError(str(payload))
        path.write_bytes(response.content);record=dict(request=request,actual_url=response.url,date=datetime.now(timezone.utc).isoformat(),sha256=digest(path))
        meta.write_text(json.dumps(record,indent=2)+'\n');requests_log.append(record);inputs.extend([path,meta]);return payload
    mapped=get('mapping','get_string_ids',dict(identifiers='\r'.join(query),species=9606,limit=1,echo_query=1))
    mapping=pd.DataFrame(mapped);mapping['mapping_source']='UniProt_accession';mapping['matched_query']=mapping.queryItem
    for uid in sorted(set(query)-set(mapping.queryItem)):
        symbols=component_symbols.get(uid,set())
        if len(symbols)!=1:continue
        symbol=next(iter(symbols))
        alternative=get('mapping_fallback_'+uid,'get_string_ids',dict(identifiers=symbol,species=9606,limit=10,echo_query=1))
        exact=[r for r in alternative if r['preferredName']==symbol and int(r['ncbiTaxonId'])==9606]
        if len({r['stringId'] for r in exact})==1:
            row=exact[0];row['matched_query']=symbol;row['queryItem']=uid;row['mapping_source']='unique_ChEMBL_component_gene_symbol'
            mapping=pd.concat([mapping,pd.DataFrame([row])],ignore_index=True)
    mapping.to_csv(out/'STRING_mapping.tsv',sep='\t',index=False)
    assert not mapping.queryItem.duplicated().any()
    missing=sorted(set(query)-set(mapping.queryItem));(out/'STRING_unmapped.json').write_text(json.dumps(missing,indent=2)+'\n')
    ids=sorted(set(mapping.stringId));network=get('network','network',dict(identifiers='\r'.join(ids),species=9606,required_score=400,add_nodes=0,network_type='functional'))
    net=pd.DataFrame(network);net=net[net.score>.4].copy();net.to_csv(out/'STRING_network.tsv',sep='\t',index=False)
    reverse={sid:set(g.queryItem) for sid,g in mapping.groupby('stringId')}
    links=[]
    for _,edge in net.iterrows():
        for a,b in [(edge.stringId_A,edge.stringId_B),(edge.stringId_B,edge.stringId_A)]:
            for candidate in reverse.get(a,set())&uids:
                for target in reverse.get(b,set())&targets:
                    for _,drug in rel[rel.uniprot==target].iterrows():
                        links.append(dict(candidate_uniprot=candidate,interacting_target_uniprot=target,molecule_chembl_id=drug.molecule_chembl_id,
                                          drug_target_type=drug.target_type,drug_target_relation=drug.relation,STRING_score=edge.score,
                                          evidence='functional_network_association_not_direct_drug_binding'))
    indirect=pd.DataFrame(links).drop_duplicates();indirect.to_csv(out/'candidate_network_drug_associations.tsv',sep='\t',index=False)
    outputs=[p for p in out.iterdir() if p.name!='run.json']
    record=dict(status='ChEMBL_STRING_public_reconstruction_computed',original_DrugBank_snapshot=False,STRING_version=version,
                indication_rule='ChEMBL psoriasis MeSH D011565 max_phase_for_ind=4; not a current regulatory approval claim',
                candidate_count=len(candidates),phase4_parent_molecules=run['phase4_parent_count'],mechanisms=run['mechanism_rows'],unresolved_mechanisms=len(unresolved),
                single_protein_candidate_rows=int((rel.candidate&(rel.relation=='single_protein_target')).sum()),
                complex_family_candidate_rows=int((rel.candidate&(rel.relation!='single_protein_target')).sum()),
                STRING_edges=len(net),network_drug_association_rows=len(indirect),STRING_unmapped=missing,
                inputs={str(p.relative_to(ROOT)):digest(p) for p in sorted(set(inputs))},script_sha256=digest(Path(__file__)),outputs={p.name:digest(p) for p in outputs})
    (out/'run.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps({k:v for k,v in record.items() if k not in ['inputs','outputs']},indent=2))

if __name__=='__main__':main()
