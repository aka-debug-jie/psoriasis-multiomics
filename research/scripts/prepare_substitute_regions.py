"""Use version-matched gene identities for substituted and author-reference targets."""
from pathlib import Path
import json
import pandas as pd
from audit_mr_selection import gene_locations,digest
from gene_identity import identity_resolver
ROOT=Path(__file__).resolve().parents[2]

def main():
    out=ROOT/'research/data/processed/coloc_substitute_registry';out.mkdir(parents=True,exist_ok=True)
    follow=ROOT/'research/results/mr_fenland_substitute/followup_proteins.tsv'
    author=ROOT/'research/references/author_main_table1.tsv'
    gtf=ROOT/'research/data/metadata/gencode.v19.annotation.gtf.gz'
    geneinfo=ROOT/'research/data/metadata/Homo_sapiens.gene_info.gz'
    loc=gene_locations(gtf);loc=loc[loc.gene_chr.isin(['chr'+str(i) for i in range(1,23)])].copy()
    loc.gene_id=loc.gene_id.str.split('.').str[0];byid=loc.set_index('gene_id')
    resolve,_=identity_resolver(geneinfo,loc)
    genes={}
    for _,r in pd.read_csv(follow,sep='\t').iterrows():
        genes[r.gene_id]=dict(gene_id=r.gene_id,symbol=byid.loc[r.gene_id,'gene'],available_followup=True,author_reference=False)
    for target in pd.read_csv(author,sep='\t').target:
        name={'IL23':'IL23A','IL12':'IL12A'}.get(target,target);gid=resolve(name)
        assert gid is not None,('Unresolved author gene',name)
        if gid not in genes:genes[gid]=dict(gene_id=gid,symbol=byid.loc[gid,'gene'],available_followup=False,author_reference=True)
        else:genes[gid]['author_reference']=True
    reg=pd.DataFrame(genes.values()).sort_values('gene_id');reg.to_csv(out/'gene_registry.tsv',sep='\t',index=False)
    targets=ROOT/'research/data/processed/Fenland_Olink/target_registry.tsv'
    t=pd.read_csv(targets,sep='\t');wanted=set(t.loc[t.ensembl_gene_id.isin(genes),'id'])
    raw=ROOT/'research/data/raw/Fenland_Olink/Summary.statistics.cis.Olink.proteins.Fenland.20220107.txt.gz'
    parts=[]
    for chunk in pd.read_csv(raw,sep='\t',chunksize=200000):
        keep=chunk[chunk.id.isin(wanted)]
        if len(keep):parts.append(keep)
    full=pd.concat(parts,ignore_index=True);full.to_csv(out/'Fenland_full_regions.tsv.gz',sep='\t',index=False)
    record=dict(status='registered',genes=len(reg),substitute_followup=int(reg.available_followup.sum()),author_reference=int(reg.author_reference.sum()),Fenland_rows=len(full),
                identity='Unique GENCODE v19 primary gene; current/legacy NCBI symbols resolved through unique GeneID. MMP12 uses GRCh37 source ENSG00000110347, not current ENSG00000262406.',
                inputs={str(p.relative_to(ROOT)):digest(p) for p in [follow,author,gtf,geneinfo,targets,raw]},script_sha256=digest(Path(__file__)),outputs={p.name:digest(p) for p in out.iterdir() if p.name!='run.json'})
    (out/'run.json').write_text(json.dumps(record,indent=2)+'\n');print(record['genes'],record['Fenland_rows'])

if __name__=='__main__':main()
