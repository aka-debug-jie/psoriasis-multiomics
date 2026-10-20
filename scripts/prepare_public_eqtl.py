"""Prepare licensed GTEx/eQTL Catalogue examples on a documented GRCh37 bridge."""
from pathlib import Path
import argparse,json,hashlib,subprocess,gzip
import pandas as pd
from pyliftover import LiftOver
p=argparse.ArgumentParser();p.add_argument('--course',type=Path,required=True);p.add_argument('--records',type=Path,required=True);p.add_argument('--chain',type=Path,required=True);p.add_argument('--mvp',type=Path,required=True);p.add_argument('--ld',type=Path,required=True);p.add_argument('--plink',type=Path,required=True);p.add_argument('--smr',type=Path,required=True);a=p.parse_args()
folder=a.course/'data/smr';folder.mkdir(exist_ok=True)
lo=LiftOver(str(a.chain));complement=str.maketrans('ACGT','TGCA')
genes={'ENSG00000164308':'ERAP2','ENSG00000107798':'LIPA'};datasets={};exclusions=[];flist=[]
for gene,symbol in genes.items():
 d=pd.read_csv(a.records/(gene+'.tsv'),sep='\t');meta=json.loads((a.records/(gene+'.json')).read_text());tss=meta['start']if meta['strand']==1 else meta['end']
 tss37=lo.convert_coordinate('chr'+str(meta['seq_region_name']),tss-1)
 assert len(tss37)==1
 records=[]
 for r in d.itertuples():
  mapped=lo.convert_coordinate('chr'+str(r.chromosome),int(r.position)-1)
  if r.type!='SNP' or len(r.ref)!=1 or len(r.alt)!=1 or len(mapped)!=1 or not str(r.rsid).startswith('rs'):
   exclusions.append({'gene':gene,'variant':r.variant,'reason':'not a uniquely mapped biallelic SNP with rsID'});continue
  chrom,pos,strand,_=mapped[0]
  alt,ref=(r.alt.translate(complement),r.ref.translate(complement))if strand=='-'else(r.alt,r.ref)
  records.append({'Chr':chrom.removeprefix('chr'),'SNP':r.rsid,'BP':int(pos)+1,'A1':alt,'A2':ref,'Freq':r.ac/r.an,'b':r.beta,'se':r.se,'p':r.pvalue})
 frame=pd.DataFrame(records).sort_values(['Chr','BP','SNP']).drop_duplicates()
 conflicts=frame.SNP.duplicated(keep=False)
 for rs in frame.loc[conflicts,'SNP']:exclusions.append({'gene':gene,'variant':rs,'reason':'duplicate rsID'})
 frame=frame[~conflicts];assert len(frame)>1000 and frame.p.min()<5e-8
 frame.to_csv(folder/(gene+'.esd'),sep='\t',index=False);datasets[gene]=frame
 flist.append({'Chr':tss37[0][0].removeprefix('chr'),'Probe':gene,'GeneticDistance':0,'Probe_bp':int(tss37[0][1])+1,'Gene':symbol,'Orientation':'+'if meta['strand']==1 else'-','ESD':'data/smr/'+gene+'.esd'})
pd.DataFrame(flist).to_csv(folder/'eqtl_flist.txt',sep='\t',index=False)
(folder/'probes.txt').write_text('\n'.join(genes)+'\n')
ids=set().union(*(set(d.SNP)for d in datasets.values()));(a.records/'extract_snps.txt').write_text('\n'.join(sorted(ids))+'\n')
parts=[]
for chunk in pd.read_csv(a.mvp,sep='\t',chunksize=300000):
 sub=chunk[chunk.SNP.isin(ids)]
 if len(sub):parts.append(sub)
mvp=pd.concat(parts,ignore_index=True).drop_duplicates();assert not mvp.SNP.duplicated().any()
mvp.to_csv(folder/'mvp_regional.ma',sep='\t',index=False)
with(a.records/'plink_prepare.log').open('w')as log:
 subprocess.run([str(a.plink.resolve()),'--bfile',str(a.ld.resolve()),'--extract',str((a.records/'extract_snps.txt').resolve()),'--make-bed','--out',str((folder/'EUR_two_regions').resolve())],stdout=log,stderr=subprocess.STDOUT,check=True)
with(a.records/'make_besd.log').open('w')as log:
 subprocess.run([str(a.smr.resolve()),'--eqtl-flist','data/smr/eqtl_flist.txt','--make-besd','--out','data/smr/public_eqtl'],cwd=a.course,stdout=log,stderr=subprocess.STDOUT,check=True)
# Remove the superseded, unlicensed eQTLGen inputs only from this public staging copy.
for name in ['ENSG00000118503.esd','ENSG00000181523.esd','source_order_eqtl.besd','source_order_eqtl.epi','source_order_eqtl.esi','research_reference.tsv']:
 (folder/name).unlink(missing_ok=True)
coloc=a.course/'data/coloc';regions=pd.read_csv(coloc/'regions.tsv',sep='\t');regions=regions[regions.gene=='IL7R']
eqtl=datasets['ENSG00000164308'];index=eqtl.loc[eqtl.p.idxmin()];rows=[]
for window in (500,1000):
 start=int(index.BP)-window*1000;end=int(index.BP)+window*1000
 selected=eqtl[(eqtl.BP>=start)&(eqtl.BP<=end)]
 shared=selected.merge(mvp,on='SNP',suffixes=('_x','_y'))
 # Preserve both allele pairs. The course explicitly harmonizes and excludes ambiguous alleles.
 out=pd.DataFrame({'SNP':shared.SNP,'position':shared.BP,'A1':shared.A1_x,'A2':shared.A2_x,'bx':shared.b_x,'sx':shared.se_x,'EA_y':shared.A1_y,'OA_y':shared.A2_y,'by':shared.b_y,'sy':shared.se_y})
 def compatible(r):
  pair=r.A1+r.A2;other=r.EA_y+r.OA_y
  return set(pair)not in ({'A','T'},{'C','G'}) and (pair==other or pair==other[::-1] or pair==other.translate(complement) or pair==other.translate(complement)[::-1])
 keep=out.apply(compatible,axis=1)
 for rs in out.loc[~keep,'SNP']:exclusions.append({'gene':'ENSG00000164308','variant':rs,'reason':'coloc ambiguous or incompatible alleles','window_kb':window})
 out=out[keep]
 out.to_csv(coloc/f'ERAP2_{window}kb.tsv.gz',sep='\t',index=False,compression={'method':'gzip','mtime':0})
 rows.append({'gene':'ERAP2','molecule':'mRNA','window_kb':window,'index_snp':index.SNP,'requested_start':start,'requested_end':end,'source_min':eqtl.BP.min(),'source_max':eqtl.BP.max(),'extends_source':start<eqtl.BP.min()or end>eqtl.BP.max(),'shared_variants':len(out)})
for file in coloc.glob('SGSH_*.tsv.gz'):file.unlink()
pd.concat([regions,pd.DataFrame(rows)]).to_csv(coloc/'regions.tsv',sep='\t',index=False)
meta={'study':'eQTL Catalogue reprocessed GTEx v8 whole blood','study_id':'QTS000015','dataset_id':'QTD000356','sample_size':670,'license':'CC-BY-4.0','license_url':'https://www.ebi.ac.uk/eqtl/License/','source_url':'https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/QTD000356/QTD000356.all.tsv.gz','source_build':'GRCh38','analysis_build':'GRCh37','effect_allele':'ALT; complemented on negative-strand liftover','effect_allele_frequency':'ac/an','chain_sha256':hashlib.sha256(a.chain.read_bytes()).hexdigest(),'selection':'Teaching examples chosen using exposure eQTL strength only; not a candidate discovery analysis','genes':[{'gene_id':g,'symbol':genes[g],'variants':len(d),'raw_sha256':hashlib.sha256((a.records/(g+'.tsv')).read_bytes()).hexdigest()}for g,d in datasets.items()]}
(folder/'public_example.json').write_text(json.dumps(meta,indent=2)+'\n');(coloc/'mrna_source.json').write_text(json.dumps(meta,indent=2)+'\n')
(a.records/'mapping_exclusions.json').write_text(json.dumps(exclusions,indent=2))
print('Prepared GTEx examples:',[(genes[g],len(d))for g,d in datasets.items()],flush=True)
