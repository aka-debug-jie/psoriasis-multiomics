"""Retrieve exact frozen instrument positions from indexed FinnGen R12 GWAS."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
import argparse,hashlib,json,os
from datetime import datetime,timezone
import pandas as pd
import requests
import pysam
from audit_mr_selection import digest

ROOT=Path(__file__).resolve().parents[2]
CODES=['K11_CD_STRICT2','K11_IBD_STRICT','K11_UC_STRICT2','NAFLD','M13_RHEUMA','I9_HEARTFAIL','C_STROKE','T2D','C3_OTHER_SKIN_EXALLC','C3_MELANOMA_SKIN_EXALLC','I9_CORATHER','I9_HYPTENS','E4_THYROIDITAUTOIM','E4_OBESITY','M13_SJOGREN','SLE_FG','M13_GOUT','E4_GRAVES_STRICT','L12_PSORIASIS','M13_PSORIARTH']

def main(include_proxies=False):
    os.environ['CURL_CA_BUNDLE']='/etc/ssl/certs/ca-certificates.crt'
    os.environ['SSL_CERT_FILE']='/etc/ssl/certs/ca-certificates.crt'
    out=ROOT/'research/data/processed/clinical_finngen';out.mkdir(parents=True,exist_ok=True)
    manifest_path=ROOT/'research/data/metadata/finngen_R12_manifest.tsv'
    manifest=pd.read_csv(manifest_path,sep='\t')
    follow_path=ROOT/'research/results/mr_fenland_substitute/followup_proteins.tsv'
    iv_path=ROOT/'research/results/mr_fenland_substitute/mvp_cohort_wald.tsv'
    coordinate_path=ROOT/'research/data/processed/substitute_mr/finngen_outcome.tsv'
    follow=pd.read_csv(follow_path,sep='\t');iv=pd.read_csv(iv_path,sep='\t')
    coord=pd.read_csv(coordinate_path,sep='\t')[['rsid','chr','pos']]
    iv=iv[iv.uniprot.isin(follow.uniprot)]
    requested=iv[['rsid']].drop_duplicates().merge(coord,on='rsid',how='left',validate='one_to_one')
    requested['coordinate_status']=requested.pos.notna().map({True:'exact_rsid_GRCh38_from_R12_psoriasis',False:'no_GRCh38_position_no_proxy'})
    requested.to_csv(out/'requested_instruments.tsv',sep='\t',index=False,na_rep='NA')
    wanted=requested.dropna(subset=['chr','pos']).copy();wanted.chr=wanted.chr.astype(int).astype(str);wanted.pos=wanted.pos.astype(int)
    codes=CODES+(['I9_CHD','E4_HYPERLIPNAS'] if include_proxies else [])
    panel=manifest[manifest.phenocode.isin(codes)].copy();assert len(panel)==len(codes)
    panel['scope']='public_R12_reconstruction_original_panel_accessions_unreported'
    panel['proxy_requires_authorization']=panel.phenocode.isin(['I9_CHD','E4_HYPERLIPNAS'])
    panel.to_csv(out/'phenotype_manifest.tsv',sep='\t',index=False)
    request_hash=digest(out/'requested_instruments.tsv')
    def one(item):
        code=item['phenocode'];url=item['path_https'];directory=out/code;directory.mkdir(exist_ok=True)
        recpath=directory/'retrieval.json';table=directory/'regional_rows.tsv';idxpath=directory/'source.gz.tbi'
        if recpath.exists():
            old=json.loads(recpath.read_text())
            if old.get('status')=='downloaded' and old.get('request_sha256')==request_hash and old.get('table_sha256')==digest(table):return old
        attempts=[]
        for attempt in range(2):
            try:
                head=requests.head(url,timeout=60);head.raise_for_status()
                index=requests.get(url+'.tbi',timeout=60);index.raise_for_status();idxpath.write_bytes(index.content)
                handle=pysam.TabixFile(url,index=str(idxpath))
                header=handle.header[-1].lstrip('#').split('\t');rows=[];queries=[]
                for chrom,pos in sorted(set(zip(wanted.chr,wanted.pos))):
                    found=list(handle.fetch(chrom,pos-1,pos))
                    queries.append(dict(chr=chrom,pos=pos,records=len(found)))
                    rows.extend(line.split('\t') for line in found)
                handle.close()
                after=requests.head(url,timeout=60);after.raise_for_status()
                assert head.headers.get('ETag')==after.headers.get('ETag'),'Source changed during extraction'
                d=pd.DataFrame(rows,columns=header).drop_duplicates();d.to_csv(table,sep='\t',index=False)
                r=dict(status='downloaded',phenocode=code,url=url,index_url=url+'.tbi',accessed_utc=datetime.now(timezone.utc).isoformat(),
                       request_sha256=request_hash,source_etag=head.headers.get('ETag'),source_generation=head.headers.get('x-goog-generation'),source_hash_header=head.headers.get('x-goog-hash'),
                       source_bytes=head.headers.get('Content-Length'),index_sha256=digest(idxpath),table_sha256=digest(table),rows=len(d),
                       queries=queries,method='pysam 0.23.3 / HTSlib remote BGZF+Tabix, 0-based half-open exact point queries',attempts=attempts)
                recpath.write_text(json.dumps(r,indent=2)+'\n');return r
            except Exception as exc:attempts.append(str(exc))
        r=dict(status='retrieval_failed',phenocode=code,url=url,attempts=attempts);recpath.write_text(json.dumps(r,indent=2)+'\n');return r
    results=[]
    with ThreadPoolExecutor(max_workers=3) as executor:
        jobs=[executor.submit(one,item) for item in panel.to_dict('records')]
        for job in as_completed(jobs):
            r=job.result();results.append(r);print(r['phenocode'],r['status'],r.get('rows'),flush=True)
    record=dict(status='complete' if all(r['status']=='downloaded' for r in results) else 'partial_retrieval',pysam=pysam.__version__,
                script_sha256=digest(Path(__file__)),inputs={str(p.relative_to(ROOT)):digest(p) for p in [manifest_path,follow_path,iv_path,coordinate_path]},
                requested=len(requested),coordinate_available=len(wanted),outcomes=[dict(phenocode=r['phenocode'],status=r['status']) for r in results])
    (out/'run.json').write_text(json.dumps(record,indent=2)+'\n')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--include-authorized-proxies',action='store_true');a=p.parse_args();main(a.include_authorized_proxies)
