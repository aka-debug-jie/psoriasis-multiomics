"""Read only the public dbGaP metadata and first gzip block via checked HTTP ranges."""
from pathlib import Path
import hashlib,json,zlib,csv,io
import requests

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'research/records/closure_20261017/mvp_publication'
URL='https://ftp.ncbi.nlm.nih.gov/dbgap/studies/phs002453/analyses/GIA/phs002453.MVP_R4.1000G_AGR.GIA.PheCodes_DermatologicalSystem.analysis-PI.MULTI.tar'


def get_range(start,length):
    with requests.get(URL,headers={'Range':f'bytes={start}-{start+length-1}'},stream=True,timeout=(20,60)) as r:
        if r.status_code!=206:raise RuntimeError(f'Expected HTTP206, received {r.status_code}; full archive download refused')
        assert r.headers['Content-Range'].startswith(f'bytes {start}-{start+length-1}/')
        b=r.raw.read(length+1);assert len(b)==length
        return b,dict(range=r.headers['Content-Range'],etag=r.headers.get('ETag'),last_modified=r.headers.get('Last-Modified'))


def main():
    OUT.mkdir(parents=True,exist_ok=True);records=[]
    for kind,offset in [('metadata',49126383104),('summary',166804524544)]:
        header,hmeta=get_range(offset,512)
        name=header[:100].rstrip(b'\0').decode();size=int(header[124:136].strip(b'\0 '),8)
        assert '696_41' in name and '.EUR.GIA' in name,name
        # TAR checksum treats its checksum bytes as spaces.
        assert sum(header[:148])+8*32+sum(header[156:])==int(header[148:156].strip(b'\0 '),8)
        length=size if kind=='metadata' else 65536
        payload,meta=get_range(offset+512,length)
        (OUT/(kind+'_tar_header.bin')).write_bytes(header)
        (OUT/(kind+'_range.bin')).write_bytes(payload)
        if kind=='metadata':(OUT/'pha006164_metadata.txt').write_bytes(payload)
        else:
            decoded=zlib.decompressobj(16+zlib.MAX_WBITS).decompress(payload,2*1024*1024)
            lines=decoded.decode().splitlines()[:3];assert len(lines)==3
            (OUT/'original_first_records.tsv').write_text('\n'.join(lines)+'\n')
        records.append(dict(kind=kind,tar_header_name=name,tar_header_name_may_be_truncated=len(name)==100,archive_member_size=size,tar_header_offset=offset,
           body_offset=offset+512,downloaded_body_bytes=length,range_sha256=hashlib.sha256(payload).hexdigest(),
           tar_header_checksum_valid=True,HTTP=meta,full_member_downloaded=kind=='metadata',full_summary_GZIP_CRC_verified=False))
    record=dict(status='public_source_metadata_and_first_two_records_observed',archive_url=URL,records=records,
      source_scope='No full tar or 592MB summary download. Gzip first-block decompression succeeds; whole-member CRC/hash is not claimed.',
      interpretation='Binary-trait source publishes OR, CI and P; no native beta/SE columns in observed header',
      files={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in OUT.iterdir() if p.is_file() and p.name!='source.json'})
    (OUT/'source.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record,indent=2))


if __name__=='__main__':main()
