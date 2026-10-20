"""Independent verification of cell identity, embeddings, provenance and replay."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]


def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda:stream.read(1024*1024),b''):
            h.update(chunk)
    return h.hexdigest()


def main():
    out=ROOT/'research/results/sc_embedding'
    run=json.loads((out/'run.json').read_text())
    result=run['result']
    assert result['status']=='embedding_completed_annotation_pending_review',result
    assert result['formal_candidate_list_used'] is False
    for p,h in run['inputs'].items():
        assert digest(Path(p))==h,p
    for p,h in run['outputs'].items():
        assert digest(out/p)==h,p
    prep=ROOT/'research/data/processed/sc_analysis'
    merged=json.loads((prep/'run.json').read_text())
    for p,h in merged['output_sha256'].items():
        assert digest(prep/p)==h,p
    assert digest(ROOT/'research/scripts/prepare_sc_matrix.py')==merged['script_sha256']
    for entry in merged['input_sha256']:
        assert digest(ROOT/entry['path'])==entry['sha256'],entry['path']
    source=pd.read_csv(prep/'cells.tsv',sep='\t')
    qc=pd.read_csv(ROOT/'research/results/sc_qc/cell_qc.tsv.gz',sep='\t')
    qc=qc[qc.keep].reset_index(drop=True)
    ids=qc.gsm+'_'+qc.barcode
    assert source.cell_id.tolist()==ids.tolist()
    cells=pd.read_csv(out/'cell_embeddings.tsv',sep='\t')
    assert cells.cell_id.tolist()==source.cell_id.tolist()
    assert not cells.cell_id.duplicated().any()
    for column in ['gsm','donor_id','tissue_type','n_counts','n_features']:
        assert cells[column].tolist()==source[column].tolist(),column
    assert np.allclose(cells.percent_mito,source.percent_mito,rtol=1e-8,atol=1e-10)
    assert np.isfinite(cells[['UMAP_1','UMAP_2']].to_numpy()).all()
    pca=pd.read_csv(out/'cell_pca.tsv',sep='\t')
    assert pca.cell_id.tolist()==source.cell_id.tolist()
    assert pca.shape==(len(cells),31)
    assert np.isfinite(pca.iloc[:,1:].to_numpy()).all()
    features=pd.read_csv(prep/'features.tsv',sep='\t')
    variable=pd.read_csv(out/'variable_features.tsv',sep='\t')
    assert not features.gene_id.duplicated().any()
    assert set(variable.gene_id).issubset(features.gene_id)
    assert len(variable)==2000 and not variable.gene_id.duplicated().any()
    assert len(cells)==result['cells']==merged['cells']==len(qc)
    assert len(features)==result['genes']==merged['genes_after_global_zero_filter']
    assert cells.cluster.nunique()==result['clusters']
    counts=pd.read_csv(out/'cluster_donor_counts.tsv',sep='\t')
    assert int(counts.Freq.sum())==len(cells)
    expected=cells.groupby(['cluster','donor_id','tissue_type']).size().sort_index()
    observed=counts[counts.Freq>0].set_index(['cluster','donor_id','tissue_type']).Freq.sort_index()
    assert expected.index.equals(observed.index)
    assert (expected.values==observed.values).all()
    approved=json.loads((ROOT/'research/config/reconstruction_defaults.json').read_text())['single_cell']
    for key,value in approved.items():
        assert result['settings'][key]==value,key
    assert result['future_workers']==1
    previous=json.loads((ROOT/'research/records/optimization_20261017/sc_first_run_tables.json').read_text())
    replay={p:digest(out/p)==h for p,h in previous.items()}
    assert all(replay.values()),replay
    report=dict(status='passed',cells=len(cells),genes=len(features),clusters=int(cells.cluster.nunique()),
                donors=int(cells.donor_id.nunique()),sample_groups=cells.groupby('tissue_type').size().to_dict(),
                independent_cell_identity_check=True,input_output_hashes_checked=True,
                two_full_computations_identical_tables=replay,full_single_cell_reproduction=False)
    (ROOT/'research/records/sc_embedding_validation.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    main()
