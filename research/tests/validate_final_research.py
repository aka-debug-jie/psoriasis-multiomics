"""Validate the unified candidate set, adjusted families, graph and preserved biological evidence."""
from pathlib import Path
import argparse,csv,json,hashlib,math,re,subprocess
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from validate_release import png_info

ROOT=Path(__file__).resolve().parents[2]


def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
    return h.hexdigest()


def close(a,b):np.testing.assert_allclose(a,b,rtol=1e-8,atol=1e-10,equal_nan=True)


def bh(p,n=None):
    p=np.asarray(p,float);n=n or len(p);good=np.isfinite(p);out=np.full(len(p),np.nan)
    indices=np.flatnonzero(good);indices=indices[np.argsort(p[good])];running=1.
    for rank,index in reversed(list(enumerate(indices,1))):
        running=min(running,p[index]*n/rank);out[index]=running
    return out


def main(output='research/final',repeat='research/records/closure_20261017/fresh_final_complete'):
    out=ROOT/output;data=out/'data';record=json.loads((out/'run.json').read_text());cfg=json.loads((ROOT/'research/config/primary_analysis.json').read_text())
    for n,h in record['inputs'].items():assert sha(ROOT/n)==h,n
    for n,h in record['outputs'].items():assert sha(out/n)==h,n
    read=lambda n:pd.read_csv(data/(n+'.tsv'),sep='\t')
    catalog=read('table_catalog')
    for r in catalog.itertuples():
        p=data/r.file;assert sha(p)==r.sha256
        with p.open(newline='') as f:
            rows=csv.reader(f,delimiter='\t');header=next(rows);count=sum(1 for _ in rows)
        assert count==r.rows and len(header)==r.columns and len(set(header))==len(header)
    full=read('primary_all_proteins');candidates=read('primary_candidates');e=read('integrated_evidence');uids=set(candidates.uniprot)
    expected=(full.p<.05/2709)&((full.finngen_p<.05)|(full.stuart_p<.05))
    assert len(full)==2709 and set(full.loc[expected,'uniprot'])==uids==set(e.uniprot)
    assert not full.uniprot.duplicated().any() and not e.gene_id.duplicated().any()
    assert not {'passes_frozen_bonferroni','externally_supported','original_p','original_followup'}&set(e.columns)
    source=pd.read_csv(ROOT/'research/results/covariance_resolution/all_proteins_comparison.tsv',sep='\t')
    source=source[(source.branch==cfg['branch'])&(source.scenario==cfg['MVP_scenario'])].set_index('uniprot')
    for field in ('beta','se','p','finngen_beta','finngen_p','stuart_beta','stuart_p'):
        close(e[field],source.loc[e.uniprot,field])
    opportunity=full.p.notna()&(full.finngen_p.notna()|full.stuart_p.notna())
    assert np.array_equal(full.background_primary,opportunity)
    mapping=read('enrichment_background_mapping');assert set(mapping.loc[mapping.background_primary,'uniprot'])==set(full.loc[opportunity,'uniprot'])
    checked_terms=0;enrichment_counts={}
    for scope in ('primary','discovery_sensitivity'):
        enrichment_counts[scope]={}
        for library in ('biological_process','molecular_function','cellular_component','KEGG'):
            d=read('enrichment/'+scope+'/'+library);ps=[]
            for r in d.itertuples():
                k,q,K,N=map(int,[r.overlap,r.query_annotated,r.term_background,r.background_annotated])
                ps.append(fisher_exact([[k,q-k],[K-k,N-K-q+k]],alternative='greater').pvalue)
            close(d.p_value,ps);close(d.fdr_bh,bh(ps));assert (d.overlap==0).any()
            assert len(d)==record['enrichment'][scope][library]['tested_terms']
            assert d.term_background.between(*cfg['enrichment_term_size']).all()
            checked_terms+=len(d);enrichment_counts[scope][library]=int((d.fdr_bh<.05).sum())
    cl=read('clinical_primary');assert len(cl)==len(uids)*22 and not cl.duplicated(['uniprot','phenocode']).any()
    close(cl.BH_candidate_family_exploratory,bh(cl.p));close(cl.MRE_BH_candidate_family_exploratory,bh(cl.p_MRE))
    old=pd.read_csv(ROOT/'research/results/covariance_resolution/clinical_original_17.tsv',sep='\t')
    old=old[old.branch==cfg['branch']];close(old.BH_fixed_374,bh(old.p,374));old=old.set_index(['uniprot','outcome'])
    ix=pd.MultiIndex.from_frame(cl[['uniprot','phenocode']]);close(cl.p,old.loc[ix,'p']);close(cl.BH_original_family_374,old.loc[ix,'BH_fixed_374'])
    legacy=pd.read_csv(ROOT/'research/results/evidence_substitute/all_registered_evidence.tsv',sep='\t').set_index('gene_id')
    for field in [c for c in e if c.startswith(('bulk_','cell_Wilcoxon_','donor_pseudobulk_','SMR_','HEIDI_','mRNA_coloc_','protein_coloc_','author_'))]:
        a=e[field].reset_index(drop=True);b=legacy.loc[e.gene_id,field].reset_index(drop=True)
        if pd.api.types.is_numeric_dtype(a) and pd.api.types.is_numeric_dtype(b):close(a,b)
        else:pd.testing.assert_series_equal(a,b,check_names=False,check_dtype=False)
    assert e[['UKB_PGS_psoriasis','UKB_PGS_PsA','reconstructed_total_score']].isna().all().all()
    assert e.loc[e.gene=='DDX58','bulk_tested'].all()
    assert e.loc[e.gene=='IL12B','donor_pseudobulk_any_support'].isna().all()
    multi=read('multi_instrument_methods');original_multi=pd.read_csv(ROOT/'research/results/multi_instrument_substitute/method_results.tsv',sep='\t');original_multi=original_multi[original_multi.uniprot.isin(uids)]
    pd.testing.assert_frame_equal(multi.drop(columns='interpretation').reset_index(drop=True),original_multi.reset_index(drop=True),check_dtype=False,rtol=1e-8,atol=1e-10)
    steiger=read('steiger_selected_instruments');assert set(steiger.uniprot)==uids
    assert (steiger.groupby(['instrument_id','rsid']).size()==8).all()
    assert steiger.forward_direction.eq(steiger.R2_exposure>steiger.R2_outcome).all()
    matrix=read('evidence_matrix');assert matrix.loc[matrix.status.isin(['not_evaluable','partial_no_observed_support']),'support'].isna().all()
    assert matrix.loc[matrix.domain.isin(['SMR / HEIDI','Colocalization']),'MVP_statistics_used'].all()
    relations=read('drugs/all_drug_relations');assert np.array_equal(relations.candidate,relations.uniprot.isin(uids))
    net=read('drugs/STRING_induced_network');nodes=read('drugs/STRING_mapping');ids=set(nodes.stringId)
    rawpath=ROOT/'research/data/raw/functional/string_substitute/network.json';raw=pd.DataFrame(json.loads(rawpath.read_text()))
    expect=raw[raw.stringId_A.isin(ids)&raw.stringId_B.isin(ids)&(raw.score>.4)]
    net['ncbiTaxonId']=net.ncbiTaxonId.astype(str);expect=expect.copy();expect['ncbiTaxonId']=expect.ncbiTaxonId.astype(str)
    keys=['stringId_A','stringId_B'];pd.testing.assert_frame_equal(net.sort_values(keys).reset_index(drop=True),expect.sort_values(keys).reset_index(drop=True),check_dtype=False,rtol=1e-8,atol=1e-10)
    assert set(read('drugs/STRING_unmapped').uniprot)=={'P11473'}
    control=read('MVP_free_control_all');native=pd.read_csv(ROOT/'research/results/covariance_resolution/all_outcome_results.tsv',sep='\t')
    native_mre=pd.read_csv(ROOT/'research/results/covariance_resolution/heterogeneity/all_outcome_MRE.tsv',sep='\t')
    for (branch,model),g in control.groupby(['branch','model']):
        src=native if model=='fixed_effect' else native_mre;pcol='p' if model=='fixed_effect' else 'p_MRE'
        f=src[(src.branch==branch)&(src.outcome=='finngen')].set_index('uniprot');s=src[(src.branch==branch)&(src.outcome=='stuart')].set_index('uniprot')
        assert len(g)==2709;close(g.FinnGen_p,f.loc[g.uniprot,pcol]);close(g.Stuart_p,s.loc[g.uniprot,pcol])
        rule=(g.FinnGen_p<.05/2709)&(g.Stuart_p<.05)&(g.FinnGen_beta*g.Stuart_beta>0)
        assert set(g.loc[rule,'uniprot'])==set(g.loc[g.control_followup.astype('boolean').fillna(False),'uniprot'])
    metadata_comparison=json.loads((ROOT/'research/records/closure_20261017/mvp_publication/catalog_comparison.json').read_text())
    assert all(all(r['field_matches'].values()) for r in metadata_comparison['rows'])
    assert metadata_comparison['MVP_convention_resolved'] is False
    for n,h in metadata_comparison['inputs'].items():assert sha(ROOT/n)==h
    # A genuinely new output directory was used for this comparison, not --resume.
    fresh=ROOT/repeat;repeat_run=json.loads((fresh/'run.json').read_text());assert record['outputs']==repeat_run['outputs']
    for n,h in repeat_run['outputs'].items():assert sha(fresh/n)==h
    figure=json.loads((out/'figures_manifest.json').read_text());assert len(figure['figures'])==4
    for n,h in figure['inputs'].items():assert sha(ROOT/n)==h,n
    for r in figure['panel_tables'].values():assert sha(out/r['file'])==r['sha256']
    for name,r in figure['figures'].items():
        for n,h in r['outputs'].items():assert sha(out/n)==h
        pdf=out/'figures'/(name+'.pdf');info=subprocess.check_output(['pdfinfo',str(pdf)],text=True)
        assert re.search(r'Pages:\s+1\b',info)
        w,h=map(float,re.search(r'Page size:\s+([\d.]+) x ([\d.]+)',info).groups())
        assert w==math.floor(72*r['width_in']+1e-8) and h==math.floor(72*r['height_in']+1e-8)
        fonts=subprocess.check_output(['pdffonts',str(pdf)],text=True).splitlines()[2:];assert fonts and all(x.split()[-5]=='yes' for x in fonts)
        png=png_info(out/'figures'/(name+'.png'));assert png['width']==round(r['width_in']*600) and png['height']==round(r['height_in']*600)
        assert png.get('unit')==1 and abs(png['xppm']*.0254-600)<.1
    report=dict(status='passed',primary_candidates=len(e),opportunity_background=int(opportunity.sum()),enrichment_terms_independently_checked=checked_terms,
        enrichment_significant=enrichment_counts,clinical_pairs=len(cl),BH374=int((cl.BH_original_family_374<.05).sum()),BH_candidate_exploratory=int((cl.BH_candidate_family_exploratory<.05).sum()),
        legacy_flags_explicit=True,source_biological_statistics_preserved=True,missingness_preserved=True,MVP_free_control_verified=True,
        original_source_records_matched=2,source_convention_resolved=False,induced_STRING_edges=len(net),
        multi_instrument_rows_preserved=len(multi),steiger_selected_grid_rows=len(steiger),
        fresh_directory_tables_identical=len(record['outputs']),result_tables=len(catalog),fresh_output_root=repeat,publication_figures=4,
        run_sha256=sha(out/'run.json'),validator_sha256=sha(Path(__file__)))
    (ROOT/'research/records/closure_20261017/validation.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',default='research/final');p.add_argument('--repeat',default='research/records/closure_20261017/fresh_final_complete');a=p.parse_args();main(a.output,a.repeat)
