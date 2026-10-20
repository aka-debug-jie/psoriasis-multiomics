#!/usr/bin/env python3
"""Post-hoc method extension: one GWAS observation per variant, signed-LD IVW."""
from pathlib import Path
from datetime import datetime, timezone
import argparse, hashlib, json, subprocess
import platform, scipy
import numpy as np
import pandas as pd
from mr_covariance import allele_sign, deduplicate_sources, freeze_ld_subset, givw, p_equivalent_se, bh_fixed

ROOT = Path(__file__).resolve().parents[2]
CONFIG = ROOT/'research/config/covariance_resolution.json'
CFG = json.loads(CONFIG.read_text())
OUT = ROOT/CFG['output']
INPUT = ROOT/CFG['input']


def sha(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as f:
        for block in iter(lambda: f.read(4*1024*1024), b''): h.update(block)
    return h.hexdigest()


def write(d, name):
    d.to_csv(OUT/name, sep='\t', index=False, na_rep='NA', float_format='%.17g')


def read(p): return pd.read_csv(p, sep='\t')


def prepare():
    OUT.mkdir(parents=True, exist_ok=True)
    x = read(INPUT/'frozen_instruments.tsv'); x = x[x.eligibility == 'eligible'].copy()
    x['instrument_id'] = x.cohort+'::'+x.uniprot
    ids = sorted({s for aliases in x.rsid for s in aliases.split(',')})
    (OUT/'requested_rsids.txt').write_text('\n'.join(ids)+'\n')
    ref = ROOT/CFG['LD_prefix']; plink = ROOT/CFG['plink']; prefix = OUT/'LD_subset'
    commands = [[str(plink), '--bfile', str(ref), '--extract', str(OUT/'requested_rsids.txt'),
                 '--keep-allele-order', '--make-bed', '--out', str(prefix), '--threads', '4'],
                [str(plink), '--bfile', str(prefix), '--keep-allele-order', '--r', 'square',
                 '--out', str(OUT/'LD_signed'), '--threads', '4']]
    for j, command in enumerate(commands):
        with (OUT/f'LD_command_{j}.log').open('w') as f: subprocess.run(command, stdout=f, stderr=subprocess.STDOUT, check=True)
    bim = pd.read_csv(str(prefix)+'.bim', sep=r'\s+', names=['chr','rsid','cm','pos','A1','A2'])
    assert not bim.rsid.duplicated().any()
    bim['LD_index'] = np.arange(len(bim))
    bim['variant_key'] = [f'{c}:{p}:{"/".join(sorted((a,b)))}' for c,p,a,b in zip(bim.chr,bim.pos,bim.A1,bim.A2)]
    write(bim, 'LD_variant_order.tsv')
    r = np.loadtxt(OUT/'LD_signed.ld'); assert r.shape == (len(bim), len(bim))
    byid = bim.set_index('rsid'); records = []
    for _, row in x.iterrows():
        d = row.to_dict(); aliases = str(row.rsid).split(','); found = byid.loc[byid.index.intersection(aliases)]
        d['source_rsid'] = row.rsid
        if found.variant_key.nunique() > 1:
            d.update(selection_status='ambiguous_aliases_multiple_reference_variants', variant_key='unresolved:'+row.instrument_id,
                     LD_index=-1, LD_sign=np.nan)
        elif len(found):
            b = found.sort_index().iloc[0]; d['rsid'] = b.name
            sign = allele_sign(row.EA, row.OA, b.A1, b.A2)
            if str(row.chr) != str(b.chr): sign = np.nan
            d.update(selection_status='resolved_reference_variant' if np.isfinite(sign) else 'reference_allele_or_chromosome_unresolved',
                     variant_key=b.variant_key, LD_index=int(b.LD_index), LD_sign=sign)
        else:
            d.update(rsid=sorted(aliases)[0], selection_status='reference_absent_fixed_rsid',
                     variant_key='rsid:'+sorted(aliases)[0], LD_index=-1, LD_sign=np.nan)
        records.append(d)
    source = pd.DataFrame(records); unique = deduplicate_sources(source)
    write(source, 'exposure_identity.tsv'); write(unique, 'unique_variant_sources.tsv')
    plan = []
    for branch, cutoff in CFG['branches'].items():
        for uid, g in unique.groupby('uniprot', sort=True):
            g = g.sort_values(['F_statistic','cohort','assay','rsid'], ascending=[False,True,True,True]).copy()
            if len(g) == 1: keep, reason = [0], 'single_unique_variant'
            elif (g.LD_index < 0).any() or not np.isfinite(g.LD_sign).all():
                keep, reason = [0], 'single_IV_missing_or_unresolved_LD'
            else:
                ix = g.LD_index.to_numpy(int); s = g.LD_sign.to_numpy()
                keep, reason = freeze_ld_subset(r[np.ix_(ix, ix)]*np.outer(s,s), cutoff, CFG['max_LD_condition_number'])
            g['branch'] = branch; g['selected'] = False; g.iloc[keep, g.columns.get_loc('selected')] = True
            g['LD_selection_reason'] = reason; g['requested_unique_variants'] = len(g)
            plan.append(g)
    plan = pd.concat(plan, ignore_index=True); write(plan, 'frozen_analysis_plan.tsv')
    inputs = [CONFIG, INPUT/'frozen_instruments.tsv', INPUT/'selection.json', Path(__file__), Path(__file__).with_name('mr_covariance.py'), plink]
    inputs += [Path(str(ref)+ext) for ext in ('.bed','.bim','.fam')]
    outputs = [OUT/n for n in ['requested_rsids.txt','LD_subset.bed','LD_subset.bim','LD_subset.fam','LD_signed.ld','LD_variant_order.tsv','exposure_identity.tsv','unique_variant_sources.tsv','frozen_analysis_plan.tsv']]
    record = dict(status='exposure_and_reference_only_plan_frozen_before_outcome_reads', created_utc=datetime.now(timezone.utc).isoformat(),
                  proteins=int(x.uniprot.nunique()), source_rows=len(x), unique_variant_rows=len(unique), LD_reference_variants=len(bim),
                  reference_matching='unique rsID plus biallelic allele/chromosome check; physical reference key collapses aliases; cross-build positions are not directly compared',
                  matrix_alleles='PLINK A1 preserved by --keep-allele-order, then signs transformed to selected exposure EA',
                  commands=commands, inputs={str(p.relative_to(ROOT)):sha(p) for p in inputs}, outputs={p.name:sha(p) for p in outputs})
    (OUT/'selection.json').write_text(json.dumps(record, indent=2)+'\n')
    print(json.dumps({k:v for k,v in record.items() if k not in ('inputs','outputs','commands')}, indent=2))


def fit_outcome(plan, r, y, label, scenario, compatibility, audit_rows, stress):
    rows = []
    assert not y.rsid.duplicated().any()
    y = y[['rsid','EA','OA','beta','se','p_reported']].set_index('rsid')
    for (branch, uid), g in plan[plan.selected].groupby(['branch','uniprot'], sort=True):
        d = g.reset_index(drop=True).rename(columns={k:k+'_x' for k in ('EA','OA','beta')})
        observed = y.reindex(g.rsid).reset_index(drop=True).rename(columns={k:k+'_y' for k in ('EA','OA','beta')})
        d = pd.concat([d, observed], axis=1)
        d['sign'] = [allele_sign(a,b,c,e) for a,b,c,e in zip(d.EA_x,d.OA_x,d.EA_y,d.OA_y)]
        d['y_aligned'] = d.beta_y*d.sign; d['sy_used'] = d.se
        d['source_ci_audit'] = d.rsid.map(compatibility)
        if scenario == 'P_equivalent':
            d['sy_used'] = [p_equivalent_se(b,p) for b,p in zip(d.beta_y,d.p_reported)]
        valid = np.isfinite(d.y_aligned) & np.isfinite(d.sy_used) & (d.sy_used > 0)
        valid &= ~d.selection_status.str.startswith('ambiguous_aliases')
        if scenario == 'CI_compatible_only': valid &= d.source_ci_audit == 'consistent'
        d['exclusion_reason'] = np.where(valid,'',np.where(~np.isfinite(d.y_aligned),'outcome_missing_or_alleles_unresolved',
            np.where(~np.isfinite(d.sy_used)|(d.sy_used<=0),'scenario_SE_undefined',
            np.where(d.selection_status.str.startswith('ambiguous_aliases'),'variant_identity_unresolved','CI_not_confirmed_compatible'))))
        d['included'] = valid; d['outcome'] = label; d['scenario'] = scenario
        audit_rows.append(d)
        v = d[valid].copy(); rec = dict(branch=branch, outcome=label, scenario=scenario, uniprot=uid, gene=g.gene.iloc[0],
              frozen_variants=len(g), observed_variants=len(v), LD_selection_reason=g.LD_selection_reason.iloc[0],
              status='not_estimable', coverage='full' if len(v)==len(g) else 'partial_or_absent')
        if len(v):
            assert not v.variant_key.duplicated().any() and not v.rsid.duplicated().any()
            if len(v) == 1: rr = np.eye(1)
            else:
                ix = v.LD_index.to_numpy(int); signs = v.LD_sign.to_numpy(); assert (ix>=0).all() and np.isfinite(signs).all()
                rr = r[np.ix_(ix, ix)]*np.outer(signs, signs)
            args = [v.beta_x.to_numpy(), v.y_aligned.to_numpy(), v.se_exposure.to_numpy(), v.sy_used.to_numpy(), rr]
            rec.update(givw(*args), status='computed', variants=';'.join(v.rsid), source_cohorts=';'.join(v.cohort))
            if label == 'mvp' and branch == 'unique_variant_givw' and scenario == 'CI_Wald':
                for xx in CFG['exposure_error_correlations']:
                    for xy in CFG['exposure_outcome_error_correlations']:
                        a = givw(*args, rho_xx=xx, rho_xy=xy)
                        stress.append(dict(uniprot=uid, gene=rec['gene'], rho_xx=xx, rho_xy=xy, **a))
        rows.append(rec)
    return pd.DataFrame(rows)


def analyse():
    frozen = json.loads((OUT/'selection.json').read_text())
    for name, h in frozen['inputs'].items(): assert sha(ROOT/name)==h, name
    for name, h in frozen['outputs'].items(): assert sha(OUT/name)==h, name
    plan = read(OUT/'frozen_analysis_plan.tsv'); r = np.loadtxt(OUT/'LD_signed.ld')
    sel = json.loads((INPUT/'selection.json').read_text()); denominator = sel['fixed_testing_denominator']
    assert plan.uniprot.nunique() == denominator
    precision_path = ROOT/'research/results/mvp_precision_audit/mvp_precision_audit.tsv'
    # The audit file name is fixed by the precision runner, not inferred from results.
    precision = read(precision_path)
    compatibility = precision.groupby('original_rsid').status.agg(lambda s:s.iloc[0] if s.nunique()==1 else 'conflicting_source_rows').to_dict()
    results = []; audits = []; stress = []; inputs = [CONFIG, OUT/'selection.json', precision_path]
    inputs += [ROOT/p for p in frozen['inputs']] + [OUT/p for p in frozen['outputs']]
    for label in ('mvp','finngen','stuart'):
        p = INPUT/(label+'_outcome.tsv'); inputs.append(p); y = read(p)
        for scenario in (CFG['mvp_scenarios'] if label=='mvp' else ['reported_SE']):
            results.append(fit_outcome(plan, r, y, label, scenario, compatibility, audits, stress))
    res = pd.concat(results, ignore_index=True); res['passes_fixed_bonferroni'] = (res.p < .05/denominator).astype('boolean').mask(res.p.isna(),pd.NA)
    write(res, 'all_outcome_results.tsv'); write(pd.concat(audits, ignore_index=True), 'harmonization_and_coverage.tsv')
    write(pd.DataFrame(stress), 'delta_error_correlation_grid.tsv')
    comparisons = []
    old = read(ROOT/'research/results/mr_fenland_substitute/all_proteins_comparison.tsv')
    for (branch, scenario), d in res[res.outcome=='mvp'].groupby(['branch','scenario'], sort=True):
        d = d.copy()
        for outcome in ('finngen','stuart'):
            e = res[(res.branch==branch)&(res.outcome==outcome)][['uniprot','beta','p']]
            d = d.merge(e.rename(columns={c:outcome+'_'+c for c in ['beta','p']}), on='uniprot', how='left', validate='one_to_one')
        def below(p): return (p<.05).astype('boolean').mask(p.isna(),pd.NA)
        d['external_nominal'] = below(d.finngen_p)|below(d.stuart_p)
        d['followup'] = d.passes_fixed_bonferroni & d.external_nominal
        def agrees(other):
            return (d.beta*other>0).astype('boolean').mask(d.beta.isna()|other.isna(),pd.NA)
        d['external_concordant'] = ((below(d.finngen_p)&agrees(d.finngen_beta))|(below(d.stuart_p)&agrees(d.stuart_beta))).mask(d.beta.isna(),pd.NA)
        d['followup_concordant'] = d.followup & d.external_concordant
        d = d.merge(old[['uniprot','p','passes_frozen_bonferroni','externally_supported']].rename(columns={'p':'original_p'}), on='uniprot', how='left')
        d['original_followup'] = d.passes_frozen_bonferroni.fillna(False)&d.externally_supported.fillna(False)
        comparisons.append(d)
    compare = pd.concat(comparisons, ignore_index=True); write(compare,'all_proteins_comparison.tsv')
    write(compare[compare.original_followup], 'original_17_comparison.tsv')
    write(compare[compare.followup.fillna(False)], 'corrected_followup.tsv')
    # Fixed original 17 x 22 clinical family; this is a comparison, not new-set enrichment.
    clinical_dir = ROOT/'research/data/processed/clinical_finngen'; panel = read(clinical_dir/'phenotype_manifest.tsv')
    ids = set(old.loc[old.passes_frozen_bonferroni&old.externally_supported,'uniprot'])
    clinical = []; clinical_audit = []; inputs += [clinical_dir/'phenotype_manifest.tsv']
    for code in panel.phenocode:
        p = clinical_dir/code/'regional_rows.tsv'; inputs.append(p); y = read(p)
        y['rsid'] = y.rsids.str.replace(';',',',regex=False).str.split(','); y = y.explode('rsid')
        y = y.rename(columns={'alt':'EA','ref':'OA','sebeta':'se','pval':'p_reported'})[['rsid','EA','OA','beta','se','p_reported']].drop_duplicates()
        y = y[~y.rsid.duplicated(keep=False)]
        clinical.append(fit_outcome(plan[plan.uniprot.isin(ids)],r,y,code,'reported_SE',{},clinical_audit,[]))
    clinical = pd.concat(clinical,ignore_index=True)
    clinical['BH_fixed_374'] = clinical.groupby('branch',sort=False).p.transform(lambda p:bh_fixed(p,374))
    write(clinical, 'clinical_original_17.tsv'); write(pd.concat(clinical_audit,ignore_index=True),'clinical_harmonization.tsv')
    inputs += [ROOT/'research/results/mr_fenland_substitute/all_proteins_comparison.tsv', INPUT/'selection.json',Path(__file__),Path(__file__).with_name('mr_covariance.py')]
    summary=[]
    for (branch,scenario), d in compare.groupby(['branch','scenario'],sort=True):
        summary.append(dict(branch=branch,scenario=scenario,estimable=int(d.p.notna().sum()),discovery=int(d.passes_fixed_bonferroni.sum()),
                            followup=int(d.followup.sum()),concordant_followup=int(d.followup_concordant.sum()),
                            original_17_retained=int((d.original_followup & d.followup.fillna(False)).sum())))
    record=dict(status='unique_variant_signed_LD_extension_computed',date=datetime.now(timezone.utc).isoformat(),
                primary_results_preserved=True,denominator=denominator,summary=summary,
                software=dict(Python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__),
                clinical=[dict(branch=b,computed=int(d.p.notna().sum()),BH05=int((d.BH_fixed_374<.05).sum())) for b,d in clinical.groupby('branch')],
                limitations=['Post hoc method extension, not recovered author implementation','NOME first-order gIVW and delta SE are separate',
                 'Reference LD and common protein effect scale are assumptions','Correlation grid is not measured sample overlap',
                 'MVP P-equivalent SE is only a sensitivity; source CI/P convention remains unconfirmed',
                 'Original 17 clinical family retained for comparison; downstream expression and enrichment not redefined as validation of a new list'],
                inputs={str(p.relative_to(ROOT)):sha(p) for p in inputs},script_sha256=sha(Path(__file__)),
                outputs={p.name:sha(p) for p in OUT.glob('*.tsv')})
    (OUT/'run.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps({k:record[k] for k in ('status','denominator','summary','clinical')},indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('--prepare', action='store_true'); a = p.parse_args()
    prepare() if a.prepare else analyse()
