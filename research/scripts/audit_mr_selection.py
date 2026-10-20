"""Audit index substitutions, F strength and gene-based MHC location; no filtering."""
from pathlib import Path
from datetime import datetime, timezone
import gzip
import hashlib
import json
import re
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]


def mhc_flags(chromosome, start, end, strand):
    tss = start if strand == '+' else end
    return dict(mhc_overlap=chromosome == 'chr6' and start <= 33000000 and end >= 25000000,
                mhc_contained=chromosome == 'chr6' and start >= 25000000 and end <= 33000000,
                mhc_tss=chromosome == 'chr6' and 25000000 <= tss <= 33000000)


def gene_locations(path):
    rows = []
    with gzip.open(path, 'rt') as stream:
        for line in stream:
            if line.startswith('#'):
                continue
            columns = line.rstrip().split('\t')
            if columns[2] != 'gene':
                continue
            attrs = dict(re.findall(r'(\w+) "([^"]*)"', columns[8]))
            start, end = int(columns[3]), int(columns[4])
            rows.append(dict(gene=attrs['gene_name'], gene_id=attrs['gene_id'], gene_chr=columns[0],
                             gene_start=start, gene_end=end, strand=columns[6],
                             **mhc_flags(columns[0], start, end, columns[6])))
    return pd.DataFrame(rows)


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def main():
    inputs = [ROOT/'research/data/processed/exposures/eldjarn_all_cis_unselected.tsv',
              ROOT/'research/results/mr_available/cohort_wald.tsv',
              ROOT/'research/data/metadata/gencode.v19.annotation.gtf.gz']
    x, selected = [pd.read_csv(p, sep='\t') for p in inputs[:2]]
    x['exposure_row_id'] = np.arange(1, len(x)+1)
    x['F_statistic'] = (x.beta/x.se_reconstructed)**2
    group = ['cohort','uniprot']
    tops = x.groupby(group).neg_log10_p.max().rename('strongest_source_neg_log10_p')
    audit = selected.merge(tops, on=group, validate='many_to_one')
    audit['selected_weaker_than_source_max'] = audit.neg_log10_p < audit.strongest_source_neg_log10_p
    audit['F_independent'] = (audit['beta.exposure']/audit.se_reconstructed)**2
    if not np.allclose(audit.F_independent, audit.F_statistic, rtol=1e-8, atol=1e-10):
        raise ValueError('Stored F statistics do not agree with independent calculation')
    strongest = x.merge(tops, on=group, validate='many_to_one')
    strongest = strongest[strongest.neg_log10_p == strongest.strongest_source_neg_log10_p]
    strongest['selected_in_current_MR'] = strongest.exposure_row_id.isin(selected.exposure_row_id)
    locations = gene_locations(inputs[2])
    counts = locations.groupby('gene').size()
    locations['gene_mapping_status'] = np.where(locations.gene.map(counts) == 1, 'unique_symbol', 'multiple_gene_locations')
    mhc = x[['cohort','uniprot','gene']].drop_duplicates().merge(locations, on='gene', how='left')
    mhc['gene_mapping_status'] = mhc.gene_mapping_status.fillna('unmapped_symbol')
    # Neither variant positions nor ambiguous/missing gene mapping define exclusion.
    mhc['used_to_filter_current_MR'] = False
    out = ROOT/'research/results/mr_selection_audit'; out.mkdir(parents=True, exist_ok=True)
    tables = {'selected_instrument_diagnostics.tsv':audit, 'strongest_source_rows.tsv':strongest,
              'gene_MHC_locations_GRCh37.tsv':mhc}
    for name, frame in tables.items():
        frame.to_csv(out/name, sep='\t', index=False)
    unique_selected = mhc.merge(selected[['cohort','uniprot']].drop_duplicates(),on=group)
    mhc_proteins = unique_selected.loc[(unique_selected.gene_mapping_status == 'unique_symbol') &
                                      (unique_selected.mhc_overlap == True), 'uniprot'].nunique()
    report = dict(status='diagnostics_completed_no_formal_filtering', cohort_estimates=len(audit),
                  weaker_than_original_strongest=int(audit.selected_weaker_than_source_max.sum()),
                  F_minimum=float(audit.F_independent.min()), F_at_or_below_10=int((audit.F_independent<=10).sum()),
                  selected_proteins_with_unique_gene_MHC_overlap=int(mhc_proteins),
                  source_proteins=int(x.uniprot.nunique()), selected_proteins=int(selected.uniprot.nunique()),
                  paper_denominators=[2598,2647], paper_reported_p=1.89e-5,
                  observed_coverage_bonferroni_p=.05/selected.uniprot.nunique(),
                  formal_denominator_frozen=False, steiger_filter_applied=False,
                  notes=['F uses reconstructed exposure SE', 'GENCODE v19 GRCh37 positions are gene intervals, not instrument coordinates',
                         'MHC overlap, containment and TSS flags retained separately; no new exclusions',
                         'Index replacement diagnosis does not imply author selection rule recovered'],
                  finished_utc=datetime.now(timezone.utc).isoformat(),
                  inputs={str(p.relative_to(ROOT)):digest(p) for p in inputs},
                  script_sha256=digest(Path(__file__)),outputs={name:digest(out/name) for name in tables})
    (out/'run.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ('inputs','outputs','notes')},indent=2))


if __name__ == '__main__':
    main()
