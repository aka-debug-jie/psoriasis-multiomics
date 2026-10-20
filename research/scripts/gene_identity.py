"""Resolve unique current/legacy symbols to a primary-assembly GENCODE gene."""
import csv
import gzip
from collections import defaultdict

def identity_resolver(gene_info,locations):
    symbols=defaultdict(set);aliases=defaultdict(set)
    with gzip.open(gene_info,'rt') as f:
        for r in csv.DictReader(f,delimiter='\t'):
            if r['#tax_id']!='9606':continue
            symbols[r['Symbol']].add(r['GeneID'])
            for a in r['Synonyms'].split('|'):
                if a!='-':aliases[a].add(r['GeneID'])
    def entrez(name):
        ids=symbols.get(name,set()) or aliases.get(name,set())
        return next(iter(ids)) if len(ids)==1 else None
    by_entrez=defaultdict(set)
    for _,r in locations.iterrows():
        eid=entrez(r.gene)
        if eid:by_entrez[eid].add(r.gene_id)
    def resolve(name):
        eid=entrez(name)
        ids=by_entrez.get(eid,set())
        return next(iter(ids)) if len(ids)==1 else None
    return resolve,entrez
