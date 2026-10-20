import gzip
import importlib.util
from pathlib import Path
import tempfile
import unittest
import pandas as pd

PATH=Path(__file__).resolve().parents[1]/'scripts/gene_identity.py'
spec=importlib.util.spec_from_file_location('identity',PATH);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)

class GeneIdentityTests(unittest.TestCase):
    def test_old_names_join_but_distinct_genes_do_not(self):
        with tempfile.TemporaryDirectory() as directory:
            p=Path(directory)/'genes.gz'
            with gzip.open(p,'wt') as f:
                f.write('#tax_id\tGeneID\tSymbol\tSynonyms\n9606\t1\tDNPH1\tRCL|AMBIG\n9606\t2\tMICA\tAMBIG\n9606\t3\tMICB\t-\n')
            locations=pd.DataFrame({'gene':['RCL','MICA','MICB'],'gene_id':['ENSG1','ENSG2','ENSG3']})
            resolve,_=module.identity_resolver(p,locations)
            self.assertEqual(resolve('DNPH1'),resolve('RCL'))
            self.assertNotEqual(resolve('MICA'),resolve('MICB'))
            self.assertIsNone(resolve('AMBIG'))
            self.assertIsNone(resolve('UNREPORTED'))

    def test_multiple_reference_locations_remain_unresolved(self):
        with tempfile.TemporaryDirectory() as directory:
            p=Path(directory)/'genes.gz'
            with gzip.open(p,'wt') as f:f.write('#tax_id\tGeneID\tSymbol\tSynonyms\n9606\t1\tGENE\tOLD\n')
            locations=pd.DataFrame({'gene':['GENE','OLD'],'gene_id':['ENSG1','ENSG2']})
            resolve,_=module.identity_resolver(p,locations)
            self.assertIsNone(resolve('GENE'))

if __name__=='__main__':unittest.main()
