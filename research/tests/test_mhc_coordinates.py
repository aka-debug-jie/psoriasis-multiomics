import unittest
from research.scripts.audit_mr_selection import mhc_flags

class MHCCoordinateTests(unittest.TestCase):
    def test_boundary_spanning_gene_and_strand(self):
        self.assertEqual(mhc_flags('chr6',24999999,25000001,'+'),
                         dict(mhc_overlap=True,mhc_contained=False,mhc_tss=False))
        self.assertTrue(mhc_flags('chr6',24999999,25000001,'-')['mhc_tss'])
        self.assertFalse(any(mhc_flags('chr1',25000000,33000000,'+').values()))
        self.assertTrue(all(mhc_flags('chr6',25000000,33000000,'+').values()))
        self.assertFalse(any(mhc_flags('chr6',33000001,33000003,'+').values()))

if __name__ == '__main__':
    unittest.main()
