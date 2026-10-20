from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
import numpy as np
import pandas as pd
from research.scripts.alternative_gwas_lookup import orient, compare


class AlternativeLookupTests(unittest.TestCase):
    def test_orientation_and_known_strand(self):
        self.assertEqual(orient("A", "C", "C", "A", -.2), (.2, "swapped"))
        self.assertEqual(orient("A", "C", "T", "G", .2), (.2, "strand_aligned"))
        self.assertTrue(np.isnan(orient("A", "C", "T", "G", .2, False)[0]))
        self.assertTrue(np.isnan(orient("A", "T", "A", "T", .2)[0]))

    def test_duplicate_effects_are_not_independent(self):
        baseline = pd.DataFrame([dict(rsid="rs1", EA="A", OA="C", beta=.1, se=.01, p=1e-8)])
        row = dict(rsid="rs1", chr="1", pos=100, EA="A", OA="C", beta=.2, se=.02, p=.01)
        with TemporaryDirectory() as d:
            result = compare(pd.DataFrame([row, row]), baseline, "identical", Path(d))
            self.assertEqual(result["matched_variants"], 1)
            result = compare(pd.DataFrame([row, dict(row, beta=-.2)]), baseline, "conflicting", Path(d))
            self.assertEqual(result["matched_variants"], 0)
            self.assertEqual(result["ambiguous_variants_excluded"], 1)


if __name__ == "__main__":
    unittest.main()
