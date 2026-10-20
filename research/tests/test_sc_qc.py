import unittest
import numpy as np
from research.scripts.sc_qc import qc_mask


class QCThresholdTests(unittest.TestCase):
    def test_inclusive_paper_boundaries(self):
        settings = dict(min_genes=500, max_genes=6000, min_counts=1000,
                        max_counts=30000, max_percent_mito=10)
        actual = qc_mask(np.array([500, 6000, 499, 6001, 500, 500, 500]),
                         np.array([1000, 30000, 1000, 30000, 999, 30001, 1000]),
                         np.array([10, 10, 0, 0, 0, 0, 10.01]), settings)
        np.testing.assert_array_equal(actual, [True, True, False, False, False, False, False])


if __name__ == "__main__":
    unittest.main()
