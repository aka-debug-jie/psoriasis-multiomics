import importlib.util
import unittest
from decimal import Decimal
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "audit_mvp_ci_precision.py"
SPEC = importlib.util.spec_from_file_location("audit_mvp_ci_precision", SCRIPT)
audit = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(audit)


class MvpPrecisionTest(unittest.TestCase):
    def test_extreme_tail_probability_is_not_cancelled_to_zero(self):
        from scipy.stats import norm
        row={'rsid':'rs-tail','odds_ratio':'1.19','ci_lower':'1.145','ci_upper':'1.237','p_value':'1e-18'}
        result=audit.audit_row(row,1,set())
        import math
        z=math.log(1.19)/((math.log(1.237)-math.log(1.145))/(2*norm.ppf(.975)))
        expected=2*norm.sf(abs(z));observed=float(result['p_ci_point_estimate'])
        self.assertGreater(observed,0)
        self.assertAlmostEqual(observed/expected,1,places=12)
    def test_rounding_precision_can_explain_coarse_fields(self):
        row = {
            "rsid": "rs-rounding",
            "odds_ratio": "1.10",
            "ci_lower": "1.00",
            "ci_upper": "1.21",
            "p_value": "0.05",
        }
        result = audit.audit_row(row, 1, {"rs-rounding"})
        self.assertEqual(result["status"], "consistent")
        self.assertEqual(result["currently_used_rsid"], "true")
        self.assertEqual(result["p_rounding_intersects_ci_feasible"], "true")

    def test_cd8a_source_precision_does_not_explain_p_value(self):
        row = {
            "rsid": "rs3020726",
            "odds_ratio": "1.051",
            "ci_lower": "1.017",
            "ci_upper": "1.085",
            "p_value": "0.001605",
        }
        result = audit.audit_row(row, 2, set())
        self.assertEqual(result["status"], "inconsistent")
        self.assertEqual(result["p_rounding_intersects_ci_feasible"], "false")

    def test_scientific_notation_and_na_are_preserved_and_classified(self):
        self.assertEqual(audit.half_ulp("1.25e-3"), Decimal("0.000005"))
        result = audit.audit_row({
            "rsid": "rs-na",
            "odds_ratio": "1.10e0",
            "ci_lower": "NA",
            "ci_upper": "1.20E+0",
            "p_value": "1.5e-2",
        }, 3, set())
        self.assertEqual(result["original_ci_lower"], "NA")
        self.assertEqual(result["status"], "unresolved")
        self.assertEqual(result["reason"], "missing_or_nonfinite_OR_or_CI")


if __name__ == "__main__":
    unittest.main()
