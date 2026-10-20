import math
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from research.scripts.prepare_mvp_smr import transform,deduplicate

class MVPToSMRTests(unittest.TestCase):
    def row(self,**kwargs):
        r=dict(rsid='rs1',effect_allele='A',other_allele='T',effect_allele_frequency='.00001',odds_ratio='2',standard_error='#NA',ci_lower='1',ci_upper='4',p_value='0',n='20')
        r.update(kwargs);return r
    def test_effect_CI_and_palindrome_preservation(self):
        v,origin,zero=transform(self.row())
        self.assertEqual(v[:3],['rs1','A','T']);self.assertTrue(zero)
        self.assertAlmostEqual(float(v[4]),math.log(2));self.assertAlmostEqual(float(v[5]),math.log(4)/(2*1.959963984540054))
        self.assertEqual(float(v[6]),1e-300);self.assertEqual(float(v[7]),20)
        v,origin,_=transform(self.row(standard_error='.2',ci_lower='#NA',ci_upper='#NA'))
        self.assertEqual(origin,'reported_SE');self.assertEqual(float(v[5]),.2)
        v,_,clamped=transform(self.row(p_value='1e-310'))
        self.assertEqual(float(v[6]),1e-310);self.assertFalse(clamped)
    def test_invalid_values_excluded(self):
        for changes in [dict(odds_ratio='0'),dict(effect_allele='AA'),dict(effect_allele_frequency='1'),dict(n='nan'),dict(ci_lower='5'),dict(p_value='1.1'),dict(rsid='chr1:20')]:
            self.assertIsNone(transform(self.row(**changes))[0])
    def test_conflicting_group_not_double_counted(self):
        with TemporaryDirectory() as d:
            p=Path(d);(p/'sorted').write_text('rs1\ta\nrs1\ta\nrs2\tb\nrs2\tb\nrs2\tc\nrs3\td\n')
            c=deduplicate(p/'sorted',p/'out',p/'conflicts')
            self.assertEqual(c['output_rows'],2);self.assertEqual(c['identical_duplicates_collapsed'],1)
            self.assertEqual(c['conflicting_rows_removed'],3);self.assertNotIn('rs2', (p/'out').read_text())

if __name__=='__main__':unittest.main()
