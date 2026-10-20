import sys, unittest
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import norm
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from mr_covariance import givw, deduplicate_sources, freeze_ld_subset, p_equivalent_se, bh_fixed


class CovarianceTests(unittest.TestCase):
    def setUp(self):
        self.x=np.array([.4,-.3]); self.y=np.array([.08,-.07])
        self.sx=np.array([.01,.02]); self.sy=np.array([.03,.04]); self.r=np.array([[1,.35],[.35,1]])

    def test_single_variant_matches_wald_and_delta(self):
        d=givw([.4],[.08],[.01],[.03],[[1]])
        self.assertAlmostEqual(d['beta'],.2)
        self.assertAlmostEqual(d['se'],.03/.4)
        self.assertAlmostEqual(d['se_delta']**2,.03**2/.4**2+.08**2*.01**2/.4**4)

    def test_duplicate_source_cannot_add_precision(self):
        d=pd.DataFrame(dict(uniprot=['P','P'],variant_key=['v','w'],F_statistic=[12,20],cohort=['B','A'],assay=['b','a'],rsid=['rs1','rs2']))
        a=deduplicate_sources(d).reset_index(drop=True)
        b=deduplicate_sources(pd.concat([d,d,d]).sample(frac=1,random_state=4)).reset_index(drop=True)
        pd.testing.assert_frame_equal(a,b)

    def test_allele_flip_requires_ld_sign_flip(self):
        a=givw(self.x,self.y,self.sx,self.sy,self.r)
        signs=np.array([-1,1]); b=givw(self.x*signs,self.y*signs,self.sx,self.sy,self.r*np.outer(signs,signs))
        for key in ('beta','se','se_delta','p','Q'): self.assertAlmostEqual(a[key],b[key],places=13)

    def test_independent_ld_matches_standard_ivw(self):
        a=givw(self.x,self.y,self.sx,self.sy,np.eye(2))
        w=1/self.sy**2; b=np.sum(w*self.x*self.y)/np.sum(w*self.x**2)
        self.assertAlmostEqual(a['beta'],b); self.assertAlmostEqual(a['se'],np.sqrt(1/np.sum(w*self.x**2)))

    def test_delta_gradient_against_finite_differences(self):
        a=givw(self.x,self.y,self.sx,self.sy,self.r)
        def f(x,y): return givw(x,y,self.sx,self.sy,self.r)['beta']
        eps=1e-6; gx=[];gy=[]
        for i in range(2):
            delta=np.eye(2)[i]*eps
            gx.append((f(self.x+delta,self.y)-f(self.x-delta,self.y))/(2*eps))
            gy.append((f(self.x,self.y+delta)-f(self.x,self.y-delta))/(2*eps))
        gx=np.array(gx);gy=np.array(gy)
        expected=gx@np.diag(self.sx**2)@gx+gy@(np.outer(self.sy,self.sy)*self.r)@gy
        self.assertAlmostEqual(a['se_delta']**2,expected,places=11)

    def test_missing_and_perfect_ld_do_not_use_pseudoinverse(self):
        self.assertEqual(freeze_ld_subset(np.ones((2,2)),.95)[0],[0])
        self.assertEqual(freeze_ld_subset([[1,np.nan],[np.nan,1]],.95),([0],'single_IV_missing_LD'))
        with self.assertRaises(np.linalg.LinAlgError): givw(self.x,self.y,self.sx,self.sy,np.ones((2,2)))

    def test_invalid_joint_covariance_rejected(self):
        with self.assertRaisesRegex(ValueError,'not_positive_semidefinite'):
            givw(self.x,self.y,self.sx,self.sy,self.r,rho_xx=0,rho_xy=1.)

    def test_fixed_family_and_missing_p_equivalence(self):
        np.testing.assert_allclose(bh_fixed([.001,np.nan,.02],4),[.004,np.nan,.04],equal_nan=True)
        for b,p in [(0,.1),(.2,0),(.2,1)]: self.assertTrue(np.isnan(p_equivalent_se(b,p)))
        se=p_equivalent_se(.2,.003)
        self.assertAlmostEqual(2*norm.sf(.2/se),.003)

    def test_gls_simulation_has_analytic_variance(self):
        rng=np.random.default_rng(20261017); s=np.outer(self.sy,self.sy)*self.r
        y=rng.multivariate_normal(self.x*.2,s,size=30000)
        w=np.linalg.solve(s,self.x); estimates=y@w/(self.x@w)
        expected=1/np.sqrt(self.x@w)
        self.assertLess(abs(np.std(estimates,ddof=1)/expected-1),.02)
        self.assertLess(abs(np.mean(np.abs(estimates-.2)<norm.ppf(.975)*expected)-.95),.005)


if __name__=='__main__': unittest.main()
