"""Fixed-support estimator identities; not tests tuned to legacy peak counts."""
import sys,unittest,math
from pathlib import Path
if len(sys.argv)>2 and sys.argv[1]=='--numpy-root':
    sys.path.insert(0,str(Path(sys.argv[2]).resolve()));del sys.argv[1:3]
import numpy as np
from common_profile import ProfileGeometry,summarize_fit

def asymmetric():
    x=np.arange(41)-20
    return np.exp(math.log(1000)+.1*x+np.where(x<0,-.025*x*x,-.08*x*x))

class ProfileTests(unittest.TestCase):
    def test_exact_asymmetric_mean_is_recovered(self):
        g=ProfileGeometry(41,np);r,data=summarize_fit(asymmetric(),g)
        self.assertEqual(r['status'],'interior_maximum');self.assertEqual(r['hinge'],20)
        self.assertAlmostEqual(r['peak_position'],20.625,places=10)
        self.assertTrue(np.allclose([r[k] for k in ('b0','b1','b2','b3')],[math.log(1000),.1,-.025,-.08],rtol=0,atol=1e-11))
        self.assertLess(r['rss'],1e-24)

    def test_amplitude_units_only_change_intercept(self):
        g=ProfileGeometry(41,np);reference,_=summarize_fit(asymmetric(),g)
        for factor in (1e-4,.01,10000):
            r,_=summarize_fit(asymmetric()*factor,g)
            self.assertEqual(r['hinge'],reference['hinge'])
            self.assertAlmostEqual(r['peak_position'],reference['peak_position'],places=10)
            self.assertAlmostEqual(r['b0']-reference['b0'],math.log(factor),places=11)
            for k in ('b1','b2','b3'):self.assertAlmostEqual(r[k],reference[k],places=11)

    def test_every_profile_matches_independent_svd_prediction(self):
        rng=np.random.default_rng(7);y=asymmetric()*np.exp(rng.normal(0,.1,41))
        g=ProfileGeometry(41,np);z,beta,rss,curves=g.fit(y)
        for j,h in enumerate(g.hinges):
            x=np.arange(41)-h;X=np.column_stack([np.ones(41),x,np.where(x<0,x*x,0),np.where(x>0,x*x,0)])
            b=np.linalg.lstsq(X,z,rcond=None)[0];prediction=X@b
            self.assertTrue(np.allclose(prediction,curves[j],rtol=0,atol=2e-10))
            self.assertLess(abs(np.sum((z-prediction)**2)-rss[j]),1e-10)
            self.assertTrue(np.allclose(X@beta[j],curves[j],rtol=0,atol=2e-10))

    def test_selected_fit_reflects_with_the_data(self):
        g=ProfileGeometry(41,np);rng=np.random.default_rng(21)
        for _ in range(10):
            y=asymmetric()*np.exp(rng.normal(0,.1,41));a,_=summarize_fit(y,g);b,_=summarize_fit(y[::-1],g)
            self.assertEqual(a['status'],b['status']);self.assertEqual(a['hinge']+b['hinge'],40)
            self.assertAlmostEqual(a['peak_position']+b['peak_position'],40,places=9)
            self.assertAlmostEqual(a['rss'],b['rss'],places=10)

    def test_constant_data_does_not_choose_a_hinge_or_peak(self):
        for n in (5,9,41,99):
            g=ProfileGeometry(n,np);r,data=summarize_fit(np.full(n,240.),g)
            self.assertEqual(r['status'],'constant_data');self.assertEqual(r['minimizers'],n-2)
            self.assertIsNone(r['hinge']);self.assertIsNone(r['peak_position']);self.assertEqual(r['rss'],0)
            self.assertTrue(np.array_equal(data[1][:,1:],np.zeros((n-2,3))))

    def test_geometry_is_reused_without_mutating_operators(self):
        g=ProfileGeometry(41,np);before=g.A.copy();q=g.Q.copy()
        summarize_fit(asymmetric(),g);summarize_fit(asymmetric()[::-1],g)
        self.assertTrue(np.array_equal(g.A,before));self.assertTrue(np.array_equal(g.Q,q))

    def test_unsupported_inputs_fail_explicitly(self):
        with self.assertRaises(ValueError):ProfileGeometry(4,np)
        g=ProfileGeometry(5,np)
        for y in ([1,2,0,4,5],[1,2,-1,4,5],[1,2,float('nan'),4,5],[1,2,float('inf'),4,5]):
            with self.assertRaises(ValueError):g.fit(y)

if __name__=='__main__':unittest.main()
