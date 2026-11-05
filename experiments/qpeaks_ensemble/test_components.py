"""Analytic model/Jacobian/nesting checks; no component-count thresholds."""
import sys,unittest,math
from pathlib import Path
while len(sys.argv)>2 and sys.argv[1] in ('--numpy-root','--scipy-root'):
    sys.path.insert(0,str(Path(sys.argv[2]).resolve()));del sys.argv[1:3]
import numpy as np
from scipy.optimize import lsq_linear,least_squares
from compare_components import Geometry,compare,fit_start

class ComponentTests(unittest.TestCase):
    def setUp(self):self.g=Geometry(29,np)

    def test_raw_addition_and_component_permutation(self):
        g=self.g;p=np.array([2.,1.,3.,4.,1.,-2.,4.,5.]);theta=p*np.r_[g.scale,g.scale]
        mean,jac=g.model_and_jacobian(theta)
        raw=np.exp(g.X@p[:4])+np.exp(g.X@p[4:])
        self.assertTrue(np.allclose(np.exp(mean),raw,rtol=1e-14,atol=0))
        swapped=np.r_[theta[4:],theta[:4]];other,j2=g.model_and_jacobian(swapped)
        self.assertTrue(np.array_equal(mean,other))
        self.assertTrue(np.allclose(jac,np.column_stack([j2[:,4:],j2[:,:4]]),rtol=0,atol=1e-15))

    def test_analytic_jacobian_matches_finite_differences(self):
        g=self.g;theta=np.array([2.,1.,3.,4.,1.,-2.,4.,5.])*np.r_[g.scale,g.scale]
        _,J=g.model_and_jacobian(theta)
        for k in range(8):
            step=np.finfo(float).eps**(1/3)*max(1.,abs(theta[k]));delta=np.zeros(8);delta[k]=step
            numerical=(g.model_and_jacobian(theta+delta)[0]-g.model_and_jacobian(theta-delta)[0])/(2*step)
            self.assertTrue(np.allclose(J[:,k],numerical,rtol=2e-7,atol=2e-10))

    def test_different_hinges_preserve_label_symmetry_and_jacobian(self):
        a=Geometry(29,np,8);b=Geometry(29,np,19)
        theta=np.array([2.,1.,3.,4.,1.,-2.,4.,5.])*np.r_[a.scale,b.scale]
        mean,J=a.model_and_jacobian(theta,b)
        swapped,other=b.model_and_jacobian(np.r_[theta[4:],theta[:4]],a)
        self.assertTrue(np.array_equal(mean,swapped))
        self.assertTrue(np.allclose(J,np.column_stack([other[:,4:],other[:,:4]]),rtol=0,atol=1e-15))
        for k in range(8):
            step=np.finfo(float).eps**(1/3)*max(1.,abs(theta[k]));delta=np.zeros(8);delta[k]=step
            numerical=(a.model_and_jacobian(theta+delta,b)[0]-a.model_and_jacobian(theta-delta,b)[0])/(2*step)
            self.assertTrue(np.allclose(J[:,k],numerical,rtol=2e-7,atol=2e-10))

    def test_gaussian_initial_curves_translate_between_hinges(self):
        a=Geometry(29,np,8);b=Geometry(29,np,19);x=np.arange(29)
        y=1000*np.exp(-(x-10)**2/12.5)+500*np.exp(-(x-18)**2/12.5)
        z,_=a.prepare(y);one=a.one(z,lsq_linear)
        seeds=dict(a.starts(z,one,b));theta=seeds['split_14']
        left=a.gaussian_seed(z,0,14);right=a.gaussian_seed(z,14,29)
        self.assertTrue(np.allclose(a.G@theta[:4],a.X@left,rtol=0,atol=1e-12))
        self.assertTrue(np.allclose(b.G@theta[4:],a.X@right,rtol=0,atol=1e-12))

    def test_duplicate_components_embed_single_and_have_rank_four(self):
        g=self.g;p=np.array([2.,.5,3.,4.]);half=p.copy();half[0]-=math.log(2)
        theta=np.r_[half*g.scale,half*g.scale];mean,J=g.model_and_jacobian(theta)
        self.assertTrue(np.allclose(mean,g.X@p,rtol=0,atol=1e-14))
        self.assertEqual(np.linalg.matrix_rank(J),4)

    def test_single_convex_fit_recovers_feasible_analytic_profile(self):
        g=self.g;truth=np.array([2.,1.,3.,4.]);z=g.X@truth
        r=g.one(z,lsq_linear)
        self.assertTrue(r['success']);self.assertLess(r['rss'],1e-25)
        self.assertTrue(np.allclose(r['theta']/g.scale,truth,rtol=0,atol=1e-12))

    def test_linear_bound_solution_satisfies_kkt(self):
        g=self.g;z=g.u*g.u
        r=g.one(z,lsq_linear);gradient=g.G.T@(r['curve']-z)
        self.assertTrue(r['success'])
        self.assertTrue(np.all(r['theta'][2:]>=0))
        self.assertLess(np.max(np.abs(gradient[:2])),1e-10)
        self.assertGreaterEqual(float(gradient[2:].min()),-1e-10)

    def test_convex_reference_matches_independent_active_set_enumeration(self):
        from audit_components import single_active_set_rss
        for h in (1,14,27):
            g=Geometry(29,np,h)
            for beta in ([2,1,3,4],[2,-1,-3,-4],[2,.2,0,0],[2,-.4,3,-2]):
                y=np.exp(g.X@np.array(beta));z,_=g.prepare(y);r=g.one(z,lsq_linear)
                self.assertTrue(r['success'])
                self.assertLess(abs(r['rss']-single_active_set_rss(y,h,np)),1e-10)

    def test_noiseless_two_gaussians_are_recovered_without_truth_starts(self):
        g=self.g;x=np.arange(g.n)
        y=1000*np.exp(-(x-10)**2/12.5)+500*np.exp(-(x-18)**2/12.5)
        r=compare(y,g,lsq_linear,least_squares)
        self.assertTrue(r['best']['success']);self.assertLess(r['best']['rss'],1e-16)
        apices=sorted(p['apex'] for p in r['components'])
        self.assertTrue(np.allclose(apices,[10,18],rtol=0,atol=1e-6))
        self.assertEqual(len(g.starts(r['z'],r['one'])),g.n-4)

    def test_component_export_ids_match_parameter_blocks(self):
        from types import SimpleNamespace
        g=self.g
        def gaussian(a,H):
            d=(g.hinge-a)/g.span;c=g.span**2/12.5
            return np.array([math.log(H)-c*d*d,-2*c*d,c,c])*g.scale
        # Deliberately store the right component first. Export IDs must retain
        # this order; sorting just the descriptions breaks curve-file identity.
        theta=np.r_[gaussian(18,500),gaussian(10,1000)]
        y=np.exp(g.model_and_jacobian(theta)[0])
        _,anchor=g.prepare(y)
        centered_theta=theta.copy();centered_theta[[0,4]]-=anchor*g.scale[0]
        def fixed_solver(fun,x0,**kwargs):
            return SimpleNamespace(x=centered_theta.copy(),success=True,status=1,nfev=1,optimality=0.,active_mask=np.zeros(8))
        r=compare(y,g,lsq_linear,fixed_solver)
        self.assertAlmostEqual(r['components'][0]['apex'],18,places=10)
        self.assertAlmostEqual(r['components'][1]['apex'],10,places=10)
        for j,p in enumerate(r['components']):
            delta=np.arange(g.n)-g.hinge
            predicted=p['log_intercept']+p['slope_per_sample']*delta+np.where(delta<0,p['left_curvature']*delta*delta,p['right_curvature']*delta*delta)
            self.assertTrue(np.allclose(predicted,g.G@r['best']['theta'][4*j:4*j+4]+r['anchor'],rtol=0,atol=1e-12))

    def test_constant_profile_does_not_imply_two_identifiable_components(self):
        for n in (9,29,99):
            g=Geometry(n,np);r=compare(np.full(n,240.),g,lsq_linear,least_squares)
            self.assertLessEqual(r['best']['rss'],1e-28)
            self.assertEqual(r['rank'],4)
            self.assertEqual(len(r['starts']),1)
            self.assertTrue(all(p['apex'] is None for p in r['components']))

    def test_intensity_scaling_only_shifts_log_predictions(self):
        g=self.g;x=np.arange(g.n);y=1000*np.exp(-(x-10)**2/12.5)+500*np.exp(-(x-18)**2/12.5)
        a=compare(y,g,lsq_linear,least_squares);b=compare(y*10000,g,lsq_linear,least_squares)
        self.assertLess(abs(a['best']['rss']-b['best']['rss']),1e-14)
        self.assertTrue(np.allclose(a['curve_two']+a['anchor']+math.log(10000),b['curve_two']+b['anchor'],rtol=0,atol=1e-7))

    def test_failed_optimizer_is_explicit_and_not_a_zero_fit(self):
        g=self.g;z,anchor=g.prepare(np.arange(1,g.n+1));one=g.one(z,lsq_linear);p=g.starts(z,one)[0][1]
        def fail(*args,**kwargs):raise RuntimeError('deliberate failure')
        r=fit_start(g,z,p,'failure',fail)
        self.assertFalse(r['success']);self.assertIsNone(r['rss']);self.assertIn('deliberate failure',r['error'])

    def test_unsupported_observations_and_saturated_reference(self):
        with self.assertRaises(ValueError):Geometry(8,np)
        for bad in (0.,-1.,float('nan'),float('inf')):
            y=np.ones(self.g.n);y[2]=bad
            with self.assertRaises(ValueError):self.g.prepare(y)

if __name__=='__main__':unittest.main()
