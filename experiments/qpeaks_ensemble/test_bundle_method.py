"""Structural and algebraic checks of the proposed two-stage reference."""
import sys,unittest
from pathlib import Path
if len(sys.argv)>2 and sys.argv[1]=='--numpy-root':
    sys.path.insert(0,str(Path(sys.argv[2]).resolve()));del sys.argv[1:3]
import numpy as np
from bundle_method import basins,topology,MomentGeometry,bundle_support,geometric_peak

def field(slopes,scale=2,block=0,start=2,offset=0):
    return [dict(id=offset+j,center=start+j,scale=scale,block=block,b1=d) for j,d in enumerate(slopes)]

class BundleMethodTests(unittest.TestCase):
    def test_maximum_corridor_and_complete_assignment(self):
        r,owners=basins(list(range(7)),[2,1,-1,-2,1,2,3])
        self.assertEqual([(x['root'],x['positions']) for x in r if x['kind']=='maximum'],[((1,2),[0,1,2,3])])
        self.assertEqual(r[owners[-1]]['kind'],'right_boundary_exit');self.assertEqual(len(owners),7)
    def test_flat_plateau_and_minimum_are_distinct(self):
        r,o=basins(list(range(7)),[1,0,0,-1,-1,0,1])
        self.assertEqual([x['root'] for x in r if x['kind']=='maximum'],[(1,2)])
        self.assertEqual(r[o[5]]['kind'],'stationary_minimum')
        r,o=basins(list(range(5)),[0]*5);self.assertFalse(any(x['kind']=='maximum' for x in r))
    def test_same_sign_stationary_shoulder_continues(self):
        r,o=basins(list(range(6)),[1,0,0,1,-1,-2])
        self.assertEqual(len(r),1);self.assertEqual(r[0]['positions'],list(range(6)))
    def test_reflection_including_stationary_runs(self):
        d=[1,0,1,-1,-2,0,2,3,-1];n=len(d)
        a,_=basins(list(range(n)),d);b,_=basins(list(range(n)),[-x for x in d[::-1]])
        canon=lambda rs:sorted((tuple(r['root']),tuple(r['positions'])) for r in rs if r['kind']=='maximum')
        reflected=[dict(r,root=None if r['root'] is None else tuple(n-1-v for v in r['root'][::-1]),positions=sorted(n-1-v for v in r['positions'])) for r in b]
        self.assertEqual(canon(a),canon(reflected))
    def test_scale_merge_retained_and_no_legacy_annotation_needed(self):
        a=field([1,-1,-1,1,1,-1,-1]);b=field([1,1,-1,-1,-1],scale=3,start=3,offset=len(a))
        r,e,g,m=topology(a+b)
        self.assertEqual(len(m),12);self.assertEqual(len(e),2);self.assertEqual(len(g),3)
        self.assertEqual(sum(x['contains_merge'] for x in g),1)
        self.assertTrue(all(x['native_fits']>0 for x in g))
    def test_diagonal_root_motion_does_not_require_equal_root_centers(self):
        a=field([1,-1,-1,-1,-1]);b=field([1,1,-1],scale=3,start=3,offset=len(a))
        r,e,g,m=topology(a+b);self.assertEqual(len(g),1)
        self.assertNotEqual(r[g[0]['basins'][0]]['root'],r[g[0]['basins'][1]]['root'])
    def test_prefix_fit_matches_direct_svd_for_every_hinge(self):
        z=np.random.default_rng(451).normal(size=31);geometry=MomentGeometry(np)
        fits,winners=geometry.fit(z,list(range(1,30)))
        for r in fits:
            h=r['hinge'];t=np.arange(len(z))-h;X=np.column_stack([np.ones(len(z)),t,np.where(t<0,t*t,0),np.where(t>0,t*t,0)])
            p=np.linalg.lstsq(X,z,rcond=None)[0];prediction=X@p;got=X@np.asarray(r['coefficients'])
            self.assertTrue(np.allclose(got,prediction,rtol=0,atol=2e-10))
            self.assertAlmostEqual(r['rss'],float(np.sum((prediction-z)**2)),places=11)
        self.assertEqual(len(geometry.cache),29)
    def test_ensemble_support_and_exact_piecewise_recovery(self):
        a=field([1,1,-1,-1,-1],scale=2);r,e,g,m=topology(a)
        left,right,hinges=bundle_support(g[0],r,a);self.assertEqual((left,right),(0,8))
        x=np.arange(9);h=4;t=x-h;z=6+.2*t+np.where(t<0,-.02,-.08)*t*t
        fits,winners=MomentGeometry(np).fit(z,hinges)
        best=min(fits,key=lambda r:r['rss']);self.assertEqual(best['hinge'],h)
        self.assertLess(best['rss'],1e-24)
    def test_log_offset_preserves_topology_and_changes_only_fit_intercept(self):
        z=np.array([1.,2.,3.,4.,3.,2.,1.]);geom=MomentGeometry(np)
        a,_=geom.fit(z,[2,3,4]);b,_=geom.fit(z+8,[2,3,4])
        for x,y in zip(a,b):
            self.assertTrue(np.allclose(np.asarray(y['coefficients'])-x['coefficients'],[8,0,0,0],atol=1e-12,rtol=0))
            self.assertAlmostEqual(x['rss'],y['rss'],places=12)
    def test_hinge_is_not_assumed_to_be_apex(self):
        r=geometric_peak([6.,.2,-.02,-.08],4,9)
        self.assertAlmostEqual(r['apex'],5.25)
        self.assertTrue(r['above_both_support_edges'])
        self.assertIsNone(geometric_peak([6.,0.,0.,0.],4,9)['apex'])
    def test_invalid_and_duplicate_inputs_are_explicit_errors(self):
        with self.assertRaises(ValueError):basins([0,2],[1.,-1.])
        with self.assertRaises(ValueError):basins([0],[float('nan')])
        a=field([1,-1]);a[1]['id']=a[0]['id']
        with self.assertRaises(ValueError):topology(a)

if __name__=='__main__':unittest.main()
