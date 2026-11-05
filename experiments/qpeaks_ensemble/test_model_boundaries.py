"""Analytic boundary-root regressions; frozen v1 expectations stay unchanged."""
import unittest,sys,math
from pathlib import Path
if len(sys.argv)>2 and sys.argv[1]=='--numpy-root':
    sys.path.insert(0,str(Path(sys.argv[2]).resolve()));del sys.argv[1:3]
import numpy as np
from bundle_method import topology,bundle_support,MomentGeometry,geometric_peak

def native(ident,center,scale,coeffs):
    return dict(id=ident,block=0,center=center,scale=scale,**dict(zip(('b0','b1','b2','b3'),coeffs)))
def values(p,c,n):
    return [p[0]+p[1]*(i-c)+(p[2] if i<c else p[3])*(i-c)**2 for i in range(n)]

class BoundaryTests(unittest.TestCase):
    def test_single_native_regression_retained_and_fitted(self):
        p=[6.,.2,-.3,-.1];z=values(p,2,5);rs=[native(0,2,2,p)]
        self.assertEqual(len(topology(rs)[2]),0) # Frozen v1 comparator, not a new acceptance expectation.
        bs,edges,groups,membership=topology(rs,model_boundaries=True,log_values=z)
        self.assertEqual(len(groups),1);self.assertEqual(len(membership),1)
        self.assertEqual(bs[0]['kind'],'boundary_model_maximum')
        lo,hi,h=bundle_support(groups[0],bs,rs);self.assertEqual((lo,hi,h),(0,4,[2]))
        fits,_=MomentGeometry(np).fit(z,h);self.assertAlmostEqual(geometric_peak(fits[0]['coefficients'],2,5)['apex'],3.)
    def test_reflection_of_single_window(self):
        p=[6.,.2,-.3,-.1];q=[6.,-.2,-.1,-.3];z=values(p,2,5)
        a=topology([native(0,2,2,p)],model_boundaries=True,log_values=z)[0]
        b=topology([native(0,2,2,q)],model_boundaries=True,log_values=z[::-1])[0]
        self.assertAlmostEqual(a[0]['root'][0]+b[0]['root'][0],4.)
    def test_zero_center_slope_is_not_a_flat_signal(self):
        p=[6.,0.,-.3,-.1];rs=[native(0,2,2,p)]
        bs,_,g,_=topology(rs,model_boundaries=True,log_values=values(p,2,5))
        self.assertEqual(len(g),1);self.assertEqual(len(bs[0]['model_anchors']),1)
    def test_exact_constant_rejects_roundoff_curvature_without_epsilon(self):
        rs=[native(0,2,2,[4.,0.,-1e-14,-1e-14])]
        bs,_,g,m=topology(rs,model_boundaries=True,log_values=[4.]*5)
        self.assertEqual(len(g),0);self.assertEqual(len(m),1)
    def test_monotone_model_is_not_closed(self):
        p=[6.,1.,0.,0.];rs=[native(0,2,2,p)]
        self.assertEqual(len(topology(rs,model_boundaries=True,log_values=values(p,2,5))[2]),0)
    def test_outside_observed_apex_is_not_extrapolated(self):
        p=[6.,1.,-.1,-.1];rs=[native(0,2,2,p)]
        self.assertEqual(len(topology(rs,model_boundaries=True,log_values=values(p,2,5))[2]),0)
    def test_grid_edge_is_not_observed_window_edge(self):
        n=9;apex=.5;z=[6.-.1*(i-apex)**2 for i in range(n)]
        rs=[native(j,c,2,[z[c],-.2*(c-apex),-.1,-.1]) for j,c in enumerate(range(2,7))]
        bs,_,g,m=topology(rs,model_boundaries=True,log_values=z)
        self.assertEqual(len(g),1);self.assertEqual(len(m),len(rs));self.assertAlmostEqual(bs[g[0]['basins'][0]]['root'][0],apex)
    def test_separate_stationary_endpoint_models_do_not_get_one_owner(self):
        z=[0.,1.,0.,1.,0.,1.,0.]
        rs=[native(0,2,2,[24/35,0.,-1/7,-1/7]),native(1,3,2,[11/35,0.,1/7,1/7]),native(2,4,2,[24/35,0.,-1/7,-1/7])]
        bs,_,g,m=topology(rs,model_boundaries=True,log_values=z)
        self.assertEqual(len(g),2);self.assertEqual(sum(len(x['candidate_ids']) for x in bs),3)
        self.assertNotEqual(m[0],m[2]);self.assertNotEqual(m[1],m[0]);self.assertNotEqual(m[1],m[2])
    def test_all_scales_preserve_unambiguous_continuation(self):
        n=9;a=3.25;z=[6.-.1*(i-a)**2 for i in range(n)];rs=[]
        for s in range(2,5):
            for c in range(s,n-s):rs.append(native(len(rs),c,s,[z[c],-.2*(c-a),-.1,-.1]))
        bs,e,g,m=topology(rs,model_boundaries=True,log_values=z)
        self.assertEqual(len(g),1);self.assertEqual(len(g[0]['boundary_model_nodes']),1)
        self.assertEqual(len(m),len(rs));self.assertFalse(g[0]['contains_merge'])
    def test_missing_native_data_are_explicit_errors(self):
        rs=[native(0,2,2,[6.,1.,-1.,-1.])]
        with self.assertRaises(ValueError):topology(rs,model_boundaries=True)
        with self.assertRaises(ValueError):topology(rs,model_boundaries=True,log_values=[1.,2.,math.nan,2.,1.])
    def test_log_offset_does_not_change_boundary_geometry(self):
        p=[6.,.2,-.3,-.1];z=values(p,2,5)
        a=topology([native(0,2,2,p)],model_boundaries=True,log_values=z)
        q=[p[0]+10,*p[1:]];b=topology([native(0,2,2,q)],model_boundaries=True,log_values=[v+10 for v in z])
        self.assertEqual(a,b)
    def test_nonpositive_separator_is_not_bridged(self):
        p=[6.,.2,-.3,-.1];z=values(p,2,5)+[math.nan]+values(p,2,5)
        a=native(0,2,2,p);b=native(1,8,2,p);b['block']=1
        bs,e,g,m=topology([a,b],model_boundaries=True,log_values=z)
        self.assertEqual(len(g),2);self.assertEqual(e,{})
        for group in g:
            lo,hi,_=bundle_support(group,bs,[a,b]);self.assertFalse(lo<=5<=hi)

if __name__=='__main__':unittest.main()
