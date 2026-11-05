"""Independent structural/matrix checks for the exploratory ensemble analysis."""
import sys,unittest
from pathlib import Path
if len(sys.argv)>2 and sys.argv[1]=='--numpy-root':
    sys.path.insert(0,str(Path(sys.argv[2]).resolve()));del sys.argv[1:3]
import numpy as np
from analyse import lineage,operator,correlated_reference

def field(points):
    a=np.zeros((len(points),12))
    for i,(center,scale) in enumerate(points):a[i,0]=center;a[i,1]=scale;a[i,2]=1;a[i,3]=center+.2;a[i,6]=0
    return a

class AnalysisTests(unittest.TestCase):
    def test_lineage_keeps_two_parents_at_merge(self):
        a=field([(4,2),(6,2),(4,3),(5,3),(6,3)])
        runs,edges=lineage(a,np)
        self.assertEqual(len(runs),3);self.assertEqual(sum(r['parents']>1 for r in runs),1)
        self.assertEqual(sum(edges.values()),2)
        self.assertEqual(sorted(r['members'] for r in runs),[1,1,3])

    def test_lineage_keeps_two_children_at_split(self):
        a=field([(4,2),(5,2),(6,2),(4,3),(6,3)])
        runs,edges=lineage(a,np)
        self.assertEqual(len(runs),3);self.assertEqual(sum(r['children']>1 for r in runs),1)
        self.assertEqual(sum(edges.values()),2)

    def test_lineage_empty_and_reflected_fields(self):
        self.assertEqual(lineage(np.empty((0,12)),np),([],{}))
        a=field([(4,2),(6,2),(4,3),(5,3),(6,3)]);b=a.copy();b[:,0]=12-a[:,0];b[:,3]=12-a[:,3]
        ar,ae=lineage(a,np);br,be=lineage(b,np)
        self.assertEqual(sorted((r['scale'],r['members'],r['parents'],r['children']) for r in ar),sorted((r['scale'],r['members'],r['parents'],r['children']) for r in br))
        self.assertEqual(sorted(ae.values()),sorted(be.values()))

    def test_nonoverlapping_kernels_have_zero_cross_covariance(self):
        a=operator(15,3,2,np);b=operator(15,11,2,np);c=operator(15,4,2,np)
        self.assertTrue(np.array_equal(a@b.T,np.zeros((4,4))))
        self.assertGreater(np.linalg.norm(a@c.T),0)
        stacked=np.vstack([a,c]);cov=stacked@stacked.T
        self.assertTrue(np.allclose(cov[:4,4:],a@c.T,rtol=0,atol=1e-14))

    def test_correlated_ensemble_recovers_original_data_reference(self):
        for n in (9,21):
            r=correlated_reference(n,np)
            self.assertEqual(r['operator_rank'],n)
            self.assertGreater(r['stacked_covariance_nullity'],0)
            self.assertLess(r['GLS_vs_direct_OLS_operator_max_abs'],1e-10)
            self.assertTrue(all(v>1 for v in r['naive_actual_to_reported_standard_error_ratio']))

if __name__=='__main__':unittest.main()
