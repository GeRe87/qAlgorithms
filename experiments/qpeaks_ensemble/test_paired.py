"""Benchmark contract checks: exact grid accounting and localization semantics."""
import sys,unittest,tempfile,struct,itertools
from pathlib import Path
if len(sys.argv)>2 and sys.argv[1]=='--numpy-root':
    sys.path.insert(0,str(Path(sys.argv[2]).resolve()));del sys.argv[1:3]
from compare_bundle_cascade import fields,match_positions
from bundle_method import topology,bundle_support

class PairedTests(unittest.TestCase):
    def test_matching_is_one_to_one_with_extra_proposals(self):
        r=match_positions([1,4],[0,1.1,4.2,7])
        self.assertEqual(len(r['pairs']),2);self.assertEqual(r['unmatched_outputs'],2)
        self.assertAlmostEqual(r['squared_error'],.05)
    def test_far_match_is_not_claimed_detection_recall(self):
        r=match_positions([0],[1000]);self.assertEqual(r['squared_error'],1000000.)
        self.assertEqual(r['pairs'],[(0.,1000.)])
    def test_matching_under_detection_and_empty_sets(self):
        r=match_positions([1,4,8],[4.1]);self.assertEqual(r['unmatched_reference'],2)
        self.assertAlmostEqual(r['squared_error'],.01)
        self.assertEqual(match_positions([],[])['pairs'],[])
        self.assertEqual(match_positions([1],[])['unmatched_reference'],1)
    def test_matching_matches_exhaustive_small_subset_optimum(self):
        small=[0.,4.,8.];large=[-1.,.5,3.,6.,9.]
        expected=min(sum((a-b)**2 for a,b in zip(small,subset)) for subset in itertools.combinations(large,len(small)))
        self.assertEqual(match_positions(small,large)['squared_error'],expected)
        self.assertEqual(match_positions(large,small)['squared_error'],expected)
    def payload(self):
        return b'QFLD0001'+struct.pack('<4I',1,0,5,1)+struct.pack('<5f',1,2,3,2,1)+struct.pack('<3I',0,5,1)+struct.pack('<HH4d',2,2,3,0,-1,-1)
    def load(self,b):
        with tempfile.TemporaryDirectory() as t:
            p=Path(t)/'fields.bin';p.write_bytes(b);return list(fields(p))
    def test_native_binary_allocation_and_truncation_accounting(self):
        b=self.payload();r=self.load(b)
        self.assertEqual(r[0][0],(1,0));self.assertEqual(len(r[0][2]),1)
        for k in range(len(b)):
            if k==8:continue
            with self.assertRaises(ValueError):self.load(b[:k])
    def test_native_duplicate_identity_fails(self):
        b=self.payload()
        with self.assertRaises(ValueError):self.load(b+b[8:])
    def test_indexed_support_is_identical_to_original_scan(self):
        candidates=[dict(id=i,block=0,center=i+2,scale=2,b1=(1 if i<3 else -1)) for i in range(6)]
        bs,_,groups,_=topology(candidates);index={r['id']:r for r in candidates}
        self.assertEqual(bundle_support(groups[0],bs,candidates),bundle_support(groups[0],bs,candidates,index))

if __name__=='__main__':unittest.main()
