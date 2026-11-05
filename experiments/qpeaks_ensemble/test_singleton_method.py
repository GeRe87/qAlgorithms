"""Conservative base-case correction; multi-model families remain frozen v1."""
import sys,unittest,math
from pathlib import Path
from bundle_method import topology,bundle_support

def reg(i,c,s,p,block=0):return dict(id=i,center=c,scale=s,block=block,**dict(zip(('b0','b1','b2','b3'),p)))
def singleton(rs,z):return topology(rs,model_boundaries=True,boundary_scope='singleton_blocks',log_values=z)
def group_members(bs,groups):return sorted(tuple(sorted(i for b in g['basins'] for i in bs[b]['candidate_ids'])) for g in groups)

class SingletonTests(unittest.TestCase):
    def test_single_available_model_can_form_a_bundle(self):
        p=[6.,.2,-.3,-.1];z=[p[0]+p[1]*(i-2)+(p[2] if i<2 else p[3])*(i-2)**2 for i in range(5)]
        rs=[reg(0,2,2,p)];bs,e,g,m=singleton(rs,z)
        self.assertEqual(len(g),1);self.assertEqual(g[0]['native_fits'],1);self.assertEqual(bs[g[0]['basins'][0]]['kind'],'boundary_model_maximum')
    def test_constant_and_monotone_singletons_are_not_forced_into_peaks(self):
        for p,z in (([4.,0.,-1e-14,-1e-14],[4.]*5),([4.,1.,0.,0.],[2.,3.,4.,5.,6.])):
            self.assertEqual(len(singleton([reg(0,2,2,p)],z)[2]),0)
    def test_all_multimodel_group_members_and_supports_are_unchanged(self):
        a=3.25;n=9;z=[6.-.1*(i-a)**2 for i in range(n)];rs=[]
        for s in range(2,5):
            for c in range(s,n-s):rs.append(reg(len(rs),c,s,[z[c],-.2*(c-a),-.1,-.1]))
        old=topology(rs);new=singleton(rs,z)
        self.assertEqual(group_members(old[0],old[2]),group_members(new[0],new[2]))
        self.assertEqual([bundle_support(g,old[0],rs) for g in old[2]],[bundle_support(g,new[0],rs) for g in new[2]])
        self.assertFalse(any(b['kind']=='boundary_model_maximum' for b in new[0]))
    def test_independent_blocks_do_not_change_each_other(self):
        z=[5.,5.5,6.,5.5,5.]+[math.nan]+[4.,5.,6.,6.,5.,4.]
        rs=[reg(0,2,2,[6.,.2,-.3,-.1]),reg(1,8,2,[6.,1.,-.5,-.5],1),reg(2,9,2,[6.,-1.,-.5,-.5],1)]
        old=topology(rs);new=singleton(rs,z)
        old_groups=[g for g in old[2] if old[0][g['basins'][0]]['block']==1];new_groups=[g for g in new[2] if new[0][g['basins'][0]]['block']==1]
        self.assertEqual(group_members(old[0],old_groups),group_members(new[0],new_groups))
        self.assertEqual(len(new[2]),len(old[2])+1)
    def test_two_model_boundary_case_remains_an_explicit_limit(self):
        # This version fixes the singleton base case, not all censored-center
        # geometry. Do not hide this remaining limitation with a fitted rule.
        z=[6.-.1*(i-.5)**2 for i in range(6)]
        rs=[reg(i,c,2,[z[c],-.2*(c-.5),-.1,-.1]) for i,c in enumerate((2,3))]
        self.assertEqual(len(topology(rs,model_boundaries=True,log_values=z)[2]),1)
        self.assertEqual(len(singleton(rs,z)[2]),0)

if __name__=='__main__':unittest.main()
