"""Structural tests for the experimental collector; no new peak decision rule."""
import csv,math,struct,subprocess,tempfile,unittest,sys
from pathlib import Path
EXE=Path(sys.argv.pop(1)).resolve() if len(sys.argv)>1 else None

def write_inputs(path,records):
    with path.open('wb') as f:
        f.write(b'QENS0001')
        for scenario,x,y in records:
            n=len(y);f.write(struct.pack('<4I',scenario,0,n,0))
            f.write(struct.pack('<'+'f'*n,*x));f.write(struct.pack('<'+'f'*n,*y))
            f.write(struct.pack('<'+'H'*n,*range(1,n+1)))

def table(folder,name):
    with (folder/name).open(newline='',encoding='utf-8') as f:return list(csv.DictReader(f))

class CollectorTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if EXE is None:raise RuntimeError('Pass collector executable as first argument.')
        cls.temp=tempfile.TemporaryDirectory(prefix='qpeaks-ensemble-')
        cls.root=Path(cls.temp.name);cls.out=cls.root/'out';cls.input=cls.root/'input.bin'
        y=[math.exp(3-0.02*(i-4.3)**2) for i in range(21)]
        split=[0,0,1,2,3,4,5,0,6,7,8,0]
        cls.records=[(0,list(range(21)),y),(1,list(range(21)),y[::-1]),
                     (2,[i/60 for i in range(21)],y),
                     (3,list(range(12)),split),(4,list(range(9)),[240]*9)]
        write_inputs(cls.input,cls.records)
        p=subprocess.run([str(EXE),str(cls.input),str(cls.out)],capture_output=True,text=True)
        if p.returncode:raise AssertionError(p.stdout+p.stderr)
        cls.traces=table(cls.out,'traces.csv');cls.candidates=table(cls.out,'candidates.csv')
        cls.edges=table(cls.out,'edges.csv');cls.components=table(cls.out,'components.csv')

    @classmethod
    def tearDownClass(cls):cls.temp.cleanup()

    def test_complete_enumeration_and_point_accounting(self):
        rows={int(r['scenario']):r for r in self.traces}
        self.assertEqual(len(rows),5)
        for scenario in (0,1,2):
            self.assertEqual(int(rows[scenario]['all_regressions']),sum(21-2*s for s in range(2,11)))
            self.assertEqual(int(rows[scenario]['positive_blocks']),1)
        self.assertEqual(int(rows[3]['all_regressions']),1)
        self.assertEqual(int(rows[3]['short_positive_samples']),3)
        self.assertEqual(int(rows[3]['nonpositive_samples']),4)
        for r in self.traces:
            rs=[c for c in self.candidates if c['scenario']==r['scenario']]
            self.assertEqual(len(rs),int(r['all_regressions']))
            self.assertEqual([int(c['id']) for c in rs],list(range(len(rs))))
        self.assertTrue(any(r['all_above_production_floor']=='0' for r in self.candidates if r['scenario']=='0'))

    def test_edges_and_components_match_independent_grid_bfs(self):
        for trace in self.traces:
            scenario=trace['scenario'];rs=[r for r in self.candidates if r['scenario']==scenario]
            grid={(int(r['center']),int(r['scale'])):int(r['id']) for r in rs if r['has_interior_maximum']=='1'}
            expected=set();neighbors={i:set() for i in grid.values()}
            for (center,scale),i in grid.items():
                for axis,key in [('center',(center+1,scale)),('scale',(center,scale+1))]:
                    if key in grid:
                        j=grid[key];expected.add((i,j,axis));neighbors[i].add(j);neighbors[j].add(i)
            actual={(int(r['from_id']),int(r['to_id']),r['axis']) for r in self.edges if r['scenario']==scenario}
            self.assertEqual(expected,actual)
            unseen=set(neighbors);groups=[]
            while unseen:
                work=[min(unseen)];group=set()
                while work:
                    i=work.pop()
                    if i in group:continue
                    group.add(i);unseen.discard(i);work.extend(neighbors[i]-group)
                groups.append(group)
            labels={int(r['id']):int(r['component']) for r in rs if r['has_interior_maximum']=='1'}
            self.assertEqual(len(groups),int(trace['graph_components']))
            for group in groups:self.assertEqual({labels[i] for i in group},{min(group)})
            self.assertEqual(sum(int(r['nodes']) for r in self.components if r['scenario']==scenario),len(grid))

    def test_axis_units_do_not_change_candidate_field(self):
        a=[r for r in self.candidates if r['scenario']=='0'];b=[r for r in self.candidates if r['scenario']=='2']
        for x,y in zip(a,b):
            self.assertEqual({k:v for k,v in x.items() if k!='scenario'},{k:v for k,v in y.items() if k!='scenario'})

    def test_reflection_of_nonboundary_gaussian(self):
        a={(int(r['center']),int(r['scale'])):r for r in self.candidates if r['scenario']=='0'}
        b={(20-int(r['center']),int(r['scale'])):r for r in self.candidates if r['scenario']=='1'}
        self.assertEqual(set(a),set(b))
        for key,x in a.items():
            y=b[key];self.assertEqual(x['has_interior_maximum'],y['has_interior_maximum'])
            if x['has_interior_maximum']=='1':
                # The existing kernel starts several sums in float; bound only
                # this diagnostic coordinate comparison at float resolution.
                self.assertLess(abs(float(x['apex'])+float(y['apex'])-20),128*2**-23*21)

    def test_exact_flatness_is_recorded_not_hidden_by_a_tolerance(self):
        flat=[r for r in self.candidates if r['scenario']=='4']
        self.assertTrue(flat);self.assertTrue(all(r['constant_log_window']=='1' for r in flat))
        self.assertTrue(all(math.isfinite(float(r[k])) for r in flat for k in ('b0','b1','b2','b3','full_rss')))
        # Do not assert that floating arithmetic produces zero stationary points.
        # Their presence is exactly the degeneracy this collector must expose.

    def test_truncation_and_nonfinite_inputs_fail_explicitly(self):
        bad=self.root/'truncated.bin';bad.write_bytes(self.input.read_bytes()[:-1])
        p=subprocess.run([str(EXE),str(bad),str(self.root/'truncated')],capture_output=True,text=True)
        self.assertNotEqual(p.returncode,0);self.assertIn('Truncated',p.stderr)
        invalid=self.root/'invalid.bin';write_inputs(invalid,[(9,list(range(5)),[1,2,float('nan'),4,5])])
        p=subprocess.run([str(EXE),str(invalid),str(self.root/'invalid')],capture_output=True,text=True)
        self.assertNotEqual(p.returncode,0);self.assertIn('finite',p.stderr)

if __name__=='__main__':unittest.main()
