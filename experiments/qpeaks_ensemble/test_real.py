"""Real-input accounting, exact window binding, and native export contracts."""
import csv,struct,sys,tempfile,unittest,subprocess
from pathlib import Path
from real_comparison import read_eics,select_windows,rows

executable=None
if len(sys.argv)>2 and sys.argv[1]=='--exporter':
    executable=Path(sys.argv[2]).resolve();del sys.argv[1:3]

def record(ident=0,area=(2.,3.),df=(1,2)):
    n=len(area)
    return struct.pack('<3I',1,ident,n)+struct.pack('<'+'f'*n,*range(n))+struct.pack('<'+'f'*n,*area)+struct.pack('<'+'f'*n,*([100.]*n))+struct.pack('<'+'H'*n,*df)+struct.pack('<'+'I'*n,*range(n))

class RealTests(unittest.TestCase):
    def setUp(self):
        self.t=tempfile.TemporaryDirectory();self.folder=Path(self.t.name)
    def tearDown(self):self.t.cleanup()
    def load(self,data):
        p=self.folder/'eics.bin';p.write_bytes(data);return read_eics(p)
    def test_binary_exact_roundtrip(self):
        d=self.load(b'QEIC0001'+record())
        self.assertEqual(d[(1,0)]['area'],(2.,3.));self.assertEqual(d[(1,0)]['scan'],(0,1))
    def test_truncations_never_yield_partial_success(self):
        b=b'QEIC0001'+record()
        for end in range(len(b)):
            if end==8:continue # An empty polarity/blank file is a valid export.
            with self.assertRaises(ValueError):self.load(b[:end])
    def test_duplicate_and_invalid_observations_rejected(self):
        for payload in (record()+record(),record(area=(0.,3.)),record(area=(float('nan'),3.)),record(df=(1,1))):
            with self.assertRaises(ValueError):self.load(b'QEIC0001'+payload)
    def fixture(self):
        eics={}
        for p,label in ((1,'positive'),(2,'negative')):
            eics[(p,0)]={'rt':tuple(0.125*i for i in range(40))}
            fields=['featureID','binID','idxPeakStart','length','isUnresolved','lowerRT','upperRT']
            with (self.folder/f'x_{label}_area_timing_audit.csv').open('w',newline='') as f:
                w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
                for ident,(n,k) in enumerate(((9,0),(11,0),(10,1),(12,1))):
                    w.writerow(dict(featureID=ident,binID=0,idxPeakStart=ident,length=n,isUnresolved=k,lowerRT=.125*ident,upperRT=.125*(ident+n-1)))
            with (self.folder/f'x_{label}_features.csv').open('w',newline='') as f:
                w=csv.writer(f);w.writerow(['ID','binID','isUnresolved']);w.writerows((i,0,k) for i,k in enumerate((0,0,1,1)))
        return eics
    def test_declared_sampling_has_complete_binding_and_deduplicates(self):
        e=self.fixture();all_rows,chosen,provenance=select_windows(e,self.folder)
        self.assertEqual((len(all_rows),len(chosen),len(provenance)),(8,4,4))
        self.assertTrue(all(len(r['reasons'])==2 for r in chosen))
    def test_fresh_rt_mismatch_rejected(self):
        e=self.fixture();e[(1,0)]['rt']=tuple(v+1 for v in e[(1,0)]['rt'])
        with self.assertRaisesRegex(ValueError,'does not match'):select_windows(e,self.folder)
    def test_duplicate_audit_identity_rejected(self):
        e=self.fixture();p=self.folder/'x_positive_area_timing_audit.csv';s=p.read_text();lines=s.splitlines();lines[-1]=lines[-2];p.write_text('\n'.join(lines)+'\n')
        with self.assertRaisesRegex(ValueError,'one-to-one'):select_windows(e,self.folder)
    @unittest.skipUnless(executable,'Native exporter not supplied')
    def test_native_polarity_blank_and_unsupported_contracts(self):
        sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'test/integration'))
        from test_mzml_cli import write_fixture
        for name,kwargs,scans,success in (
            ('positive',{},[65,0],True),('mixed',{'polarity':'mixed'},[65,65],True),
            ('blank',{'blank':True},[8,0],True),('unsupported',{'ms_level':2},[0,0],False)):
            path=self.folder/(name+'.mzML');out=self.folder/name;write_fixture(path,**kwargs)
            r=subprocess.run([str(executable),str(path),str(out)],capture_output=True,text=True)
            self.assertEqual(r.returncode==0,success,r.stdout+r.stderr)
            if success:
                inventory=rows(out/'polarities.csv');self.assertEqual([int(x['selected_ms1_profile_scans']) for x in inventory],scans)
                e=read_eics(out/'eics.bin');self.assertEqual(len(e),sum(int(x['eics']) for x in inventory))
                self.assertEqual(bool(e),name!='blank')
            else:self.assertIn('no MS1 spectra with a known polarity',r.stderr)

if __name__=='__main__':unittest.main()
