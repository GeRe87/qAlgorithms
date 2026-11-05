"""Bind all native records, reflection graphs, and fresh real EIC observations."""
import argparse,csv,json,sys,itertools,hashlib,re
from pathlib import Path
from collections import defaultdict
from real_comparison import read_eics

def rows(p):
    with p.open(encoding='utf-8',newline='') as f:yield from csv.DictReader(f)

def run(production_collection,centered_collection,production,centered,scenario_metadata,eic_binary,out):
    out.mkdir(parents=True,exist_ok=True);counts={}
    for label,source,target in (('production',production_collection,production),('centered',centered_collection,centered)):
        count=0
        for a,b in itertools.zip_longest(rows(source/'candidates.csv'),rows(target/'native-membership.csv')):
            assert a is not None and b is not None
            assert all(a[k]==b[k] for k in ('scenario','trial','id','block','center','scale','b0','b1','b2','b3'))
            count+=1
        report=json.loads((target/'report.json').read_text());assert count==report['native_fits_accounted']
        assert report['shared_fit_errors']==0 and not report['old_cascade_used_for_grouping_or_fits']
        counts[label]=count
    identical={name:(production/name).read_bytes()==(centered/name).read_bytes() for name in ('basins.csv','edges.csv','bundles.csv','hinge-fits.csv','fit-curves.csv','fits.json')}
    assert all(identical.values()),identical
    # Native original and centered arithmetic yield exactly the same slope-sign
    # topology here. This is a dataset check, not a universal sign guarantee.
    bs=defaultdict(dict);sizes={};edges=defaultdict(set)
    obs=defaultdict(list)
    for r in rows(centered_collection/'observations.csv'):obs[int(r['scenario'])].append(r)
    for s,r in obs.items():sizes[s]=len(r)
    for r in rows(production/'basins.csv'):bs[int(r['scenario'])][int(r['id'])]=r
    for r in rows(production/'edges.csv'):edges[int(r['scenario'])].add((int(r['from_basin']),int(r['to_basin']),int(r['shared_centers'])))
    def canonical(r,n=None):
        kind=r['kind'];a,b=int(r['first_center']),int(r['last_center'])
        root=() if not r['root_left'] else (int(r['root_left']),int(r['root_right']))
        if n is not None:
            a,b=n-1-b,n-1-a;root=tuple(n-1-v for v in root[::-1])
            kind={'left_boundary_exit':'right_boundary_exit','right_boundary_exit':'left_boundary_exit'}.get(kind,kind)
        return int(r['scale']),kind,a,b,root
    scenarios=json.loads(scenario_metadata.read_text(encoding='utf-8'))['scenarios']
    reflected=[]
    for scenario in scenarios:
        if 'reflection_of' not in scenario:continue
        s=scenario['reflection_of'];other=scenario['id']
        n=sizes[s];left={i:canonical(r) for i,r in bs[s].items()};right={i:canonical(r,n) for i,r in bs[other].items()}
        assert set(left.values())==set(right.values()),s
        a={(left[u],left[v],k) for u,v,k in edges[s]};b={(right[u],right[v],k) for u,v,k in edges[other]}
        assert a==b,s
        reflected.append(dict(scenario=s,basins=len(left),edges=len(a),exact_reflected_topology=True))
    eics=read_eics(eic_binary);bindings=[]
    for item in scenarios:
        match=re.fullmatch(r'real_EIC(\d+)_observations',item['name'])
        if not match:continue
        scenario=item['id'];binid=int(match.group(1));source=obs[scenario];native=eics[(1,binid)]
        assert len(source)==len(native['area'])
        for i,r in enumerate(source):
            assert float(r['y'])==native['area'][i] and float(r['x'])==native['rt'][i] and int(r['cumulative_df'])==native['df'][i]
        bindings.append(dict(scenario=scenario,polarity=1,binID=binid,observations=len(source),fresh_original_vectors_identical=True))
    report=dict(native_fit_records_verified=counts,original_centered_outputs_identical=identical,reflection=reflected,fresh_real_trace_bindings=bindings,real_eic_binary_sha256=hashlib.sha256(eic_binary.read_bytes()).hexdigest(),physical_peak_decisions=False)
    (out/'audit.json').write_text(json.dumps(report,indent=2),encoding='utf-8');return report

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('production_collection','centered_collection','production','centered','scenario_metadata','eic_binary','output'):p.add_argument(k,type=Path)
    a=p.parse_args();print(json.dumps(run(a.production_collection,a.centered_collection,a.production,a.centered,a.scenario_metadata,a.eic_binary,a.output),indent=2))
if __name__=='__main__':main()
