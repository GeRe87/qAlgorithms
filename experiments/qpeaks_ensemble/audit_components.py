"""Independently bind inputs, starts, parameters, curves, and hinge coverage."""
import argparse,csv,json,sys,math,hashlib,itertools
from pathlib import Path
from collections import Counter,defaultdict
from analyse import read_inputs

def rows(path):
    with path.open(encoding='utf-8',newline='') as f:yield from csv.DictReader(f)

def geometry(n,h,np):
    u=(np.arange(n)-h)/(n-1)
    X=np.column_stack([np.ones(n),u,-np.where(u<0,u*u,0),-np.where(u>0,u*u,0)])
    scale=np.sqrt(np.sum(X*X,axis=0))
    return X/scale,scale

def score(y,theta,h1,h2,np):
    n=len(y);z=np.log(y);anchor=z[n//2] if n%2 else (z[n//2-1]+z[n//2])/2;z=z-anchor
    G1,_=geometry(n,h1,np);G2,_=geometry(n,h2,np);a=G1@theta[:4];b=G2@theta[4:]
    mean=np.logaddexp(a,b);residual=mean-z
    return float(np.sum(residual*residual)),mean,anchor,z,a,b

def check_close(a,b,n,np):
    # Floating arithmetic audit only; not an inferential or component threshold.
    tolerance=256*np.finfo(float).eps*n*max(1.,abs(a),abs(b))
    assert abs(a-b)<=tolerance,(a,b,tolerance)

def single_active_set_rss(y,h,np):
    # Independent convex reference: enumerate all four possible active sets
    # of the two nonnegative curvature parameters, using SVD free solves.
    y=np.asarray(y);n=len(y);z=np.log(y)
    z=z-(z[n//2] if n%2 else (z[n//2-1]+z[n//2])/2)
    G,scale=geometry(n,h,np);candidates=[]
    for active in ((),(2,),(3,),(2,3)):
        free=[i for i in range(4) if i not in active]
        p=np.zeros(4);p[free]=np.linalg.lstsq(G[:,free],z,rcond=None)[0]
        if np.any(p[2:]<0):continue
        residual=G@p-z;candidates.append(float(residual@residual))
    assert candidates
    return min(candidates)

def run(common,comparison,grid,out,np):
    out.mkdir(parents=True,exist_ok=True)
    protocol=json.loads((comparison/'protocol.json').read_text());source=read_inputs(common/'common-inputs.bin')
    assert hashlib.sha256((common/'common-inputs.bin').read_bytes()).hexdigest()==protocol['source_common_input_sha256']
    supports=json.loads((common/'supports.json').read_text());ledger=list(rows(comparison/'selection-ledger.csv'))
    assert len(ledger)==len(supports)==len(source)==6425
    assert {int(r['source_case']) for r in ledger}=={r['case'] for r in supports}
    inputs=json.loads((comparison/'comparison-inputs.json').read_text());cases={c['id']:c for c in inputs}
    selected={int(r['source_case']) for r in ledger if r['selected']=='True'}
    assert selected=={c['source_case'] for c in inputs if c['source_case'] is not None and c['scenario']<7000}
    for c in inputs:
        y=np.asarray(c['y']);assert np.isfinite(y).all() and np.all(y>0) and len(y)>8
        if c['source_case'] is not None and c['scenario']<7000:assert tuple(c['y'])==source[(c['source_case'],0)][2]
        if c['scenario']==7000:
            parent=next(p for p in inputs if p['scenario']==7 and p['trial']==0)
            assert c['y']==parent['y'][::-1]
        if str(c['scenario']) in protocol['generated_control_definitions'] and c['scenario']>=7001:
            d=protocol['generated_control_definitions'][str(c['scenario'])];x=np.arange(d['n']);mean=np.full(d['n'],float(d.get('background',0)))
            for apex,H,width in zip(d['apices'],d['heights'],d['widths']):mean+=H*np.exp(-(x-apex)**2/(2*width**2))
            if d['log_noise_sd']:mean*=np.exp(np.random.default_rng(d['seed']).normal(0,d['log_noise_sd'],d['n']))
            assert np.array_equal(mean.astype(np.float32).astype(float),y)
    summary={int(r['case']):r for r in rows(comparison/'summary.csv')};assert set(summary)==set(cases) and len(cases)==109
    starts=defaultdict(list)
    for r in rows(comparison/'starts.csv'):starts[int(r['case'])].append(r)
    components=defaultdict(dict)
    for r in rows(comparison/'components.csv'):
        i,j=int(r['case']),int(r['component']);assert j not in components[i];components[i][j]=r
    curves=defaultdict(list)
    for r in rows(comparison/'curves.csv'):curves[int(r['case'])].append(r)
    max_score_error=0.;failures=0;checked_starts=0;max_kkt=0.;single_reference_max_error=0.;apex_errors=[];geometry_by_group=defaultdict(list)
    for case,c in cases.items():
        row=summary[case];n=len(c['y']);h=(n-1)/2;y=np.asarray(c['y']);G,scale=geometry(n,h,np)
        rs=starts[case];names=[r['start'] for r in rs];assert len(names)==len(set(names))==int(row['starts'])
        assert 'nested_single_candidate' in names
        if not np.all(y==y[0]):
            assert 'balanced_single_embedding' in names
            assert all('split_'+str(k) in names for k in range(3,n-2))
        else:assert names==['nested_single_candidate']
        for r in rs:
            if not r['rss']:
                assert r['success']=='False' and r['error'];failures+=1;continue
            theta=np.array([float(r['theta_'+str(i)]) for i in range(8)])
            assert np.isfinite(theta).all() and np.all(theta[[2,3,6,7]]>=0)
            rss,mean,anchor,z,a,b=score(y,theta,h,h,np);reported=float(r['rss']);check_close(rss,reported,n,np)
            max_score_error=max(max_score_error,abs(rss-reported));checked_starts+=1
            if r['success']=='False':failures+=1
        best=next(r for r in rs if r['start']==row['best_start']);finite=[float(r['rss']) for r in rs if r['rss']]
        assert float(best['rss'])==min(finite)==float(row['two_rss'])
        assert (best['success']=='True')==(row['status']=='best_found_converged')
        assert sum(r['success']=='False' for r in rs)==int(row['unsuccessful_starts'])
        theta=np.array([float(best['theta_'+str(i)]) for i in range(8)])
        rss,mean,anchor,z,a,b=score(y,theta,h,h,np)
        embedded=next(r for r in rs if r['start']=='nested_single_candidate')
        p=np.array([float(embedded['theta_'+str(i)]) for i in range(4)]);p[0]+=math.log(2)*scale[0]
        one=G@p;one_rss=float(np.sum((one-z)**2));check_close(one_rss,float(row['one_rss']),n,np)
        independent_one=single_active_set_rss(y,h,np);check_close(independent_one,one_rss,n,np);single_reference_max_error=max(single_reference_max_error,abs(independent_one-one_rss))
        grad=G.T@(one-z);kkt=max(float(np.max(np.abs(grad[:2]))),max((abs(grad[i]) if p[i]>0 else max(0.,-grad[i])) for i in (2,3)))
        assert kkt<=64*protocol['refinement_tolerance']*max(1.,float(np.linalg.norm(z)))
        max_kkt=max(max_kkt,kkt)
        assert float(row['two_rss'])<=float(row['one_rss'])+256*np.finfo(float).eps*n*max(1.,one_rss)
        assert set(components[case])=={0,1}
        positions=[]
        for j,r in components[case].items():
            p=theta[4*j:4*j+4]/scale
            expected=[p[0]+anchor,p[1]/(n-1),-p[2]/(n-1)**2,-p[3]/(n-1)**2]
            actual=[float(r[k]) for k in ('log_intercept','slope_per_sample','left_curvature','right_curvature')]
            assert np.allclose(actual,expected,rtol=0,atol=1e-12*max(1.,max(map(abs,expected))))
            if r['apex']:positions.append(float(r['apex']))
        geometry_by_group[c['scenario']].append(sum(r['interior_apex']=='True' for r in components[case].values()))
        if c['scenario']==7 and len(positions)==2:
            error=np.sort(positions)-np.sort(c['truth']);apex_errors.extend(map(float,error))
        if c['trial']==0:
            cr=curves[case];assert len(cr)==n and [int(r['index']) for r in cr]==list(range(n))
            for i,r in enumerate(cr):
                assert float(r['observed'])==y[i]
                for key,value in [('one_log_median',one[i]+anchor),('two_log_median',mean[i]+anchor),('component1_log_median',a[i]+anchor),('component2_log_median',b[i]+anchor)]:check_close(float(r[key]),float(value),n,np)
        else:assert case not in curves
    gridreport=json.loads((grid/'report.json').read_text());pairs=defaultdict(dict);gridstarts=defaultdict(list)
    for r in rows(grid/'pairs.csv'):
        key=int(r['case']);pair=float(r['hinge1']),float(r['hinge2']);assert pair not in pairs[key];pairs[key][pair]=r
    for r in rows(grid/'grid-starts.csv'):gridstarts[int(r['case'])].append(r)
    single_hinges=defaultdict(dict)
    for r in rows(grid/'single-hinges.csv'):
        i,h=int(r['case']),float(r['hinge']);assert h not in single_hinges[i];single_hinges[i][h]=r
    checked_grid=0
    for r in gridreport['cases']:
        c=cases[r['case']];y=np.asarray(c['y']);n=len(y);hs=r['hinges']
        expected=set(itertools.combinations_with_replacement(hs,2));assert set(pairs[c['id']])==expected and len(expected)==r['pairs_processed']==r['pairs_expected']
        assert set(single_hinges[c['id']])==set(hs)
        for h,q in single_hinges[c['id']].items():
            assert q['success']=='True';independent=single_active_set_rss(y,h,np);check_close(independent,float(q['rss']),n,np);single_reference_max_error=max(single_reference_max_error,abs(independent-float(q['rss'])))
        check_close(min(float(q['rss']) for q in single_hinges[c['id']].values()),r['profiled_one_rss'],n,np)
        rs=gridstarts[c['id']];assert len(rs)==r['start_candidates']+1
        assert sum(q['success']=='False' for q in rs if q['start']!='refine_grid_best')==r['unsuccessful_starts']
        for q in rs:
            if not q['rss']:assert q['success']=='False' and q['error'];continue
            theta=np.array([float(q['theta_'+str(i)]) for i in range(8)])
            value,*_=score(y,theta,float(q['hinge1']),float(q['hinge2']),np);check_close(value,float(q['rss']),n,np);checked_grid+=1
        check_close(r['fixed_one_rss'],float(summary[c['id']]['one_rss']),n,np);check_close(r['fixed_two_rss'],float(summary[c['id']]['two_rss']),n,np)
        assert r['best_found_two_rss']<=r['fixed_two_rss']+1e-12 and r['profiled_one_rss']<=r['fixed_one_rss']+1e-12
        best_recorded=min([float(q['rss']) for q in rs if q['rss']]+[r['fixed_two_rss']]);check_close(best_recorded,r['best_found_two_rss'],n,np)
    report=dict(compared_cases=len(cases),source_cases_accounted=len(ledger),source_selection_reasons=dict(Counter(r['reason'] for r in ledger)),start_records_checked=checked_starts,unsuccessful_starts_explicit=failures,maximum_independent_rss_difference=max_score_error,maximum_single_reference_kkt_residual=max_kkt,independent_bounded_single_max_rss_difference=single_reference_max_error,independent_bounded_single_profiles=109+sum(len(v) for v in single_hinges.values()),hinge_pairs_verified=sum(len(p) for p in pairs.values()),grid_start_records_checked=checked_grid,all_component_ids_and_exported_curves_bound=True,overlap_reference_component_apex_rmse=float(np.sqrt(np.mean(np.square(apex_errors)))) if apex_errors else None,overlap_apex_error_values=len(apex_errors),geometric_interior_apices_by_scenario={str(k):dict(Counter(v)) for k,v in geometry_by_group.items()},statistical_model_count_decision=False)
    (out/'audit.json').write_text(json.dumps(report,indent=2,allow_nan=False),encoding='utf-8');return report

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('common',type=Path);p.add_argument('comparison',type=Path);p.add_argument('grid',type=Path);p.add_argument('output',type=Path);p.add_argument('--numpy-root',type=Path);a=p.parse_args()
    if a.numpy_root:sys.path.insert(0,str(a.numpy_root.resolve()))
    import numpy as np
    print(json.dumps(run(a.common,a.comparison,a.grid,a.output,np),indent=2))
if __name__=='__main__':main()
