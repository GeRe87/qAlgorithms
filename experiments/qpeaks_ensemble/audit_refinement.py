"""Integrity and independent arithmetic checks for the second experiment."""
import argparse,csv,json,hashlib,sys,re,math
from pathlib import Path
from collections import Counter,defaultdict
from analyse import read_inputs,design

def rows(path):
    with path.open(encoding='utf-8',newline='') as f:yield from csv.DictReader(f)

def run(study,refinement,header,np):
    original={(int(r['scenario']),int(r['trial'])):r for r in rows(study/'release/traces.csv')}
    current={(int(r['scenario']),int(r['trial'])):r for r in rows(refinement/'centered-release/traces.csv')}
    assert set(current)==set(original) and len(current)==5220
    changed=[]
    for k,r in current.items():
        for f in ('n','truth_count','positive_blocks','short_positive_samples','nonpositive_samples','all_regressions','constant_log_windows','current_accepted_windows','strict_output_count','strict_error'):
            assert r[f]==original[k][f],(k,f)
        if r['interior_maxima']!=original[k]['interior_maxima']:changed.append(dict(scenario=k[0],trial=k[1],before=int(original[k]['interior_maxima']),centered=int(r['interior_maxima'])))
    assert current[(3000,0)]['interior_maxima']=='0' and current[(3000,0)]['graph_components']=='0'
    build=[]
    for p in sorted((refinement/'centered-release').glob('*.csv')):
        q=refinement/'centered-debug-checked'/p.name
        a=hashlib.sha256(p.read_bytes()).hexdigest();b=hashlib.sha256(q.read_bytes()).hexdigest();assert a==b
        build.append(dict(file=p.name,sha256=a,bytes=p.stat().st_size,debug_checked_identical=True))
    # Reconstruct the exact float log observations exported by the native code.
    logs=defaultdict(list);input_records=read_inputs(study/'inputs.bin')
    for r in rows(refinement/'centered-release/observations.csv'):
        k=int(r['scenario']),int(r['trial']);i=int(r['index']);assert i==len(logs[k])
        logs[k].append(float(r['log_y_float']) if r['log_y_float'] else np.nan)
    blocks={}
    for k,z in logs.items():
        y=input_records[k][2];segments=[];start=0
        while start<len(y):
            while start<len(y) and y[start]<=0:start+=1
            end=start
            while end<len(y) and y[end]>0:end+=1
            if end-start>=5:segments.append((start,end))
            if end==len(y):break
            start=end+1
        blocks[k]=segments
        assert len(z)==len(y)
    content=header.read_text(encoding='utf-8');body=content.split('const MatInverse qalgo_matInverse[] = {',1)[1].split('};',1)[0]
    values=[[float(v) for v in m.split(',')] for m in re.findall(r'\{([^{}]+)\}',body)]
    mats=[]
    for A,B,C,D,E,F in values:mats.append(np.array([[A,0,B,B],[0,C,D,-D],[B,D,E,F],[B,-D,F,E]]))
    detailed_counts=Counter();native_normalized_bound_max=0.;qr_discrepancy_max=0.;verified=0;selected=[]
    eps=np.finfo(float).eps
    # All detailed synthetic/control windows; first/last center per scale in
    # each real positive block. This samples arithmetic, not statistical data.
    for r in rows(refinement/'centered-release/candidates.csv'):
        k=int(r['scenario']),int(r['trial']);detailed_counts[k]+=1
        center=int(r['center']);scale=int(r['scale']);block=int(r['block']);start,end=blocks[k][block]
        if 2000<=k[0]<3000 and center not in (start+scale,end-scale-1):continue
        z=np.asarray(logs[k]);window=z[center-scale:center+scale+1]
        length=end-start;anchor=z[start+length//2] if length%2 else (z[start+length//2-1]+z[start+length//2])/2
        centered=window-anchor;X=design(scale,np);M=mats[scale]
        moment=X.T@centered;predicted=M@moment;predicted[0]+=anchor
        actual=np.array([float(r[f'b{i}']) for i in range(4)])
        n=len(window);gamma=2*n*eps/(1-2*n*eps);gamma8=8*eps/(1-8*eps)
        bound=2*np.abs(M)@(gamma*(np.abs(X).T@np.abs(centered)))+2*gamma8*np.abs(M)@np.abs(moment)
        bound+=eps*(np.abs(predicted)+np.array([abs(anchor),0,0,0]))
        difference=np.abs(actual-predicted)
        assert np.all(difference<=bound),(k,center,scale,difference,bound)
        relative=np.divide(difference,bound,out=np.zeros(4),where=bound>0)
        native_normalized_bound_max=max(native_normalized_bound_max,float(relative.max()))
        Q,R=np.linalg.qr(X,mode='reduced');independent=np.linalg.solve(R,Q.T@centered);independent[0]+=anchor
        discrepancy=float(np.max(np.abs(actual-independent)));qr_discrepancy_max=max(qr_discrepancy_max,discrepancy)
        selected.append(dict(scenario=k[0],trial=k[1],center=center,scale=scale,native_vs_direct_table_max=float(difference.max()),rounding_bound_max=float(bound.max()),native_vs_QR_max=discrepancy))
        verified+=1
    for k,r in current.items():
        if r['detail']=='1':assert detailed_counts[k]==int(r['all_regressions'])
    # Complete binding of every common-profile result to its serialized support.
    common=refinement/'common-profile';inputs=read_inputs(common/'common-inputs.bin');supports=json.loads((common/'supports.json').read_text())
    results=list(rows(common/'case-results.csv'));assert len(inputs)==len(supports)==len(results)==6425
    expected_profiles={};expected_curves={}
    for info,r in zip(supports,results):
        case=info['case'];assert int(r['case'])==case
        truth,x,y=inputs[(case,0)];assert len(y)==int(r['n'])==info['end']-info['start']+1
        if info['scenario']<4000:
            original_key=info['scenario'],info['trial'];source=input_records[original_key]
            assert tuple(y)==source[2][info['start']:info['end']+1]
        count=sum(v<=0 for v in y);assert count==int(r['nonpositive'])
        if count:assert r['status']=='nonpositive_observations' and not r['peak_position']
        elif len(y)<5:assert r['status']=='too_short'
        elif info['trial']==0:
            expected_profiles[case]=len(y)-2
            if r['status']!='multiple_equal_minimizers':expected_curves[case]=len(y)
    actual_profiles=Counter(int(r['case']) for r in rows(common/'hinge-profiles.csv'))
    actual_curves=Counter(int(r['case']) for r in rows(common/'fit-curves.csv'))
    assert actual_profiles==expected_profiles and actual_curves==expected_curves
    report=dict(records=5220,all_production_comparisons_unchanged=True,build_comparison=build,maximum_classification_changes=changed,independently_checked_windows=verified,arithmetic_max_fraction_of_derived_rounding_bound=native_normalized_bound_max,max_coefficient_discrepancy_vs_independent_QR=qr_discrepancy_max,QR_discrepancy_interpretation='Reported, not hidden by an acceptance tolerance. Existing inverse-table rounding remains in the centered reference.',common_profile_cases=6425,all_common_supports_and_statuses_verified=True,detailed_profile_rows=sum(actual_profiles.values()),detailed_curve_rows=sum(actual_curves.values()),production_changed=False)
    (refinement/'audit.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    with (refinement/'arithmetic-checks.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(selected[0]));w.writeheader();w.writerows(selected)
    return report

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('study',type=Path);p.add_argument('refinement',type=Path);p.add_argument('--matrix-header',type=Path,default=Path(__file__).resolve().parents[2]/'external/qalgorithms_matinverse.h');p.add_argument('--numpy-root',type=Path);a=p.parse_args()
    if a.numpy_root:sys.path.insert(0,str(a.numpy_root.resolve()))
    import numpy as np
    print(json.dumps(run(a.study,a.refinement,a.matrix_header,np),indent=2))
if __name__=='__main__':main()
