"""Independently recompute every real-window fit objective and bind original data."""
import argparse,csv,json,sys,math,hashlib
from collections import defaultdict,Counter
from pathlib import Path
from real_comparison import read_eics,select_windows,rows
from audit_components import score,geometry,single_active_set_rss,check_close

def run(export,debug,features,comparison,out,np):
    out.mkdir(parents=True,exist_ok=True)
    for name in ('eics.bin','inventory.csv','polarities.csv'):
        assert (export/name).read_bytes()==(debug/name).read_bytes(),name
    eics=read_eics(export/'eics.bin');protocol=json.loads((comparison/'protocol.json').read_text())
    assert hashlib.sha256((export/'eics.bin').read_bytes()).hexdigest()==protocol['eic_binary_sha256']
    annotation,selected,provenance=select_windows(eics,features,protocol['historical_positive_bins'])
    assert provenance==protocol['feature_provenance']
    ledger=rows(comparison/'feature-selection-ledger.csv')
    assert len(ledger)==len(annotation)==protocol['all_prior_feature_windows_verified']
    for a,b in zip(annotation,ledger):
        assert all(str(a[k])==b[k] for k in ('polarity','featureID','binID','start','length','old_collective'))
    cases=json.loads((comparison/'comparison-inputs.json').read_text())
    assert len(cases)==len(selected)==protocol['selected_windows']
    summary={int(r['case']):r for r in rows(comparison/'summary.csv')}
    assert len(summary)==len(cases)
    starts=defaultdict(list);curves=defaultdict(list);components=defaultdict(list)
    for r in rows(comparison/'starts.csv'):starts[int(r['case'])].append(r)
    for r in rows(comparison/'curves.csv'):curves[int(r['case'])].append(r)
    for r in rows(comparison/'components.csv'):components[int(r['case'])].append(r)
    checked=0;failed=0;max_error=0.;max_single_error=0.
    for case,(c,s) in enumerate(zip(cases,selected)):
        assert c['case']==case and all(c[k]==s[k] for k in s)
        d=eics[(c['polarity'],c['binID'])];begin=c['start'];n=c['length'];h=(n-1)/2
        for ck,dk in (('y','area'),('rt','rt'),('mz','mz'),('scan','scan')):assert tuple(c[ck])==d[dk][begin:begin+n]
        row=summary[case];y=np.asarray(c['y']);rs=starts[case]
        assert row['status']!='case_error',row
        assert len(rs)==int(row['starts'])==len({r['start'] for r in rs})
        assert {'nested_single_candidate','balanced_single_embedding',*('split_'+str(k) for k in range(3,n-2))}<=set(r['start'] for r in rs)
        for r in rs:
            if r['rss']:
                theta=np.array([float(r['theta_'+str(i)]) for i in range(8)])
                assert np.isfinite(theta).all() and np.all(theta[[2,3,6,7]]>=0)
                rss,*_=score(y,theta,h,h,np);reported=float(r['rss']);check_close(rss,reported,n,np)
                max_error=max(max_error,abs(rss-reported));checked+=1
            else:assert r['success']=='False' and r['error']
        failed+=sum(r['success']=='False' for r in rs)
        assert int(row['unsuccessful_starts'])==sum(r['success']=='False' for r in rs)
        best=next(r for r in rs if r['start']==row['best_start'])
        assert float(best['rss'])==min(float(r['rss']) for r in rs if r['rss'])==float(row['two_rss'])
        assert (best['success']=='True')==(row['status']=='best_found_converged')
        theta=np.array([float(best['theta_'+str(i)]) for i in range(8)])
        rss,mean,anchor,z,a,b=score(y,theta,h,h,np);G,scale=geometry(n,h,np)
        nested=next(r for r in rs if r['start']=='nested_single_candidate')
        p=np.array([float(nested['theta_'+str(i)]) for i in range(4)]);p[0]+=math.log(2)*scale[0]
        one=G@p;single=float(np.sum((one-z)**2));independent=single_active_set_rss(y,h,np)
        check_close(single,float(row['one_rss']),n,np);check_close(single,independent,n,np)
        max_single_error=max(max_single_error,abs(single-independent));assert rss<=single+256*np.finfo(float).eps*n*max(1.,single)
        cr=curves[case];assert len(cr)==n and [int(r['index']) for r in cr]==list(range(n))
        for i,r in enumerate(cr):
            assert float(r['observed'])==y[i] and float(r['rt'])==c['rt'][i]
            for key,value in (('one_log_median',one[i]+anchor),('two_log_median',mean[i]+anchor),('one_log_residual',one[i]-z[i]),('two_log_residual',mean[i]-z[i])):check_close(float(r[key]),float(value),n,np)
        assert {int(r['component']) for r in components[case]}=={0,1}
        for r in components[case]:
            j=int(r['component']);pp=theta[4*j:4*j+4]/scale
            expected=[pp[0]+anchor,pp[1]/(n-1),-pp[2]/(n-1)**2,-pp[3]/(n-1)**2]
            for key,value in zip(('log_intercept','slope_per_sample','left_curvature','right_curvature'),expected):check_close(float(r[key]),float(value),n,np)
    report=json.loads((comparison/'report.json').read_text());assert failed==report['unsuccessful_starts']
    assert dict(Counter(r['status'] for r in summary.values()))==report['status_counts']
    result=dict(debug_release_all_export_files_identical=True,fresh_eics=len(eics),all_prior_windows_bound=len(annotation),real_windows_audited=len(cases),start_scores_independently_recomputed=checked,unsuccessful_starts_retained=failed,max_score_error=max_error,independent_convex_single_max_error=max_single_error,physical_peak_count_validated=False)
    (out/'audit.json').write_text(json.dumps(result,indent=2),encoding='utf-8');return result

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('export','debug','features','comparison','output'):p.add_argument(name,type=Path)
    p.add_argument('--numpy-root',type=Path);a=p.parse_args()
    if a.numpy_root:sys.path.insert(0,str(a.numpy_root.resolve()))
    import numpy as np
    print(json.dumps(run(a.export,a.debug,a.features,a.comparison,a.output,np),indent=2))
if __name__=='__main__':main()
