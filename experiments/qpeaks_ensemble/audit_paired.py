"""Independent field, output, full hinge-score and best-curve verification."""
import argparse,csv,json,math,struct,hashlib,sys,re
from pathlib import Path
from collections import defaultdict
from functools import lru_cache
from analyse import read_inputs
from compare_bundle_cascade import fields,rows
from bundle_method import topology,bundle_support,geometric_peak
from real_comparison import read_eics

def run(inputs,release,debug,comparison,old_collection,fresh,out,np):
    out.mkdir(parents=True,exist_ok=True)
    comparison_report=json.loads((comparison/'report.json').read_text())
    model_boundaries=comparison_report.get('topology_version','center_slope_v1')!='center_slope_v1'
    boundary_scope='singleton_blocks' if comparison_report.get('topology_version')=='model_singleton_v3' else 'all_edges'
    assert (release/'fields.bin').read_bytes()==(debug/'fields.bin').read_bytes()
    for name in ('cascade-runs.csv','cascade-peaks.csv','field-timing.csv'):
        a=list(rows(release/name));b=list(rows(debug/name));assert len(a)==len(b)
        ignored={'seconds','kernel_seconds'}
        assert [{k:v for k,v in r.items() if k not in ignored} for r in a]==[{k:v for k,v in r.items() if k not in ignored} for r in b],name
    source=read_inputs(inputs/'inputs.bin');cases={(c['scenario'],c['trial']):c for c in json.loads((inputs/'cases.json').read_text())}
    policy=json.loads((release/'field-policy.json').read_text());assert policy==json.loads((debug/'field-policy.json').read_text())
    threshold=policy['minimum_intensity'];oldrows=defaultdict(list);bad_prefix={}
    for key,(truth,x,y) in source.items():
        cumulative=[0]
        for v in y:cumulative.append(cumulative[-1]+int(v<=threshold))
        bad_prefix[key]=cumulative
    for r in rows(old_collection/'candidates.csv'):
        key=int(r['scenario']),int(r['trial'])
        if key not in source:continue
        oldrows[key] # Keep an explicit empty reference if the domain removes all fits.
        center,scale=int(r['center']),int(r['scale']);bad=bad_prefix[key]
        if bad[center+scale+1]!=bad[center-scale]:continue
        oldrows[key].append((center,scale,*[float(r[k]) for k in ('b0','b1','b2','b3')]))
    oldhashes={};oldcounts={}
    for key,rs in oldrows.items():
        digest=hashlib.sha256()
        for r in sorted(rs):digest.update(struct.pack('<HH4d',*r))
        oldhashes[key]=digest;oldcounts[key]=len(rs)
    del oldrows
    scored=defaultdict(list)
    with (comparison/'hinge-scores.bin').open('rb') as f:
        assert f.read(8)==b'QSCR0001'
        while True:
            b=f.read(32)
            if not b:break
            assert len(b)==32;s,t,g,h,rss,fast=struct.unpack('<4I2d',b);scored[(s,t,g)].append((h,rss,fast))
    fits=defaultdict(dict)
    for r in json.loads((comparison/'bundle-fits.json').read_text()):
        key=r['scenario'],r['trial'];assert r['bundle'] not in fits[key];fits[key][r['bundle']]=r
    summaries={(c['scenario'],c['trial']):c for c in json.loads((comparison/'case-results.json').read_text())}
    @lru_cache(maxsize=4096)
    def reference(n,h):
        t=(np.arange(n)-h)/(n-1);X=np.column_stack([np.ones(n),t,np.where(t<0,t*t,0),np.where(t>0,t*t,0)])
        norm=np.linalg.norm(X,axis=0);G=X/norm;U,s,_=np.linalg.svd(G,full_matrices=False)
        return U,float(s[0]/s[-1])
    checked=0;max_score_error=0.;max_curve_error=0.;old_verified=0;native_count=0;real_verified=[];fresh_eics=read_eics(fresh)
    for identity,z,candidates,intervals in fields(release/'fields.bin'):
        meta=cases[identity];truth,x,y=source[identity];bs,edges,groups,membership=topology(candidates,model_boundaries=model_boundaries,log_values=z,boundary_scope=boundary_scope);by_id={r['id']:r for r in candidates}
        native_count+=len(candidates);assert len(candidates)==summaries[identity]['native_fits']
        if identity in oldhashes:
            digest=hashlib.sha256()
            for r in sorted(candidates,key=lambda r:(r['center'],r['scale'])):digest.update(struct.pack('<HH4d',r['center'],r['scale'],r['b0'],r['b1'],r['b2'],r['b3']))
            assert digest.digest()==oldhashes[identity].digest() and len(candidates)==oldcounts[identity];old_verified+=1
        assert set(fits[identity])=={g['id'] for g in groups}
        for g in groups:
            record=fits[identity][g['id']];assert record['status']!='explicit_fit_error'
            left,right,hinges=bundle_support(g,bs,candidates,by_id);assert (left,right)==(record['left'],record['right'])
            n=right-left+1;local=np.asarray(z[left:right+1]);anchor=local[n//2] if n%2 else (local[n//2-1]+local[n//2])/2;w=local-anchor
            scores=scored[(*identity,g['id'])];assert len(scores)==len(hinges) and {h-left for h,rss,fast in scores}==set(hinges)
            assert min(rss for h,rss,fast in scores)==record['best_rss']
            for h,rss,fast in scores:
                Q,condition=reference(n,h-left);residual=w-Q@(Q.T@w);independent=float(residual@residual)
                error=abs(independent-rss);tolerance=2048*np.finfo(float).eps*n*max(1.,rss,independent)
                assert error<=tolerance,(identity,g['id'],h,error,tolerance)
                max_score_error=max(max_score_error,error);checked+=1
            h=record['hinge_global']-left;p=record['coefficients'];u=np.arange(n)-h
            curve=p[0]+p[1]*u+np.where(u<0,p[2],p[3])*u*u
            Q,condition=reference(n,h);independent=Q@(Q.T@w)+anchor
            curve_error=float(np.max(np.abs(curve-independent)));max_curve_error=max(max_curve_error,curve_error)
            assert curve_error<=4096*np.finfo(float).eps*n*max(1.,float(np.linalg.norm(local)))*condition**2
            if record['numeric_minimizers']==1:
                peak=geometric_peak(p,h,n);eligible=peak['apex'] is not None and peak['above_both_support_edges'] is True
                assert eligible==record['eligible_geometric_proposal']
                assert record['apex']==(None if peak['apex'] is None else left+peak['apex'])
                assert not record['statistically_accepted']
        closed=[r for r in fits[identity].values() if r['eligible_geometric_proposal']]
        assert len(closed)==summaries[identity]['geometry_only_proposals']
        match=re.fullmatch(r'real_EIC(\d+)_observations',meta['name'])
        key=(meta['polarity'],meta['binID']) if 'polarity' in meta else (1,int(match.group(1))) if match else None
        if key is not None:
            data=fresh_eics[key];assert tuple(x)==data['rt'] and tuple(y)==data['area']
            real_verified.append(dict(scenario=identity[0],polarity=key[0],binID=key[1],n=len(y)))
        print(f'audited={identity} groups={len(groups)} scores_total={checked}',flush=True)
    report=json.loads((comparison/'report.json').read_text());assert checked==report['hinge_scores'] and native_count==report['native_regressions_accounted']
    result=dict(native_field_policy=policy,debug_release_native_fields_identical=True,debug_release_all_cascade_outputs_identical_except_timing=True,cases=len(source),all_native_regressions=native_count,prior_native_detailed_fields_exactly_bound=old_verified,all_hinge_scores_independently_recomputed=checked,all_best_curves_independently_recomputed=sum(len(v) for v in fits.values()),max_hinge_rss_error=max_score_error,max_best_curve_log_error=max_curve_error,real_vectors_bound_to_fresh_mzml=real_verified,statistical_ensemble_peak_decision_validated=False)
    result['topology_version']=comparison_report.get('topology_version','center_slope_v1')
    result['verified_hashes']={name:hashlib.sha256(path.read_bytes()).hexdigest() for name,path in [('inputs',inputs/'inputs.bin'),('fields',release/'fields.bin'),('hinge_scores',comparison/'hinge-scores.bin'),('bundle_fits',comparison/'bundle-fits.json'),('case_results',comparison/'case-results.json')]}
    (out/'audit.json').write_text(json.dumps(result,indent=2),encoding='utf-8');return result

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('inputs','release','debug','comparison','old_collection','fresh','output'):p.add_argument(k,type=Path)
    p.add_argument('--numpy-root',type=Path);a=p.parse_args()
    if a.numpy_root:sys.path.insert(0,str(a.numpy_root.resolve()))
    import numpy as np
    print(json.dumps(run(a.inputs,a.release,a.debug,a.comparison,a.old_collection,a.fresh,a.output,np),indent=2))
