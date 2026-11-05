"""Exhaustive discrete hinge-pair coverage with finite local multistarts.
This assesses hinge sensitivity, not a certified global nonlinear optimum.
"""
import argparse,csv,json,math,sys,itertools,time
from pathlib import Path
from compare_components import Geometry,fit_start

def rows(path):
    with path.open(encoding='utf-8',newline='') as f:yield from csv.DictReader(f)

def run(comparison,out,np,scipy):
    from scipy.optimize import lsq_linear,least_squares
    out.mkdir(parents=True,exist_ok=True)
    cases=json.loads((comparison/'comparison-inputs.json').read_text(encoding='utf-8'))
    summaries={int(r['case']):r for r in rows(comparison/'summary.csv')}
    starts=list(rows(comparison/'starts.csv'));reports=[]
    sf=(out/'grid-starts.csv').open('w',newline='',encoding='utf-8');sw=csv.writer(sf)
    sw.writerow(['case','hinge1','hinge2','start','rss','success','status','nfev','error',*['theta_'+str(i) for i in range(8)]])
    pf=(out/'pairs.csv').open('w',newline='',encoding='utf-8');pw=csv.writer(pf)
    pw.writerow(['case','hinge1','hinge2','attempts','unsuccessful','best_rss','best_success','best_start'])
    nf=(out/'single-hinges.csv').open('w',newline='',encoding='utf-8');nw=csv.writer(nf);nw.writerow(['case','hinge','rss','success','status'])
    cf=(out/'curves.csv').open('w',newline='',encoding='utf-8');cw=csv.writer(cf);cw.writerow(['case','index','observed','single_log_median','two_log_median'])
    for c in [next(c for c in cases if c['scenario']==s and c['trial']==0) for s in (7,9)]:
        begin=time.perf_counter();n=len(c['y']);mid=(n-1)/2;hinges=sorted(set(map(float,range(1,n-1)))|{mid})
        geometries={h:Geometry(n,np,h) for h in hinges};g0=geometries[mid];z,anchor=g0.prepare(c['y'])
        nulls={h:g.one(z,lsq_linear) for h,g in geometries.items()}
        if not all(r['success'] for r in nulls.values()):raise RuntimeError('A convex hinge reference did not converge.')
        for h,r in nulls.items():nw.writerow([c['id'],h,r['rss'],r['success'],r['status']])
        null_h=min(nulls,key=lambda h:nulls[h]['rss']);null=nulls[null_h]
        reference=summaries[c['id']]
        refrow=next(r for r in starts if int(r['case'])==c['id'] and r['start']==reference['best_start'])
        reference_theta=np.array([float(refrow['theta_'+str(i)]) for i in range(8)])
        reference_curves=[g0.G@reference_theta[:4],g0.G@reference_theta[4:]]
        best=dict(rss=float(reference['two_rss']),theta=reference_theta.copy(),success=reference['status']=='best_found_converged',start='fixed_hinge_reference',h1=mid,h2=mid,status='reference')
        pairs=list(itertools.combinations_with_replacement(hinges,2));failures=0;attempts=0
        for h1,h2 in pairs:
            g1,g2=geometries[h1],geometries[h2];fits=[]
            if h1==h2:
                p=nulls[h1]['theta'].copy();p[0]-=math.log(2)*g1.scale[0];theta=np.r_[p,p]
                curve,_=g1.model_and_jacobian(theta,g2)
                fits.append(dict(start='nested_single',rss=float((curve-z)@(curve-z)),success=True,status='exact_embedding',nfev=0,error='',theta=theta))
            # Two warm label assignments from the preceding conditional fit.
            for swapped in range(1 if h1==h2 else 2):
                a=g1.one(reference_curves[swapped],lsq_linear);b=g2.one(reference_curves[1-swapped],lsq_linear)
                if not a['success'] or not b['success']:raise RuntimeError('A warm component projection failed.')
                fits.append(fit_start(g1,z,np.r_[a['theta'],b['theta']],f'warm_{swapped}',least_squares,other=g2))
            # One geometrically chosen split (and swapped assignment), strictly
            # a starting-value budget; not a boundary or distance for detection.
            split=min(n-3,max(3,int(math.floor((h1+h2)/2))+1))
            a=g0.gaussian_seed(z,0,split);b=g0.gaussian_seed(z,split,n)
            def translate(p,target):
                d=(target.hinge-mid)/(n-1);b0,b1,curv,_=p
                return np.array([b0+b1*d-curv*d*d,b1-2*curv*d,curv,curv])*target.scale
            for swapped in range(1 if h1==h2 else 2):
                p,q=(a,b) if not swapped else (b,a)
                fits.append(fit_start(g1,z,np.r_[translate(p,g1),translate(q,g2)],f'geometric_split_{split}_swap_{swapped}',least_squares,other=g2))
            finite=[r for r in fits if r['rss'] is not None];local=min(finite,key=lambda r:r['rss'])
            attempts+=len(fits);failed=sum(not r['success'] for r in fits);failures+=failed
            pw.writerow([c['id'],h1,h2,len(fits),failed,local['rss'],local['success'],local['start']])
            for r in fits:sw.writerow([c['id'],h1,h2,r['start'],r['rss'],r['success'],r['status'],r['nfev'],r['error'],*(['']*8 if r['theta'] is None else list(map(float,r['theta'])))])
            if local['rss']<best['rss']:best=dict(local,h1=h1,h2=h2)
        g1,g2=geometries[best['h1']],geometries[best['h2']]
        refined=fit_start(g1,z,best['theta'],'refine_grid_best',least_squares,other=g2,tight=True)
        sw.writerow([c['id'],best['h1'],best['h2'],refined['start'],refined['rss'],refined['success'],refined['status'],refined['nfev'],refined['error'],*(['']*8 if refined['theta'] is None else list(map(float,refined['theta'])))])
        if refined['rss'] is not None and refined['rss']<=best['rss']:best=dict(refined,h1=best['h1'],h2=best['h2'])
        curve,J=g1.model_and_jacobian(best['theta'],g2);singular=np.linalg.svd(J,compute_uv=False)
        for i in range(n):cw.writerow([c['id'],i,c['y'][i],null['curve'][i]+anchor,curve[i]+anchor])
        components=[g1.component(best['theta'][:4],anchor),g2.component(best['theta'][4:],anchor)]
        reports.append(dict(case=c['id'],scenario=c['scenario'],n=n,hinges=hinges,pairs_expected=len(hinges)*(len(hinges)+1)//2,pairs_processed=len(pairs),start_candidates=attempts,unsuccessful_starts=failures,fixed_one_rss=float(reference['one_rss']),fixed_two_rss=float(reference['two_rss']),profiled_one_rss=null['rss'],profiled_one_hinge=null_h,best_found_two_rss=best['rss'],best_hinges=[best['h1'],best['h2']],best_success=best['success'],best_start=best['start'],singular_values=list(map(float,singular)),components=components,fixed_deviance=float(reference['best_found_deviance']),profiled_best_found_deviance=n*math.log(null['rss']/best['rss']) if null['rss']>0 and best['rss']>0 else None,seconds=time.perf_counter()-begin))
        print(f"scenario={c['scenario']} pairs={len(pairs)} starts={attempts} best_success={best['success']} seconds={reports[-1]['seconds']:.3f}",flush=True)
    for f in (sf,pf,nf,cf):f.close()
    report=dict(cases=reports,all_discrete_pairs_covered=True,global_continuous_optimum_certified=False,statistical_model_count_decision=False,protocol='Hinges=all interior integer positions plus the fixed midpoint reference. Every unordered pair including duplicates. Warm projections in both label assignments and a geometric-midpoint Gaussian split. These are finite local starts, not a global certificate.')
    (out/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False),encoding='utf-8');return report

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('comparison',type=Path);p.add_argument('output',type=Path);p.add_argument('--numpy-root',type=Path);p.add_argument('--scipy-root',type=Path);a=p.parse_args()
    for path in (a.numpy_root,a.scipy_root):
        if path:sys.path.insert(0,str(path.resolve()))
    import numpy as np,scipy
    print(json.dumps(run(a.comparison,a.output,np,scipy),indent=2))
if __name__=='__main__':main()
