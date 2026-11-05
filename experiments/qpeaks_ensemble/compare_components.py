"""One vs two additive concave log-profile components on fixed supports.
Conditional reference fits, not a calibrated component-count test or detector.
"""
import argparse,csv,json,math,sys,time,hashlib
from pathlib import Path
from collections import Counter,defaultdict
from analyse import read_inputs

GENERATED_CONTROLS = {
    '7000': {'kind':'paired_reflection_of_overlap_trial_zero_on_same_support'},
    '7001': {'n':29,'apices':[10,18],'heights':[1000,500],'widths':[2.5,2.5],'log_noise_sd':0},
    '7002': {'n':41,'apices':[19,21],'heights':[1000,500],'widths':[3,3],'log_noise_sd':.1,'seed':20261105+9001},
    '7003': {'n':41,'apices':[20],'heights':[1000],'widths':[3],'background':240,'log_noise_sd':.1,'seed':20261105+9002},
}

class Geometry:
    def __init__(self,n,np,hinge=None):
        self.n=n;self.np=np;self.span=n-1
        self.hinge=(n-1)/2 if hinge is None else float(hinge)
        if n<=8:raise ValueError('Two four-coefficient components require n>8 for this unsaturated reference.')
        u=(np.arange(n)-self.hinge)/self.span
        self.u=u
        self.X=np.column_stack([np.ones(n),u,-np.where(u<0,u*u,0),-np.where(u>0,u*u,0)])
        self.scale=np.linalg.norm(self.X,axis=0);self.G=self.X/self.scale
        self.lower=np.array([-np.inf,-np.inf,0.,0.]);self.upper=np.full(4,np.inf)
        self.tol=np.finfo(float).eps**.5;self.tight=np.finfo(float).eps**(2/3)

    def prepare(self,y):
        np=self.np;y=np.asarray(y,dtype=float)
        if len(y)!=self.n or not np.isfinite(y).all() or np.any(y<=0):raise ValueError('Finite positive observations on the complete support are required.')
        z=np.log(y);anchor=float(z[self.n//2]) if self.n%2 else float((z[self.n//2-1]+z[self.n//2])/2)
        return z-anchor,anchor

    def one(self,z,solver):
        np=self.np
        if np.all(z==z[0]):
            p=np.zeros(4);p[0]=z[0]*self.scale[0];return dict(theta=p,rss=0.,curve=z.copy(),success=True,status='analytic_constant',optimality=0.,active=[0,0,-1,-1])
        r=solver(self.G,z,bounds=(self.lower,self.upper),method='bvls',tol=self.tight)
        residual=self.G@r.x-z
        return dict(theta=r.x,rss=float(residual@residual),curve=self.G@r.x,success=bool(r.success),status=str(r.status),optimality=float(r.optimality),active=list(map(int,r.active_mask)))

    def model_and_jacobian(self,theta,other=None):
        np=self.np;g2=self.G if other is None else other.G
        with np.errstate(over='raise',invalid='raise',divide='raise'):
            q1=self.G@theta[:4];q2=g2@theta[4:]
            mean=np.logaddexp(q1,q2);w1=np.exp(q1-mean);w2=np.exp(q2-mean)
            jac=np.column_stack([w1[:,None]*self.G,w2[:,None]*g2])
        if not np.isfinite(mean).all() or not np.isfinite(jac).all():raise FloatingPointError('Nonfinite model or Jacobian.')
        return mean,jac

    def gaussian_seed(self,z,first,last):
        np=self.np;u=self.u[first:last];X=np.column_stack([np.ones(len(u)),u,-u*u])
        coef=np.linalg.lstsq(X,z[first:last],rcond=None)[0]
        if coef[2]<0:
            coef[:2]=np.linalg.lstsq(X[:,:2],z[first:last],rcond=None)[0];coef[2]=0
        return np.array([coef[0],coef[1],coef[2],coef[2]])

    def starts(self,z,one,other=None):
        np=self.np;other=self if other is None else other
        result=[]
        if other.hinge==self.hinge:
            p=one['theta'].copy();p[0]-=math.log(2)*self.scale[0]
            result.append(('balanced_single_embedding',np.r_[p,p]))
        for split in range(3,self.n-2):
            a=self.gaussian_seed(z,0,split);b=self.gaussian_seed(z,split,self.n)
            # The Gaussian seed is translated exactly to either hinge; the
            # final two-component fit is not restricted to Gaussian curvature.
            def translated(p,target):
                d=(target.hinge-self.hinge)/self.span;b0,b1,c,_=p
                return np.array([b0+b1*d-c*d*d,b1-2*c*d,c,c])*target.scale
            result.append((f'split_{split}',np.r_[translated(a,self),translated(b,other)]))
            if other.hinge!=self.hinge:result.append((f'split_{split}_swapped',np.r_[translated(b,self),translated(a,other)]))
        return result

    def component(self,theta,anchor):
        b0,b1,cL,cR=map(float,theta/self.scale)
        c=cL if b1<0 else cR
        apex=self.hinge+self.span*b1/(2*c) if c>0 else None
        if apex is not None and not math.isfinite(apex):apex=None
        return dict(log_intercept=b0+anchor,slope_per_sample=b1/self.span,left_curvature=-cL/self.span**2,right_curvature=-cR/self.span**2,apex=apex,interior_apex=apex is not None and 0<apex<self.n-1,zero_curvature_boundary=(cL==0 or cR==0))


def fit_start(geometry,z,theta,name,least_squares,other=None,tight=False):
    np=geometry.np;other=geometry if other is None else other
    cache={}
    def evaluate(p):
        if 'p' not in cache or not np.array_equal(cache['p'],p):
            mean,jac=geometry.model_and_jacobian(p,other);cache.update(p=p.copy(),r=mean-z,j=jac)
        return cache['r'],cache['j']
    try:
        tol=geometry.tight if tight else geometry.tol
        r=least_squares(lambda p:evaluate(p)[0].copy(),theta,jac=lambda p:evaluate(p)[1].copy(),bounds=(np.r_[geometry.lower,other.lower],np.r_[geometry.upper,other.upper]),method='trf',loss='linear',x_scale='jac',ftol=tol,xtol=tol,gtol=tol,max_nfev=(400 if tight else 100)*8)
        residual,jac=evaluate(r.x);score=float(residual@residual)
        if not math.isfinite(score):raise FloatingPointError('Nonfinite residual sum.')
        return dict(start=name,rss=score,success=bool(r.success),status=str(r.status),nfev=int(r.nfev),optimality=float(r.optimality),theta=r.x.copy(),active=list(map(int,r.active_mask)),error='')
    except Exception as e:
        return dict(start=name,rss=None,success=False,status='exception',nfev=0,optimality=None,theta=None,active=[],error=str(e))


def compare(y,geometry,linear_solver,nonlinear_solver,refine=True):
    np=geometry.np;z,anchor=geometry.prepare(y);one=geometry.one(z,linear_solver)
    if not one['success']:raise RuntimeError('Convex single-component reference did not converge.')
    starts=geometry.starts(z,one);fits=[]
    # Include the exact nested single-component representation as a candidate.
    name,p=starts[0];curve,jac=geometry.model_and_jacobian(p)
    fits.append(dict(start='nested_single_candidate',rss=float((curve-z)@(curve-z)),success=True,status='exact_embedding',nfev=0,optimality=None,theta=p.copy(),active=[],error=''))
    # Constant data have a proven RSS lower bound of zero. No optimizer is
    # needed to invent a two-peak interpretation of this boundary case.
    if not np.all(z==z[0]):
        for name,p in starts:fits.append(fit_start(geometry,z,p,name,nonlinear_solver))
    finite=[f for f in fits if f['rss'] is not None]
    best=min(finite,key=lambda f:f['rss'])
    if refine and best['status']!='exact_embedding':
        refinement=fit_start(geometry,z,best['theta'],'refine_best',nonlinear_solver,tight=True);fits.append(refinement)
        if refinement['rss'] is not None and refinement['rss']<=best['rss']:best=refinement
    curve,jac=geometry.model_and_jacobian(best['theta']);singular=np.linalg.svd(jac,compute_uv=False)
    rank_tolerance=max(jac.shape)*np.finfo(float).eps*singular[0]
    rank=int(sum(singular>rank_tolerance))
    condition=float(singular[0]/singular[-1]) if singular[-1]>0 else None
    if condition is not None and not math.isfinite(condition):condition=None
    components=[geometry.component(best['theta'][:4],anchor),geometry.component(best['theta'][4:],anchor)]
    # Preserve optimizer block IDs across parameter and curve exports. Label
    # permutation is a model symmetry, not an implicit reordering of one file.
    deviance=geometry.n*math.log(one['rss']/best['rss']) if one['rss']>0 and best['rss']>0 else None
    return dict(one=one,best=best,starts=fits,anchor=anchor,z=z,curve_two=curve,jacobian=jac,singular=singular,rank=rank,condition=condition,components=components,deviance=deviance)


def make_protocol(common,np,trials):
    records=read_inputs(common/'common-inputs.bin');supports=json.loads((common/'supports.json').read_text(encoding='utf-8'));cases=[];ledger=[]
    for info in supports:
        case=info['case'];truth,x,y=records[(case,0)];scenario=info['scenario'];trial=info['trial'];variant=info['variant']
        pick=(scenario in (4,4010) and variant=='whole_input' and trial<trials) or (scenario==7 and variant=='fixed_24_52' and trial<trials) or (trial==0 and scenario in (8,9,3000,3001,4001,4002,4003,1004,4011) and variant=='whole_input')
        eligible=len(y)>8 and all(v>0 for v in y)
        reason='selected' if pick and eligible else ('nonpositive_log_input' if any(v<=0 for v in y) else 'n_not_greater_than_8' if len(y)<=8 else 'outside_predeclared_first_trials_protocol')
        ledger.append(dict(source_case=case,scenario=scenario,trial=trial,variant=variant,selected=bool(pick and eligible),reason=reason))
        if pick and eligible:
            sigma=.05 if scenario==7 else .1 if scenario in (4,1004,4010,4011) else None
            cases.append(dict(id=len(cases),name=info['name'],scenario=scenario,trial=trial,variant=variant,source_case=case,origin=info['start'],truth=[p-info['start'] for p in info['truth']],sigma=sigma,y=list(y)))
    source=next(c for c in cases if c['scenario']==7 and c['trial']==0)
    cases.append(dict(id=len(cases),name='overlap_reflected_same_support',scenario=7000,trial=0,variant='generated_paired_reflection',source_case=source['source_case'],origin=0,truth=sorted(len(source['y'])-1-p for p in source['truth']),sigma=.05,y=source['y'][::-1]))
    def gaussian(x,a,H,s):return H*np.exp(-(x-a)**2/(2*s*s))
    x=np.arange(29,dtype=float)
    y=(gaussian(x,10,1000,2.5)+gaussian(x,18,500,2.5)).astype(np.float32)
    cases.append(dict(id=len(cases),name='two_gaussians_noiseless',scenario=7001,trial=0,variant='generated',source_case=None,origin=0,truth=[10,18],sigma=None,y=list(map(float,y))))
    x=np.arange(41,dtype=float);mean=gaussian(x,19,1000,3)+gaussian(x,21,500,3)
    rng=np.random.default_rng(20261105+9001)
    y=(mean*np.exp(rng.normal(0,.1,41))).astype(np.float32)
    cases.append(dict(id=len(cases),name='close_pair_log_noise',scenario=7002,trial=0,variant='generated',source_case=None,origin=0,truth=[19,21],sigma=.1,y=list(map(float,y))))
    rng=np.random.default_rng(20261105+9002)
    mean=gaussian(x,20,1000,3)+240
    y=(mean*np.exp(rng.normal(0,.1,41))).astype(np.float32)
    cases.append(dict(id=len(cases),name='one_peak_plus_constant_background',scenario=7003,trial=0,variant='generated_misspecification_control',source_case=None,origin=0,truth=[20],sigma=.1,y=list(map(float,y))))
    return cases,ledger


def run(common,out,np,scipy,plots=True,trials=32,pilot=False):
    from scipy.optimize import lsq_linear,least_squares
    out.mkdir(parents=True,exist_ok=True);cases,ledger=make_protocol(common,np,trials)
    if pilot:cases=[next(c for c in cases if c['scenario']==s and c['trial']==0) for s in (4,4010,7)]
    with (out/'selection-ledger.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(ledger[0]));w.writeheader();w.writerows(ledger)
    (out/'comparison-inputs.json').write_text(json.dumps(cases,indent=2,allow_nan=False),encoding='utf-8')
    protocol=dict(trials_per_main_generator=trials,pilot=pilot,component_family='q=b0+b1*u-cL*u^2(left)-cR*u^2(right), cL,cR>=0; total intensity is exp(q1)+exp(q2)',observation_model='log(y)=logaddexp(q1,q2)+iid N(0,sigma^2); noise applies to the total intensity',hinge='fixed at support midpoint in both models; conditional reference, no full hinge search',starts='balanced exact single embedding plus Gaussian seeds from EVERY split with >=3 points per half; all final objectives use all original observations once',numerical_tolerance=float(np.finfo(float).eps**.5),refinement_tolerance=float(np.finfo(float).eps**(2/3)),nfev_budget_per_start=800,nfev_refinement_budget=3200,minimum_peak_distance=None,minimum_component_amplitude=None,peak_count_decision=False,known_shape_model_assumption=True,scipy=scipy.__version__,numpy=np.__version__,generated_close_pair_seed=20261105+9001,generated_background_seed=20261105+9002,constant_strategy='analytic zero-RSS lower-bound control; no optimizer multistart needed',generated_control_definitions=GENERATED_CONTROLS,source_common_input_sha256=hashlib.sha256((common/'common-inputs.bin').read_bytes()).hexdigest())
    (out/'protocol.json').write_text(json.dumps(protocol,indent=2),encoding='utf-8')
    summary=[];cache={};allfits={};start_file=(out/'starts.csv').open('w',newline='',encoding='utf-8');sw=csv.writer(start_file)
    sw.writerow(['case','start','rss','success','status','nfev','optimality','error',*['theta_'+str(i) for i in range(8)]])
    curve_file=(out/'curves.csv').open('w',newline='',encoding='utf-8');cw=csv.writer(curve_file);cw.writerow(['case','index','observed','one_log_median','two_log_median','component1_log_median','component2_log_median'])
    param_file=(out/'components.csv').open('w',newline='',encoding='utf-8');pw=csv.writer(param_file);pw.writerow(['case','component','log_intercept','slope_per_sample','left_curvature','right_curvature','apex','interior_apex','zero_curvature_boundary'])
    for c in cases:
        started=time.perf_counter();n=len(c['y'])
        if n not in cache:cache[n]=Geometry(n,np)
        g=cache[n]
        try:
            result=compare(c['y'],g,lsq_linear,least_squares);one=result['one'];best=result['best'];allfits[c['id']]=result
            failed=sum(not f['success'] for f in result['starts']);finite=[f['rss'] for f in result['starts'] if f['rss'] is not None]
            row=dict(case=c['id'],name=c['name'],scenario=c['scenario'],trial=c['trial'],n=n,source_case=c['source_case'],status='best_found_converged' if best['success'] else 'best_found_unconverged',one_rss=one['rss'],two_rss=best['rss'],gain=one['rss']-best['rss'],best_found_deviance=result['deviance'],known_sigma_gain=(one['rss']-best['rss'])/c['sigma']**2 if c['sigma'] else None,starts=len(result['starts']),unsuccessful_starts=failed,best_start=best['start'],min_finite_start_rss=min(finite),max_finite_start_rss=max(finite),jacobian_rank=result['rank'],jacobian_condition=result['condition'],singular_values=json.dumps(list(map(float,result['singular']))),seconds=time.perf_counter()-started,error='')
            for f in result['starts']:
                sw.writerow([c['id'],f['start'],f['rss'],f['success'],f['status'],f['nfev'],f['optimality'],f['error'],*(['']*8 if f['theta'] is None else list(map(float,f['theta'])))])
            for j,p in enumerate(result['components']):pw.writerow([c['id'],j,*p.values()])
            if c['trial']==0:
                q1=g.G@best['theta'][:4]+result['anchor'];q2=g.G@best['theta'][4:]+result['anchor']
                for i in range(n):cw.writerow([c['id'],i,c['y'][i],one['curve'][i]+result['anchor'],result['curve_two'][i]+result['anchor'],q1[i],q2[i]])
        except Exception as e:
            row=dict(case=c['id'],name=c['name'],scenario=c['scenario'],trial=c['trial'],n=n,source_case=c['source_case'],status='case_error',one_rss=None,two_rss=None,gain=None,best_found_deviance=None,known_sigma_gain=None,starts=0,unsuccessful_starts=0,best_start=None,min_finite_start_rss=None,max_finite_start_rss=None,jacobian_rank=None,jacobian_condition=None,singular_values=None,seconds=time.perf_counter()-started,error=str(e))
        summary.append(row)
        print(f"case={c['id']} scenario={c['scenario']} trial={c['trial']} status={row['status']} starts={row['starts']} seconds={row['seconds']:.3f}",flush=True)
    start_file.close();curve_file.close();param_file.close()
    with (out/'summary.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(summary[0]));w.writeheader();w.writerows(summary)
    aggregate=[]
    for scenario in sorted({r['scenario'] for r in summary}):
        rs=[r for r in summary if r['scenario']==scenario];ok=[r for r in rs if r['two_rss'] is not None]
        aggregate.append(dict(scenario=scenario,cases=len(rs),errors=len(rs)-len(ok),best_unconverged=sum(r['status']=='best_found_unconverged' for r in rs),unsuccessful_starts=sum(r['unsuccessful_starts'] for r in rs),median_one_rss=float(np.median([r['one_rss'] for r in ok])) if ok else None,median_two_rss=float(np.median([r['two_rss'] for r in ok])) if ok else None,median_best_found_deviance=float(np.median([r['best_found_deviance'] for r in ok if r['best_found_deviance'] is not None])) if any(r['best_found_deviance'] is not None for r in ok) else None))
    report=dict(cases=len(cases),status_counts=dict(Counter(r['status'] for r in summary)),aggregates=aggregate,automatic_peak_count_decisions=False,global_two_component_optimum_certified=False,nominal_chi_square_pvalues_reported=False,production_changed=False)
    (out/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False),encoding='utf-8')
    if plots:
        import matplotlib
        matplotlib.use('Agg');import matplotlib.pyplot as plt
        fig,axes=plt.subplots(2,2,figsize=(12,7),layout='constrained')
        for ax,scenario,title in zip(axes.flat,(4,4010,7,9),('Ein Gaußpeak','Ein asymmetrischer Peak','Zwei überlappende Peaks','pt_04: mehr als zwei Strukturen')):
            selected=next((c for c in cases if c['scenario']==scenario and c['trial']==0),None)
            if selected is None:ax.set_visible(False);continue
            c=selected;r=allfits.get(c['id']);ax.plot(range(len(c['y'])),c['y'],label='Messwerte',color='#23546e')
            if r is not None:
                ax.plot(np.exp(r['one']['curve']+r['anchor']),label='Ein Profil',color='#b14c30');ax.plot(np.exp(r['curve_two']+r['anchor']),label='Summe zweier Profile',color='#408553')
            ax.set_title(title);ax.set_xlabel('Index im festgelegten Bereich');ax.set_ylabel('Intensität');ax.legend(fontsize=8)
        fig.suptitle('Ein- gegen Zweikomponentenmodell: gemeinsamer Datenbereich, gleiche Komponentenfamilie\nKleinere Residuen sind keine kalibrierte Entscheidung für zwei Peaks.',fontsize=11)
        fig.savefig(out/'model-comparison.png',dpi=150);fig.savefig(out/'model-comparison.pdf');plt.close(fig)
        fig,axes=plt.subplots(1,3,figsize=(14,4),layout='constrained')
        for ax,scenario,title in zip(axes,(3000,7002,7003),('Konstant: Komponenten nicht eindeutig','Zwei nahe Komponenten, ein Signalmaximum','Ein Peak mit zusätzlichem Hintergrund')):
            selected=next((c for c in cases if c['scenario']==scenario and c['trial']==0),None)
            if selected is None:ax.set_visible(False);continue
            c=selected;r=allfits.get(c['id']);ax.plot(c['y'],label='Messwerte',color='#23546e')
            if r is not None:
                g=cache[len(c['y'])];p=r['best']['theta']
                ax.plot(np.exp(r['one']['curve']+r['anchor']),label='Ein Profil',color='#b14c30')
                ax.plot(np.exp(r['curve_two']+r['anchor']),label='Summe zweier Profile',color='#408553')
                ax.plot(np.exp(g.G@p[:4]+r['anchor']),ls=':',color='#888888',label='Einzelkomponenten')
                ax.plot(np.exp(g.G@p[4:]+r['anchor']),ls=':',color='#888888')
            ax.set_title(title,fontsize=10);ax.set_xlabel('Index im festgelegten Bereich');ax.set_ylabel('Intensität');ax.legend(fontsize=7)
        fig.suptitle('Kontrollen: Zwei Modellkomponenten bedeuten nicht automatisch zwei physische Peaks',fontsize=11)
        fig.savefig(out/'identifiability-controls.png',dpi=150);fig.savefig(out/'identifiability-controls.pdf');plt.close(fig)

    return report


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('common',type=Path);p.add_argument('output',type=Path);p.add_argument('--numpy-root',type=Path);p.add_argument('--scipy-root',type=Path);p.add_argument('--trials',type=int,default=32);p.add_argument('--pilot',action='store_true');p.add_argument('--no-plots',action='store_true');a=p.parse_args()
    for path in (a.numpy_root,a.scipy_root):
        if path:sys.path.insert(0,str(path.resolve()))
    import numpy as np,scipy
    if a.trials<1 or a.trials>300:raise ValueError('Protocol trial budget must be in 1..300.')
    report=run(a.common,a.output,np,scipy,not a.no_plots,a.trials,a.pilot);print(json.dumps(report,indent=2))
if __name__=='__main__':main()
