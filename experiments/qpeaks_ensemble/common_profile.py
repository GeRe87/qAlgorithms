"""Common-support piecewise-quadratic log-profile reference, not a peak detector.
All geometry-dependent QR operators are cached before fitting each support size.
No fit votes, automatic support choice, p-values, areas or production changes.
"""
import argparse,csv,json,math,struct,sys
from pathlib import Path
from collections import Counter,defaultdict
from analyse import read_inputs

class ProfileGeometry:
    def __init__(self,n,np):
        if n<5:raise ValueError('Four coefficients and positive residual df require n>=5.')
        self.n=n;self.np=np;self.span=n-1;self.hinges=np.arange(1,n-1)
        operators=[];projectors=[]
        for h in self.hinges:
            t=np.arange(n,dtype=float)-h
            X=np.column_stack([np.ones(n),t,np.where(t<0,t*t,0),np.where(t>0,t*t,0)])
            norm=np.linalg.norm(X,axis=0)
            scaled=X/norm
            Q,R=np.linalg.qr(scaled,mode='reduced')
            normalized_operator=np.linalg.solve(R,Q.T)
            # Check the dimensionless basis. Undoing column units can greatly
            # amplify off-diagonal identity errors for a one-point tail.
            rounding_bound=128*n*np.finfo(float).eps*np.linalg.norm(normalized_operator,2)*np.linalg.norm(scaled,2)
            assert np.linalg.norm(normalized_operator@scaled-np.eye(4))<=rounding_bound
            A=normalized_operator/norm[:,None]
            operators.append(A);projectors.append(Q)
        self.A=np.asarray(operators);self.Q=np.asarray(projectors)

    def fit(self,y):
        np=self.np;y=np.asarray(y,dtype=float)
        if len(y)!=self.n or not np.isfinite(y).all() or np.any(y<=0):raise ValueError('Common log support requires finite positive observations.')
        z=np.log(y)
        anchor=float(z[self.n//2]) if self.n%2 else float((z[self.n//2-1]+z[self.n//2])/2)
        centered=z-anchor
        beta=np.einsum('hkn,n->hk',self.A,centered)
        projection=np.einsum('hnk,hk->hn',self.Q,np.einsum('hnk,n->hk',self.Q,centered))
        residual=centered[None,:]-projection;rss=np.einsum('hn,hn->h',residual,residual)
        beta[:,0]+=anchor
        assert np.isfinite(beta).all() and np.isfinite(rss).all()
        return z,beta,rss,projection+anchor

    def apex(self,beta,hinge):
        b0,b1,b2,b3=map(float,beta);left=-hinge;right=self.n-1-hinge
        root=None
        if b1<0 and b2<0:
            value=-b1/(2*b2)
            if left<value<0:root=value
        elif b1>0 and b3<0:
            value=-b1/(2*b3)
            if 0<value<right:root=value
        elif b1==0 and b2<0 and b3<0:root=0.0
        return None if root is None else hinge+root


def summarize_fit(y,geometry):
    np=geometry.np;z,beta,rss,curves=geometry.fit(y)
    minima=np.flatnonzero(rss==rss.min())
    result=dict(status='interior_maximum',n=len(y),nonpositive=0,minimizers=len(minima),hinge=None,peak_position=None,rss=float(rss.min()),b0=None,b1=None,b2=None,b3=None)
    index=None
    if np.all(z==z[0]):result['status']='constant_data'
    elif len(minima)!=1:result['status']='multiple_equal_minimizers'
    else:
        index=int(minima[0]);h=int(geometry.hinges[index]);position=geometry.apex(beta[index],h)
        result.update(hinge=h,peak_position=position,**dict(zip(('b0','b1','b2','b3'),map(float,beta[index]))))
        if position is None:result['status']='no_interior_maximum'
    return result,(z,beta,rss,curves,index)


def cases_from_study(folder,np):
    records=read_inputs(folder/'inputs.bin');meta=json.loads((folder/'scenarios.json').read_text(encoding='utf-8'))
    names={s['id']:s['name'] for s in meta['scenarios']};cases=[]
    def add(scenario,trial,truth,x,y,start,end,variant,name,known_hinge=None):
        cases.append(dict(case=len(cases),scenario=scenario,trial=trial,name=name,variant=variant,start=start,end=end,truth=list(truth),x=list(x[start:end+1]),y=list(y[start:end+1]),known_hinge_global=known_hinge))
    for (scenario,trial),(truth,x,y) in records.items():
        known=(len(y)-1)/2 if scenario in (4,1004,3001) else None
        add(scenario,trial,truth,x,y,0,len(y)-1,'whole_input',names[scenario],known)
    # These are declared simulation-reference supports, not data-chosen regions.
    for (scenario,trial),(truth,x,y) in records.items():
        if scenario in (5,1005):add(scenario,trial,truth,x,y,18,22,'fixed_18_22',names[scenario],20)
        if scenario in (7,1007):add(scenario,trial,truth,x,y,24,52,'fixed_24_52',names[scenario])
    t=np.arange(41,dtype=float)-20
    logmean=math.log(1000)+.1*t+np.where(t<0,-.025*t*t,-.08*t*t)
    rng=np.random.default_rng(20261105)
    for trial in range(300):
        y=np.exp(logmean+rng.normal(0,.1,41)).astype(np.float32).astype(float)
        add(4010,trial,[20.625],range(41),y,0,40,'whole_input','asymmetric_log_gaussian',20)
        add(4011,trial,[19.375],range(41),y[::-1],0,40,'whole_input','asymmetric_reflected',20)
    for scenario,factor in ((4001,1.0),(4002,.01),(4003,1e4)):
        y=(np.exp(logmean)*factor).astype(np.float32).astype(float)
        add(scenario,0,[20.625],range(41),y,0,40,'whole_input','asymmetric_noiseless_unit_'+str(factor),20)
    return cases


def run(folder,output,np,plots=True):
    output.mkdir(parents=True,exist_ok=True);cases=cases_from_study(folder,np)
    # Serialize every actual fitting support, including explicitly unsupported
    # cases, so no selected or transformed input has to be inferred later.
    with (output/'common-inputs.bin').open('wb') as f:
        f.write(b'QENS0001')
        for c in cases:
            n=len(c['y']);truth=[p-c['start'] for p in c['truth'] if c['start']<=p<=c['end']]
            f.write(struct.pack('<4I',c['case'],0,n,len(truth)))
            for values in (truth,c['x'],c['y']):f.write(struct.pack('<'+'f'*len(values),*values))
            f.write(struct.pack('<'+'H'*n,*range(1,n+1)))
    (output/'supports.json').write_text(json.dumps([{k:v for k,v in c.items() if k not in ('x','y')} for c in cases],indent=2),encoding='utf-8')
    (output/'generated-model.json').write_text(json.dumps(dict(seed=20261105,trials=300,log_noise_sd=.1,n=41,hinge=20,b0=math.log(1000),b1_per_sample=.1,b2_per_sample2=-.025,b3_per_sample2=-.08,reflections_are_paired=True,noiseless_intensity_factors=[1,.01,1e4],interpretation='Specified model controls, not new acceptance thresholds. Fixed slice protocols use known synthetic layouts; they are not automatic support selection.'),indent=2),encoding='utf-8')
    cache={};results=[];metric_rows=[];plot_cases={}
    profile_file=(output/'hinge-profiles.csv').open('w',newline='',encoding='utf-8');pw=csv.writer(profile_file)
    pw.writerow(['case','hinge','rss','b0','b1','b2','b3','geometric_apex_index'])
    curve_file=(output/'fit-curves.csv').open('w',newline='',encoding='utf-8');cw=csv.writer(curve_file)
    cw.writerow(['case','original_index','x','observed','fitted_log_median','fitted_median','log_residual'])
    for c in cases:
        y=np.asarray(c['y'],dtype=float);n=len(y);base={k:c[k] for k in ('case','scenario','trial','name','variant','start','end')}
        nonpositive=int(sum(y<=0))
        if not np.isfinite(y).all():raise ValueError('Nonfinite observations cannot be omitted.')
        if nonpositive or n<5:
            results.append(dict(base,status='nonpositive_observations' if nonpositive else 'too_short',n=n,nonpositive=nonpositive,minimizers=0,hinge=None,peak_position=None,rss=None,b0=None,b1=None,b2=None,b3=None));continue
        if n not in cache:cache[n]=ProfileGeometry(n,np)
        geometry=cache[n];result,data=summarize_fit(y,geometry);z,beta,rss,curves,index=data
        if result['hinge'] is not None:result['hinge']+=c['start']
        if result['peak_position'] is not None:result['peak_position']+=c['start']
        results.append(dict(base,**result))
        if c['trial']==0:
            for j,h in enumerate(geometry.hinges):
                position=geometry.apex(beta[j],int(h))
                pw.writerow([c['case'],c['start']+int(h),float(rss[j]),*map(float,beta[j]),'' if position is None else c['start']+position])
            if index is not None or result['status']=='constant_data':
                fitted=curves[index] if index is not None else np.full(n,z[0])
                for i,pred in enumerate(fitted):cw.writerow([c['case'],c['start']+i,c['x'][i],y[i],float(pred),float(np.exp(pred)),float(z[i]-pred)])
            if (c['scenario'],c['variant']) in ((4010,'whole_input'),(7,'fixed_24_52'),(9,'whole_input'),(3000,'whole_input')):plot_cases[(c['scenario'],c['variant'])]=(c,result,data)
        if len(c['truth'])==1 and c['scenario'] in (4,1004,5,1005,3001,4010,4011):
            truth=c['truth'][0];metric_rows.append(dict(scenario=c['scenario'],variant=c['variant'],method='selected_hinge_four_coefficients',error=None if result['peak_position'] is None else result['peak_position']-truth))
            h=c['known_hinge_global']
            if h is not None and c['start']<h<c['end']:
                j=int(h-c['start'])-1;position=geometry.apex(beta[j],int(h-c['start']))
                metric_rows.append(dict(scenario=c['scenario'],variant=c['variant'],method='fixed_hinge_four_coefficients',error=None if position is None else c['start']+position-truth))
            t=(np.arange(n)-(n-1)/2)/(n-1);X=np.column_stack([np.ones(n),t,t*t]);coef=np.linalg.lstsq(X,z,rcond=None)[0]
            position=(n-1)/2-(n-1)*coef[1]/(2*coef[2]) if coef[2]<0 else None
            if position is not None and not 0<position<n-1:position=None
            metric_rows.append(dict(scenario=c['scenario'],variant=c['variant'],method='single_quadratic_three_coefficients_comparator',error=None if position is None else c['start']+position-truth))
    profile_file.close();curve_file.close()
    with (output/'case-results.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(results[0]));w.writeheader();w.writerows(results)
    metrics=[]
    for key in sorted({(r['scenario'],r['variant'],r['method']) for r in metric_rows}):
        rows=[r for r in metric_rows if (r['scenario'],r['variant'],r['method'])==key];errors=[r['error'] for r in rows if r['error'] is not None]
        metrics.append(dict(scenario=key[0],variant=key[1],method=key[2],cases=len(rows),geometric_apex_available=len(errors),bias=float(np.mean(errors)) if errors else None,rmse=float(np.sqrt(np.mean(np.square(errors)))) if errors else None))
    scope_counts=[]
    for scenario,variant in sorted({(r['scenario'],r['variant']) for r in results}):
        rows=[r for r in results if (r['scenario'],r['variant'])==(scenario,variant)]
        scope_counts.append(dict(scenario=scenario,variant=variant,cases=len(rows),status_counts=dict(Counter(r['status'] for r in rows)),geometric_apex_available=sum(r['peak_position'] is not None for r in rows)))
    report=dict(cases=len(cases),results=len(results),scope_counts=scope_counts,status_counts=dict(Counter(r['status'] for r in results)),geometry_cache_sizes=sorted(cache),production_changed=False,automatic_support_selection=False,statistical_peak_decisions=False,calibrated_intervals_reported=False,metrics=metrics,model='One shared piecewise-quadratic LOG mean, discrete hinge on all full-rank interior sample positions, unrestricted coefficients; classify the global LS winner geometrically. Predictions are fitted lognormal medians, not raw-intensity means.',covariance_scope='Fixed-support original observations counted once. No independence of overlapping fits assumed. Profile hinge selection and support selection do not inherit fixed-model coverage.')
    (output/'metrics.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    if plots:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig,axes=plt.subplots(2,2,figsize=(12,7),layout='constrained')
        for ax,key,title in zip(axes.flat,[(4010,'whole_input'),(7,'fixed_24_52'),(9,'whole_input'),(3000,'whole_input')],['Asymmetrischer Einzelpeak','Überlappung: Einzelmodell reicht nicht','pt_04: drei Strukturen, ein Modell','Konstante Kontrolle']):
            c,result,data=plot_cases[key];z,beta,rss,curves,index=data;x=np.arange(c['start'],c['end']+1)
            ax.plot(x,c['y'],label='Messwerte',color='#23546e')
            if index is not None:ax.plot(x,np.exp(curves[index]),label='Gemeinsames Einzelmodell',color='#b14c30')
            elif result['status']=='constant_data':ax.plot(x,c['y'],ls='--',label='Konstante, keine Peakform',color='#b14c30')
            ax.set_title(title);ax.set_xlabel('Messpunktindex');ax.set_ylabel('Intensität');ax.legend(fontsize=8)
        fig.suptitle('Referenzschätzung auf vorgegebenem Bereich – keine automatische Peakfreigabe\nAlle Übergangsstellen geprüft; jeder Messpunkt einmal. Kurven zeigen gefittete Medianintensitäten.',fontsize=11)
        fig.savefig(output/'common-profile.png',dpi=150);fig.savefig(output/'common-profile.pdf');plt.close(fig)
    return report

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('study',type=Path);p.add_argument('output',type=Path);p.add_argument('--numpy-root',type=Path);p.add_argument('--no-plots',action='store_true');a=p.parse_args()
    if a.numpy_root:sys.path.insert(0,str(a.numpy_root.resolve()))
    import numpy as np
    r=run(a.study,a.output,np,not a.no_plots);print(json.dumps(dict(cases=r['cases'],status_counts=r['status_counts'],metrics=[m for m in r['metrics'] if m['scenario'] in (4,4010,4011)]),indent=2))
if __name__=='__main__':main()
