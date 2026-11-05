"""Account for every detailed native fit; pilot ensemble-derived shared fits.

The pilot fits every tied largest bundle per input, not every bundle. Native
fields are copied without consulting current_accepted or the old feature lists.
"""
import argparse,csv,json,hashlib,itertools,sys
from pathlib import Path
from collections import defaultdict,Counter
from bundle_method import topology,bundle_support,MomentGeometry,geometric_peak

def rows(path):
    with path.open(encoding='utf-8',newline='') as f:yield from csv.DictReader(f)

def run(collection,out,np,observation_file=None):
    out.mkdir(parents=True,exist_ok=True)
    observations=defaultdict(list)
    observation_file=observation_file or collection/'observations.csv'
    source_observations=list(rows(collection/'observations.csv'))
    supplied_observations=list(rows(observation_file))
    assert len(source_observations)==len(supplied_observations)
    for a,b in zip(source_observations,supplied_observations):
        assert all(a[k]==b[k] for k in ('scenario','trial','index','x','y','cumulative_df'))
    for r in supplied_observations:observations[(int(r['scenario']),int(r['trial']))].append(r)
    traces={(int(r['scenario']),int(r['trial'])):r for r in rows(collection/'traces.csv') if r['detail']=='1'}
    assert set(traces)==set(observations)
    geometry=MomentGeometry(np);summaries=[];all_fits=[];figdata={}
    files={name:(out/(name+'.csv')).open('w',newline='',encoding='utf-8') for name in ('native-membership','basins','edges','bundles','hinge-fits','fit-curves')}
    writers={k:csv.writer(f) for k,f in files.items()}
    writers['native-membership'].writerow(['scenario','trial','id','block','center','scale','b0','b1','b2','b3','basin','bundle'])
    writers['basins'].writerow(['scenario','trial','id','block','scale','kind','root_left','root_right','first_center','last_center','fits'])
    writers['edges'].writerow(['scenario','trial','from_basin','to_basin','shared_centers'])
    writers['bundles'].writerow(['scenario','trial','bundle','basins','native_fits','contains_merge','contains_split','fit_in_pilot','reason'])
    writers['hinge-fits'].writerow(['scenario','trial','bundle','hinge_global','rss','fast_rss','b0','b1','b2','b3'])
    writers['fit-curves'].writerow(['scenario','trial','bundle','index','original_x','observed_y','native_log_y','fitted_log'])
    seen=set();max_coefficient_score_error=0.;max_prefix_vs_svd_prediction_error=0.
    key=lambda r:(int(r['scenario']),int(r['trial']))
    for identity,iterator in itertools.groupby(rows(collection/'candidates.csv'),key):
        assert identity not in seen;seen.add(identity)
        candidates=list(iterator);obs=observations[identity];trace=traces[identity]
        assert len(candidates)==int(trace['all_regressions'])
        assert [int(r['index']) for r in obs]==list(range(len(obs)))
        bs,edges,groups,membership=topology(candidates);by_basin={b:g['id'] for g in groups for b in g['basins']}
        for r in candidates:
            owner=membership[int(r['id'])]
            writers['native-membership'].writerow([*identity,*[r[k] for k in ('id','block','center','scale','b0','b1','b2','b3')],owner,by_basin.get(owner,'')])
        for b in bs:
            writers['basins'].writerow([*identity,b['id'],b['block'],b['scale'],b['kind'],*(['',''] if b['root'] is None else b['root']),min(b['positions']),max(b['positions']),len(b['candidate_ids'])])
        for (a,b),count in sorted(edges.items()):writers['edges'].writerow([*identity,a,b,count])
        most=max((g['native_fits'] for g in groups),default=0);fit_count=0;errors=0
        for g in groups:
            selected=g['native_fits']==most
            writers['bundles'].writerow([*identity,g['id'],len(g['basins']),g['native_fits'],g['contains_merge'],g['contains_split'],selected,'all_tied_largest_bundles_pilot' if selected else 'not_fitted_in_this_pilot'])
            if not selected:continue
            record=dict(scenario=identity[0],trial=identity[1],bundle=g['id'],native_fits=g['native_fits'],contains_merge=g['contains_merge'],contains_split=g['contains_split'],statistical_peak_release=False)
            try:
                left,right,hinges=bundle_support(g,bs,candidates);subset=obs[left:right+1]
                if len(subset)!=right-left+1 or any(not r['log_y_float'] for r in subset):raise ValueError('Missing native log observation in ensemble support.')
                z=np.array([float(r['log_y_float']) for r in subset]);fits,winners=geometry.fit(z,hinges)
                for r in fits:
                    writers['hinge-fits'].writerow([*identity,g['id'],left+r['hinge'],r['rss'],r['fast_rss'],*r['coefficients']])
                    max_coefficient_score_error=max(max_coefficient_score_error,abs(r['rss']-r['fast_rss']))
                    # Independent direct normalized SVD reference, every hinge.
                    t=(np.arange(len(z))-r['hinge'])/(len(z)-1);X=np.column_stack([np.ones(len(z)),t,np.where(t<0,t*t,0),np.where(t>0,t*t,0)])
                    norm=np.linalg.norm(X,axis=0);Q=X/norm
                    pred=Q@np.linalg.lstsq(Q,z,rcond=None)[0]
                    p=r['coefficients'];u=np.arange(len(z))-r['hinge'];got=p[0]+p[1]*u+np.where(u<0,p[2],p[3])*u*u
                    error=float(np.max(np.abs(pred-got)));max_prefix_vs_svd_prediction_error=max(max_prefix_vs_svd_prediction_error,error)
                    tolerance=4096*np.finfo(float).eps*len(z)*max(1.,float(np.linalg.norm(z)))*float(np.linalg.cond(Q))**2
                    assert error<=tolerance,(identity,g['id'],r['hinge'],error,tolerance)
                collective=g['contains_merge'] or g['contains_split']
                best=min(fits,key=lambda r:r['rss']);record.update(status='collective_diagnostic_fit' if collective else 'conditional_single_bundle_fit',left=left,right=right,hinge_proposals=len(hinges),equal_numeric_minimizers=len(winners),selected_hinge=None if len(winners)!=1 else left+best['hinge'],rss=best['rss'])
                if len(winners)==1:
                    peak=geometric_peak(best['coefficients'],best['hinge'],len(z))
                    record.update(geometric_apex=None if peak['apex'] is None else left+peak['apex'],fitted_log_height=peak['log_height'],above_both_support_edges=peak['above_both_support_edges'])
                    p=best['coefficients'];u=np.arange(len(z))-best['hinge'];curve=p[0]+p[1]*u+np.where(u<0,p[2],p[3])*u*u
                    for i,(r,v) in enumerate(zip(subset,curve)):writers['fit-curves'].writerow([*identity,g['id'],left+i,r['x'],r['y'],r['log_y_float'],v])
                    if identity not in figdata:figdata[identity]=(obs,bs,by_basin,(left,right,curve,collective))
                fit_count+=1
            except Exception as e:
                record.update(status='explicit_fit_error',error=str(e));errors+=1
            all_fits.append(record)
        figdata.setdefault(identity,(obs,bs,by_basin,None))
        indegrees=Counter(b for a,b in edges)
        summaries.append(dict(scenario=identity[0],trial=identity[1],n=len(obs),native_fits=len(candidates),maximum_basins=sum(b['kind']=='maximum' for b in bs),nonmaximum_basins=sum(b['kind']!='maximum' for b in bs),scale_edges=len(edges),topological_bundles=len(groups),merge_basins=sum(count>1 for count in indegrees.values()),fit_pilot_completed=fit_count,fit_pilot_errors=errors))
        print(json.dumps(summaries[-1]),flush=True)
    assert seen==set(traces)
    for f in files.values():f.close()
    (out/'fits.json').write_text(json.dumps(all_fits,indent=2),encoding='utf-8')
    report=dict(detailed_inputs=len(seen),native_fits_accounted=sum(r['native_fits'] for r in summaries),scope='All detailed native trial-zero inputs; other native trials are outside this topology pilot, not silently failed inputs.',summaries=summaries,shared_fit_pilot='All tied largest bundles by member count per detailed input, including junction bundles; this is a diagnostic workload selection, not a peak acceptance rule.',shared_fit_results=len(all_fits),shared_fit_errors=sum(r['status']=='explicit_fit_error' for r in all_fits),max_fast_rss_vs_explicit_difference=max_coefficient_score_error,max_prefix_vs_independent_svd_prediction_error=max_prefix_vs_svd_prediction_error,geometry_cache_entries=len(geometry.cache),old_cascade_used_for_grouping_or_fits=False,existing_feature_windows_used=False,physical_peak_count_decision=False,numerical_slope_signs_certified=False,source={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (collection/'candidates.csv',collection/'observations.csv',collection/'traces.csv')})
    report['source']['fitting_observations_sha256']=hashlib.sha256(observation_file.read_bytes()).hexdigest()
    report['fitting_observations_bound_to_native_input']=True
    (out/'report.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    import matplotlib
    matplotlib.use('Agg');import matplotlib.pyplot as plt
    examples=[k for k in ((4,0),(7,0),(8,0),(9,0),(2000,0),(3000,0)) if k in figdata]
    fig,axes=plt.subplots(len(examples),2,figsize=(13,3*len(examples)),layout='constrained',squeeze=False)
    for ax,k in zip(axes,examples):
        obs,bs,by_basin,fit=figdata[k];x=[int(r['index']) for r in obs];y=[float(r['y']) for r in obs]
        ax[0].plot(x,y,label='Originale Beobachtungen',color='#23546e')
        if fit:
            left,right,curve,collective=fit;ax[0].plot(range(left,right+1),np.exp(curve),label='Kollektiv-Diagnosefit' if collective else 'Fit aus Pilotbündel',color='#b14c30')
        ax[0].set_title(f"Szenario {k[0]}: bedingter Fit, keine Peakfreigabe");ax[0].set_xlabel('Messpunktindex');ax[0].set_ylabel('Intensität');ax[0].legend(fontsize=8)
        modes=[b for b in bs if b['kind']=='maximum']
        ax[1].scatter([sum(b['root'])/2 for b in modes],[b['scale'] for b in modes],c=[by_basin[b['id']]%20 for b in modes],cmap='tab20',s=8)
        ax[1].set_xlim(0,len(obs)-1);ax[1].set_xlabel('Maximum-Korridor, Mittelpunkt nur zur Darstellung');ax[1].set_ylabel('Skala');ax[1].set_title('Bündel aus nativen Steigungen; Farben wiederholen sich')
    fig.suptitle('Methodenprototyp: native Running-Regressionen → Bündelgraph → gemeinsamer Piecewise-Fit\nAlle Rauschzweige bleiben erhalten; größte Bündel nur als deklarierte Fit-Stichprobe.',fontsize=11)
    fig.savefig(out/'bundle-method.png',dpi=140);plt.close(fig)
    return report

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('collection',type=Path);p.add_argument('output',type=Path);p.add_argument('--numpy-root',type=Path);p.add_argument('--observations',type=Path);a=p.parse_args()
    if a.numpy_root:sys.path.insert(0,str(a.numpy_root.resolve()))
    import numpy as np
    report=run(a.collection,a.output,np,a.observations);print(json.dumps({k:v for k,v in report.items() if k not in ('summaries','source')},indent=2))
    if report['shared_fit_errors']:raise SystemExit(1)
if __name__=='__main__':main()
