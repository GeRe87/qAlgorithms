"""Paired geometry-only ensemble ablation versus the unchanged native cascade.

All bundles are fitted. A complete geometric proposal is NOT a statistically
accepted peak. No fitted threshold, oracle selection or distance cutoff.
"""
import argparse,csv,json,struct,sys,time,math,hashlib
from pathlib import Path
from collections import defaultdict,Counter
from analyse import read_inputs
from bundle_method import topology,bundle_support,MomentGeometry,geometric_peak

def rows(path):
    with path.open(encoding='utf-8',newline='') as f:yield from csv.DictReader(f)

def fields(path):
    with path.open('rb') as f:
        if f.read(8)!=b'QFLD0001':raise ValueError('Unknown native field format.')
        seen=set()
        def read(n):
            data=f.read(n)
            if len(data)!=n:raise ValueError('Truncated native field file.')
            return data
        while True:
            first=f.read(4)
            if not first:break
            if len(first)!=4:raise ValueError('Truncated native field header.')
            s,t,n,blocks=struct.unpack('<4I',first+read(12));identity=s,t
            if identity in seen or n<1 or n>65535 or blocks>n:raise ValueError('Invalid/duplicate native field record.')
            seen.add(identity);z=struct.unpack('<'+'f'*n,read(4*n));candidates=[];intervals=[]
            for b in range(blocks):
                start,length,count=struct.unpack('<3I',read(12));maximum=min(215,(length-1)//2)
                expected=sum(length-2*k for k in range(2,maximum+1))
                if length<5 or start+length>n or expected!=count:raise ValueError('Native field grid accounting mismatch.')
                if intervals and start<=intervals[-1][1]:raise ValueError('Native positive blocks overlap.')
                intervals.append((start,start+length-1));keys=set()
                for _ in range(count):
                    center,scale,b0,b1,b2,b3=struct.unpack('<HH4d',read(36))
                    key=center,scale
                    if key in keys or not 2<=scale<=maximum or not start+scale<=center<start+length-scale or not all(math.isfinite(v) for v in (b0,b1,b2,b3)):raise ValueError('Invalid native regression.')
                    keys.add(key);candidates.append(dict(id=len(candidates),block=b,center=center,scale=scale,b0=b0,b1=b1,b2=b2,b3=b3))
            yield identity,z,candidates,intervals

def match_positions(truth,predictions):
    """Minimum squared-distance one-to-one subset match, without a radius.

    It matches all points in the smaller set; therefore this is a localization
    diagnostic, never a recall/precision definition. Unmatched counts stay visible.
    """
    truth=sorted(map(float,truth));predictions=sorted(map(float,predictions));swap=len(truth)>len(predictions)
    small,large=(predictions,truth) if swap else (truth,predictions);m,n=len(small),len(large)
    cost=[[math.inf]*(n+1) for _ in range(m+1)];take=[[False]*(n+1) for _ in range(m+1)];cost[0]=[0.]*(n+1)
    for i in range(1,m+1):
        for j in range(i,n+1):
            matched=cost[i-1][j-1]+(small[i-1]-large[j-1])**2
            if matched<cost[i][j-1]:cost[i][j]=matched;take[i][j]=True
            else:cost[i][j]=cost[i][j-1]
    i,j=m,n;pairs=[]
    while i:
        if take[i][j]:
            pairs.append((large[j-1],small[i-1]) if swap else (small[i-1],large[j-1]));i-=1
        j-=1
    pairs.reverse()
    return dict(pairs=pairs,errors=[p-t for t,p in pairs],unmatched_reference=len(truth)-m,unmatched_outputs=len(predictions)-m,squared_error=cost[m][n])

def run(inputs,native,out,np,topology_version='center_slope_v1',estimator='refit'):
    if topology_version not in ('center_slope_v1','model_boundary_v2','model_singleton_v3'):raise ValueError('Unknown topology version.')
    if estimator not in ('refit','native-representative'):raise ValueError('Unknown estimator.')
    if estimator=='native-representative':
        from native_representative import native_bundle
    out.mkdir(parents=True,exist_ok=True)
    source=read_inputs(inputs/'inputs.bin');cases={(c['scenario'],c['trial']):c for c in json.loads((inputs/'cases.json').read_text())}
    assert set(source)==set(cases)
    cascade={(int(r['scenario']),int(r['trial']),r['policy']):r for r in rows(native/'cascade-runs.csv')}
    timings={(int(r['scenario']),int(r['trial'])):r for r in rows(native/'field-timing.csv')}
    peakrows=defaultdict(list)
    for r in rows(native/'cascade-peaks.csv'):peakrows[(int(r['scenario']),int(r['trial']),r['policy'])].append(r)
    geometry=MomentGeometry(np);results=[];bundle_results=[];seen=set();timing_rows=[]
    score_file=(out/'hinge-scores.bin').open('wb') if estimator=='refit' else None
    if score_file is not None:score_file.write(b'QSCR0001')
    score_count=0;max_fast_error=0.
    for identity,z,candidates,intervals in fields(native/'fields.bin'):
        assert identity not in seen;seen.add(identity);s,t=identity;c=cases[identity];truth,x,y=source[identity]
        assert len(z)==len(y) and len(candidates)==int(timings[identity]['regressions'])
        assert all(math.isfinite(v) if yy>0 else math.isnan(v) for yy,v in zip(y,z))
        began=time.perf_counter();bs,edges,groups,membership=topology(candidates,model_boundaries=topology_version!='center_slope_v1',log_values=z,boundary_scope='singleton_blocks' if topology_version=='model_singleton_v3' else 'all_edges');by_id={int(r['id']):r for r in candidates}
        assert len(membership)==len(candidates)
        proposals=[];errors=0;local_records=[]
        for g in groups:
            record=dict(scenario=s,trial=t,bundle=g['id'],native_fits=g['native_fits'],contains_merge=g['contains_merge'],contains_split=g['contains_split'],eligible_geometric_proposal=False,statistically_accepted=False)
            if topology_version!='center_slope_v1':record['boundary_model_nodes']=len(g['boundary_model_nodes'])
            try:
                left,right,hinges=bundle_support(g,bs,candidates,by_id)
                if estimator=='native-representative':
                    members=[by_id[i] for basin in g['basins'] for i in bs[basin]['candidate_ids']]
                    best,rss,ties,coefficients,peak=native_bundle(z[left:right+1],left,members)
                    center,scale=best['center'],best['scale'];apex_origin=center-scale
                    record.update(left=left,right=right,numeric_minimizers=ties,best_rss=rss)
                    if ties==1:
                        record.update(hinge_global=center,coefficients=coefficients,native_candidate_id=best['id'],
                                      native_left=apex_origin,native_right=center+scale,
                                      extrapolated_points=(right-left+1)-(2*scale+1))
                else:
                    local=np.array(z[left:right+1]);fits,winners=geometry.fit(local,hinges)
                    for fit in fits:
                        score_file.write(struct.pack('<4I2d',s,t,g['id'],left+fit['hinge'],fit['rss'],fit['fast_rss']));score_count+=1
                        max_fast_error=max(max_fast_error,abs(fit['rss']-fit['fast_rss']))
                    best=min(fits,key=lambda r:r['rss']);ties=len(winners);apex_origin=left
                    record.update(left=left,right=right,hinge_proposals=len(hinges),numeric_minimizers=ties,best_rss=best['rss'],hinge_global=left+best['hinge'],coefficients=best['coefficients'])
                    peak=geometric_peak(best['coefficients'],best['hinge'],len(local)) if ties==1 else None
                if ties!=1:record['status']='multiple_equal_minimizers'
                else:
                    apex=peak['apex']
                    complete=apex is not None and peak['above_both_support_edges'] is True
                    record.update(status='closed_geometric_proposal' if complete else 'no_closed_geometric_peak',eligible_geometric_proposal=complete,apex=None if apex is None else apex_origin+apex,log_height=peak['log_height'],above_both_support_edges=peak['above_both_support_edges'])
                    if complete:proposals.append(record)
            except Exception as e:
                errors+=1;record.update(status='explicit_fit_error',error=str(e))
            local_records.append(record)
        elapsed=time.perf_counter()-began;bundle_results.extend(local_records)
        base=dict(scenario=s,trial=t,name=c['name'],reference_kind=c['reference_kind'],reference_count=None if c['reference_kind']=='unknown_real' else len(truth),native_fits=len(candidates),atomic_bundles=len(groups),all_bundle_fits=len(local_records),fit_errors=errors,geometry_only_proposals=len(proposals),junction_proposals=sum(r['contains_merge'] or r['contains_split'] for r in proposals),ensemble_status='complete_ablation' if not errors else 'incomplete_fit_errors',ensemble_alarm=None if errors else bool(proposals),ensemble_proposal_apices=[r['apex'] for r in proposals])
        for policy in ('strict','collective'):
            key=s,t,policy;runrow=cascade[key];pr=peakrows[key];failed=bool(runrow['error']) or int(runrow['status'])<0
            assert len(pr)==(0 if failed else int(runrow['description_count']))
            base[policy+'_status']='explicit_error' if failed else 'complete'
            base[policy+'_descriptions']=None if failed else len(pr)
            base[policy+'_collectives']=None if failed else int(runrow['collective_count'])
            base[policy+'_alarm']=None if failed else bool(pr)
            base[policy+'_positions_axis']=[float(r['position_axis']) for r in pr]
            base[policy+'_seconds']=float(runrow['seconds'])
            if c['reference_kind']!='unknown_real':
                assert tuple(x)==tuple(map(float,range(len(x))))
                base[policy+'_localization']=None if failed else match_positions(truth,[float(r['position_axis']) for r in pr if r['is_collective']=='0'])
        if c['reference_kind']!='unknown_real':base['ensemble_proposal_localization']=match_positions(truth,[r['apex'] for r in proposals])
        results.append(base);timing_rows.append(dict(scenario=s,trial=t,n=len(y),native_field_seconds=float(timings[identity]['kernel_seconds']),python_all_bundle_stage_seconds_including_score_logging=elapsed,strict_seconds=base['strict_seconds'],collective_seconds=base['collective_seconds']))
        print(f"case={len(results)} id={identity} fits={len(local_records)} proposals={len(proposals)} errors={errors} python_seconds={elapsed:.4f}",flush=True)
    if score_file is not None:score_file.close()
    assert seen==set(source)==set(timings)
    (out/'case-results.json').write_text(json.dumps(results,indent=2),encoding='utf-8')
    (out/'bundle-fits.json').write_text(json.dumps(bundle_results,indent=2),encoding='utf-8')
    with (out/'timing.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(timing_rows[0]));w.writeheader();w.writerows(timing_rows)
    grouped=[]
    for scenario in sorted({r['scenario'] for r in results}):
        rs=[r for r in results if r['scenario']==scenario];entry=dict(scenario=scenario,name=rs[0]['name'],reference_kind=rs[0]['reference_kind'],cases=len(rs),reference_count=rs[0]['reference_count'],median_atomic_bundles=float(np.median([r['atomic_bundles'] for r in rs])),median_geometry_proposals=float(np.median([r['geometry_only_proposals'] for r in rs])))
        for method in ('strict','collective','ensemble'):
            alarms=[r[method+'_alarm'] for r in rs];entry[method]=dict(alarm=sum(a is True for a in alarms),no_alarm=sum(a is False for a in alarms),error=sum(a is None for a in alarms))
        paired=[r for r in rs if r['strict_alarm'] is not None and r['ensemble_alarm'] is not None]
        entry['paired_screening']={f'cascade_{a}_ensemble_{b}':sum(r['strict_alarm']==a and r['ensemble_alarm']==b for r in paired) for a in (False,True) for b in (False,True)}
        both=[r for r in rs if r['reference_count']==1 and r['reference_kind']=='generated_reference' and r['strict_descriptions']==1 and r['strict_collectives']==0 and r['geometry_only_proposals']==1 and r['junction_proposals']==0 and not r['fit_errors']]
        entry['paired_single_output_localization']=dict(cases=len(both),strict_rmse=None if not both else math.sqrt(sum(r['strict_localization']['squared_error'] for r in both)/len(both)),ensemble_rmse=None if not both else math.sqrt(sum(r['ensemble_proposal_localization']['squared_error'] for r in both)/len(both)),conditioning='Both have exactly one complete noncollective output/proposal on the same generated single-peak inputs; not unconditional detector accuracy.')
        grouped.append(entry)
    report=dict(cases=len(results),native_regressions_accounted=sum(r['native_fits'] for r in results),all_bundles_fitted=len(bundle_results),fit_errors=sum(r['fit_errors'] for r in results),hinge_scores=score_count,max_fast_rss_difference=max_fast_error,groups=grouped,ensemble_variant='Geometry-only all-bundle ablation; no statistical acceptance and no oracle selection.',proposal_rule='Unique numeric best hinge, strict interior local maximum, fitted height strictly above both fitted support edges. All such atomic-bundle proposals retained, including merge/split bundles flagged as junctions. These are analytic shape conditions, not calibrated significance.',ensemble_final_physical_peak_counts_defined=False,screening_event='At least one complete geometric proposal versus at least one native output description. Errors are third outcomes, never zero detections.',localization_metric='Minimum-distance one-to-one subset assignment with no tolerance. Unmatched counts remain explicit. It is not recall; junction fits and extra proposals can be mixtures/duplicates.',areas_compared=False,whole_real_file_feature_lists_regenerated=False,timing_interpretation='Native cascade executes optimized C++; ensemble uses the same native kernel plus unoptimized Python graph/fit reference and score logging. Cold/warm cache depends on fixed case order. No algorithmic speed winner can be inferred.',inputs_sha256=hashlib.sha256((inputs/'inputs.bin').read_bytes()).hexdigest(),native_fields_sha256=hashlib.sha256((native/'fields.bin').read_bytes()).hexdigest())
    report['native_field_policy']=json.loads((native/'field-policy.json').read_text())
    report['geometry_cache_entries']=len(geometry.cache)
    if estimator=='native-representative':
        report['estimator']=estimator
        report['all_bundles_evaluated']=report.pop('all_bundles_fitted')
        report['ensemble_variant']='Geometry-only selection from existing native regressions. No coefficients fitted.'
        report['proposal_rule']='Unique minimum log RSS across the common bundle support, including any extrapolated points. Geometric apex and edge checks use the selected model’s original native window.'
        report['timing_interpretation']='Native cascade executes optimized C++; ensemble uses the same native kernel plus Python grouping and scoring. Timings do not predict C++ performance.'
    report['topology_version']=topology_version
    (out/'report.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    return report

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('inputs','native','output'):p.add_argument(k,type=Path)
    p.add_argument('--numpy-root',type=Path);p.add_argument('--topology',choices=('center_slope_v1','model_boundary_v2','model_singleton_v3'),default='center_slope_v1')
    p.add_argument('--estimator',choices=('refit','native-representative'),default='refit');a=p.parse_args()
    if a.numpy_root:sys.path.insert(0,str(a.numpy_root.resolve()))
    import numpy as np
    r=run(a.inputs,a.native,a.output,np,a.topology,a.estimator);print(json.dumps({k:v for k,v in r.items() if k!='groups'},indent=2))
    if r['fit_errors']:raise SystemExit(1)
