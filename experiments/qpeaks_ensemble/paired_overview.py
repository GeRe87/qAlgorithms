"""Verify the matched-domain control and display paired descriptive results."""
import argparse,csv,json,math,itertools,struct,sys
from pathlib import Path
from collections import Counter
from analyse import read_inputs
from compare_bundle_cascade import fields,rows
from bundle_method import geometric_peak

def run(inputs,positive_native,matched_native,positive,matched,out,np):
    out.mkdir(parents=True,exist_ok=True);source=read_inputs(inputs/'inputs.bin')
    for name in ('cascade-runs.csv','cascade-peaks.csv'):
        a=[{k:v for k,v in r.items() if k!='seconds'} for r in rows(positive_native/name)]
        b=[{k:v for k,v in r.items() if k!='seconds'} for r in rows(matched_native/name)]
        assert a==b,name
    floor=json.loads((matched_native/'field-policy.json').read_text())['minimum_intensity']
    count=0;full_total=0;matched_total=0;single_windows=[]
    for a,b in itertools.zip_longest(fields(positive_native/'fields.bin'),fields(matched_native/'fields.bin')):
        assert a is not None and b is not None and a[0]==b[0]
        identity,z,full,_=a;_,zz,restricted,_=b;assert struct.pack('<'+'f'*len(z),*z)==struct.pack('<'+'f'*len(zz),*zz)
        y=source[identity][2];bad=[0]
        for v in y:bad.append(bad[-1]+int(v<=floor))
        key=lambda r:(r['center'],r['scale'])
        value=lambda r:tuple(r[k] for k in ('b0','b1','b2','b3'))
        expected={key(r):value(r) for r in full if bad[r['center']+r['scale']+1]==bad[r['center']-r['scale']]}
        actual={key(r):value(r) for r in restricted};assert expected==actual,identity
        if len(restricted)==1:
            r=restricted[0];scale=r['scale'];peak=geometric_peak(value(r),scale,2*scale+1)
            if peak['apex'] is not None and peak['above_both_support_edges'] is True:
                single_windows.append(dict(scenario=identity[0],trial=identity[1],native_center=r['center'],scale=scale,native_coefficients=list(value(r)),native_geometric_apex=r['center']-scale+peak['apex']))
        count+=1;full_total+=len(full);matched_total+=len(restricted)
    assert count==len(source)
    pos={(r['scenario'],r['trial']):r for r in json.loads((positive/'case-results.json').read_text())}
    mat={(r['scenario'],r['trial']):r for r in json.loads((matched/'case-results.json').read_text())};assert set(pos)==set(mat)==set(source)
    common=[]
    for identity,a in pos.items():
        b=mat[identity];truth=source[identity][0]
        if a['reference_kind']!='generated_reference' or len(truth)!=1:continue
        if not(a['strict_descriptions']==1 and a['strict_collectives']==0 and a['geometry_only_proposals']==b['geometry_only_proposals']==1 and a['junction_proposals']==b['junction_proposals']==0 and not a['fit_errors'] and not b['fit_errors']):continue
        target=truth[0];common.append(dict(scenario=identity[0],trial=identity[1],reference=target,cascade_error=a['strict_positions_axis'][0]-target,positive_error=a['ensemble_proposal_apices'][0]-target,matched_error=b['ensemble_proposal_apices'][0]-target))
    groups=[]
    for s in sorted({k[0] for k in source}):
        rs=[r for k,r in pos.items() if k[0]==s];ms=[mat[(r['scenario'],r['trial'])] for r in rs]
        entry=dict(scenario=s,name=rs[0]['name'],cases=len(rs),reference_kind=rs[0]['reference_kind'],reference_count=rs[0]['reference_count'])
        for name,records,alarmkey in (('cascade_strict',rs,'strict_alarm'),('cascade_collective',rs,'collective_alarm'),('ensemble_positive',rs,'ensemble_alarm'),('ensemble_matched',ms,'ensemble_alarm')):
            values=[r[alarmkey] for r in records];entry[name]=dict(alarm=sum(v is True for v in values),no_alarm=sum(v is False for v in values),error=sum(v is None for v in values))
        entry['counts']=dict(cascade_strict={str(k):v for k,v in Counter(r['strict_descriptions'] for r in rs).items()},ensemble_positive={str(k):v for k,v in Counter(r['geometry_only_proposals'] for r in rs).items()},ensemble_matched={str(k):v for k,v in Counter(r['geometry_only_proposals'] for r in ms).items()})
        paired=[r for r in common if r['scenario']==s]
        entry['common_single_output_localization']=dict(cases=len(paired),**{name:None if not paired else math.sqrt(sum(r[name+'_error']**2 for r in paired)/len(paired)) for name in ('cascade','positive','matched')})
        exact=[(a,b) for a,b in zip(rs,ms) if a['reference_kind']=='generated_reference' and a['reference_count'] and a['strict_descriptions']==a['reference_count'] and b['geometry_only_proposals']==a['reference_count'] and b['junction_proposals']==0 and not b['fit_errors']]
        positions=sum(a['reference_count'] for a,b in exact)
        entry['matched_exact_output_localization']=dict(cases=len(exact),reference_positions=positions,cascade_rmse=None if not positions else math.sqrt(sum(a['strict_localization']['squared_error'] for a,b in exact)/positions),ensemble_matched_rmse=None if not positions else math.sqrt(sum(b['ensemble_proposal_localization']['squared_error'] for a,b in exact)/positions),conditioning='Both output the reference count, no ensemble junction proposals; subset matching without a tolerance. Not unconditional recognition accuracy.')
        groups.append(entry)
    missed_single=[r for r in single_windows if len(source[(r['scenario'],r['trial'])][0])==1 and mat[(r['scenario'],r['trial'])]['ensemble_alarm'] is False]
    report=dict(matched_control_cases_verified=count,positive_native_windows=full_total,matched_native_windows=matched_total,matched_windows_are_exact_unchanged_kernel_subset=True,cascade_outputs_identical_across_domain_runs=True,existing_cascade_floor=floor,groups=groups,common_single_output_cases=common,single_center_missed_generated_peak_counterexamples=missed_single,localization_condition='Same generated single-peak inputs where cascade and both ensemble variants each produce exactly one complete nonjunction output/proposal. Conditional localization, not end-to-end recognition.',statistical_ensemble_detector_completed=False)
    (out/'overview.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    import matplotlib
    matplotlib.use('Agg');import matplotlib.pyplot as plt
    g={r['name']:r for r in groups};colors=('#23546e','#b14c30','#408553')
    fig,axes=plt.subplots(2,2,figsize=(14,10),layout='constrained')
    names=['flat_log','flat_raw','zero_clipped_raw_24','zero_clipped_raw_2400'];labels=['Log-Rauschen','Roh-Rauschen','Nullbeschnitten 24','Nullbeschnitten 2400'];x=np.arange(4)
    for j,(key,label) in enumerate((('cascade_strict','Kaskade strikt'),('ensemble_positive','Ensemble: alle positiven Punkte'),('ensemble_matched','Ensemble: gleicher Intensitätsboden'))):
        a=[g[n][key]['alarm'] for n in names];e=[g[n][key]['error'] for n in names];xx=x+(j-1)*.25
        axes[0,0].bar(xx,a,width=.24,label=label,color=colors[j]);axes[0,0].bar(xx,e,bottom=a,width=.24,color='none',edgecolor='#8b2222',hatch='///',label='Verarbeitungsfehler' if j==0 else None)
    axes[0,0].set_xticks(x,labels,fontsize=8);axes[0,0].set_ylabel('Eingaben mit mindestens einer Ausgabe / einem Vorschlag');axes[0,0].set_title('Rauschkontrollen: jeweils 32 Eingaben; Fehler separat schraffiert');axes[0,0].legend(fontsize=7);axes[0,0].set_ylim(0,35)
    names=['gaussian_wide','gaussian_sparse'];x=np.arange(2)
    for j,(key,label) in enumerate((('cascade','Kaskade'),('positive','Ensemble: alle positiven Punkte'),('matched','Ensemble: gleicher Intensitätsboden'))):
        axes[0,1].bar(x+(j-1)*.25,[float('nan') if g[n]['common_single_output_localization'][key] is None else g[n]['common_single_output_localization'][key] for n in names],width=.24,color=colors[j],label=label)
    for i,name in enumerate(names):
        if not g[name]['common_single_output_localization']['cases']:
            axes[0,1].text(i,.05,'Kein gemeinsamer Lagevergleich:\nEnsemble bei gleichem Boden\nohne Vorschlag',ha='center',fontsize=8,transform=axes[0,1].get_xaxis_transform())
    axes[0,1].set_xlim(-.6,1.6)
    axes[0,1].set_xticks(x,[f"Breiter Einzelpeak (n={g[names[0]]['common_single_output_localization']['cases']})",f"Schmaler Einzelpeak (n={g[names[1]]['common_single_output_localization']['cases']})"]);axes[0,1].set_ylabel('RMSE der Lage, Messpunktindex');axes[0,1].set_title('Lage: identische Teilmenge mit genau einer Ausgabe aller drei Verfahren');axes[0,1].legend(fontsize=7)
    for ax,name,title in ((axes[1,0],'three_gaussians','Drei getrennte Komponenten'),(axes[1,1],'overlapping_pair','Zwei überlappende Komponenten')):
        entry=g[name]
        upper=max(int(k) for histogram in entry['counts'].values() for k in histogram if k!='None')
        for j,key in enumerate(('cascade_strict','ensemble_positive','ensemble_matched')):
            hist=entry['counts'][key];xs=list(range(upper+1));ys=[hist.get(str(k),0) for k in xs]
            ax.plot(xs,ys,'o-',label=('Kaskaden-Beschreibungen','Rohe Ensemble-Vorschläge','Ensemble-Vorschläge, gleicher Boden')[j],color=colors[j])
        ax.axvline(entry['reference_count'],color='#555555',ls='--',lw=.8,label='Generierte Komponentenzahl');ax.set_xlabel('Anzahl Ausgaben bzw. ungeprüfter geometrischer Vorschläge');ax.set_ylabel('Eingaben');ax.set_title(title+' – jeweils 32 Eingaben');ax.legend(fontsize=7)
    fig.suptitle('Gepaarter Pilot: unveränderte Testkaskade gegen geometrischen Ensemble-Rohansatz\nEnsemble-Vorschläge sind keine statistisch bestätigten Peaks. Lagegewinn und Erkennungsleistung getrennt bewerten.',fontsize=12)
    fig.savefig(out/'paired-comparison.png',dpi=145);fig.savefig(out/'paired-comparison.pdf');plt.close(fig)
    return report

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('inputs','positive_native','matched_native','positive','matched','output'):p.add_argument(k,type=Path)
    p.add_argument('--numpy-root',type=Path);a=p.parse_args()
    if a.numpy_root:sys.path.insert(0,str(a.numpy_root.resolve()))
    import numpy as np
    r=run(a.inputs,a.positive_native,a.matched_native,a.positive,a.matched,a.output,np);print(json.dumps({k:v for k,v in r.items() if k not in ('groups','common_single_output_cases')},indent=2))
