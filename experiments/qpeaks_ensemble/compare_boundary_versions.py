"""Frozen v1/v2 comparison with distinct randomized trial intervals."""
import argparse,json,sys,math,hashlib,struct
from pathlib import Path
from collections import Counter

def records(folder):return {(r['scenario'],r['trial']):r for r in json.loads((folder/'case-results.json').read_text())}
def metrics(rs1,rs2,revision_label='v2'):
    assert set(rs1)==set(rs2);groups=[]
    for s in sorted({k[0] for k in rs1}):
        pairs=[(rs1[k],rs2[k]) for k in sorted(rs1) if k[0]==s];a0=pairs[0][0]
        g=dict(scenario=s,name=a0['name'],cases=len(pairs),reference_kind=a0['reference_kind'],reference_count=a0['reference_count'])
        for label,items,key in (('cascade',[a for a,b in pairs],'strict_alarm'),('v1',[a for a,b in pairs],'ensemble_alarm'),(revision_label,[b for a,b in pairs],'ensemble_alarm')):
            alarm=[r[key] for r in items];counts=[r['strict_descriptions'] if label=='cascade' else r['geometry_only_proposals'] for r in items]
            g[label]=dict(alarm=sum(x is True for x in alarm),no_alarm=sum(x is False for x in alarm),error=sum(x is None for x in alarm),counts={str(k):v for k,v in Counter(counts).items()})
        exact=[(a,b) for a,b in pairs if a['reference_kind']=='generated_reference' and a['reference_count'] and a['strict_descriptions']==a['reference_count'] and b['geometry_only_proposals']==a['reference_count'] and b['junction_proposals']==0 and not b['fit_errors']]
        n=sum(a['reference_count'] for a,b in exact)
        g['cascade_revision_localization']=dict(cases=len(exact),positions=n,cascade_rmse=None if not n else math.sqrt(sum(a['strict_localization']['squared_error'] for a,b in exact)/n),revision_rmse=None if not n else math.sqrt(sum(b['ensemble_proposal_localization']['squared_error'] for a,b in exact)/n),conditioning='Both have the reference output count, no revised-model junctions. Conditional localization; not end-to-end accuracy.')
        common=[(a,b) for a,b in exact if a['geometry_only_proposals']==a['reference_count'] and a['junction_proposals']==0 and not a['fit_errors']]
        n=sum(a['reference_count'] for a,b in common)
        g['common_three_way_localization']=dict(cases=len(common),positions=n,cascade_rmse=None if not n else math.sqrt(sum(a['strict_localization']['squared_error'] for a,b in common)/n),v1_rmse=None if not n else math.sqrt(sum(a['ensemble_proposal_localization']['squared_error'] for a,b in common)/n),revision_rmse=None if not n else math.sqrt(sum(b['ensemble_proposal_localization']['squared_error'] for a,b in common)/n))
        groups.append(g)
    changes=[dict(scenario=k[0],trial=k[1],reference_kind=rs1[k]['reference_kind'],reference_count=rs1[k]['reference_count'],v1=rs1[k]['ensemble_alarm'],revision=rs2[k]['ensemble_alarm'],v1_proposals=rs1[k]['geometry_only_proposals'],revision_proposals=rs2[k]['geometry_only_proposals']) for k in sorted(rs1) if rs1[k]['ensemble_alarm']!=rs2[k]['ensemble_alarm']]
    return dict(groups=groups,changed_screening_outcomes=changes)

def run(reference_input,reference_v1,reference_v2,new_input,new_v1,new_v2,freeze,out,np,freeze_source=None):
    out.mkdir(parents=True,exist_ok=True);frozen=json.loads(freeze.read_text())
    for inp,comparison in ((reference_input,reference_v1),(reference_input,reference_v2),(new_input,new_v1),(new_input,new_v2)):
        digest=hashlib.sha256((inp/'inputs.bin').read_bytes()).hexdigest()
        assert json.loads((comparison/'report.json').read_text())['inputs_sha256']==digest
        assert json.loads((inp/'protocol.json').read_text())['input_sha256']==digest
    frozen_source=freeze_source or Path(__file__).parent
    for name,digest in frozen['method_files'].items():assert hashlib.sha256((frozen_source/name).read_bytes()).hexdigest()==digest
    revision_label=frozen['method'].rsplit('_',1)[-1]
    old=records(reference_v1);old2=records(reference_v2);new=records(new_v1);new2=records(new_v2)
    meta={(c['scenario'],c['trial']):c for c in json.loads((new_input/'cases.json').read_text())}
    additional={k for k,c in meta.items() if c['evaluation_role']=='trial_interval'}
    randomized={k[0] for k in additional};previous={k for k in old if k[0] in randomized}
    assert previous.isdisjoint(additional) and len(previous)==len(additional)==256
    lo,hi=frozen['new_trial_interval'];assert all(lo<=k[1]<hi for k in additional)
    for first,last in frozen.get('development_trial_intervals',[]):assert all(not first<=k[1]<last for k in additional)
    fixed=set(new)-additional;assert fixed<=set(old) and len(fixed)==24
    for label,one,two in (('reference',old,old2),('additional',new,new2)):
        assert set(one)==set(two)
        for k,a in one.items():
            b=two[k]
            assert all(a.get(n)==b.get(n) for n in ('name','reference_kind','reference_count','strict_descriptions','strict_collectives','strict_alarm','strict_positions_axis','collective_descriptions','collective_collectives','collective_alarm'))
    result=dict(revision_label=revision_label,method_source_frozen_before_additional_evaluation=True,method_freeze=frozen,previous_randomized_cases=len(previous),additional_randomized_cases=len(additional),trial_intervals_disjoint=True,fixed_diagnostics_reused=len(fixed),unique_input_identities=len(set(old)|set(new)),reference=metrics(old,old2,revision_label),additional=metrics({k:new[k] for k in additional},{k:new2[k] for k in additional},revision_label),fixed_diagnostics=metrics({k:new[k] for k in fixed},{k:new2[k] for k in fixed},revision_label),statistical_peak_acceptance=False,runtime_speed_claim=False)
    if revision_label=='v3':
        signature=lambda r:(r['left'],r['right'],r['hinge_global'],r['best_rss'],tuple(r['coefficients']),r['status'],r['contains_merge'],r['contains_split'])
        checks=[]
        for label,baseline,revised in (('reference',reference_v1,reference_v2),('additional',new_v1,new_v2)):
            a=json.loads((baseline/'bundle-fits.json').read_text());b=json.loads((revised/'bundle-fits.json').read_text())
            original=Counter((r['scenario'],r['trial'],signature(r)) for r in a)
            updated=Counter((r['scenario'],r['trial'],signature(r)) for r in b)
            assert not(original-updated),'An existing fitted bundle changed under the singleton-only correction.'
            def normalized_scores(folder,records):
                geometry={(r['scenario'],r['trial'],r['bundle']):(r['left'],r['right']) for r in records};values=Counter()
                with (folder/'hinge-scores.bin').open('rb') as f:
                    assert f.read(8)==b'QSCR0001'
                    while True:
                        raw=f.read(32)
                        if not raw:break
                        assert len(raw)==32
                        s,t,g,h,rss,fast=struct.unpack('<4I2d',raw);left,right=geometry[(s,t,g)]
                        values[(s,t,left,right,h,rss,fast)]+=1
                return values
            old_scores=normalized_scores(baseline,a);new_scores=normalized_scores(revised,b)
            assert not(old_scores-new_scores),'An existing hinge-score changed under the singleton-only correction.'
            checks.append(dict(cohort=label,existing_bundle_fits_preserved=len(a),new_singleton_bundle_fits=len(b)-len(a),existing_hinge_scores_preserved=sum(old_scores.values())))
        result['all_existing_fitted_bundles_preserved']=checks
    (out/'comparison.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    import matplotlib
    matplotlib.use('Agg');import matplotlib.pyplot as plt
    a={r['name']:r for r in result['additional']['groups']};ref={r['name']:r for r in result['reference']['groups']}
    fig,axes=plt.subplots(2,2,figsize=(14,9),layout='constrained');colors=['#23546e','#b14c30','#408553']
    for j,(method,title) in enumerate((('cascade','Kaskade'),('v1','Ensemble v1'),(revision_label,'Ensemble '+revision_label))):
        axes[0,0].bar(np.arange(2)+(j-1)*.25,[ref['gaussian_sparse'][method]['alarm'],a['gaussian_sparse'][method]['alarm']],width=.24,color=colors[j],label=title)
    axes[0,0].set_xticks([0,1],['Bisherige Trials 0–31',f'Zusätzliche Trials {lo}–{hi-1}']);axes[0,0].set_ylim(0,35);axes[0,0].set_ylabel('Eingaben mit mindestens einem Vorschlag / einer Ausgabe');axes[0,0].set_title('Schmale Peaks: jeweils 32 Eingaben');axes[0,0].legend(fontsize=8)
    names=['gaussian_wide','gaussian_sparse','three_gaussians'];x=np.arange(3)
    for j,(method,label) in enumerate((('cascade_rmse','Kaskade'),('revision_rmse','Ensemble '+revision_label))):
        axes[0,1].bar(x+(j-.5)*.3,[float('nan') if a[n]['cascade_revision_localization'][method] is None else a[n]['cascade_revision_localization'][method] for n in names],width=.29,label=label,color=colors[j*2])
    axes[0,1].set_xticks(x,[f"{label}\nn={a[n]['cascade_revision_localization']['cases']}" for n,label in zip(names,['Breit','Schmal','Drei Komponenten'])]);axes[0,1].set_ylabel('Lage-RMSE, Messpunktindex');axes[0,1].set_title('Zusätzliche Trials: gleiche Fälle, passende Anzahl, keine Junctions');axes[0,1].legend(fontsize=8)
    noise=['flat_log','flat_raw','zero_clipped_raw_24','zero_clipped_raw_2400'];x=np.arange(4)
    for j,(method,label) in enumerate((('cascade','Kaskade'),('v1','Ensemble v1'),(revision_label,'Ensemble '+revision_label))):
        alarms=[a[n][method]['alarm'] for n in noise];errors=[a[n][method]['error'] for n in noise];xx=x+(j-1)*.25
        axes[1,0].bar(xx,alarms,width=.24,color=colors[j],label=label);axes[1,0].bar(xx,errors,bottom=alarms,width=.24,facecolor='none',edgecolor='#8b2222',hatch='///',label='Fehler' if j==0 else None)
    axes[1,0].set_xticks(x,['Log-Rauschen','Roh-Rauschen','Beschnitten 24','Beschnitten 2400'],fontsize=8);axes[1,0].set_ylim(0,35);axes[1,0].set_title('Zusätzliche Rauschkontrollen: jeweils 32 Eingaben');axes[1,0].set_ylabel('Eingaben mit mindestens einer Ausgabe / einem Vorschlag');axes[1,0].legend(fontsize=8)
    overlap=a['overlapping_pair'];upper=max(int(k) for m in ('cascade','v1',revision_label) for k in overlap[m]['counts'] if k!='None');xs=list(range(upper+1))
    for j,(m,label) in enumerate((('cascade','Kaskaden-Ausgaben'),('v1','v1-Vorschläge'),(revision_label,revision_label+'-Vorschläge'))):axes[1,1].plot(xs,[overlap[m]['counts'].get(str(i),0) for i in xs],'o-',label=label,color=colors[j])
    axes[1,1].axvline(2,color='#555555',ls='--',lw=.8);axes[1,1].set_title('Zusätzliche Überlappungen: zwei generierte Komponenten');axes[1,1].set_xlabel('Anzahl Ausgaben bzw. ungeprüfter Vorschläge');axes[1,1].set_ylabel('Eingaben');axes[1,1].legend(fontsize=8)
    fig.suptitle(f"{'Einzelmodell-Basisfall' if revision_label=='v3' else 'Modellgestützte Randknoten'}: {revision_label} gegen eingefrorene v1 und unveränderte Kaskade\nGleiche native Fenster und Intensitätsboden; Ensemble-Geometrie weiterhin ohne statistische Peakfreigabe.",fontsize=12)
    fig.savefig(out/'model-boundary-comparison.png',dpi=145);fig.savefig(out/'model-boundary-comparison.pdf');plt.close(fig)
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('reference_input','reference_v1','reference_v2','new_input','new_v1','new_v2','freeze','output'):p.add_argument(key,type=Path)
    p.add_argument('--numpy-root',type=Path);p.add_argument('--freeze-source',type=Path);a=p.parse_args()
    if a.numpy_root:sys.path.insert(0,str(a.numpy_root.resolve()))
    import numpy as np
    r=run(a.reference_input,a.reference_v1,a.reference_v2,a.new_input,a.new_v1,a.new_v2,a.freeze,a.output,np,a.freeze_source)
    print(json.dumps({k:v for k,v in r.items() if k not in ('reference','additional','fixed_diagnostics')},indent=2))
