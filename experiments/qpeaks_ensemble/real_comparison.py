"""Apply the fixed-support comparison to freshly exported real mzML EICs.
Prior production labels specify sampling, not truth or new acceptance rules.
"""
import argparse,csv,json,sys,struct,math,hashlib,time
from pathlib import Path
from collections import Counter
from compare_components import Geometry,compare

def rows(path):
    with path.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))

def read_eics(path):
    b=path.read_bytes()
    if b[:8]!=b'QEIC0001':raise ValueError('Unknown EIC format.')
    offset=8;result={}
    while offset<len(b):
        if offset+12>len(b):raise ValueError('Truncated EIC header.')
        polarity,ident,n=struct.unpack_from('<3I',b,offset);offset+=12
        if polarity not in (1,2) or n<1 or n>65535:raise ValueError('Invalid EIC header.')
        if offset+18*n>len(b):raise ValueError('Truncated EIC vectors.')
        data={}
        for name in ('rt','area','mz'):
            data[name]=struct.unpack_from('<'+'f'*n,b,offset);offset+=4*n
        data['df']=struct.unpack_from('<'+'H'*n,b,offset);offset+=2*n
        data['scan']=struct.unpack_from('<'+'I'*n,b,offset);offset+=4*n
        key=(polarity,ident)
        if key in result:raise ValueError('Duplicate EIC identity.')
        if data['df']!=tuple(range(1,n+1)):raise ValueError('Unsupported measured-point counts.')
        if not all(math.isfinite(v) and v>0 for v in data['area']):raise ValueError('Invalid EIC intensities.')
        if not all(math.isfinite(v) for v in (*data['rt'],*data['mz'])):raise ValueError('Nonfinite EIC coordinates.')
        result[key]=data
    if offset!=len(b):raise ValueError('EIC byte accounting mismatch.')
    return result

def select_windows(eics,feature_dir,historical_bins=()):
    annotations=[];provenance=[]
    for polarity,label in ((1,'positive'),(2,'negative')):
        matches=list(feature_dir.glob('*_'+label+'_area_timing_audit.csv'))
        if len(matches)!=1:raise ValueError('Expected exactly one timing audit per polarity.')
        audit=matches[0];feature=audit.with_name(audit.name.replace('_area_timing_audit.csv','_features.csv'))
        feature_rows=rows(feature);actual={int(r['ID']):r for r in feature_rows};timing=rows(audit)
        if len(actual)!=len(feature_rows) or len(actual)!=len(timing) or {int(r['featureID']) for r in timing}!=set(actual):
            raise ValueError('Feature/audit identities are not one-to-one.')
        provenance.extend([dict(file=p.name,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in (audit,feature)])
        for r in timing:
            ident=int(r['featureID']);old=actual[ident];binid=int(r['binID']);start=int(r['idxPeakStart']);n=int(r['length'])
            data=eics[(polarity,binid)]
            if start<0 or n<1 or start+n>len(data['rt']):raise ValueError('Feature window outside fresh EIC.')
            if old['binID']!=r['binID'] or old['isUnresolved']!=r['isUnresolved']:raise ValueError('Feature metadata mismatch.')
            # The existing audit serialized float RTs with nine significant
            # digits. Compare that serialization, not an invented RT tolerance.
            if float(format(data['rt'][start],'.9g'))!=float(r['lowerRT']) or float(format(data['rt'][start+n-1],'.9g'))!=float(r['upperRT']):
                raise ValueError('Fresh EIC does not match the saved feature window.')
            annotations.append(dict(polarity=polarity,featureID=ident,binID=binid,start=start,length=n,old_collective=int(r['isUnresolved']),reasons=[]))
    selected={}
    def choose(r,reason):
        key=(r['polarity'],r['binID'],r['start'],r['length'])
        if key not in selected:selected[key]=dict(r,reasons=[])
        selected[key]['reasons'].append(reason)
    for polarity in (1,2):
        for collective in (0,1):
            eligible=sorted((r for r in annotations if r['polarity']==polarity and r['old_collective']==collective and r['length']>8),key=lambda r:(r['length'],r['binID'],r['featureID']))
            if not eligible:raise ValueError('Requested real sampling stratum is empty.')
            choose(eligible[len(eligible)//2],f'length_median_old_class_{collective}')
            choose(eligible[-1],f'length_maximum_old_class_{collective}')
    for binid in historical_bins:
        eligible=sorted((r for r in annotations if r['polarity']==1 and r['binID']==binid and r['old_collective']==1 and r['length']>8),key=lambda r:r['featureID'])
        if not eligible:raise ValueError('A historical conflict EIC has no eligible prior collective window.')
        choose(eligible[0],'historical_conflict_first_collective_feature')
    return annotations,list(selected.values()),provenance

def diagnostics(residual,prediction,np):
    residual=np.asarray(residual);prediction=np.asarray(prediction)
    def correlation(a,b):
        a=a-a.mean();b=b-b.mean();den=float(np.linalg.norm(a)*np.linalg.norm(b))
        return float(a@b/den) if den else None
    return dict(rss=float(residual@residual),mean_log_residual=float(residual.mean()),lag1_adjacent_observation_correlation=correlation(residual[:-1],residual[1:]),squared_residual_vs_fitted_log_correlation=correlation(residual*residual,prediction))

def run(export,features,mzml,out,np,scipy,historical_bins=()):
    from scipy.optimize import lsq_linear,least_squares
    out.mkdir(parents=True,exist_ok=True);eics=read_eics(export/'eics.bin')
    inventory=rows(export/'inventory.csv');native=rows(export/'polarities.csv')
    if len(eics)!=len(inventory) or {(int(r['polarity']),int(r['binID'])) for r in inventory}!=set(eics):raise ValueError('Incomplete EIC inventory.')
    if len(native)!=2 or {int(r['polarity']) for r in native}!={1,2}:raise ValueError('Incomplete polarity inventory.')
    for r in native:
        selected=[d for key,d in eics.items() if key[0]==int(r['polarity'])]
        if len(selected)!=int(r['eics']) or sum(len(d['area']) for d in selected)!=int(r['eic_observations']):raise ValueError('Polarity accounting mismatch.')
    for r in inventory:
        d=eics[(int(r['polarity']),int(r['binID']))]
        if len(d['area'])!=int(r['length']):raise ValueError('Inventory length mismatch.')
    annotations,selected,provenance=select_windows(eics,features,historical_bins)
    keys={(r['polarity'],r['binID'],r['start'],r['length']):i for i,r in enumerate(selected)}
    with (out/'feature-selection-ledger.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.writer(f);w.writerow(['polarity','featureID','binID','start','length','old_collective','selected_case','reason'])
        for r in annotations:
            key=r['polarity'],r['binID'],r['start'],r['length'];case=keys.get(key)
            w.writerow([r[k] for k in ('polarity','featureID','binID','start','length','old_collective')]+['' if case is None else case,'n_not_greater_than_8' if r['length']<=8 else 'selected_prior_window' if case is not None else 'outside_declared_real_sampling_protocol'])
    cases=[]
    for i,r in enumerate(selected):
        d=eics[(r['polarity'],r['binID'])];s=r['start'];end=s+r['length'];rt=d['rt'][s:end]
        if any(rt[j+1]<=rt[j] for j in range(len(rt)-1)):raise ValueError('Selected real RTs must increase.')
        cases.append(dict(r,case=i,rt=list(rt),y=list(d['area'][s:end]),mz=list(d['mz'][s:end]),scan=list(d['scan'][s:end])))
    (out/'comparison-inputs.json').write_text(json.dumps(cases,indent=2),encoding='utf-8')
    protocol=dict(mzml_name=mzml.name,mzml_bytes=mzml.stat().st_size,mzml_sha256=hashlib.sha256(mzml.read_bytes()).hexdigest(),eic_binary_sha256=hashlib.sha256((export/'eics.bin').read_bytes()).hexdigest(),fresh_eic_records=len(eics),native_polarity_counts=native,all_prior_feature_windows_verified=len(annotations),feature_provenance=provenance,selected_windows=len(cases),sampling='Per polarity/old class: median and maximum length among n>8; plus first prior collective window in six historical positive conflict EICs. Deduplicated, chosen before new fits.',prior_labels_are_truth=False,mean_model='Same fixed-midpoint concave log-profile family as the preceding one-versus-two reference; raw intensities add; iid log errors are a working assumption, not an instrument calibration.',index_model=True,rt_resampled=False,area_corrected=False,peak_count_release=False,scipy=scipy.__version__,numpy=np.__version__)
    protocol['historical_positive_bins']=list(historical_bins)
    protocol['sampling']='Per polarity/old class: median and maximum length among n>8; plus first prior collective window in each explicitly supplied historical positive conflict EIC. Deduplicated, chosen before new fits.'
    (out/'protocol.json').write_text(json.dumps(protocol,indent=2),encoding='utf-8')
    summaries=[];curves={};cache={}
    sf=(out/'starts.csv').open('w',newline='',encoding='utf-8');sw=csv.writer(sf)
    sw.writerow(['case','start','rss','success','status','nfev','error',*['theta_'+str(i) for i in range(8)]])
    cf=(out/'curves.csv').open('w',newline='',encoding='utf-8');cw=csv.writer(cf)
    cw.writerow(['case','index','rt','observed','one_log_median','two_log_median','one_log_residual','two_log_residual'])
    pf=(out/'components.csv').open('w',newline='',encoding='utf-8');pw=csv.writer(pf)
    pw.writerow(['case','component','log_intercept','slope_per_sample','left_curvature','right_curvature','apex_index','interior_apex'])
    for c in cases:
        began=time.perf_counter();n=len(c['y'])
        if n not in cache:cache[n]=Geometry(n,np)
        g=cache[n];base={k:c[k] for k in ('case','polarity','featureID','binID','start','length','old_collective')}
        try:
            r=compare(c['y'],g,lsq_linear,least_squares);one=r['one'];best=r['best']
            d1=diagnostics(one['curve']-r['z'],one['curve'],np);d2=diagnostics(r['curve_two']-r['z'],r['curve_two'],np)
            dt=np.diff(c['rt']);cv=float(np.std(dt)/np.mean(dt))
            row=dict(base,status='best_found_converged' if best['success'] else 'best_found_unconverged',one_rss=one['rss'],two_rss=best['rss'],best_start=best['start'],starts=len(r['starts']),unsuccessful_starts=sum(not s['success'] for s in r['starts']),one_lag1=d1['lag1_adjacent_observation_correlation'],two_lag1=d2['lag1_adjacent_observation_correlation'],one_squared_residual_correlation=d1['squared_residual_vs_fitted_log_correlation'],two_squared_residual_correlation=d2['squared_residual_vs_fitted_log_correlation'],jacobian_rank=r['rank'],jacobian_condition=r['condition'],rt_interval_cv=cv,seconds=time.perf_counter()-began,error='')
            for s in r['starts']:sw.writerow([c['case'],s['start'],s['rss'],s['success'],s['status'],s['nfev'],s['error'],*(['']*8 if s['theta'] is None else list(map(float,s['theta'])))])
            for j,p in enumerate(r['components']):pw.writerow([c['case'],j,*[p[k] for k in ('log_intercept','slope_per_sample','left_curvature','right_curvature','apex','interior_apex')]])
            onecurve=one['curve']+r['anchor'];twocurve=r['curve_two']+r['anchor'];curves[c['case']]=(onecurve,twocurve)
            for i in range(n):cw.writerow([c['case'],i,c['rt'][i],c['y'][i],onecurve[i],twocurve[i],one['curve'][i]-r['z'][i],r['curve_two'][i]-r['z'][i]])
        except Exception as e:
            row=dict(base,status='case_error',one_rss=None,two_rss=None,best_start=None,starts=0,unsuccessful_starts=0,one_lag1=None,two_lag1=None,one_squared_residual_correlation=None,two_squared_residual_correlation=None,jacobian_rank=None,jacobian_condition=None,rt_interval_cv=None,seconds=time.perf_counter()-began,error=str(e))
        summaries.append(row)
        print(f"case={c['case']} polarity={c['polarity']} EIC={c['binID']} n={n} status={row['status']} seconds={row['seconds']:.3f}",flush=True)
    for f in (sf,cf,pf):f.close()
    with (out/'summary.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(summaries[0]));w.writeheader();w.writerows(summaries)
    report=dict(fresh_eics=len(eics),prior_feature_windows_bound=len(annotations),compared_windows=len(cases),status_counts=dict(Counter(r['status'] for r in summaries)),unsuccessful_starts=sum(r['unsuccessful_starts'] for r in summaries),physical_truth_known=False,automatic_peak_count_decision=False,residual_diagnostics='Descriptive adjacent-observation correlations; systematic mean-model error can cause them. They do not identify instrument noise correlation or a validated false-alarm rate.')
    (out/'report.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    import matplotlib
    matplotlib.use('Agg');import matplotlib.pyplot as plt
    # One resolved and one collective median-length reference in each polarity.
    examples=[next(c for c in cases if c['polarity']==p and f'length_median_old_class_{k}' in c['reasons']) for p in (1,2) for k in (0,1)]
    fig,axes=plt.subplots(4,2,figsize=(13,11),layout='constrained')
    for ax,c in zip(axes,examples):
        if c['case'] not in curves:
            for a in ax:a.text(.5,.5,f"case {c['case']}: fit error; see summary.csv",ha='center',transform=a.transAxes)
            continue
        one,two=curves[c['case']];rt=c['rt'];y=np.asarray(c['y']);z=np.log(y)
        ax[0].plot(rt,y,label='Reale Messwerte',color='#23546e');ax[0].plot(rt,np.exp(one),label='Ein Profil',color='#b14c30');ax[0].plot(rt,np.exp(two),label='Zwei Profile',color='#408553')
        ax[0].set_title(f"{'Positiv' if c['polarity']==1 else 'Negativ'}, EIC {c['binID']}, altes Kollektivlabel {c['old_collective']}");ax[0].set_ylabel('Intensität');ax[0].set_xlabel('Originale Retentionszeit');ax[0].legend(fontsize=7)
        ax[1].plot(rt,one-z,label='Ein Profil',color='#b14c30');ax[1].plot(rt,two-z,label='Zwei Profile',color='#408553');ax[1].axhline(0,color='#777777',lw=.6);ax[1].set_xlabel('Originale Retentionszeit');ax[1].set_ylabel('Log-Residuum');ax[1].legend(fontsize=7)
    fig.suptitle('Reale mzML: Anpassung und verbleibende Residuen auf vorab ausgewählten Fenstern\nFit weiterhin auf Messpunktindex; Original-RT nur zur Zuordnung. Alte Labels sind keine Referenzwahrheit.',fontsize=11)
    fig.savefig(out/'real-comparison.png',dpi=145);fig.savefig(out/'real-comparison.pdf');plt.close(fig)
    return report

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('export',type=Path);p.add_argument('features',type=Path);p.add_argument('mzml',type=Path);p.add_argument('output',type=Path);p.add_argument('--numpy-root',type=Path);p.add_argument('--scipy-root',type=Path);p.add_argument('--historical-positive-bin',type=int,action='append',default=[]);a=p.parse_args()
    for path in (a.numpy_root,a.scipy_root):
        if path:sys.path.insert(0,str(path.resolve()))
    import numpy as np,scipy
    print(json.dumps(run(a.export,a.features,a.mzml,a.output,np,scipy,a.historical_positive_bin),indent=2))
if __name__=='__main__':main()
