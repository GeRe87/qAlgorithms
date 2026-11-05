"""Freeze an untuned paired workload from prior simulations and fresh EICs."""
import argparse,json,hashlib,struct,re
from pathlib import Path
from collections import Counter
from analyse import read_inputs
from real_comparison import read_eics

def prepare(study,eic_binary,out,trials,trial_start=0):
    if trials<1 or trial_start<0:raise ValueError('Trial workload must be positive with nonnegative start.')
    out.mkdir(parents=True,exist_ok=True);source=read_inputs(study/'inputs.bin')
    meta=json.loads((study/'scenarios.json').read_text(encoding='utf-8'));scenarios={s['id']:s for s in meta['scenarios']}
    selected=[];ledger=[]
    for (s,t),(truth,x,y) in source.items():
        m=scenarios[s];randomized=int(m['trials'])>1
        take=trial_start<=t<trial_start+trials if randomized else True
        ledger.append(dict(scenario=s,trial=t,selected=take,reason=('declared_trial_interval' if randomized else 'fixed_control') if take else 'outside_declared_pilot_workload'))
        if take:
            m=scenarios[s];name=m['name'];kind='unknown_real' if m.get('truth_unknown') else 'legacy_expectation' if name.startswith('pt_') or name.startswith('pt_03_') or name.startswith('pt_04_') else 'generated_reference'
            selected.append(dict(scenario=s,trial=t,name=name,truth=list(truth),x=list(x),y=list(y),reference_kind=kind,source='preceding_frozen_study',evaluation_role='trial_interval' if randomized else 'fixed_real_diagnostic' if m.get('truth_unknown') else 'fixed_control'))
    fresh=read_eics(eic_binary);extra=[];nextid=max(scenarios)+100
    for polarity in (1,2):
        eligible=sorted(((k,d) for k,d in fresh.items() if k[0]==polarity and len(d['area'])>8),key=lambda kd:(len(kd[1]['area']),kd[0][1]))
        if not eligible:raise ValueError('Requested real polarity has no eligible EIC.')
        for reason,(key,d) in (('median_length',eligible[len(eligible)//2]),('maximum_length',min(eligible,key=lambda kd:(-len(kd[1]['area']),kd[0][1])))):
            existing=next((c for c in selected if tuple(c['x'])==d['rt'] and tuple(c['y'])==d['area']),None)
            if existing:
                extra.append(dict(polarity=polarity,binID=key[1],reason=reason,duplicate_of=[existing['scenario'],existing['trial']]))
                continue
            selected.append(dict(scenario=nextid,trial=0,name=f'real_polarity_{polarity}_{reason}',truth=[],x=list(d['rt']),y=list(d['area']),reference_kind='unknown_real',source='fresh_mzml_EIC',polarity=polarity,binID=key[1],evaluation_role='fixed_real_diagnostic'))
            extra.append(dict(polarity=polarity,binID=key[1],reason=reason,case=[nextid,0]));nextid+=1
    with (out/'inputs.bin').open('wb') as f:
        f.write(b'QENS0001')
        for c in selected:
            n=len(c['y']);f.write(struct.pack('<4I',c['scenario'],c['trial'],n,len(c['truth'])))
            for values in (c['truth'],c['x'],c['y']):f.write(struct.pack('<'+'f'*len(values),*values))
            f.write(struct.pack('<'+'H'*n,*range(1,n+1)))
    descriptions=[{k:v for k,v in c.items() if k not in ('x','y')} for c in selected]
    (out/'cases.json').write_text(json.dumps(descriptions,indent=2),encoding='utf-8')
    (out/'source-selection-ledger.json').write_text(json.dumps(ledger,indent=2),encoding='utf-8')
    protocol=dict(selected_cases=len(selected),trials_per_existing_scenario_cap=trials,source_records_accounted=len(source),source_outside_workload=sum(not r['selected'] for r in ledger),extra_real_selection=extra,reference_counts=dict(Counter(c['reference_kind'] for c in selected)),input_sha256=hashlib.sha256((out/'inputs.bin').read_bytes()).hexdigest(),source_sha256=hashlib.sha256((study/'inputs.bin').read_bytes()).hexdigest(),fresh_eics_sha256=hashlib.sha256(eic_binary.read_bytes()).hexdigest(),algorithm_parameters_tuned=False,workload_rule='First declared trials of each fixed scenario; all shorter scenarios; per real polarity median and maximum complete-EIC length among n>8, lowest binID breaks maximum ties. Exact x/y duplicates retained once.',statistical_ensemble_acceptance_defined=False)
    protocol['trial_start']=trial_start
    protocol['randomized_records_selected']=sum(c['evaluation_role']=='trial_interval' for c in selected)
    protocol['fixed_records_selected']=len(selected)-protocol['randomized_records_selected']
    protocol['workload_rule']='Declared trial interval in each repeated scenario; all fixed controls reused and identified; same deterministic real EIC selection. Workload is not an algorithm threshold.'
    (out/'protocol.json').write_text(json.dumps(protocol,indent=2),encoding='utf-8');return protocol

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('study','eic_binary','output'):p.add_argument(key,type=Path)
    p.add_argument('--trials',type=int,default=32);p.add_argument('--trial-start',type=int,default=0);a=p.parse_args();print(json.dumps(prepare(a.study,a.eic_binary,a.output,a.trials,a.trial_start),indent=2))
