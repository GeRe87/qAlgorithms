"""Reproducible topology audit and correlated-ensemble reference calculation.
No ensemble peak detector or significance cutoff is defined here.
"""
import argparse,csv,json,math,sys,hashlib
from pathlib import Path
from collections import defaultdict,Counter

def design(scale,np):
    x=np.arange(-scale,scale+1,dtype=float)
    return np.column_stack([np.ones_like(x),x,np.where(x<0,x*x,0),np.where(x>0,x*x,0)])

def operator(n,center,scale,np):
    X=design(scale,np);Q,R=np.linalg.qr(X,mode='reduced')
    A=np.zeros((4,n));A[:,center-scale:center+scale+1]=np.linalg.solve(R,Q.T)
    return A

def correlated_reference(n,np):
    # Predetermined ensemble of ALL local windows, without geometric selection.
    t=(np.arange(n)-(n-1)/2)/((n-1)/2)
    X=np.column_stack([np.ones(n),t,t*t]);blocks=[];information=np.zeros((3,3));rhs=np.zeros((3,n))
    for scale in range(2,(n-1)//2+1):
        for center in range(scale,n-scale):
            A=operator(n,center,scale,np);B=A@X
            W=np.linalg.inv(A@A.T)
            information+=B.T@W@B;rhs+=B.T@W@A;blocks.append(A)
    A=np.vstack(blocks);_,singular,Vh=np.linalg.svd(A,full_matrices=False)
    rank=int(sum(singular>max(A.shape)*np.finfo(float).eps*singular[0]))
    P=Vh[:rank].T@Vh[:rank]
    cov_gls=np.linalg.inv(X.T@P@X);K_gls=cov_gls@X.T@P
    cov_direct=np.linalg.inv(X.T@X);K_direct=cov_direct@X.T
    cov_naive_reported=np.linalg.inv(information);K_naive=cov_naive_reported@rhs
    cov_naive_actual=K_naive@K_naive.T
    difference=float(np.max(np.abs(K_gls-K_direct)))
    assert difference<1e-10
    assert rank==n
    assert np.linalg.eigvalsh(cov_naive_actual-cov_direct).min()>-1e-10
    ratio=np.sqrt(np.diag(cov_naive_actual)/np.diag(cov_naive_reported))
    assert np.all(ratio>1)
    return dict(samples=n,local_models=len(blocks),stacked_coefficients=int(A.shape[0]),operator_rank=rank,stacked_covariance_rank=rank,stacked_covariance_nullity=int(A.shape[0]-rank),GLS_vs_direct_OLS_operator_max_abs=difference,naive_actual_to_reported_standard_error_ratio=ratio.tolist(),naive_actual_variance_over_direct_OLS=(np.diag(cov_naive_actual)/np.diag(cov_direct)).tolist(),assumptions='common global quadratic log-mean; iid unit-variance errors; predetermined windows, no data selection; linear reference, not a new peak estimator')

def read_inputs(path):
    import struct
    b=path.read_bytes();assert b[:8]==b'QENS0001';offset=8;records={}
    while offset<len(b):
        scenario,trial,n,t=struct.unpack_from('<4I',b,offset);offset+=16
        truth=struct.unpack_from('<'+'f'*t,b,offset);offset+=4*t
        x=struct.unpack_from('<'+'f'*n,b,offset);offset+=4*n
        y=struct.unpack_from('<'+'f'*n,b,offset);offset+=4*n
        df=struct.unpack_from('<'+'H'*n,b,offset);offset+=2*n
        assert all(math.isfinite(v) and v>=0 for v in y)
        assert df==tuple(range(1,n+1))
        key=(scenario,trial);assert key not in records
        records[key]=(truth,x,y)
    assert offset==len(b)
    return records

def csvrows(path):
    with path.open(newline='',encoding='utf-8') as f:yield from csv.DictReader(f)

def lineage(arr,np):
    # Contract only existing horizontal grid edges. Keep every vertical edge
    # between runs, including multiple parents/children; no peak assignment.
    rows=sorted((int(i) for i in np.flatnonzero(arr[:,2]==1)),key=lambda i:(arr[i,1],arr[i,0]))
    groups=[]
    for i in rows:
        if not groups or arr[i,1]!=arr[groups[-1][-1],1] or arr[i,0]!=arr[groups[-1][-1],0]+1:groups.append([])
        groups[-1].append(i)
    node_run={i:r for r,g in enumerate(groups) for i in g}
    grid={(int(arr[i,0]),int(arr[i,1])):i for i in rows};edges=Counter()
    for i in rows:
        j=grid.get((int(arr[i,0]),int(arr[i,1])+1))
        if j is not None:edges[(node_run[i],node_run[j])]+=1
    parents=defaultdict(set);children=defaultdict(set)
    for a,b in edges:children[a].add(b);parents[b].add(a)
    result=[]
    for r,g in enumerate(groups):
        a=arr[g];component=int(a[0,6]);assert np.all(a[:,6]==component)
        result.append(dict(run=r,scale=int(a[0,1]),first_center=int(a[0,0]),last_center=int(a[-1,0]),members=len(g),min_apex=float(a[:,3].min()),max_apex=float(a[:,3].max()),display_mean_apex=float(a[:,3].mean()),current_accepted=int(sum(a[:,5])),constant_members=int(sum(a[:,4])),component=component,parents=len(parents[r]),children=len(children[r])))
    assert sum(r['members'] for r in result)==len(rows)
    return result,edges

def run(folder,output,np,plots,collection_name='release'):
    output.mkdir(parents=True,exist_ok=True)
    meta=json.loads((folder/'scenarios.json').read_text(encoding='utf-8'));records=read_inputs(folder/'inputs.bin')
    assert len(records)==meta['records']
    assert hashlib.sha256((folder/'inputs.bin').read_bytes()).hexdigest()==meta['input_sha256']
    traces={(int(r['scenario']),int(r['trial'])):r for r in csvrows(folder/collection_name/'traces.csv')}
    assert set(traces)==set(records)
    byscenario={s['id']:s for s in meta['scenarios']}
    arrays={k:np.empty((int(r['all_regressions']),12)) for k,r in traces.items() if r['detail']=='1'}
    positions=Counter()
    for r in csvrows(folder/collection_name/'candidates.csv'):
        key=(int(r['scenario']),int(r['trial']));i=int(r['id']);assert i==positions[key];positions[key]+=1
        arrays[key][i]=[int(r['center']),int(r['scale']),int(r['has_interior_maximum']),float(r['apex']) if r['apex'] else np.nan,int(r['constant_log_window']),int(r['current_accepted']),int(r['component']) if r['component'] else -1,float(r['full_rss']),float(r['b1']),float(r['b2']),float(r['b3']),int(r['all_above_production_floor'])]
    assert set(positions)<=set(arrays)
    group_sizes=defaultdict(Counter);group_ranges=defaultdict(dict);max_scales=defaultdict(int)
    for r in csvrows(folder/collection_name/'components.csv'):
        key=(int(r['scenario']),int(r['trial']));component=int(r['component'])
        assert component not in group_ranges[key]
        group_sizes[key][component]=int(r['nodes'])
        group_ranges[key][component]=(float(r['min_apex']),float(r['max_apex']))
        max_scales[key]=max(max_scales[key],int(r['scales']))
    edge_counts=Counter()
    for r in csvrows(folder/collection_name/'edges.csv'):
        key=(int(r['scenario']),int(r['trial']));a,b=int(r['from_id']),int(r['to_id']);arr=arrays[key]
        assert arr[a,2]==arr[b,2]==1 and arr[a,6]==arr[b,6]
        dc,ds=arr[b,:2]-arr[a,:2]
        assert (dc==1 and ds==0 and r['axis']=='center') or (dc==0 and ds==1 and r['axis']=='scale')
        edge_counts[key]+=1
    for key,r in traces.items():
        truth,x,y=records[key];assert int(r['n'])==len(y) and int(r['truth_count'])==len(truth)
        runs=[];length=0
        for v in y:
            if v>0:length+=1
            else:
                if length:runs.append(length)
                length=0
        if length:runs.append(length)
        expected=sum(sum(n-2*s for s in range(2,min(215,(n-1)//2)+1)) for n in runs if n>=5)
        assert expected==int(r['all_regressions'])
        assert sum(n for n in runs if n<5)==int(r['short_positive_samples'])
        assert sum(v<=0 for v in y)==int(r['nonpositive_samples'])
        assert len(group_sizes[key])==int(r['graph_components'])
        assert sum(group_sizes[key].values())==int(r['interior_maxima'])
        if r['strict_error']:assert r['strict_output_count']=='-1'
        else:assert int(r['strict_output_count'])>=0
        if key in arrays:
            arr=arrays[key];assert positions[key]==len(arr)
            assert int(sum(arr[:,2]))==int(r['interior_maxima'])
            assert int(sum(arr[:,4]))==int(r['constant_log_windows'])
            assert int(sum(arr[:,5]))==int(r['current_accepted_windows'])
            assert np.isfinite(arr[:,[0,1,2,4,5,6,7,8,9,10,11]]).all()
            assert np.isfinite(arr[arr[:,2]==1,3]).all()
            assert Counter(arr[arr[:,2]==1,6].astype(int))==group_sizes[key]
            assert edge_counts[key]==int(r['graph_edges'])
    lineages={};lineage_summary=[]
    with (output/'runs.csv').open('w',newline='',encoding='utf-8') as f, (output/'lineage-edges.csv').open('w',newline='',encoding='utf-8') as ef:
        writer=None;ew=csv.writer(ef);ew.writerow(['scenario','trial','from_run','to_run','shared_centers'])
        for key,arr in arrays.items():
            runs,le=lineage(arr,np);lineages[key]=(runs,le)
            for r in runs:
                row=dict(scenario=key[0],trial=key[1],**r)
                if writer is None:writer=csv.DictWriter(f,fieldnames=list(row));writer.writeheader()
                writer.writerow(row)
            for (a,b),weight in sorted(le.items()):
                assert runs[b]['scale']==runs[a]['scale']+1
                ew.writerow([*key,a,b,weight])
            lineage_summary.append(dict(scenario=key[0],trial=key[1],runs=len(runs),vertical_run_edges=len(le),births=sum(r['parents']==0 for r in runs),merge_nodes=sum(r['parents']>1 for r in runs),split_nodes=sum(r['children']>1 for r in runs),deaths=sum(r['children']==0 for r in runs)))
    reflection=[]
    for s in meta['scenarios']:
        if 'reflection_of' not in s:continue
        a=arrays[(s['reflection_of'],0)];b=arrays[(s['id'],0)];n=s['samples']
        lookup={(int(row[0]),int(row[1])):row for row in b};mismatches=0;max_error=0
        for row in a:
            other=lookup[(n-1-int(row[0]),int(row[1]))]
            if row[2]!=other[2]:mismatches+=1
            elif row[2]:max_error=max(max_error,abs(row[3]+other[3]-(n-1)))
        reflection.append(dict(scenario=s['name'],candidate_count=len(a),maximum_classification_mismatches=mismatches,max_reflected_apex_sum_error=max_error))
    aggregated=[]
    for s in meta['scenarios']:
        rs=[r for (scenario,trial),r in traces.items() if scenario==s['id']]
        assert len(rs)==s['trials']
        aggregated.append(dict(scenario=s['id'],name=s['name'],trials=len(rs),traces_with_any_component=sum(int(r['graph_components'])>0 for r in rs),median_graph_components=float(np.median([int(r['graph_components']) for r in rs])),median_largest_scale_count=float(np.median([max_scales[(s['id'],int(r['trial']))] for r in rs])),strict_errors=sum(bool(r['strict_error']) for r in rs),strict_positive=sum(not r['strict_error'] and int(r['strict_output_count'])>0 for r in rs),strict_zero=sum(not r['strict_error'] and r['strict_output_count']=='0' for r in rs)))
    references=[correlated_reference(n,np) for n in (9,21,41)]
    report=dict(records=len(records),all_records_and_detailed_graphs_verified=True,production_policy_changed=False,new_peak_decisions_made=False,scenarios=aggregated,reflection=reflection,lineage=lineage_summary,correlated_ensemble_reference=references,constant_control=traces[(3000,0)],coordinate_system='Graph/apex coordinates use sample indices. Original x is retained for current-cascade comparison and observations; physical RT was not regularized.')
    (output/'analysis.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    if plots:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        for ids,name in [([8,9,3000,7],'controls'),([4,5,6,3001],'signals'),([2000,2001,2002,2003],'real-eics')]:
            ids=[i for i in ids if (i,0) in records]
            if not ids:continue
            fig,axes=plt.subplots(len(ids),2,figsize=(13,3.1*len(ids)),layout='constrained',squeeze=False)
            for row,scenario in enumerate(ids):
                key=(scenario,0);truth,x,y=records[key];arr=arrays[key];mask=arr[:,2]==1;accepted=arr[:,5]==1
                axes[row,0].plot(range(len(y)),y,color='#23546e',lw=1)
                for position in truth:axes[row,0].axvline(position,color='#9c3e31',lw=.8,ls='--')
                axes[row,0].set_title(byscenario[scenario]['name']);axes[row,0].set_ylabel('Intensity');axes[row,0].set_xlabel('Sample index')
                # Points are fits, not independent observations. Marker opacity
                # and size are display settings, never acceptance thresholds.
                axes[row,1].scatter(arr[mask,3],arr[mask,1],c=arr[mask,6]%20,cmap='tab20',s=18 if sum(mask)<200 else 2,alpha=.8 if sum(mask)<200 else .35,rasterized=True)
                axes[row,1].scatter(arr[accepted,3],arr[accepted,1],s=9,facecolors='none',edgecolors='black',linewidths=.4,rasterized=True)
                axes[row,1].set_xlabel('Fitted apex (sample index)');axes[row,1].set_ylabel('Scale (half-window)')
                axes[row,1].set_title(f"{int(sum(mask))} max-bearing fits; {int(traces[key]['graph_components'])} components")
                axes[row,1].set_xlim(-.5,len(y)-.5)
            fig.suptitle('Diagnostic graph only: colours are adjacency components; black rings pass current local cascade.\nConnected components are NOT inferred peaks. No significance filter applied to coloured fits.',fontsize=11)
            fig.savefig(output/f'{name}.png',dpi=150);fig.savefig(output/f'{name}.pdf');plt.close(fig)
        fig,axes=plt.subplots(2,2,figsize=(13,7),layout='constrained')
        for row,scenario in enumerate((4,7)):
            truth,x,y=records[(scenario,0)];runs,le=lineages[(scenario,0)]
            axes[row,0].plot(range(len(y)),y,color='#23546e')
            for position in truth:axes[row,0].axvline(position,color='#a33b30',ls='--',lw=.8)
            axes[row,0].set_title('Einzelpeak' if scenario==4 else 'Zwei überlappende Peaks')
            axes[row,0].set_xlabel('Messpunktindex');axes[row,0].set_ylabel('Intensität')
            for a,b in le:
                axes[row,1].plot([runs[a]['display_mean_apex'],runs[b]['display_mean_apex']],[runs[a]['scale'],runs[b]['scale']],color='#b9c2c7',lw=.7)
            for r in runs:
                axes[row,1].plot([r['min_apex'],r['max_apex']],[r['scale'],r['scale']],color='#23546e',lw=1.2)
                axes[row,1].plot(r['display_mean_apex'],r['scale'],'o',color='#ba4e2e' if r['parents']>1 else '#23546e',ms=4 if r['parents']>1 else 2)
            axes[row,1].set_xlim(-.5,len(y)-.5);axes[row,1].set_xlabel('Gefitteter Apex (Messpunktindex)');axes[row,1].set_ylabel('Skala (halbe Fensterbreite)')
        fig.suptitle('Experimentelle Bündelstruktur: Verzweigungen bleiben erhalten\nStrecken zeigen die Apexspanne der Fits, keine Konfidenzintervalle. Orange: Zusammenführung im Nachbarschaftsgraphen.',fontsize=11)
        fig.savefig(output/'lineage-example.png',dpi=160);fig.savefig(output/'lineage-example.pdf');plt.close(fig)
    return report

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('folder',type=Path);parser.add_argument('output',type=Path);parser.add_argument('--numpy-root',type=Path);parser.add_argument('--collection-name',default='release');parser.add_argument('--no-plots',action='store_true');args=parser.parse_args()
    if args.numpy_root:sys.path.insert(0,str(args.numpy_root.resolve()))
    import numpy as np
    report=run(args.folder,args.output,np,not args.no_plots,args.collection_name)
    print(json.dumps(dict(records=report['records'],verified=True,references=report['correlated_ensemble_reference'],selected_scenarios=[s for s in report['scenarios'] if s['scenario'] in (0,1,2,3,4,5,6,7,8,9,3000,3001)]),indent=2))
if __name__=='__main__':main()
