"""Check window-root anchors, v1 graph retention and reflected model geometry."""
import argparse,json,math,sys,csv
from pathlib import Path
from collections import defaultdict
from compare_bundle_cascade import fields
from bundle_method import topology

def canonical(bs,edges,n=None):
    def key(r):
        positions=tuple(r['positions']);kind=r['kind']
        if n is not None:
            positions=tuple(sorted(n-1-c for c in positions));kind={'left_boundary_exit':'right_boundary_exit','right_boundary_exit':'left_boundary_exit'}.get(kind,kind)
        return r['scale'],kind,positions
    keys={r['id']:key(r) for r in bs}
    return set(keys.values()),{(keys[a],keys[b],weight) for (a,b),weight in edges.items()},keys

def verify(candidates,z,scope='all_edges'):
    old,old_edges,old_groups,_=topology(candidates);bs,edges,groups,membership=topology(candidates,model_boundaries=True,log_values=z,boundary_scope=scope)
    assert len(membership)==len(candidates)==sum(len(r['candidate_ids']) for r in bs)
    signature=lambda r:(r['block'],r['scale'],tuple(r['positions']))
    new_modes={signature(r):r['id'] for r in bs if r['kind']=='maximum'}
    assert set(new_modes)=={signature(r) for r in old if r['kind']=='maximum'}
    mapping={r['id']:new_modes[signature(r)] for r in old if r['kind']=='maximum'};inverse={v:k for k,v in mapping.items()}
    induced={(inverse[a],inverse[b]):w for (a,b),w in edges.items() if a in inverse and b in inverse}
    assert induced==old_edges
    index={r['id']:r for r in candidates};bounds={}
    block_counts=defaultdict(int)
    for r in candidates:
        block_counts[r['block']]+=1
        key=r['block'],r['scale'];c=r['center']
        if key not in bounds:bounds[key]=[c,c]
        else:bounds[key]=[min(bounds[key][0],c),max(bounds[key][1],c)]
    anchors=0
    for b in bs:
        if b['kind']!='boundary_model_maximum':continue
        if scope=='singleton_blocks':assert block_counts[b['block']]==1
        for a in b['model_anchors']:
            r=index[a['candidate_id']];c=r['center'];s=r['scale'];local=z[c-s:c+s+1]
            assert len(local)==2*s+1 and all(math.isfinite(v) for v in local) and any(v!=local[0] for v in local)
            assert r['b1']-2*s*r['b2']>0 and r['b1']+2*s*r['b3']<0
            root=(-r['b1']/(2*r['b2']) if r['b1']<0 else -r['b1']/(2*r['b3']) if r['b1']>0 else 0.)
            assert -s<root<s
            # Compare the independent local-root formula at floating precision;
            # never turn a fitted root span into a statistical confidence interval.
            assert abs(c+root-a['apex'])<=8*sys.float_info.epsilon*max(1.,abs(c+root))
            lo,hi=bounds[(r['block'],s)]
            if a['side']=='left':assert c==lo and a['apex']<=lo
            else:assert c==hi and a['apex']>=hi
            anchors+=1
    blocks=sorted({r['block'] for r in candidates});reverse_blocks={b:blocks[-1-i] for i,b in enumerate(blocks)}
    reflected=[dict(r,center=len(z)-1-r['center'],block=reverse_blocks[r['block']],b1=-r['b1'],b2=r['b3'],b3=r['b2']) for r in candidates]
    rb,re,_,_=topology(reflected,model_boundaries=True,log_values=z[::-1],boundary_scope=scope)
    assert canonical(bs,edges)[:2]==canonical(rb,re,len(z))[:2]
    if scope=='singleton_blocks':
        def unchanged(nodes,gs):
            return {tuple(sorted(i for b in g['basins'] for i in nodes[b]['candidate_ids'])):(g['contains_merge'],g['contains_split']) for g in gs if block_counts[nodes[g['basins'][0]]['block']]>1}
        assert unchanged(old,old_groups)==unchanged(bs,groups)
    return bs,edges,anchors

def run(native,out,scope='all_edges'):
    out.mkdir(parents=True,exist_ok=True);cases=0;anchors=0;native_pairs={};max_native_root_error=0.
    af=(out/'model-anchors.csv').open('w',newline='',encoding='utf-8');aw=csv.writer(af)
    aw.writerow(['scenario','trial','basin','block','scale','candidate','center','source_kind','side','apex_index','observed_start','observed_end','left_slope','right_slope','b0','b1','b2','b3'])
    for identity,z,candidates,intervals in fields(native/'fields.bin'):
        bs,edges,count=verify(candidates,z,scope);cases+=1;anchors+=count
        for b in bs:
            for a in b.get('model_anchors',[]):
                r=candidates[a['candidate_id']]
                aw.writerow([*identity,b['id'],b['block'],b['scale'],a['candidate_id'],a['center'],b['source_kind'],a['side'],a['apex'],*a['data_support'],a['left_slope'],a['right_slope'],*[r[k] for k in ('b0','b1','b2','b3')]])
        if identity[1]==0 and (0<=identity[0]<10 or 1000<=identity[0]<1010):native_pairs[identity[0]]=(len(z),bs,edges)
        print(f'verified case={identity} native_model_anchors={count}',flush=True)
    pairs=[]
    for s in range(10):
        if s not in native_pairs or s+1000 not in native_pairs:continue
        n,a,ae=native_pairs[s];_,b,be=native_pairs[s+1000]
        ac=canonical(a,ae);bc=canonical(b,be,n);assert ac[:2]==bc[:2],s
        roots_a={ac[2][r['id']]:r['root'] for r in a if r['kind']=='boundary_model_maximum'}
        roots_b={bc[2][r['id']]:tuple(n-1-v for v in r['root'][::-1]) for r in b if r['kind']=='boundary_model_maximum'}
        assert set(roots_a)==set(roots_b)
        error=max((abs(x-y) for k in roots_a for x,y in zip(roots_a[k],roots_b[k])),default=0.)
        max_native_root_error=max(max_native_root_error,error);pairs.append(dict(scenario=s,exact_structural_reflection=True,max_native_apex_reflection_error=error))
    af.close()
    result=dict(boundary_scope=scope,cases=cases,model_anchors_verified=anchors,every_native_fit_has_one_owner=True,all_v1_maximum_nodes_and_edges_preserved=True,all_v1_multimodel_bundles_preserved=True if scope=='singleton_blocks' else None,reflected_coefficient_geometry_verified_cases=cases,actual_native_reflection_pairs=pairs,max_actual_native_boundary_apex_reflection_error=max_native_root_error,numerical_sign_certificate=False,statistical_peak_acceptance=False)
    (out/'audit.json').write_text(json.dumps(result,indent=2),encoding='utf-8');return result

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('native',type=Path);p.add_argument('output',type=Path);p.add_argument('--scope',choices=('all_edges','singleton_blocks'),default='all_edges');a=p.parse_args();print(json.dumps(run(a.native,a.output,a.scope),indent=2))
