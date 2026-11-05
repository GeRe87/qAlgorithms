"""Research reference: native slope-field basins, scale DAG, ensemble-derived fit.

Topology proposals are not statistical peak detections. No cascade annotations
are read. Exact floating signs are used; numerical sign certification is pending.
"""
import math
from collections import defaultdict


def basins(centers, slopes):
    """Partition one contiguous center grid by ascent in its local slope field.

    A +,- pair is a maximum corridor, not an interpolated measurement. Exact
    zero plateaus preserve their topology; boundary exits remain explicit.
    """
    n=len(centers)
    if not n or len(slopes)!=n or any(centers[i+1]!=centers[i]+1 for i in range(n-1)):
        raise ValueError('A nonempty contiguous center grid is required.')
    if not all(math.isfinite(d) for d in slopes):raise ValueError('Nonfinite slopes.')
    directions=[int(d>0)-int(d<0) for d in slopes]
    terminal={};records=[]
    def add(indices,kind,root):
        ident=len(records);records.append(dict(kind=kind,root=root,positions=[]))
        for i in indices:terminal[i]=ident
        return ident
    i=0
    while i<n:
        if directions[i]:i+=1;continue
        j=i+1
        while j<n and not directions[j]:j+=1
        left=directions[i-1] if i else None;right=directions[j] if j<n else None
        if left==1 and right==-1:add(range(i,j),'maximum',(centers[i],centers[j-1]))
        elif left is not None and left==right:
            for k in range(i,j):directions[k]=left
        else:
            kind='stationary_minimum' if left==-1 and right==1 else 'flat_or_boundary_stationary'
            add(range(i,j),kind,None)
        i=j
    for i in range(n-1):
        if directions[i]==1 and directions[i+1]==-1:
            add((i,i+1),'maximum',(centers[i],centers[i+1]))
    exits={}
    for first in range(n):
        path=[];i=first
        while i not in terminal:
            path.append(i);step=directions[i]
            if step==0:raise AssertionError('Unclassified stationary center.')
            nxt=i+step
            if nxt<0 or nxt>=n:
                side='left' if nxt<0 else 'right'
                if side not in exits:exits[side]=add((),side+'_boundary_exit',None)
                terminal[i]=exits[side];break
            i=nxt
        owner=terminal[i]
        for p in path:terminal[p]=owner
    for i,c in enumerate(centers):records[terminal[i]]['positions'].append(c)
    assert sum(len(r['positions']) for r in records)==n
    return records,[terminal[i] for i in range(n)]


def topology(candidates):
    """Every native fit belongs to one basin; retain every maximum DAG edge."""
    grids=defaultdict(list)
    for r in candidates:grids[(int(r['block']),int(r['scale']))].append(r)
    result=[];membership={};lookup={}
    for (block,scale),rs in sorted(grids.items()):
        rs=sorted(rs,key=lambda r:int(r['center']));centers=[int(r['center']) for r in rs]
        records,owners=basins(centers,[float(r['b1']) for r in rs]);offset=len(result)
        for i,r in enumerate(records):
            result.append(dict(r,id=offset+i,block=block,scale=scale,candidate_ids=[]))
        for r,owner in zip(rs,owners):
            ident=int(r['id'])
            if ident in membership:raise ValueError('Duplicate native fit ID.')
            owner+=offset;membership[ident]=owner;result[owner]['candidate_ids'].append(ident)
            if result[owner]['kind']=='maximum':lookup[(block,scale,int(r['center']))]=owner
    edges=defaultdict(int)
    for (block,scale,center),owner in lookup.items():
        nxt=lookup.get((block,scale+1,center))
        if nxt is not None:edges[(owner,nxt)]+=1
    incoming=defaultdict(set);outgoing=defaultdict(set)
    for a,b in edges:outgoing[a].add(b);incoming[b].add(a)
    parents={r['id']:r['id'] for r in result if r['kind']=='maximum'}
    def root(i):
        while parents[i]!=i:parents[i]=parents[parents[i]];i=parents[i]
        return i
    # Contract only unique continuations. Merges and splits are never erased.
    for a,b in edges:
        if len(outgoing[a])==len(incoming[b])==1:
            aa,bb=root(a),root(b);parents[max(aa,bb)]=min(aa,bb)
    bundles=defaultdict(list)
    for i in parents:bundles[root(i)].append(i)
    groups=[dict(id=k,basins=sorted(v)) for k,v in sorted(bundles.items())]
    for g in groups:
        members=set(g['basins'])
        g['incoming']=[list(e) for e in edges if e[1] in members and e[0] not in members]
        g['outgoing']=[list(e) for e in edges if e[0] in members and e[1] not in members]
        g['native_fits']=sum(len(result[i]['candidate_ids']) for i in members)
        g['contains_merge']=any(len(incoming[i])>1 for i in members)
        g['contains_split']=any(len(outgoing[i])>1 for i in members)
    assert len(membership)==len(candidates)
    return result,dict(edges),groups,membership


class MomentGeometry:
    """Cached geometry, prefix-moment coefficients, explicit reference scoring.

    QR builds the data-independent inverse once for each (n,h). This extends
    the same four-column native model to asymmetric supports. Explicit curve
    RSS is intentionally retained here to audit cancellation in fast RSS.
    """
    def __init__(self,np):self.np=np;self.cache={}
    def geometry(self,n,h):
        key=n,h
        if key not in self.cache:
            np=self.np;t=(np.arange(n)-h)/(n-1)
            X=np.column_stack([np.ones(n),t,np.where(t<0,t*t,0),np.where(t>0,t*t,0)])
            norm=np.linalg.norm(X,axis=0)
            if n<5 or not 0<h<n-1 or np.any(norm==0):raise ValueError('Unsupported geometry.')
            G=X/norm;_,R=np.linalg.qr(G,mode='reduced');inv=np.linalg.solve(R,np.eye(4))
            self.cache[key]=(inv@inv.T,norm)
        return self.cache[key]
    def fit(self,z,hinges):
        np=self.np;z=np.asarray(z,dtype=float);n=len(z)
        if n<5 or not np.isfinite(z).all():raise ValueError('Invalid original log observations.')
        if not hinges:raise ValueError('No ensemble hinge proposals.')
        anchor=float(z[n//2] if n%2 else (z[n//2-1]+z[n//2])/2);w=z-anchor
        x=np.arange(n);prefix=[np.r_[0.,np.cumsum(w*x**k)] for k in (0,1,2)]
        records=[]
        for h in sorted(set(hinges)):
            if int(h)!=h:raise ValueError('Reference uses enumerated integer hinges.')
            h=int(h);M,norm=self.geometry(n,h);span=n-1
            p0=prefix[0][-1];p1=(prefix[1][-1]-h*p0)/span
            left=(prefix[2][h]-2*h*prefix[1][h]+h*h*prefix[0][h])/span**2
            total=(prefix[2][-1]-2*h*prefix[1][-1]+h*h*p0)/span**2
            rhs=np.array([p0,p1,left,total-left]);beta=(M@(rhs/norm))/norm
            t=(x-h)/span;curve=beta[0]+beta[1]*t+np.where(t<0,beta[2],beta[3])*t*t
            residual=curve-w;rss=float(residual@residual);fast_rss=float(w@w-beta@rhs)
            actual=beta.copy();actual[0]+=anchor;actual[1]/=span;actual[2:]/=span**2
            records.append(dict(hinge=h,rss=rss,fast_rss=fast_rss,coefficients=actual.tolist()))
        best=min(r['rss'] for r in records);winners=[r['hinge'] for r in records if r['rss']==best]
        return records,winners


def bundle_support(group,basin_rows,candidates):
    """Support union and hinge proposals come from the native ensemble itself."""
    ids={i for b in group['basins'] for i in basin_rows[b]['candidate_ids']}
    rs=[r for r in candidates if int(r['id']) in ids]
    intervals=sorted((int(r['center'])-int(r['scale']),int(r['center'])+int(r['scale'])) for r in rs)
    left,right=intervals[0]
    for lo,hi in intervals[1:]:
        if lo>right+1:raise ValueError('Disconnected support must remain separate.')
        right=max(right,hi)
    hinges=sorted({int(r['center'])-left for r in rs})
    return left,right,hinges


def geometric_peak(coefficients,hinge,n):
    """Describe a strict interior local maximum; never assert peak detection."""
    b0,b1,b2,b3=coefficients;root=None
    if b1<0 and b2<0:
        a=-b1/(2*b2)
        if -hinge<a<0:root=a
    elif b1>0 and b3<0:
        a=-b1/(2*b3)
        if 0<a<n-1-hinge:root=a
    elif b1==0 and b2<0 and b3<0:root=0.
    if root is None:return dict(apex=None,log_height=None,above_both_support_edges=None)
    value=lambda t:b0+b1*t+(b2 if t<0 else b3)*t*t
    height=value(root)
    return dict(apex=hinge+root,log_height=height,above_both_support_edges=height>max(value(-hinge),value(n-1-hinge)))
