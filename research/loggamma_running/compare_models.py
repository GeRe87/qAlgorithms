"""Synthetic peak-family comparison (optional research deps: numpy, scipy).

Not a mass-accuracy validation; the instrument physics are not represented.
Run: pip install numpy scipy ; python compare_models.py
"""
import numpy as np
from scipy.integrate import quad
from scipy.special import betaln, kve, ndtr
from linear_loggamma import Geometry, fit_pooled

LAMBDA = .8
Q = 9.
P = Q / LAMBDA**2


def truth_log(kind, z, height):
    if kind == 'gaussian':
        return np.log(height) - Q*z*z/2
    if kind == 'loggamma+':
        return np.log(height) + P*(LAMBDA*z+1-np.exp(LAMBDA*z))
    if kind == 'loggamma-':
        return np.log(height) + P*(-LAMBDA*z+1-np.exp(-LAMBDA*z))
    if kind == 'logcosh':
        return np.log(height) - P*np.log(np.cosh(LAMBDA*z))
    if kind == 'hyperbolic':
        return np.log(height) - (Q/LAMBDA**2)*(np.cosh(LAMBDA*z)-1) - (3/LAMBDA**2)*(np.sinh(LAMBDA*z)-LAMBDA*z)
    if kind == 'piecewise':
        return np.log(height) - .5*np.where(z < 0, 7, 12)*z*z
    raise ValueError(kind)


def profiles_for(kind, sigma=.015, seed=42):
    rng = np.random.default_rng(seed)
    positions = np.linspace(-.4, .4, 8)
    heights = 700 + 1100*rng.uniform(size=8)
    profiles=[]
    for mu, height in zip(positions, heights):
        phase = rng.uniform(-.04,.04)
        x = mu + phase + np.linspace(-.65,.65,31)
        logy = truth_log(kind,x-mu,height) + rng.normal(0,sigma,len(x))
        anchor = x[np.argmax(logy)]  # observed maximum, NOT true apex
        profiles.append((x-anchor,logy,anchor,mu,height))
    return profiles


def pooled_ols(profiles, local_functions, shared_functions):
    local_grams=[]
    common_schur=None
    common_rhs=None
    for z,y,*_ in profiles:
        L = np.column_stack([f(z) for f in local_functions])
        Z = np.column_stack([f(z) for f in shared_functions])
        gram = L.T@L
        h = L.T@Z
        rhs = L.T@y
        gin_h = np.linalg.solve(gram,h)
        gin_y = np.linalg.solve(gram,rhs)
        sj = Z.T@Z - h.T@gin_h
        rj = Z.T@y - h.T@gin_y
        common_schur = sj if common_schur is None else common_schur+sj
        common_rhs = rj if common_rhs is None else common_rhs+rj
        local_grams.append((gram,h,rhs))
    shared = np.linalg.solve(common_schur,common_rhs)
    local = [np.linalg.solve(gram,rhs-h@shared) for gram,h,rhs in local_grams]
    residual = np.concatenate([y - np.column_stack([f(z) for f in local_functions])@b -
                               np.column_stack([f(z) for f in shared_functions])@shared
                               for (z,y,*_),b in zip(profiles,local)])
    return np.array(local), shared, float(np.sqrt(np.mean(residual**2)))


def fit_kind(data, candidate):
    k=LAMBDA
    if candidate in ('loggamma+','loggamma-'):
        s = 1 if candidate[-1]=='+' else -1
        entries=[]
        for z,y,anchor,*_ in data:
            g=Geometry.from_offsets(z,k,s)
            entries.append((g,g.moments(y),anchor))
        fit=fit_pooled(entries)
        return (np.array([p.apex_global for p in fit.peaks]),
                np.exp([p.log_area_local for p in fit.peaks]),
                np.sqrt(fit.residual_sum_squares/sum(g.n for g,_,_ in entries)))
    if candidate == 'logcosh':
        loc, shared, rms = pooled_ols(data,[lambda z: np.ones_like(z),lambda z:z],
                                  [lambda z:np.log(np.cosh(k*z))])
        c=-shared[0]
        if c<=0 or np.any(np.abs(loc[:,1]) >= c*k):
            raise ValueError('Nonintegrable log-cosh fit')
        mus=np.array([anchor+np.arctanh(b/(c*k))/k for (a,b),(_,_,anchor,*_) in zip(loc,data)])
        logarea=np.array([a+(c-1)*np.log(2)-np.log(k)+betaln((c+b/k)/2,(c-b/k)/2) for a,b in loc])
        return mus,np.exp(logarea),rms
    if candidate == 'hyperbolic':
        loc, shared, rms = pooled_ols(data,[lambda z:np.ones_like(z),lambda z:z],
             [lambda z:(np.cosh(k*z)-1)/k**2,lambda z:(np.sinh(k*z)-k*z)/k**2])
        q,r=-shared
        if q<=abs(r):
            raise ValueError('Nonintegrable hyperbolic fit')
        Qh,R=q/k**2,r/k**2
        mus=[]; areas=[]
        for (a,b),(_,_,anchor,*_) in zip(loc,data):
            v=R+b/k
            t=(v+np.sqrt(v*v+Qh*Qh-R*R))/(Qh+R)
            mus.append(anchor+np.log(t)/k)
            nu=b/k+R
            lam=np.sqrt(Qh*Qh-R*R)
            logarea=a+Qh+np.log(2/k)+.5*nu*np.log((Qh-R)/(Qh+R))+np.log(kve(nu,lam))-lam
            areas.append(np.exp(logarea))
        return np.array(mus),np.array(areas),rms
    if candidate == 'piecewise':
        loc, shared, rms=pooled_ols(data,[lambda z:np.ones_like(z),lambda z:z],
                [lambda z:np.where(z<0,.5*z*z,0),lambda z:np.where(z>=0,.5*z*z,0)])
        qleft,qright=-shared
        if qleft<=0 or qright<=0:
            raise ValueError('Nonintegrable piecewise fit')
        mus=[];areas=[]
        for (a,b),(_,_,anchor,*_) in zip(loc,data):
            mus.append(anchor+b/(qleft if b<0 else qright))
            area=np.exp(a)*(np.exp(b*b/(2*qleft))*np.sqrt(2*np.pi/qleft)*ndtr(-b/np.sqrt(qleft))+
                          np.exp(b*b/(2*qright))*np.sqrt(2*np.pi/qright)*ndtr(b/np.sqrt(qright)))
            areas.append(area)
        return np.array(mus),np.array(areas),rms
    raise ValueError(candidate)


def main():
    kinds=['gaussian','loggamma+','loggamma-','logcosh','hyperbolic','piecewise']
    methods=['loggamma+','loggamma-','logcosh','hyperbolic','piecewise']
    print('true_shape,fit_model,rms_log,apex_mae_u,mean_abs_area_error_percent,valid')
    for kind in kinds:
        data=profiles_for(kind)
        true_apex=np.array([d[3] for d in data])
        true_area=np.array([quad(lambda z:np.exp(truth_log(kind,z,d[4])), -12,12)[0] for d in data])
        for model in methods:
            try:
                apex,area,rms=fit_kind(data,model)
                valid=np.isfinite(apex)&np.isfinite(area)
                mae=np.mean(np.abs(apex[valid]-true_apex[valid])) if valid.any() else np.nan
                error=100*np.mean(np.abs(area[valid]/true_area[valid]-1)) if valid.any() else np.nan
                print(f'{kind},{model},{rms:.5f},{mae:.5f},{error:.3f},{valid.sum()}')
            except Exception as error:
                print(f'{kind},{model},INVALID,{error}')


if __name__ == '__main__':
    main()
