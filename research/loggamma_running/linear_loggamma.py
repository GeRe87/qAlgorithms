"""Research-only pooled log-Gamma peak fitting with cached running moments.

No imports from qAlgorithms and no changes to the production pipeline.
Coordinates z are dimensionless local native-axis offsets. One lambda and
orientation are shared within a pooled fit.
"""
from __future__ import annotations

from dataclasses import dataclass
from math import exp, fsum, isfinite, lgamma, log
from typing import Sequence


@dataclass(frozen=True)
class Moments:
    y0: float
    yz: float
    ye: float
    yy: float


@dataclass(frozen=True)
class Geometry:
    z: tuple[float, ...]
    e: tuple[float, ...]
    lam: float
    sign: int
    ginv00: float
    ginv01: float
    ginv11: float
    hz0: float
    hz1: float
    schur: float
    zsum: float
    zesum: float
    esum: float
    n: int

    @classmethod
    def from_offsets(cls, offsets: Sequence[float], lam: float, sign: int) -> Geometry:
        if not isfinite(lam) or lam <= 0 or sign not in (-1, 1):
            raise ValueError('Need lambda > 0 and sign in {-1,+1}')
        z = tuple(float(v) for v in offsets)
        if len(z) < 4 or not all(isfinite(v) for v in z):
            raise ValueError('At least four finite coordinates are required')
        e = tuple(exp(sign * lam * v) for v in z)
        n = len(z)
        esum = fsum(e)
        eesum = fsum(v * v for v in e)
        zsum = fsum(z)
        zesum = fsum(zi * ei for zi, ei in zip(z, e))
        zzsum = fsum(v * v for v in z)
        det = n * eesum - esum * esum
        if det <= 1e-13 * n * eesum:
            raise ValueError('Unidentifiable local intercept/exponential basis')
        g00, g01, g11 = eesum/det, -esum/det, n/det
        hz0 = g00*zsum + g01*zesum
        hz1 = g01*zsum + g11*zesum
        schur = zzsum - zsum*hz0 - zesum*hz1
        if schur <= 1e-11 * max(zzsum, 1.):
            raise ValueError('Shared log-Gamma slope not identifiable')
        return cls(z,e,lam,sign,g00,g01,g11,hz0,hz1,schur,zsum,zesum,esum,n)

    @classmethod
    def uniform(cls, n: int, step: float, lam: float, sign: int) -> Geometry:
        if not isfinite(step) or step <= 0:
            raise ValueError('Need positive finite step')
        return cls.from_offsets([(i-(n-1)/2)*step for i in range(n)],lam,sign)

    def moments(self, y: Sequence[float]) -> Moments:
        if len(y) != self.n or not all(isfinite(float(v)) for v in y):
            raise ValueError('Log-intensity vector has wrong length/nonfinite values')
        return Moments(fsum(y),fsum(z*v for z,v in zip(self.z,y)),
                       fsum(e*v for e,v in zip(self.e,y)),fsum(v*v for v in y))

    def slide_right(self, m: Moments, outgoing: float, incoming: float) -> Moments:
        """Constant-time update for a uniform window advanced one index right.

        Prefer for sign=+1. For sign=-1 use slide_left to avoid
        amplification of roundoff in the exponential-weighted sum.
        """
        self._require_uniform()
        h = self.z[1]-self.z[0]
        return Moments(m.y0-outgoing+incoming,
                       m.yz-self.z[0]*outgoing-h*(m.y0-outgoing)+self.z[-1]*incoming,
                       exp(-self.sign*self.lam*h)*(m.ye-self.e[0]*outgoing)+self.e[-1]*incoming,
                       m.yy-outgoing*outgoing+incoming*incoming)

    def slide_left(self, m: Moments, outgoing: float, incoming: float) -> Moments:
        """Constant-time update for a uniform window advanced one index left.

        Prefer for sign=-1.
        """
        self._require_uniform()
        h = self.z[1]-self.z[0]
        return Moments(m.y0-outgoing+incoming,
                       m.yz-self.z[-1]*outgoing+h*(m.y0-outgoing)+self.z[0]*incoming,
                       exp(self.sign*self.lam*h)*(m.ye-self.e[-1]*outgoing)+self.e[0]*incoming,
                       m.yy-outgoing*outgoing+incoming*incoming)

    def _require_uniform(self) -> None:
        h = self.z[1]-self.z[0]
        if any(abs((self.z[i]-self.z[i-1])-h)>1e-10*max(1,abs(h)) for i in range(2,self.n)):
            raise ValueError('Running update requires uniform window geometry')


@dataclass(frozen=True)
class PeakEstimate:
    intercept: float
    exponential_coefficient: float
    apex_local: float
    apex_global: float
    log_area_local: float
    curvature_at_apex: float


@dataclass(frozen=True)
class FoldingFit:
    shared_slope: float
    peaks: tuple[PeakEstimate, ...]
    residual_sum_squares: float
    degrees_of_freedom: int
    mse: float


def fit_pooled(profiles: Sequence[tuple[Geometry,Moments,float]]) -> FoldingFit:
    """Pooled OLS via local 2x2 eliminations and one scalar Schur step.

    Each entry is (cached_geometry, running_moments, native_axis_anchor).
    Different profile lengths are allowed, with shared lambda and sign.
    Unconstrained OLS coefficients not defining an integrable peak raise.
    """
    if not profiles:
        raise ValueError('Need at least one profile')
    lam=profiles[0][0].lam
    sign=profiles[0][0].sign
    if any(g.lam!=lam or g.sign!=sign for g,_,_ in profiles):
        raise ValueError('One common lambda and orientation per pooled fit')
    rhs=0.0
    schur=0.0
    for g,m,_ in profiles:
        p0=g.ginv00*m.y0+g.ginv01*m.ye
        p1=g.ginv01*m.y0+g.ginv11*m.ye
        rhs += m.yz-g.zsum*p0-g.zesum*p1
        schur += g.schur
    d=rhs/schur
    if not isfinite(d) or sign*d<=0:
        raise ValueError('OLS slope does not give an integrable peak')
    estimates=[]
    rss=0.0
    for g,m,anchor in profiles:
        ty0=m.y0-d*g.zsum
        tye=m.ye-d*g.zesum
        a=g.ginv00*ty0+g.ginv01*tye
        c=g.ginv01*ty0+g.ginv11*tye
        if not isfinite(c) or c>=0:
            raise ValueError('OLS exponential coefficient must be negative')
        shape=d/(sign*lam)
        zmax=log(-d/(sign*lam*c))/(sign*lam)
        logarea=a-log(lam)+lgamma(shape)-shape*log(-c)
        estimates.append(PeakEstimate(a,c,zmax,anchor+zmax,logarea,-sign*d*lam))
        rss +=m.yy-a*m.y0-c*m.ye-d*m.yz
    dof=sum(g.n for g,_,_ in profiles)-2*len(profiles)-1
    if dof<=0:
        raise ValueError('Not enough residual degrees of freedom')
    if rss< -1e-7*sum(m.yy for _,m,_ in profiles):
        raise ArithmeticError('Severe loss of precision in residual computation')
    rss=max(rss,0.0)
    return FoldingFit(d,tuple(estimates),rss,dof,rss/dof)


def log_peak_height(z: float, apex: float, log_height: float,
                    shape: float, lam: float, sign: int = 1) -> float:
    """Generate an exactly shifted synthetic log-Gamma peak."""
    t=sign*lam*(z-apex)
    return log_height+shape*(t+1-exp(t))
