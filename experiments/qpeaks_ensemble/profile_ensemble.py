"""Model averaging on native windows, with no new fit of the observations."""

from functools import lru_cache
import math

import numpy as np


def arrays(models):
    return np.array([[m[k] for k in ("center", "scale", "b0", "b1", "b2", "b3")]
                     for m in models], dtype=float)


def terms(models, x):
    center, scale, b0, b1, b2, b3 = models.T[:, :, None]
    d = x - center
    u = d / scale
    taper = np.maximum(0., 1. - u*u)
    w = taper*taper
    dw = -4*u*taper/scale
    curvature = np.where(d < 0, b2, b3)
    q = b0 + d*(b1 + d*curvature)
    return w, dw, q, b1 + 2*d*curvature


def evaluate(models, x):
    x = np.atleast_1d(x)
    sums = np.zeros((4, len(x)))
    for start in range(0, len(models), 128):
        w, dw, q, dq = terms(models[start:start+128], x)
        sums += np.array([w.sum(axis=0), (w*q).sum(axis=0),
                          dw.sum(axis=0), (dw*q + w*dq).sum(axis=0)])
    den, num, dden, dnum = sums
    value = np.divide(num, den, out=np.full(len(x), np.nan), where=den > 0)
    slope = np.divide(dnum - value*dden, den, out=np.full(len(x), np.nan), where=den > 0)
    return value, slope, den


def crossings(slope, first, direction):
    """Keep strict sign transitions, including crossings through zero runs."""
    nonzero = np.flatnonzero(slope != 0)
    return [(first + int(a), first + int(b)) for a, b in zip(nonzero, nonzero[1:])
            if direction*slope[a] > 0 > direction*slope[b]]


def describe(models):
    from scipy.optimize import brentq
    models = arrays(models)
    c, s = models.T[:2]
    first, last = int(min(c-s))+1, int(max(c+s))-1
    grid = np.arange(first, last+1)
    log_profile, slope, coverage = evaluate(models, grid)
    if np.any(coverage <= 0):
        return "disconnected_supported_grid", []
    maxima = crossings(slope, first, 1)
    if not maxima:
        return "no_bracketed_maximum", []
    # Between adjacent maxima, use the lowest evaluated scan as the separator.
    separators = []
    for (_, a), (b, _) in zip(maxima, maxima[1:]):
        values = log_profile[a-first:b-first+1]
        lowest = np.flatnonzero(values == values.min())
        if len(lowest) != 1:
            return "ambiguous_valley_boundary", []
        separators.append(a + int(lowest[0]))
    boundaries = [first, *separators, last]
    output = []
    for (a, b), left, right in zip(maxima, boundaries, boundaries[1:]):
        local = models[(c+s > a) & (c-s < b)]
        apex = brentq(lambda x: evaluate(local, x)[1][0], a, b)
        height = float(np.exp(evaluate(local, apex)[0][0]))
        intensities = np.exp(log_profile[left-first:right-first+1])
        output.append(dict(apex=apex, intensity=height, left=left, right=right,
                           **shape(models, apex, height, left, right, intensities)))
    return "multiple_maxima" if len(output) > 1 else "one_maximum", output


def shape(models, apex, height, left, right, intensities):
    from scipy.optimize import brentq
    if not left < apex < right:
        return dict(area_scan=None, asymmetry=None, shape_status="apex_at_boundary")
    lo, hi = intensities[[0, -1]]
    chord = np.linspace(lo, hi, len(intensities))
    excess = intensities - chord
    area = float(excess.sum())  # Unit-grid trapezoid, with zero endpoint excess.
    def relative_height(x):
        baseline = lo + (hi-lo)*(x-left)/(right-left)
        return math.exp(evaluate(models, x)[0][0]) - baseline
    half = (height - (lo + (hi-lo)*(apex-left)/(right-left)))/2
    if half <= 0:
        return dict(area_scan=None, asymmetry=None, shape_status="nonpositive_peak_above_chord")
    # Nearest observed brackets of the half-height crossings.
    lower = np.flatnonzero((np.arange(left, right+1) < apex) & (excess <= half))[-1] + left
    upper = np.flatnonzero((np.arange(left, right+1) > apex) & (excess <= half))[0] + left
    hleft = brentq(lambda x: relative_height(x)-half, lower, min(lower+1, apex))
    hright = brentq(lambda x: relative_height(x)-half, max(upper-1, apex), upper)
    return dict(area_scan=area, asymmetry=(hright-apex)/(apex-hleft),
                half_left=hleft, half_right=hright, shape_status="chord_diagnostic")


@lru_cache(maxsize=None)
def kernel(scale):
    d = np.arange(-scale, scale+1, dtype=float)/scale
    design = np.column_stack([np.ones_like(d), d, np.minimum(d, 0)**2, np.maximum(d, 0)**2])
    return np.linalg.pinv(design)


def influence(models, x, first, length):
    """Log-profile and slope Jacobians conditional on fixed members."""
    w, dw, _, _ = terms(models, np.array([x]))
    w, dw = w[:, 0], dw[:, 0]
    total = w.sum()
    if total <= 0:
        raise ValueError("No profile support at the evaluation point.")
    value, slope = np.zeros((2, length))
    for j in np.flatnonzero(w):
        center, scale = map(int, models[j, :2])
        u = (x-center)/scale
        basis = np.array([1., u, min(u, 0)**2, max(u, 0)**2])
        derivative = np.array([0., 1., 2*min(u, 0), 2*max(u, 0)])/scale
        window = slice(center-scale-first, center+scale-first+1)
        a = basis @ kernel(scale)
        weight = w[j]/total
        value[window] += weight*a
        slope[window] += weight*(derivative @ kernel(scale)) + (dw[j]-weight*dw.sum())/total*a
    return value, slope


def conditional_gradients(models, peak):
    """Apex, chord area and width ratio derivatives w.r.t. original log scans."""
    models = arrays(models)
    c, s = models.T[:2]
    first, last = int(min(c-s)), int(max(c+s))
    def at(x):
        return influence(models, x, first, last-first+1)
    apex = peak["apex"]
    if np.any((apex == c) & (models[:, 4] != models[:, 5]) & (np.abs(apex-c) < s)):
        raise ValueError("Curvature may be discontinuous at an active hinge.")
    if np.any((apex == c-s) | (apex == c+s)):
        raise ValueError("Apex lies at a native window boundary.")
    w, dw, q, dq = (column[:, 0] for column in terms(models, np.array([apex])))
    u = (apex-c)/s
    ddw = np.where(np.abs(u) < 1, (12*u*u-4)/(s*s), 0.)
    ddq = 2*np.where(u < 0, models[:, 4], models[:, 5])
    value, slope, den = evaluate(models, apex)
    curvature = (sum(ddw*(q-value[0]) + 2*dw*(dq-slope[0]) + w*ddq)/den[0])
    if curvature >= 0 or not np.isfinite(curvature):
        raise ValueError("Apex curvature does not support linear propagation.")
    value, derivative = at(apex)
    da = -derivative/curvature
    result = dict(apex=da, intensity=peak["intensity"]*value)
    if peak["area_scan"] is None:
        return result
    left, right = peak["left"], peak["right"]
    grid = np.arange(left, right+1)
    intensities = np.exp(evaluate(models, grid)[0])
    area = sum(intensity*at(x)[0] for x, intensity in zip(grid[1:-1], intensities[1:-1]))
    dlo = intensities[0]*at(left)[0]
    dhi = intensities[-1]*at(right)[0]
    result["area_scan"] = area - (right-left-1)/2*(dlo+dhi)
    chord_slope = (intensities[-1]-intensities[0])/(right-left)
    def dexcess(x, intensity):
        fraction = (x-left)/(right-left)
        return intensity*at(x)[0] - ((1-fraction)*dlo + fraction*dhi)
    dhalf = (dexcess(apex, peak["intensity"]) - chord_slope*da)/2
    def crossing_gradient(x):
        log_value, log_slope, _ = evaluate(models, x)
        intensity = math.exp(log_value[0])
        slope = intensity*log_slope[0] - chord_slope
        if slope == 0:
            raise ValueError("Half-height crossing is tangent.")
        return (dhalf-dexcess(x, intensity))/slope
    l, r = peak["half_left"], peak["half_right"]
    result["asymmetry"] = ((crossing_gradient(r)-da)
                           - peak["asymmetry"]*(da-crossing_gradient(l)))/(apex-l)
    return result
