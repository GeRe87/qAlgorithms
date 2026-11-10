"""Independent dense-OLS and numerical-integral references (standard library)."""
import math
import random
import unittest

from linear_loggamma import Geometry, fit_pooled, log_peak_height


def dense_ols(rows, targets):
    """Independent full normal equations, dense Gaussian elimination with pivots."""
    p = len(rows[0])
    gram = [[math.fsum(row[i] * row[j] for row in rows)
             for j in range(p)] for i in range(p)]
    rhs = [math.fsum(row[i] * y for row, y in zip(rows, targets))
           for i in range(p)]
    for i in range(p):
        pivot = max(range(i, p), key=lambda j: abs(gram[j][i]))
        gram[pivot], gram[i] = gram[i], gram[pivot]
        rhs[pivot], rhs[i] = rhs[i], rhs[pivot]
        if abs(gram[i][i]) < 1e-13:
            raise ValueError('Dense OLS design is singular')
        for j in range(i + 1, p):
            factor = gram[j][i] / gram[i][i]
            for k in range(i, p):
                gram[j][k] -= factor * gram[i][k]
            rhs[j] -= factor * rhs[i]
    beta = [0.0] * p
    for i in reversed(range(p)):
        beta[i] = (rhs[i] - math.fsum(gram[i][j] * beta[j]
                                   for j in range(i + 1, p))) / gram[i][i]
    rss = math.fsum((y-math.fsum(v * b for v, b in zip(row, beta)))**2
                    for row, y in zip(rows, targets))
    return beta, rss


class IndependentValidationTests(unittest.TestCase):
    def test_pooled_schur_matches_independent_full_ols_with_noise(self):
        rng = random.Random(193)
        for sign in (1, -1):
            lam = 0.63
            parts = []
            rows, targets = [], []
            ns = (19, 26, 33, 42)
            m = len(ns)
            for j, n in enumerate(ns):
                offsets = sorted(-0.7+1.4*i/(n-1)+rng.uniform(-0.009, 0.009)
                                 for i in range(n))
                geom = Geometry.from_offsets(offsets, lam, sign)
                mode = rng.uniform(-0.16, 0.16)
                y = [log_peak_height(z, mode, math.log(600+150*j), 11., lam, sign)
                     + rng.gauss(0, 0.015) for z in geom.z]
                parts.append((geom, geom.moments(y), j*2.))
                for z, e, yi in zip(geom.z, geom.e, y):
                    row = [0.] * (2*m+1)
                    row[2*j], row[2*j+1], row[-1] = 1., e, z
                    rows.append(row)
                    targets.append(yi)
            fitted = fit_pooled(parts)
            expected, reference_rss = dense_ols(rows, targets)
            self.assertAlmostEqual(fitted.shared_slope, expected[-1], delta=5e-7)
            for j, peak in enumerate(fitted.peaks):
                self.assertAlmostEqual(peak.intercept, expected[2*j], delta=5e-7)
                self.assertAlmostEqual(peak.exponential_coefficient, expected[2*j+1], delta=5e-7)
            self.assertAlmostEqual(fitted.residual_sum_squares, reference_rss, delta=1e-7)

    def test_analytic_area_against_independent_simpson_quadrature(self):
        lam, shape, mu, height = 0.7, 8., 0.3, 1500.
        for sign in (1, -1):
            geom = Geometry.uniform(55, 0.032, lam, sign)
            log_h = math.log(height)
            y = [log_peak_height(z, mu, log_h, shape, lam, sign) for z in geom.z]
            fit = fit_pooled([(geom, geom.moments(y), 0.)])
            analytic = math.exp(fit.peaks[0].log_area_local)
            lo, hi, n = mu-8.0, mu+8.0, 20000
            dx = (hi-lo)/n
            numerical = 0.
            for i in range(n+1):
                x = lo+dx*i
                w = (1 if i==0 or i==n else 4 if i%2 else 2)
                numerical += w * math.exp(log_peak_height(x,mu,log_h,shape,lam,sign))
            numerical *= dx/3
            self.assertAlmostEqual(analytic/numerical, 1.0, delta=2e-8)


if __name__ == '__main__':
    unittest.main()
