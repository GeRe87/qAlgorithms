import math
import unittest

from scipy.special import kv

from area_audit import STEP, fitted_areas


class AreaAuditTests(unittest.TestCase):
    def test_symmetric_and_asymmetric_integrable_shapes(self):
        for slope, curvature, skew in [(0., 2., 0.), (1., 2., .4), (20., 2., -.4), (-20., 2., .4)]:
            with self.subTest(slope=slope, skew=skew):
                scale = 4
                intercept = 3.
                grid = [(i - scale) * STEP for i in range(2 * scale + 1)]
                raw = [math.exp(intercept + slope * u - curvature * math.cosh(u) - skew * math.sinh(u))
                       for u in grid]
                observed, sampled, window, whole, tails, error = fitted_areas(raw, scale, slope, curvature, skew)
                p, q = (curvature + skew) / 2, (curvature - skew) / 2
                exact = 2 * math.exp(intercept) * (q / p) ** (slope / 2) * kv(slope, 2 * math.sqrt(p * q))
                self.assertAlmostEqual(whole / exact, 1, places=10)
                self.assertAlmostEqual(sampled / observed, 1, places=13)
                self.assertAlmostEqual(window / whole + tails, 1, places=12)
                self.assertLess(error, 1e-8)

    def test_constant_profile_has_only_window_area(self):
        observed, sampled, window, whole, tails, error = fitted_areas([10.] * 9, 4, 0., 0., 0.)
        self.assertAlmostEqual(observed, 28)
        self.assertAlmostEqual(sampled, observed)
        self.assertAlmostEqual(window, observed)
        self.assertEqual((whole, tails, error), (None, None, None))

    def test_rejects_missing_or_zero_observations(self):
        for raw in [[1.] * 8, [1.] * 8 + [0.]]:
            with self.assertRaises(ValueError):
                fitted_areas(raw, 4, 0., 2., 0.)


if __name__ == "__main__":
    unittest.main()
