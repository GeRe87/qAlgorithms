import unittest

import numpy as np

from audit_profile_corridors import corridor_kind, union_corridors, valley
from profile_ensemble import conditional_gradients, describe, evaluate, arrays, influence
from study_profile_uncertainty import field


class ProfileEnsembleTests(unittest.TestCase):
    def test_slope_includes_change_in_weights(self):
        models = arrays([dict(center=4, scale=3, b0=5., b1=.2, b2=-.1, b3=-.2),
                         dict(center=6, scale=4, b0=6., b1=-.1, b2=-.2, b3=-.1)])
        x, step = 4.5, 1e-5
        values = evaluate(models, [x-step, x+step])[0]
        self.assertAlmostEqual(evaluate(models, x)[1][0], (values[1]-values[0])/(2*step), places=8)
        self.assertAlmostEqual(evaluate(models, 8)[0][0], 5.4)

    def test_reflected_profile_preserves_area_and_inverts_width_ratio(self):
        models = field(5. -.05*(np.arange(21)-10.3)**2, range(2, 6))
        _, (peak,) = describe(models)
        reflected = [dict(m, center=20-m["center"], b1=-m["b1"], b2=m["b3"], b3=m["b2"])
                     for m in reversed(models)]
        _, (reverse,) = describe(reflected)
        self.assertAlmostEqual(peak["apex"] + reverse["apex"], 20.)
        self.assertAlmostEqual(peak["area_scan"], reverse["area_scan"])
        self.assertAlmostEqual(peak["asymmetry"] * reverse["asymmetry"], 1.)

    def test_overlapping_windows_propagate_shared_scan_errors(self):
        x = np.arange(21)
        y = 5. -.04*(x-10.3)**2
        scales = range(2, 6)
        models = field(y, scales)
        _, (peak,) = describe(models)
        gradients = conditional_gradients(models, peak)
        direction = np.random.default_rng(42).normal(size=len(y))
        step = 1e-5
        perturbed = [describe(field(y + sign*step*direction, scales))[1][0] for sign in (-1, 1)]
        for metric, gradient in gradients.items():
            measured = (perturbed[1][metric]-perturbed[0][metric])/(2*step)
            np.testing.assert_allclose(gradient @ direction, measured, rtol=2e-6, atol=1e-7)
        value, slope = influence(arrays(models), peak["apex"], 0, len(y))
        self.assertAlmostEqual(value.sum(), 1.)
        self.assertAlmostEqual(slope.sum(), 0.)


    def test_corridor_attribution_does_not_fill_gaps(self):
        own = [(2, 3), (5, 6)]
        block = union_corridors([*own, (3, 3.5), (7, 8)])
        self.assertEqual(corridor_kind(3, own, block), "own_corridor")
        self.assertEqual(corridor_kind(3.5, own, block), "other_corridor")
        self.assertEqual(corridor_kind(4, own, block), "between_own_corridors")
        self.assertEqual(corridor_kind(9, own, block), "outside_corridors")

    def test_valley_uncertainty_includes_shared_scan_effects(self):
        y = 5. + .04*(np.arange(11)-5)**2
        members = field(y, (3,))
        left, right = dict(apex=4.5, right=5), dict(apex=5.5, left=5)
        sigma = .05
        result = valley(members, left, right, sigma)
        np.testing.assert_allclose(result["log_drops"], [.01, .01])
        step = 1e-5
        numerical = []
        for direction in np.eye(len(y)):
            scores = [valley(field(y+sign*step*direction, (3,)), left, right)["log_drops"]
                      for sign in (-1, 1)]
            numerical.append((np.array(scores[1])-scores[0])/(2*step))
        np.testing.assert_allclose(result["conditional_log_drop_se"],
                                   sigma*np.linalg.norm(numerical, axis=0), rtol=1e-8)


if __name__ == "__main__":
    unittest.main()
