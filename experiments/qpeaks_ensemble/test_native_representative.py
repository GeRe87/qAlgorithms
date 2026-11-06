import unittest

from native_representative import native_bundle, representative


class NativeRepresentativeTests(unittest.TestCase):
    def test_preserves_coefficients_and_native_window(self):
        narrow = dict(id=0, center=3, scale=2, b0=2., b1=.5, b2=-.25, b3=-.25)
        wide = dict(narrow, id=1, scale=3, b0=1.)
        z = [2. + .5*x - .25*x*x for x in range(-3, 4)]
        model, rss, ties, beta, peak = native_bundle(z, 0, [wide, narrow])
        self.assertIs(model, narrow)
        self.assertEqual((rss, ties, beta), (0., 1, [2., .5, -.25, -.25]))
        self.assertEqual(peak['apex'], 3.)
        self.assertTrue(peak['above_both_support_edges'])

    def test_points_outside_native_window_count_in_score(self):
        model = dict(center=3, scale=2, b0=0., b1=0., b2=0., b3=0.)
        self.assertEqual(representative([2., 0., 0., 0., 0., 0., 3.], 0, [model])[1], 13.)

    def test_exact_ties_are_explicit(self):
        model = dict(center=2, scale=2, b0=1., b1=0., b2=-1., b3=-1.)
        _, _, ties, _, peak = native_bundle([0.]*5, 0, [model, dict(model)])
        self.assertEqual(ties, 2)
        self.assertIsNone(peak)

    def test_reflection_preserves_score_and_mirrors_apex(self):
        model = dict(center=2, scale=2, b0=4., b1=.5, b2=-.5, b3=-.25)
        reflected = dict(model, center=3, b1=-.5, b2=-.25, b3=-.5)
        z = [0., 1., 3., 2., 1., 0.]
        _, rss, _, _, peak = native_bundle(z, 0, [model])
        _, reverse_rss, _, _, reverse_peak = native_bundle(z[::-1], 0, [reflected])
        self.assertEqual(rss, reverse_rss)
        self.assertEqual(peak['apex'] + 1 + reverse_peak['apex'], 5.)


if __name__ == '__main__':
    unittest.main()
