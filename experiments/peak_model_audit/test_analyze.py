"""Independent identities and simulation controls for the real-data diagnostic."""
import unittest
import numpy as np
from analyze import geometry


class ResidualReferenceTests(unittest.TestCase):
    def test_projection_and_kernel_basis(self):
        for scale in (4, 8):
            design, q, r, operator, lag, dw = geometry(scale)
            np.testing.assert_allclose(q @ r, design, atol=2e-14)
            np.testing.assert_allclose(operator @ design, 0, atol=3e-13)
            np.testing.assert_allclose(operator @ operator, operator, atol=2e-14)
            self.assertAlmostEqual(np.trace(operator), len(design) - 4, places=12)
            self.assertLess(float(np.median(lag)), -.1)
            self.assertTrue(np.all(np.isfinite(dw)))

    def test_projection_noise_is_not_raw_noise(self):
        design, q, r, operator, reference_lag, _ = geometry(4)
        noise = np.random.default_rng(29136).normal(size=(5000, 9))
        # Independent per-draw least-squares residuals, compared to the cached
        # projection operation and its independently seeded reference envelope.
        beta = np.linalg.lstsq(design, noise.T, rcond=None)[0]
        residual = (noise.T - design @ beta).T
        np.testing.assert_allclose(residual, noise @ operator, atol=3e-14)
        lag = np.sum(residual[:, :-1] * residual[:, 1:], axis=1) / np.sum(residual ** 2, axis=1)
        lo, hi = np.quantile(reference_lag, [.025, .975])
        fraction = np.mean((lag < lo) | (lag > hi))
        self.assertTrue(.03 < fraction < .075)
        raw_lag = np.sum(noise[:, :-1] * noise[:, 1:], axis=1) / np.sum(noise ** 2, axis=1)
        self.assertGreater(abs(np.median(lag) - np.median(raw_lag)), .15)


if __name__ == "__main__":
    unittest.main()
