"""Select an existing native regression by log RSS on common bundle support."""

import numpy as np

from bundle_method import geometric_peak


def representative(observations, left, candidates):
    """Return the native model, its common-support score, and exact tie count.

    Scoring can extrapolate. The selected coefficients and geometric peak
    support remain those of the native regression.
    """
    observations = np.asarray(observations, dtype=float)
    positions = np.arange(left, left + len(observations))
    best = None
    best_rss = float("inf")
    ties = 0
    # Bound the working arrays even for long EICs with many candidates.
    for start in range(0, len(candidates), 128):
        batch = candidates[start:start + 128]
        center, b0, b1, b2, b3 = np.array([
            [row[key] for key in ("center", "b0", "b1", "b2", "b3")]
            for row in batch
        ]).T[:, :, None]
        delta = positions - center
        residual = (b0 - observations) + delta * (b1 + delta * np.where(delta < 0, b2, b3))
        squared = residual * residual
        # Pair the two ends so reflection preserves the summation order.
        half = len(observations) // 2
        scores = (squared[:, :half] + squared[:, :len(observations)-half-1:-1]).sum(axis=1)
        if len(observations) % 2:
            scores += squared[:, half]
        winner = int(scores.argmin())
        rss = float(scores[winner])
        count = int(np.count_nonzero(scores == rss))
        if rss < best_rss:
            best, best_rss, ties = batch[winner], rss, count
        elif rss == best_rss:
            ties += count
    return best, best_rss, ties


def native_bundle(observations, left, candidates):
    best, rss, ties = representative(observations, left, candidates)
    center, scale = best["center"], best["scale"]
    coefficients = [best[key] for key in ("b0", "b1", "b2", "b3")]
    peak = geometric_peak(coefficients, scale, 2 * scale + 1) if ties == 1 else None
    return best, rss, ties, coefficients, peak
