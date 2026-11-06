"""Probe valley contrasts on overlapping Gaussian signals with known log noise."""

import argparse
from collections import Counter
import json
from pathlib import Path

import numpy as np

from audit_profile_corridors import valley
from bundle_method import topology
from profile_ensemble import describe
from study_profile_uncertainty import field


def run(trials, output):
    rng = np.random.default_rng(251311)
    x = np.arange(41)-20
    width, sigma = 2.5, .05
    report = []
    for separation in (0., 4., 6., 9.):
        intensity = sum(300*np.exp(-.5*((x-sign*separation/2)/width)**2) for sign in (-1, 1))
        # A constant baseline makes the log observations finite on the tails.
        truth = np.log(20 + intensity)
        counts = Counter()
        cases = []
        for trial in range(trials):
            y = truth + rng.normal(0., sigma, len(x))
            models = field(y, range(2, 13))
            basins, _, groups, _ = topology(models)
            result = []
            for group in groups:
                ids = [i for b in group["basins"] for i in basins[b]["candidate_ids"]]
                members = [models[i] for i in ids]
                status, peaks = describe(members)
                valleys = [valley(members, a, b, sigma) for a, b in zip(peaks, peaks[1:])]
                result.append(dict(status=status, apices=[p["apex"] for p in peaks], valleys=valleys))
            pairs = [v for g in result for v in g["valleys"]]
            counts["cases_with_multiple_branches"] += sum(bool(g["apices"]) for g in result) > 1
            counts["cases_with_within_branch_valley"] += bool(pairs)
            # Nominal pointwise contrasts are diagnostics after searching for maxima.
            counts["cases_with_both_drops_above_1_96_se"] += any(
                all(d > 1.96*se for d, se in zip(v["log_drops"], v["conditional_log_drop_se"]))
                for v in pairs)
            cases.append(dict(trial=trial, branches=result))
        report.append(dict(separation_scans=separation, gaussian_width_scans=width,
                           known_sigma_log=sigma, trials=trials,
                           underlying_profile_modes=1 if separation < 2*width else 2,
                           counts=counts, cases=cases))
    output.write_text(json.dumps(report, allow_nan=False), encoding="utf-8")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--trials", type=int, default=64)
    args = parser.parse_args()
    run(args.trials, args.output)
