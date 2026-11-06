"""Conditional coverage on exact piecewise peaks with known iid log noise."""

import argparse
from collections import Counter
import json
from pathlib import Path

import numpy as np

from bundle_method import topology
from profile_ensemble import conditional_gradients, describe, kernel


def field(y, scales):
    models = []
    for scale in scales:
        for center in range(scale, len(y)-scale):
            beta = kernel(scale) @ y[center-scale:center+scale+1]
            b0, b1, b2, b3 = beta / np.array([1., scale, scale**2, scale**2])
            models.append(dict(center=center, scale=scale, b0=b0, b1=b1, b2=b2, b3=b3,
                               block=0, id=len(models)))
    return models


def run(trials, output):
    rng = np.random.default_rng(698093)
    report = []
    sigma = 0.05
    for name, right_curvature in (("symmetric", -.045), ("asymmetric", -.025)):
        truth_model = dict(center=14, scale=14, b0=6., b1=.05, b2=-.045, b3=right_curvature)
        x = np.arange(29)-14
        truth = 6. + .05*x + x*x*np.where(x<0, -.045, right_curvature)
        _, (target,) = describe([truth_model])
        metrics = ("apex", "intensity", "area_scan", "asymmetry")
        values = {m: [] for m in metrics}
        status = Counter()
        auto_counts = Counter()
        auto_matches = []
        for _ in range(trials):
            y = truth + rng.normal(0., sigma, len(truth))
            models = field(y, range(2, 10))
            # Oracle membership, fixed across trials and independent of the noise.
            members = [m for m in models if abs(target["apex"]-m["center"]) < m["scale"]]
            outcome, peaks = describe(members)
            status[outcome] += 1
            if len(peaks) == 1:
                peak = peaks[0]
                try:
                    gradients = conditional_gradients(members, peak)
                except ValueError as error:
                    status[str(error)] += 1
                else:
                    for m, gradient in gradients.items():
                        se = sigma*np.linalg.norm(gradient)
                        values[m].append((peak[m]-target[m], se))
            basins, _, groups, _ = topology(models)
            auto = [p for g in groups for p in describe(
                [models[i] for b in g["basins"] for i in basins[b]["candidate_ids"]])[1]]
            auto_counts[len(auto)] += 1
            if len(auto) == 1:
                auto_matches.append(auto[0]["apex"]-target["apex"])
        summary = {}
        for m, pairs in values.items():
            if pairs:
                errors, se = np.array(pairs).T
                summary[m] = dict(intervals=len(pairs), bias=float(errors.mean()),
                                  rmse=float(np.sqrt(np.mean(errors**2))),
                                  mean_se=float(se.mean()),
                                  covered=int(np.count_nonzero(np.abs(errors) <= 1.96*se)))
        report.append(dict(name=name, target={m: target[m] for m in metrics}, sigma_log=sigma,
                           trials=trials, oracle_status=status, conditional_metrics=summary,
                           automatic_proposal_counts=auto_counts,
                           automatic_single_apex_rmse=(float(np.sqrt(np.mean(np.square(auto_matches))))
                                                       if auto_matches else None)))
    output.write_text(json.dumps(report, indent=2), encoding="utf-8")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--trials", type=int, default=128)
    args = parser.parse_args()
    run(args.trials, args.output)
