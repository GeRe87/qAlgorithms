"""Audit matched TripleTOF windows using dimensionless index coordinates."""

import argparse
import csv
import json
import math
import warnings
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy import integrate


STEP = 0.35
RELATIVE_AREA_TOLERANCE = 2e-8
BASE_EVIDENCE_HEAD = "69f0bc653252b09df9c254b33fd87884eb5ad583"


def quadrature(function, lower, upper):
    with warnings.catch_warnings():
        warnings.simplefilter("error", integrate.IntegrationWarning)
        value, error = integrate.quad(function, lower, upper, epsabs=1e-12, epsrel=1e-10, limit=200)
    if not math.isfinite(value) or value < 0:
        raise ArithmeticError("Invalid quadrature result")
    return value, error


def fitted_areas(raw, scale, slope, curvature, skew):
    """Reconstruct the scalar raw LS refit and integrate that fixed shape."""
    if len(raw) != 2 * scale + 1 or not all(math.isfinite(y) and y > 0 for y in raw):
        raise ValueError("Expected the positive observations of the matched window")

    def log_shape(u):
        return slope * u - curvature * math.cosh(u) - skew * math.sinh(u)

    grid = [(i - scale) * STEP for i in range(len(raw))]
    shapes = [log_shape(u) for u in grid]
    maximum = max(shapes)
    weights = [math.exp(value - maximum) for value in shapes]
    amplitude = math.fsum(y * w for y, w in zip(raw, weights)) / math.fsum(w * w for w in weights)
    observed = STEP * (math.fsum(raw[1:-1]) + .5 * (raw[0] + raw[-1]))
    sampled = STEP * amplitude * (math.fsum(weights[1:-1]) + .5 * (weights[0] + weights[-1]))
    lower, upper = grid[0], grid[-1]

    if curvature <= abs(skew):
        # The observed-window integral remains meaningful for this fit.
        window, _ = quadrature(lambda u: math.exp(log_shape(u) - maximum), lower, upper)
        return observed, sampled, amplitude * window, None, None, None

    p = .5 * (curvature + skew)
    q = .5 * (curvature - skew)
    apex = math.asinh(slope / (2 * math.sqrt(p) * math.sqrt(q))) + .5 * (math.log(q) - math.log(p))
    log_p, log_q = math.log(p) + apex, math.log(q) - apex
    p_at_apex, q_at_apex = math.exp(log_p), math.exp(log_q)
    width = 1 / math.sqrt(p_at_apex + q_at_apex)

    def normalized(t):
        offset = width * t
        try:
            exponent = (slope * offset + p_at_apex + q_at_apex
                        - math.exp(log_p + offset) - math.exp(log_q - offset))
        except OverflowError:
            return 0.0
        return math.exp(exponent)

    def integrate_piece(lo, hi):
        if lo < 0 < hi:
            left, left_error = quadrature(normalized, lo, 0)
            right, right_error = quadrature(normalized, 0, hi)
            return left + right, left_error + right_error
        return quadrature(normalized, lo, hi)

    lo, hi = (lower - apex) / width, (upper - apex) / width
    left, left_error = integrate_piece(-math.inf, lo)
    middle, middle_error = integrate_piece(lo, hi)
    right, right_error = integrate_piece(hi, math.inf)
    total = left + middle + right
    log_multiplier = math.log(amplitude) + log_shape(apex) - maximum + math.log(width)
    whole = math.exp(log_multiplier + math.log(total))
    window = math.exp(log_multiplier + math.log(middle)) if middle > 0 else 0.0
    return observed, sampled, window, whole, (left + right) / total, (left_error + middle_error + right_error) / total


def audit(records, windows):
    with windows.open(newline="") as stream:
        fitted = {(int(row["scan"]), int(row["center"]), int(row["scale"]))
                  for row in csv.DictReader(stream)}
    rows = []
    with records.open() as stream:
        for line in stream:
            record = json.loads(line)
            if record["type"] != "window" or record["reason"] != "eligible":
                continue
            key = (record["scan"], record["center"], record["scale"])
            if key not in fitted:
                continue
            _, slope, curvature, skew = record["coefficients"]
            observed, sampled, window, whole, tails, error = fitted_areas(
                record["raw"], key[2], slope, curvature, skew)
            analytic = record["area_u"]
            difference = abs(whole / analytic - 1) if whole is not None and analytic is not None and analytic > 0 else None
            if difference is not None and difference > RELATIVE_AREA_TOLERANCE:
                raise ArithmeticError(f"Whole-line area mismatch at {key}: {difference:.3g}")
            rows.append(dict(scan=key[0], center=key[1], scale=key[2], cohort=record["cohort"],
                             observed_trapezoid_u=observed, sampled_fit_trapezoid_u=sampled,
                             fitted_window_u=window, fitted_whole_quadrature_u=whole,
                             fitted_whole_bessel_u=analytic, tail_fraction=tails,
                             quadrature_relative_error_estimate=error, whole_relative_difference=difference,
                             window_vs_observed=window / observed - 1,
                             sampled_fit_vs_observed=sampled / observed - 1))
    if not rows:
        raise ValueError("No matched eligible windows")
    return rows


def summarize(rows):
    groups = defaultdict(list)
    for row in rows:
        groups[(row["cohort"], row["scale"])].append(row)
    metrics = ("window_vs_observed", "sampled_fit_vs_observed", "tail_fraction",
               "whole_relative_difference", "quadrature_relative_error_estimate")
    output = {}
    for (cohort, scale), group in sorted(groups.items()):
        summary = {"records": len(group)}
        for key in metrics:
            values = [row[key] for row in group if row[key] is not None]
            summary[key] = {"count": len(values), "p10_p50_p90": np.quantile(values, [.1, .5, .9]).tolist(),
                            "max": max(values)} if values else {"count": 0}
        output[f"{cohort}_n{2 * scale + 1}"] = summary
    return dict(records=len(rows), distinct_windows=len({(r["scan"], r["center"], r["scale"]) for r in rows}),
                whole_comparisons=sum(r["whole_relative_difference"] is not None for r in rows),
                max_whole_relative_difference=max((r["whole_relative_difference"] for r in rows
                                                   if r["whole_relative_difference"] is not None), default=None),
                groups=output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("evidence", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    if (args.evidence / "head.txt").read_text().strip() != BASE_EVIDENCE_HEAD:
        raise ValueError("Unexpected segmentation evidence commit")
    args.output.mkdir(parents=True, exist_ok=True)
    summaries = {}
    for injection, expected in [(1, 933), (19, 1094)]:
        result = audit(args.evidence / f"audit{injection}.jsonl", args.evidence / f"windows{injection}.csv")
        if len(result) != expected:
            raise ValueError(f"Injection {injection} has {len(result)} matched records, expected {expected}")
        with (args.output / f"areas{injection}.csv").open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=result[0].keys())
            writer.writeheader()
            writer.writerows(result)
        summaries[str(injection)] = summarize(result)
    report = dict(source_head=BASE_EVIDENCE_HEAD, step=STEP, relative_area_tolerance=RELATIVE_AREA_TOLERANCE,
                  input_manifest=json.loads((args.evidence / "input-manifest.json").read_text()), injections=summaries)
    (args.output / "summary.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
