"""Paired central-window fits against known whole-line target areas."""

import argparse
import csv
import itertools
import json
import math
from pathlib import Path
import subprocess

import numpy as np
from scipy.integrate import quad

TARGET = 1000.
SEED = 74821
REPLICATES = 40


def gaussian(u, center, width):
    return np.exp(-.5 * ((u - center) / width) ** 2) / (width * math.sqrt(2 * math.pi))


def hyper_shape(u):
    return np.exp(.7 * u - 4 * np.cosh(u) - .4 * np.sinh(u))


HYPER_NORM = quad(hyper_shape, -20, 20, epsabs=1e-13)[0]


def target(u, shape):
    if shape == "hyperbolic":
        return TARGET * hyper_shape(u) / HYPER_NORM
    if shape == "gaussian":
        return TARGET * gaussian(u, 0, .5)
    return TARGET * (.8 * gaussian(u, 0, .45) + .2 * gaussian(u, .65, .9))


def trials():
    rng = np.random.default_rng(SEED)
    for shape, scale, step, noise, background, overlap in itertools.product(
            ("hyperbolic", "gaussian", "asymmetric_mixture"), (4, 8), (.175, .35),
            ("none", "iid_04", "iid_12", "ar1_12"), (0., .02), (0., .15)):
        u = step * np.arange(-scale, scale + 1)
        signal = target(u, shape)
        mean = signal + background * float(target(np.array([0.]), shape)[0])
        mean = mean + overlap * TARGET * gaussian(u, .9, .4)
        sigma = {"none": 0., "iid_04": .04, "iid_12": .12, "ar1_12": .12}[noise]
        for repeat in range(REPLICATES if sigma else 1):
            z = rng.normal(size=u.size)
            if noise == "ar1_12":
                for i in range(1, u.size):
                    z[i] = .65 * z[i - 1] + math.sqrt(1 - .65 ** 2) * z[i]
            raw = (mean * np.exp(sigma * z - .5 * sigma ** 2)).astype(np.float32)
            yield dict(shape=shape, samples=u.size, step=step, noise=noise,
                       background=background, overlap=overlap, replicate=repeat), raw


def fit_metrics(model, coefficients, half_width, noisy):
    intercept, slope, left, right, area, se = coefficients
    finite = math.isfinite(area) and area > 0
    if not finite:
        return dict(area=None, relative_error=None, tail_fraction=None, native_se=None, se_scale=None, covered=None)

    def profile(u):
        if model == "hyperbolic":
            return math.exp(intercept + slope * u - left * math.cosh(u) - right * math.sinh(u))
        return math.exp(intercept + slope * u + (left if u < 0 else right) * u * u)

    window = quad(profile, -half_width, half_width, points=[0], epsabs=1e-8)[0]
    tail = 1 - window / area
    if not -1e-7 <= tail <= 1:
        raise ArithmeticError(f"Inconsistent window and total area for {model}: {tail}")
    valid_se = noisy and math.isfinite(se) and se > 0
    log_scale = model == "piecewise_quadratic"
    distance = abs(math.log(area / TARGET)) if log_scale else abs(area - TARGET)
    return dict(area=area, relative_error=area / TARGET - 1, tail_fraction=max(0., tail),
                native_se=se if valid_se else None, se_scale=("log_area" if log_scale else "area") if valid_se else None,
                covered=distance <= 1.96 * se if valid_se else None)


def summarize(rows):
    errors = [r["relative_error"] for r in rows if r["relative_error"] is not None]
    covered = [r["covered"] for r in rows if r["covered"] is not None]
    tails = [r["tail_fraction"] for r in rows if r["tail_fraction"] is not None]
    return dict(attempted=len(rows), finite=len(errors), coverage_count=len(covered),
                coverage=float(np.mean(covered)) if covered else None,
                mean_bias=float(np.mean(errors)) if errors else None,
                median_absolute_error=float(np.median(np.abs(errors))) if errors else None,
                error_p10_p50_p90=np.quantile(errors, [.1, .5, .9]).tolist() if errors else None,
                tail_p50_p90=np.quantile(tails, [.5, .9]).tolist() if tails else None,
                fraction_absolute_error_gt_10pct=float(np.mean(np.abs(errors) > .1)) if errors else None)


def run(executable, output):
    cases = list(trials())
    payload = "".join(f'{(case["samples"]-1)//2} {case["step"]} ' + " ".join(map(str, raw)) + "\n"
                      for case, raw in cases)
    result = subprocess.run([str(executable.resolve())], input=payload, capture_output=True, text=True, check=True)
    lines = result.stdout.splitlines()
    if len(lines) != len(cases):
        raise RuntimeError("Native output does not match the input trial count")
    rows, pairs = [], []
    for (case, _), line in zip(cases, lines):
        values = list(map(float, line.split()))
        if len(values) != 12:
            raise RuntimeError("Invalid native fit record")
        pair = []
        for model, coefficients in (("piecewise_quadratic", values[:6]), ("hyperbolic", values[6:])):
            metrics = fit_metrics(model, coefficients, (case["samples"] - 1) / 2 * case["step"], case["noise"] != "none")
            rows.append(dict(**case, model=model, **metrics))
            pair.append(metrics["relative_error"])
            matched_model = {"gaussian": "piecewise_quadratic", "hyperbolic": "hyperbolic"}.get(case["shape"])
            if model == matched_model and case["noise"] == "none" and case["background"] == case["overlap"] == 0:
                if metrics["relative_error"] is None or abs(metrics["relative_error"]) > 1e-5:
                    raise AssertionError(f"Native matched-shape recovery failed for {case}")
        pairs.append(dict(**case, piecewise_error=pair[0], hyperbolic_error=pair[1]))
    output.mkdir(parents=True, exist_ok=True)
    with (output / "fits.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    groups = {}
    for row in rows:
        key = (row["shape"], row["samples"], row["step"], row["noise"], row["background"], row["overlap"], row["model"])
        groups.setdefault(key, []).append(row)
    summary = [dict(zip(("shape", "samples", "step", "noise", "background", "overlap", "model"), key),
                    **summarize(group)) for key, group in groups.items()]
    diagnostics = {}
    for model in ("piecewise_quadratic", "hyperbolic"):
        subset = [r for r in rows if r["model"] == model and r["tail_fraction"] is not None and r["noise"] != "none"]
        diagnostics[model] = {label: summarize([r for r in subset if lo <= r["tail_fraction"] < hi])
                              for label, lo, hi in [("below_5pct", 0, .05), ("5_to_20pct", .05, .2),
                                                    ("20pct_or_more", .2, 1.01)]}
    both = [p for p in pairs if p["piecewise_error"] is not None and p["hyperbolic_error"] is not None]
    report = dict(seed=SEED, target_area=TARGET, noisy_replicates=REPLICATES, paired_windows=len(cases),
                  both_finite=len(both), diagnostics=diagnostics, conditions=summary,
                  hyperbolic_lower_absolute_error=sum(abs(p["hyperbolic_error"]) < abs(p["piecewise_error"]) for p in both))
    (output / "summary.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    run(args.executable, args.output)
