"""Compare fixed hyperbolic fits in nested windows at the same TripleTOF anchor."""

import argparse
import csv
import json
import math
from collections import Counter
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.stats import spearmanr


STEP = .35
SOURCE_HEADS = {
    "segments": "69f0bc653252b09df9c254b33fd87884eb5ad583",
    "areas": "35fa8de09bb30ed0ee8ad94de8f45a1d845c888c",
}
FEATURES = ("tail9", "log_rms9", "edge9")


def key(record):
    return int(record["scan"]), int(record["center"]), record["cohort"], int(record["scale"])


def positive(value):
    return value is not None and math.isfinite(value) and value > 0


def diagnostics(record):
    intercept, slope, curvature, skew = record["coefficients"]
    scale = record["scale"]
    grid = (np.arange(len(record["raw"])) - scale) * STEP
    shape = slope * grid - curvature * np.cosh(grid) - skew * np.sinh(grid)
    raw = np.asarray(record["raw"])
    log_rms = float(np.sqrt(np.mean((np.log(raw) - intercept - shape) ** 2)))
    edge = float(np.mean(raw[[0, -1]]) / np.max(raw))

    def log_shape(u):
        return slope * u - curvature * math.cosh(u) - skew * math.sinh(u)

    maximum = float(np.max(shape))
    weights = np.exp(shape - maximum)
    amplitude = (raw @ weights) / (weights @ weights)
    common_area = amplitude * quad(lambda u: math.exp(log_shape(u) - maximum),
                                   -4 * STEP, 4 * STEP, epsrel=1e-10)[0]
    return log_rms, edge, float(common_area)


def pair_records(segments, areas, injection):
    records = {}
    for line in (segments / f"audit{injection}.jsonl").read_text().splitlines():
        record = json.loads(line)
        if record["type"] == "window" and record["scale"] in (4, 8):
            records[key(record)] = record
    with (areas / f"areas{injection}.csv").open(newline="") as stream:
        fits = {key(row): row for row in csv.DictReader(stream)}

    rows = []
    for scan, center, cohort, scale in sorted(records):
        short_key, long_key = (scan, center, cohort, 4), (scan, center, cohort, 8)
        if scale != 4 or long_key not in records:
            continue
        short, long = records[short_key], records[long_key]
        row = dict(injection=injection, scan=scan, center=center, cohort=cohort,
                   short_reason=short["reason"], long_reason=long["reason"],
                   state="missing_matched_fit", area9=None, area17=None, change=None,
                   abs_log_change=None, common_change=None, tail9=None, tail17=None,
                   log_rms9=None, log_rms17=None, edge9=None, edge17=None)
        rows.append(row)
        if short_key not in fits or long_key not in fits:
            continue
        for field in ("raw", "mz"):
            if short[field] != long[field][4:-4]:
                raise ValueError(f"Windows are not nested at {injection}, {short_key}, {field}")
        common_areas = []
        for n, record, index in ((9, short, short_key), (17, long, long_key)):
            fit = fits[index]
            area = record["area_u"]
            csv_area = float(fit["fitted_whole_bessel_u"]) if fit["fitted_whole_bessel_u"] else None
            if area != csv_area:
                raise ValueError(f"Area evidence differs at {injection}, {index}")
            row[f"area{n}"] = area
            row[f"tail{n}"] = float(fit["tail_fraction"]) if fit["tail_fraction"] else None
            row[f"log_rms{n}"], row[f"edge{n}"], common_area = diagnostics(record)
            common_areas.append(common_area)
        row["common_change"] = common_areas[1] / common_areas[0] - 1
        has_short, has_long = positive(row["area9"]), positive(row["area17"])
        row["state"] = {(True, True): "both_finite", (True, False): "short_only",
                        (False, True): "long_only", (False, False): "neither"}[has_short, has_long]
        if has_short and has_long:
            row["change"] = row["area17"] / row["area9"] - 1
            row["abs_log_change"] = abs(math.log(row["area17"]) - math.log(row["area9"]))
    return rows


def correlation(x, y):
    if len(x) < 3 or len(set(x)) < 2 or len(set(y)) < 2:
        return None
    return float(spearmanr(x, y).statistic)


def describe(rows):
    finite = [row for row in rows if row["state"] == "both_finite"]
    common = [abs(row["common_change"]) for row in rows if row["common_change"] is not None]
    return dict(anchors=len(rows), states=dict(Counter(row["state"] for row in rows)),
                signed_change_p10_p50_p90=np.quantile([r["change"] for r in finite], [.1, .5, .9]).tolist() if finite else None,
                above_10_percent=sum(abs(r["change"]) > .10 for r in finite),
                median_abs_common_change=float(np.median(common)) if common else None)


def evaluate(training, evaluation):
    def arrays(rows):
        complete = [r for r in rows if r["state"] == "both_finite" and all(r[f] is not None for f in FEATURES)]
        return np.array([[r[f] for f in FEATURES] for r in complete]), np.array([r["abs_log_change"] for r in complete])

    train_x, train_y = arrays(training)
    test_x, test_y = arrays(evaluation)
    if min(len(train_y), len(test_y)) <= len(FEATURES):
        return dict(training=len(train_y), evaluation=len(test_y), status="insufficient_data")
    mean, std = train_x.mean(axis=0), train_x.std(axis=0)
    std[std == 0] = 1
    model = np.linalg.lstsq((train_x - mean) / std, train_y - train_y.mean(), rcond=None)[0]
    combined = (test_x - mean) / std @ model
    return dict(training=len(train_y), evaluation=len(test_y),
                individual_spearman={feature: correlation(test_x[:, i], test_y)
                                     for i, feature in enumerate(FEATURES)},
                combined_spearman=correlation(combined, test_y))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for arg in ("segments", "areas", "output"):
        parser.add_argument(arg, type=Path)
    args = parser.parse_args()
    for source, expected in SOURCE_HEADS.items():
        if (getattr(args, source) / "head.txt").read_text().strip() != expected:
            raise ValueError(f"Unexpected {source} evidence commit")
    rows = [row for injection in (1, 19) for row in pair_records(args.segments, args.areas, injection)]
    unique = {}
    for row in rows:
        index = row["injection"], row["scan"], row["center"]
        previous = unique.setdefault(index, row)
        if any(previous[field] != value for field, value in row.items() if field != "cohort"):
            raise ValueError(f"Cohorts disagree at {index}")
    groups = {injection: [r for r in unique.values() if r["injection"] == injection] for injection in (1, 19)}
    summary = dict(sources=SOURCE_HEADS, cohort_records=len(rows),
                   pooled=describe(list(unique.values())),
                   by_cohort={cohort: describe([r for r in rows if r["cohort"] == cohort])
                              for cohort in sorted({r["cohort"] for r in rows})},
                   injections={i: describe(group) for i, group in groups.items()},
                   train1_test19=evaluate(groups[1], groups[19]),
                   train19_test1=evaluate(groups[19], groups[1]))
    args.output.mkdir(parents=True, exist_ok=True)
    with (args.output / "pairs.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    (args.output / "summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
