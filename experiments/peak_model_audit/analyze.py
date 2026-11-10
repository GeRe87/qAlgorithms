"""Exploratory audit: QR-projected iid reference, never a raw-noise diagnosis."""
import argparse
import collections
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

STEP = 0.35
SEED = 617281


def geometry(scale, simulations=10000):
    u = STEP * np.arange(-scale, scale + 1, dtype=float)
    design = np.column_stack((np.ones(u.size), u, -np.cosh(u), -np.sinh(u)))
    q, r = np.linalg.qr(design, mode="reduced")
    residual_operator = np.eye(u.size) - q @ q.T
    noise = np.random.default_rng(SEED + scale).normal(size=(simulations, u.size))
    projected = noise @ residual_operator
    energy = np.sum(projected * projected, axis=1)
    lag = np.sum(projected[:, :-1] * projected[:, 1:], axis=1) / energy
    dw = np.sum(np.diff(projected, axis=1) ** 2, axis=1) / energy
    return design, q, r, residual_operator, lag, dw


def audit(path, output, source):
    cache = {scale: geometry(scale) for scale in (4, 8)}
    groups = {}
    scans = collections.Counter()
    scan_samples = collections.Counter()
    rows = []
    widths = {}
    maximum_coefficient_error = 0.0
    for line in path.read_text().splitlines():
        record = json.loads(line)
        if record["type"] == "scan":
            scans[record["status"]] += 1
            for field in ("points", "zeros", "negatives", "nonfinite"):
                scan_samples[field] += record.get(field, 0)
            assert record.get("candidates", 0) == 0, "Exploratory audit must not grant detections"
            continue
        key = record["cohort"] + "_n" + str(2 * record["scale"] + 1)
        group = groups.setdefault(key, {
            "attempted": 0, "reasons": collections.Counter(), "statuses": collections.Counter(),
            "eligible": 0, "integrable": 0, "apex_supported": 0, "refit_valid": 0,
            "area_se_available": 0, "numeric_zero_energy": 0, "outside_projection_reference": 0,
            "lag_values": [], "relative_area_se": [], "scan_outside": collections.defaultdict(list),
        })
        group["attempted"] += 1
        group["reasons"][record["reason"]] += 1
        if record["reason"] != "eligible":
            continue
        group["eligible"] += 1
        assert record["decision"] == 0, "Data-selected diagnostics cannot become prespecified decisions"
        group["statuses"][record["inference_status"]] += 1
        for field in ("integrable", "apex_supported", "refit_valid"):
            group[field] += record[field]
        y = np.asarray(record["logs"], dtype=float)
        design, q, r, residual_operator, reference_lag, reference_dw = cache[record["scale"]]
        beta = np.linalg.solve(r, q.T @ y)
        actual = np.asarray(record["coefficients"])
        error = float(np.max(np.abs(beta - actual) / (1 + np.abs(beta))))
        maximum_coefficient_error = max(maximum_coefficient_error, error)
        if error > 2e-7:
            raise AssertionError(f"Kernel vs independent QR coefficient error {error}")
        residual = residual_operator @ y
        rss = float(residual @ residual)
        floor = 128 * np.finfo(float).eps ** 2 * max(1, float(y @ y))
        row = {"scan": record["scan"], "cohort": record["cohort"], "center": record["center"],
               "n": y.size, "rss": rss, "spacing_cv": record["spacing_cv"],
               "status": record["inference_status"], "integrable": record["integrable"],
               "apex_supported": record["apex_supported"], "area_u": record["area_u"],
               "area_se_u": record["area_se_u"], "lag1": None, "dw": None}
        if rss <= floor:
            group["numeric_zero_energy"] += 1
        else:
            lag = float(residual[:-1] @ residual[1:] / rss)
            dw = float(np.diff(residual) @ np.diff(residual) / rss)
            lo, hi = np.quantile(reference_lag, [.025, .975])
            outside = bool(lag < lo or lag > hi)
            group["lag_values"].append(lag)
            group["outside_projection_reference"] += outside
            group["scan_outside"][record["scan"]].append(outside)
            row.update(lag1=lag, dw=dw)
        area, se = record["area_u"], record["area_se_u"]
        if area is not None and area > 0:
            widths.setdefault((record["scan"], record["cohort"], record["center"]), {})[record["scale"]] = area
            if se is not None and se >= 0:
                group["area_se_available"] += 1
                group["relative_area_se"].append(se / area)
        rows.append(row)
    for key, group in groups.items():
        for name in ("lag_values", "relative_area_se"):
            values = group.pop(name)
            group[name + "_median"] = float(np.median(values)) if values else None
        fractions = [np.mean(values) for values in group.pop("scan_outside").values()]
        group["per_scan_outside_fraction_median"] = float(np.median(fractions)) if fractions else None
    reference = {}
    for scale, (_, _, _, operator, lag, dw) in cache.items():
        reference[str(2 * scale + 1)] = {
            "residual_df": int(round(np.trace(operator))),
            "lag_95_reference": np.quantile(lag, [.025, .975]).tolist(),
            "lag_median": float(np.median(lag)),
            "dw_95_reference": np.quantile(dw, [.025, .975]).tolist(),
        }
    sensitivity = {}
    for cohort in ("fixed_index", "observed_maximum"):
        ratios = [value[8] / value[4] for key, value in widths.items()
                  if key[1] == cohort and 4 in value and 8 in value]
        sensitivity[cohort] = {"paired_windows": len(ratios),
                               "area17_over_area9_quantiles": np.quantile(ratios, [.1, .5, .9]).tolist() if ratios else []}
    summary = {
        "source": str(source.name), "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "source_md5": hashlib.md5(source.read_bytes()).hexdigest(), "numpy_version": np.__version__,
        "step": STEP, "seed": SEED, "simulations_per_geometry": 10000,
        "interpretation": "descriptive fixed-window projected-iid reference, not selection-adjusted p-values or raw-noise correlation",
        "scans": dict(scans), "scan_samples": dict(scan_samples), "cohorts": groups, "projected_iid_reference": reference,
        "width_sensitivity": sensitivity, "max_kernel_qr_relative_coefficient_error": maximum_coefficient_error,
    }
    output.mkdir(parents=True, exist_ok=True)
    (output / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    if rows:
        with (output / "residuals.csv").open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
            writer.writeheader()
            writer.writerows(rows)
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("records", type=Path)
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    summary = audit(args.records, args.output, args.source)
    print(json.dumps({"scans": summary["scans"], "cohorts": {
        key: {field: value for field, value in group.items() if field not in ("reasons", "statuses")}
        for key, group in summary["cohorts"].items()}}, indent=2))
