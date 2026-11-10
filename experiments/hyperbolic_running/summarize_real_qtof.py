"""Summarize paired real-QTOF window fits. Descriptive, not mass accuracy."""
import argparse
import csv
import json
import math
import statistics
from collections import defaultdict
from pathlib import Path


def finite(x):
    try:
        value = float(x)
        return value if math.isfinite(value) else None
    except (ValueError, TypeError):
        return None


def mean(values):
    v = [x for x in values if x is not None]
    return statistics.fmean(v) if v else None


def summarise(input_csv):
    with Path(input_csv).open(encoding="utf-8", newline="") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        raise ValueError("No eligible measured profile windows")
    groups = defaultdict(list)
    for r in rows:
        groups[(int(r["half_width"]), int(r["trace_count"]))].append(r)
    out = {
        "dataset_kind": "experimental profile-mode TripleTOF, no known concentration/area truth",
        "area_metric": "complete model-tail area in dimensionless u, not calibrated m/z area",
        "sampling": "deterministic time-stratified spectra, no resampling",
        "total_evaluated_trace_windows": len(rows),
        "unique_scans": len({r["scan"] for r in rows}),
        "groups": [],
    }
    for (h, m), group in sorted(groups.items()):
        def values(field):
            return [finite(r[field]) for r in group]
        raw_pre_old=values("old_raw_rss_pre")
        raw_post_old=values("old_raw_rss_post")
        raw_pre_new=values("new_raw_rss_pre")
        raw_post_new=values("new_raw_rss_post")
        a_old=values("old_area_u_post")
        a_new=values("new_area_u_post")
        ratios = [new/old for old, new in zip(a_old, a_new)
                  if old is not None and new is not None and old > 0 and new > 0]
        oldwins=0
        newwins=0
        for row in group:
            aicold=finite(row["old_AICc"])
            aicnew=finite(row["new_AICc"])
            if aicold is not None and aicnew is not None:
                oldwins+=aicold < aicnew
                newwins+=aicnew < aicold
        out["groups"].append({
            "half_window_samples": h,
            "pooled_traces": m,
            "rows": len(group),
            "old_integrable_fraction": mean(values("old_integrable")),
            "new_integrable_fraction": mean(values("new_integrable")),
            "valid_amplitude_refit_fraction": mean(values("raw_refit_valid")),
            "old_mean_raw_rss_before": mean(raw_pre_old),
            "old_mean_raw_rss_after": mean(raw_post_old),
            "new_mean_raw_rss_before": mean(raw_pre_new),
            "new_mean_raw_rss_after": mean(raw_post_new),
            "mean_log_rss_old": mean(values("old_log_rss")),
            "mean_log_rss_new": mean(values("new_log_rss")),
            "aicc_preferred_old": oldwins,
            "aicc_preferred_new": newwins,
            "median_new_to_old_postfit_full_area_ratio": statistics.median(ratios) if ratios else None,
            "count_finite_model_area_pairs": len(ratios),
        })
    return out


def main():
    p=argparse.ArgumentParser()
    p.add_argument("csv", type=Path)
    p.add_argument("output", type=Path)
    args=p.parse_args()
    results=summarise(args.csv)
    args.output.write_text(json.dumps(results, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
